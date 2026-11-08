#!/usr/bin/env python3
"""Zn-Br project structure generator (charter section 2).
Generates single-species and complex POSCARs per the v1.0 charter.
Usage: python3 generate_structures.py <project_calc_root>"""
import math, os, sys, csv
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
BOX = 32.0                      # uniform box (rev 2026-10-29: 44->32 A for memory; ~7 A vacuum)
DISTS = [3.0, 3.4, 3.8]         # N...Br approach distances (charter section 2.3)
PP_ORDER = ['N', 'C', 'H', 'Br', 'O']   # species order convention per case subset

CATIONS = {
    'C10TAB': 'CCCCCCCCCC[N+](C)(C)C',
    'C12TAB': 'CCCCCCCCCCCC[N+](C)(C)C',
    'C14TAB': 'CCCCCCCCCCCCCC[N+](C)(C)C',
    'MEPBr':  'CC[N+]1(CCCC1)C',
    'TPABr':  'CCC[N+](CCC)(CCC)CCC',
    'EMImBr': 'CCn1cc[n+](C)c1',
}
EXPECT_NATOMS = {'C10TAB': 44, 'C12TAB': 50, 'C14TAB': 56, 'MEPBr': 24, 'TPABr': 41, 'EMImBr': 19}


def write_poscar(path, title, species_counts, cart):
    """cart: list of (symbol, [x,y,z]); species_counts: list of (sym, n) in POSCAR order."""
    with open(path, 'w', newline='\n') as f:
        f.write(title + '\n1.0\n')
        for v in ('%.1f 0.0 0.0' % BOX, '0.0 %.1f 0.0' % BOX, '0.0 0.0 %.1f' % BOX):
            f.write('  ' + v + '\n')
        f.write('  ' + ' '.join(s for s, _ in species_counts) + '\n')
        f.write('  ' + ' '.join(str(n) for _, n in species_counts) + '\n')
        f.write('Direct\n')
        for sym, p in cart:
            fx = (p[0] / BOX) % 1.0
            fy = (p[1] / BOX) % 1.0
            fz = (p[2] / BOX) % 1.0
            f.write(f'  {fx:.10f}  {fy:.10f}  {fz:.10f}\n')


def group_cart(cart):
    """Group cart (list of (sym, pos)) into species_counts ordered by PP_ORDER subset."""
    order = [s for s in PP_ORDER if any(sym == s for sym, _ in cart)]
    counts = [(s, sum(1 for sym, _ in cart if sym == s)) for s in order]
    return order, counts


def centered(cart):
    """Return cart shifted so molecular centroid sits at box center (22,22,22)."""
    n = len(cart)
    c = [sum(p[i] for _, p in cart) / n for i in range(3)]
    return [(s, [p[i] - c[i] + BOX / 2 for i in range(3)]) for s, p in cart]


def write_xyz(path, title, cart):
    with open(path, 'w', newline='\n') as f:
        f.write(f'{len(cart)}\n{title}\n')
        for s, p in cart:
            f.write(f'{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n')


def cation_geometry(smiles, natoms_expect):
    """ETKDG 24 conformers -> MMFF94 -> lowest. Returns (symbols, positions, n_pos)."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f'SMILES parse failed: {smiles}'
    mol = Chem.AddHs(Chem.Mol(mol), False)
    assert len(mol.GetAtoms()) == natoms_expect, f'atom count {len(mol.GetAtoms())} != {natoms_expect}'
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xf00d
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=24, params=params)
    assert len(cids) >= 20, f'only {len(cids)} conformers'
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    best = min(range(len(cids)), key=lambda i: res[i][1])
    conf = mol.GetConformer(best)
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    pos = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    # N index: the charged N
    n_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge() == 1 and a.GetSymbol() == 'N'][0]
    return sym, pos, n_idx


def anion_geometry(kind):
    """Return (symbols, positions) for the anion/neutral per charter section 2.2."""
    if kind == 'Br':
        return ['Br'], [[0.0, 0.0, 0.0]]
    if kind == 'Br3':
        # linear, 2.55 A spacing; index 0 = end Br facing N+
        return ['Br'] * 3, [[0.0, 0.0, 0.0], [2.55, 0.0, 0.0], [5.10, 0.0, 0.0]]
    if kind == 'Br2':
        return ['Br'] * 2, [[0.0, 0.0, 0.0], [2.28, 0.0, 0.0]]
    if kind == 'H2O':
        # O at origin, H at 0.97 A, angle 104.5 deg
        ang = math.radians(104.5 / 2)
        return ['O', 'H', 'H'], [[0.0, 0.0, 0.0],
                                 [0.97 * math.cos(ang), 0.97 * math.sin(ang), 0.0],
                                 [0.97 * math.cos(ang), -0.97 * math.sin(ang), 0.0]]
    if kind == 'Br5':
        # V-shaped Br5- per charter 2.2: (Br2)2...Br- = central Br- (vertex, faces N+)
        # with two Br2 arms. Middle contact M...X = 3.1 A, terminal Br-Br bond
        # X-Y = 2.4 A, V opening angle 115 deg; arms in the y-z plane, opening
        # along +z (away from the cation when placed).
        half = math.radians(115.0) / 2
        s, c = math.sin(half), math.cos(half)
        P = []
        for sgn in (1.0, -1.0):
            x = [0.0, sgn * 3.1 * s, 3.1 * c]
            y = [0.0, x[1] + sgn * 2.4 * s, x[2] + 2.4 * c]
            P.append(x)
            P.append(y)
        # order: M(vertex, contact atom), X1, Y1, X2, Y2
        return ['Br'] * 5, [[0.0, 0.0, 0.0]] + P
    raise KeyError(kind)


def place_complex(cat_sym, cat_pos, n_idx, an_kind, dist, h_flip=False):
    """Place anion so contact atom A sits at N + dhat*dist for a collision-free direction dhat.
    Deterministic: sample 60 golden-spiral directions, score by min anion-cation contact, keep best."""
    Np = cat_pos[n_idx]

    def local_anion():
        """Anion in local frame: A at origin, main axis +z (away from cation), V/H plane in z-y."""
        if an_kind == 'Br':
            return [('Br', [0., 0., 0.])]
        if an_kind == 'Br3':
            return [('Br', [0., 0., 0.]), ('Br', [0., 0., 2.55]), ('Br', [0., 0., 5.10])]
        if an_kind == 'Br2':
            return [('Br', [0., 0., 0.]), ('Br', [0., 0., 2.28])]
        if an_kind == 'Br5':
            # V-shaped Br5-: vertex Br- at origin (contact atom, faces N+),
            # two Br2 arms in the z-y plane opening along +z:
            # M...X = 3.1 (middle contact), X-Y = 2.4 (terminal bond), angle 115 deg
            half = math.radians(115.0) / 2
            s5, c5 = math.sin(half), math.cos(half)
            out = [('Br', [0., 0., 0.])]
            for sgn in (1.0, -1.0):
                x = [0., sgn * 3.1 * s5, 3.1 * c5]
                y = [0., x[1] + sgn * 2.4 * s5, x[2] + 2.4 * c5]
                out.append(('Br', x))
                out.append(('Br', y))
            return out
        if an_kind == 'H2O':
            s = math.sin(math.radians(104.5 / 2)); c = math.cos(math.radians(104.5 / 2))
            sign = -1.0 if h_flip else 1.0   # H's toward cation if flipped
            return [('O', [0., 0., 0.]), ('H', [0.97 * s, 0., 0.97 * c * sign]), ('H', [-0.97 * s, 0., 0.97 * c * sign])]
        raise KeyError(an_kind)

    loc = local_anion()
    extra = 0.8 if (an_kind == 'H2O' and h_flip) else 0.0   # H-toward: keep H tips clear of cation H shell

    def world_for(dhat):
        ref = [0.0, 0.0, 1.0] if abs(dhat[2]) < 0.9 else [1.0, 0.0, 0.0]
        e1 = [ref[k] - sum(ref[j] * dhat[j] for j in range(3)) * dhat[k] for k in range(3)]
        e1l = math.sqrt(sum(x * x for x in e1)); e1 = [x / e1l for x in e1]
        e2 = [dhat[1] * e1[2] - dhat[2] * e1[1], dhat[2] * e1[0] - dhat[0] * e1[2], dhat[0] * e1[1] - dhat[1] * e1[0]]
        origin = [Np[k] + dhat[k] * (dist + extra) for k in range(3)]
        return [(s, [origin[k] + p[0] * e1[k] + p[1] * e2[k] + p[2] * dhat[k] for k in range(3)]) for s, p in loc]

    # golden-spiral direction sampling (deterministic)
    NC = 60
    golden = math.pi * (3.0 - math.sqrt(5.0))
    best, best_score = None, -1e9
    for k in range(NC):
        zz = 1.0 - 2.0 * (k + 0.5) / NC
        rr = math.sqrt(max(0.0, 1.0 - zz * zz))
        cand = [rr * math.cos(golden * k), rr * math.sin(golden * k), zz]
        w = world_for(cand)
        score = 9e9
        for _, ap in w:
            for cp in cat_pos:
                score = min(score, math.dist(ap, cp))
        if score > best_score:
            best_score, best = score, w
    return [s for s, _ in best], [p for _, p in best]


manifest = []
def emit(root_dir, title, cart, note='', charge=0):
    order, counts = group_cart(cart)
    cart_g = []
    for s in order:
        cart_g += [(sym, p) for sym, p in cart if sym == s]
    cart_c = centered(cart_g)
    os.makedirs(root_dir, exist_ok=True)
    write_poscar(os.path.join(root_dir, 'POSCAR'), title, counts, cart_c)
    write_xyz(os.path.join(root_dir, 'init.xyz'), title, cart_c)
    with open(os.path.join(root_dir, 'species_order.txt'), 'w', newline='\n') as f:
        f.write(' '.join(order) + '\n')
    with open(os.path.join(root_dir, 'charge.txt'), 'w', newline='\n') as f:
        f.write(f'{charge}\n')
    manifest.append({'case': os.path.relpath(root_dir, ROOT).replace('\\', '/'), 'natoms': len(cart_c),
                     'species': ' '.join(order), 'charge': charge, 'note': note})


# ---- 1. cations (best conformer each) ----
cat_geo = {}
for name, smi in CATIONS.items():
    sym, pos, n_idx = cation_geometry(smi, EXPECT_NATOMS[name])
    cat_geo[name] = (sym, pos, n_idx)
    cart = centered(list(zip(sym, pos)))
    emit(os.path.join(ROOT, '01_single', name), f'{name} cation best-conf', cart,
         'best of 24 ETKDG/MMFF94', charge=1)

# ---- 2. anion singles ----
for kind in ['Br', 'Br3', 'Br5', 'Br2', 'H2O']:
    sym, pos = anion_geometry(kind)
    cart = centered(list(zip(sym, pos)))
    d = 'H2O' if kind == 'H2O' else kind
    # charged anions: Br-, Br3-, Br5- carry -1 (NELECT handled by setup_case.sh)
    emit(os.path.join(ROOT, '01_single', d), f'{kind} single', cart, 'manual internal coords',
         charge=(-1 if kind in ('Br', 'Br3', 'Br5') else 0))

# ---- 3. complexes: 3 x 5 main + 4 controls, 3 confs each ----
MATRIX = {
    'C10TAB': ['Br', 'Br3', 'Br5', 'Br2', 'H2O'],
    'C12TAB': ['Br', 'Br3', 'Br5', 'Br2', 'H2O'],
    'C14TAB': ['Br', 'Br3', 'Br5', 'Br2', 'H2O'],
    'MEPBr':  ['Br3'],
    'TPABr':  ['Br3'],
    'EMImBr': ['Br3', 'H2O'],
}
for cat, anions in MATRIX.items():
    sym, pos, n_idx = cat_geo[cat]
    for an in anions:
        for ci in range(3):
            if an == 'H2O':
                dist = DISTS[ci] if ci < 2 else DISTS[0]
                h_flip = (ci == 2)
                note = f'H2O {"O-toward" if not h_flip else "H-toward"} d={dist + (0.8 if h_flip else 0.0)}'
            else:
                dist = DISTS[ci]
                h_flip = False
                note = f'N...Br={dist}'
            asym, apos = place_complex(sym, pos, n_idx, an, dist, h_flip)
            cart = centered(list(zip(sym, pos)) + list(zip(asym, apos)))
            cname = f'{cat}_{an}_c{ci+1}'
            # net charge = cation(+1) + anion: Br/Br3/Br5 are -1 -> 0; Br2/H2O neutral -> +1
            cplx_charge = 0 if an in ('Br', 'Br3', 'Br5') else 1
            emit(os.path.join(ROOT, '02_complex', cname), cname, cart, note, charge=cplx_charge)

# ---- manifest ----
try:
    _csvf = open(os.path.join(ROOT, '00_structures', 'cases_manifest.csv'), 'w', newline='\n', encoding='utf-8')
except PermissionError:
    _csvf = None
    print('WARN: cases_manifest.csv locked (open in editor); CSV not refreshed')
if _csvf:
    w = csv.DictWriter(_csvf, fieldnames=['case', 'natoms', 'species', 'charge', 'note'])
    w.writeheader()
    for m in manifest:
        w.writerow(m)
    _csvf.close()
print(f'generated {len(manifest)} cases; box={BOX} A uniform')
