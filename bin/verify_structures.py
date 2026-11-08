#!/usr/bin/env python3
"""Spot-check generated Zn-Br POSCARs: atom counts, internal distances, centering."""
import math, os, sys, csv

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
checked = 0
fails = []
EXPECT_SINGLE = {'C10TAB': 44, 'C12TAB': 50, 'C14TAB': 56, 'MEPBr': 24, 'TPABr': 41,
                 'EMImBr': 19, 'Br': 1, 'Br3': 3, 'Br5': 5, 'Br2': 2, 'H2O': 3}

def read_poscar(p):
    """Return (species, cartesian coords, box length).

    VASP convention: line 2 = universal scale (1.0 for our files), lines 3-5 =
    lattice vectors, line 6 = species, line 7 = counts, line 8 = Direct.
    """
    with open(p) as f:
        lines = f.read().splitlines()
    scale = float(lines[1].split()[0])
    lat = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    syms = lines[5].split()
    nums = [int(x) for x in lines[6].split()]
    species = []
    for s, n in zip(syms, nums):
        species += [s] * n
    cart = []
    for ln in lines[8:8 + sum(nums)]:
        f = [float(x) for x in ln.split()[:3]]
        cart.append([sum(f[k] * lat[k][j] for k in range(3)) * scale for j in range(3)])
    box = scale * lat[0][0]
    return species, cart, box


def check(name, species, cart, box):
    global checked
    checked += 1
    n = len(species)
    # centering
    cx = sum(p[0] for p in cart) / n
    centre = box / 2.0
    for k in range(3):
        c = sum(p[k] for p in cart) / n
        if abs(c - centre) > 1e-6:
            fails.append(f'{name}: centroid off {c:.4f} on axis {k}')
    nbr = species.count('Br')

    def dist(a, b):
        return math.dist(cart[a], cart[b])

    if '_c' not in name:
        if n != EXPECT_SINGLE.get(name, -1):
            fails.append(f'{name}: natoms {n} != {EXPECT_SINGLE.get(name)}')
        return
    # complexes: format CAT_AN_cK
    base = name.rsplit('_c', 1)[0]
    cat, an = base.split('_', 1)
    if an == 'Br':  # single Br anion: N...Br = 3.0/3.4/3.8
        iBr = len(species) - 1  # Br last in PP_ORDER
        dN = [dist(i, iBr) for i, s in enumerate(species) if s == 'N']
        d = min(dN)
        if not any(abs(x - t) < 0.02 for x in dN for t in (3.0, 3.4, 3.8)):
            fails.append(f'{name}: N...Br {d:.3f} not in {3.0,3.4,3.8}')
    elif an == 'Br3':
        idxs = [i for i, s in enumerate(species) if s == 'Br']
        d_ab = dist(idxs[0], idxs[1]); d_bc = dist(idxs[1], idxs[2])
        if abs(d_ab - 2.55) > 0.02 or abs(d_bc - 2.55) > 0.02:
            fails.append(f'{name}: Br3 spacing {d_ab:.3f}/{d_bc:.3f}')
        dNs = [dist(i, idxs[0]) for i, s in enumerate(species) if s == 'N']
        if not any(abs(x - t) < 0.02 for x in dNs for t in (3.0, 3.4, 3.8)):
            fails.append(f'{name}: N...Br3-end {min(dNs):.3f}')
    elif an == 'Br5':
        idxs = [i for i, s in enumerate(species) if s == 'Br']
        if len(idxs) != 5:
            fails.append(f'{name}: Br5 needs 5 Br, got {len(idxs)}')
        else:
            # order in file: M(vertex, faces N), X1, Y1, X2, Y2
            M, X1, Y1, X2, Y2 = idxs
            for lbl, dd, exp in (('M-X1', dist(M, X1), 3.10), ('X1-Y1', dist(X1, Y1), 2.40),
                                 ('M-X2', dist(M, X2), 3.10), ('X2-Y2', dist(X2, Y2), 2.40)):
                if abs(dd - exp) > 0.02:
                    fails.append(f'{name}: Br5 {lbl} {dd:.3f} != {exp}')
            dX1X2 = dist(X1, X2)
            if not (4.5 < dX1X2 < 6.0):
                fails.append(f'{name}: Br5 X1-X2 {dX1X2:.3f} out of V range')
        dNs = [dist(i, idxs[0]) for i, s in enumerate(species) if s == 'N']
        if not any(abs(x - t) < 0.02 for x in dNs for t in (3.0, 3.4, 3.8)):
            fails.append(f'{name}: N...Br5-vertex {min(dNs):.3f}')
    elif an == 'Br2':
        idxs = [i for i, s in enumerate(species) if s == 'Br']
        d = dist(idxs[0], idxs[1])
        if abs(d - 2.28) > 0.02:
            fails.append(f'{name}: Br2 bond {d:.3f}')
    elif an == 'H2O':
        idxsO = species.index('O')
        wh = [i for i, s in enumerate(species) if s == 'H' and abs(dist(idxsO, i) - 0.97) < 0.02]
        if len(wh) != 2:
            fails.append(f'{name}: water H count {len(wh)} != 2')
        else:
            v1 = [cart[wh[0]][k]-cart[idxsO][k] for k in range(3)]
            v2 = [cart[wh[1]][k]-cart[idxsO][k] for k in range(3)]
            dot = sum(a*b for a, b in zip(v1, v2))
            ang = math.degrees(math.acos(max(-1, min(1, dot / (0.97*0.97)))))
            if abs(ang - 104.5) > 1.0:
                fails.append(f'{name}: HOH angle {ang:.1f}')

    # global cation-anion min contact (anion = all Br + water O/H for H2O cases)
    anion_set = {i for i, s in enumerate(species) if s == 'Br'}
    if an == 'H2O':
        iO = species.index('O')
        anion_set |= {iO} | {i for i, s in enumerate(species) if s == 'H' and dist(iO, i) < 1.1}
    if anion_set:
        mn = min(dist(i, j) for i in anion_set for j in range(len(species)) if j not in anion_set)
        if mn < 1.5:
            fails.append(f'{name}: cation-anion min contact {mn:.3f} < 1.5')

# singles
for d in sorted(os.listdir(os.path.join(ROOT, '01_single'))):
    p = os.path.join(ROOT, '01_single', d, 'POSCAR')
    sp, c, box = read_poscar(p)
    if len(sp) != EXPECT_SINGLE[d]:
        fails.append(f'01_single/{d}: natoms {len(sp)} != {EXPECT_SINGLE[d]}')
    check(d, sp, c, box)
    if d == 'Br3':
        dd = math.dist(c[0], c[1])
        if abs(dd - 2.55) > 0.02: fails.append(f'Br3 spacing {dd:.3f}')
    if d == 'Br5':
        # order: M(vertex), X1, Y1, X2, Y2
        for lbl, dd, exp in (('M-X1', math.dist(c[0], c[1]), 3.10),
                             ('X1-Y1', math.dist(c[1], c[2]), 2.40),
                             ('M-X2', math.dist(c[0], c[3]), 3.10),
                             ('X2-Y2', math.dist(c[3], c[4]), 2.40)):
            if abs(dd - exp) > 0.02: fails.append(f'Br5 single {lbl} {dd:.3f}')
    if d == 'H2O':
        iO = sp.index('O')
        hs = [i for i, s in enumerate(sp) if s == 'H']
        ang = math.degrees(math.acos(max(-1, min(1, sum((c[hs[0]][k]-c[iO][k])*(c[hs[1]][k]-c[iO][k]) for k in range(3)) / (0.97*0.97)))))
        if abs(ang - 104.5) > 1.0: fails.append(f'H2O angle {ang:.1f}')

# complexes
for d in sorted(os.listdir(os.path.join(ROOT, '02_complex'))):
    p = os.path.join(ROOT, '02_complex', d, 'POSCAR')
    sp, c, box = read_poscar(p)
    check(d, sp, c, box)

with open(os.path.join(ROOT, '00_structures', 'cases_manifest.csv'), encoding='utf-8-sig') as f:
    rows = list(csv.DictReader(f))
print(f'checked {checked} POSCARs; manifest rows {len(rows)}')
if fails:
    print(f'FAIL ({len(fails)}):')
    for x in fails[:20]:
        print(' ', x)
    sys.exit(1)
print('ALL GEOMETRY CHECKS PASSED')
