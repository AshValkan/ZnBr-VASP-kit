#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Export verified case directories (01_single/, 02_complex/) into kit layout:
00_structures/<name>.vasp + 00_structures/manifest.json.

Run after generate_structures.py + verify_structures.py:
    python3 bin/export_kit_manifest.py <calc_root>
"""
import csv, json, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'

CATION_ANIONS = {   # complex matrix (must mirror generate_structures.py)
    'C10TAB': ['Br', 'Br3', 'Br5', 'Br2', 'H2O'],
    'C12TAB': ['Br', 'Br3', 'Br5', 'Br2', 'H2O'],
    'C14TAB': ['Br', 'Br3', 'Br5', 'Br2', 'H2O'],
    'MEPBr':  ['Br3'],
    'TPABr':  ['Br3'],
    'EMImBr': ['Br3', 'H2O'],
}
CHARGED_ANIONS = ('Br', 'Br3', 'Br5')     # carry -1
KIND_ANION = ('Br', 'Br3', 'Br5')         # kind=anion; Br2/H2O kind=molecule


def parse_case(d):
    pos = os.path.join(d, 'POSCAR')
    lines = open(pos, encoding='utf-8').read().splitlines()
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    charge = int(open(os.path.join(d, 'charge.txt')).read().strip())
    return species, counts, charge, lines


def main():
    out_dir = os.path.join(ROOT, '00_structures')
    os.makedirs(out_dir, exist_ok=True)
    systems = []
    rows = []
    man_csv = os.path.join(out_dir, 'cases_manifest.csv')
    if os.path.exists(man_csv):
        rows = list(csv.DictReader(open(man_csv, encoding='utf-8-sig')))
    notes = {r['case']: r['note'] for r in rows}

    for sub, kind_base in (('01_single', 'single'), ('02_complex', 'complex')):
        base = os.path.join(ROOT, sub)
        if not os.path.isdir(base):
            continue
        for name in sorted(os.listdir(base)):
            d = os.path.join(base, name)
            if not os.path.isdir(d):
                continue
            species, counts, charge, lines = parse_case(d)
            natoms = sum(counts)
            # classify
            if sub == '01_single':
                if name in CATION_ANIONS:
                    kind, cat, an = 'cation', name, None
                elif name in KIND_ANION:
                    kind, cat, an = 'anion', None, name
                else:
                    kind, cat, an = 'molecule', None, name
            else:
                cat, rest = name.split('_', 1)
                an, conf = rest.rsplit('_c', 1)
                kind = 'complex'
            # key atoms (0-based indices into grouped POSCAR order)
            idx = []
            for s, c in zip(species, counts):
                idx += [s] * c
            key = {
                'N': [i for i, e in enumerate(idx) if e == 'N'],
                'Br': [i for i, e in enumerate(idx) if e == 'Br'],
                'O': [i for i, e in enumerate(idx) if e == 'O'],
            }
            # write <name>.vasp (same content as POSCAR)
            with open(os.path.join(out_dir, name + '.vasp'), 'w', newline='\n') as f:
                f.write('\n'.join(lines) + '\n')
            systems.append({
                'name': name, 'kind': kind, 'cation': cat, 'anion': an,
                'charge': charge, 'natoms': natoms, 'box': 44.0,
                'species': species, 'counts': counts,
                'key': {'N': key['N'] if kind != 'anion' else [],
                        'Br': key['Br'], 'O': key['O']},
                'file': f'00_structures/{name}.vasp',
                'note': notes.get(f'{sub}/{name}', ''),
            })
    man = {'box': 44.0, 'systems': systems}
    with open(os.path.join(out_dir, 'manifest.json'), 'w', encoding='utf-8', newline='\n') as f:
        json.dump(man, f, ensure_ascii=False, indent=1)
    kinds = {}
    for s in systems:
        kinds[s['kind']] = kinds.get(s['kind'], 0) + 1
    print(f'exported {len(systems)} systems -> 00_structures/manifest.json  {kinds}')


if __name__ == '__main__':
    main()
