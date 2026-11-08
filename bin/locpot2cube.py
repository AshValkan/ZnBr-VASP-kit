#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert LOCPOT to a Gaussian cube file (for VESTA / Multiwfn).

    python3 bin/locpot2cube.py 02_calc/sp_vac/C12TAB [stride]

The cube is subsampled (default stride 3, ~0.3-0.5 A) to keep the file small.
Potential is written in eV.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridio

BOHR = 1.0 / 0.529177210903
ZVAL = {"H": 1, "C": 4, "N": 5, "O": 6, "Br": 7}


def main():
    if len(sys.argv) < 2:
        print("usage: locpot2cube.py <calc_dir> [stride]")
        return
    d = sys.argv[1]
    stride = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    dims, pot, meta = gridio.read_grid(os.path.join(d, "LOCPOT"), stride)
    lat = np.array(meta["lattice"])
    nx, ny, nz = dims
    sx, sy, sz = nx // stride, ny // stride, nz // stride

    elems = []
    for s, c in zip(meta["species"], meta["counts"]):
        elems += [s] * c
    coords = np.array(meta["coords"])

    out = os.path.join(d, "LOCPOT.cube")
    with open(out, "w", encoding="utf-8") as f:
        f.write("LOCPOT %s  potential(eV) stride=%d\n" % (os.path.basename(d.rstrip("/")), stride))
        f.write("outside\n")
        f.write("%5d %12.6f %12.6f %12.6f\n" % (len(elems), 0.0, 0.0, 0.0))
        for k in range(3):
            vec = lat[k] / sx if k == 0 else (lat[k] / sy if k == 1 else lat[k] / sz)
            f.write("%5d %12.6f %12.6f %12.6f\n" %
                    ((sx, sy, sz)[k], vec[0] * BOHR, vec[1] * BOHR, vec[2] * BOHR))
        for e, r in zip(elems, coords):
            f.write("%5d %12.6f %12.6f %12.6f %12.6f\n" %
                    (ZVAL.get(e, 1), float(r[0]) * BOHR, float(r[1]) * BOHR,
                     float(r[2]) * BOHR, 0.0))
        for iz in range(sz):
            for iy in range(sy):
                row = pot[iz, iy, :]
                for ix in range(0, sx, 6):
                    f.write("".join("%13.5E" % v for v in row[ix:ix + 6]) + "\n")
    print("written", out, "grid", pot.shape)


if __name__ == "__main__":
    main()
