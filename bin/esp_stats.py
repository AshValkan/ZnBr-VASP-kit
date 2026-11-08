#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Electrostatic-potential statistics on a density isosurface.

    python3 bin/esp_stats.py 02_calc/sp_vac/C12TAB

Reads LOCPOT + CHGCAR, samples the rho = 0.001 e/A^3 surface, subtracts the
vacuum level (mean potential where rho < 1e-6) and reports V_S statistics.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridio

RHO0 = 1.0e-3          # e/A^3 isosurface
TOL = 0.25             # relative window around RHO0


def main():
    if len(sys.argv) < 2:
        print("usage: esp_stats.py <calc_dir> [stride]")
        return
    d = sys.argv[1]
    stride = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    dims, pot, meta = gridio.read_grid(os.path.join(d, "LOCPOT"), stride)
    _, chg, _ = gridio.read_grid(os.path.join(d, "CHGCAR"), stride)
    vol = gridio.volume(meta["lattice"])
    rho = chg / vol

    vac = pot[rho < 1e-6]
    ref = float(vac.mean()) if vac.size else 0.0
    v = pot - ref

    sel = np.abs(rho - RHO0) < TOL * RHO0
    if sel.sum() < 100:
        sel = np.abs(rho - RHO0) < 3 * TOL * RHO0
    vs = v[sel]
    out = {
        "dir": os.path.basename(os.path.normpath(d)),
        "points": int(sel.sum()),
        "ref_eV": ref,
        "Vmax": float(vs.max()),
        "Vmin": float(vs.min()),
        "Vmean": float(vs.mean()),
        "sigma": float(vs.std()),
        "frac_pos_1eV": float((vs > 1.0).mean()),
    }
    out_path = os.path.join(d, "esp_stats.csv")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("quantity,value\n")
        for k, val in out.items():
            f.write("%s,%s\n" % (k, val if isinstance(val, str) else "%.6f" % val))
    print("V_S max = %.3f eV   V_S min = %.3f eV   sigma = %.3f eV   "
          "area>1eV = %.1f%%" % (out["Vmax"], out["Vmin"], out["sigma"],
                                 100 * out["frac_pos_1eV"]))
    print("written", out_path)


if __name__ == "__main__":
    main()
