#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarise Bader charges from ACF.dat.

    python3 bin/bader_summary.py 02_calc/sp_vac/C12TAB

Prints per-element net charges and writes bader.csv next to ACF.dat.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
import json

ZVAL = {"H": 1, "C": 4, "N": 5, "O": 6, "Br": 7}


def main():
    if len(sys.argv) < 2:
        print("usage: bader_summary.py <calc_dir>")
        return
    d = sys.argv[1]
    acf = os.path.join(d, "ACF.dat")
    if not os.path.exists(acf):
        print("ACF.dat not found in", d)
        return
    pos = common.read_poscar(os.path.join(d, "POSCAR"))
    elems = pos["elems"]
    rows = []
    for line in open(acf, encoding="utf-8", errors="replace"):
        p = line.split()
        if len(p) >= 5 and p[0].isdigit():
            rows.append((int(p[0]), float(p[4])))
    n = min(len(rows), len(elems))
    per = {}
    with open(os.path.join(d, "bader.csv"), "w", encoding="utf-8") as f:
        f.write("index,element,bader_charge,net_charge\n")
        for i in range(n):
            idx, q = rows[i]
            e = elems[i]
            net = ZVAL.get(e, 0) - q
            per.setdefault(e, []).append(net)
            f.write("%d,%s,%.4f,%.4f\n" % (idx, e, q, net))
    print("element   n   mean net charge (e)")
    for e in ("N", "Br", "O", "C"):
        if e in per:
            v = per[e]
            print("%-8s %3d   %+.4f" % (e, len(v), sum(v) / len(v)))
    print("written", os.path.join(d, "bader.csv"))


if __name__ == "__main__":
    main()
