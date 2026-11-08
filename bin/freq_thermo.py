#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Harmonic thermochemistry from a VASP frequency run.

    python3 bin/freq_thermo.py [--all]

Reads every 02_calc/freq/<system>, parses the frequencies from OUTCAR, and
writes 04_results/thermo.csv (E, ZPE, S_vib, S_trans, S_rot, G) plus dG for
complexes when the monomers are available.
"""

import argparse
import csv
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

HC = 1.239841984e-4      # eV per cm^-1
KB = 8.617333262e-5      # eV/K
KB_SI = 1.380649e-23
H_SI = 6.62607015e-34
AMU = 1.66053906660e-27
P0 = 101325.0
T = 298.15


def freqs(path):
    out = []
    for line in open(path, encoding="utf-8", errors="replace"):
        if "cm-1" in line and ("f  =" in line or "f/i=" in line):
            m = re.search(r"([-+]?\d+\.\d+)\s+cm-1", line)
            if m:
                out.append(float(m.group(1)))
    return out


def thermo_one(d, cfg):
    e = common.outcar_energy(os.path.join(d, "OUTCAR"))
    if not e:
        return None
    fr = [f for f in freqs(os.path.join(d, "OUTCAR")) if f > 0]
    n_imag = len([f for f in freqs(os.path.join(d, "OUTCAR")) if f < 0])
    if not fr:
        return None
    pos = common.read_poscar(os.path.join(d, "CONTCAR")
                             if os.path.exists(os.path.join(d, "CONTCAR"))
                             else os.path.join(d, "POSCAR"))
    m = np.array([cfg["masses"].get(x, 12.0) for x in pos["elems"]])
    xyz = np.array(pos["coords"])
    com = (m[:, None] * xyz).sum(0) / m.sum()
    r = xyz - com
    I = np.zeros((3, 3))
    for mi, ri in zip(m, r):
        I += mi * (np.dot(ri, ri) * np.eye(3) - np.outer(ri, ri))
    I *= AMU * 1e-20
    w = np.linalg.eigvalsh(I)
    w = np.clip(w, 1e-45, None)

    zpe = 0.5 * HC * sum(fr)
    x = np.array(fr) * HC / (KB * T)
    s_vib = KB * np.sum(x / (np.exp(x) - 1.0) - np.log(1.0 - np.exp(-x)))

    m_tot = m.sum() * AMU
    s_trans = KB_SI * (np.log((2 * np.pi * m_tot * KB_SI * T / H_SI ** 2) ** 1.5
                              * KB_SI * T / P0) + 2.5) / 1.602176634e-19
    s_rot = KB_SI * (np.log(np.sqrt(np.pi * w[0] * w[1] * w[2])
                            * (8 * np.pi ** 2 * KB_SI * T / H_SI ** 2) ** 1.5) + 1.5) \
        / 1.602176634e-19
    s_tot = s_vib + s_trans + s_rot
    g = e["toten"] + zpe - T * s_tot
    return {"system": os.path.basename(os.path.normpath(d)), "E": e["toten"],
            "ZPE": zpe, "S_vib": s_vib, "S_trans": s_trans, "S_rot": s_rot,
            "S_tot": s_tot, "G": g, "n_imag": n_imag, "nmodes": len(fr)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    cfg = common.load_config()
    root = os.path.join(common.KIT, "02_calc", "freq")
    if not os.path.isdir(root):
        print("no 02_calc/freq")
        return
    rows = []
    for name in sorted(os.listdir(root)):
        r = thermo_one(os.path.join(root, name), cfg)
        if r:
            rows.append(r)
            print("%-22s G=%10.4f  ZPE=%.3f  S=%.4f eV/K  imag=%d"
                  % (r["system"], r["G"], r["ZPE"], r["S_tot"], r["n_imag"]))
    out = os.path.join(common.KIT, "04_results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "thermo.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["system"])
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.6f" % v if isinstance(v, float) else v)
                        for k, v in r.items()})
    print("written 04_results/thermo.csv")

    G = {r["system"]: r["G"] for r in rows}
    import json
    man = json.load(open(os.path.join(common.KIT, "00_structures", "manifest.json"),
                         encoding="utf-8"))
    dg = []
    for rec in man["systems"]:
        if rec["kind"] != "complex":
            continue
        n, c, a = rec["name"], rec["cation"], rec["anion"]
        if n in G and c in G and a in G:
            dg.append((n, G[n] - G[c] - G[a]))
    if dg:
        with open(os.path.join(out, "dg.csv"), "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["system", "dG_eV"])
            for n, v in dg:
                w.writerow([n, "%.4f" % v])
        print("written 04_results/dg.csv (%d reactions)" % len(dg))


if __name__ == "__main__":
    main()
