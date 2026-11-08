#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Optional figures from 04_results (needs matplotlib)."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    print("matplotlib not available - CSV results are complete, skipping figures")
    sys.exit(0)

import csv


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    res = os.path.join(common.KIT, "04_results")
    fig = os.path.join(res, "figures")
    os.makedirs(fig, exist_ok=True)

    rows = [r for r in read_csv(os.path.join(res, "eb.csv")) if r.get("best") == "1"]
    if rows:
        cats = sorted({r["cation"] for r in rows}, key=lambda c: ["C10TAB", "C12TAB", "C14TAB", "MEP", "TPA", "EMIm"].index(c) if c in ["C10TAB", "C12TAB", "C14TAB", "MEP", "TPA", "EMIm"] else 99)
        ans = ["Br", "Br3", "Br5", "Br2", "H2O"]
        w = 0.15
        fig1, ax = plt.subplots(figsize=(8, 4.2))
        for k, a in enumerate(ans):
            vals = []
            for c in cats:
                v = [r for r in rows if r["cation"] == c and r["anion"] == a]
                vals.append(float(v[0]["eb_solv"] or v[0]["eb_vac"]) if v and (v[0]["eb_solv"] or v[0]["eb_vac"]) else 0.0)
            ax.bar([i + k * w for i in range(len(cats))], vals, width=w, label=a)
        ax.set_xticks([i + 2 * w for i in range(len(cats))])
        ax.set_xticklabels(cats)
        ax.set_ylabel("Eb (eV)")
        ax.axhline(0, color="k", lw=0.6)
        ax.legend(ncol=5, fontsize=8)
        fig1.tight_layout()
        fig1.savefig(os.path.join(fig, "eb.png"), dpi=200)

    ts = read_csv(os.path.join(common.KIT, "04_results", "aimd_timeseries.csv"))
    if not ts and os.path.isdir(os.path.join(common.KIT, "02_calc", "aimd")):
        for name in os.listdir(os.path.join(common.KIT, "02_calc", "aimd")):
            p = os.path.join(common.KIT, "02_calc", "aimd", name, "aimd_timeseries.csv")
            ts = read_csv(p)
            if ts:
                break
    if ts:
        fig2, ax = plt.subplots(figsize=(6, 3.6))
        ax.plot([float(r["t_ps"]) for r in ts], [float(r["disp_A"]) for r in ts])
        ax.set_xlabel("t (ps)")
        ax.set_ylabel("|dr| of Br3- (A)")
        fig2.tight_layout()
        fig2.savefig(os.path.join(fig, "aimd_displacement.png"), dpi=200)
    print("figures in", fig)


if __name__ == "__main__":
    main()
