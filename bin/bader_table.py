#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Consolidate Bader charges for the LAECHG=True subset.

    python3 bin/bader_table.py

Reads  02_calc/sp_vac/<case>/bader.csv    (per-atom net charge)
       04_results/eb.csv                 (binding energies)
Writes 04_results/bader/bader_summary.{csv,md}

Only the three cases that were run with LAECHG=.TRUE. have AECCAR files, so this
is deliberately a small table: C10TAB_Br3_c1 / C12TAB_Br3_c2 / C14TAB_Br3_c1.
"""

import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

KIT = common.KIT
CASES = ["C10TAB_Br3_c1", "C12TAB_Br3_c2", "C14TAB_Br3_c1"]
OUT = os.path.join(KIT, "04_results", "bader")
LABEL = {"C10TAB_Br3_c1": "C10TAB", "C12TAB_Br3_c2": "C12TAB", "C14TAB_Br3_c1": "C14TAB"}


def main():
    eb = {}
    p = os.path.join(KIT, "04_results", "eb.csv")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("eb_vac"):
                    eb[r["system"]] = float(r["eb_vac"])

    rows = []
    for c in CASES:
        p = os.path.join(KIT, "02_calc", "sp_vac", c, "bader.csv")
        if not os.path.exists(p):
            print("missing", p)
            continue
        agg = defaultdict(lambda: [0, 0.0])
        with open(p, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                e = r["element"]
                q = float(r["net_charge"])
                agg[e][0] += 1
                agg[e][1] += q
        anion = agg.get("Br", [0, 0.0])
        row = {
            "case": c, "label": LABEL.get(c, c),
            "eb_vac": eb.get(c),
            "n_Br": anion[0], "q_Br_total": anion[1],
            "q_Br_mean": (anion[1] / anion[0]) if anion[0] else None,
            "q_N": agg.get("N", [0, 0.0])[1] or None,
            "q_C_total": agg.get("C", [0, 0.0])[1],
            "q_H_total": agg.get("H", [0, 0.0])[1],
        }
        rows.append(row)

    if not rows:
        return

    os.makedirs(OUT, exist_ok=True)
    csv_path = os.path.join(OUT, "bader_summary.csv")
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.4f" % v if isinstance(v, float) else v) for k, v in r.items()})

    lines = ["# Bader 电荷（LAECHG=.TRUE. 的三个案例）\n",
             "净电荷 = ZVAL − Bader 电荷。单位 e。\n",
             "| 体系 | Eb (eV) | Br₃⁻ 总电荷 | 单个 Br 平均 | N |",
             "|---|---|---|---|---|"]
    for r in rows:
        lines.append("| %s | %.3f | %+.3f | %+.3f | %+.3f |" % (
            r["label"], r["eb_vac"], r["q_Br_total"], r["q_Br_mean"], r["q_N"]))
    lines.append("")
    lines.append("对照组：Br₃⁻ 的总电荷越负 = 从阳离子拿走的电子越少 = 结合越弱。")

    md_path = os.path.join(OUT, "bader_summary.md")
    with open(md_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

    print("written", csv_path)
    print("written", md_path)
    print()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
