#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Consolidate per-case esp_stats.csv into one table.

    python3 bin/esp_table.py

Reads  02_calc/sp_vac/<case>/esp_stats.csv
       04_results/eb.csv                     (binding energies)
Writes 04_results/esp/esp_summary.csv
       04_results/esp/esp_summary.md         (matrices, Chinese captions)

Only the fields that survive scrutiny are carried forward: Vmax / Vmin / Vmean /
sigma.  frac_pos_1eV is deliberately dropped -- it is a one-sided threshold that
reads 0 for every cation and therefore carries no information.
"""

import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

KIT = common.KIT
SP = os.path.join(KIT, "02_calc", "sp_vac")
OUT = os.path.join(KIT, "04_results", "esp")

ANION_ORDER = ["Br", "Br3", "Br5", "Br2", "H2O"]
ANION_LABEL = {"Br": "Br⁻", "Br3": "Br₃⁻", "Br5": "Br₅⁻", "Br2": "Br₂", "H2O": "H₂O"}
CATION_ORDER = ["C10TAB", "C12TAB", "C14TAB"]


def read_esp(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.reader(f):
            if len(row) == 2 and row[0] != "quantity":
                try:
                    out[row[0]] = float(row[1])
                except ValueError:
                    out[row[0]] = row[1]
    return out


def split_name(name):
    """C12TAB_Br3_c2 -> (C12TAB, Br3, c2); C12TAB -> (C12TAB, '', '')"""
    m = re.match(r"^(C\d+TAB)(?:_(Br3|Br5|Br2|Br|H2O))?(?:_(c\d+))?$", name)
    if not m:
        return name, "", ""
    return m.group(1), m.group(2) or "", m.group(3) or ""


def main():
    eb = {}
    eb_path = os.path.join(KIT, "04_results", "eb.csv")
    if os.path.exists(eb_path):
        with open(eb_path, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("eb_vac"):
                    eb[r["system"]] = float(r["eb_vac"])

    rows = []
    for name in sorted(os.listdir(SP)):
        p = os.path.join(SP, name, "esp_stats.csv")
        if not os.path.exists(p):
            continue
        d = read_esp(p)
        cat, an, conf = split_name(name)
        rows.append({
            "system": name, "cation": cat, "anion": an, "conf": conf,
            "points": int(d.get("points", 0)),
            "Vmax": d.get("Vmax"), "Vmin": d.get("Vmin"),
            "Vmean": d.get("Vmean"), "sigma": d.get("sigma"),
            "eb_vac": eb.get(name),
        })

    if not rows:
        print("no esp_stats.csv found")
        return

    os.makedirs(OUT, exist_ok=True)
    csv_path = os.path.join(OUT, "esp_summary.csv")
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.4f" % v if isinstance(v, float) else v)
                        for k, v in r.items()})

    # 矩阵：sigma
    mat = {}
    for r in rows:
        mat[(r["cation"], r["anion"])] = r["sigma"]

    lines = []
    lines.append("# 静电势统计（ESP）\n")
    lines.append("等值面 ρ = 1.0×10⁻³ e/Å³，已扣除真空能级。单位 eV。\n")
    lines.append("## 表面电势标准差 σ\n")
    lines.append("| 体系 | " + " | ".join(ANION_LABEL[a] for a in ANION_ORDER) + " | 单体 |")
    lines.append("|---" * 7 + "|")
    for c in CATION_ORDER:
        cells = []
        for a in ANION_ORDER:
            v = mat.get((c, a))
            cells.append("%.3f" % v if v is not None else "—")
        mono = mat.get((c, ""))
        cells.append("%.3f" % mono if mono is not None else "—")
        lines.append("| %s | %s |" % (c, " | ".join(cells)))
    lines.append("")
    lines.append("## 逐体系明细\n")
    lines.append("| 体系 | 采样点 | V_max | V_min | V_mean | σ | Eb (eV) |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda x: x["system"]):
        lines.append("| %s | %d | %.3f | %.3f | %.3f | %.3f | %s |" % (
            r["system"], r["points"], r["Vmax"], r["Vmin"], r["Vmean"], r["sigma"],
            ("%.3f" % r["eb_vac"]) if r["eb_vac"] is not None else "—"))

    md_path = os.path.join(OUT, "esp_summary.md")
    with open(md_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

    print("cases: %d" % len(rows))
    print("written", csv_path)
    print("written", md_path)
    print()
    print("\n".join(lines[5:11]))


if __name__ == "__main__":
    main()
