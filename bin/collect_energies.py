#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Collect OUTCAR energies, compute binding / exchange energies.

    python3 bin/collect_energies.py
writes 04_results/{energies_raw,eb,exchange}.csv and eb.md
"""

import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

TASKS = ["opt", "opt_sol", "sp_vac", "sp_sol"]
CATIONS = ["C10TAB", "C12TAB", "C14TAB", "MEPBr", "TPABr", "EMImBr"]
ANIONS = ["Br", "Br3", "Br5", "Br2", "H2O"]


def gather():
    E = {}
    for task in TASKS:
        td = os.path.join(common.KIT, "02_calc", task)
        if not os.path.isdir(td):
            continue
        for name in os.listdir(td):
            d = os.path.join(td, name)
            e = common.outcar_energy(os.path.join(d, "OUTCAR"))
            # only converged relaxations (or finished single points) count:
            # a running job's OUTCAR already holds energies that are still drifting
            if not e:
                continue
            if task.startswith("sp"):
                if e["finished"]:
                    E.setdefault(name, {})[task] = e
            elif e["converged"]:
                E.setdefault(name, {})[task] = e
    return E


def level(E, name, solv):
    if name not in E:
        return None, None
    for k in (["sp_sol", "opt_sol"] if solv else ["sp_vac", "opt"]):
        if k in E[name]:
            return E[name][k]["toten"], k
    return None, None


def fnum(x, nd=3):
    return "" if x is None else ("%.*f" % (nd, x))


def main():
    out = os.path.join(common.KIT, "04_results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(common.KIT, "00_structures", "manifest.json"),
              encoding="utf-8") as f:
        man = json.load(f)
    E = gather()

    with open(os.path.join(out, "energies_raw.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["system", "task", "toten_eV", "sigma0_eV", "finished"])
        for name in sorted(E):
            for task in TASKS:
                if task in E[name]:
                    e = E[name][task]
                    w.writerow([name, task, "%.6f" % e["toten"],
                                fnum(e["sigma0"], 6), int(e["finished"])])

    rows = []
    for rec in man["systems"]:
        if rec["kind"] != "complex":
            continue
        name, cat, an = rec["name"], rec["cation"], rec["anion"]
        row = {"system": name, "cation": cat, "anion": an,
               "conf": name[-2:], "eb_vac": None, "eb_solv": None,
               "level_vac": "", "level_solv": ""}
        for solv, key, lk in ((False, "eb_vac", "level_vac"),
                              (True, "eb_solv", "level_solv")):
            ec, lc = level(E, name, solv)
            ea, _ = level(E, cat, solv)
            eb, _ = level(E, an, solv)
            if None not in (ec, ea, eb):
                row[key] = ec - ea - eb
                row[lk] = lc
        rows.append(row)

    best = {}
    for r in rows:
        key = (r["cation"], r["anion"])
        val = r["eb_solv"] if r["eb_solv"] is not None else r["eb_vac"]
        if val is None:
            continue
        if key not in best or val < best[key][0]:
            best[key] = (val, r["system"])
    for r in rows:
        r["best"] = int(best.get((r["cation"], r["anion"]), (None, ""))[1] == r["system"])

    with open(os.path.join(out, "eb.csv"), "w", newline="", encoding="utf-8") as f:
        cols = ["system", "cation", "anion", "conf", "eb_vac", "eb_solv",
                "level_vac", "level_solv", "best"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            rr = dict(r)
            for k in ("eb_vac", "eb_solv"):
                rr[k] = fnum(rr[k], 4)
            w.writerow(rr)

    ex = []
    for cat in CATIONS:
        rec = {"cation": cat}
        for solv, suf in ((False, "vac"), (True, "solv")):
            def get(n):
                return level(E, n, solv)[0]
            e_cb, e_cb3 = get(cat + "_Br_c1"), get(cat + "_Br3_c1")
            e_b, e_b3 = get("Br"), get("Br3")
            if None not in (e_cb, e_cb3, e_b, e_b3):
                # exchange reaction: complex_Br + Br3-  ->  complex_Br3 + Br-
                # dE_ex > 0 => the Br- complex is the lower-energy one
                rec["dE_ex_" + suf] = (e_cb3 + e_b) - (e_cb + e_b3)
            e_h2o = get(cat + "_H2O_c1")
            if None not in (e_cb3, e_h2o, e_b3):
                eb_b3 = e_cb3 - get(cat) - e_b3
                eb_h2o = e_h2o - get(cat) - get("H2O")
                rec["dE_H2O_" + suf] = eb_b3 - eb_h2o
        ex.append(rec)

    with open(os.path.join(out, "exchange.csv"), "w", newline="",
              encoding="utf-8") as f:
        cols = ["cation", "dE_ex_vac", "dE_ex_solv", "dE_H2O_vac", "dE_H2O_solv"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in ex:
            w.writerow({k: (fnum(r.get(k), 3) if k != "cation" else r[k]) for k in cols})

    lines = ["# Binding energies (eV)", "",
             "Best conformer per cation/anion; vac = vacuum, solv = VASPsol.", "",
             "| cation | " + " | ".join(ANIONS) + " |",
             "|---|" + "---|" * len(ANIONS)]
    table = {}
    for r in rows:
        if r["best"]:
            table[(r["cation"], r["anion"])] = r
    for cat in CATIONS:
        cells = []
        for an in ANIONS:
            r = table.get((cat, an))
            if not r:
                cells.append("")
            else:
                v = r["eb_solv"] if r["eb_solv"] is not None else r["eb_vac"]
                cells.append(fnum(v, 2))
        lines.append("| %s | %s |" % (cat, " | ".join(cells)))
    lines += ["", "# Exchange energies (eV)", "",
              "| cation | dE_ex vac | dE_ex solv | dE_H2O vac | dE_H2O solv |",
              "|---|---|---|---|---|"]
    for r in ex:
        lines.append("| %s | %s | %s | %s | %s |" % (
            r["cation"], fnum(r.get("dE_ex_vac"), 2), fnum(r.get("dE_ex_solv"), 2),
            fnum(r.get("dE_H2O_vac"), 2), fnum(r.get("dE_H2O_solv"), 2)))
    lines.append("")
    lines.append("dE_ex < 0: Br3- binds preferentially over Br-; "
                 "dE_H2O < 0: Br3- preferred over water.")
    with open(os.path.join(out, "eb.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))

    print("wrote 04_results/{energies_raw,eb,exchange}.csv and eb.md")
    print("\n".join(lines[:12]))


if __name__ == "__main__":
    main()
