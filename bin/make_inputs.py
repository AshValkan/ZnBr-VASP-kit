#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create VASP run directories from 00_structures + manifest.json.

    python3 bin/make_inputs.py --tasks opt,opt_sol
    python3 bin/make_inputs.py --tasks sp_vac,freq --systems C12TAB_Br3_c1

Adapted from ZnBr_VASP_kit: NELECT/ISPIN written automatically from the
manifest charge (charged anions/cations, odd-electron -> ISPIN=2).
"""

import argparse
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
import vasp_inputs


def potcar(d, cfg, species, warn=True):
    src = cfg.get("potcar_dir", "")
    out = os.path.join(d, "POTCAR")
    if not src or not os.path.isdir(src):
        if warn:
            print("  ! potcar_dir not set/found -> POTCAR skipped in", d)
        return False
    chunks = []
    for e in species:
        name = cfg["potcar_map"].get(e, e)
        p = os.path.join(src, name, "POTCAR")
        if not os.path.exists(p):
            print("  ! missing POTCAR for", e, "at", p)
            return False
        with open(p, encoding="utf-8", errors="replace") as f:
            chunks.append(f.read())
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("".join(chunks))
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", default="opt,opt_sol")
    ap.add_argument("--systems", default="all")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-potcar", action="store_true")
    args = ap.parse_args()

    cfg = common.load_config()
    struct_dir = os.path.join(common.KIT, "00_structures")
    with open(os.path.join(struct_dir, "manifest.json"), encoding="utf-8") as f:
        man = json.load(f)
    recs = {s["name"]: s for s in man["systems"]}
    systems = list(recs)
    if args.systems != "all":
        systems = [s for s in systems if s in args.systems.split(",")]

    tasks = [t.strip() for t in args.tasks.split(",") if t.strip()]
    total = 0
    for task in tasks:
        solv = task.endswith("_sol")
        for name in systems:
            rec = recs[name]
            d = os.path.join(common.KIT, "02_calc", task, name)
            if os.path.exists(os.path.join(d, "OUTCAR")) and not args.force:
                continue
            os.makedirs(d, exist_ok=True)
            src = os.path.join(struct_dir, name + ".vasp")
            dst = os.path.join(d, "POSCAR")
            if not os.path.exists(dst) or args.force:
                shutil.copyfile(src, dst)
            vasp_inputs.write_inputs(d, cfg, name, task, solv,
                                     charge=int(rec.get("charge", 0)),
                                     poscar_path=src)
            if not args.no_potcar:
                potcar(d, cfg, rec["species"])
            total += 1
            print("created", os.path.relpath(d, common.KIT))
    print("\n%d directories ready under 02_calc/" % total)


if __name__ == "__main__":
    main()
