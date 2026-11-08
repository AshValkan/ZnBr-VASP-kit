#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Environment doctor: run before first use on a new machine.

    python3 bin/doctor.py

Checks config paths, POTCAR availability, binaries, python deps and the
structure manifest.  Exit code 0 = ready, 1 = something needs fixing.
"""
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

FAILS, WARNS, OKS = [], [], []


def check(name, ok, hint="", warn_only=False):
    if ok:
        OKS.append(name)
        print("  [PASS] %s" % name)
    else:
        (WARNS if warn_only else FAILS).append(name)
        tag = "WARN" if warn_only else "FAIL"
        print("  [%s] %s%s" % (tag, name, ("  <- " + hint) if hint else ""))


def main():
    cfg = common.load_config()
    print("== ZnBr VASP kit doctor ==")

    # 1. manifest + structures
    man = os.path.join(common.KIT, "00_structures", "manifest.json")
    check("manifest.json", os.path.exists(man))
    n_vasp = len(glob_vasp()) if os.path.isdir(os.path.join(common.KIT, "00_structures")) else 0
    n_sys = 0
    if os.path.exists(man):
        systems = json.load(open(man, encoding="utf-8"))["systems"]
        n_sys = len(systems)
        missing = [s["name"] for s in systems
                   if not os.path.exists(os.path.join(common.KIT, "00_structures", s["name"] + ".vasp"))]
        check("structures %d/%d .vasp present" % (n_vasp, n_sys),
              n_vasp >= n_sys and not missing,
              "missing: %s" % missing[:5] if missing else "")

    # 2. binaries
    vasp = os.path.expanduser(cfg.get("vasp_std", ""))
    check("vasp binary: %s" % vasp, os.path.exists(vasp),
          "edit config.json -> vasp_std")
    mpi = cfg.get("mpi_cmd", "").split()[0] if cfg.get("mpi_cmd") else ""
    check("mpi launcher: %s" % mpi, bool(mpi) and (os.path.exists(mpi) or shutil.which(mpi)),
          "edit config.json -> mpi_cmd")

    # 3. POTCAR
    pp = os.path.expanduser(cfg.get("potcar_dir", ""))
    if os.path.isdir(pp):
        need = set()
        if os.path.exists(man):
            for s in systems:
                need.update(s["species"])
        missing_pp = [e for e in sorted(need)
                      if not os.path.exists(os.path.join(pp, cfg.get("potcar_map", {}).get(e, e), "POTCAR"))]
        check("POTCAR for %s (in %s)" % (", ".join(sorted(need)), pp), not missing_pp,
              "missing elements: %s" % missing_pp)
    else:
        check("potcar_dir: %s" % pp, False, "edit config.json -> potcar_dir")

    # 4. scheduler sanity
    blocks = cfg.get("core_blocks") or []
    nproc = os.cpu_count() or 1
    if cfg.get("scheduler") == "direct":
        ok_blocks = all(max(_rng(b)) < nproc for b in blocks) if blocks else False
        check("core_blocks within %d cpus: %s" % (nproc, blocks), ok_blocks or not blocks,
              "core_blocks exceed cpu count")
        check("max_jobs=%s, allow_overlap=%s" % (cfg.get("max_jobs"), cfg.get("allow_overlap")), True)

    # 5. python deps
    try:
        import numpy  # noqa
        check("numpy", True)
    except ImportError:
        check("numpy", False, "pip install numpy  (needed by analysis scripts)", warn_only=True)
    try:
        import rdkit  # noqa
        check("rdkit (structure regeneration only)", True)
    except ImportError:
        check("rdkit (optional)", False, "pip install rdkit  (only needed to regenerate structures)", warn_only=True)

    # 6. run dirs
    has_runs = os.path.isdir(os.path.join(common.KIT, "02_calc"))
    check("02_calc run dirs %s" % ("(built)" if has_runs else "(not built yet)"), True,
          "" if has_runs else "run: python3 bin/make_inputs.py --tasks opt,opt_sol")

    print("\n%d PASS, %d WARN, %d FAIL" % (len(OKS), len(WARNS), len(FAILS)))
    if FAILS:
        print("NOT READY - fix the FAIL items above.")
        return 1
    print("READY.  quickstart:  bash run.sh   ->   python3 bin/submit.py --tasks opt,opt_sol")
    return 0


def glob_vasp():
    import glob
    return glob.glob(os.path.join(common.KIT, "00_structures", "*.vasp"))


def _rng(s):
    out = []
    for part in str(s).split(","):
        if "-" in part:
            a, b = part.split("-")
            out += [int(a), int(b)]
        elif part.strip().isdigit():
            out.append(int(part))
    return out or [0]


if __name__ == "__main__":
    sys.exit(main())
