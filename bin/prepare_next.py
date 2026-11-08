#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""After a geometry optimisation finishes, feed its CONTCAR to the dependent jobs.

    python3 bin/prepare_next.py            # opt -> sp_vac, freq ; opt_sol -> sp_sol, aimd
"""

import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

DEPENDENTS = {
    "opt": ["sp_vac", "freq"],
    "opt_sol": ["sp_sol", "aimd"],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    root = os.path.join(common.KIT, "02_calc")
    n = 0
    for src_task, dst_tasks in DEPENDENTS.items():
        sd = os.path.join(root, src_task)
        if not os.path.isdir(sd):
            continue
        for name in sorted(os.listdir(sd)):
            d = os.path.join(sd, name)
            contcar = os.path.join(d, "CONTCAR")
            if not os.path.exists(contcar) or os.path.getsize(contcar) < 100:
                continue
            if not common.is_done(d):
                print("skip (opt not finished):", name)
                continue
            for dst_task in dst_tasks:
                td = os.path.join(root, dst_task, name)
                if not os.path.isdir(td):
                    continue
                pos = os.path.join(td, "POSCAR")
                if os.path.exists(os.path.join(td, ".from_opt")) and not args.force:
                    continue
                shutil.copyfile(contcar, pos)
                with open(os.path.join(td, ".from_opt"), "w") as f:
                    f.write(src_task + "\n")
                # remove stale submit marker so it can be queued again
                p = os.path.join(td, ".submitted")
                if os.path.exists(p):
                    os.remove(p)
                print("copied CONTCAR -> %s/%s" % (dst_task, name))
                n += 1
    print("\n%d dependent directories updated" % n)


if __name__ == "__main__":
    main()
