#!/usr/bin/env python3
"""Rewrite job.sh for every existing case directory from the current template.

job.sh is a pure derivative of config.json (vasp binary, np, mpi_cmd, launcher),
so regenerating it never touches INCAR / POSCAR / POTCAR / results.
"""
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import common          # noqa: E402
import vasp_inputs     # noqa: E402

ROOT = os.path.dirname(HERE)
cfg = common.load_config()


def running_dirs():
    """Directories that currently hold a live vasp_std / prterun process.

    job.sh is being read by bash while the job runs; rewriting it under a running
    shell can make bash resume at a stale byte offset and re-execute a line.
    """
    live = set()
    for name in ("vasp_std", "prterun", "mpirun"):
        try:
            out = subprocess.run(["pgrep", "-x", name],
                                 capture_output=True, text=True).stdout
        except Exception:
            continue
        for pid in out.split():
            try:
                live.add(os.path.realpath("/proc/%s/cwd" % pid))
            except OSError:
                pass
    return live


LIVE = running_dirs()

n = 0
changed = 0
skipped = 0
for d in sorted(glob.glob(os.path.join(ROOT, "02_calc", "*", "*"))):
    if not os.path.isdir(d) or not os.path.exists(os.path.join(d, "INCAR")):
        continue
    n += 1
    if os.path.realpath(d) in LIVE:
        skipped += 1
        print("skip (job running): %s" % os.path.relpath(d, ROOT))
        continue
    name = os.path.basename(d)
    task = "opt"
    tfile = os.path.join(d, "TASK")
    if os.path.exists(tfile):
        parts = open(tfile, encoding="utf-8").read().split()
        if parts:
            task = parts[0]
    new = vasp_inputs.job_script(cfg, name, task)
    p = os.path.join(d, "job.sh")
    old = ""
    if os.path.exists(p):
        old = open(p, encoding="utf-8").read()
    if old != new:
        with open(p, "w", encoding="utf-8", newline="\n") as f:
            f.write(new)
        os.chmod(p, 0o755)
        changed += 1

print("case dirs: %d, job.sh rewritten: %d, skipped running: %d" % (n, changed, skipped))

sample = os.path.join(ROOT, "02_calc", "opt", "C10TAB_H2O_c1", "job.sh")
if os.path.exists(sample):
    print("--- sample ---")
    print(open(sample, encoding="utf-8").read())
