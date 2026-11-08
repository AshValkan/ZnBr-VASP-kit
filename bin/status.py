#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""One-screen status for every calculation directory.

    python3 bin/status.py            # full table
    python3 bin/status.py --brief    # compact block for the poller / automation

Case states come from common.run_state(): done / partial / fresh.  Anything
outside the conformer policy is reported separately and never counted as work.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

MPI_NAMES = ("vasp_std", "prterun", "mpirun")


def ionic_steps(d):
    p = os.path.join(d, "OSZICAR")
    if not os.path.exists(p):
        return 0
    n = 0
    for line in open(p, encoding="utf-8", errors="replace"):
        if "F=" in line and "E0=" in line:
            n += 1
    return n


def _pids(name):
    try:
        import subprocess
        out = subprocess.run(["pgrep", "-x", name], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True, timeout=10).stdout
        return [int(x) for x in out.split()]
    except Exception:
        return []


def _expand(mask):
    s = set()
    for part in mask.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-")
            s.update(range(int(a), int(b) + 1))
        else:
            s.add(int(part))
    return s


def _mask_of(pid):
    try:
        for ln in open("/proc/%d/status" % pid):
            if ln.startswith("Cpus_allowed_list"):
                return ln.split()[1]
    except OSError:
        pass
    return ""


def running_jobs():
    """[{rel, ionic, cores, ranks}] for VASP jobs currently running under 02_calc."""
    jobs = {}
    for name in MPI_NAMES:
        for pid in _pids(name):
            try:
                cwd = os.readlink("/proc/%d/cwd" % pid)
            except OSError:
                continue
            if "/02_calc/" in cwd:
                jobs.setdefault(cwd, [])
                if name == "vasp_std":
                    jobs[cwd].append(pid)
    out = []
    for cwd, ps in sorted(jobs.items()):
        rel = cwd.split("/02_calc/")[-1]
        cores = set()
        for p in ps:
            cores |= {c // 2 for c in _expand(_mask_of(p))}
        out.append({"rel": rel, "ionic": ionic_steps(cwd),
                    "cores": sorted(cores), "ranks": len(ps)})
    return out


def survey(cfg):
    """Count in-scope cases by state; return counts, manual list, skipped, rows."""
    root = os.path.join(common.KIT, "02_calc")
    counts = {"done": 0, "partial": 0, "fresh": 0}
    manual, skipped, rows = [], 0, []
    if not os.path.isdir(root):
        return counts, manual, skipped, rows
    for task in sorted(os.listdir(root)):
        td = os.path.join(root, task)
        if not os.path.isdir(td):
            continue
        for name in sorted(os.listdir(td)):
            d = os.path.join(td, name)
            if not os.path.isdir(d):
                continue
            rows.append((task, name, d))
            if not common.conformer_allowed(cfg, name):
                skipped += 1
                continue
            st = common.run_state(d)
            counts[st] += 1
            if st == "partial":
                manual.append("%s/%s" % (task, name))
    return counts, manual, skipped, rows


def brief(jobs, counts, manual, skipped):
    for j in jobs:
        print("job %-28s ionic=%-4d ranks=%-2d phys_cores=%s"
              % (j["rel"], j["ionic"], j["ranks"], j["cores"]))
        if j["ionic"] >= 60:
            print("ZNBR_SLOW (%s at %d ionic steps, flat-PES drift -> "
                  "consider IBRION=1 or EDIFFG -0.08)" % (j["rel"], j["ionic"]))

    overlap = []
    for i in range(len(jobs)):
        for k in range(i + 1, len(jobs)):
            c = set(jobs[i]["cores"]) & set(jobs[k]["cores"])
            if c:
                overlap.append((jobs[i]["rel"], jobs[k]["rel"], sorted(c)))
    if overlap:
        for a, b, c in overlap:
            print("CPU_OVERLAP %s <-> %s share physical cores %s" % (a, b, c))
    else:
        used = sum(len(j["cores"]) for j in jobs)
        print("cpu binding: OK (%d jobs, %d physical cores, no overlap)"
              % (len(jobs), used))

    print("scope %d: done %d | partial %d (manual) | fresh %d   [%d c2/c3 skipped]"
          % (sum(counts.values()), counts["done"], counts["partial"],
             counts["fresh"], skipped))
    if manual:
        shown = ", ".join(manual[:8])
        more = "" if len(manual) <= 8 else " ... +%d" % (len(manual) - 8)
        print("manual: " + shown + more)

    if not jobs:
        if counts["fresh"] == 0 and counts["partial"] > 0:
            print("ZNBR_MANUAL_ONLY (auto queue drained; %d half-run case(s) need "
                  "manual handling)" % counts["partial"])
        elif counts["fresh"] == 0:
            print("ZNBR_QUEUE_EMPTY (all in-scope cases finished)")
        else:
            print("ZNBR_STALLED (fresh cases pending but nothing running)")


def table(cfg, jobs, rows):
    live = {j["rel"] for j in jobs}
    print("%-9s %-24s %-9s %6s %14s" % ("task", "system", "state", "steps", "E (eV)"))
    print("-" * 68)
    for task, name, d in rows:
        if not common.conformer_allowed(cfg, name):
            continue
        st = common.run_state(d)
        e = common.outcar_energy(os.path.join(d, "OUTCAR"))
        rel = "%s/%s" % (task, name)
        mark = "running" if rel in live else st
        en = "%.4f" % e["toten"] if e and e.get("toten") is not None else "-"
        print("%-9s %-24s %-9s %6d %14s" % (task, name, mark, ionic_steps(d), en))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--brief", action="store_true",
                    help="compact output for the poller / automation")
    args = ap.parse_args()

    cfg = common.load_config()
    jobs = running_jobs()
    counts, manual, skipped, rows = survey(cfg)

    if args.brief:
        brief(jobs, counts, manual, skipped)
    else:
        table(cfg, jobs, rows)
        print("-" * 68)
        print("scope %d   done %d   partial %d (manual)   fresh %d   skipped %d"
              % (sum(counts.values()), counts["done"], counts["partial"],
                 counts["fresh"], skipped))
        print("running %d   %s"
              % (len(jobs), ", ".join(j["rel"] for j in jobs) or "-"))


if __name__ == "__main__":
    main()
