#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Submit ready calculations, honouring a global queue limit.

    python3 bin/submit.py --tasks opt,opt_sol --dry-run
    python3 bin/submit.py --tasks sp_vac,freq --max 30

Adds a "direct" scheduler for standalone hosts (no Slurm): picks a free
taskset core block per job, launches it in the background with JOB_CORES
set, and records the block in .submitted.
"""

import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common


def _ranges_to_set(s):
    out = set()
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a), int(b) + 1))
        elif part.strip().isdigit():
            out.add(int(part))
    return out


def vasp_rank_masks():
    """CPU affinity masks of running vasp_std ranks (taskset exec()s into
    prterun, so /proc is the only reliable source)."""
    import glob
    masks = []
    for p in glob.glob("/proc/[0-9]*"):
        try:
            with open(os.path.join(p, "comm"), encoding="ascii", errors="ignore") as f:
                if f.read().strip() != "vasp_std":
                    continue
            with open(os.path.join(p, "status"), encoding="ascii", errors="ignore") as f:
                for ln in f:
                    if ln.startswith("Cpus_allowed_list"):
                        masks.append(_ranges_to_set(ln.split("\t")[-1].strip()))
                        break
        except Exception:
            continue
    return masks


def block_usage(cfg):
    """[(block, nranks)] -- ranks whose affinity mask intersects each block."""
    usage = []
    masks = vasp_rank_masks()
    for b in cfg.get("core_blocks") or []:
        s = _ranges_to_set(b)
        usage.append((b, sum(1 for m in masks if m & s)))
    return usage


def free_blocks(cfg, extra_busy=()):
    """Blocks with zero vasp ranks; extra_busy = blocks claimed this run."""
    if not (cfg.get("core_blocks") or []):
        return [None]
    claimed = set(extra_busy)
    return [b for b, n in block_usage(cfg) if n == 0 and b not in claimed]


def launch_direct(cfg, d, block, dry):
    env = dict(os.environ)
    if block:
        env["JOB_CORES"] = block
    cmd = ["bash", "job.sh"]
    if dry:
        print("[dry] cd %s && JOB_CORES=%s bash job.sh (background)"
              % (os.path.relpath(d, common.KIT), block))
        return "dry"
    log = open(os.path.join(d, "submit.log"), "a", encoding="utf-8")
    if cfg.get("launcher") == "systemd":
        import time as _t
        unit = "znbr-%s-%d" % (os.path.basename(d), int(_t.time()))
        cmd2 = ["sudo", "-n", "/usr/bin/systemd-run", "--unit", unit,
                "--uid=%d" % os.getuid(), "--gid=%d" % os.getgid(),
                "--working-directory", d,
                "--setenv=HOME=%s" % os.path.expanduser("~"),
                "--setenv=PATH=/usr/local/bin:/usr/local/sbin:/usr/bin:/bin"]
        if block:
            cmd2.append("--setenv=JOB_CORES=%s" % block)
        cmd2 += ["bash", "job.sh"]
        r = subprocess.run(cmd2, stdout=log, stderr=log, env=env)
        return "unit %s" % unit if r.returncode == 0 else "systemd-run FAILED rc=%d" % r.returncode
    p = subprocess.Popen(cmd, cwd=d, env=env, stdout=log, stderr=log,
                         start_new_session=True)
    return "pid %d" % p.pid


def mem_available_mb():
    """MemAvailable in MB (0 on failure)."""
    try:
        with open("/proc/meminfo") as f:
            for ln in f:
                if ln.startswith("MemAvailable:"):
                    return int(ln.split()[1]) // 1024
    except Exception:
        pass
    return 0


def queue_key(name):
    """Order: singles first, then one conformer of every complex (c1), then c2, c3.

    Running all c1 first yields a complete Eb matrix early; extra conformers
    only refine it, so they can wait behind the full coverage set.

    """
    m = re.search(r"_c(\d+)$", name)
    return (int(m.group(1)) if m else 0, name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", default="opt,opt_sol")
    ap.add_argument("--max", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--include-partial", action="store_true",
                    help="also (re)submit cases that already have run artefacts")
    args = ap.parse_args()

    cfg = common.load_config()
    sched = cfg.get("scheduler", "direct").lower()
    if sched == "local":      # alias kept for the upstream kit / smoke test
        sched = "direct"
    if sched == "local":      # alias kept for the upstream kit / smoke test
        sched = "direct"
    max_jobs = args.max if args.max is not None else int(cfg["max_jobs"])
    active = common.active_jobs()
    room = max(0, max_jobs - active)
    print("queue: %d running, %d slots free (max %d)" % (active, room, max_jobs))
    if sched == "direct" and room > 0:
        fb = free_blocks(cfg)
        print("core blocks: free=%s" % (fb or ["<none>"]))

    tasks = [t.strip() for t in args.tasks.split(",") if t.strip()]
    submitted = 0
    locally_assigned = []   # track blocks picked this run (async start race)

    # pre-pass: conformer policy + case state (fresh / partial / done)
    wl = cfg.get("conformer_whitelist") or []
    dropped = 0
    state_of = {}
    n_fresh = n_partial = n_done = 0
    manual = []
    for task in tasks:
        td = os.path.join(common.KIT, "02_calc", task)
        if not os.path.isdir(td):
            continue
        for name in os.listdir(td):
            d = os.path.join(td, name)
            if not os.path.isdir(d):
                continue
            if not common.conformer_allowed(cfg, name):
                dropped += 1
                continue
            st = common.run_state(d)
            state_of[d] = st
            if st == "done":
                n_done += 1
            elif st == "partial":
                n_partial += 1
                manual.append(os.path.relpath(d, common.KIT))
            else:
                n_fresh += 1
    print("conformer policy: c1 for all systems, c2/c3 only for %d key cases (%d dirs skipped)"
          % (len(wl), dropped))
    print("case states: fresh %d | partial %d (left for manual) | done %d"
          % (n_fresh, n_partial, n_done))

    for task in tasks:
        td = os.path.join(common.KIT, "02_calc", task)
        if not os.path.isdir(td):
            continue
        for name in sorted(os.listdir(td), key=queue_key):
            d = os.path.join(td, name)
            if not os.path.isdir(d):
                continue
            if not common.conformer_allowed(cfg, name):
                continue
            st = state_of.get(d) or common.run_state(d)
            if st == "done" and not args.force:
                continue
            if st == "partial" and not (args.force or args.include_partial):
                # already started: never touched automatically, operator decides
                continue
            if not os.path.exists(os.path.join(d, "POTCAR")):
                print("skip (no POTCAR):", os.path.relpath(d, common.KIT))
                continue
            if submitted >= room:
                print("queue limit reached; stopping")
                return
            avail = mem_available_mb()
            floor = int(cfg.get("mem_floor_mb", 6000))
            if avail and avail < floor:
                print("memory floor reached (avail %d MB < %d MB); stopping" % (avail, floor))
                return
            if sched == "direct":
                usage = dict(block_usage(cfg))
                for b in locally_assigned:
                    usage[b] = usage.get(b, 0) + int(cfg.get("nprocs", 1))
                free = [b for b in usage if usage[b] == 0]
                if free:
                    block = free[0]
                elif usage and cfg.get("allow_overlap", False):
                    # run freely: share the least-loaded block instead of stalling
                    block = min(usage, key=lambda b: usage[b])
                elif usage:
                    print("no free core block; stopping")
                    return
                else:
                    block = None
                locally_assigned.append(block)
                tag = launch_direct(cfg, d, block, args.dry_run)
                info = "block %s %s" % (block, tag)
            else:
                cmd = {"slurm": ["sbatch", "job.sh"],
                       "pbs": ["qsub", "job.sh"]}.get(sched, ["bash", "job.sh"])
                if args.dry_run:
                    print("[dry] cd %s && %s" % (os.path.relpath(d, common.KIT), " ".join(cmd)))
                    info = "dry"
                else:
                    p = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True)
                    out = (p.stdout or "").strip()
                    print("%-40s %s" % (os.path.relpath(d, common.KIT), out))
                    info = out
            if not args.dry_run:
                with open(os.path.join(d, ".submitted"), "w") as f:
                    f.write(info + "\n")
            print("%-40s %s" % (os.path.relpath(d, common.KIT), info))
            submitted += 1
    print("\n%d jobs submitted" % submitted)
    if manual:
        print("%d case(s) already started -> left for manual handling "
              "(--include-partial to override):" % len(manual))
        for m in manual[:15]:
            print("  " + m)
        if len(manual) > 15:
            print("  ... and %d more" % (len(manual) - 15))


if __name__ == "__main__":
    main()
