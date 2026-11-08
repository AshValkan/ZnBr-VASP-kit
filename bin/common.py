#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared helpers: config, POSCAR I/O, OUTCAR parsing.  Stdlib only.

Runs on a standalone Linux host: direct scheduler (MPI + taskset core blocks,
no Slurm) plus NELECT/ISPIN auto handling.
"""

import json
import os
import re
import subprocess
import sys

KIT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DEFAULT_CONFIG = {
    "vasp_std": "~/calc/vasp_std",
    "vasp_gam": "~/calc/vasp_std",   # fallback: Gamma-only runs may use vasp_std
    "mpi_cmd": "/usr/local/bin/mpirun --bind-to none",
    "potcar_dir": "~/calc/VASP-MXenes-alloy/PP",
    "potcar_map": {"C": "C", "H": "H", "N": "N", "Br": "Br", "O": "O"},
    "scheduler": "direct",           # direct | slurm | pbs
    "partition": "normal",
    "nprocs": 5,
    "nodes": 1,
    "max_jobs": 4,
    "core_blocks": ["0-9", "10-19", "20-29", "30-39"],
    "modules": [],
    "walltime": {"opt": "48:00:00", "opt_sol": "48:00:00",
                 "sp": "08:00:00", "freq": "96:00:00", "aimd": "192:00:00"},
    "encut": 500,
    "sigma": 0.05,
    "ivdw": 12,
    "ediff": 1e-5,
    "ediffg": -0.02,
    "nsw": 300,
    "ncore": 1,
    "solvent": {"EB_K": 78.54, "LAMBDA_D_k": 1.0e4, "SIGMA_K": 0.6,
                "TAU": 0.0, "NC_K": 1.0e-4},
    "aimd": {"temperature": 298, "steps": 10000, "potim": 1.0},
    "masses": {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "Br": 79.904},
    "zval": {"H": 1, "C": 4, "N": 5, "O": 6, "Br": 7},
}


def load_config(path=None):
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))
    path = path or os.path.join(KIT, "config.json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            user = json.load(f)
        for k, v in user.items():
            if isinstance(v, dict) and isinstance(cfg.get(k), dict):
                cfg[k].update(v)
            else:
                cfg[k] = v
    if not os.path.isabs(cfg["potcar_dir"]):
        cfg["potcar_dir"] = os.path.expanduser(cfg["potcar_dir"])
    for k in ("vasp_std", "vasp_gam"):
        if cfg[k].startswith("~"):
            cfg[k] = os.path.expanduser(cfg[k])
    return cfg


# ---------------------------------------------------------------- POSCAR ----
def read_poscar(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    scale = float(lines[1].split()[0])
    lat = [[float(x) * scale for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    i = 5
    if lines[i].strip() and not lines[i].strip()[0].isdigit():
        species = lines[i].split()
        i += 1
    else:
        species = ["X"]
    counts = [int(x) for x in lines[i].split()]
    i += 1
    if lines[i].strip().lower().startswith("s"):
        i += 1
    mode = lines[i].strip().lower()
    i += 1
    direct = mode.startswith("d")
    coords = []
    for k in range(sum(counts)):
        parts = lines[i + k].split()
        v = [float(parts[0]), float(parts[1]), float(parts[2])]
        if direct:
            v = [v[0] * lat[0][0] + v[1] * lat[1][0] + v[2] * lat[2][0],
                 v[0] * lat[0][1] + v[1] * lat[1][1] + v[2] * lat[2][1],
                 v[0] * lat[0][2] + v[1] * lat[1][2] + v[2] * lat[2][2]]
        coords.append(v)
    elems = []
    for s, c in zip(species, counts):
        elems += [s] * c
    return {"scale": scale, "lattice": lat, "species": species,
            "counts": counts, "coords": coords, "elems": elems,
            "comment": lines[0]}


def write_poscar(path, comment, lat, elems, coords):
    species, counts = [], []
    for e in elems:
        if not species or species[-1] != e:
            species.append(e)
            counts.append(1)
        else:
            counts[-1] += 1
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(comment + "\n1.0\n")
        for v in lat:
            f.write("  %.10f %.10f %.10f\n" % tuple(v))
        f.write(" ".join(species) + "\n")
        f.write(" ".join(str(c) for c in counts) + "\nDirect\n")
        for v in coords:
            f.write("  %.10f %.10f %.10f\n" % tuple(v))


# ---------------------------------------------------- charge / electrons ----
def nelect_of(poscar_path, cfg, charge):
    """NELECT = sum(n_species * ZVAL) - system charge."""
    pos = read_poscar(poscar_path)
    z = cfg["zval"]
    total = 0
    for s, c in zip(pos["species"], pos["counts"]):
        if s not in z:
            raise KeyError("ZVAL missing for element %r" % s)
        total += c * z[s]
    return int(round(total - charge))


# ---------------------------------------------------------------- OUTCAR ----
_FLOAT = r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?"


def outcar_energy(path):
    """Return dict(toten, sigma0, finished, converged) or None if unreadable."""
    if not os.path.exists(path):
        return None
    toten = None
    sigma0 = None
    finished = False
    converged = False
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if "free  energy   TOTEN" in line:
                    m = re.search(_FLOAT, line.split("=")[-1])
                    if m:
                        toten = float(m.group())
                elif "energy(sigma->0)" in line:
                    m = re.search(_FLOAT, line.split("=")[-1])
                    if m:
                        sigma0 = float(m.group())
                elif "reached required accuracy" in line:
                    converged = True
                elif "General timing and accounting" in line:
                    finished = True
    except OSError:
        return None
    if toten is None:
        return None
    return {"toten": toten, "sigma0": sigma0,
            "finished": finished, "converged": converged}


def is_done(d):
    e = outcar_energy(os.path.join(d, "OUTCAR"))
    return bool(e and e["finished"])


# ---------------------------------------------------------------------------
# case state: one definition shared by submit.py, status.py and the poller
# ---------------------------------------------------------------------------
RUN_MARKERS = ("OUTCAR", "OSZICAR", "CONTCAR", "vasp.out", "vasprun.xml",
               "EDD", "DOSCAR", "EIGENVAL", "XDATCAR", "CHG", "CHGCAR", "WAVECAR")


def run_state(d):
    """Classify a relaxation case directory: 'done' / 'partial' / 'fresh'.

    done    : OUTCAR reports 'reached required accuracy'
    partial : started at least once (submission marker or any run artefact)
    fresh   : input files only, never launched

    Only 'fresh' cases may be submitted automatically; 'partial' ones are left
    to the operator so an unattended run can never redo work in progress.
    (Single-point tasks never print the convergence line, so they always read
    as 'partial' -- submit.py only queues opt/opt_sol, which is intended.)
    """
    e = outcar_energy(os.path.join(d, "OUTCAR"))
    if e and e["converged"]:
        return "done"
    if os.path.exists(os.path.join(d, ".submitted")):
        return "partial"
    for f in RUN_MARKERS:
        p = os.path.join(d, f)
        try:
            if os.path.exists(p) and os.path.getsize(p) > 0:
                return "partial"
        except OSError:
            pass
    return "fresh"


def conformer_allowed(cfg, name):
    """Conformer policy: c1 for every complex, c2/c3 only for whitelisted cases."""
    m = re.search(r"_c(\d+)$", name)
    if not m:
        return True
    if int(m.group(1)) <= int(cfg.get("conformers_default", 1)):
        return True
    return name in set(cfg.get("conformer_whitelist") or [])


def active_jobs():
    """Number of VASP jobs actually running (not inferred from markers)."""
    try:
        out = subprocess.run(["pgrep", "-x", "prterun"],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             text=True, timeout=10).stdout
        return len([x for x in out.split() if x.strip()])
    except Exception:
        return 0


def run(cmd, cwd=None, dry=False):
    if dry:
        print("[dry] " + " ".join(cmd))
        return 0
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    if p.stdout:
        sys.stdout.write(p.stdout)
    return p.returncode
