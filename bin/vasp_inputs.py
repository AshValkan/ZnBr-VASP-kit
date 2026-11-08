#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""INCAR / KPOINTS / job-script builders.

Adds NELECT/ISPIN auto handling (charged anion / cation systems) and a
"direct" scheduler template for standalone hosts (MPI --bind-to none +
taskset core block + per-rank pinning, no Slurm).
"""

import os

KPOINTS_GAMMA = """Gamma-point only
0
Gamma
1 1 1
0 0 0
"""


def _fmt(v):
    if isinstance(v, bool):
        return ".TRUE." if v else ".FALSE."
    if isinstance(v, float):
        s = ("%.6g" % v)
        if "e" in s or "E" in s:
            m, e = s.replace("E", "e").split("e")
            return "%sE%s" % (m, e.replace("+", ""))
        return s
    return str(v)


def _block(d):
    return "\n".join("%-10s = %s" % (k, _fmt(v)) for k, v in d.items())


def incar_opt(cfg, name, solv=False):
    d = {
        "SYSTEM": name,
        "PREC": "Accurate",
        "ENCUT": int(cfg["encut"]),
        "EDIFF": "1E-5",
        "EDIFFG": cfg["ediffg"],
        "IBRION": 2,
        "ISIF": 2,
        "NSW": int(cfg["nsw"]),
        "ISYM": 0,
        "ISPIN": 1,
        "ALGO": "Normal",
        "LREAL": cfg.get("lreal", False),
        "ISMEAR": 0,
        "SIGMA": cfg["sigma"],
        "IVDW": int(cfg["ivdw"]),
        "NCORE": int(cfg["ncore"]),
        "LWAVE": False,
        "LCHARG": False,
    }
    if solv:
        s = cfg["solvent"]
        d.update({"LSOL": True, "EB_K": s["EB_K"], "LAMBDA_D_k": s["LAMBDA_D_k"],
                  "SIGMA_K": s["SIGMA_K"], "TAU": s["TAU"], "NC_K": s["NC_K"]})
    return _block(d) + "\n"


def _common(cfg, name):
    return {
        "SYSTEM": name,
        "PREC": "Accurate",
        "ENCUT": int(cfg["encut"]),
        "ISYM": 0,
        "ISPIN": 1,
        "LREAL": cfg.get("lreal", False),
        "ISMEAR": 0,
        "SIGMA": cfg["sigma"],
        "IVDW": int(cfg["ivdw"]),
        "NCORE": int(cfg["ncore"]),
        "LWAVE": False,
    }


def incar_sp(cfg, name, solv=False, locpot=True):
    d = _common(cfg, name)
    d.update({"IBRION": -1, "NSW": 0, "EDIFF": "1E-6",
              "LCHARG": True, "LAECHG": True,
              "LVTOT": bool(locpot), "LVHAR": False})
    if solv:
        s = cfg["solvent"]
        d.update({"LSOL": True, "EB_K": s["EB_K"], "LAMBDA_D_k": s["LAMBDA_D_k"],
                  "SIGMA_K": s["SIGMA_K"], "TAU": s["TAU"], "NC_K": s["NC_K"]})
    return _block(d) + "\n"


def incar_freq(cfg, name, solv=False):
    d = _common(cfg, name)
    d.update({"IBRION": 5, "NFREE": 2, "POTIM": 0.015, "EDIFF": "1E-6",
              "EDIFFG": -0.005, "ADDGRID": True, "NWRITE": 0,
              "LCHARG": False, "LVTOT": False})
    if solv:
        s = cfg["solvent"]
        d.update({"LSOL": True, "EB_K": s["EB_K"], "LAMBDA_D_k": s["LAMBDA_D_k"],
                  "SIGMA_K": s["SIGMA_K"], "TAU": s["TAU"], "NC_K": s["NC_K"]})
    return _block(d) + "\n"


def incar_aimd(cfg, name, solv=False):
    a = cfg["aimd"]
    d = {
        "SYSTEM": name,
        "PREC": "Accurate",
        "ENCUT": int(cfg["encut"]),
        "ISYM": 0,
        "ISPIN": 1,
        "LREAL": cfg.get("lreal", False),
        "ISMEAR": 0,
        "SIGMA": cfg["sigma"],
        "IVDW": int(cfg["ivdw"]),
        "NCORE": int(cfg["ncore"]),
        "LWAVE": False,
        "LCHARG": False,
        "IBRION": 0,
        "NSW": int(a["steps"]),
        "POTIM": float(a["potim"]),
        "TEBEG": int(a["temperature"]),
        "TEEND": int(a["temperature"]),
        "MDALGO": 2,
        "SMASS": 0,
        "EDIFF": "1E-4",
        "ALGO": "Fast",
        "NBLOCK": 10,
        "KBLOCK": 100,
        "NWRITE": 0,
    }
    if solv:
        s = cfg["solvent"]
        d.update({"LSOL": True, "EB_K": s["EB_K"], "LAMBDA_D_k": s["LAMBDA_D_k"],
                  "SIGMA_K": s["SIGMA_K"], "TAU": s["TAU"], "NC_K": s["NC_K"]})
    return _block(d) + "\n"


SLURM = """#!/bin/bash
#SBATCH -J {name}
#SBATCH -N {nodes}
#SBATCH --ntasks-per-node {nprocs}
#SBATCH -p {partition}
#SBATCH -t {walltime}
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
{modules}

cd $SLURM_SUBMIT_DIR
echo "start $(date)  $(pwd)"
srun {vasp} > vasp.out 2>&1
echo "stop  $(date)"
"""

PBS = """#!/bin/bash
#PBS -N {name}
#PBS -l nodes={nodes}:ppn={nprocs}
#PBS -l walltime={walltime}
#PBS -j oe
{modules}

cd $PBS_O_WORKDIR
echo "start $(date)  $(pwd)"
mpirun -np {nprocs} {vasp} > vasp.out 2>&1
echo "stop  $(date)"
"""

def job_script(cfg, name, task):
    """Direct template needs the mpirun line built here; slurm/pbs as in kit."""
    sched = cfg.get("scheduler", "direct").lower()
    if sched == "local":      # alias kept for the upstream kit / smoke test
        sched = "direct"
    if sched == "local":      # alias kept for the upstream kit / smoke test
        sched = "direct"
    vasp = cfg["vasp_gam"]   # every task here is Gamma-only
    if sched == "direct":
        mpi = cfg.get("mpi_cmd", "mpirun --bind-to none")
        vasp_sh = cfg["vasp_gam"]
        # NOTE: mpirun must NOT do the binding itself.  PRRTE builds its map from
        # the whole-machine topology and ignores the taskset mask, so a job
        # restricted to 10-19 still lands on cores 0-9 -- two concurrent jobs then
        # share the same physical cores.  Instead: taskset confines every rank to
        # the block, and pin_ranks.sh spreads one rank per physical core inside it.
        return """#!/bin/bash
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=$HOME/calc/pmix-5.0.1/lib:/usr/local/lib:$LD_LIBRARY_PATH
export PATH=/usr/local/bin:$PATH
export OMP_NUM_THREADS=1
echo "start $(date)  cores=${{JOB_CORES:-none}}  supervisor=$$  $(pwd)"

# Supervised loop: VASP has been observed to die right after an ionic step
# (twice, no OOM / no systemd event).  Restarting it blindly used to throw away
# the whole run, so every attempt first rewinds to the last CONTCAR.
attempt=0
max_attempts=30
while [ $attempt -lt $max_attempts ]; do
    attempt=$((attempt+1))
    if [ -s CONTCAR ] && ! grep -q 'reached required accuracy' OUTCAR 2>/dev/null; then
        nat=$(awk 'NR==7{{s=0; for(i=1;i<=NF;i++) s+=$i; print s}}' POSCAR)
        lines=$(wc -l < CONTCAR)
        if [ -n "$nat" ] && [ "$lines" -ge $((nat+8)) ]; then
            cp -f CONTCAR POSCAR
            echo "attempt $attempt: POSCAR <- CONTCAR ($lines lines, $nat atoms)"
        fi
    fi
    [ $attempt -gt 1 ] && [ -f vasp.out ] && mv -f vasp.out "vasp.out.attempt$((attempt-1))"
    echo "attempt $attempt: launch $(date)"
    if [ -n "$JOB_CORES" ]; then
        taskset -c "$JOB_CORES" {mpi} -np {np} {vasp} > vasp.out 2>&1 &
        MPI_PID=$!
        timeout 180 bash "$HOME/calc/ZnBr/bin/pin_ranks.sh" "$JOB_CORES" {np} > pin.log 2>&1 &
        wait $MPI_PID
        rc=$?
    else
        {mpi} -np {np} {vasp} > vasp.out 2>&1
        rc=$?
    fi
    echo "attempt $attempt: vasp exit rc=$rc $(date)"
    if grep -q 'reached required accuracy' OUTCAR 2>/dev/null; then
        echo "converged after $attempt attempt(s)"
        break
    fi
    sleep 5
done
echo "stop  $(date)  attempts=$attempt  supervisor=$$"
""".format(mpi=mpi, np=int(cfg["nprocs"]), vasp=vasp_sh)
    w = cfg["walltime"].get(task, "24:00:00")
    modules = "\n".join(cfg.get("modules", []))
    nprocs = int(cfg["nprocs"])
    if sched == "slurm":
        s = SLURM.format(name=name[:32], nodes=int(cfg["nodes"]), nprocs=nprocs,
                         partition=cfg["partition"], walltime=w,
                         modules=modules, vasp=vasp)
    elif sched == "pbs":
        s = PBS.format(name=name[:15], nodes=int(cfg["nodes"]), nprocs=nprocs,
                       walltime=w, modules=modules, vasp=vasp)
    else:
        raise ValueError("unknown scheduler: %s" % sched)
    return s


def write_nelect(d, cfg, poscar_path, charge):
    """Append NELECT to INCAR; flip ISPIN=2 when the electron count is odd."""
    nelect = __import__("common").nelect_of(poscar_path, cfg, charge)
    inc = os.path.join(d, "INCAR")
    txt = open(inc, encoding="utf-8").read()
    txt += "NELECT    = %d\n" % nelect
    if nelect % 2 == 1:
        txt = txt.replace("ISPIN      = 1", "ISPIN      = 2")
    with open(inc, "w", encoding="utf-8", newline="\n") as f:
        f.write(txt)
    return nelect


def write_inputs(d, cfg, name, task, solv, charge=0, poscar_path=None):
    os.makedirs(d, exist_ok=True)
    if task.startswith("opt"):
        inc = incar_opt(cfg, name, solv)
    elif task.startswith("sp"):
        inc = incar_sp(cfg, name, solv, locpot=True)
    elif task.startswith("freq"):
        inc = incar_freq(cfg, name, solv)
    elif task.startswith("aimd"):
        inc = incar_aimd(cfg, name, solv)
    else:
        raise ValueError(task)
    with open(os.path.join(d, "INCAR"), "w", encoding="utf-8", newline="\n") as f:
        f.write(inc)
    if poscar_path is not None:
        # always write NELECT (even for charge=0): explicit value is auditable,
        # and odd-electron systems (e.g. MEPBr_Br3, 69 e) must get ISPIN=2.
        write_nelect(d, cfg, poscar_path, charge)
    with open(os.path.join(d, "KPOINTS"), "w", encoding="utf-8", newline="\n") as f:
        f.write(KPOINTS_GAMMA)
    js = job_script(cfg, name, task)
    with open(os.path.join(d, "job.sh"), "w", encoding="utf-8", newline="\n") as f:
        f.write(js)
    os.chmod(os.path.join(d, "job.sh"), 0o755)
    # task metadata used by the analysis scripts
    with open(os.path.join(d, "TASK"), "w", encoding="utf-8") as f:
        f.write("%s %s\n" % (task, "solv" if solv else "vac"))
