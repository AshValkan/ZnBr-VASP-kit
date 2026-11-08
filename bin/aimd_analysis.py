#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AIMD analysis: confinement of Br3-, RDF, potential-energy fluctuation.

    python3 bin/aimd_analysis.py 02_calc/aimd/C12TAB_Br3_c1 [--skip 1000]

Outputs (next to XDATCAR):
  aimd_timeseries.csv  t(ps), msd(A^2), d_br_cat(A), epot(eV)
  rdf_NBr.csv, rdf_BrBr.csv
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common


def read_xdatcar(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    elems = []
    for s, c in zip(species, counts):
        elems += [s] * c
    natoms = len(elems)
    frames = []
    i = 7
    while i < len(lines):
        if lines[i].startswith("Direct configuration"):
            arr = np.array([[float(x) for x in lines[i + 1 + k].split()[:3]]
                            for k in range(natoms)])
            frames.append(arr @ lat)
            i += 1 + natoms
        else:
            i += 1
    return lat, elems, frames


def read_epot(path):
    import re
    out = []
    if not os.path.exists(path):
        return out
    pat = re.compile(r"E0=s*([-+]?\d+\.\d+E[-+]?\d+)")
    for line in open(path, encoding="utf-8", errors="replace"):
        if "E0=" in line:
            m = pat.search(line)
            if m:
                out.append(float(m.group(1)))
    return out


def mic(d, lat):
    """Minimum image in a general (possibly triclinic) cell."""
    inv = np.linalg.inv(np.asarray(lat))
    f = np.asarray(d) @ inv
    f = f - np.round(f)
    return f @ np.asarray(lat)


def rdf(pairs, lat, n_center, n_frames, n_partner, rmax=8.0, nbins=80):
    hist = np.zeros(nbins)
    for d in pairs:
        r = np.linalg.norm(mic(d, lat))
        if r < rmax:
            k = int(r / rmax * nbins)
            if k < nbins:
                hist[k] += 1
    edges = np.linspace(0, rmax, nbins + 1)
    r = 0.5 * (edges[1:] + edges[:-1])
    vol = abs(np.linalg.det(lat))
    rho = n_partner / vol
    shell = 4 * np.pi * r ** 2 * (rmax / nbins)
    g = hist / (n_frames * n_center * shell * rho + 1e-12)
    return r, g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--skip", type=int, default=0,
                    help="frames to discard as equilibration (default 10%%)")
    args = ap.parse_args()
    d = args.dir
    lat, elems, frames = read_xdatcar(os.path.join(d, "XDATCAR"))
    n = len(frames)
    if n < 10:
        print("not enough frames:", n)
        return
    skip = args.skip or max(1, n // 10)
    epot = read_epot(os.path.join(d, "OSZICAR"))
    br = [i for i, e in enumerate(elems) if e == "Br"]
    nidx = [i for i, e in enumerate(elems) if e == "N"]
    cat = [i for i in range(len(elems)) if i not in br]

    com_br = np.array([f[br].mean(0) for f in frames])
    com_cat = np.array([f[cat].mean(0) for f in frames])
    # unwrap the Br3- centre of mass relative to the cation
    rel = np.zeros_like(com_br)
    rel[0] = mic(com_br[0] - com_cat[0], lat)
    for i in range(1, n):
        step = mic((com_br[i] - com_cat[i]) - (com_br[i - 1] - com_cat[i - 1]), lat)
        rel[i] = rel[i - 1] + step
    d0 = rel[skip].copy()
    disp = np.linalg.norm(rel - d0, axis=1)
    msd = disp ** 2

    potim = 1.0
    inc = os.path.join(d, "INCAR")
    if os.path.exists(inc):
        for line in open(inc, encoding="utf-8", errors="replace"):
            if line.strip().upper().startswith("POTIM"):
                try:
                    potim = float(line.split("=")[1].split()[0])
                except (IndexError, ValueError):
                    pass
    t = np.arange(n) * potim * 1e-3  # ps
    out = os.path.join(d, "aimd_timeseries.csv")
    with open(out, "w", encoding="utf-8") as f:
        f.write("frame,t_ps,msd_A2,disp_A,d_br_cat_A,epot_eV\n")
        for i in range(n):
            e = epot[i] if i < len(epot) else float("nan")
            f.write("%d,%.4f,%.4f,%.4f,%.4f,%.6f\n"
                    % (i, t[i], msd[i], disp[i],
                       np.linalg.norm(mic(com_br[i] - com_cat[i], lat)), e))

    pairs_nbr, pairs_brbr = [], []
    sampled = list(range(skip, n, max(1, (n - skip) // 200)))
    for i in sampled:
        fr = frames[i]
        for a in nidx:
            for b in br:
                pairs_nbr.append(fr[a] - fr[b])
        for a in br:
            for b in br:
                if a != b:
                    pairs_brbr.append(fr[a] - fr[b])
    r1, g1 = rdf(pairs_nbr, lat, max(1, len(nidx)), len(sampled), len(br))
    r2, g2 = rdf(pairs_brbr, lat, len(br), len(sampled), len(br))
    np.savetxt(os.path.join(d, "rdf_NBr.csv"),
               np.column_stack([r1, g1]), delimiter=",", header="r_A,g_NBr", comments="")
    np.savetxt(os.path.join(d, "rdf_BrBr.csv"),
               np.column_stack([r2, g2]), delimiter=",", header="r_A,g_BrBr", comments="")

    e = np.array(epot[-n + skip:]) if len(epot) >= n else np.array(epot[skip:])
    print("frames           %d  (equilibration %d)" % (n, skip))
    print("max |dr|         %.2f A" % disp[skip:].max())
    print("final |dr|       %.2f A" % disp[-1])
    print("mean d(Br-cat)   %.2f A" % np.mean([np.linalg.norm(mic(com_br[i] - com_cat[i], lat))
                                              for i in range(skip, n)]))
    if e.size:
        print("std(E_pot)       %.4f eV" % e.std())
    print("written", out)


if __name__ == "__main__":
    main()
