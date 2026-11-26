#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Local electrostatic-potential descriptors on the rho = 1e-3 e/A^3 surface.

    python3 bin/esp_local.py <calc_dir> [stride] [R_head]

The global sigma from esp_stats.py turned out to be useless for the chain-length
question: it is a whole-surface variance, dominated by the cation as a whole, and
does not move when the head-group charge is diluted (C10 1.037 / C12 0.984 /
C14 1.069 eV while Eb changes by 0.46 eV).

This script samples instead at *specific places*, so the descriptors can be
compared against where the anion actually sits:

  head_all      mean V_S over surface points within R_head of the quaternary N
  head_hemi     same, restricted to the hemisphere facing the anion
  anion_site    V at the anion centroid (mean over a small sphere)
  anion_surf    mean V_S over surface points within 2.5 A of any Br

Units: eV, already referenced to the vacuum level.
"""

import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
import gridio

RHO0 = 1.0e-3
TOL = 0.25


def grid_cart(dims, lattice, stride):
    """Cartesian coordinates for the subsampled grid, shape (nz,ny,nx,3)."""
    nx, ny, nz = dims
    nsx = (nx + stride - 1) // stride
    nsy = (ny + stride - 1) // stride
    nsz = (nz + stride - 1) // stride
    fx = (np.arange(nsx) * stride) / float(nx)
    fy = (np.arange(nsy) * stride) / float(ny)
    fz = (np.arange(nsz) * stride) / float(nz)
    a1, a2, a3 = (np.array(v, dtype=float) for v in lattice[:3])
    X, Y, Z = np.meshgrid(fx, fy, fz, indexing="xy")  # (nsz,nsy,nsx)
    cart = (X[..., None] * a1 + Y[..., None] * a2 + Z[..., None] * a3)
    return cart


def split_species(meta):
    """meta['coords'] is a flat list; rebuild per-element lists."""
    species = meta["species"]
    counts = meta["counts"]
    coords = meta["coords"]
    out = {}
    k = 0
    for s, c in zip(species, counts):
        out.setdefault(s, []).extend(coords[k:k + c])
        k += c
    return out


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    d = sys.argv[1]
    stride = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    r_head = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0

    dims, pot, meta = gridio.read_grid(os.path.join(d, "LOCPOT"), stride)
    _, chg, _ = gridio.read_grid(os.path.join(d, "CHGCAR"), stride)
    lat = np.array(meta["lattice"], dtype=float)
    vol = abs(np.linalg.det(lat))
    rho = chg / vol

    cart = grid_cart(dims, lat, stride)

    vac = pot[rho < 1e-6]
    ref = float(vac.mean()) if vac.size else 0.0
    v = pot - ref

    sel = np.abs(rho - RHO0) < TOL * RHO0
    if sel.sum() < 100:
        sel = np.abs(rho - RHO0) < 3 * TOL * RHO0

    atoms = split_species(meta)
    n_pos = np.array(atoms.get("N", [np.nan, np.nan, np.nan])[0], dtype=float)
    br = np.array(atoms.get("Br", []), dtype=float)

    pts = cart[sel]
    vals = v[sel]
    out = {"dir": os.path.basename(d.rstrip("/")), "points": int(sel.sum()),
           "ref_eV": ref}

    # --- head region ---
    if not np.isnan(n_pos).any():
        dist = np.linalg.norm(pts - n_pos, axis=1)
        m = dist < r_head
        out["head_all"] = float(vals[m].mean()) if m.sum() else float("nan")
        out["head_n"] = int(m.sum())
        if len(br):
            axis = br.mean(axis=0) - n_pos
            axis /= np.linalg.norm(axis)
            proj = (pts - n_pos) @ axis
            m2 = m & (proj > 0)
            out["head_hemi"] = float(vals[m2].mean()) if m2.sum() else float("nan")
            out["head_hemi_n"] = int(m2.sum())
    # --- anion ---
    if len(br):
        c = br.mean(axis=0)
        near = np.linalg.norm(cart - c, axis=-1) < 1.5
        out["anion_site"] = float(v[near].mean()) if near.sum() else float("nan")
        dist_br = np.linalg.norm(pts[:, None, :] - br[None, :, :], axis=2).min(axis=1)
        m3 = dist_br < 2.5
        out["anion_surf"] = float(vals[m3].mean()) if m3.sum() else float("nan")
        out["anion_surf_n"] = int(m3.sum())

    path = os.path.join(d, "esp_local.csv")
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value"])
        for k, val in out.items():
            w.writerow([k, ("%.4f" % val) if isinstance(val, float) else val])

    print("%-20s pts=%-5d head=%-8s hemi=%-8s site=%-8s surf=%-8s" % (
        out["dir"], out["points"],
        ("%.3f" % out["head_all"]) if "head_all" in out else "—",
        ("%.3f" % out["head_hemi"]) if "head_hemi" in out else "—",
        ("%.3f" % out["anion_site"]) if "anion_site" in out else "—",
        ("%.3f" % out["anion_surf"]) if "anion_surf" in out else "—"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
