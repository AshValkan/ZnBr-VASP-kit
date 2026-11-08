#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reader for VASP CHGCAR / LOCPOT / AECCAR style grid files."""

import numpy as np


def read_grid(path, stride=1, nmax=None):
    """Return (dims, data[nz,ny,nx], meta).  data is subsampled by stride."""
    with open(path, encoding="utf-8", errors="replace") as f:
        comment = f.readline()
        scale = float(f.readline().split()[0])
        lat = [[float(x) * scale for x in f.readline().split()[:3]] for _ in range(3)]
        species = f.readline().split()
        counts = [int(x) for x in f.readline().split()]
        line = f.readline()
        while line.strip() == "":
            line = f.readline()
        direct = line.strip().lower().startswith("d")
        natoms = sum(counts)
        coords = []
        for _ in range(natoms):
            v = [float(x) for x in f.readline().split()[:3]]
            if direct:
                v = [v[0] * lat[0][0] + v[1] * lat[1][0] + v[2] * lat[2][0],
                     v[0] * lat[0][1] + v[1] * lat[1][1] + v[2] * lat[2][1],
                     v[0] * lat[0][2] + v[1] * lat[1][2] + v[2] * lat[2][2]]
            coords.append(v)
        line = f.readline()
        while line.strip() == "":
            line = f.readline()
        nx, ny, nz = [int(x) for x in line.split()[:3]]

        nsx = (nx + stride - 1) // stride
        nsy = (ny + stride - 1) // stride
        nsz = (nz + stride - 1) // stride
        data = np.zeros((nsz, nsy, nsx), dtype=np.float64)
        i = 0
        total = nx * ny * nz
        for tok in _tokens(f, total):
            ix = i % nx
            iy = (i // nx) % ny
            iz = i // (nx * ny)
            if ix % stride == 0 and iy % stride == 0 and iz % stride == 0:
                data[iz // stride, iy // stride, ix // stride] = tok
            i += 1
        return (nx, ny, nz), data, {"lattice": lat, "species": species,
                                    "counts": counts, "coords": coords}


def _tokens(f, total):
    n = 0
    for line in f:
        for tok in line.split():
            if n >= total:
                return
            yield float(tok)
            n += 1


def volume(lat):
    a = np.array(lat[0])
    b = np.array(lat[1])
    c = np.array(lat[2])
    return abs(np.dot(a, np.cross(b, c)))
