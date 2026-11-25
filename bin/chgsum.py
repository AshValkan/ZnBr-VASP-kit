#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sum two VASP charge files (typically AECCAR0 + AECCAR2) into one CHGCAR.

Equivalent to VTST's chgsum.pl for the common case of two files with the same
grid.  Usage:

    python3 bin/chgsum.py AECCAR0 AECCAR2 > CHGCAR_sum

Why this exists: the bader binary needs a *total* charge reference, which is
the sum of the core (AECCAR0) and valence (AECCAR2) augmentation charges.  The
upstream Perl helper is no longer downloadable, and this does the same job in
~60 lines with no extra dependency.
"""

import sys


def _is_int(tok):
    try:
        int(tok)
        return True
    except ValueError:
        return False


def read_charge(path):
    """Return (header_lines, (nx, ny, nz), values, tail_lines).

    Parses the VASP CHGCAR layout: comment, scale, 3 lattice vectors, species,
    counts, coordinate mode, N coordinate lines, blank, grid dims, data.
    """
    with open(path, "r", errors="replace") as f:
        lines = f.readlines()

    if len(lines) < 9:
        raise ValueError("%s: too short to be a CHGCAR/AECCAR file" % path)

    # line 6 (0-based) holds the atom counts
    counts_line = 6
    try:
        natoms = sum(int(x) for x in lines[counts_line].split())
    except ValueError:
        raise ValueError("%s: cannot parse atom counts on line 7" % path)

    grid_line = counts_line + 2 + natoms          # mode line + coords
    while grid_line < len(lines) and not lines[grid_line].strip():
        grid_line += 1

    toks = lines[grid_line].split()
    if len(toks) != 3 or not all(_is_int(t) for t in toks):
        raise ValueError("%s: cannot locate grid dimensions (line %d: %r)"
                         % (path, grid_line + 1, lines[grid_line][:60]))
    ngrid = tuple(int(t) for t in toks)

    values = []
    for raw in lines[grid_line + 1:]:
        toks = raw.split()
        if not toks:
            break                                  # end of the data block
        values.extend(float(t) for t in toks)

    header = lines[:grid_line + 1]
    return header, ngrid, values, lines[grid_line + 1 + len(values):]


def main(argv):
    if len(argv) < 3:
        print("usage: chgsum.py FILE1 FILE2 [...] > OUT", file=sys.stderr)
        return 2

    header = None
    ngrid = None
    acc = None
    for path in argv[1:]:
        h, g, v, _ = read_charge(path)
        if ngrid is None:
            header, ngrid, acc = h, g, list(v)
        else:
            if g != ngrid:
                print("grid mismatch: %s has %s, expected %s" % (path, g, ngrid),
                      file=sys.stderr)
                return 1
            if len(v) != len(acc):
                print("length mismatch in %s: %d vs %d" % (path, len(v), len(acc)),
                      file=sys.stderr)
                return 1
            for i, x in enumerate(v):
                acc[i] += x

    sys.stdout.writelines(header)
    for i in range(0, len(acc), 5):
        sys.stdout.write("".join(" %0.11E" % x for x in acc[i:i + 5]) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
