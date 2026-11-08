#!/bin/bash
# pin_ranks.sh BLOCK NP
#
# Pin each of the NP VASP ranks of the job running in the current directory to
# one physical core inside BLOCK.  This box enumerates one physical core as two
# consecutive logical CPUs (core c -> cpus 2c, 2c+1), so BLOCK "10-19" yields
# pairs 10-11 / 12-13 / 14-15 / 16-17 / 18-19.
#
# Called in the background by job.sh right after mpirun starts: mpirun itself is
# launched with --bind-to none because PRRTE ignores the taskset mask.
set -u

BLOCK="${1:-}"
NP="${2:-5}"
HERE="$(pwd)"
DEADLINE=$(( $(date +%s) + 150 ))

if [ -z "$BLOCK" ]; then
    echo "usage: pin_ranks.sh BLOCK NP"
    exit 2
fi

mapfile -t CPUS < <(python3 -c '
import sys
out = []
for part in sys.argv[1].split(","):
    part = part.strip()
    if not part:
        continue
    if "-" in part:
        a, b = part.split("-")
        out.extend(range(int(a), int(b) + 1))
    else:
        out.append(int(part))
print("\n".join(str(c) for c in sorted(set(out))))
' "$BLOCK")

if [ ${#CPUS[@]} -eq 0 ]; then
    echo "bad block: $BLOCK"
    exit 2
fi

PAIRS=()
i=0
while [ $i -lt ${#CPUS[@]} ]; do
    if [ $((i + 1)) -lt ${#CPUS[@]} ]; then
        PAIRS+=("${CPUS[$i]}-${CPUS[$((i + 1))]}")
    else
        PAIRS+=("${CPUS[$i]}")
    fi
    i=$((i + 2))
done

echo "block=$BLOCK np=$NP pairs=${PAIRS[*]} dir=$HERE"

while [ "$(date +%s)" -lt "$DEADLINE" ]; do
    PIDS=()
    for p in $(pgrep -x vasp_std 2>/dev/null); do
        d=$(readlink "/proc/$p/cwd" 2>/dev/null) || continue
        [ "$d" = "$HERE" ] && PIDS+=("$p")
    done
    if [ ${#PIDS[@]} -ge "$NP" ]; then
        mapfile -t SORTED < <(printf '%s\n' "${PIDS[@]}" | sort -n)
        k=0
        for p in "${SORTED[@]}"; do
            [ $k -ge ${#PAIRS[@]} ] && break
            if taskset -pc "${PAIRS[$k]}" "$p" >/dev/null 2>&1; then
                echo "rank $p -> ${PAIRS[$k]}"
            fi
            k=$((k + 1))
        done
        echo "pinned $k/$NP ranks"
        exit 0
    fi
    sleep 2
done

echo "timeout: only found ${#PIDS[@]}/$NP ranks, left unbound inside $BLOCK"
exit 1
