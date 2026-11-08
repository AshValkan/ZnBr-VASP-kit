#!/bin/bash
# ZnBr pipeline automation (cron every 10 min):
#   memory guard -> submit ready jobs -> collect energies
cd "$(dirname "$0")" || exit 1
AVAIL=$(free -m | awk 'NR==2{print $7}')
echo "=== autorun $(date) avail=${AVAIL}M ===" >> queue.log
if [ "${AVAIL:-0}" -lt 6000 ]; then
    echo "memory guard: only ${AVAIL}M available (<6000M), skip submit" >> queue.log
else
    /usr/bin/python3 bin/submit.py >> queue.log 2>&1
fi
/usr/bin/python3 bin/collect_energies.py >> queue.log 2>&1 || true
