#!/bin/bash
# ZnBr pipeline automation (cron every 10 min):
#   prepare dependents -> submit ready jobs -> collect energies
# 内存守卫已移除（用户要求）：swap 足够，接受内存压力换吞吐。
# 任务链：opt 收敛 -> sp_vac（ESP：LVTOT/LAECHG，产出 LOCPOT+CHGCAR）
cd "$(dirname "$0")" || exit 1
echo "=== autorun $(date) avail=$(free -m | awk 'NR==2{print $7}')M ===" >> queue.log
/usr/bin/python3 bin/prepare_next.py >> queue.log 2>&1 || true
/usr/bin/python3 bin/submit.py >> queue.log 2>&1
/usr/bin/python3 bin/collect_energies.py >> queue.log 2>&1 || true
