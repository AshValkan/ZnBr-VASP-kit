#!/bin/bash
# run_bader.sh  —— 对已完成的 sp_vac 案例做 Bader 电荷分析
#
# 用法：
#   bash bin/run_bader.sh 02_calc/sp_vac/C10TAB_Br3_c1 [更多目录...]
#   bash bin/run_bader.sh --all-bader      # 自动处理 _hold_bader 里已跑完的三个
#
# 依赖：~/calc/bin/bader（Henkelman 组 1.05，已编译）、bin/chgsum.py
# 产出：CHGCAR_sum、ACF.dat、BCF.dat、bader.out、bader.csv
#
# 前提：该案例必须是 LAECHG=.TRUE. 跑的，目录里要有 AECCAR0 / AECCAR2。
set -u

KIT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$HOME/calc/bin:$PATH"

if ! command -v bader >/dev/null 2>&1; then
    echo "ERROR: bader 未安装（期望 ~/calc/bin/bader）" >&2
    exit 1
fi

if [ "${1:-}" = "--all-bader" ]; then
    # 只收实际存在、且确实带 AECCAR 的目录；_hold_bader 里的已归位
    set --
    for c in C10TAB_Br3_c1 C12TAB_Br3_c2 C14TAB_Br3_c1; do
        for base in "$KIT/02_calc/sp_vac" "$KIT/_hold_bader"; do
            if [ -f "$base/$c/AECCAR2" ]; then
                set -- "$@" "$base/$c"
            fi
        done
    done
fi

if [ $# -eq 0 ]; then
    sed -n '2,12p' "$0"
    exit 2
fi

fail=0
for d in "$@"; do
    echo "===== $(basename "$d") ====="
    if [ ! -d "$d" ]; then
        echo "  skip: 目录不存在"; continue
    fi

    missing=""
    for f in AECCAR0 AECCAR2 CHGCAR; do
        [ -s "$d/$f" ] || missing="$missing $f"
    done
    if [ -n "$missing" ]; then
        echo "  skip: 缺少$missing（该案例可能不是 LAECHG=.TRUE. 跑的）"
        fail=$((fail + 1))
        continue
    fi

    # 已有结果就跳过：一次重算要 15–20 分钟，重复跑纯属浪费
    if [ -s "$d/ACF.dat" ] && [ "$d/ACF.dat" -nt "$d/AECCAR2" ]; then
        echo "  skip: ACF.dat 已存在（$(stat -c '%y' "$d/ACF.dat" | cut -c1-16)）；要重算请先删掉它"
        continue
    fi

    echo "  chgsum: AECCAR0 + AECCAR2 -> CHGCAR_sum"
    python3 "$KIT/bin/chgsum.py" "$d/AECCAR0" "$d/AECCAR2" > "$d/CHGCAR_sum" || {
        echo "  chgsum 失败"; fail=$((fail + 1)); continue; }

    echo "  bader: CHGCAR -ref CHGCAR_sum"
    ( cd "$d" && bader CHGCAR -ref CHGCAR_sum > bader.out 2>&1 )
    if [ ! -s "$d/ACF.dat" ]; then
        echo "  bader 失败，见 $d/bader.out"
        tail -3 "$d/bader.out" 2>/dev/null | sed 's/^/    /'
        fail=$((fail + 1)); continue
    fi

    echo "  汇总:"
    python3 "$KIT/bin/bader_summary.py" "$d" 2>&1 | sed 's/^/    /'
    rm -f "$d/CHGCAR_sum"      # 中间产物，几 GB，用完即删
done

echo
if [ "$fail" -gt 0 ]; then
    echo "完成，$fail 个案例有问题"
    exit 1
fi
echo "全部完成"
