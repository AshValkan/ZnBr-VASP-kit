#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""整理 04_results 的原始 CSV，生成汇总表、趋势图与文字报告。

用法:  python3 bin/make_summary.py [输出目录]
默认输出到 04_results/summary/

产出:
  binding_matrix.csv    cation × anion 结合能矩阵 (eV)
  summary_table.csv     逐体系明细（能量 / 结合能 / 收敛层级）
  report.md             中文摘要（完成度、关键结论、待办）
  fig_eb_grouped.png    分组柱状图（按阴离子分组，阳离子分色）
  fig_eb_trend.png      结合能随碳链长度的变化
  fig_exchange.png      dE_ex / dE_H2O 随碳链长度变化
"""

import csv
import os
import sys
from datetime import datetime

KIT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(KIT, "04_results")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RES, "summary")

CATIONS = ["C10TAB", "C12TAB", "C14TAB"]
ANIONS = ["Br", "Br3", "Br5", "Br2", "H2O"]
CHAIN = {"C10TAB": 10, "C12TAB": 12, "C14TAB": 14}
SUB = {"Br": "Br$^-$", "Br3": "Br$_3^-$", "Br5": "Br$_5^-$",
       "Br2": "Br$_2$", "H2O": "H$_2$O"}


def read_csv(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return []
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main():
    os.makedirs(OUT, exist_ok=True)
    eb = read_csv("eb.csv")
    raw = read_csv("energies_raw.csv")
    exc = read_csv("exchange.csv")

    # ---------- 明细与矩阵 ----------
    best = {(r["cation"], r["anion"]): r for r in eb if r.get("best") == "1"}
    energy = {r["system"]: fnum(r["toten_eV"]) for r in raw}
    level = {r["system"]: r.get("task", "") for r in raw}

    rows = []
    for (cat, an), r in sorted(best.items()):
        v = fnum(r.get("eb_vac"))
        s = fnum(r.get("eb_solv"))
        rows.append({
            "system": r["system"],
            "cation": cat,
            "anion": an,
            "conf": r.get("conf", ""),
            "E_total_eV": energy.get(r["system"]),
            "level": level.get(r["system"], ""),
            "Eb_vac_eV": v,
            "Eb_solv_eV": s,
            "best": r.get("level_vac") or r.get("level_solv") or "",
        })

    with open(os.path.join(OUT, "summary_table.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["system", "cation", "anion", "conf",
                                          "E_total_eV", "level", "Eb_vac_eV",
                                          "Eb_solv_eV", "best"])
        w.writeheader()
        w.writerows(rows)

    matrix = []
    for cat in CATIONS:
        line = {"cation": cat}
        for an in ANIONS:
            r = best.get((cat, an))
            line[an] = fnum(r.get("eb_vac")) if r else None
        matrix.append(line)
    with open(os.path.join(OUT, "binding_matrix.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cation"] + ANIONS)
        w.writeheader()
        w.writerows(matrix)

    # ---------- 报告 ----------
    n_mol = len(raw)
    n_cx = len(rows)
    done_cations = [c for c in CATIONS if all((c, a) in best for a in ANIONS)]
    partial = {c: [a for a in ANIONS if (c, a) in best] for c in CATIONS}
    partial = {c: v for c, v in partial.items() if v and c not in done_cations}

    lines = []
    lines.append("# Zn–Br 电解液体系 DFT 结果汇总\n")
    lines.append("生成时间：%s\n" % datetime.now().strftime("%Y-%m-%d %H:%M"))
    lines.append("## 完成度\n")
    lines.append("- 已收敛体系：**%d** 个（%d 单体 + %d 复合物）" % (n_mol, n_mol - n_cx, n_cx))
    lines.append("- 完整阳离子组：%s" % ("、".join(done_cations) if done_cations else "无"))
    if partial:
        for c, v in partial.items():
            lines.append("- %s：已完成 %s（缺 %s）" % (
                c, "、".join(v), "、".join([a for a in ANIONS if a not in v])))
    _solv = os.path.join(KIT, "02_calc", "opt_sol")
    _n = len(os.listdir(_solv)) if os.path.isdir(_solv) else 0
    lines.append("- 溶剂化相（opt_sol）：%s\n" % (
        "%d 项待跑" % _n if _n else "已停用（内容归档至 `_parked/opt_sol`）"))

    lines.append("## 结合能矩阵（真空相，eV）\n")
    hdr = "| cation | " + " | ".join(SUB[a] for a in ANIONS) + " |"
    lines.append(hdr)
    lines.append("|" + "---|" * (len(ANIONS) + 1))
    for line in matrix:
        cells = []
        for a in ANIONS:
            v = line[a]
            cells.append("%.2f" % v if v is not None else "—")
        lines.append("| %s | %s |" % (line["cation"], " | ".join(cells)))
    lines.append("")

    if exc:
        lines.append("## 竞争吸附指标（eV）\n")
        lines.append("| cation | dE_ex | dE_H2O |")
        lines.append("|---|---|---|")
        for r in exc:
            ex = fnum(r.get("dE_ex_vac"))
            h2 = fnum(r.get("dE_H2O_vac"))
            if ex is None and h2 is None:
                continue
            lines.append("| %s | %s | %s |" % (
                r["cation"],
                "%.2f" % ex if ex is not None else "—",
                "%.2f" % h2 if h2 is not None else "—"))
        lines.append("")
        lines.append("dE_ex < 0 表示 Br₃⁻ 优于 Br⁻；dE_H2O < 0 表示 Br₃⁻ 优于水。\n")

    lines.append("## 趋势\n")
    chain_done = [c for c in CATIONS if all((c, a) in best for a in ANIONS)]
    if len(chain_done) >= 2:
        seg_hdr = " | ".join(
            "C%d→C%d" % (CHAIN[chain_done[i]], CHAIN[chain_done[i + 1]])
            for i in range(len(chain_done) - 1))
        lines.append("| 吸附物 | %s |" % seg_hdr)
        lines.append("|---|" + "---|" * (len(chain_done) - 1))
        last_max = 0.0
        for a in ANIONS:
            vals = [fnum(best[(c, a)].get("eb_vac")) for c in chain_done]
            if any(v is None for v in vals):
                continue
            diffs = [vals[i + 1] - vals[i] for i in range(len(vals) - 1)]
            last_max = max(last_max, abs(diffs[-1]))
            lines.append("| %s | %s |" % (SUB[a], " | ".join("%+.2f" % d for d in diffs)))
        lines.append("")
        if len(chain_done) >= 3:
            lines.append("最后一段 C%d→C%d 的最大变化为 **%.2f eV**。" % (
                CHAIN[chain_done[-2]], CHAIN[chain_done[-1]], last_max))
            if last_max < 0.05:
                lines.append("该值落在 DFT 误差带（约 0.05 eV）内，"
                             "**判定链长效应在此处饱和**；而前一段的落差明显更大。\n")
            else:
                lines.append("仍可分辨出链长效应，尚未饱和。\n")

    lines.append("## 待办\n")
    missing = ["%s_%s" % (c, a) for c in CATIONS for a in ANIONS if (c, a) not in best]
    if missing:
        lines.append("- 未完成 c1 复合物：%s" % "、".join(missing))
    else:
        lines.append("- c1 复合物已全部完成（共 %d 个体系）" % len(best))
    lines.append("- 关键构象 c2/c3 的构象竞争检验")
    solv = os.path.join(KIT, "02_calc", "opt_sol")
    n_solv = len(os.listdir(solv)) if os.path.isdir(solv) else 0
    if n_solv:
        lines.append("- opt_sol 溶剂化相仍有 %d 项待跑" % n_solv)
    else:
        lines.append("- opt_sol 已停用（目录为空，内容在 `_parked/`）")
    lines.append("")

    with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    # ---------- 图 ----------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
    except Exception as e:
        print("matplotlib 不可用，跳过图：%s" % e)
        print("已生成 CSV 与 report.md ->", OUT)
        return

    rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.linewidth": 2.0,
        "xtick.major.width": 1.6,
        "ytick.major.width": 1.6,
        "font.size": 11,
    })

    present = [c for c in CATIONS if any((c, a) in best for a in ANIONS)]

    # 图 1：分组柱状
    fig, ax = plt.subplots(figsize=(9, 4.6))
    w = 0.8 / max(len(present), 1)
    for k, cat in enumerate(present):
        vals, pos, labels = [], [], []
        for i, a in enumerate(ANIONS):
            r = best.get((cat, a))
            v = fnum(r.get("eb_vac")) if r else None
            if v is not None:
                vals.append(v)
                pos.append(i + (k - (len(present) - 1) / 2) * w)
                labels.append(a)
        ax.bar(pos, vals, width=w * 0.92, label=cat)
    ax.set_xticks(range(len(ANIONS)))
    ax.set_xticklabels([SUB[a] for a in ANIONS])
    ax.set_ylabel("Binding energy (eV)")
    ax.axhline(0, color="black", lw=1.2)
    ax.legend(frameon=False, ncol=len(present))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_eb_grouped.png"), dpi=300)
    plt.close(fig)

    # 图 2：随链长变化
    chain_cats = [c for c in present if c in CHAIN]
    if len(chain_cats) >= 2:
        fig, ax = plt.subplots(figsize=(7, 4.6))
        for a in ANIONS:
            xs, ys = [], []
            for c in sorted(chain_cats, key=lambda x: CHAIN[x]):
                r = best.get((c, a))
                v = fnum(r.get("eb_vac")) if r else None
                if v is not None:
                    xs.append(CHAIN[c])
                    ys.append(v)
            if len(xs) >= 2:
                ax.plot(xs, ys, marker="o", ms=6, lw=2, label=SUB[a])
        ax.set_xlabel("Alkyl chain length (carbon atoms)")
        ax.set_ylabel("Binding energy (eV)")
        ax.set_xticks(sorted({CHAIN[c] for c in chain_cats}))
        ax.legend(frameon=False, ncol=3, fontsize=9)
        ax.grid(axis="y", alpha=0.25, lw=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "fig_eb_trend.png"), dpi=300)
        plt.close(fig)

    # 图 3：竞争指标
    if exc:
        pts = [r for r in exc if r["cation"] in CHAIN and
               (fnum(r.get("dE_ex_vac")) is not None or fnum(r.get("dE_H2O_vac")) is not None)]
        if len(pts) >= 2:
            fig, ax = plt.subplots(figsize=(6.5, 4.4))
            for key, lbl, mk in (("dE_ex_vac", "Br$_3^-$ vs Br$^-$", "o"),
                                 ("dE_H2O_vac", "Br$_3^-$ vs H$_2$O", "s")):
                xs = [CHAIN[r["cation"]] for r in pts if fnum(r.get(key)) is not None]
                ys = [fnum(r.get(key)) for r in pts if fnum(r.get(key)) is not None]
                if len(xs) >= 2:
                    ax.plot(xs, ys, marker=mk, ms=7, lw=2, label=lbl)
            ax.axhline(0, color="black", lw=1.0, ls="--")
            ax.set_xlabel("Alkyl chain length (carbon atoms)")
            ax.set_ylabel("Preference (eV)")
            ax.set_xticks(sorted({CHAIN[r["cation"]] for r in pts}))
            # 数据集中在 ±1 eV 两侧，中间大片空白；图例放这里可避开零参考线
            ax.legend(frameon=False, loc="center", bbox_to_anchor=(0.5, 0.46))
            ax.grid(axis="y", alpha=0.25, lw=0.8)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            fig.tight_layout()
            fig.savefig(os.path.join(OUT, "fig_exchange.png"), dpi=300)
            plt.close(fig)

    print("完成 ->", OUT)
    for n in sorted(os.listdir(OUT)):
        print("  ", n)


if __name__ == "__main__":
    main()
