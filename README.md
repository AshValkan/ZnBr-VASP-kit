# Zn–Br VASP Kit

水系锌溴电池理论计算工具链（VASP · PBE + D3(BJ) · VASPsol）。
config 驱动、manifest 驱动：68 个体系（11 单体 + 57 复合物构象）从结构到能量到分析一条流水线。

## 目录

```
README.md          本文件（快速上手）
config.json        全局配置：路径 / 调度器 / 所有 INCAR 参数（ENCUT、IVDW、VASPsol、AIMD…）
run.sh             一键助手：make_inputs + submit --dry-run
bin/               全部 python 工具（见下表）
docs/checklist.md  提交前自检清单
docs/charter.md    立项章程（参数依据 §3 / 任务定义 §4 / INCAR 模板 §10）
00_structures/     68 个初始结构 <name>.vasp + manifest.json（含 charge/natoms/关键原子索引）
                   + cases_manifest.csv（生成清单）
02_calc/           运行目录（make_inputs 生成）：02_calc/<task>/<name>/
results/           汇总输出（collect_energies / plot_results 生成）
```

## bin/ 工具一览

| 脚本 | 作用 |
|---|---|
| `generate_structures.py` | 从 SMILES + 内坐标生成 68 案例（RDKit ETKDG/MMFF94；黄金螺旋选向防碰撞） |
| `verify_structures.py` | 几何校验：内距 / 角度 / N···Br 目标距 / 居中 / 最小接触 |
| `export_kit_manifest.py` | 校验过的案例 → `00_structures/*.vasp` + `manifest.json` |
| `make_inputs.py` | 建运行目录：POSCAR + POTCAR 拼接 + INCAR/KPOINTS + 作业脚本 |
| `doctor.py` | 新机器环境体检：路径/赝势/二进制/依赖/结构清单一键自检 |
| `submit.py` | 批量提交：三态队列（只提交 fresh）、并发上限按实际进程数、`--dry-run` 预览、direct 模式核块分配 |
| `status.py` | 一屏进度（done / partial / fresh + 运行中）；`--brief` 供轮询脚本调用 |
| `regen_jobsh.py` | 按当前模板重写全部 `job.sh`（自动跳过正在运行的目录） |
| `pin_ranks.sh` | 由 `job.sh` 调用：把每个 rank 钉到核块内一对连续 CPU |
| `prepare_next.py` | 阶段链：opt CONTCAR → sp / freq / aimd |
| `collect_energies.py` | Eb（真空+溶剂）+ 交换能 + 水竞争 + `results/eb.md` |
| `freq_thermo.py` | 频率 → 零点/熵/自由能校正 |
| `aimd_analysis.py` | MSD / RDF / 能量漂移 |
| `esp_stats.py` `locpot2cube.py` | LOCPOT → cube → 静电势统计 |
| `bader_summary.py` | Bader 电荷汇总 |
| `plot_results.py` | 结果出图 |

## 快速上手

```bash
# 0) 开箱两步（新机器）
#    a. 编辑 config.json 的三个路径: vasp_std / mpi_cmd / potcar_dir
#    b. 环境体检（9 项自检，不跑 VASP）:
python3 bin/doctor.py

# 1) 生成 T1 输入（opt = 真空, opt_sol = VASPsol）并预览队列
bash run.sh                       # make_inputs + submit --dry-run

# 2) 确认无误后正式提交
python3 bin/submit.py --tasks opt,opt_sol
python3 bin/status.py             # 盯进度

# 3) 收能量
python3 bin/collect_energies.py   # → results/eb.csv, eb.md

# 4) 阶段链（freq/aimd 等）
python3 bin/prepare_next.py --from opt --to freq
python3 bin/submit.py --tasks freq
```

### 无人值守（可选）

```bash
crontab -e   # 加一行（每10分钟自动续队 + 自动收能量）:
*/10 * * * * /bin/bash /path/to/kit/autorun.sh
```

### 调度行为（direct 模式）

- 优先使用空闲核心块；`allow_overlap=true` 时若无空闲块则选**最闲块共享**（按运行中 rank 数均衡），`false` 则等待
- `max_jobs` 限制本 kit 并发作业数；`.submitted` 标记幂等防重，重启后自动恢复
- 串行保守配置：`max_jobs=1, allow_overlap=false`（内存受限机器推荐）；放开：`max_jobs=4~8, allow_overlap=true`

### 队列范围（构象策略）

```json
"conformers_default": 1,
"conformer_whitelist": ["C12TAB_Br3_c2", "C12TAB_Br3_c3", "..."]
```

- 默认**只跑每个复合物的 c1**（构象 1），其余 c2/c3 不进队列
- `conformer_whitelist` 里的案例名允许跑到更高的构象号，用于关键结论的稳健性检验
- 结构仍然全部生成（68 个 POSCAR），只是不提交计算；后续要补算把名字加进白名单即可
- 同一套规则由 `bin/submit.py`（提交）与 `poll_status.sh`（进度统计）共用，
  避免"被跳过的目录永远算作待跑"导致队列清空判定失效

### 半途任务不自动接管

`bin/submit.py` 把每个案例分成三态：

| 状态 | 判据 | 自动队列行为 |
|---|---|---|
| `done` | OUTCAR 含 `reached required accuracy` | 跳过 |
| `partial` | 有 `.submitted` 标记，或目录里已有 OUTCAR/OSZICAR/CONTCAR/vasp.out/CHG 等运行产物 | **不碰**，留给人工 |
| `fresh` | 只有 POSCAR/POTCAR/INCAR/KPOINTS/job.sh | 正常提交 |

也就是说，**跑过一半的案例永远不会被自动队列重跑或覆盖**；要接管需显式
`python3 bin/submit.py --include-partial`（保留原目录与 CONTCAR，从断点续算）。
并发上限按**实际运行的 prterun 进程数**统计，而不是数 `.submitted` 标记，
这样崩溃残留的标记不会把队列卡死。

## 关键自动化

- **NELECT / ISPIN**：按 `manifest.json` 的 `charge` 自动写 `NELECT = Σ(n·ZVAL) − charge`，
  电子数奇数自动 `ISPIN=2`（阳离子 +1、Br⁻/Br₃⁻/Br₅⁻ −1、阳离子+Br₂/H₂O 复合物 +1）。
  结论只建立在 ΔE > 0.2 eV 之上（带电胞误差预算，见章程 §3.6）。
- **scheduler**：`config.json → scheduler`。
  - `direct`（默认，单机无 Slurm）：扫描运行中 `taskset` CPU 集，自动把作业
    分配到 `core_blocks` 中与所有占用块不相交的块，后台 setsid 启动；
  - `slurm` / `pbs`：模板化 sbatch/qsub（换集群时切换）。
- **防重**：`.submitted` 标记 + 队列上限 `max_jobs`；重跑安全。
- **从头再生成结构**：`python3 bin/generate_structures.py . && python3 bin/verify_structures.py .`
  然后 `python3 bin/export_kit_manifest.py .`（需要 RDKit）。

## 依赖

- 运行：numpy（esp/aimd 分析需要）、VASP 6.4.2（含 VASPsol，判定看 OUTCAR 回读 `LSOL=T` + `Solvation Ediel_sol`）
- 结构再生成（可选）：rdkit
- 赝势：PBE PAW，N/C/H/Br/O，`config.json → potcar_dir` 指向按元素分目录的 PP 根

## 适配记录（相对上游 ZnBr_VASP_kit）

1. 补齐电荷/自旋自动处理（上游缺失，本体系科学正确性的关键）
2. 新增 `direct` 调度器（taskset 核块隔离 + 竞态防护）
3. `vasp_gam` → `vasp_std`；阳离子命名 MEP/TPA → MEPBr/TPABr
4. 结构矩阵扩为 68 案例（上游 43 构象），全部过几何校验
