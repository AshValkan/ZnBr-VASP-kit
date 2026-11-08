# 季铵盐–多溴阴离子体系理论计算方案（VASP 版）

> **版本** v1.0　|　**日期** 2026-10-26　|　**计算平台** VASP 6.x（PAW/PBE + DFT-D3 + VASPsol）+ 开源后处理
> **适用体系** 水系静态 Zn–Br₂ 电池正极：C₁₀/C₁₂/C₁₄TAB⁺ 及对照 BCA（MEPBr⁺、TPABr⁺、EMImBr⁺）
> **关联文献** ①双咪唑 BCA（Angew. Chem. 2025，anie202509293）；②疏水有机溴化物 M0–M6（Small Methods 2025，smtd70112）；③疏水季铵碘化物 C₈–C₁₈I₃（Angew. Chem. 2026，anie74000）

---

## 0. 结论先行

1. **本方案保留原 DMol3 方案的三个核心计算**（阳离子–多溴结合能、阳离子–水结合能、静电势 ESP），全部迁移到 VASP 的 PAW/PBE + D3(BJ) + 隐式溶剂框架，并补充**交换能判据、显式水微溶剂化、Bader/DDEC6 电荷、ΔG 热化学、AIMD 限域动力学**五项，使"静电锚定 + 疏水屏障"从定性图像变成可量化、可复核的证据链。
2. **关键方法学改动**：①用"Br₃⁻/Br⁻ 交换能"替代直接比较两个带电单体的绝对结合能，从源头消除周期带电胞的背景电荷误差；②用 Bader/DDEC6 替代 Mulliken/Hirshfeld（VASP 无原生 Mulliken，且 Mulliken 在平面波基下无定义）；③ESP 用 LOCPOT + 真空参考 + Multiwfn 定量表面分析，出图统一色标与等值面；④ΔG 用"全 Hessian + 头基部分 Hessian"两级策略控制机时。
3. **预期结果**：Eb(C₁₂TAB⁺–Br₃⁻) 约 −0.5 ~ −0.7 eV，比 EMImBr⁺ 更负 ≥ 0.15 eV；Eb(C₁₂TAB⁺–H₂O) 仅 −0.1 ~ −0.25 eV；C₁₂/C₁₄ 的 N⁺ 区域正电势最强；AIMD 中 Br₃⁻ 质心位移 8–12 Å，落在 EMIm⁺（>15 Å）与 TPA⁺（<6 Å）之间。
4. **一条必须提前想清楚的科学逻辑**：计算很可能给出 C₁₄ 结合最强，但实验最优是 C₁₂。论文不能把"最强结合"直接等同于"最优 BCA"，必须引入"结合强度–疏水性–离子传输"的权衡描述符（见 §1.6），否则 C₁₄ > C₁₂ 的结果会与"C₁₂ 最优"的结论自相矛盾。

---

## 1. 目标与论证逻辑

### 1.1 一句话思路

用长链季铵溴（C₁₀/C₁₂/C₁₄TAB）作为溴络合剂，与 Br₃⁻/Br₅⁻ 预络合形成疏水离子对/固体配合物做自支撑正极；季铵阳离子靠**静电作用锚定多溴阴离子**，长烷基链提供**疏水屏障**，二者协同抑制多溴溶解与穿梭、保护锌负极。计算要证明：**C₁₂TAB⁺ 抓得住多溴、同时排斥水**。

### 1.2 论证链

~~~
C12TAB+ 与 Br3- 结合强（强于与 H2O 结合）
        │
        ├─→ 优先锚定多溴，而不是被水溶剂化
        │
        ├─→ 多溴不溶解、不穿梭到负极
        │
        ├─→ 锌负极不被腐蚀、库伦效率高
        │
        └─→ 电池循环稳定
~~~

三个计算分别承担链条上的一环，缺一不可：

| 计算 | 回答的问题 | 在链条中的位置 |
|---|---|---|
| ① 阳离子–多溴结合能 Eb | 抓得牢不牢？ | 锚定能力 |
| ② 阳离子–水结合能 Eb(H₂O) | 更愿意抓多溴还是抓水？ | 选择性 |
| ③ 静电势 ESP | 正电荷在哪？为什么能静电锚定？ | 机理解释 |

### 1.3 三个核心计算的任务定义

| 编号 | 计算 | 体系 | 输出量 | 判据（越负越强） |
|---|---|---|---|---|
| C1 | 阴阳离子结合能 | 阳离子 + Br⁻/Br₃⁻/Br₅⁻/Br₂ | Eb（真空、溶剂化两套） | C₁₂ 比 EMIm⁺ 更负；Br₃⁻/Br₅⁻ 明显强于 Br₂ |
| C2 | 与水的结合能 | 阳离子 + H₂O | Eb(H₂O) | 明显弱于 Eb(Br₃⁻)，差值 > 0.2 eV |
| C3 | 静电势 | 阳离子（及复合物） | V_S,max、V_S,min、正电势面积占比、ESP 图 | C₁₂/C₁₄ 的 N⁺ 区最正；烷基链接近中性 |

### 1.4 预期结果（方案预设，待 VASP 验证）

> 下表的 −0.55 / −0.25 / −0.15 eV 是前期 DMol3/COSMO 口径的参考值（C₁₀–C₁₄ 为本项目新体系，EMIm⁺ 等对照值参考三篇文献）。VASP 的绝对数值预计偏移 ±0.1–0.2 eV，**以排序和差值为准，不以绝对值为准**。

| 计算 | 参考值（DMol3 口径） | VASP 预期区间 | 说明 |
|---|---|---|---|
| Eb(C₁₂TAB⁺–Br₃⁻) | −0.55 eV | −0.5 ~ −0.75 eV | 强于 EMImBr⁺ |
| Eb(EMImBr⁺–Br₃⁻) | −0.25 eV | −0.2 ~ −0.4 eV | 对照下限 |
| Eb(C₁₀TAB⁺–Br₃⁻) | — | 比 C₁₂ 弱 0.05–0.15 eV | 链短、范德华弱 |
| Eb(C₁₄TAB⁺–Br₃⁻) | — | 与 C₁₂ 相当或略强 | 链长、范德华强 |
| Eb(C₁₂TAB⁺–H₂O) | −0.15 eV | −0.1 ~ −0.25 eV | 远弱于 Br₃⁻ |
| Eb(C₁₄TAB⁺–H₂O) | — | 三者中最弱 | 疏水性最强 |
| ESP（C₁₂/C₁₄ 的 N⁺ 区） | 深蓝、强正电势 | V_S,max 最高组 | 静电锚定能力最强 |

### 1.5 与原 DMol3 方案的差异（为什么必须改）

| 项目 | 原方案（DMol3） | 本方案（VASP） | 原因 |
|---|---|---|---|
| 基组 | DNP+ 数值基 | PAW + 平面波 500 eV | 无 BSSE、可扩展到表界面；但绝对能量不可直接对比 |
| 色散 | Grimme 2006（D2） | D3(BJ) 为主，D2 对照 | D3(BJ) 对极性/带电体系更稳健；D2 用于与文献对表 |
| 溶剂 | COSMO（ε=78.54） | VASPsol（ε=78.54，线性化 PB） | 物理同类、数值不同；需重新标定 |
| 电荷 | Mulliken/Hirshfeld | Bader / DDEC6 | VASP 无原生 Mulliken；Hirshfeld 需外挂 |
| 带电胞 | 分子轨道法无背景电荷 | 周期胞有均匀背景 | 必须做盒子收敛 + 交换能判据（§5.2） |
| 热化学 | 未做 / 经验 | 全 Hessian + 部分 Hessian | VASP 频率贵，需分级策略 |
| 动力学 | AIMD（文献有） | AIMD 298 K / 1 fs / 10 ps | 保留并量化限域指标 |

### 1.6 如何解释"C₁₂ 最优"而不是 C₁₄（必读）

若计算显示 C₁₄ 结合最强（很可能），需要三段论收口：

1. **结合强度不是唯一指标**。C₁₄ 结合过强 → Br₃⁻ 被"锁死"，界面电荷转移与 Br₂ 释放动力学变差；AIMD 中表现为 Br₃⁻ 位移过小（<6 Å）、停留时间过长。
2. **疏水性与润湿性的权衡**。C₁₄ 水结合最弱、表面能最低，电解液润湿差、离子电导下降；C₁₂ 处于"够疏水 + 可润湿"的窗口。
3. **二维判据而非单一标量**。以 |Eb(Br₃⁻)| 为横轴、|Eb(H₂O)| 为纵轴做象限图：理想 BCA 应落在"结合强（右）+ 水结合弱（下）"的右下象限。C₁₄ 很可能落在"最强结合但水合最弱"的极端，C₁₂ 落在"强结合 + 适度水合"的平衡点；再以气泡大小叠加 AIMD 位移，即可在一张图里同时表达结合、疏水、动力学三个维度。若一定要给标量，用 Φ = |Eb(Br₃⁻)| − |Eb(H₂O)| − λ·d（d 为 AIMD 中 Br₃⁻ 的质心位移，λ 由实验标定，单位 eV/Å），但论文正文建议用象限图，直观且不必人为定 λ。

---

## 2. 计算体系与模型构建

### 2.1 阳离子清单

| 缩写 | 名称 | 分子式 | SMILES（仅用于生成初猜） | 电荷 | 原子数（含 H） |
|---|---|---|---|---|---|
| C₁₀TAB⁺ | 十烷基三甲基铵 | C₁₃H₃₀N⁺ | CCCCCCCCCC[N+](C)(C)C | +1 | 44 |
| C₁₂TAB⁺ | 十二烷基三甲基铵 | C₁₅H₃₄N⁺ | CCCCCCCCCCCC[N+](C)(C)C | +1 | 50 |
| C₁₄TAB⁺ | 十四烷基三甲基铵 | C₁₇H₃₈N⁺ | CCCCCCCCCCCCCC[N+](C)(C)C | +1 | 56 |
| MEPBr⁺ | 1-乙基-1-甲基吡咯烷鎓 | C₇H₁₆N⁺ | CC[N+]1(CCCC1)C | +1 | 24 |
| TPABr⁺ | 四丙基铵 | C₁₂H₂₈N⁺ | CCC[N+](CCC)(CCC)CCC | +1 | 41 |
| EMImBr⁺ | 1-乙基-3-甲基咪唑鎓 | C₆H₁₁N₂⁺ | CCn1ccn(C)c1 | +1 | 18 |

### 2.2 阴离子 / 中性分子

| 物种 | 初始构型 | 关键内坐标 | 电荷 |
|---|---|---|---|
| Br⁻ | 单原子 | — | −1 |
| Br₃⁻ | 线性对称 | Br–Br ≈ 2.55 Å（优化后 2.5–2.7 Å） | −1 |
| Br₅⁻ | V 型（Br₂···Br₃⁻） | 端基 Br–Br 2.4 Å，中间 Br···Br 3.1 Å | −1 |
| Br₂ | 线性 | Br–Br = 2.28 Å | 0 |
| H₂O | 弯曲 | O–H 0.97 Å，∠HOH 104.5° | 0 |

> **体系边界**：只做纯溴体系，不含 Cl⁻、不涉及 Br⁺、不涉及四电子过程（与 Zn–I 四电子体系区分）。Br₅⁻ 在文献中常以 Br₂·Br₃⁻ 形式存在，优化后若解离为 Br₂ + Br₃⁻，需在结果中如实说明并按"弱结合复合物"处理。

### 2.3 复合物矩阵（主矩阵 15 个 + 扩展）

| 阳离子 \ 阴离子 | Br⁻ | Br₃⁻ | Br₅⁻ | Br₂ | H₂O |
|---|---|---|---|---|---|
| C₁₀TAB⁺ | ✓ | ✓ | ✓ | ✓ | ✓ |
| C₁₂TAB⁺ | ✓ | ✓✓ | ✓ | ✓ | ✓✓ |
| C₁₄TAB⁺ | ✓ | ✓ | ✓ | ✓ | ✓ |
| MEPBr⁺ | — | ✓ | — | — | — |
| TPABr⁺ | — | ✓ | — | — | — |
| EMImBr⁺ | — | ✓✓ | — | — | ✓ |

- ✓✓ = 重点体系，需做全 Hessian / AIMD / 显式水。
- 每个复合物至少 3 个初始构型：N⁺ 正对 Br₃⁻ 端基 Br（3.0 / 3.4 / 3.8 Å 三种距离），取优化后能量最低者。
- 水复合物两种构型：①水 O 朝向 N⁺；②水 H 朝向 Br₃⁻（若有阴离子共存）。两者都要算，取最低。
- 备选 2:1 复合物（[C₁₂TAB]₂(Br₃)₂ 或 C₁₂TAB⁺–(Br₃⁻)₂）用于讨论高载量下的配位模式，非必做。

### 2.4 建模流程

1. **生成初猜**：SMILES → Open Babel（obabel -:"SMILES" --gen3d -O init.xyz）→ RDKit ETKDG 生成 ≥20 个构象 → MMFF94 排序 → 取前 3 个。
2. **阴离子**：按 §2.2 内坐标手工搭建（线性/V 型），比 SMILES 可靠。
3. **复合物**：用 Python（ASE / pymatgen）把阳离子与阴离子按指定距离和取向拼装，N⁺ 指向 Br 端基。
4. **居中**：所有结构置于立方盒中心，分子质心在 (0.5, 0.5, 0.5)。
5. **命名**：<cation>_<anion>_<conf>，如 C12TAB_Br3_c1；目录见 §9.3。

### 2.5 盒子与边界条件

| 体系 | 推荐盒子 | 真空层 | 备注 |
|---|---|---|---|
| 阴离子 / H₂O / Br₂ | 25 × 25 × 25 Å | ≥ 10 Å | 仅用于单体能量 |
| C₁₀ / C₁₂ 及复合物 | 40 × 40 × 40 Å | ≥ 10 Å | 基准 |
| C₁₄ 及复合物 | 44 × 44 × 44 Å | ≥ 9 Å | 分子长约 20 Å，必须加大 |
| AIMD（全部） | 40–44 Å | — | 与静态一致 |

> **一致性优先**：同一张结合能表里的所有物种必须用**同一盒子、同一截断、同一溶剂设置**。若 C₁₄ 用 44 Å、C₁₂ 用 40 Å，两者 Eb 之差会混入 0.01–0.02 eV 的盒子误差，虽然不大，但审稿人会问。建议统一 44 Å，机时增加约 40%，可接受。

---

## 3. 计算参数（VASP）

### 3.1 参数总表

| 参数 | 取值 | 说明 |
|---|---|---|
| 代码 | VASP 6.x | 需确认编译了 VASPsol（LSOL 标签可读） |
| 泛函 | PBE（GGA） | 与三篇文献一致 |
| 赝势 | PAW，标准 PBE（Br、N、C、O、H） | 标量相对论已含；无需全电子 |
| ENCUT | 500 eV | 测试 400 / 500 / 600，ΔE < 5 meV/atom |
| EDIFF | 1×10⁻⁵ eV（几何）；1×10⁻⁶（频率/单点） | 电子步收敛 |
| EDIFFG | −0.02 eV/Å（几何）；−0.005（频率） | 力收敛；比文献 2×10⁻³ Ha/Å（≈0.10 eV/Å）更严 |
| IBRION / ISIF | 2 / 2 | 固定晶胞，只弛豫离子 |
| NSW | 300 | 长链柔性大，留足步数 |
| ISYM | 0 | 分子体系关对称 |
| ISPIN | 1 | 全部闭壳层；若出现 Br 自由基再开 2 |
| ISMEAR / SIGMA | 0 / 0.05 eV | 分子用 Gaussian 展宽 |
| LREAL | .FALSE. | 小体系精确 |
| ALGO | Normal（静态）/ Fast（AIMD） | |
| NCORE | 8 | 按机器调整 |
| IVDW | 12（D3-BJ）；对照 1（D2） | 见 §3.2 |
| 溶剂 | VASPsol，EB_K=78.54 | 见 §3.3 |
| k 点 | Γ 点（1×1×1） | 44 Å 盒子，Γ 已足够；测试 2×2×2 验证 |
| 偶极修正 | 不做（IDIPOL=0） | 0D 周期体系不需要 |

### 3.2 色散校正：D3(BJ) 为主，D2 对照

- **主方案 IVDW=12（D3-BJ）**：对带电、极性、氢键体系更稳健，长链烷基–多溴的范德华描述更好。
- **对照 IVDW=1（D2）**：三篇文献用的是 Grimme 2006（即 D2）。为与文献数值对话，对 5 个关键体系（C₁₂–Br₃⁻、C₁₀–Br₃⁻、C₁₄–Br₃⁻、EMIm–Br₃⁻、TPA–Br₃⁻）各补一套 D2 单点。
- **判据**：两套色散的 Eb 排序必须一致；若排序翻转，说明结论对色散方案敏感，需在文中说明并选 D3(BJ) 为主。

### 3.3 隐式溶剂：VASPsol

~~~
LSOL       = .TRUE.
EB_K       = 78.54      ! 水的介电常数
LAMBDA_D_k = 1.0E4      ! 纯水：设大以关闭离子屏蔽；模拟 1 m 电解液可取 3.0
SIGMA_K    = 0.6        ! 空腔宽度
TAU        = 0.0        ! 表面张力项
NC_K       = 1.0E-4     ! 空腔截断
~~~

- VASPsol 解线性化 Poisson–Boltzmann，物理上与 COSMO 同类（连续介质 + 空腔），但**数值不同**，不能指望复现 −0.55 eV 这种具体数字。
- **纯水 vs 电解液**：实验电解液是 1 m ZnSO₄ + 1 m ZnBr₂，不是纯水。生产计算用纯水（λ_D 大），另做一组 λ_D = 3 Å 的敏感性测试，说明"离子屏蔽不改变排序"。
- **备选**：Environ（SCCS 模型）或"真空 + 显式水簇 + 经验校正"。若 VASPsol 编译有问题，先落到"真空 + 3–5 个显式水"的方案，结论仍可成立（但绝对值不可与文献直接比）。

### 3.4 带电周期胞的处理（关键）

阳离子 +1、阴离子 −1，VASP 自动加均匀背景电荷。背景带来两项误差：

1. **单胞自能**：E_corr ≈ q²·α_M / (2εL)，α_M = 2.8373（立方），L 为盒长。
   - L = 44 Å、ε = 78.54：E_corr ≈ 0.006 eV，可忽略。
   - L = 44 Å、ε = 1（真空）：E_corr ≈ 0.5 eV，**不可忽略**。
2. **Eb 中的叠加**：Eb 减去两个带电单体，误差约为 2E_corr。溶剂化下 ≈ 0.01 eV，真空下 ≈ 1 eV。

**应对策略（按优先级）**：

- **策略 A（推荐）**：所有生产计算都开隐式溶剂，带电误差被 ε=78.5 压低两个数量级；再用"交换能"（§5.2）做最终判据，交换能两侧带电胞数目相同，误差精确抵消。
- **策略 B**：真空计算只用于"本征相互作用"趋势，报告时注明未做背景修正，或按 Makov–Payne 加回 q²α/(2L)。
- **策略 C**：盒子收敛测试（36 / 40 / 44 Å），要求 Eb 漂移 < 0.03 eV。

### 3.5 与三篇文献参数的对照

| 项目 | 双咪唑 BCA（Angew 2025） | 疏水有机溴 M0–M6（Small Methods 2025） | 疏水季铵碘 C₈–C₁₈I₃（Angew 2026） | 本方案（VASP） |
|---|---|---|---|---|
| 代码 | DMol3 | DMol3 | DMol3（Materials Studio 2020） | VASP 6.x |
| 泛函 | PBE | PBE | PBE | PBE |
| 基组 | DNP+ | DNP+ | DNP | PAW + 500 eV |
| 色散 | Grimme（D2） | Grimme（D2） | **未说明**（疑未加） | D3(BJ)，D2 对照 |
| 溶剂 | COSMO，ε=78.54 | COSMO，ε=78.54 | COSMO，ε=78.54 | VASPsol，ε=78.54 |
| 轨道截断 | 5.0 Å | 5.0 Å | — | 平面波无对应 |
| 能量收敛 | 1×10⁻⁵ Ha | 1×10⁻⁵ Ha | 1×10⁻⁵ Ha（SCF 1×10⁻⁶） | 1×10⁻⁵ eV |
| 力收敛 | 2×10⁻³ Ha/Å（≈0.10 eV/Å） | 同左 | 2×10⁻³ Ha/Å | 0.02 eV/Å（更严） |
| 位移收敛 | 5×10⁻³ Å | 同左 | — | 以力为准 |
| 动力学 | AIMD 298 K，1 fs，10 ps | — | — | AIMD 298 K，1 fs，10 ps |

> **注意**：三篇文献的色散处理并不统一（第三篇未写 Grimme），所以它们之间、它们与本方案之间的绝对 Eb 都不可直接比较。**只有同一套参数内部算出的排序才有意义**，这是本方案做方法标定的原因。

### 3.6 收敛性与方法标定清单（T0 任务）

| 测试项 | 变量 | 体系 | 通过标准 |
|---|---|---|---|
| 截断能 | 400 / 500 / 600 eV | C₁₂–Br₃⁻ | ΔEb < 0.03 eV |
| 盒子 | 36 / 40 / 44 Å | C₁₂–Br₃⁻ | ΔEb < 0.03 eV |
| k 点 | Γ vs 2×2×2 | C₁₂–Br₃⁻ | ΔE < 1 meV |
| 展宽 | SIGMA 0.02 / 0.05 / 0.10 | C₁₂–Br₃⁻ | ΔE < 5 meV |
| 色散 | D2 vs D3(BJ) | 5 个关键体系 | 排序一致 |
| 溶剂 | 真空 / VASPsol | C₁₂–Br₃⁻、EMIm–Br₃⁻ | 报告两套，排序一致 |
| 背景电荷 | 加/不加 Makov–Payne | C₁₂–Br₃⁻ | ΔEb < 0.02 eV（溶剂化下） |

**方法标定**：用同一套 VASP 参数重算 EMIm⁺–Br₃⁻、TPA⁺–Br₃⁻，与文献/前期的 −0.25 eV 量级对比，记录系统偏差 δ = E_VASP − E_DMol3。若 δ 稳定（例如 −0.15 eV），可在文中说明"绝对值整体偏移 δ，排序不变"，**不要**做逐体系的经验平移。

---

## 4. 计算任务分解

### T0　方法标定与收敛（3–5 天）

- **目的**：确定一套可复现的参数，回答"VASP 的数字和 DMol3 差多少、排序是否一致"。
- **体系**：C₁₂–Br₃⁻、EMIm–Br₃⁻。
- **步骤**：按 §3.6 逐项跑；每个测试只做单点或短优化。
- **输出**：收敛性表格、δ 值、最终参数集（写入 INCAR.template）。
- **判据**：六项测试全部通过；δ 的符号和量级可解释（平面波基 + D3BJ + VASPsol vs 数值基 + D2 + COSMO）。
- **机时**：约 300–500 核时。

### T1　几何优化与真空结合能（1 周）

- **目的**：拿到所有单体和复合物的可靠基态结构，给出本征结合能排序。
- **体系**：§2.3 主矩阵 15 个 + 4 个对照 = 19 个复合物，加 11 个单体。
- **步骤**：
  1. 单体优化（阳离子、Br⁻、Br₃⁻、Br₅⁻、Br₂、H₂O）。
  2. 复合物优化：每个 3 个初始构型 → 取最低。
  3. 单点（IBRION=-1，LVTOT=.TRUE.，LAECHG=.TRUE.）拿 LOCPOT/CHGCAR。
  4. 计算 Eb = E(complex) − E(cation) − E(anion)。
- **输出**：Eb_vac.csv、优化结构（CONTCAR）、能量汇总表。
- **判据**：
  - C₁₂–Br₃⁻ 比 EMIm–Br₃⁻ 更负 ≥ 0.15 eV；
  - Br₃⁻/Br₅⁻ 比 Br₂ 更负 ≥ 0.15 eV；
  - C₁₀ < C₁₂ ≤ C₁₄（结合强度递增）。
- **注意**：Br⁻ 可能因电荷密度高而结合更强，这是正常的；**不要**用 Br⁻ 的数据否定 Br₃⁻ 的选择性，选择性要看 T2 的溶剂化竞争和 §5.2 的交换能。
- **机时**：约 1,500–3,000 核时。

### T2　溶剂化结合能与水竞争（1 周）

- **目的**：回答"阳离子更愿意抓多溴还是抓水"。
- **体系**：C₁₀/C₁₂/C₁₄–Br₃⁻、C₁₀/C₁₂/C₁₄–H₂O、EMIm–Br₃⁻、EMIm–H₂O、MEP/TPA–Br₃⁻。
- **步骤**：
  1. 全部加 VASPsol 重新优化（溶剂会改变构象，不能拿真空结构直接单点）。
  2. 计算 Eb_solv。
  3. 水复合物两种构型都算，取最低。
  4. 显式水微溶剂化：C₁₂–Br₃⁻–(H₂O)ₙ，n = 0–4，看第一溶剂化壳对 Eb 的影响。
- **输出**：Eb_solv.csv、水竞争柱状图、显式水能量曲线。
- **判据**：|Eb(Br₃⁻)| − |Eb(H₂O)| > 0.2 eV；C₁₄–H₂O 最弱。
- **机时**：约 1,500–2,500 核时。

### T3　静电势与原子电荷（3–5 天）

- **目的**：给出"正电荷在 N⁺、烷基链疏水"的图像，并量化。
- **体系**：C₁₀/C₁₂/C₁₄TAB⁺、MEP⁺、TPA⁺、EMIm⁺；可选 C₁₂–Br₃⁻ 复合物。
- **步骤**：
  1. 单点 + LVTOT=.TRUE. 得 LOCPOT；LAECHG=.TRUE. 得 AECCAR0/2。
  2. LOCPOT → cube（VASPKIT / vasp2cube）。
  3. Multiwfn 主功能 12：在 0.001 e/Å³ 等值面上统计 V_S,max、V_S,min、平均、方差、正电势面积占比。
  4. Bader：chgsum.pl AECCAR0 AECCAR2 → bader CHGCAR -ref CHGCAR_sum；可选 DDEC6（Chargemol）。
  5. 出图：VESTA 或 Multiwfn 导出，统一色标、统一等值面、统一取向（N⁺ 朝上）。
- **输出**：ESP 图（三阳离子并排）、esp_stats.csv、bader.csv。
- **判据**：C₁₂/C₁₄ 的 V_S,max 高于 C₁₀；烷基链区 V_S ≈ 0；N⁺ 的 Bader 电荷与 ESP 趋势一致。
- **注意**：
  - VASP 的 LOCPOT 含 XC 势，严格说不是纯静电势；在 0.001 e/Å³ 面外 XC 已衰减，定性/半定量可用。若要严格 ESP，用 Multiwfn 从密度 cube 重算。
  - 带电体系的绝对电势有任意常数，**只比较同一电荷态、同一盒子下的相对值**；报告 V_S,max − V_S,min 和排序，不报绝对值。
  - Bader 电荷在 C₁₀/C₁₂/C₁₄ 之间的差异可能只有 0.01–0.02 e，**别过度解读**；主证据用 ESP，Bader 作辅证。
- **机时**：约 300–600 核时。

### T4　热化学与 ΔG（1 周）

- **目的**：把电子能量升级为自由能，回答"络合是否自发"。
- **体系**：
  - 全 Hessian（≤25 原子）：Br₃⁻、Br₂、H₂O、EMIm⁺、MEP⁺、TPA⁺、C₁₀TAB⁺（验证用）。
  - 部分 Hessian（头基 ~20 原子）：C₁₂TAB⁺、C₁₄TAB⁺ 及其 Br₃⁻/H₂O 复合物。
- **步骤**：
  1. IBRION=5, NFREE=2, POTIM=0.015, EDIFF=1E-6, ADDGRID=.TRUE.。
  2. 提取频率 → ZPE、298.15 K 下的振动熵/平动转动熵。
  3. ΔG = ΔE + ΔZPE − TΔS；VASPsol 的溶剂化自由能已含在 ΔE 中。
  4. 虚频检查：出现虚频则微调结构重算。
- **输出**：thermo.csv、ΔG 表、ΔG 台阶图。
- **判据**：ΔG(C₁₂–Br₃⁻) < 0；ΔG(C₁₂–H₂O) 更正（更不利）；部分 Hessian 与全 Hessian 在 C₁₀ 上的偏差 < 0.03 eV。
- **注意**：长链的甲基/亚甲基转动是低频模式，谐振近似会高估熵；建议对 < 50 cm⁻¹ 的模式按自由转子处理，或直接引用"ΔG 的绝对值不确定度 ±0.1 eV，排序可信"。
- **机时**：约 800–1,500 核时（频率是机时大户，务必先做 C₁₀ 的收敛测试）。

### T5　AIMD 限域动力学（1–2 周）

- **目的**：证明"多溴被限制在阳离子周围 ~10 Å，不会跑掉"，并与 EMIm/TPA 对比。
- **体系**：C₁₀/C₁₂/C₁₄–Br₃⁻、EMIm–Br₃⁻、TPA–Br₃⁻（5 个）。
- **参数**：NVT，298 K（筛选可用 330 K 加速），Nosé–Hoover（MDALGO=2, SMASS=0），1 fs，10 ps（筛选 5 ps），Γ 点，EDIFF=1E-4。
- **步骤**：
  1. 从 T2 的溶剂化优化结构出发。
  2. 先 1 ps 平衡（丢弃），再 9 ps 采样。
  3. 分析（见 §5.6）。
- **输出**：aimd_stats.csv、MSD 曲线、RDF、代表性 snapshots、势能涨落曲线。
- **判据**：
  - C₁₂ 的 Br₃⁻ 质心最大位移 8–12 Å；
  - EMIm > 15 Å，TPA < 6 Å；
  - C₁₀ > C₁₂ > C₁₄（链越长越受限）；
  - Br–Br 内坐标 RDF 第一峰不消失（Br₃⁻ 未解离）。
- **注意**：
  - 隐式溶剂 + MD 可能数值不稳；若不稳，改用真空 AIMD 做"位阻对比"（结论只谈相对限域，不谈溶剂化）。
  - 显式水 AIMD（C₁₂–Br₃⁻ + 40 H₂O，10 ps）作为加分项，机时约 3–5 倍。
- **机时**：约 8,000–15,000 核时（最大头，先跑 5 ps 筛选）。

### T6　可选加分项（按机时决定）

| 项目 | 内容 | 价值 |
|---|---|---|
| 界面模型 | 石墨（0001）6×6 超胞 + C₁₂–Br₃⁻ 复合物，算吸附能 | 支撑"界面富集/疏水屏障" |
| Zn 负极 | Zn(0001) + 阳离子吸附，算吸附能 | 解释负极保护 |
| Br₂ 释放 | 复合物 → 阳离子 + Br₂ 的解离能垒（CI-NEB） | 解释可逆性/动力学 |
| 显式电解液 | C₁₂–Br₃⁻ + Zn²⁺/SO₄²⁻ + 水 | 更接近真实电解液，但极贵 |

---

## 5. 数据处理与判据

### 5.1 结合能

**Eb = E(complex) − E(cation) − E(anion)**

- 符号约定：**负值 = 放热 = 结合有利**。
- 所有项必须同盒子、同截断、同溶剂。
- 2:1 复合物按化学计量比扩展：Eb = E(complex) − 2E(cation) − E(anion)。

### 5.2 交换能（推荐作为选择性判据）

直接比较两个带电单体的绝对能量，会引入背景电荷误差。用**电荷守恒的交换反应**可以精确抵消：

**Br₃⁻/Br⁻ 选择性**

ΔE_ex = [E(C₁₂–Br₃⁻) + E(Br⁻)] − [E(C₁₂–Br⁻) + E(Br₃⁻)]

（相当于反应 [C₁₂TAB⁺–Br⁻] + Br₃⁻ → [C₁₂TAB⁺–Br₃⁻] + Br⁻ 的反应能）

- 两侧各有一个 −1 带电胞 + 一个中性复合物胞，背景误差同号同大小，相减抵消。
- **ΔE_ex < 0 → 优先结合 Br₃⁻**，这正是"选择性锚定多溴"的热力学证据。
- 代数上 ΔE_ex ≡ Eb(Br₃⁻) − Eb(Br⁻)，所以它同时也是两种阴离子结合能之差；用交换式写是为了让背景误差显式抵消。

**水竞争**（无电荷守恒的交换式，改用直接比较）

ΔE_H2O = Eb(Br₃⁻) − Eb(H₂O)

- 溶剂化下带电胞误差仅 ~0.01 eV，可直接比。
- 判据：ΔE_H2O < −0.2 eV → 明显偏好 Br₃⁻。

**一致性检验**：交换式与直接差减在代数上恒等，但数值上只有在盒子/截断/溶剂完全一致时才严格相等。若两者差 > 0.02 eV，说明某项设置不统一，应回查。

### 5.3 自由能

**ΔG = ΔE_elec + ΔZPE − TΔS**

- T = 298.15 K。
- VASPsol 的溶剂化自由能已含在 E_elec 中；它主要是静电（+空腔）项，不含溶质–溶剂色散，D3 只补溶质内部色散，这部分系统性偏差计入 §5.7 的误差预算。
- 谐振近似 + 理想气体平动转动；低频模式（<50 cm⁻¹）按自由转子处理。
- 不确定度：±0.1 eV；**只报告 > 0.15 eV 的差异为显著**。

### 5.4 ESP 定量指标

| 指标 | 定义 | 用途 |
|---|---|---|
| V_S,max | 表面静电势最大值（eV） | N⁺ 区正电势强度 |
| V_S,min | 表面静电势最小值 | 烷基链/阴离子区负电势 |
| σ² | 表面电势方差 | 电荷分离程度 |
| A⁺ | V_S > +1.0 eV 的面积占比 | "正电势区"大小 |
| V_S(N⁺) | N⁺ 原子外 1.2 Å 处的电势 | 更局域的锚定位点强度 |

- 等值面统一用 0.001 e/Å³；色标统一 −1.0 ~ +3.0 eV；三张图同一取向、同一比例。
- 报告"排序"而非绝对值（带电体系绝对电势有任意常数）。

### 5.5 原子电荷

- **Bader（AIM）**：主用。chgsum.pl AECCAR0 AECCAR2 → bader CHGCAR -ref CHGCAR_sum。
- **DDEC6**（Chargemol）：可选，更化学直观。
- **不要用 Mulliken**：平面波基下无定义；Hirshfeld 需外挂，且对阴离子偏差大。
- 报告 N、Br、C(α) 的电荷；预期 N 的电荷差异 < 0.03 e，**只作辅证**。

### 5.6 AIMD 指标

| 指标 | 算法 | 物理含义 |
|---|---|---|
| 势能涨落 | std(E_pot) over 9 ps | 结构刚性/柔性 |
| Br₃⁻ 质心位移 | r_CM(t) 相对初值的最大位移（去质心漂移） | 限域程度 |
| N⁺···Br 距离分布 | g(r) 第一峰位置/高度 | 锚定强度 |
| Br–Br RDF | 分子内第一峰是否保留 | Br₃⁻ 是否解离 |
| 停留时间 τ | 距离 < 5 Å 的连续时间自相关 | 交换动力学 |
| 扩散系数 D | MSD 斜率 / 6 | 迁移率 |
| 水配位数 | 第一壳层水分子数（显式水） | 疏水屏障 |

### 5.7 误差预算

| 来源 | 量级 | 处理 |
|---|---|---|
| 泛函（PBE） | 0.1–0.2 eV | 只比排序；必要时 PBE0 单点验证 2 个体系 |
| 色散 | 0.05–0.1 eV | D2/D3 双算 |
| 溶剂模型 | 0.1–0.2 eV | 真空/隐式/显式水三套对照 |
| 背景电荷 | 0.01 eV（溶剂化） | 交换能 + 盒子收敛 |
| 热化学 | 0.1 eV | 部分 Hessian 验证 |
| 构象 | 0.05–0.15 eV | 多构象 Boltzmann 平均 |
| **总计** | **±0.2 eV** | 结论只建立在 > 0.2 eV 的差异上 |

---

## 6. 结果表模板

### 6.1 结合能（eV）

| 阳离子 | Br⁻ | Br₃⁻ | Br₅⁻ | Br₂ | H₂O |
|---|---|---|---|---|---|
| C₁₀TAB⁺ | | | | | |
| C₁₂TAB⁺ | | | | | |
| C₁₄TAB⁺ | | | | | |
| MEPBr⁺ | — | | — | — | — |
| TPABr⁺ | — | | — | — | — |
| EMImBr⁺ | — | | — | — | |

> 每格填 Eb(vac) / Eb(solv) 两个数。

### 6.2 交换能与选择性（eV）

| 体系 | ΔE_ex(Br₃⁻/Br⁻) | ΔE_H2O = Eb(Br₃⁻) − Eb(H₂O) | 结论 |
|---|---|---|---|
| C₁₀TAB⁺ | | | |
| C₁₂TAB⁺ | | | |
| C₁₄TAB⁺ | | | |
| EMImBr⁺ | | | |

### 6.3 ESP 与电荷

| 体系 | V_S,max (eV) | V_S,min (eV) | σ² | A⁺ (%) | q_Bader(N) |
|---|---|---|---|---|---|
| C₁₀TAB⁺ | | | | | |
| C₁₂TAB⁺ | | | | | |
| C₁₄TAB⁺ | | | | | |
| MEPBr⁺ | | | | | |
| TPABr⁺ | | | | | |
| EMImBr⁺ | | | | | |

### 6.4 热化学与动力学

| 体系 | ΔG(Br₃⁻) (eV) | ΔG(H₂O) (eV) | std(E_pot) (eV) | 最大位移 (Å) | τ (ps) | D (cm²/s) |
|---|---|---|---|---|---|---|
| C₁₀TAB⁺ | | | | | | |
| C₁₂TAB⁺ | | | | | | |
| C₁₄TAB⁺ | | | | | | |
| EMImBr⁺ | | | | | | |
| TPABr⁺ | | | | | | |

---

## 7. 论文论点与计算图的对应

| 论文论点 | 支撑计算 | 建议图 | 关键判据 |
|---|---|---|---|
| C₁₂TAB⁺ 能牢牢抓住多溴 | Eb、ΔG | 柱状图（6 个阳离子 × 5 个阴离子） | C₁₂ 比 EMIm 更负 ≥ 0.15 eV |
| 优先结合多溴而非水 | Eb(H₂O)、ΔE_H2O | 双柱对比 + 差值标注 | 差值 > 0.2 eV |
| N⁺ 是静电锚定位点 | ESP、Bader | ESP 三图并排 + V_S,max 标注 | N⁺ 区最正 |
| 长链提供疏水屏障 | ESP 链区、C₁₄–H₂O | ESP 侧视图 + 水结合能趋势 | 链区 ≈ 0；C₁₄–H₂O 最弱 |
| 多溴被限域、不穿梭 | AIMD | MSD + RDF + snapshots | C₁₂ 位移 8–12 Å |
| C₁₂ 综合最优 | 象限图 + Φ | Eb(Br₃⁻)–Eb(H₂O) 象限图（气泡 = AIMD 位移） | 结合–疏水–动力学平衡 |

---

## 8. 风险与备选方案

| 风险 | 触发条件 | 备选方案 |
|---|---|---|
| VASPsol 编译失败 | LSOL 报 unknown tag | 改用 Environ；或真空 + 显式水簇；或只做真空趋势 |
| 带电胞不收敛 | SCF 振荡、能量随盒子漂移 | 加大盒子、加大 SIGMA、先用 ALGO=All 预收敛 |
| 隐式溶剂下结构跑飞 | 优化中分子散架 | 先真空优化 → 再溶剂化优化；减小 POTIM |
| C₁₄ 结合最强、与 C₁₂ 最优矛盾 | T1 结果 | 用 §1.6 的象限图/Φ 判据；补 AIMD 动力学惩罚 |
| 频率计算太贵 | C₁₂/C₁₄ 全 Hessian 超预算 | 部分 Hessian；或只对 Br₃⁻/H₂O 做全 Hessian，阳离子用刚性近似 |
| AIMD 太贵 | 10 ps × 5 体系超机时 | 5 ps 筛选 + 只对 C₁₂ 做 10 ps；提高温度到 330 K |
| Bader 电荷差异太小 | Δq < 0.02 e | 以 ESP 为主证据，Bader 只作辅证 |
| 与文献数值差太多 | δ > 0.3 eV | 说明基组/溶剂/色散差异；只比排序；补 D2 对照 |

---

## 9. 进度、机时与交付

### 9.1 六周计划

| 周 | 任务 | 里程碑 |
|---|---|---|
| W1 | T0 收敛 + 建模 + 单体优化 | 参数定稿，19 个复合物初构 |
| W2 | T1 几何优化 + 真空 Eb | Eb 表初稿 |
| W3 | T2 溶剂化 + 水竞争 | 选择性结论 |
| W4 | T3 ESP/电荷 + T4 热化学 | ESP 图、ΔG 表 |
| W5 | T5 AIMD（5 ps 筛选） | MSD/RDF 初稿 |
| W6 | 补算、出图、写作 | 图表终稿 + 数据归档 |

### 9.2 机时估算（128 核，仅供参考）

| 任务 | 体系数 | 单体系 | 合计 |
|---|---|---|---|
| 单体 + 复合物优化 | ~30 | 2–4 h | 60–120 h |
| 溶剂化重优化 | ~15 | 2–4 h | 30–60 h |
| 频率（全/部分 Hessian） | 8 | 10–30 h | 100–200 h |
| AIMD 5 ps | 5 | 25–40 h | 125–200 h |
| AIMD 10 ps（C₁₂） | 1 | 50–80 h | 50–80 h |
| **合计** | | | **约 400–650 h（≈ 5–8 万核时）** |

> 先跑 T0 的基准算例，用实测的 SCF/步时间替换上表，再决定 AIMD 体系数。

### 9.3 目录结构

~~~
calc/
├── 00_structures/          # 初始 xyz/POSCAR
├── 01_single/              # 单体
│   ├── C12TAB/
│   ├── Br3/
│   └── ...
├── 02_complex/             # 复合物
│   ├── C12TAB_Br3_c1/
│   └── ...
├── 03_freq/
├── 04_aimd/
├── 05_esp/
├── 06_bader/
├── scripts/
│   ├── run_vasp.sh
│   ├── collect_energy.py
│   └── analysis_aimd.py
└── results/
    ├── eb.csv
    ├── esp_stats.csv
    └── figures/
~~~

### 9.4 交付物

1. 结合能表（真空 + 溶剂化，含 D2/D3 对照）。
2. 交换能/选择性表 + 水竞争图。
3. ESP 图（三阳离子并排，统一色标）+ V_S 统计表。
4. Bader/DDEC6 电荷表。
5. ΔG 表 + 台阶图。
6. AIMD MSD/RDF/snapshots + 指标表。
7. 优化结构（CONTCAR）、INCAR/KPOINTS、后处理脚本、原始 OUTCAR 归档。

---

## 10. 附录

### A. INCAR 模板

**A1 几何优化（含隐式溶剂）**

~~~
SYSTEM = C12TAB-Br3 opt
PREC   = Accurate
ENCUT  = 500
EDIFF  = 1E-5
EDIFFG = -0.02
IBRION = 2
ISIF   = 2
NSW    = 300
ISYM   = 0
ISPIN  = 1
ALGO   = Normal
LREAL  = .FALSE.
ISMEAR = 0
SIGMA  = 0.05
IVDW   = 12
NCORE  = 8
LWAVE  = .FALSE.
LCHARG = .TRUE.
LSOL       = .TRUE.
EB_K       = 78.54
LAMBDA_D_k = 1.0E4
SIGMA_K    = 0.6
TAU        = 0.0
NC_K       = 1.0E-4
~~~

**A2 单点 + LOCPOT + AECCAR**

~~~
IBRION = -1
NSW    = 0
EDIFF  = 1E-6
LVTOT  = .TRUE.
LCHARG = .TRUE.
LAECHG = .TRUE.
~~~

**A3 频率**

~~~
IBRION  = 5
NFREE   = 2
POTIM   = 0.015
EDIFF   = 1E-6
EDIFFG  = -0.005
ISYM    = 0
ADDGRID = .TRUE.
NWRITE  = 0
~~~

**A4 AIMD**

~~~
IBRION = 0
NSW    = 10000
POTIM  = 1.0
TEBEG  = 298
TEEND  = 298
MDALGO = 2
SMASS  = 0
ISYM   = 0
EDIFF  = 1E-4
ALGO   = Fast
NCORE  = 8
NBLOCK = 10
KBLOCK = 100
~~~

### B. KPOINTS

~~~
Gamma only
0
Gamma
1 1 1
0 0 0
~~~

### C. 后处理命令

~~~
# Bader
chgsum.pl AECCAR0 AECCAR2
bader CHGCAR -ref CHGCAR_sum

# LOCPOT/CHGCAR -> cube（任选其一）
vaspkit          # 功能 5xx：格点数据转换
vasp2cube LOCPOT

# ESP 表面统计（Multiwfn 主功能 12）
Multiwfn C12TAB.cube
~~~

### D. 最小可行集（机时紧张时先做这 8 个）

1. C₁₂TAB⁺–Br₃⁻（核心）
2. C₁₂TAB⁺–Br⁻（交换能）
3. C₁₂TAB⁺–H₂O（水竞争）
4. C₁₀TAB⁺–Br₃⁻（链长对照）
5. C₁₄TAB⁺–Br₃⁻（链长对照）
6. EMImBr⁺–Br₃⁻（对照）
7. TPABr⁺–Br₃⁻（对照）
8. C₁₂TAB⁺ 单体（ESP）

这 8 个算完就能给出"结合能排序 + 水竞争 + ESP"三条主结论；AIMD 和 ΔG 作为加强项后补。

### E. 常见坑

1. **盒子不够大**：C₁₄ 链长约 20 Å，25 Å 盒子会让分子周期性镜像自相互作用，Eb 偏负 0.1–0.3 eV。
2. **溶剂化后没重新优化**：拿真空结构直接单点算溶剂化能量，构象不对，误差可达 0.2 eV。
3. **带电单体的背景电荷**：真空下不修正会引入 ~1 eV 误差；务必用溶剂化或交换能。
4. **ESP 参考点不一致**：不同体系用不同的真空参考，图好看但不可比；统一盒子、统一参考。
5. **频率虚频**：长链低频模式多，优化不严容易出虚频；先 EDIFFG=-0.005 收紧再算频率。
6. **AIMD 热浴**：用 TEBEG=298 直接起步，前 1 ps 仍可能温度漂移，平衡段一律丢弃；不要从 0 K 缓慢升温再计入采样。
7. **VASPsol 的 LAMBDA_D_k**：纯水设 1e4，模拟电解液设 3.0；设错会让溶剂化能差 0.5 eV 以上。
8. **D2 与 D3 混用**：同一张表里不能一半 D2 一半 D3；对照要成对出现。

---

## 11. 参考资料

**体系与实验背景**

1. Y. Wu, C. Xu, C. Lei, W. Ma, W. Zhang, X. He, X. Liang. Optimizing Bromine Complexation and Kinetics: A Bisimidazole Strategy for High-Performance Zn–Br Static Batteries. *Angew. Chem. Int. Ed.* **2025**, e202509293.（SI：anie202509293）
2. J. Jin, C. Liu, C. Lei, C. Xu, H. Wang, J. Li, X. Liang. Efficient Aqueous Static Zinc–Bromine Batteries Enabled by Hydrophobic Organic Bromide Design. *Small Methods* **2025**.（SI：smtd70112）
3. X. Bao, J. Shen, W. Yuan, Y. Wang, X. Li, X. Ru, Z. Shen, Y. Wang, L. Ma, N. Zhang. Hydrophobic Quaternary Ammonium Iodides for Highly Stable Four-Electron Zinc–Iodine Batteries. *Angew. Chem. Int. Ed.* **2026**.（SI：anie74000）

**计算方法**

4. G. Kresse, J. Furthmüller. Efficient iterative schemes for ab initio total-energy calculations using a plane-wave basis set. *Phys. Rev. B* **1996**, 54, 11169.
5. J. P. Perdew, K. Burke, M. Ernzerhof. Generalized Gradient Approximation Made Simple. *Phys. Rev. Lett.* **1996**, 77, 3865.
6. S. Grimme. Semiempirical GGA-type density functional constructed with a long-range dispersion correction. *J. Comput. Chem.* **2006**, 27, 1787.（D2）
7. S. Grimme, J. Antony, S. Ehrlich, H. Krieg. A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H–Pu. *J. Chem. Phys.* **2010**, 132, 154104；S. Grimme, S. Ehrlich, L. Goerigk. *J. Comput. Chem.* **2011**, 32, 1456.（D3-BJ）
8. K. Mathew, R. Sundararaman, K. Letchworth-Weaver, T. A. Arias, R. G. Hennig. Implicit solvation model for density-functional study of nanocrystal surfaces and reaction pathways. *J. Chem. Phys.* **2014**, 140, 084106.（VASPsol）
9. O. Andreussi, I. Dabo, N. Marzari. Revised self-consistent continuum solvation in electronic-structure calculations. *J. Chem. Phys.* **2012**, 136, 064102.（Environ/SCCS，备选）
10. W. Tang, E. Sanville, G. Henkelman. A grid-based Bader analysis algorithm without lattice bias. *J. Phys.: Condens. Matter* **2009**, 21, 084204.
11. T. Lu, F. Chen. Multiwfn: A multifunctional wavefunction analyzer. *J. Comput. Chem.* **2012**, 33, 580.
12. K. Momma, F. Izumi. VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data. *J. Appl. Crystallogr.* **2011**, 44, 1272.

---

*本方案基于用户提供的两张方案图与三篇文献（含 SI）整理；具体数值以实际 VASP 计算结果为准。*
