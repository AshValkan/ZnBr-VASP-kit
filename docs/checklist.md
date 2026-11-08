# 计算检查清单

## A. 提交前（一次性）

- [x] VASP 能识别 VASPsol（G0 已验证：OUTCAR 回读 LSOL=T + Ediel_sol 输出行）
- [ ] POTCAR 路径正确（config.json potcar_dir），且每个元素目录存在（C/H/N/O/Br）
- [ ] POTCAR 顺序 = POSCAR 元素顺序（N C H Br O），可用 grep TITEL POTCAR 核对
- [ ] vasp_std 在服务器可用（~/calc/vasp_std，Gamma-only 用 vasp_std 跑）
- [ ] 核块配置与机器一致（4 块 × 5 物理核，MXene NEB 共享机器时 submit.py 自动避开占用的块）
- [ ] 磁盘：68 体系 × 多任务，建议 ≥ 200 GB（WAVECAR 已关，主要占用是 CHGCAR/LOCPOT）

## B. 方法标定（T0，必须先做）

- [ ] C12TAB-Br3 截断能 400/500/600，Eb 漂移 < 0.03 eV
- [ ] 盒子 36/40/44 A，Eb 漂移 < 0.03 eV
- [ ] k 点 Gamma vs 2x2x2，能量差 < 1 meV
- [ ] D2 与 D3(BJ) 排序一致
- [ ] 真空 vs VASPsol 的 Eb 都记录，排序一致
- [ ] 与文献 EMIm/TPA 的 -0.25 eV 量级比对，记录系统偏差 delta

## C. 电荷与自旋（本项目特有）

- [ ] INCAR 的 NELECT = Σ(n×ZVAL) − charge（make_inputs.py 自动写入，抽查 2-3 例）
  - 阳离子单体 +1（C12TAB⁺ 价电子 98）；Br⁻/Br₃⁻/Br₅⁻ 单体 −1（8/22/29）
  - 复合物：Br/Br₃/Br₅ 类 0；**Br₂/H₂O 类 +1**
- [ ] 奇电子体系（Br₅⁻ 29e）INCAR 里 ISPIN=2（make_inputs.py 自动翻转）
- [ ] Eb 差减的电荷态自洽：E(complex) − E(cation) − E(anion) 各项 NELECT 相加一致

## D. 几何优化

- [ ] OUTCAR 末尾有 "reached required accuracy"
- [ ] 最大力 < 0.02 eV/A（EDIFFG）
- [ ] 没有原子飞出盒子（检查 CONTCAR 分数坐标在 0-1）
- [ ] 溶剂化优化确实重新做了（不是真空结构直接单点）
- [ ] 三个初始构象的最终能量差 < 0.1 eV，否则说明存在多个局部极小，需在文中说明取最低

## E. 单点 / Bader / ESP

- [ ] sp 任务 OUTCAR 有 "General timing"
- [ ] LOCPOT 与 CHGCAR 网格维数一致（同一目录、同一任务）
- [ ] AECCAR0 + AECCAR2 求和后再做 Bader
- [ ] ACF.dat 行数 = 原子数
- [ ] 净电荷总和 ≈ 体系总电荷（见 C 节），偏差 > 0.05 e 要查
- [ ] ESP：真空参考取远离分子的区域；不同体系用同一盒子、同一等值面

## F. 频率 / 热化学

- [ ] 虚频数 = 0（线性分子/过渡态除外）
- [ ] 频率数目 = 3N-6（非线性）或 3N-5（线性）
- [ ] ZPE 量级合理（每原子 ~0.1-0.2 eV）
- [ ] ΔG 的绝对值不确定度按 ±0.1 eV 报告，只强调排序

## G. AIMD

- [ ] 温度在前 1 ps 后稳定在 298 ± 20 K
- [ ] 势能无持续漂移（有漂移说明失衡或步长过大）
- [ ] 丢掉平衡段后再统计 MSD/RDF
- [ ] Br-Br 分子内 RDF 第一峰仍在（Br3- 未解离）
