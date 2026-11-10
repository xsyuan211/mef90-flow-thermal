# 施工排量–页岩可压性匹配：研究基础、代码接入可行性与技术路线

**日期**：2026-11-10
**基线版本**：`HEAD=027fa86`（CI #121 收口，闸门 4 = 14/14、24 项已落地、0 项在册待做）；**B0 落地版**见第 11 节（基线 `3a819ab`）
**依据**：公开文献与软件发布信息（2026-11-10 检索）、项目源码静态阅读、CI #111–#121 真实 diagnostics artifact、纯 `numpy` 离线标度计算
**约束**：零本地编译、零本地求解器运行。本文所有数值均为离线标度/网格元数据读取，不是求解器结果。
**状态**：阶段 A 研究设计稿 + **阶段 B0 已落地接线、待 CI 验收**（新增 2 个 deck + V1–V8 判据 + CI 硬闸门，仍未改任何 Fortran）；引用线索见第 9 节，投稿前需二次核对卷期与 DOI。

---

## 0. 结论先行

| 问题 | 结论 |
|---|---|
| 现有研究是否已回答"排量–可压性匹配"？ | **没有形成统一结论**。文献对"Q↑ 是否增加缝网复杂度"给出相反答案，分歧源于各自固定了不同的约束（相同注入时间 vs 相同注入体积、是否含孔眼摩阻/限流、是否显式 DFN、是否滤失主导）。这本身就是可研究空白。 |
| 接入现有代码是否可行？ | **可行，且比初判更省力**。阶段 B 原则上**不需要新增 Fortran**：`flow.source` + `source_scaling` 已能给出常排量/线性变排量，`FractureToughness` / `internalLength` / `toughnessAnisotropyMatrix` / `damage_type` 已覆盖可压性，`Gate 4.5` 能量分解与 EXO 输出已覆盖目标量。需要的是**新 deck + 新判据 + CI 接线**。 |
| 最大的技术障碍是什么？ | **不是代码，是相似律**：现 deck 的 `k_min=5e-11 m²` 让 `T=c_v t/H²≈1`（1 s 跑完固结），而现场页岩 1 h 注入 `T≈3e-8~1.6e-5`。要复现机制必须把 `k` 降 **6–8 个量级**并联合调整 `S`，否则压力瞬间扩散、注入压力上不去。 |
| 平台能证明什么 / 不能证明什么？ | 能证明**耗散分区（黏性↔韧性）与孔弹时间尺度耦合下的排量响应、能量分配、损伤形态**；**不能**独立复现滤失主导区、井筒/孔眼摩阻、多簇分流、支撑剂、DFN、地应力张量、现场产能递减。 |
| 建议下一步 | **阶段 B0 已落地**：源项量纲定论 `[1/s]`、相似校准双 deck、V1–V8 判据与 CI 硬闸门（第 11 节），待首轮 CI 验收后决定是否启动 B1 主实验。**仍待用户确认**：是否提供区块实测参数、B2 之后是否允许新增 Fortran。 |

### 对上一轮（2026-11-10 首版能力盘点）的三处更正

| # | 上一轮表述 | 静态核查后的事实 | 证据 |
|---|---|---|---|
| 1 | "`vDefHF.F90` 不在当前 CI 构建链" | **错**。CI 每轮都编译 `vDef vDefP vDefUpa vDefBT vDefHF`；但**从不运行**，且 `vDefHF` 链接不含 `MEF90_FLOW_OBJS`（无孔隙渗流/滤失） | `.github/workflows/build.yml:389`；`vDef/Makefile:58-62` |
| 2 | "注入源项通路从未跑过" | **不精确**。已有 `terzaghi1D_source.yaml` 的 `source: 1.0e-6`（Gate 4.4 验证 `sourceSec` 接线与 `W>0`）；缺的是**排量扫描与损伤–渗流耦合下的注入驱动** | `vDef/data/terzaghi1D_source.yaml:143`；`vDef/vDef.F90:598-608` |
| 3 | "IV 类强滤失 ⇒ 单靠排量失效" | **需撤回**。上一轮的 `C_L=2e-5~6e-4 m/√s` 是假设值；由 `k,S,μ,Δp` 物理反算的 Carter 系数为 `1.1e-7~2.6e-6 m/√s`，`Cdim` 全部 ≪1，**四类页岩 1 h 内全在存储主导区**。现场高滤失来自 DFN 宏观导流，不是基质渗流 | 本文第 4.2 节离线计算 |

> 第 3 条改变了研究范围的边界，但**强化而非削弱**了本研究的必要性：它说明"排量–可压性匹配"在致密页岩里主要是**耗散机制 + 孔弹时间尺度**的问题，恰好是本平台能严格控制的两个自由度；而"滤失主导"必须交给 DFN/商业软件。

### B0 实施后的两处补充更正（2026-11-10 同日）

| # | 本文原表述 | 实施核查后的事实 | 证据 |
|---|---|---|---|
| 4 | 第 3.2 节：源项「量纲不一致（B0 必须先校准）」，注释 `[Pa/s]` 与 bag 注册 `1/s` 悬而未决 | **已定论：`f` 是 `[1/s]`**。`-div(A∇p)=S dp/dt+f` 中 `S dp/dt` 为 `[1/s]`；bag 注册 `Volumetric source/sink rate (1/s)`；`MEF90FlowEnergy` 计算 `W=∫f·p dΩ`。旧 deck 注释 `[Pa/s]` 是错的。`f/S` 才是压力爬升率 `[Pa/s]`，`Q_total=f·Area` 是**单位厚度体积流量 `[m²/s]`**，换算 3D 需乘缝高 `H_f` | `m_MEF90_FlowCtx.F90:791`；`m_MEF90_FlowAssembly.F90:499-615`；`check_gate45_energy.py:35-70` |
| 5 | 第 5.3 节「用 PyYAML 校验单行差异契约」 | **CI runner 无 numpy / scipy / PyYAML**（现有全部 `check_gate4*.py` 的 docstring 均明确声明）。B0 判据改用纯标准库 `code_lines()` 剥离注释后逐代码行比较；本机 PyYAML 只用于提交前离线复核，不进 CI | `check_gate4*.py` docstring；`check_gate415_rate_fracability.py` 的 `code_lines` / `code_diff` / `block_values` |

> 第 4 条把第 4.3 节风险 R2 从「待标定」降级为「已解析 + 可证伪」：B0 不靠注释换算，而是由真实 EXO 的 `pmax` 斜率反解 `S_eff=f/slope` 与 deck 的 `S_num` 比对，并显式拒绝 `slope≈f` 的 `[Pa/s]` 解释（两者相差 `1/S = 6.67e9` 倍）。第 5 条是工程约束更正，不影响物理结论。

---

## 1. 问题界定：把"匹配"写成可证伪的约束满足问题

不把"匹配"定义成"排量越大越好"或一条拟合曲线，而是：给定储层可压性向量 $\mathbf{m}=(E',\nu,G_c,\ell,\mathcal{A},k,S,\Delta\sigma,\sigma_{weak})$，求可行排量区间 $Q^\*$ 使下列约束同时成立：

| 编号 | 约束 | 物理含义 | 本平台能否判定 |
|---|---|---|---|
| C1 | $p_{net}\ge \Delta\sigma+\sigma_{weak}$ | 净压力足以激活天然裂缝/弱面 | **能**（相对量纲下） |
| C2 | $\sigma_{h}+p_{net}+p_{perf}\le p_{lim}$ | 不超泵压/套管限压 | **部分**（$p_{perf}$ 不在模型内） |
| C3 | $V_{leak}/V_{inj}\le \eta_{max}$ | 液体效率 / 可加砂 | **受限**（无独立 Carter 系数，见 4.2） |
| C4 | $w\ge 3d_p$、缝高不穿层 | 支撑剂可进入、缝高受控 | **不能**（无支撑剂、无层状地应力） |

因此本平台的产出应是**"存在性 + 窗口宽度"的机制图版**（可压性 → $\mathcal{K}$ 与 $T$ → 可行 Q 区间），而不是现场施工参数表。现场侧由阶段 C 用施工曲线/微地震/产剖闭合。

---

## 2. 相关研究成果梳理

### 2.1 五条研究线与代表作

| 线 | 代表成果 | 核心贡献 | 与本项目的关系 |
|---|---|---|---|
| L1 变分/相场断裂源头 | Francfort & Marigo (1998)；Bourdin, Francfort & Marigo (2008), *J. Elasticity* 91 | Griffith 准则改写为全局能量最小化，裂纹路径无需预设 | `vDef` / `m_DefMech` 的理论母体 |
| L1b 相场水力压裂 | **Bourdin, Chukwudozie & Yoshioka (2012), SPE 159154-MS**, SPE ATCE, pp.1442–1452, DOI `10.2118/159154-MS`；Chukwudozie, Bourdin & Yoshioka (2019), *CMAME* 347:957–982, DOI `10.1016/j.cma.2018.12.037`；Yoshioka & Bourdin (2016), 与油藏模拟器耦合；Mikelić/Wheeler/Wick 系列；Miehe & Mauthe 多孔介质相场 | 把水力加载写成压力做功 / 统一裂缝–孔隙流动方程；用解析解验证；可处理分叉、合并、转向 | **本项目上游参考实现**。注意 2012 那篇明确假设：不渗透、无孔隙、不可压流体、**无滤失**、缝内压力均匀且仅由总注入体积决定 |
| L1c 综述与开放问题 | **Heider (2021), *Eng. Fract. Mech.* 253:107881**, DOI `10.1016/j.engfracmech.2021.107881` | 系统梳理裂缝宽度估计、驱动力、数值实现；列出开放问题 | 权威背书三条空白：①受控实验/现场条件验证不足；②**低黏流体 + 高排量下的振荡孔压与阶跃扩展**现有 PFM 多为平滑稳态、研究不足；③AI/数据驱动与 PFM 结合是未来方向 |
| L2 传播机制分区 | Savitski & Detournay (2002), *IJSS*；Garagash & Detournay；Bunger, Detournay & Garagash (2005) 滤失–韧性；Adachi & Detournay (2008) 滤失–黏性；**Detournay (2016), *Annu. Rev. Fluid Mech.* 48:311–339** | 四类渐近顶点（存储/滤失 × 黏性/韧性），参数空间轨迹；无量纲韧性 $\mathcal{K}$ 与滤失 $\mathcal{C}$ | 提供**可证伪的标度锚**，是 Gate 4.15 判据的理论依据 |
| L2b 给定流量的相场孔弹 | **Santillán, Juanes & Cueto-Felgueroso (2018), *JGR: Solid Earth* 123(3):2127–2155**, DOI `10.1002/2017JB014740` | 给定注入流量的相场孔弹模型；发现**三种排量模式**（低于阈值不扩展 / 中等流量达到稳态缝长 / 高流量分叉）；用量纲分析塌缩结果；做了网格与阈值无关性检验 | **最直接的方法学先例**，几乎就是本研究的目标，差别在：他们聚焦滤失–韧性主导与五点井网，本项目聚焦可压性参数场 + 硬判据验证范式 |
| L3 可压性评价 | **Rickman, Mullen, Petre, Grieser & Kundert (2008), SPE 115258**；Sondergeld et al. (2010) 矿物法；Jin et al. (2014) 能量耗散法；Yuan et al. (2013) 引入断裂韧性；Sharma, Chopra & Lines (2019) 用应变能密度 + 韧性替代脆性指数 | 从 E–ν 脆性指数 → 矿物/成岩/天然裂缝综合 FI → 能量/韧性物理量 | 支持本文把"可压性"落到 $G_c,\ell,E',\mathcal{A},\Delta\sigma$，而不是无量纲经验 BI |
| L3b 对脆性指数的批判 | 教科书式批评（无守恒律、动态/静态模量差、忽略地应力）；Sharma et al. (2019) 指出 E 与 ν 正相关区块里 BI 失去分辨力 | 脆性指数只是相关性，不是物理量 | 本研究用能量/韧性参数直接接入相场，避开该缺陷，可作为论文的一个立论点 |
| L4 完井侧：限流与孔眼 | Crump & Conway (1988), *JPT* 40(8):1041–1048；El-Rabba, Shah & Lord (1999), *SPE Prod. & Facil.* 14(1):63–71, DOI `10.2118/54533-PA`；Long & Xu (2017), *SPE J.* 22(2):645–659, DOI `10.2118/185173-PA`；Lecampion & Desroches (2015)；XLE 极限限流设计（长庆油田，2023 Springer 章节） | 孔眼摩阻是簇效率**最关键**因素；孔眼冲蚀会削弱限流；高排量利于均匀进液但加速冲蚀 | 支撑上一轮"Q↑ 的增压主要来自孔眼摩阻"的归因；同时说明**本平台不含该通道**，故 Q 在产能模型里必须与孔数/孔径/簇数做交互项 |
| L5 排量效应的实验/数值证据 | 见 2.2 分歧表 | 见 2.2 | 研究空白的直接证据 |
| L5b 压降反演 | **Nolte (1979, 1986)** G 函数；Barree et al. (2009)；Wang & Sharma 系列 DFIT；CPIT 恒压注入测试（2022） | 由停泵压降反演闭合应力、$C_L$、液体效率、渗透率 | 阶段 C 的现场入口：$C_L$ 与可压性由施工曲线反演，再喂给阶段 B |
| L5c 数据驱动 | 机器学习压裂参数优化、88 口井产能预测（本项目中石油算法大赛任务） | 从施工参数到产量轨迹 | 阶段 C：把 $\mathcal{K}$、$\mathcal{L}$、$p_{perf}/p_{net}$ 作为**物理交互特征**入模 |

### 2.2 关键分歧：文献对"Q↑ 是否有利于复杂缝网"结论相反

| 来源 | 结论 | 隐含固定的约束 |
|---|---|---|
| Beugelsdijk et al. (2000)；Guo et al. (2014) | **低排量**更易与天然裂缝作用、SRV 更大 | 天然裂缝发育、相同注入体积、含 DFN |
| Li & Einstein (2019) AE 实验 | **降低排量**产生更复杂微裂纹网络、剪切机制占比更低 | 室内小尺度、无围压或低围压 |
| Gil et al. (2011) | 低排量 → 剪切破坏；高排量 → 拉伸破坏 | DFN 数值模型 |
| Riahi & Damjanac (2013) | **给定注入体积**：低排量增大受影响 DFN 面积；HF 形成后高排量反而不利于剪切激活 | 显式 DFN、相同注入体积 |
| Frontiers (2021) Global CZM；MDPI *Processes* (2026) 14(17):2808 | **增大排量**提高净压/破裂压力、增加分支与改造面积（3→6 m³/min：面积 +23.4%、分支 6→11、破裂压力 +19.9%） | 相同注入时间、2D 均质、无孔眼摩阻；后者明确声明"配置特定、非网格无关普适值" |
| 相场双裂缝数值 (2024, *RSM*) | 排量对**峰值注入压力影响很小**，但加快扩展、加长路径 | 饱和多孔弹性、KGD 验证 |
| 相场双孔定向压裂 (2026, *EFM*) | 高排量 + 高黏度缩短穿透时间但**提高破裂压力**、陡化压力–时间斜率 | 含预制裂缝、双孔 |
| Santillán et al. (2018) | 存在**排量阈值**：低于阈值不扩展 / 中等达稳态 / 高分叉 | 五点井网、滤失–韧性主导 |

**归因**：分歧不是"谁对谁错"，而是三组被混淆的控制量——
① **固定注入时间 vs 固定注入体积**（前者 Q↑ 等价于多注入能量，后者 Q↑ 等价于更短时间）；
② **储层净压 vs 井底压力**（孔眼摩阻/限流把 Q 的平方项加进 BHTP，而净压在黏性主导下与 Q 无关）；
③ **基质渗流 vs DFN 宏观导流**（决定滤失是否成为独立约束）。

**本研究可贡献点**：用无量纲相似 + 硬判据把这三组控制量逐一隔离，给出"在什么约束下哪个结论成立"的判别图版；并用相场能量分解（弹性/表面/内聚/注入功）解释机制，而不是只报缝长或面积。这正是 Heider (2021) 指出的"验证严格度不足"与"低黏高排量非稳态"两处空白的交集。

### 2.3 核心标度关系（用于判据，不用于现场配产）

$$
E'=\frac{E}{1-\nu^2},\quad \mu'=12\mu,\quad q=\frac{Q}{H_f},\quad G_c=\frac{K_{Ic}^2}{E'},\quad \sigma_c\sim\sqrt{\frac{E'G_c}{\ell}}
$$

$$
\mathcal{K}=\frac{4K_{Ic}}{\sqrt{\pi}}\left(\frac{1}{3qE'^3\mu'}\right)^{1/4},\qquad
\mathcal{C}=2C_L\left(\frac{E't}{12\mu q^3}\right)^{1/6},\qquad
T=\frac{c_v t}{H^2},\ c_v=\frac{k}{\mu S}
$$

- $\mathcal{K}\propto q^{-1/4}$：**加排量把系统推向黏性主导**，与直觉相反；
- 黏性主导 $p_{net}\sim(E'^2\mu'/t)^{1/3}$ 与 Q **无关**；韧性主导 $p_{net}\propto Q^{-1/3}$ 反而下降；
- 孔眼摩阻 $p_{perf}=\dfrac{\rho}{2C_d^2}\left(\dfrac{Q}{N_pA_p}\right)^2$ 是 Q 的**平方**项；
- 相场峰值应力 $\sigma_c\sim\sqrt{E'G_c/\ell}$ 把可压性参数与起裂门槛绑定，是无量纲化必须用的压力尺度。

---

## 3. 现有代码能力静态核查（零编译、`HEAD=027fa86`）

### 3.1 两个主程序的互补性（本次最重要的结构性发现）

| | `vDef/vDef.F90`（CI 生产路径） | `vDef/vDefHF.F90`（编译但从不运行） |
|---|---|---|
| 链接库 | `MEF90_OBJS + DEFMECH + HEATXFER + FLOW`（`vDef/Makefile:36-38`） | 仅 `MEF90_OBJS + DEFMECH + HEATXFER`，**无 `MEF90_FLOW_OBJS`**（`vDef/Makefile:58-60`） |
| 压力 | 全场 `SectionReal`/Vec，Biot 孔弹耦合 | 标量乘子：`VecScale(CrackPressure, CrackPressureSave(I3))`（`vDefHF.F90:408`） |
| 注入控制 | 给定 `flow.source`（1/s 体积源），压力是响应 | **体积控制注入**：secant 迭代令 `CrackVolume = time(step)`（`vDefHF.F90:387-427`） |
| 滤失 | 有（基质 Darcy + 损伤增强渗透率） | 无（Bourdin 2012 的不渗透假设） |
| CI 状态 | 12 个 2D deck + 1 个 3D deck 实跑，16 组 verdict | `build.yml:389` 编译，**从未执行** |

**判断**：`vDefHF` 的 secant 体积控制适用于"不渗透 + 无限导流"的解析基准（Sneddon / KGD 验证），**不适合作为排量–可压性主模型**——在有滤失的孔弹系统里，质量守恒是「注入 = 缝内储存 + 滤失 + 基质储存」，强行令 $V_{crack}=Qt$ 会破坏滤失守恒。反过来，**给定 `f_source` 才是带滤失系统的物理正确注入方式**。因此阶段 B 的主模型应选 `vDef + flow.source`，把 `vDefHF` 留作「无滤失极限」的解析锚与交叉校验。

> 这修正了上一轮"需把 `CrackVolumeControlled` 接进 `vDef`"的建议：**不必接，也不该接**。

### 3.2 排量入口

| 入口 | 位置 | 量纲 / 行为 | 可用性 |
|---|---|---|---|
| `flow.source`（cell-set） | `m_Flow/m_MEF90_FlowCtx.F90:791-792` | 注册说明 `Volumetric source/sink rate (1/s)` | **可直接用**：改 deck 数值即扫 Q |
| `source_scaling: constant / linear / file / null` | `FlowCtx.F90:722-723`；`FlowAssembly.F90:869-878` | `constant` 用 `SetSourceCst`；`linear` 再 `VecScale(source,timeVal)` ⇒ $Q\propto t$；`file` 从 EXO `source_offset` 逐步读场 | **可直接用**：变排量程序无需改 Fortran（`file` 需先离线生成带源项场的 EXO，见风险 R3） |
| `source` 参与残差 | `FlowAssembly.F90:171-180, 217-227` | Gate 4.4 已把符号修正为 `G - F`（正 deck 源 = 注入） | **已验证**：`terzaghi1D_source.yaml` + Gate 4.5 `W>0` |
| `source_offset` EXO 输出 | `FlowCtx.F90:711-712`；`FlowAssembly.F90:635-637` | 逐单元源场可离线读回 | **可直接用**做体积守恒判据 |
| `flow.viscosity` | deck `:57` | `μ` 显式进入 $k/\mu$ | **可直接用**扫 $\mu'$ |

⚠️ **量纲不一致（B0 必须先校准）**：方程头写 `-div(k/μ ∇p) = S dp/dt + f_source`（`FlowCtx.F90:16`），bag 注册为 `1/s`，而 `terzaghi1D_source.yaml:143` 注释写 `[Pa/s]`。判据 R1 的体积守恒必须先用一次实测把 `Q_total = Σ(f_i V_i)` 与 EXO 数值对上，不能凭注释直接换算。

### 3.3 可压性入口

| 参数 | 键名 / 位置 | 默认 / 现值 | 可用性 |
|---|---|---|---|
| 断裂能 $G_c$ | `FractureToughness`，`MEF90/m_MEF90_Materials.F90:505,600` | 默认 1.0 N/m | **可直接扫**（由 $K_{Ic}^2/E'$ 反算，见 4.2 表 A） |
| 相场长度 $\ell$ | `internalLength`，`:527,619`；deck `:98` = 0.2 m | 0.2 m | **可直接扫**；注意 $\sigma_c\propto\ell^{-1/2}$，$\ell/H$ 必须固定 |
| 韧性各向异性 | `toughnessAnisotropyMatrix`，`:507 (2D,3 分量) / :602 (3D,6 分量)` | 单位阵 | **可直接用**表征页岩层理（文献 L3/L5 的层理角度效应） |
| 弹性 | `hookeslaw_YoungsModulus` / `hookeslaw_PoissonRatio` / `hookeslaw_tensor` / `hookeslaw_planeStress`，`:520-522,615-616` | deck `:91` 起 | **可直接扫** $E,\nu$ |
| 损伤律 | `damage_type`（AT1/AT2/LinSoft…）、`damage_AT1exp_b`、`damage_LinSoft_k`，`DefMechCtx.F90:943-945` | deck 用 AT2 | **可直接扫**；换律会改变 $\sigma_c$ 的 O(1) 系数，需同律比较 |
| 塑性/剪切驱动 | `plasticity_type`、`drivingForce_type`、Drucker–Prager 模块 `m_MEF90_DefMechDrivingForceDruckerPrager.F90:75-76,185-186` | — | **可用**：对应文献"低排量→剪切、高排量→拉伸"的机制隔离 |
| 损伤→渗透率 | `k(v)=k_{min}+(k_{max}-k_{min})\max(v,v_{res})`，`m_MEF90_Flow.F90:209-210,327`（含 3D 孪生） | 20× 增强 | **可直接扫**裂缝导流；也是本平台唯一能间接代表 DFN 的旋钮 |
| 孔弹 | `biot_coeff`、`storativity`、`perm_min/max` | deck `:57-59` | **可直接扫** $T,\alpha$ |

### 3.4 输出与目标量（均已存在，无需新增 Fortran）

| 量 | 位置 | 用途 |
|---|---|---|
| 压力 / 位移 / 损伤 EXO | `disp_offset=3, damage_offset=2, source_offset=2`（deck `:74-76`） | 场级判据、逐位回归 |
| Darcy 耗散能 `E` + 注入功 `W` | `Gate 4.5: step N t T E ... W ...`，`vDef.F90:604-608` | 能量守恒、排量单调性 |
| 弹性能 / 表面能 / 内聚能 / 外力功 / 总机械能 | `vDef.F90:789-814` | **能量分配**：黏性 vs 韧性耗散的直接证据 |
| `Gate 4.12` 步间压降、`Gate 4.1` Biot 探针 | `vDef.F90` | 注入压力演化、孔弹耦合活性 |
| 3D `Displacement_Z`、`coordz` | CI #121 `out_dim3.gen`（15,106,908 B，11 步，10 条 `iterations 1`） | 3D 一致性通道已打通 |

### 3.5 从未被闸门覆盖的通路（Gate 4.15 的立项依据）

1. `source ≠ 0` 只在符号/接线层验证过（Gate 4.4/4.5），**从未用于排量扫描，也从未与损伤演化耦合**；
2. 所有 deck 的 `FractureToughness` / `internalLength` / `toughnessAnisotropyMatrix` 都是默认值或防回归置零，**从未作为自变量扫描过**；
3. `vDefHF` 编译但从未运行，`CrackVolumeControlled` 与 `WorkControlled` 路径零覆盖；
4. `source_scaling: linear / file` 从未实跑；
5. 3D 通道（CI #121 已绿）目前只跑固结，未跑注入。

> 与 Gate 4.14 之前「3D 编译孪生存在但从未执行」是同一类验证空白，处置范式可直接复用：新 deck + 可证伪判据 + CI 硬闸门 + 延迟 `exit 1` + 生产零扰动证明。

---

## 4. 可行性评估

### 4.1 能力 / 缺口矩阵

| 维度 | 分级 | 依据与代价 |
|---|---|---|
| 常排量 Q 扫描 | 🟢 **零代码改动** | 复制 `terzaghi1D_source.yaml`，只改 `source` 数值；PyYAML 可验「单行差异」契约 |
| 变排量程序 $Q(t)$ | 🟢 **零代码改动**（`linear`）/ 🟡 `file` 需离线 EXO 生成器 | `linear` 给 $Q\propto t$；`file` 需 CDF2 写能力（风险 R3） |
| 可压性扫描 $G_c,\ell,E,\nu,\mathcal{A}$ | 🟢 **零代码改动** | Materials bag 键已注册，deck 可直接写 |
| 损伤律 / 剪切驱动机制隔离 | 🟢 **零代码改动** | `damage_type`、`drivingForce_type` 已注册 |
| 能量分配与体积守恒判据 | 🟢 **零代码改动** | Gate 4.5 + EXO source/damage 场足够 |
| 2D / 3D 一致性 | 🟢 **零代码改动** | CI #121 已打通 3D 诊断运行（`-np 1` + 直接 LU，1800 s） |
| 孔弹时间尺度 $T$ 相似 | 🟡 **需 deck 参数大改**（$k$ 降 6–8 量级、$S$ 调 0.33×） | 见 4.2；可能触发求解器条件数变化，需复核 KSP 收敛（#120 教训） |
| 滤失主导 $\mathcal{C}\gtrsim1$ | 🔴 **做不到（物理上不独立）** | $C_L$ 由 $k,S,\mu,\Delta p$ 涌现，与 $T$ 共用 $k/(\mu S)$；只能靠等效 $k$ 放大代表 DFN，不能独立设 $C_L$。若要真 DFN，需新增 Fortran/网格，超出阶段 B |
| 孔眼摩阻 / 限流 / 多簇分流 | 🔴 **模型外** | 无井筒、无孔眼、无簇；只能作为**离线后处理**叠加 $p_{perf}(Q)$，并在论文里显式声明是外部关联 |
| 支撑剂 / 缝宽 $w\ge3d_p$ / 缝高穿层 | 🔴 **模型外** | 无支撑剂、无层状地应力、无 $K_{Ic}$ 分层 |
| 绝对应力 / 现场 BHTP | 🟡 **只能相对量纲** | $\sigma_c\sim1.3$–4.5 MPa（表 A），远小于现场 $\sigma_h\sim50$ MPa；$\sigma_h$ 只能作为 offset 后处理 |
| 现场产能递减 | 🔴 **模型外** | 交阶段 C 的 88 口井数据驱动线 |

### 4.2 无量纲相似：现场 → 平台秒级实验的定量映射

**设定**：现场 $H_f=H=50$ m、$t_f=3600$ s、$\mu=10^{-3}$ Pa·s、$\Delta p=3$ MPa、$S_0=4.6\times10^{-10}$ Pa⁻¹；数值 $H_n=10$ m、$\ell_n=0.2$ m、$t_n=1$ s。四类页岩参数与上一轮首版一致。

**表 A：由 $K_{Ic}$ 反算的相场材料参数**

| 类型 | $E'$ (GPa) | $G_c=K_{Ic}^2/E'$ (N/m) | $\sigma_c\sim\sqrt{E'G_c/\ell}$ (MPa) | $c_v$ (m²/s) | $T_{field}=c_vt/H^2$ |
|---|---|---|---|---|---|
| I 高脆低韧/致密 | 42.0 | 8.56 | 1.34 | 2.17e−8 | 3.13e−8 |
| II 中等脆性 | 32.0 | 37.8 | 2.46 | 1.09e−7 | 1.57e−7 |
| III 低脆高韧/富黏土 | 21.7 | 184.3 | 4.47 | 2.17e−7 | 3.13e−7 |
| IV 天然裂缝发育（等效 $k$） | 37.1 | 17.2 | 1.79 | 1.09e−5 | 1.57e−5 |

**表 B：无量纲分区（$C_L$ 由 $C_L\simeq\Delta p\sqrt{kS/(\pi\mu)}$ 反算）**

| 类型 | $C_L$ (m/√s) | $\mathcal{K}$ @4→20 m³/min | $\mathcal{C}$ @4→20 | 耗散分区 | 存储分区 |
|---|---|---|---|---|---|
| I | 1.15e−7 | 0.175 → 0.117 | 0.0030 → 0.0014 | 黏性主导 | 存储 |
| II | 2.57e−7 | 0.394 → 0.264 | 0.0065 → 0.0029 | 过渡→黏性 | 存储 |
| III | 3.63e−7 | 0.959 → 0.641 | 0.0086 → 0.0038 | 过渡 | 存储 |
| IV | 2.57e−6 | 0.256 → 0.171 | 0.0664 → 0.0297 | 黏性主导 | 存储 |

$\mathcal{K}$ 全部落在 **0.117–0.959**：无一进入韧性主导（$\mathcal{K}>3$），且 $\mathcal{K}\propto Q^{-1/4}$ 使加排量更黏性主导——与上一轮趋势一致，现在有文献公式背书。$\mathcal{C}\ll1$ 说明**在纯基质渗流假设下滤失约束根本不激活**。

**表 C：固定 $(T,\mathcal{K},\mathcal{C},\ell/H)$ 反解数值 deck 参数**

| 类型 | $k_{num}$ (m²) | $S_{num}$ (Pa⁻¹) | $k_{num}/k_{deck}$ | $S_{num}/S_{deck}$ |
|---|---|---|---|---|
| I | 4.70e−19 | 1.50e−10 | 9.4e−9 | 0.326 |
| II | 2.35e−18 | 1.50e−10 | 4.7e−8 | 0.326 |
| III | 4.70e−18 | 1.50e−10 | 9.4e−8 | 0.326 |
| IV | 2.35e−16 | 1.50e−10 | 4.7e−6 | 0.326 |

**三条硬结论**：

1. $\mathcal{K}$ **不含** $H$ 与 $t$ ⇒ 要复现机制必须 $q_{num}=q_{field}$，即 $Q_{num}=Q_{field}\cdot H_n/H_f=Q_{field}/5$；总排量按几何比缩放，不能照抄现场数字。
2. $T\propto k/(\mu S)$ 与 $C_L\propto\sqrt{kS/\mu}$ **共用 $(k,S)$** ⇒ 滤失不是独立旋钮；必须联立反解 $(k_{num},S_{num})$，且 $(T,\mathcal{C})$ 只有两个自由度。这是平台的结构性限制。
3. 现 deck $k_{min}=5\times10^{-11}$ m² 比相似所需高 **6–8 个量级**；若照搬 51 步/1 s，压力会在第一步内扩散完，注入压力与损伤都起不来。**阶段 B0 必须先做这个尺度校准**。

### 4.3 风险清单

| ID | 风险 | 触发条件 | 缓解 |
|---|---|---|---|
| R1 | `k` 骤降 6–8 量级使矩阵条件数恶化，KSP 不收敛 | 阶段 B0 | 复用 CI #121 的单 rank 直接 LU 诊断配置；`-flow_ksp_monitor` + `-flow_ksp_converged_reason`；`timeout` 与 `s415_rc` 必须参与裁决（#120 教训） |
| R2 | 源项量纲注释与 bag 注册不一致，体积守恒判据算错 | R1 判据 | B0 先做一次 $Q_{total}$ ↔ EXO 实测标定，冻结换算常数 |
| R3 | `source_scaling: file` 需要离线生成带源项场的 EXO，`parse_cdf2.py` 现为只读 | 变排量程序 | 阶段 B1 先用 `constant` + `linear`；`file` 留到 B2，届时新增纯 Python CDF2 writer（不触碰 Fortran） |
| R4 | 相场 $\ell$ 与网格 $h$ 强耦合，$\ell/H$ 固定但 $h/\ell$ 未收敛 | 可压性扫描 | 固定网格、只改 $\ell$ 时须做 $h/\ell\le0.5$ 检查；把网格敏感性作为 R6 负例之一 |
| R5 | 新 deck 扰动生产判据 | 每次新增 | 沿用 4.1b/4.9/4.14 范式：生产零扰动、15 组 verdict 逐字节一致、`out_main.gen` 物理数值逐位不变 |
| R6 | 把「单 rank 直接 LU」误读成并行规模证据 | 论文写作 | 明确写成「控制求解代价的诊断配置」，与 `RESEARCH_POSITIONING.md` 质疑 3 保持一致 |
| R7 | 把离线标度当成求解器验证结论 | 论文写作 | 所有标度结论标注「待 Gate 4.15 数值实验检验」，不得写成已验证 |

---

## 5. 技术路线

### 5.1 三阶段总览

| 阶段 | 输入 | 动作 | 产出 | 退出条件 | 是否需改 Fortran |
|---|---|---|---|---|---|
| **A 离线标度**（本轮已完成首版） | 文献公式 + 区块参数 | 修正 $C_L$、构造 $\mathcal{K}$/$\mathcal{C}$/$T$ 图版与相似映射、排出敏感性顺序 | 本文 + 参数映射表 + 88 口井物理交互特征清单 | 图版可离线复算、参数区间与代码可达区间对齐 | 否 |
| **B0 校准冒烟（已落地接线，待 CI 验收）** | 表 C 的 $(k,S)$ | **2 个新 2D deck**：零源对照 `terzaghi415_scale.yaml`（4 行差异）+ 常排量注入 `terzaghi415_inject.yaml`（1 行差异，11 步）；V1–V8 硬判据 + CI 延迟硬裁决 | 冻结量纲常数 `[1/s]`、`Q_total=1.5e-4 m²/s`、`T=1.57e-7`、`S_eff`；求解器配置、EXO 槽位 | 体积守恒闭合、零源冻结、注入线性、无 fatal、KSP 每步收敛、生产零扰动 | 否 |
| **B1 主实验 / Gate 4.15** | B0 校准结果 | $Q\times G_c\times k\times\mu$ 正交矩阵；新判据 R1–R8；CI 硬闸门 | 16+1 组 verdict、机制分区图版、能量分配曲线 | 全部 R 判据 PASS、生产零扰动 | 否 |
| **B2 扩展** | B1 绿后 | 各向异性韧性（层理角）、`linear` 变排量、3D tetra4 一致性；可选 `file` 程序 | 层理效应图版、3D 对照 | 同上 + 3D 通道零扰动 | 否（`file` 需 Python writer） |
| **C 现场与数据驱动** | 施工曲线、DFIT/G 函数、微地震、产剖、88 口井 | 反演 $C_L$ 与可压性；把 $\mathcal{K},\mathcal{L},p_{perf}/p_{net}$ 作为交互特征入产量模型；统计检验 | 排量优选表、论文验证段 | 现场数据闭合 | 否（本平台外） |

### 5.2 Gate 4.15 判据设计（可证伪，沿用现有硬闸门范式）

| 判据 | 检查内容 | 数据来源 | 负例设计 |
|---|---|---|---|
| R1 注入体积守恒 | $\sum_i f_iV_i\Delta t$ 与 压力储存 + 边界流出 + 损伤增强储存 之差 < tol；$Q_{total}$ 与 B0 标定常数一致 | `terzaghi415_*.log` 的 Gate 4.5 行 + EXO `source`/pressure | 篡改一个源项数值 / 去掉一步 |
| R2 排量单调性 | Q 增大 ⇒ 峰值 $p$、注入功 $W$、总表面能不减；零源对照 `W=0` | 多 deck 日志 | 打乱 Q 与 p 的对应、把零源 deck 的 W 置非零 |
| R3 可压性单调性 | $G_c$ 增大 ⇒ 起裂延迟、损伤总长/表面能占比下降（固定 Q） | EXO damage 场 + 能量行 | 反转 $G_c$ 顺序、令损伤场与对照逐位相同 |
| R4 机制分区标度 | 双对数斜率：黏性主导段 $\partial\ln p_{net}/\partial\ln Q\approx0$；提高 $G_c$ 或 $\mu$ 后斜率向 $-1/3$ 移动；$\mathcal{K}$ 与斜率单调对应 | 多 deck 峰值净压 | 用常数 $p_{net}$ 冒充全部机制、破坏 $\mathcal{K}\propto Q^{-1/4}$ |
| R5 各向异性效应 | `toughnessAnisotropyMatrix≠I` ⇒ 损伤路径偏向低韧性方向，且与各向同性对照有可测差 | EXO damage | 令两 deck 输出逐位相同、把矩阵改回单位阵却报有差异 |
| R6 网格/长度尺度 | $h/\ell\le0.5$；$\ell/H$ 固定的两组给出同量级 $\mathcal{K}$ 指标 | 网格元数据 + 静态扫描 | 用未收敛网格蒙混 |
| R7 生产零扰动 | 15 组既有 verdict 与 CI #121 逐字节一致；`out_main.gen` 物理数值逐位不变（仅 QA changeset 字段允许变） | artifact 对比 | 伪造 verdict |
| R8 静态反回归 + CI 裁决 | 源码/量纲注释、`flow.source` 注册、`source_scaling` 分支、CI 接线、`s415_rc` 参与延迟 fail、`timeout` 与 KSP monitor 存在 | 源码静态扫描 | 删任一接线（复用 #119/#120 两类失败） |

判据脚本命名 `vDef/data/check_gate415_rate_fracability.py`，自带 `--selftest`（B0 实装 1 正例 + 20 负例），CI 接线完全复用 `check_gate414_dim3.py --gate` 的模板与延迟 `exit 1` 结构。**B0 只实装准入判据 V1–V8**（运行完整性、相似参数冻结、零源存储主导、注入线性与 `S_eff` 量纲标定、注入功、KSP、静态/deck/CI 契约、生产隔离），它们是 R1/R4/R8 的前置子集；**R2/R3/R5/R6（排量单调性、可压性单调性、各向异性、网格尺度）需 B1 的 Q×Gc×k×μ 正交矩阵才有数据**，留到 B1 扩展，不在此轮伪造判据。

### 5.3 CI 与零扰动约束

- 新 deck 一律**单行差异**契约（B0 为「校准 deck 对基准 4 行白名单」+「注入 deck 对校准 deck 1 行」两级），由判据的纯标准库 `code_lines()` / `code_diff()` 校验（CI 无 PyYAML，见第 0 节补充更正 5）；
- 所有新运行放在 `Run Terzaghi` 后段，生产 12 个 2D deck 与 3D 诊断运行**保持原配置不动**；
- 捕获 `s415_rc`（mpirun 退出码）并接入延迟 hard fail，避免 #120 的「rc=124 只 echo」漏洞；
- artifact 增 `gate415_verdict.log` 与新 deck 的 `.log`/`.gen`；
- 零本地编译、零本地求解器：所有验证靠 `git diff --check`、`ast`、PyYAML（**仅本机离线复核，不进 CI**）、10 个 run block `bash -n`、判据 selftest、真实 artifact 离线回放。

### 5.4 决策点

1. **是否立项 Gate 4.15** —— 用户以「推进」默认立项，阶段 B0 已实施落地（第 11 节）；`FLOW_FIX_CHECKLIST.md` 现为「24 项已落地 + **1 项在册待验收**」。**B1 主实验是否启动**待 B0 首轮 CI 验收结果；
2. **区块实测参数**（$E,\nu,K_{Ic},C_L,k,\Delta\sigma$、孔数/孔径/簇数、Q 范围）：**仍待确认**。有则把表 A–C 重算成可直接用的排量优选表；无则继续用四类典型值；
3. 阶段 B2 之后是否允许新增 Fortran（例如把等效 DFN 通道、显式 Carter 滤失或井筒限流耦合进 `m_Flow`）—— **仍待确认**，这会突破"零代码改动"的优势，需单独评估。

---

## 6. 与作者三条研究线的接口

| 研究线 | 接口 | 可复用范式 |
|---|---|---|
| ISC + NaCl 盐析封堵（中文 SCI 手稿） | 「逐单元参数场 + 损伤→渗透率」与盐析堵塞导致 $k$ 突变同构；TG/SEM/高温岩心驱替给出 $k(v)$ 的实测标定 | 硬判据 + 能量分解 + 生产零扰动验证链 |
| CO₂–水–岩石 LBM | Gate 4.15 输出的 $k(v)$ 与损伤场可直接作为 LBM 升尺度边界条件；反过来 LBM 给出 $k_{min/max}$ 的微观依据 | 无量纲相似（$T,\mathcal{K},\mathcal{C}$） |
| UCG 多场耦合 | 排量/注入功 ↔ 气化剂注入速率与合成气品质；共用「排量程序 + 温度相关 $\mu(T)$ + 逐单元渗透率」 | `source_scaling` 变排量 + `UpdateMuFromTemperature`（已有） |
| 中石油算法大赛 88 口井 | 把 $\mathcal{K}$、$\mathcal{C}$、$p_{perf}/p_{net}$ 作为**物理交互特征**：Q 必须与孔数/孔径/簇数交互，单独入模会系统性误读（见 2.2/2.3） | 机理归因 → 特征工程 |

---

## 7. 能力边界与审稿质疑预案

| 质疑 | 回应 |
|---|---|
| 「2D / 10 m / 单 rank，凭什么谈现场排量？」 | 定位为 **verification + 机制判别**，不是现场尺度 validation；说服力来自无量纲相似与逐位回归，不是规模。3D tetra4 已由 CI #121 实跑（2635 节点、11 步、`Displacement_Z`），但明确声明单 rank 直接 LU 是诊断配置、非并行规模证据 |
| 「没有井筒/孔眼/多簇，排量结论是否失真？」 | 正是本研究的**隔离设计**：把 Q 的储层净压通道与完井摩阻通道拆开，$p_{perf}(Q)$ 作为外部关联后处理叠加；这恰是商业全耦合软件难以做到的机理归因 |
| 「滤失主导区没覆盖，结论不完整」 | 明确界定：由 $k,S,\mu,\Delta p$ 反算 Carter 系数后 $\mathcal{C}\ll1$，致密页岩在 1 h 内本质存储主导；滤失主导需 DFN 宏观导流，超出本平台，交给阶段 C / 商业软件。把边界写出来比强行外推更可信 |
| 「$\ell$ 是数值参数还是材料参数？」 | 相场长度尺度决定 $\sigma_c\sim\sqrt{E'G_c/\ell}$，不可当纯正则化参数；本文固定 $\ell/H$ 并把 $h/\ell$ 收敛列为 R6 判据 |
| 「与 Santillán et al. (2018) 有何区别？」 | 他们做给定流量下的孔弹相场 + 五点井网 + 滤失韧性分区与量纲塌缩；本研究做**可压性参数场（$G_c,\ell,\mathcal{A},\Delta\sigma$）作为自变量**、能量分配归因、以及可证伪硬判据验证体系（16 组闸门 → 17 组），并接数据驱动产能线 |
| 「结论是否与文献冲突？」 | 不冲突：2.2 的分歧由三组被混淆的控制量解释，本研究给出判别条件而非又一条经验曲线 |

---

## 8. 与 `RESEARCH_POSITIONING.md` 的关系

`RESEARCH_POSITIONING.md` 回答「与现有研究相比，平台整体定位如何」，结论是**领先在参数入口 + 验证严格度，落后在工程广度与并行规模**。本文是该定位在「施工排量–可压性匹配」这一具体问题上的展开：再次确认——

- **领先项**：参数入口齐全（`$G_c,\ell,\mathcal{A},k,S,\mu,Q(t)$` 全部可从 deck 扫）、验证严格（16 组硬判据、逐位回归、生产零扰动）、能量分解可归因；
- **落后项**：完井/井筒/DFN/支撑剂/地应力张量/现场产能全部在模型外；
- 因此**不主张替代 CMG / ECLIPSE / TOUGH2**，定位为「可验证的力学–渗流机理判别骨架 + 参数可辨识性基准 + 升尺度边界条件生成器」。

---

## 9. 引用线索（2026-11-10 检索；卷期与 DOI 投稿前需二次核对）

**已核实到出版信息的**
1. Bourdin, B., Chukwudozie, C., Yoshioka, K. (2012). *A variational approach to the numerical simulation of hydraulic fracturing*. SPE ATCE, SPE **159154-MS**, pp. 1442–1452. DOI `10.2118/159154-MS`.（**注意：不是 140927**）
2. Chukwudozie, C., Bourdin, B., Yoshioka, K. (2019). *A variational phase-field model for hydraulic fracturing in porous media*. **CMAME 347**, 957–982. DOI `10.1016/j.cma.2018.12.037`.
3. Heider, Y. (2021). *A review on phase-field modeling of hydraulic fracturing*. **Engineering Fracture Mechanics 253**, 107881. DOI `10.1016/j.engfracmech.2021.107881`.
4. Santillán, D., Juanes, R., Cueto-Felgueroso, L. (2018). *Phase Field Model of Hydraulic Fracturing in Poroelastic Media: Fracture Propagation, Arrest, and Branching Under Fluid Injection and Extraction*. **JGR: Solid Earth 123(3)**, 2127–2155. DOI `10.1002/2017JB014740`.
5. Detournay, E. (2016). *Mechanics of hydraulic fractures*. **Annual Review of Fluid Mechanics 48**, 311–339.
6. Savitski, A.A., Detournay, E. (2002). *Propagation of a penny-shaped fluid-driven fracture in an impermeable rock: asymptotic solutions*. **Int. J. Solids and Structures**（PII `S0020768302004924`）.
7. Rickman, R., Mullen, M.J., Petre, J.E., Grieser, W.V., Kundert, D. (2008). *A practical use of shale petrophysics for stimulation design optimization*. **SPE 115258**.
8. Crump, J., Conway, M. (1988). *Effects of perforation-entry friction on bottomhole treating analysis*. **JPT 40(8)**, 1041–1048.
9. El-Rabba, A., Shah, S., Lord, D. (1999). *New perforation pressure-loss correlations for limited-entry fracturing treatments*. **SPE Production & Facilities 14(1)**, 63–71. DOI `10.2118/54533-PA`.
10. Long, G., Xu, G. (2017). *The effects of perforation erosion on practical hydraulic-fracturing applications*. **SPE Journal 22(2)**, 645–659. DOI `10.2118/185173-PA`.
11. Francfort, G.A., Marigo, J.-J. (1998); Bourdin, B., Francfort, G.A., Marigo, J.-J. (2008). *The variational approach to fracture*. **Journal of Elasticity 91**, 5–148 / 1012–1030（分册页码需核对）.

**需二次核对卷期/DOI 的**
12. Yoshioka, K., Bourdin, B. (2016). *A variational hydraulic fracturing model coupled to a reservoir simulator*（ScienceDirect PII `S1365160916301484`；期刊名待核）.
13. Nolte, K.G. (1979, 1986). G 函数 / 压降分析原始 SPE 论文号待核.
14. Barree, R.D. et al. (2009). G 函数导数与天然裂缝识别，SPE 号待核.
15. Carrier, B., Granet, S. (2012). 数值水力压裂参数空间综述（$\mathcal{K},\mathcal{C}$ 公式转引自该文，需核对原式系数）.
16. Bunger, A.P., Detournay, E., Garagash, D.I. (2005). *Toughness-dominated hydraulic fracture with leak-off*，IJRMMS 卷期待核.
17. Adachi, J., Detournay, E. (2008). 滤失–黏性主导自相似解，卷期待核.
18. Mikelić, A., Wheeler, M.F., Wick, T. 系列（2015–2019）；Miehe, C., Mauthe, S. 多孔介质相场（2016）；Wheeler, Wick & Wollner (2014) 增广拉格朗日——卷期/DOI 待核.
19. 中文/工程侧：XLE 极限限流设计（长庆油田，Springer 2023 章节）、Frontiers in Earth Science (2021) Global CZM、*Processes* (2026) 14(17):2808、*Rock Mechanics* (2024) 相场双裂缝——作者与页码待核.

**公式出处**：$\mathcal{K},\mathcal{C}$ 采用 Carrier & Granet (2012) 转引形式（见检索到的 ENPC 博士论文 HAL `tel-04155433` 第 1.2 式）；$G_c=K_{Ic}^2/E'$ 与 $\sigma_c\sim\sqrt{E'G_c/\ell}$ 为平面应变 Griffith / AT2 标度，O(1) 系数 $c_w$ 省略，**用于分区判别，不用于定量配产**。

---

## 10. 首版文档交付与版本状态（阶段 A）

- 新增本文 `RATE_FRACABILITY_ROADMAP.md`；`README.md` 分支说明加一行索引；`.workbuddy/memory/2026-11-10.md` 追加研究结论。
- **该轮未改任何 Fortran / deck / CI，未做本地编译、未跑求解器**；全部计算为纯 `numpy` 离线标度与源码静态阅读。
- 基线 `HEAD=027fa86`；该轮提交为纯文档（`3a819ab`），按 CI #115/#122 的处置不为其追加提交。
- 阶段 B0 的启动决策与实施记录见第 11 节。

---

## 11. 阶段 B0 落地记录（2026-11-10，待 CI 验收）

**结论先行**：用户以「推进」默认立项 Gate 4.15，B0 已完成**全部代码与 CI 接线、零 Fortran 改动、零本地编译、零本地求解器**；判据离线自检 1 正例 + 20 负例全 PASS。**尚未取得 CI 实跑结果，因此本文所有 B0 数值仍是离线标度与解析预测，不是求解器验证结论**（对应第 4.3 节风险 R7）。

### 11.1 交付清单

| 文件 | 性质 | 内容 |
|---|---|---|
| `vDef/data/terzaghi415_scale.yaml` | 新增 | 相似校准零源对照；相对 `terzaghi1D.yaml` **仅 4 行代码差异** |
| `vDef/data/terzaghi415_inject.yaml` | 新增 | 常排量注入；相对校准 deck **仅 1 行代码差异** `source 0.0→1.5e-5` |
| `vDef/data/check_gate415_rate_fracability.py` | 新增 | V1–V8 硬判据，`--gate` 读真实 artifact / `--selftest` 1 正例 + 20 负例；纯标准库 |
| `.github/workflows/build.yml` | 修改 +88 行 | 2 个 2D 诊断 run、`g415_rc` verdict、两个 `s415_*_rc` 与 `g415_rc` 延迟 hard fail、artifact 5 个新成员 |
| `FLOW_FIX_CHECKLIST.md` | 修改 | 4.15 B0 作为 `- [ ]` 在册条目；状态总览改「14 项收口 / 15 项在册、1 项待验收」 |
| `RATE_FRACABILITY_ROADMAP.md` | 修改 | 第 0 节补充更正 4/5、5.1 B0 行、5.2 判据分期、5.3 PyYAML 更正、5.4 决策点、本节 |

### 11.2 冻结的相似与量纲常数（离线复算，非求解器结果）

| 量 | 值 | 来源 / 含义 |
|---|---|---|
| $k_{num}$ | $2.35\times10^{-18}$ m² | 表 C，II 类；$k/k_{deck}=4.7\times10^{-8}$（低 7.3 个量级） |
| $S_{num}$ | $1.50\times10^{-10}$ Pa⁻¹ | 表 C；$S/S_{deck}=0.326$ |
| $c_v=k/(\mu S)$ | $1.5667\times10^{-5}$ m²/s | $H²/c_v=6.38\times10^{6}$ s，1 s 内扩散被冻结 |
| $T=c_vt/H²$ @1 s | $1.5667\times10^{-7}$ | **与现场 II 类 $t_f=3600$ s / $H_f=50$ m 的 $1.57\times10^{-7}$ 一致** |
| 边界层 $\sqrt{c_vt}$ | 3.96 mm | 远小于网格尺度 0.2 m，压力最大值不受排水面影响 |
| $f$ | $1.5\times10^{-5}$ s⁻¹ | `[1/s]`，**不是** `[Pa/s]` |
| $f/S$ | $1.0\times10^{5}$ Pa/s | 注入压力爬升率；`pmax` 从 1e5 → 2e5 Pa |
| $Q_{total}=f\cdot A$ | $1.5\times10^{-4}$ m²/s | **单位厚度**体积流量；乘 $H_f$ 才是 m³/s |
| $W_{an}(1\,\mathrm{s})$ | 30 | $f A (p_0+f/S)$，V5 的宽松量级锚 |

B0 刻意取保守 $f$（$\Delta p\sim p_0$）使响应保持线性、无损伤、可解析验证；现场排量窗口（$Q_{num}=Q_{field}/5$）留到 B1 扫描。这一点必须在论文里写明，否则会被误读成"已在现场排量下验证"。

### 11.3 判据与 CI 设计

- **V1–V8**：完整性 / 相似参数 / 零源冻结 / 注入斜率与 `S_eff` / 注入功 / KSP / 静态与 deck 契约 / 生产隔离。V4 从真实 EXO 反解 `S_eff=f/slope` 并要求 `|S_eff/S_num−1|<10%`，同时要求 `slope/f>1e6`，从而**可证伪地排除 `[Pa/s]` 解释**。
- **零扰动**：12 个生产 2D deck 与 3D 诊断运行配置一字未改；新 run 放在 `Run Terzaghi` 后段并各自 `cp` 到独立 `.log` / `.gen`，随后仍执行原有「恢复生产 mesh」步骤；V8 检查生产日志仍在固结且无 415 前缀泄漏。
- **CI #120 教训复用**：`s415_scale_rc`、`s415_inject_rc`、`g415_rc` 三者全部接入延迟 `exit 1`，超时留下的半截 EXO 无法假绿。`timeout -k 10s 600` + `-flow_ksp_monitor` + `-flow_ksp_converged_reason` 保证 V6 可观测。
- **CI 环境差异**：判据不依赖 numpy / scipy / PyYAML；deck 差异与标量解析用 `code_lines()` / `block_values()` 纯标准库实现。

### 11.4 提交前离线复核（全部通过）

| 检查 | 结果 |
|---|---|
| `git diff --check` | OK |
| 判据 `ast.parse` | OK |
| PyYAML 解析 3 份 deck + workflow | OK（`numstep=11/11/51`、`k=2.35e-18`、`S=1.5e-10`、`source=1.5e-5/0/0`） |
| 10 个 workflow `run` block `bash -n` | 10/10 PASS |
| 判据 `--selftest` | **PASS（1 正例 + 20 负例）** |
| 真实 CDF2 IO 烟测 | `Beam1x10-tri3.gen` = num_dim 2 / 356 节点 / sets [1,10,20]；`Beam1x10-tri3_out.gen` = 51 记录 / 4 nodal vars，Pressure 为 var1 |
| deck 代码行契约 | 主→校准 4 行（`numstep`/`perm_min`/`perm_max`/`storativity`），校准→注入 1 行（`source`） |

### 11.5 首轮 CI 的主要风险与预案

| 风险 | 触发表现 | 预案 |
|---|---|---|
| 真实 `pmax` 斜率与 $f/S$ 有偏差（Dirichlet 修正、BE 离散、边界层） | V4 FAIL | 先读 `gate415_verdict.log` 的 `g415_slope_measured` / `g415_S_eff`；若为系统性 O(%) 偏差则放宽 `TOL_SLOPE` 并在文档记录实测值，**不得**改 deck 去凑解析解 |
| `-np 2` 下 11 步只打印 <10 条 `Linear solve converged` | V6 FAIL | 与 CI #121 的 3D 步进计数对齐；必要时改为 `>= N_REC-2` 并补证据，或加 `-flow_ksp_monitor` 逐迭代证据 |
| 注入功 `W` 量级与 `f·A·p` 不符（单元积分/边界层权重） | V5 FAIL | 现有量级带已放到 0.05–20 倍；若仍越界，以 Gate 4.5 / 4.8 的既有离散守恒式替换该检查 |
| 生产日志 `Gate 3.2` 首末值与合成正例假设不同 | V8 FAIL | 用真实 #121 `/tmp/terzaghi.log` 离线回放校正阈值（该判据只做隔离，不应误伤） |

> **状态纪律**：本节所有结论标注「待 Gate 4.15 数值实验检验」。首轮 CI 结果出来后，需把 `FLOW_FIX_CHECKLIST.md` 的 `- [ ]` 改为 `- [x]`、状态总览改回「0 项在册待做」，并补一轮「17 组 verdict 同轮全绿 + 生产零扰动」证据，才算 B0 收口。
