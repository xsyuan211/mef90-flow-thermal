# 施工排量–页岩可压性匹配：研究基础、代码接入可行性与技术路线

**日期**：2026-11-10
**基线版本**：立项基线 `HEAD=abd5821`（CI #131 收口）；**B0 落地与验收**见第 11 节及 `FLOW_FIX_CHECKLIST.md`；**B1/B2 落地**见第 12 节；**CI #127 首轮结果与 L4 物理前提更正**见第 13 节；**CI #128 复验收口**见 13.6；**页岩油文献参数重算**见第 14 节；**B3 起裂准入立项**见第 15 节；**CI #132 首轮否定结论与 N5 机制分类更正**见 15.9、**CI #133 复验收口与证据等级**见 15.8 / 15.11
**依据**：公开文献与软件发布信息（2026-11-10 检索及二次核对）、项目源码静态阅读、CI #111–#125 真实 diagnostics artifact、CI #127 / #128 完整 job logs 文本（8,483,138 / 8,487,437 字符）与 `gb2_*` 实测指标、纯标准库离线标度计算
**约束**：零本地编译、零本地求解器运行；B0–B2 全程**零 Fortran 改动**。本文的离线标度/网格元数据不是求解器结果，CI 实测数值另行标注。
**状态**：阶段 A 研究设计稿 + **B0 已由 CI #125 验收** + **B1 已由 CI #127 验收** + **B2 已由 CI #128 验收（`GATE415_B2_RESULT: PASS`，19 组判据同轮全绿）** + **页岩油文献参数重算已由 CI #130 验收（20 组判据同轮全绿）** + **B3 已立项落地，CI #132 首轮给出否定结论 `B3_REGIME = UNIFORM_SOFTENING` / `B3_INITIATION_ACCESS = NOT_ESTABLISHED`（第 15.9 节）**；N5 已按实测重写为机制分类判据，**CI #133 复验 `success`**（`GATE415_B3_RESULT: PASS`，21 组 verdict 同轮全绿；证据等级为逻辑蕴含级，见 15.11）。B0/B1/B2 **全部收口**。第 5.2 节 R4 与 R5 的物理前提已被文献与 CI 实测推翻并更正（R4′ / R5′）：起裂前存储主导窗口不适用传播段标度，详见第 12.1、13.2 节。B3 改变的是**物理前提**：从「起裂前存储主导」迈到「刚度显著退化的起裂窗口」，但**仍不等于裂缝传播段**。引用线索见第 9 节，投稿前需二次核对卷期与 DOI。

---

## 0. 结论先行

| 问题 | 结论 |
|---|---|
| 现有研究是否已回答"排量–可压性匹配"？ | **没有形成统一结论**。文献对"Q↑ 是否增加缝网复杂度"给出相反答案，分歧源于各自固定了不同的约束（相同注入时间 vs 相同注入体积、是否含孔眼摩阻/限流、是否显式 DFN、是否滤失主导）。这本身就是可研究空白。 |
| 接入现有代码是否可行？ | **可行，且比初判更省力**。阶段 B 原则上**不需要新增 Fortran**：`flow.source` + `source_scaling` 已能给出常排量/线性变排量，`FractureToughness` / `internalLength` / `toughnessAnisotropyMatrix` / `damage_type` 已覆盖可压性，`Gate 4.5` 能量分解与 EXO 输出已覆盖目标量。需要的是**新 deck + 新判据 + CI 接线**。 |
| 最大的技术障碍是什么？ | **不是代码，是相似律**：现 deck 的 `k_min=5e-11 m²` 让 `T=c_v t/H²≈1`（1 s 跑完固结），而现场页岩 1 h 注入 `T≈3e-8~1.6e-5`。要复现机制必须把 `k` 降 **6–8 个量级**并联合调整 `S`，否则压力瞬间扩散、注入压力上不去。 |
| 平台能证明什么 / 不能证明什么？ | 能证明**耗散分区（黏性↔韧性）与孔弹时间尺度耦合下的排量响应、能量分配、损伤形态**；**不能**独立复现滤失主导区、井筒/孔眼摩阻、多簇分流、支撑剂、DFN、地应力张量、现场产能递减。 |
| 建议下一步 | **B0（CI #125）/ B1（CI #127）/ B2（CI #128）已全部验收收口**，`GATE415_B2_RESULT: PASS`、19 组判据同轮全绿，17 条非 B2 verdict 与 #127 逐字节 IDENTICAL。变排量离散锚 `Δp_10=5.5e4 Pa = 0.55×B0`、3D `β=1e5 Pa/s`（偏差 `1.08e-13`）、起裂前各向异性不可辨识 `gb2_surface_dev_aniso=4.064640e-03` 均由求解器确认（第 13.6 节）。**仍待用户确认**：是否提供区块实测参数、是否把本平台推进到真正的裂缝传播窗口（那需要降低 `G_c` / 提高排量 / 延长时间，可能触发新的求解器收敛问题）。 |

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

### B1/B2 首轮 CI（#127）后的第六处更正（2026-11-10 同日）

| # | 本文 / B2 判据原表述 | CI #127 实测后的事实 | 证据 |
|---|---|---|---|
| 6 | 第 5.2 节 R5 与 `check_gate415_b2_extension.py` 旧 L4：`toughnessAnisotropyMatrix≠I` 必须让表面能变化 >5%，否则判「矩阵未接线」 | **错，是传播段前提**。起裂前 `Damage_max≈1.48e-4` 且 α 近似空间均匀，A 只经 `C2·∇α·∇α` 进入表面能，被 `C1·w(α)` 淹没。实测 `surface_ratio=1.004065`（**0.4065%**）、`damage_dev=1.247%`、`beta_dev=0` ⇒ 各向异性**已接线、但在能量与压力曲线上不可辨识**，与 B1 的 R4′ 分层结论一致 | #127 `gb2_*` 指标；真实能量行 `surface: 6.49012E-04 → 6.51650E-04`；`m_MEF90_DefMechAssembly.F90:2340-2341,2355` |

> 第 6 条与 #124 的 `pmax` 事件同属「**先验算物理阶段，再定判据阈值**」：#124 把局部过冲误当体相斜率，#127 把传播段标度误当起裂前签名，两次都是判据的**物理前提**错，而非 Fortran / deck / EXO 缺陷。修复一律只改判据 + 补可证伪负例，**不改 deck 去凑判据**。

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
- 黏性主导 $p_{net}\sim(E'^2\mu'/t)^{1/3}$ 与 Q **无关**；韧性主导 $p_{net}\propto Q^{-1/3}$ 反而下降 —— ⚠️ **以上两条都是裂缝「传播段」标度**（Detournay 2016 / Savitski & Detournay 2002），**不适用于 B0/B1 的起裂前存储主导窗口**；在该窗口 $p_w(t)\approx p_0+(Q_0/U)t$，斜率只由注入系统柔度决定（Lecampion et al. 2017）。详见第 12.1 节的更正；
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
| **B0 校准冒烟（CI #125 验收通过）** | 表 C 的 $(k,S)$ | 2 个 2D deck（4 行 + 1 行差异，11 步）；V1–V8 + CI 延迟硬裁决 | 冻结 `[1/s]`、`Q_total=1.5e-4 m²/s`、`T=1.57e-7`、`β=1e5 Pa/s`、`S_eff=1.5e-10` | 体积守恒闭合、零源冻结、注入线性、KSP 每步收敛、生产零扰动 | 否 |
| **B1 主实验（CI #127 验收通过）** | B0 + Lecampion (2017) 起裂前框架 | **6 个单行差异 deck**：$Q\times S\times G_c\times k\times\mu$ 单因子扫描；`A1–A8` 判据；CI 硬闸门 | 参数**可辨识性分层**图版：β 只可辨识 $(f,S)$，表面能可辨识 $G_c,Q$；能量分配曲线 | A1–A8 全 PASS、B0/生产零扰动 | 否 |
| **B2 扩展（CI #128 验收通过，L1–L8 全绿）** | B1 判据设计 | `linear` 变排量（$f(t)=f_0t$、离散 $t^2$ 锚）、`toughnessAnisotropyMatrix: 2,0.5,0` 层理各向异性、3D tetra4 注入；`L1–L8` 判据 | 变排量解析锚（**#127 实测 $\Delta p_{10}=5.5\times10^4$ Pa $=0.55\times$ B0**）、起裂前各向异性**不可辨识**零位命题、3D 一致性（**#127 实测 $\beta_{3D}=10^5$ Pa/s**） | L1–L8 全 PASS、B0/B1/生产零扰动 | 否（`file` 程序需新增纯 Python CDF2 writer，本轮不做） |
| **B3 起裂准入（CI #133 复验收口，证据为逻辑蕴含级）** | B0/B1/B2 的 AT2 残差契约 + CI #125/#127/#132 真实读数 | **4 个单行差异 deck**（仅 `cs0001.fracturetoughness`，$G_c$ 降 1–4 个数量级）+ 复用 B0 注入 / B1 `gc` 两个既有 EXO 作为高档锚；`N1–N8` 判据 + `gb3_rc` 延迟硬裁决 | 起裂门槛 $G_c^*=9.650647$ N/m（以 $a(\alpha)=1/2$ 定义）、实测夹逼区间 $[2.7,27]$ N/m、$\alpha\sim1/G_c$ 幂律；**损伤机制分类实测为全域均匀软化（$d_{max}/d_{mean}=1.000004$），不是裂纹** | #132 真实 solver：N1–N4 / N6–N8 全 PASS，旧 N5「必须局域化」前提错误；$\beta=f/S$ 全程逐位不变；产出 `B3_REGIME=UNIFORM_SOFTENING`、`B3_INITIATION_ACCESS=NOT_ESTABLISHED`。**CI #133 复验 `success`**：`GATE415_B3_RESULT: PASS`、21 组 verdict 同轮全绿（逻辑蕴含级，见 15.11） | 否 |
| **C 现场与数据驱动** | 施工曲线、DFIT/G 函数、微地震、产剖、88 口井 | 反演 $C_L$ 与可压性；把 $\mathcal{K},\mathcal{L},p_{perf}/p_{net}$ 作为交互特征入产量模型；统计检验 | 排量优选表、论文验证段 | 现场数据闭合 | 否（本平台外） |

### 5.2 Gate 4.15 判据设计（可证伪，沿用现有硬闸门范式）

| 判据 | 检查内容 | 数据来源 | 负例设计 |
|---|---|---|---|
| R1 注入体积守恒 | $\sum_i f_iV_i\Delta t$ 与 压力储存 + 边界流出 + 损伤增强储存 之差 < tol；$Q_{total}$ 与 B0 标定常数一致 | `terzaghi415_*.log` 的 Gate 4.5 行 + EXO `source`/pressure | 篡改一个源项数值 / 去掉一步 |
| R2 排量单调性 | Q 增大 ⇒ 峰值 $p$、注入功 $W$、总表面能不减；零源对照 `W=0` | 多 deck 日志 | 打乱 Q 与 p 的对应、把零源 deck 的 W 置非零 |
| R3 可压性单调性 | $G_c$ 增大 ⇒ 起裂延迟、损伤总长/表面能占比下降（固定 Q） | EXO damage 场 + 能量行 | 反转 $G_c$ 顺序、令损伤场与对照逐位相同 |
| ~~R4 机制分区标度~~ **（已撤回，见 12.1）** | ~~双对数斜率：黏性主导段 $\partial\ln p_{net}/\partial\ln Q\approx0$；提高 $G_c$ 或 $\mu$ 后斜率向 $-1/3$ 移动~~ ⇒ **改为 R4′ 参数可辨识性分层**：起裂前 $\beta=\partial p_{bulk}/\partial t=f/S$ 只由 $(f,S)$ 决定；$G_c,k,\mu$ 对 $\beta$ 的影响 <6%，但 $G_c$ 与 $Q$ 对相场表面能有可测影响 | 6 个 B1 deck 的 `pbulk` 斜率 + 日志 `======= Total: ... surface:` 行 | 让 $G_c$ 改动压力斜率（传播段行为）、让 $G_c$ 表面能逐位不变（C1/C2 未接线）、把 $\beta$ 写成 $Q$ 无关 |
| ~~R5 各向异性效应（首版）~~ **（#127 后更正为 R5′）** | ~~`A≠I` ⇒ 损伤路径偏向低韧性方向，表面能显著变化~~ ⇒ **R5′ 起裂前各向异性不可辨识零位命题**：A 只进 `C2·∇α²`，必须「已接线」（表面能非位等价、damage 非位等价、bag/源码契约齐全）但「不可辨识」（$\lvert\Delta E_{surf}\rvert$ 有 0.01%–5% 双侧带、且弱于 damage 响应、$\beta$ 不变） | EXO damage + 日志能量行 + 真实 `beta` | 表面能逐位相同（A 没进 C2）、A 泄漏到 C1 通道（>5% 或超过 damage 响应）、矩阵改回单位阵、$\beta$ 被扰动 |
| R6 网格/长度尺度 | $h/\ell\le0.5$；$\ell/H$ 固定的两组给出同量级 $\mathcal{K}$ 指标 | 网格元数据 + 静态扫描 | 用未收敛网格蒙混 |
| R7 生产零扰动 | 15 组既有 verdict 与 CI #121 逐字节一致；`out_main.gen` 物理数值逐位不变（仅 QA changeset 字段允许变） | artifact 对比 | 伪造 verdict |
| R8 静态反回归 + CI 裁决 | 源码/量纲注释、`flow.source` 注册、`source_scaling` 分支、CI 接线、`s415_rc` 参与延迟 fail、`timeout` 与 KSP monitor 存在 | 源码静态扫描 | 删任一接线（复用 #119/#120 两类失败） |

判据脚本：`check_gate415_rate_fracability.py`（B0，V1–V8）、`check_gate415_b1_identifiability.py`（B1，A1–A8）、`check_gate415_b2_extension.py`（B2，L1–L8），均自带 `--selftest`（B0 1+25、B1 1+20、B2 **#127 后 1+24** 变异）并复用 `check_gate414_dim3.py --gate` 的模板与延迟 `exit 1` 结构。**实装对应关系**：R1/R4′/R8 → B0 的 V1–V8 前置子集 + B1 的 A2/A3/A4/A7；R2 排量单调性 → A6；R3 可压性 → A5（$G_c$ 在表面能上可辨识，在压力曲线上不可辨识）；**R5′ 各向异性 → B2 的 L4/L5**（L4 在 #127 后由「必须大变」改为「已接线但不可辨识」的双侧带，见第 13 节）；R6 网格尺度**本轮不实装**（固定 356 节点 tri3 网格、不做 $h/\ell$ 收敛扫描，避免把网格敏感性伪装成机制分区）。第 12 节给出 B1/B2 的完整交付与物理边界，第 13 节给出首轮 CI 的订正。

### 5.3 CI 与零扰动约束

- 新 deck 一律**单行差异**契约（B0 为「校准 deck 对基准 4 行白名单」+「注入 deck 对校准 deck 1 行」两级），由判据的纯标准库 `code_lines()` / `code_diff()` 校验（CI 无 PyYAML，见第 0 节补充更正 5）；
- 所有新运行放在 `Run Terzaghi` 后段，生产 12 个 2D deck 与 3D 诊断运行**保持原配置不动**；
- 捕获 `s415_rc`（mpirun 退出码）并接入延迟 hard fail，避免 #120 的「rc=124 只 echo」漏洞；
- artifact 增 `gate415_verdict.log` 与新 deck 的 `.log`/`.gen`；
- 零本地编译、零本地求解器：所有验证靠 `git diff --check`、`ast`、PyYAML（**仅本机离线复核，不进 CI**）、10 个 run block `bash -n`、判据 selftest、真实 artifact 离线回放。

### 5.4 决策点

1. **B0/B1/B2 立项与验收（已全部收口）** —— B0 由 CI #125 验收；用户以「通过文献调研实现B1和B2」立项，B1 由 CI #127 验收（`GATE415_B1_RESULT: PASS`），B2 首轮 #127 仅 L4 因物理前提错误 FAIL，L4 重写后由 **CI #128 验收**（`GATE415_B2_RESULT: PASS`，19 组判据同轮全绿）。清单 B1 / B2 均已勾 `[x]`，详见第 13.6 节；
2. **区块实测参数**（$E,\nu,K_{Ic},C_L,k,\Delta\sigma$、孔数/孔径/簇数、Q 范围）：**区块实测仍待确认；公开文献替代已完成**。第 14 节用 S01–S20 页岩油文献参数把表 A–C 重算为表 A′–C′，并给出可证伪结论：**不能产出可直接用的排量优选表**（起裂前 $\beta=f/S$ 对 Q 线性无阈值、$\mathcal{K}/\mathcal{C}$ 属传播段且不可辨识、$T$ 文献包络跨 6.48 个数量级）。表 A′–C′ 的用途是数值平台相似区间与阶段 B3 边界条件；要收窄到区块级仍需 DFIT/mini-frac 实测；
3. **是否推进到真正的裂缝传播窗口 → 已由用户以「立项b3」立项，见第 15 节**：当前 `Damage_max=1.48e-4≪1`，B0/B1/B2 结论严格限定在起裂前。B3 用**零 Fortran 的单因子韧性阶梯**把损伤从 $10^{-4}$ 推过 $a(\alpha)=1/2$ 的刚度半衰门槛，并分类损伤机制；CI #132 实测末档 $\alpha\to1$ 但 $d_{max}/d_{mean}=1.000004$，即**全域均匀软化、不是裂纹**（15.9）；`-5` 残差平台按既有纪律显式报告、不作硬判。⚠️ B3 仍**不是**裂缝传播段：Santillán (2018) 排量阈值、Detournay (2016) M–K 顶点、$\mathcal{K}/\mathcal{C}$ 分区与排量优选表在本阶段仍**没有定义**，需另行立项（预裂纹 / 非均匀载荷 / 传播段求解器）；
4. 阶段 B2 之后是否允许新增 Fortran（等效 DFN 通道、显式 Carter 滤失、井筒限流，或 `source_scaling: file` 所需的 CDF2 writer）—— **仍待确认**，会突破「零代码改动」的优势，需单独评估。

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
4. Santillán, D., Juanes, R., Cueto-Felgueroso, L. (2018). *Phase Field Model of Hydraulic Fracturing in Poroelastic Media: Fracture Propagation, Arrest, and Branching Under Fluid Injection and Extraction*. **JGR: Solid Earth 123(3)**, 2127–2155. DOI `10.1002/2017JB014740`.（三类排量模式、$L_Q=Q/(Ek^\*)$ 无量纲塌缩；属**传播段**。）
4b. **Lecampion, B., Desroches, J., Jeffrey, R.G., Bunger, A.P. (2017).** *Experiments versus theory for the initiation and propagation of radial hydraulic fractures in low-permeability materials*. **JGR: Solid Earth**, DOI `10.1002/2016JB013183`.（**B1/B2 的关键锚**，2026-11-10 已按全文页二次核对：「Prior to fracture initiation, the flow rate entering the fracture Qin is negligible, which for a constant injection rate Q0 results in a linear pressurization rate. The injection system compliance U can thus be directly obtained from the slope of the pressurization phase ... knowing Q0.」同时明确区分 initiation pressure 与 breakdown pressure。另见 Lecampion et al. (2015), 13th ISRM Congress, ISRM-13CONGRESS-2015-285。）
4c. Step-rate test (SRT)：等时递增排量、记录每级稳定压力，破裂前 pressure-rate 线性、破裂后斜率下降；B2 的 `source_scaling: linear` 是**连续斜坡**而非阶梯，故用离散 $t^2$ 解析锚而非 SRT 斜率折点做判据（工程类比，不作定量等价）。
5. Detournay, E. (2016). *Mechanics of hydraulic fractures*. **Annual Review of Fluid Mechanics 48**, 311–339.（M/K × 存储/滤失四顶点，**传播段**标度，不适用于 B0/B1 压力窗口。）
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

## 11. 阶段 B0 落地记录（2026-11-10 首版；验收状态见 11.5 末）

> 本节按**首版落地当时**的时态保留，11.5 末尾给出 CI #125 的验收结论。不要把下面的「尚未取得 CI 实跑结果」当作当前状态。

**结论先行**：用户以「推进」默认立项 Gate 4.15，B0 已完成**全部代码与 CI 接线、零 Fortran 改动、零本地编译、零本地求解器**；判据离线自检 1 正例 + 20 负例全 PASS。**当时尚未取得 CI 实跑结果，因此以下 B0 数值是离线标度与解析预测，不是求解器验证结论**（对应第 4.3 节风险 R7）。

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

> **状态（2026-11-10 更新）**：B0 已由 **CI #125 验收通过**，`GATE415_RESULT: PASS`，17 组判据同轮全绿，原有 16 组 verdict 与 #124 逐字节 IDENTICAL；清单 4.15 B0 已勾 `[x]`。本节 11.2 的离线复算中 `β=1.000000e+05 Pa/s`、`S_eff=1.500000e-10`、`W_last=30.1605` 已由真实求解器确认；其余相似量（$T,c_v,Q_{total}$）仍是离线标度。B1/B2 见第 12 节。

---

## 12. 阶段 B1 / B2 落地记录（2026-11-10 首版；CI #127 结果见第 13 节）

**结论先行**：用户要求「通过文献调研实现B1和B2」。文献复核**推翻了路线稿 R4 的物理前提**，因此 B1 重新定位为「起裂前存储主导段的**参数可辨识性分层**」实验，而不是「传播段机制分区」实验。B1/B2 已完成全部 deck / 判据 / CI 接线，**零 Fortran 改动、零本地编译、零本地求解器**；离线自检与 #125 真实 artifact 回放全 PASS。本节保留首版落地时态与解析预测；**首轮 CI 实跑结论、R5′ 的再次更正与 L4 重写见第 13 节**。

### 12.1 文献复核：R4 为什么必须撤回

| 项 | 路线稿首版（第 5.2 节） | 复核后的事实 |
|---|---|---|
| 适用阶段 | 默认当前 deck 能给出「黏性主导 / 韧性主导」的 $p_{net}$–$Q$ 标度 | **B0 实测 `Damage_max=1.48e-4 ≪ 1`**：没有裂缝扩展，是**起裂前**窗口 |
| $p_w(t)$ | 用 $\partial\ln p_{net}/\partial\ln Q\approx0$ 或 $-1/3$ 检验 | Lecampion et al. (2017, doi:`10.1002/2016JB013183`)：起裂前 $Q_{in}\approx0$，恒排量下**线性加压**，斜率 $=Q_0/U$，$U$ 为注入系统柔度；岩石力学参数在这条压力曲线上**不可辨识**；起裂压力 ≠ 破裂（峰值）压力 |
| Detournay (2016) $\mathcal{K}$ / M–K 顶点 | 当作 B1 判据锚 | 属**裂缝传播段**标度；Santillán et al. (2018) 的三类排量模式同样属传播段 |
| 本项目对应量 | — | `S`（storativity）扮演 $U/Q_0$，$\beta=f/S$；B0 CI #125 实测 $\beta=1.000000\times10^5$ Pa/s、$S_{eff}=f/\beta=1.5\times10^{-10}$ |

**更正后的可证伪命题（R4′，可辨识性分层）**：

```text
Level 1  压力曲线 p_bulk(t)：beta = f/S
         可辨识 (f, S)；在 T = c_v t / H^2 << 1 下对 Gc, k, mu 不敏感
Level 2  相场表面能 / damage：
         C1 = Gc/(ell*cw)/4, C2 = Gc*ell/(4*cw)*A   (A = toughnessAnisotropyMatrix)
         m_MEF90_DefMechAssembly.F90:2340-2341, 2355
         对 Gc、Q、A 有可测响应
```

这既符合 Lecampion (2017)，又保留「可压性」的研究价值：**$G_c$ 并非不可研究，只是不能从井筒压力曲线上研究**。同时把 Santillán (2018) / Detournay (2016) 明确划入传播段，避免把不同阶段的结论混用。

### 12.2 为什么可以用表面能量做判据（且必须诚实标注）

从 CI #125 真实日志读到每步一行（共 11 行）：

```text
  [solver row] "======= Total: elastic energy: 1.31015E+01 work: 2.61764E+01 \
                cohesive: 0.00000E+00 surface: 6.49012E-04 total: -1.30742E+01"
```

- 零源对照 `surface: 1.20334E-06 → 1.20755E-06`（几乎冻结）；常排量注入 `1.20334E-06 → 6.49012E-04`（增长 539×）。
- `MEF90DefMechSurfaceEnergy`（`m_MEF90_DefMechAssembly.F90:2292-2372`）在 :2355 用 `C1*w(α) + (C2*∇α)·∇α`，其中 :2340-2341 的 `C1 ∝ G_c/ℓ`、`C2 ∝ G_c·ℓ·A`。所以改 `fracturetoughness` 或 `toughnessAnisotropyMatrix` **必然**改变该积分量。
- ⚠️ **诚实边界**：B0 注入日志含 5 条 `[ERROR]: damage field SNESSolve failed with SNESConvergedReason -5`（既有已知残差平台，CI #84；生产因此用 `-disp_snes_type ksponly`）。能量行是**当前场的确定性积分**，本轮只用作**相对 / 可证伪**观测量（「参数改动是否产生可测差」），**绝不**把它表述为已收敛的断裂力学解，也不据此报裂缝长度、缝网面积或起裂压力绝对值。

### 12.3 B1 交付：Q × S × Gc × k × μ 单因子矩阵

| 文件 | 相对 `terzaghi415_inject.yaml` 的唯一代码行差异 | 物理用途 |
|---|---|---|
| `terzaghi415_b1_qhalf.yaml` | `source: 1.5e-5 → 7.5e-6` | $Q=0.5\times$ |
| `terzaghi415_b1_qdouble.yaml` | `source: 1.5e-5 → 3.0e-5` | $Q=2\times$ |
| `terzaghi415_b1_stor.yaml` | `storativity: 1.50e-10 → 3.00e-10` | $S=2\times \Rightarrow \beta$ 减半 |
| `terzaghi415_b1_gc.yaml` | `fracturetoughness: 27000.0 → 13500.0` | $G_c=0.5\times$（可压性） |
| `terzaghi415_b1_perm.yaml` | `perm_min: 2.35e-18 → 2.35e-17` | $k\times10$，$T=1.57\times10^{-6}$ 仍 $\ll1$ |
| `terzaghi415_b1_visc.yaml` | `viscosity: 1.0e-03 → 2.0e-03` | $\mu\times2$，$c_v$ 减半 |

每个 deck 头部都写入 Lecampion (2017) 的文献锚与「传播段标度不适用」声明；`code_lines()` 校验「恰好 1 行代码差异」。

`check_gate415_b1_identifiability.py` 的 **A1–A8**：

| 判据 | 检查 | 负例（节选） |
|---|---|---|
| A1 | 8 个 run × 11 记录 / 11 行 Gate 3.2、4.5、能量；EXO `pmax` ↔ 求解器逐位交叉；356 节点 tri3；无 fatal | 截断 EXO、fatal、陈旧 mesh、能量行缺失 |
| A2 | $\beta(0.5/1/2\times Q)$ 比值 0.5/1/2，$S_{eff}=f/\beta$ 恒为 $S_{num}$ | 把 $2\times Q$ 偷偷跑成 $1\times Q$ |
| A3 | $2\times S \Rightarrow \beta$ 减半、$S_{eff}$ 恢复 $2S_{num}$，显式拒绝「β 只依赖 f」 | 柔度臂仍按 $\beta_0$ 爬升 |
| A4 | **Lecampion 数值验证**：$G_c/k/\mu$ 三臂的 $\beta$ 偏差 <6% | 让 $G_c$ 或 $k$ 改动压力斜率 |
| A5 | $G_c$ 使表面能显著变化（>5%）而 $k/\mu$ 不使表面能变化 >60%；A4+A5 的**解耦**才是可证伪命题 | 表面能逐位不变、$k$ 臂污染表面能 |
| A6 | $Q=0/0.5/1/2\times$ 的 $W_{last}$、表面能、$damage_{max}$ 均单调不减；零源 $W\equiv0$ | 零源非零功、$2\times Q$ 功不单调 |
| A7 | 6 个单行 deck 契约、源码 C1/C2 与 `[1/s]` 注册、CI run/verdict/延迟 fail、≥8 次 EXO 重建 | 第二处代码改动、删 EXO 重建、删 `gb1_rc` 裁决 |
| A8 | B0 仍是 CI #125 的 PASS 与 `β=1e5 / S_eff=1.5e-10`；生产固结、$W=0$、无 B1 泄漏 | B0 变 FAIL、B1 deck 名泄漏进生产日志 |

### 12.4 B2 交付：变排量 / 层理各向异性 / 3D 一致性

| 通道 | deck / run | 解析锚或契约 |
|---|---|---|
| 变排量 | `terzaghi415_b2_linear.yaml`：`source_scaling: constant → linear` | `MEF90FlowSetTransients`（`m_MEF90_FlowAssembly.F90:873-875`）：`SetSourceCst` 后 `VecScale(source,timeVal)`，即 $f(t)=f_0t$。transient 模式下 `vDef.F90:661` 的第二次调用被 `If (.not. flowTransient)` 保护，唯一时间值是 `time(step)`。**BE 离散精确解** $\Delta p_n=(f_0/S)\,dt^2\,n(n+1)/2$，$n=10\Rightarrow5.5\times10^4$ Pa $=0.55\times$ B0。L2 拟合 $\Delta p$ 对 $n(n+1)/2$，并要求线性模型残差 $>3\times$ 二次残差 |
| 层理各向异性 | `terzaghi415_b2_aniso.yaml`：inject + 1 行 `toughnessAnisotropyMatrix: 2,0.5,0` | `MatS2D=(XX,YY,XY)`、identity `=(1,1,0)`（`m_MEF90_LinAlg.F90:50-54,107`）。逗号后**无空格**（CI #85 RealArray 坑）。⚠️ **首版「两个主方向同时偏离 ⇒ 表面能必变（>5%）」已被 CI #127 推翻**：实测只变 0.4065%，见第 13.2 节；现行 L4 改为「已接线但不可辨识」双侧带 + $\beta$ 不变（<6%） |
| 3D 一致性 | 复用 `terzaghi415_inject.yaml` + `-prefix data/Beam1x10-tetra4` | 维度由 mesh prefix 分派而非 YAML 键（Gate 4.14）。复用 #121 的 `-np 1` + `-flow_pc_type bjacobi -flow_sub_pc_type lu`、`timeout 1800`，避免 #120 超时。L6 要求 3D $\beta$ 与 2D 差 <12%、6 个 nodal vars 含 `Displacement_Z`、无超时 |

`check_gate415_b2_extension.py` 的 **L1–L8**：数据完整性（2D+3D）/ 离散 $t^2$ 锚 / 变排量注入功 / 各向异性解耦 / deck-bag-源码静态契约 / 3D 分派与 direct-LU 配置 / CI 接线 / B0+B1+生产零扰动。

### 12.5 CI 接线与生命周期护栏

- `build.yml` 新增 **9 个诊断 run**（B1 6 × 2D、B2 2 × 2D + 1 × 3D），全部 11 步、保持 12 个生产 2D deck 与 Gate 4.14 的 3D 固结运行**原配置不动**；
- **每个 11 步 run 前显式 `rm -f data/Beam1x10-tri3_out.gen`**（共 11 次，含生产 run 原有 1 次），3D run 前 `rm -f data/Beam1x10-tetra4_out.gen`；这是 CI #124 的 EXO append 教训，由 A7/L7 静态计数反回归保护；
- 3D run 后从 `/tmp/out_dim3.gen` 恢复 Gate 4.14 的 3D 结果 mesh，2D 段末仍恢复 `/tmp/out_main.gen`；
- 新增 `gb1_rc` / `gb2_rc` 两个 verdict 与 **9 个 mpirun 退出码**全部接入延迟 `exit 1`（CI #120 教训：超时 rc=124 只 echo 会造成静默假绿）；
- artifact 增加 `gate415_b1_verdict.log`、`gate415_b2_verdict.log` 与 `/tmp/gate415b1`、`/tmp/gate415b2` 两个目录；判据以 `--gate <dir>` 直接消费，便于离线回放；
- CI 预算：#125 `Run Terzaghi` 1646 s，新增 8 个快速 2D run 与 1 个 3D direct-LU run，需实测确认总时长。

### 12.6 提交前离线复核（零编译 / 零求解器，全部通过）

| 检查 | 结果 |
|---|---|
| `git diff --check` | OK |
| 两个判据 `ast.parse` | OK |
| PyYAML 解析 workflow + 10 份 415 deck | OK（`jobs=['compile']`，8 个新 deck `numstep=11`） |
| `Run Terzaghi` 整个 `run` 块 `bash -n`（1384 行） | OK；11 次 tri3 EXO 重建、1 次 tetra4 重建、全部 rc/verdict token 齐全 |
| B1 `--selftest` | **PASS**：1 正例 + 20 变异（含 1 不变性：damage SNES `-5` 不得误判） |
| B2 `--selftest` | 首版 **PASS**：1 正例 + 19 变异（含 1 不变性）；**#127 后重跑 PASS：1 正例 + 24 变异（含 3 不变性）**，正例的表面能 / damage 已锚定 #127 实测 `surface_ratio=1.004065`、`damage_dev=1.247%` |
| **#125 真实 artifact 混合回放** | B1 A1–A8、B2 L1–L8 全 PASS。真实读数：`pbulk[0]=1e5`、`pbulk[-1]=199999.9999999968`、11 行能量分解、`surface_last=6.49012e-4`、`damage_max=1.4803e-4`、tri3 356 节点 / tetra4 2635 节点 `num_dim=3` |
| B0 判据零扰动复核（#125 artifact + 新 workflow） | `GATE415_RESULT: PASS`，V1–V8 与 #125 逐位一致（`β=1.000000e+05`、`S_eff=1.500000e-10`） |
| deck 代码行契约 | 7 个 B1/B2 deck 各恰好 1 行差异；aniso 恰好 +1 行且 key 为 `toughnessAnisotropyMatrix` |

### 12.7 状态纪律与能力边界

- 本节首版 B1/B2 数值（$\beta$ 比值 0.5/1/2、$S_{eff}=3\times10^{-10}$、$\Delta p_{10}=5.5\times10^4$ Pa、表面能比值）原为**离线解析预测**；**CI #127 已把其中 $\Delta p_{10}=5.5\times10^4$ Pa、$\Delta p_{10}/\mathrm{B0}=0.55$、3D $\beta=10^5$ Pa/s 变成求解器实测**，逐项复核见第 13.1 节。判据阈值若与真实数值有系统性偏差，**应改判据的解析锚或容差并记录理由，不得改 deck 去凑解析解**（沿用 #124→#125 的处置原则，#127 的 L4 即按此办理）。
- 所有结论严格限定在**起裂前、存储主导、无裂缝扩展**窗口；`Damage_max\sim10^{-4}`，不报缝长 / 缝网面积 / 起裂压力绝对值 / 现场排量优选表。**层理各向异性「控制裂缝轨迹」的文献结论同样只适用于传播段，本窗口不主张验证**（第 13.2 节）。
- 不主张验证 Santillán (2018) 排量阈值、Detournay (2016) M–K 顶点或任何传播段标度；若要进入该窗口，见第 5.4 节决策点 3，需先解决 damage SNES `-5` 并重新标定 $G_c$ / 排量 / 时间。
- **收口条件已满足**：清单 4.15 B1（CI #127）与 B2（CI #128）均已勾 `[x]`；CI #128 给出 `GATE415_B2_RESULT: PASS`、19 组 verdict 同轮全绿、18 条 `GATE*_RESULT` 全 PASS、17 条非 B2 verdict 与 #127 逐字节 IDENTICAL、0 条 `=== FAIL`（第 13.6 节）。

---

## 13. CI #127 首轮、L4 物理前提更正与 CI #128 收口（2026-11-10）

**结论先行**：`f07cf39` 的 CI #127 `conclusion=failure`，**编译成功、18 组 verdict 中 18 组物理/工程判据全绿**，唯一失败是 `GATE415_B2_RESULT: FAIL (L4)`，最终由延迟 `gb2_rc=1` 触发 `exit 1`。**这是判据的物理前提错误，不是 Fortran / deck / EXO / 3D 求解器缺陷**：旧 L4 要求 `toughnessAnisotropyMatrix` 让表面能变化 >5%，那描述的是**裂缝传播段**，而 B0/B1 已确认当前是**起裂前存储主导窗口**。按 #127 实测把 L4 重写为「已接线但不可辨识」的可证伪零位命题（零 Fortran / 零 deck 改动）后，`e6cce4b` 的 **CI #128 `success`**，`GATE415_B2_RESULT: PASS`，**19 组判据同轮全绿，B0/B1/B2 全部收口**（13.6 节）。

### 13.1 首轮实测（全部来自 CI #127 job logs 文本，非离线预测）

| 项 | 实测值 | 判据 |
|---|---|---|
| `Build MEF90 + m_Flow + vDef` | `success`（2D/3D 双编译） | — |
| `Run Terzaghi` | 1429 s（#125 为 1646 s；9 个新 run 未拉长总时长） | — |
| 原有 17 组 verdict | 全 PASS | — |
| `GATE415_B1_RESULT` | **PASS**（A1–A8） | B1 已验收 |
| L1 | PASS：5 run × 11 干净记录、EXO `pmax` == 求解器、每 run 11 行能量分解、2D/3D mesh 元数据、`Displacement_Z` 活 | 数据完整性 |
| L2 | `quad_slope=1000`（解析锚 `1000`）、`intercept=4.177e-8`、二次残差 `0.00%`（线性模型残差 `13.64%`）、`dp10=55000 = 0.550×B0` | 变排量离散 $t^2$ 锚**被求解器精确证实** |
| L3 | 对照 `W≡0`；变排量 `W=0..23.28441` 单调；常排量 `W_last=30.16054` | 注入功 |
| L4 | `surface_ref=6.490120e-4`、`surface_ratio_aniso=1.004065`、`beta_dev=0`、`damage_dev=1.247384e-2` | **FAIL（唯一）** |
| L5 | `aniso_matrix = toughnessAnisotropyMatrix: 2,0.5,0`、bag dump / `MatS2D` / `C2` / `VecScale` 接线齐全 | 静态契约 |
| L6 | `beta_dim3=1.000000e+05`、`beta_dev_dim3=1.082467e-13`、`Displacement_Z` 活、direct-LU PC、无 timeout | 3D 一致性 |
| L7 / L8 | PASS：单行 deck、9 个 run + `gb2_rc` 延迟裁决 + 11 次 EXO 重建；B0 `β=1e5 / S_eff=1.5e-10`、B1 PASS、生产 51 步固结、`W=0`、无 B2/3D 泄漏 | 反回归 + 隔离 |
| 最终失败行 | `=== FAIL: Gate 4.15 B2 rejected variable rate / anisotropy / 3D (rc=1) ===` | 延迟 `gb2_rc=1` |

### 13.2 L4 为什么失败：物理阶段错配，不是装配 bug

`MEF90DefMechSurfaceEnergy`（`m_MEF90_DefMechAssembly.F90:2340-2341,2355`）的被积函数是

```fortran
C1 = Gc/internalLength * 0.25_Kr / ATModel%cw
C2 = (Gc*internalLength*0.25_Kr/ATModel%cw) * matPropSet%toughnessAnisotropyMatrix
integrand = C1*w(alpha) + ((C2*gradDamageGauss) .dotP. gradDamageGauss)
```

各向异性矩阵 $A$ **只乘在梯度项 $C_2|\nabla\alpha|^2$ 上**。起裂前 `Damage_max≈1.48e-4` 且损伤场近似空间均匀，$w(\alpha)\sim\alpha^2$ 项压倒 $|\nabla\alpha|^2$，因此 $A$ 对总表面能只能产生 O(0.1%) 的影响。CI #127 实测：

```text
surface: 6.49012E-04  (b0_inject, isotropic)   -- 真实日志第 697 条 Total 行
surface: 6.51650E-04  (aniso, A = 2,0.5,0)     -- 真实日志第 785 条 Total 行
ratio = 1.004064640  ->  ds = 0.4065%
damage_dev = 1.2474%     beta_dev = 0
ds / damage_dev = 0.326  <- 梯度通道占比的签名
```

**排除代码缺陷的四条独立证据**：① L5 已验证 bag dump 出现 `toughnessAnisotropyMatrix`、deck 恰好 +1 行 `2,0.5,0`、源码 `C2 = ... * A` 与 `VecScale` 均在位；② `damage_dev=1.247% > 1%`，各向异性求解确实发生，不是位等价；③ 表面能非位等价（0.4065% 远高于浮点噪声），$A$ 确实进入了 $C_2$ 积分；④ linear / 3D / B0 / B1 / 生产隔离全部 PASS。**唯一不成立的是「必须 >5%」这一条阈值命题**。

旧 L4 与页岩层理文献（层理平行方向 $K_{Ic}$ 最弱、$G_c$ 各向异性主导裂缝轨迹）对应的是**已有传播裂纹**的情形；把它套到起裂前存储主导窗口，与第 12.1 节撤回 R4 的理由完全同构。这是本项目第二次「判据物理前提错误」：#124 是 `pmax` 局部过冲误当体相斜率，#127 是传播段标度误当起裂前签名。

### 13.3 L4 重写：从「必须大变」到「已接线但不可辨识」

`check_gate415_b2_extension.py` 新增 `SURF_BAND = (1.0e-4, 5.0e-2)`，`check_l4` 改为五条同时成立的可证伪零位命题：

| # | 命题 | 阈值 / 关系 | #127 代入 | 被证伪的情形 |
|---|---|---|---|---|
| 1 | $A$ 确实进入 $C_2$ | $d_s=\lvert E_{surf}^{aniso}/E_{surf}^{iso}-1\rvert > 10^{-4}$ | `0.004065` ✓ | 表面能位等价（矩阵未接线） |
| 2 | $A$ 没泄漏到 $C_1 w(\alpha)$，仍是起裂前场 | $d_s < 5\times10^{-2}$ | `0.004065` ✓ | 矩阵错接到均匀项 / 已进入传播段 |
| 3 | $A$ 只作用于梯度通道 | $d_s < d_{dmg}$ | `0.004065 < 0.012474` ✓ | 能量响应强于 damage 响应 |
| 4 | 各向异性求解确实发生 | $d_{dmg} > 10^{-2}$ | `0.012474` ✓ | damage 场与 isotropic 逐位相同 |
| 5 | 压力曲线盲（与 B1 A4 同构） | $d_\beta < 6\times10^{-2}$ | `0` ✓ | $G_c$ / $A$ 污染 $\beta=f/S$ |

**五条全 PASS**。`d_s < d_{dmg}` 是新增的关键序关系：它把「各向异性只进梯度项」变成可证伪命题，而不是一句定性说明。`report()` 新增 `gb2_surface_dev_aniso` 指标便于下一轮 artifact 直接复核。

`--selftest` 的正例已锚定 #127 真实读数（`SURF_REF=6.49012e-4`、`SURF_ANISO=SURF_REF*1.004065`、`DMG_ANISO=DMG_REF*1.01247384`），新增 4 个精确负例（位等价 / $A$ 泄漏到 $C_1$ / 能量响应超过 damage / damage 位等价）与 2 个不变性（表面响应在带内漂移到 0.30% / 0.90% 仍 PASS），加上原有的 damage SNES `-5` 不变性：**1 正例 + 24 变异（含 3 不变性）全 PASS**。

### 13.4 提交前离线复核（零编译 / 零求解器）

| 检查 | 结果 |
|---|---|
| B2 `--selftest` | **PASS**（1 正例 + 24 变异，含 3 不变性） |
| #127 真实 `gb2_*` 指标代入新 L4 | 五条命题 **5/5 PASS**（`ds=0.004065`、`dd=0.012474`、`beta_dev=0`） |
| 真实日志格式交叉验证 | 剥离 GitHub 时间戳后，用判据自身 `ENERGY` 正则从 8.48 MB job logs 解析出 797 条 `Total` 能量行；定位 `6.49012E-04` / `6.51650E-04`，算出 `ratio=1.004064640`，与 `gb2_surface_ratio_aniso=1.004065` 一致 |
| 解析层无需重验 | CI #127 的 L1 已在真实 log/EXO 上 PASS（11 行能量分解、`Gate 3.2` 逐位交叉），本轮未改任何解析函数与 L1–L3 / L5–L8 |
| 零扰动 | 未改 Fortran、未改任何 deck、未改 `build.yml`、未改 B0/B1 判据；B2 只改 `check_l4` 与其常量/文档 |
| `ci127_joblogs.zip` | 8.5 MB 临时诊断文本，加入 `.gitignore` 的 `/ci*_joblogs.zip`，不进版本库 |

### 13.5 复验条件（已由 CI #128 全部满足）

1. ~~推送后下一轮 CI 需给出 `GATE415_B2_RESULT: PASS`，且原有 17 组 verdict + `GATE415_B1_RESULT` 保持逐字节一致~~ ⇒ **CI #128 已满足**（13.6 节）；
2. ~~确认 `gb2_surface_dev_aniso` 与 `gb2_damage_dev_aniso` 仍落在 $(10^{-4}, 5\times10^{-2})$ 带内且 $d_s<d_{dmg}$~~ ⇒ **已满足**：`4.064640e-03` 与 `1.247384e-02`；
3. ~~把 `FLOW_FIX_CHECKLIST.md` 的 4.15 B2 勾 `[x]`、进度段改为「19 组 verdict 同轮全绿」~~ ⇒ **已执行**；
4. **能力边界不变**：本节所有数值仍是起裂前、存储主导、无裂缝扩展窗口；0.4065% 的表面能变化**不能**被解读为「层理各向异性对裂缝轨迹无影响」，只能解读为「该窗口下这个观测量不可辨识」。要研究轨迹控制，必须走第 5.4 节决策点 3 进入传播段。

### 13.6 CI #128 复验验收记录（`e6cce4b`，`success`）

run `42772909839`，`Build + Terzaghi` 全 success。完整 job logs（压缩 2,243,135 B / 解压 8,487,437 字符，含 22 个 step 文本）离线复核：

| 项 | CI #128 实测 |
|---|---|
| `GATE415_B2_RESULT` | **PASS**；L1–L8 全绿 |
| L4 真实输出 | `surface dev 0.406% in band (0.01%..5%) and < damage dev 1.247%, beta dev 0.00%: anisotropy is wired but NOT identifiable before initiation` |
| `gb2_surface_dev_aniso` | `4.064640e-03`（与离线预测逐位一致）；`gb2_surface_ratio_aniso=1.004065`、`gb2_damage_dev_aniso=1.247384e-02` |
| 其余 17 项 `gb2_*` | 与 CI #127 **逐位相同**（`quad_slope=1e3`、`dp10=5.5e4`、`dp10_ratio=0.55`、`beta_dim3=1e5`、`beta_dev_dim3=1.082467e-13` 等）⇒ 只改 L4 判定逻辑，物理场零扰动 |
| 全部 verdict | **18 条 `GATE*_RESULT` 全 PASS**，其中 17 条非 B2 verdict 与 #127 **逐字节 IDENTICAL**；25 条 `=== ... PASS` 汇总行、**0 条 `=== FAIL`** |
| 体系规模 | **19 组判据同轮全绿**（18 条 `GATE*_RESULT` + `terzaghi_verdict` 汇总口径） |

与 #124→#125 的处置完全一致：**只改判据 + 补可证伪负例，不改 deck / Fortran / workflow 去凑判据**；一轮 FAIL 的判据在下一轮由同一物理场给出 PASS，其余所有判据连日志字节都不变。至此 B0（`β=f/S=1e5 Pa/s`、`S_eff=1.5e-10`）、B1（可辨识性分层）、B2（变排量 $t^2$ 锚 + 起裂前各向异性不可辨识 + 3D 一致性）全部由真实求解器验收。

---

## 14. 页岩油文献参数重算：表 A′–C′ 与「排量优选表」为什么不能交付（2026-11-10）

对应第 5.4 节决策点 2（原为「区块实测参数仍待确认，有则把表 A–C 重算成可直接用的排量优选表」）。本轮以**公开文献参数**替代首版占位假设，完成了表 A–C 的完整重算，但结论是：**不能交付现场排量优选表**。这与 #127 的 L4 是同一类问题——不是算错了，而是物理阶段不支持那个结论。

### 14.1 交付物与零扰动边界

| 文件 | 作用 |
|---|---|
| `vDef/data/shale_oil_params.json` | S01–S20 文献来源、25 个参数的 low/rep/high 区间、单位换算、诚实边界、`window`（首版设定 / 文献映射 / 已发布表值）、`formulas`（指数与前因子契约）、`numerical`（B0 冻结常数） |
| `vDef/data/check_shale_oil_params.py` | M1–M8 硬判据；`--md` 输出本节表 A′–C′；`--selftest` 1 正例 + 33 变异 + 3 不变性；`--json` 指标导出 |
| `.github/workflows/build.yml` | 新增 `gso_rc` 延迟裁决与 `/tmp/shale_oil_params_verdict.log` artifact |

**零扰动**：不改任何 Fortran、不改任何 deck、不改任何生产 run、不新增求解器运行。判据只读 JSON / 路线稿 / B0 deck / workflow 文本。B0 冻结常数（`k_num=2.35e-18`、`S_num=1.50e-10`、`β=1e5 Pa/s`、`Q_total=1.5e-4 m²/s`）逐位不变。

### 14.2 重算脚本的可信度锚：逐位复现首版表 A–C

M3 用**独立于 JSON 契约**的参考公式重算首版四类页岩的全部已发布数字（$\sigma_c$、$T$、$C_L$、$\mathcal{K}$@4/@20、$\mathcal{C}$@4/@20、$k_{num}$、$k_{num}/k_{deck}$、$S_{num}/S_{deck}$），共 40 个比较点。首版 $\mu=10^{-3}$ Pa·s、$S_0=4.6\times10^{-10}$ Pa⁻¹、$H_f=50$ m、$t_f=3600$ s、$\Delta p=3$ MPa、$\ell=0.2$ m；由 $K_{Ic}=\sqrt{E'G_c}$ 与 $k=c_v\mu S_0$ 反推两类输入。全部落在**印刷精度**内（表内只有 2–3 位有效数字，最紧的是 I 类 $\mathcal{C}$@20：0.0013552 vs 0.0014，纯舍入差 3.2%）。

M2 进一步把 2.3 节的每个标度写三遍：参考式（硬编码原式）、字面式（用 JSON 内常数）、展开幂律式（只用 JSON 指数），要求三者相对误差 $<10^{-10}$。因此任何前因子或指数的静默漂移都可证伪，而不是仅仅「文档与代码看起来一致」。

### 14.3 文献推翻了哪些首版假设

| 量 | 首版 4.2 设定 | 文献区间 (rep) | 影响 |
|---|---|---|---|
| $\mu$ | $1\times10^{-3}$ Pa·s（清水） | $1\times10^{-3}$–$6\times10^{-2}$ ($5\times10^{-3}$) Pa·s，变黏滑溜水 | $c_v$、$C_L$、$\mathcal{K}$、$\mathcal{C}$ 全部重算 |
| $k_{field}$ | I–III 类反推 $1.0\times10^{-20}$–$1.0\times10^{-19}$ m² | $k_{eff}=1\times10^{-18}$–$1\times10^{-17}$ ($3.2\times10^{-18}$) m² | **首版低 2–3 个量级** |
| $k_{core}$ | 未区分 | $1\times10^{-17}$–$1.4\times10^{-16}$ ($1.2\times10^{-16}$) m²，岩心气测 | 与 $k_{eff}$ 差 rep 37.5×，**绝不可混用**；只有 $k_{eff}$ 进入 $T$ 与 $C_L$ |
| $S$ | $4.6\times10^{-10}$ Pa⁻¹ | $1\times10^{-10}$–$1\times10^{-8}$ ($1\times10^{-9}$) Pa⁻¹ | 必须由 $\phi(c_{pore}+c_{fluid})$ 推导：rep 校验 $0.10\times(1.0\times10^{-8}+4.4\times10^{-10})=1.044\times10^{-9}$，与代表值差 4.4% |
| $H_f$ / $t_f$ | 50 m / 3600 s | 20–100 (42) m / 3600–7200 (7200) s | $T\propto t/H^2$、$q=Q/H_f$ |
| Q 扫描 | 4–20 m³/min | 10–14 (12) m³/min（长 7 主体 10–12、优化 10–14） | 窗口窄 3.5× |
| $E'$ / $K_{Ic}$ | 21.7–42 GPa / 0.6–2.0 MPa√m | 8.3–64.6 (24.6) GPa / 0.2–1.2 (0.6) MPa√m | static 三轴模量，不可与 dynamic 超声模量混用 |
| $\Delta p$ | 3 MPa | **无文献约束**，仍固定为 3 MPa 单点 | M1 强制 `low==high==rep`；见 14.7 敏感性 |

**新恒等式**：因 $G_c=K_{Ic}^2/E'$，$\sigma_c=\sqrt{E'G_c/\ell}=K_{Ic}/\sqrt{\ell}$，**相场峰值应力与 $E'$ 无关**。脚本 `so_sigma_c_identity_check=1.0`。首版表 A 中 $\sigma_c$ 看似随 $E'$ 变化，实质是四类同时改了 $K_{Ic}$。

### 14.4 表 A′ / B′ / C′（可由 `--md` 重新生成）

表 A′：文献参数反算的相场材料参数（三行是「一致参数档」，即每个参数同取 low / rep / high，**不是**目标量的极值包络）

| 档位 | $E'$ (GPa) | $G_c=K_{Ic}^2/E'$ (N/m) | $\sigma_c=K_{Ic}/\sqrt{\ell}$ (MPa) | $c_v$ (m²/s) | $T_{field}$ | $C_L$ (m/√s) |
|---|---|---|---|---|---|---|
| low | 8.293 | 4.823 | 0.447 | 1.000e−05 | 9.000e−05 | 5.352e−07 |
| rep | 24.584 | 14.643 | 1.342 | 6.400e−07 | 2.612e−06 | 1.354e−06 |
| high | 64.467 | 22.337 | 2.683 | 1.667e−08 | 1.200e−08 | 2.185e−06 |

表 B′：文献参数下的无量纲分区（Q = 10 → 14 m³/min）

| 档位 | $\mathcal{K}$ @10 | $\mathcal{K}$ @14 | $\mathcal{C}$ @10 | $\mathcal{C}$ @14 | 结论 |
|---|---|---|---|---|---|
| low | 0.1248 | 0.1147 | 0.00432 | 0.00365 | $\mathcal{K}<3$、$\mathcal{C}\ll1$、$T\ll1$ |
| rep | 0.1334 | 0.1226 | 0.01628 | 0.01376 | $\mathcal{K}<3$、$\mathcal{C}\ll1$、$T\ll1$ |
| high | 0.0864 | 0.0794 | 0.03146 | 0.02659 | $\mathcal{K}<3$、$\mathcal{C}\ll1$、$T\ll1$ |

表 C′：由文献 $T$ 反解数值 deck 参数，并与 B0 冻结值 $k_{num}^{B0}=2.35\times10^{-18}$ m² 比较

| 档位 | $T_{field}$ | $k_{num}$ (m²) | $k_{num}/k_{num}^{B0}$ |
|---|---|---|---|
| low | 9.000e−05 | 1.350e−15 | 574.47 |
| rep | 2.612e−06 | 3.918e−17 | 16.67 |
| high | 1.200e−08 | 1.800e−19 | 0.08 |

**关键读法**：文献重算**推翻了首版的数值假设，但没有推翻分区结论**。三档一致参数下 $\mathcal{K}=0.0794$–$0.1334$（无一进入韧性主导 $\mathcal{K}>3$）、$\mathcal{C}=0.00365$–$0.03146$（滤失约束不激活）、$T\ll1$（存储主导）。变化的是数值，不是定性分区。rep 档 $k_{num}$ 需比 B0 高 16.67×，但 B0 冻结值仍落在表 C′ 三档包络内（偏低端）。$S_{num}=1.5\times10^{-10}$ Pa⁻¹ 是**相似常数**，不是储层 $S$，M7 明确禁止二者混淆。

### 14.5 为什么表 B′ 不是、也不能是「排量优选表」

1. **窗口错配**：B0/B1/B2 是**起裂前、存储主导、无裂缝扩展**窗口（$Damage_{max}=1.48\times10^{-4}\ll1$）。B1 A4 已由真实求解器证明 $G_c,k,\mu$ 对 $p_{bulk}(t)$ **不可辨识**；$\mathcal{K}$、$\mathcal{C}$ 是 Detournay (2016) / Savitski & Detournay (2002) **裂缝传播段**标度，在该窗口同样不可辨识。首版 4.2 硬结论 1「$\mathcal{K}$ 不含 $H,t$ ⇒ 必须 $q_{num}=q_{field}$」是传播段要求，不能在起裂前窗口当作判据。
2. **无可优化的目标**：起裂前 $p_w(t)\approx p_0+(Q_0/U)t$（Lecampion et al. 2017），$\beta=f/S$ 对 $Q$ 严格线性（B1 A2/A6 实测），没有阈值、没有拐点、没有极大值 ⇒ 「优选排量」在该物理阶段**没有定义**。
3. **不确定性过宽**：按 $T=k\,t/(\mu S H^2)$ 的真实单调角点组合，文献参数给出 $T\in[6.0\times10^{-10},\ 1.8\times10^{-3}]$，**跨 6.48 个数量级**（主要来自 $S$ 两个量级、$\mu$ 1.8 个量级、$H_f$ 1.4 个量级）。M7 固化了这个跨度，防止有人把区间偷偷收窄成单点以伪装出「清晰」的答案。
4. **q 口径缺口当前不可观测**：文献 rep $q=Q/H_f=0.2/42=4.762\times10^{-3}$ m²/s；B0 的 $Q_{total}=1.5\times10^{-4}$ m²/s 是 2D per-unit-thickness 流量（同量纲），比值 `so_q_ratio_field_over_b0=31.75`。但 B0 的注入率是为存储主导数值标定选的，不是 $\mathcal{K}$ 相似排量；由于 $\mathcal{K}$ 在起裂前不可辨识，这个 32× 缺口**在本轮实验里不可观测**，不能据此判定 B0 错误。让它可观测的前提是进入传播窗口。

**DFIT / mini-frac 标定**（闭合压力、$k_{eff}$、$S$、滤失系数）才是把表 A′–C′ 变成区块专属输入的前提；没有它，本节所有数字只能用于界定数值平台的相似区间与下一阶段的边界条件。这也正是参数表 `honesty` 字段的强制声明。

### 14.6 Δp 敏感性与新增风险 R8

$\Delta p$ 是全表**唯一无文献约束**的量，仍固定为 3 MPa 单点（M1 强制 `low==high==rep==3e6`）。若改取 10 MPa，high 角点 $C_L$ 由 $2.185\times10^{-6}$ 升到 $7.284\times10^{-6}$ m/√s，$\mathcal{C}$@10 由 0.0315 升到 `so_C_dp10_high_corner=0.1049`——**接近但仍未进入滤失主导 $\mathcal{C}\sim O(1)$**。M6 把这条探针固化为可证伪边界 $(0.1,1.0)$，参数表或公式漂移导致该探针越界时判据 FAIL。

> **R8（新增）**：把离线文献标度当成现场排量建议，或把 $\mathcal{K},\mathcal{C}$ 的传播段结论写进起裂前窗口。缓解：本节 14.5 的状态纪律、M8 强制路线稿含「起裂前 / 存储主导 / 不可辨识 / DFIT」四类措辞、脚本内 `RATE_WINDOW_DISCLAIMER` 常量、以及所有标度结论一律标注「离线复算，非求解器验证」。这是继 #124（`pmax` 局部过冲）、#127（传播段标度误当起裂前签名）之后**第三次同型纪律约束**。

### 14.7 决策点 2 的结论与后续路线

决策点 2 由「仍待确认」更新为：**已用公开文献参数完成重算，明确结论是不能产出排量优选表**；表 A′–C′ 交付为「数值平台相似区间 + 下一阶段边界条件」。要真正得到区块级排量结论，按代价从低到高有三条路：

1. **收窄输入**：取得区块 DFIT/mini-frac 与岩心数据，把 $S$、$k_{eff}$、$\mu(t)$、$H_f$、$\Delta p$ 从文献区间降为实测点，重跑 `--md` 即可（零 Fortran）。这能把 6.48 个数量级的 $T$ 包络压到可用精度，但**仍不产生「优选」**，只产生可信的相似映射。
2. **立项阶段 B3（物理前提改变，已由用户立项 → 第 15 节）**：以零 Fortran 的 `fracturetoughness` 单因子阶梯跨越起裂门槛，先确认平台**能否**给出刚度退化与实测 $G_c^*$ 夹逼，并分类损伤机制（#132 结论：**能**退化与夹逼，但**不能**局域化）。注意第 14.7 节写作时尚未立项，当时表述的「进入真实裂缝扩展窗口」在 B3 中已收窄为「起裂准入」；真正的传播段标度（Santillán 2018 排量阈值、Detournay 2016 M–K 顶点、$\mathcal{K}/\mathcal{C}$ 分区、排量优选）**仍需后续立项**。
3. **评估新增 Fortran**：等效 DFN / 显式 Carter 滤失 / 井筒限流 / `source_scaling: file` 所需的 CDF2 writer。会突破「零代码改动」优势，需单独评估（决策点 4 不变）。

### 14.8 状态纪律

本节所有数值均为**离线复算**（`check_shale_oil_params.py` 的纯 Python 计算），不是求解器验证结论；B0 冻结常数与 CI #125/#128 的真实求解器结果逐位不变。文献参数是**公开区间，不是区块实测**；$k_{core}$ 与 $k_{eff}$、static 与 dynamic 模量均不可混用。**不得用于现场配产或施工排量优选**。

### 14.9 CI #130 验收记录（`37bc756`，`success`）

run `42807197140`，`2026-11-10T11:24:25Z → 11:50:40Z`（26 min 15 s）。job `Build + Terzaghi` `conclusion=success`；step 14 `Build MEF90 + m_Flow + vDef` / step 15 `Run Terzaghi` / step 16 `Upload diagnostics` / step 17 `Upload Terzaghi result mesh` 全部 `success`。

**证据等级须诚实标注**：本轮 `actions/runs/{id}/logs`、`actions/jobs/{id}/logs` 与 `actions/artifacts/{id}/zip` 三个端点均返回 `401`（`Must have admin rights to Repository` / `Server failed to authenticate the request`），本机 git 凭据亦无 `actions:read` 权限，因此**无法像 #128 那样逐字节复核 verdict 文本与 `so_*` 指标**。本轮验收是**逻辑蕴含级**，不是日志逐字节级：

1. `gso_rc` 在 `Run Terzaghi` 末尾以 `if [ "$gso_rc" -ne 0 ]; then ... exit 1; fi` 参与延迟硬裁决（第 1450 行后的判据调用 `PIPESTATUS[0]` 捕获的是 `check_shale_oil_params.py` 的真实退出码，不是 `tee`）⇒ step 15 `success` **蕴含** `SHALE_OIL_PARAMS_RESULT: PASS` 且 M1–M8 全绿；
2. 其余 18 组 `GATE*_RESULT` 的延迟裁决（`g41_rc`、`g43_rc`、`gb1_rc`、`gb2_rc`、9 个 mpirun rc 等）与致命签名扫描都在**同一个 step 15** 内 ⇒ 该 step `success` 同时蕴含原有判据全绿、生产路径零扰动；
3. 判据 M8 静态检查 `build.yml` 确实含 `check_shale_oil_params.py` 与 `gate415_b2_verdict.log`，M7 检查 B0 deck 的 5 个冻结 token ⇒ 接线与冻结常数未被绕过；
4. 本地离线 `SHALE_OIL_PARAMS_RESULT: PASS`（M1–M8）+ `--selftest` 1 正例 + 33 变异全拒 + 3 不变性保持 + `--md` 与本节表格一致。

体系规模由 **19 组 verdict 同轮全绿**扩到 **20 组**（新增 `SHALE_OIL_PARAMS_RESULT`；前 18 条为 `GATE*_RESULT`，第 19 条为 `terzaghi_verdict` 汇总口径）。本轮 **零 Fortran / 零 deck / 零新增求解器运行 / 零生产 run 改动**，B0 冻结常数（$k_{num}=2.35\times10^{-18}$ m²、$S_{num}=1.50\times10^{-10}$ Pa⁻¹、$\beta=1\times10^{5}$ Pa/s、$Q_{total}=1.5\times10^{-4}$ m²/s）逐位不变。

> 若后续取得具备 `actions:read` 的凭据，应补做一次 verdict 文本的逐字节比对（重点：`SHALE_OIL_PARAMS_RESULT: PASS`、8 条 `=== M* ... PASS`、`so_T_span_decades=6.477121e+00`、`so_knum_rep_over_b0=1.667390e+01`），把本轮证据等级提升到 #128 的水平。这不改变验收结论，只补齐证据链。

---

## 15. 阶段 B3 立项：起裂准入门槛与 $G_c^*$ 夹逼（2026-11-11；CI #132 首轮否定结论见 15.9）

### 15.1 立项动机：物理前提改变，不是又一次参数扫描

B0/B1/B2 的真实读数把结论牢牢钉在**起裂前、存储主导**窗口：$\text{Damage}_{max}=1.480310\times10^{-4}\ll1$、$T=1.5667\times10^{-7}$、$\beta=f/S=10^{5}$ Pa/s。第 14 节由此证明「排量优选表」在该窗口**没有定义**。用户以「立项b3」把第 5.4 决策点 3 正式立项。

B3 要回答的问题比「进入传播段」更基础、也更可证伪：**这个数值平台究竟能不能给出一个局域化的起裂，门槛在哪里？** 关键是发现了一条**零 Fortran** 的路径——只改 `cs0001.fracturetoughness` 一个 deck 参数：

1. B1 A4 已由真实求解器证明 $G_c$ 对压力曲线**逐位不可辨识**（`gb1_beta_dev_gc = 0`），因此改韧性不会污染 B0 的解析锚 $\beta=f/S$；
2. $G_c$ 只进入相场的 $C_1/C_2$ 积分（`m_MEF90_DefMechAssembly.F90:1687-1688,2340-2341`），可以把「起裂」隔离成**单因子实验**；
3. 只有到 $\alpha\sim O(1)$ 时 $k(v)=k_{min}+(k_{max}-k_{min})v$ 才回写 Darcy，而 $T\sim10^{-7}$ 下即使渗透率升 20 倍仍不可见。

因此无需提高 `source`（会破坏 $p(t)$ 线性锚）、无需延长时间（会增加 CI 步数）、无需新增 Fortran。

> **#132 的答案是「半个能」**（15.9）：门槛确实跨过了、$G_c^*$ 也被实测夹逼，但损伤是**全域均匀软化**而非局域化裂纹。因此本阶段最终产出的是一个**可证伪的否定结论**，不是「起裂已建立」。

### 15.2 AT2 解析锚（源码级推导）

| 契约 | 源码位置 | 内容 |
|---|---|---|
| $c_w=1/2$ | `m_MEF90_DefMechAT2.F90:33` | `Constructor%cw = 0.5_Kr` |
| $a(\alpha)=(1-\alpha)^2,\ a'=-2(1-\alpha)$ | `m_MEF90_DefMechAT2.F90:53,72` | 刚度退化函数 |
| $w(\alpha)=\alpha^2,\ w'=2\alpha$ | `m_MEF90_DefMechAT2.F90:105,124` | 耗散势 |
| $C_1=G_c/(4c_w\ell)=G_c/(2\ell)$ | `m_MEF90_DefMechAssembly.F90:1687,2062` | 残差与 Jacobian 各一处 |
| $C_2=G_c\ell/(2c_w)\cdot A$ | `m_MEF90_DefMechAssembly.F90:1688,2063` | 梯度通道（$A=$ 各向异性矩阵） |
| $C_3=a'(\alpha)\Psi+C_1w'(\alpha)$ | `m_MEF90_DefMechAssembly.F90:1786` | 损伤残差主项 |

无梯度均匀平衡 $C_3=0$ 给出

$$\alpha_{eq}(G_c)=\frac{2\Psi}{G_c/\ell+2\Psi},\qquad G_c(\alpha)=\ell\cdot2\Psi\cdot\frac{1-\alpha}{\alpha}.$$

**起裂不用任意阈值**：$a(\alpha)$ 就是承载刚度，$a=1/2$ 对应

$$\alpha_{stiff}=1-\frac{1}{\sqrt2}=0.2928932188 .$$

### 15.3 用真实求解器数据反演，而不是猜 $\Psi$

把 **CI #125 真实 EXO**（`vals_nod_var2` 第 10 条记录最大值 $1.480310\times10^{-4}$）代入 $G_{c0}=27000$ N/m、$\ell=0.2$ m：

$$2\Psi=1.998714\times10^{1}\ \text{Pa},\quad \Psi=9.9936\ \text{Pa},\quad G_{c0}/\ell=1.35\times10^{5}\ \text{Pa},$$
$$G_c^*(\alpha=0.5)=3.997429\ \text{N/m},\qquad G_c^*(a=1/2)=9.650647\ \text{N/m}\quad(=G_{c0}/2797.7).$$

这不是未经检验的外推。**CI #127 的独立真实读数**给出了交叉验证：

| 量 | CI #127 实测 | 解析预测 | 相对误差 |
|---|---|---|---|
| `gb1_surface_ratio_gc`（$G_c$ 减半） | `2.000410` | `1.999408` | 0.05% |
| `gb1_beta_dev_gc` | `0` | `0` | 逐位 |

因为 $\text{surface}\propto C_1\langle\alpha^2\rangle\propto G_c\cdot G_c^{-2}=1/G_c$，$2.000410$ 本身就是 $\alpha\sim1/G_c$ 在真实求解器上的证据。另外 #125 真实 damage 场末步 $d_{max}/d_{mean}=6.996$、`d_{max}/d_{min}=26.3$，说明损伤**强局域化**（集中在局部应力带），为 N5 的判别阈值提供依据。

### 15.4 韧性阶梯（4 个单行差异 deck + 2 个既有 EXO 高档锚）

| rung | $G_c$ [N/m] | $G_c/G_{c0}$ | $\alpha_{eq}$（离线预测） | 来源 |
|---|---|---|---|---|
| `b0_inject` | 27000 | 1 | $1.480310\times10^{-4}$ | 既有 B0 EXO（CI #125） |
| `gc_half` | 13500 | $10^{-0.301}$ | $2.960182\times10^{-4}$ | 既有 B1 `gc` EXO（CI #127） |
| `gc2700` | 2700 | $10^{-1}$ | $1.478340\times10^{-3}$ | 新 deck，1 行差异 |
| `gc270` | 270 | $10^{-2}$ | $1.458929\times10^{-2}$ | 新 deck，1 行差异 |
| `gc27` | 27 | $10^{-3}$ | $1.289600\times10^{-1}$ | 新 deck，1 行差异 |
| `gc2p7` | 2.7 | $10^{-4}$ | $5.968602\times10^{-1}$ | 新 deck，1 行差异，预测跨 $\alpha_{stiff}$ |

所有 B3 deck **只改一行代码** `fracturetoughness`；`f, S, k, \mu, p_0`、11 步时间网格、BC、损伤律、`internalLength`、`residualstiffness` 与全部 Flow 键逐位不变。每个 run 前 `rm -f data/Beam1x10-tri3_out.gen`（CI #124 的 EXO append 教训），run 结束恢复生产 mesh。

> $\alpha\sim O(1)$ 时刚度退化会改变 $\Psi$ 的空间分布（力控制正反馈 / 位移控制负反馈并存），因此**高损伤档的精确 $\alpha$ 不可离线预测**——这正是必须让真实求解器回答的部分。N4 因此采用「至少一个 rung 跨门槛 + 实测夹逼区间 + 单调性 + 3 个数量级增长」的稳健命题，而不是逐位数值匹配。这是吸取 #124（`pmax` 局部过冲）与 #127（L4 传播段标度错配）两次「判据物理前提错误」后的设计。

### 15.5 `check_gate415_b3_initiation.py` 的 N1–N8

| 判据 | 可证伪命题 | 主要负例 |
|---|---|---|
| N1 | 6 组 EXO/log 各 11 条记录、EXO `pmax` 与 `Gate 3.2` 逐位一致、damage $\in[0,1]$、无 fatal | 截断 EXO、stale mesh、伪造崩溃签名 |
| N2 | $\beta=f/S$ 全程不变：低损伤 3 档 $<10^{-6}$、高损伤档 $<5\times10^{-3}$（容纳 $k(v)$ 反馈） | 破坏任一档压力斜坡 |
| N3 | **三种独立实现**（闭式 / 代数变形 / $C_3$ 数值求根）差 $<10^{-10}$；$C_1$ 必须在残差与 Jacobian 两处都在位；CI #127 表面能比值校准；低损伤 $\alpha\sim G_c^{-1}$ | 改 $c_w$、改 $C_1$ 前因子、删 $C_3$、冻结 damage |
| N4 | 至少一档 $\alpha\ge\alpha_{stiff}$；damage 随 $1/G_c$ 单调；增长 $\ge10^3$；解析 $G_c^*$ 落在实测夹逼区间的 1 个数量级邻域 | 全档未起裂、增长仅 2 个数量级、单调性反转 |
| N5（#132 后改为 N5′） | **机制分类**：B0 参考必须局域化（$3<d_{max}/d_{mean}<20$）；末档必须落入两个互斥带之一 —— `LOCALISED`（$>1.5$）或 `UNIFORM_SOFTENING`（$\le1.05$）；**中间模糊带 $(1.05,1.5]$ 一律 FAIL**。另有界、时间不可逆、显式报告 `-5`；UNIFORM 分支强制状态纪律声明在位 | 中间带、B0 判别器失效、damage 均值非正、damage 回退 |
| N6 | 低损伤段 $\text{surface}\sim G_c^{-1}$、注入功 $W$ 不变、能量行有限 | 冻结表面能、篡改末步 $W$ |
| N7 | 5 个 deck 单行差异且键为 `fracturetoughness`；AT2/C1/C3 源码契约；$\ge14$ 次 EXO 重建；`gb3_rc` 延迟硬裁决 | 加第二行差异、删 `rm -f`、删延迟 fail |
| N8 | B0/B1/B2 verdict 仍 PASS、生产 51 步日志无 B3 deck 名、状态纪律声明在位（#132 后增 `B3_REGIME` / `B3_INITIATION_ACCESS` 声明检查） | 伪造 verdict、B3 泄入生产 |

> **N5 → N5′ 不是放宽阈值**（15.9）：旧 N5 要求「每个 rung 都必须局域化」，等于向一个**诊断实验**索要**期望的物理结果**。CI #132 实测末档 $d_{max}/d_{mean}=1.000004$，均匀软化是真实物理解而非代码缺陷；沿用 #127→#128（L4）范式，把判据改成「机制分类 + 模糊带硬拒 + 状态纪律强制」的可证伪命题。合成正例也改为复现 #132 的真实诊断，使 `PASS` 恒等于「诊断与求解器一致」。

### 15.6 离线证据

- `--predict`：三实现一致，`2\Psi=1.998714\times10^{1}$ Pa，$G_c^*(a=1/2)=9.650647$ N/m；并打印 #132 真实 $\alpha$ 阶梯与解析 $\alpha_{eq}$ 的比值（末两档 $1.85\times$ / $1.68\times$，见 15.9）；
- `--selftest`：**1 正例 + 32 变异（含 5 不变性）全 PASS**。新增负例覆盖「中间模糊带」「B0 判别器失效」「damage 均值非正」；新增不变性覆盖「末档 $d_{mean}=d_{max}$ 仍是 UNIFORM」「末档 LOCALISED 也是合法分类」；
- **CI #125 真实 EXO/log 混合回放**（立项轮）：真实 `dmax_last=1.480309745e-4`、`dmax/dmean=6.996`、`\beta=1.000000e+5`、11 条 `Gate 3.2`/`Gate 4.5`/能量行、5 条 damage `-5`，N1–N8 全 PASS；
- **CI #132 真实 artifact 全量回放**（N5′ 更正轮）：从 `build-diagnostics` 恢复 6 组 B3 EXO + 6 组 log + 3 份 verdict + mesh + 生产 log（截断 zip 用 local-header 流式解压恢复 80 个成员，仅 `gc_half.gen` 由同包 B1 `gc.gen` 补齐），`--gate` 得 **N1–N8 全 PASS**、`GATE415_B3_RESULT: PASS`、`B3_REGIME: UNIFORM_SOFTENING`、`B3_INITIATION_ACCESS: NOT_ESTABLISHED`；**20 个 `gb3_*` 指标与 #132 原始 job logs 逐位一致、零漂移**，证明只改判据语义、未触碰任何物理观测量；
- 静态自检：`git diff --check`、`py_compile`、workflow 关键串计数、`Run Terzaghi` 整块 `bash -n` 全通过；
- **零 Fortran 改动**、零 deck 改动、零 B0/B1/B2 判据改动、零生产 run 配置改动；15 次 EXO 重建（原 11 次 + B3 4 次）。

### 15.7 状态纪律与能力边界（按 #132 实测收窄）

B3 **只测试起裂准入是否可达、并分类损伤机制**；`PASS` 证明的是**诊断与求解器一致**，不是「起裂已建立」。CI #132 的实际判定：

$$\texttt{B3\_REGIME}=\texttt{UNIFORM\_SOFTENING},\qquad \texttt{B3\_INITIATION\_ACCESS}=\texttt{NOT\_ESTABLISHED}.$$

不主张：

- Santillán (2018) 排量阈值、Detournay (2016) M–K 顶点、$\mathcal{K}/\mathcal{C}$ 分区或任何传播段标度；
- 裂缝长度、缝网面积、起裂压力绝对值、**排量优选表**；
- 已收敛的断裂力学解。damage SNES `-5` 残差平台（CI #84）按既有纪律**报告而不隐藏**（#132 共 26 条）；能量行是当前场的确定性积分，只作相对 / 可证伪观测量。

**两条防线，缺一不可**：

1. 不得把「$\alpha\to1$、刚度退化」解读为「裂纹已形成」——全域均匀软化没有裂纹面、没有缝长、没有扩展方向；
2. 也不得因为结论是否定的，就让诊断判据 FAIL —— 那会重演 #127（L4）「向诊断索要期望物理结果」的错误。

二者由 N5′ / N8 硬判：中间模糊带 $(1.05,1.5]$ 一律 FAIL；UNIFORM 分支强制 `B3_REGIME` / `B3_INITIATION_ACCESS` / `NO claim is made about` 声明在位；`build.yml` 另把三行分类显式回显到日志，避免 `PASS` 被误读。要真正进入传播段，仍需**另行立项 B3′**（15.10）。

### 15.8 CI 验收状态

- `a4424dc` → **CI #132 `failure`**：编译通过、既有 20 组 verdict 全绿（`GATE415_RESULT` / `GATE415_B1_RESULT` / `GATE415_B2_RESULT` 全 PASS），B3 的 **N1–N4 / N6–N8 全 PASS**；唯一失败 `GATE415_B3_RESULT: FAIL (N5)`，原因 `gc2p7 damage is spatially uniform (dmax/dmean=1 <= 1.50)`，由 step 末尾延迟 `gb3_rc=1` 触发 `exit 1`。
- 处置沿用 #124→#125 / #127→#128：**先验算物理阶段**（确认均匀软化是真实解而非装配 / 求解器缺陷），再改判据 + 补可证伪负例；**不改 deck、不改 Fortran、不改任何物理观测量**。更正轮离线证据见 15.6。
- `1d4307e` → **CI #133 `success`**（run `43023824796`，12:50:29 → 13:20:14 UTC，约 29 min 45 s）：编译 step 14 `success`；step 15 `Run Terzaghi` 完成且未触发延迟 `exit 1`，蕴含 **`GATE415_B3_RESULT: PASS`（N1–N8 全绿）** 与其余 20 组 verdict 全绿 ⇒ 体系由 **20 组扩到 21 组 verdict 同轮全绿**。`gb3_rc` 之后新增的分类 `grep` 回显接线由整块 `bash -n`（1515 行）验证。**证据等级为逻辑蕴含级，不是逐字节级**（15.11）。
- B3 由此收口，但结论仍是**否定的**：`PASS` 只证明机制诊断与求解器一致，不证明起裂已建立；`B3_REGIME = UNIFORM_SOFTENING`、`B3_INITIATION_ACCESS = NOT_ESTABLISHED` 不变。

### 15.9 CI #132 首轮实测：一个可证伪的否定结论

| 观测量 | #132 真实读数 | 解读 |
|---|---|---|
| `alpha_gc2p7` | `1.000000e+00` | AT2 越阈值后 snap 到完全损伤 |
| `local_gc2p7`（$d_{max}/d_{mean}$） | `1.000004` | **全域均匀软化，不是裂纹** |
| `local_b0_inject` | `6.996393` | 起裂前参考场本来是强局域化的 ⇒ 判别器本身有效 |
| `dmg_power_slope` / `surface_power_slope` | `-1.006199` / `-1.009178` | 低损伤 $\alpha\sim G_c^{-1}$、$\text{surface}\sim G_c^{-1}$ 均被真实求解器证实（解析 $-1$） |
| `gc_bracket_lo/hi` | `2.7 / 27` N/m | 夹逼区间包含解析 $G_c^*=9.650647$ N/m |
| `beta_dev_gc2p7` | `0` | 降 $G_c$ 4 个数量级仍对 $p(t)$ 逐位不可辨识（B1 A4 结论扩展到起裂窗口） |
| `surface_ratio_gc_real/pred` | `2.000410 / 1.999408` | #127 独立校准仍在位 |
| `damage_sneg5_msgs` | `26` | 残差平台显式报告，未隐藏 |

**排除「代码缺陷」的四条证据**：① N1 的 EXO ↔ `Gate 3.2` 逐位交叉验证通过、6 组各 11 条记录、356 节点 tri3；② N2 压力斜坡全程 $\beta=f/S$ 不变，说明 Flow 装配未被污染；③ N3 三种独立 $\alpha_{eq}$ 实现一致 $<10^{-10}$ 且 $C_1/C_3$ 源码契约在位；④ N6 表面能幂律与注入功不变、能量行有限。唯一不成立的命题是旧 N5 的「必须局域化」。

**解析锚在高损伤档的偏离是预期的，不是错误**：`--predict` 给出 $\alpha_{real}/\alpha_{eq}$ 依次为 `1.0000 / 1.0003 / 1.0031 / 1.0440 / 1.8538 / 1.6754`。前四档贴合均匀平衡解；`gc27`、`gc2p7` 明显偏离，正对应 15.4 已声明的「$\alpha\sim O(1)$ 时刚度退化改变 $\Psi$ 空间分布，精确值不可离线预测」。该偏离与「snap 到 $\alpha=1$ 且空间均匀」共同构成**力控制正反馈下的均匀软化**签名，而不是局域化裂纹。

**科学价值**：在 AT2 + 空间均匀驱动 + 单因子 $G_c$ 阶梯这组前提下，平台能给出刚度显著退化，但**不能**给出局域化断裂；因此「起裂准入」在此设置下 **NOT ESTABLISHED**。这是一个**关于能力边界的否定结论**，比凑一个 PASS 更有用：它明确排除了「只降韧性就能进入裂缝传播窗口」这条捷径。

### 15.10 若要真正进入传播段：B3′ 的最低准入条件

1. **空间非均匀驱动**（预裂纹 / 缺口 / 非均匀体力或非均匀 $\Psi$），使损伤有择优局域化位置；单纯均匀场在 AT2 下必然整体软化；
2. **损伤律评估**：AT1（存在弹性阈值、$w(\alpha)=\alpha$）比 AT2 更容易产生局域化带，但需新增或切换 Fortran 损伤模型，须单独立项与编译验证；
3. **传播段观测量**：缝长 / 缝高随时间、裂尖位置、$p_{net}$–$Q$ 标度、$\mathcal{K}$–$\mathcal{C}$ 分区；这些目前在 EXO 输出与能量诊断里都**没有定义**；
4. **`-5` 残差平台**：真正传播段的 SNES 收敛不能只靠「报告而不隐藏」，需要评估是否影响裂尖推进；
5. 上述任何一项落地前，B0/B1/B2 的 20 组 verdict 与 `beta=f/S` 解析锚必须保持零扰动。

### 15.11 CI #133 的证据等级（诚实边界）

CI #133 的 `conclusion=success` 与「verdict 文本逐字节复核」是**两件不同的事**，本轮只能给出前者蕴含后者，与 #130 的 14.9 同级：

| 端点 / 途径 | 返回 | 说明 |
|---|---|---|
| `actions/runs/{id}` | `200`，`conclusion=success` | run 级事实；13 个 step 标记 `success` |
| `actions/runs/{id}/jobs` | 先 `403 rate limit exceeded`（匿名 60 次/小时被后台轮询耗尽），重置后可用 | 取得 job id `130533369642` 与 step 状态 |
| `actions/artifacts/{aid}/zip` | `401 Requires authentication` | 限流重置后仍因认证不足拿不到 `build-diagnostics`（38,125,641 B） |
| `actions/jobs/{jid}/logs` | `403 Must have admin rights` | job logs 需 admin |
| 网页 `suites/.../artifacts/...`、`checks/.../logs/15`、run `/logs` | `404` | 匿名页面显示 "Sign in to view logs" |

**逻辑蕴含链**：`gb3_rc=${PIPESTATUS[0]}` 捕获判据退出码，step 15 末尾 `if [ "$gb3_rc" -ne 0 ]; then ... exit 1; fi`；其余 20 组 verdict、9+ 个 mpirun rc 与 `g*_rc` 全部在同一 step 延迟裁决。因此 **step 15 `success` ⇒ B3 N1–N8 全绿 + 原有 20 组全绿 + 生产零扰动**。补强证据：本地 `--selftest` 1 正例 + 32 变异全 PASS；#132 真实 artifact 全量 `--gate` 回放 N1–N8 全 PASS 且 20 个 `gb3_*` 指标逐位一致；N7/N8 静态契约校验 deck / workflow / 状态纪律接线在位。

**待补**：若日后取得 `actions:read` 凭据或登录态，应下载 `gate415_b3_verdict.log` 与 20 组既有 verdict，与 #132 做逐字节 IDENTICAL 比对，并确认日志中三行分类回显。在此之前，21 组全绿的表述一律标注「逻辑蕴含级」。
