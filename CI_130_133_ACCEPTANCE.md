# CI #130 / #133 验收结论汇编

> 整理日期：2026-11-26（仓库 HEAD `aaa15c1`）
> 用途：把两次**逻辑蕴含级**验收的事实、判据、科学结论与证据边界集中到一份可独立阅读的文档；原始记录在 `RATE_FRACABILITY_ROADMAP.md` 第 14 / 15 节与 `FLOW_FIX_CHECKLIST.md`。
> 本文不新增任何物理结论，也不修改 Fortran / deck / 判据 / workflow。

## 0. 一句话结论

| run | 主题 | 验收结论 | 科学结论 | 证据等级 |
|---|---|---|---|---|
| **#130** (`37bc756`) | 页岩油文献参数重算，表 A′–C′ | `success`，`SHALE_OIL_PARAMS_RESULT: PASS`，M1–M8 全绿，体系 19 → **20 组 verdict 同轮全绿** | **不能交付现场排量优选表**；表 A′–C′ 只能界定数值平台相似区间与下一阶段边界条件 | 逻辑蕴含级（非逐字节） |
| **#133** (`1d4307e`) | B3 起裂准入阶梯（AT2 韧性 $G_c$ 单因子） | `success`，`GATE415_B3_RESULT: PASS`，N1–N8 全绿，体系 20 → **21 组 verdict 同轮全绿** | `B3_REGIME = UNIFORM_SOFTENING`、`B3_INITIATION_ACCESS = NOT_ESTABLISHED`——**不是「已起裂」** | 逻辑蕴含级（非逐字节） |

两次验收都遵循同一条纪律：**判据向实验索要可证伪命题，而不是索要期望的物理结果**。因此 #132 首轮的 `FAIL (N5)` 与 #127 首轮的 `FAIL (L4)` 一样，是判据物理前提错误，修正后复验通过；物理结论本身没有被「调绿」。

---

## 1. 两次 run 的可核查事实

| 项 | #130 | #133 |
|---|---|---|
| run id / URL | `42807197140` · https://github.com/xsyuan211/mef90-flow-thermal/actions/runs/42807197140 | `43023824796` · https://github.com/xsyuan211/mef90-flow-thermal/actions/runs/43023824796 |
| head SHA | `37bc75638a3da45f357c283a712d8f27264cd581` | `1d4307ef4ab4d8779fc40bb61c3701e1ac7c30bd` |
| 创建 → 更新（UTC） | `2026-11-10T11:24:22Z → 11:50:41Z`（约 26 min 19 s；job step 口径 11:24:25 → 11:50:40，26 min 15 s） | `2026-11-11T12:50:29Z → 13:20:14Z`（约 29 min 45 s） |
| job / 关键 step | `Build + Terzaghi` `success`；step 14 `Build MEF90 + m_Flow + vDef`、15 `Run Terzaghi`、16/17 两组 upload 全 `success` | 编译 step 14 `success`；step 15 `Run Terzaghi` 完成且未触发延迟 `exit 1`；13 个 step 标记 `success` |
| 新增判据 | `check_shale_oil_params.py`（M1–M8，纯后处理） | `check_gate415_b3_initiation.py`（N1–N8，复用真实 EXO + 4 个新 deck） |
| 新增 verdict | `SHALE_OIL_PARAMS_RESULT`（第 20 组） | `GATE415_B3_RESULT`（第 21 组） |
| 延迟裁决变量 | `gso_rc`，`PIPESTATUS[0]` 捕获，step 15 末尾 `if [ "$gso_rc" -ne 0 ]; then ... exit 1; fi`（`build.yml:1567, 1823`） | `gb3_rc`，同型延迟裁决（`build.yml:1554, 1813`） |
| 前置失败轮 | 无 | #132 (`a4424dc`) `failure`：编译 + 20 组既有 verdict 全绿，B3 N1–N4 / N6–N8 全 PASS，唯一 `FAIL (N5)` |

---

## 2. #130：页岩油文献参数重算（表 A′–C′）

### 2.1 交付物与零扰动边界

| 文件 | 作用 |
|---|---|
| `vDef/data/shale_oil_params.json` | S01–S20 文献来源、25 个参数 low/rep/high 区间、单位换算、诚实边界、`window`（首版设定 / 文献映射 / 已发布表值）、`formulas` 契约、`numerical` B0 冻结常数 |
| `vDef/data/check_shale_oil_params.py` | M1–M8 判据；`--md` 重新生成本文表 A′–C′；`--selftest` 1 正例 + 33 变异 + 3 不变性；`--json` 指标导出 |
| `.github/workflows/build.yml` | 新增 `gso_rc` 延迟裁决与 `/tmp/shale_oil_params_verdict.log` artifact |

**零扰动**：不改 Fortran、不改任何 deck、不改生产 run、不新增求解器运行。判据只读 JSON / 路线稿 / B0 deck / workflow 文本。B0 冻结常数逐位不变：$k_{num}=2.35\times10^{-18}$ m²、$S_{num}=1.50\times10^{-10}$ Pa⁻¹、$\beta=f/S=1\times10^{5}$ Pa/s、$Q_{total}=1.5\times10^{-4}$ m²/s。

### 2.2 M1–M8 可证伪命题

| 判据 | 命题 |
|---|---|
| M1 | 参数表完整性 + 冻结常数（唯一无文献约束的 $\Delta p=3$ MPa 强制 `low==rep==high`） |
| M2 | 三实现公式恒等 + 指数契约：参考式（硬编码首版印刷公式）/ 字面式（JSON 常数）/ 幂律展开式（仅 JSON 指数），三者相对误差 $<10^{-10}$ |
| M3 | 逐位复现首版表 A–C 的 40 个已发布数字（印刷精度内，容差 4%） |
| M4 | $k_{core}$ 与 $k_{eff}$ 分离；$S$ 必须由 $\phi(c_{pore}+c_{fluid})$ 推导 |
| M5 | 标度单调性与指数符号一致 |
| M6 | 文献分区：存储主导、$\mathcal{K}<3$、$\mathcal{C}\ll1$ |
| M7 | $T$ 包络 → $(k_{num},S_{num})$ 相似映射 + B0 缺口 + $q$ 口径缺口声明 |
| M8 | 状态纪律 + 生产零扰动（路线稿含「起裂前 / 存储主导 / 不可辨识 / DFIT」四类措辞；workflow 接线在位） |

### 2.3 验收依据（为何 step `success` 等价于 M1–M8 全绿）

1. `gso_rc` 捕获 `check_shale_oil_params.py` 的**真实退出码**（不是 `tee`），并在 step 15 末尾参与延迟硬裁决 ⇒ step `success` **蕴含** `SHALE_OIL_PARAMS_RESULT: PASS`；
2. 其余 18 组 `GATE*_RESULT` 与 9 个 mpirun rc 的延迟裁决、致命签名扫描都在**同一个 step 15** ⇒ 同时蕴含原有判据全绿、生产路径零扰动；
3. M8 静态检查 workflow 确实含判据调用，M7 检查 B0 deck 的 5 个冻结 token ⇒ 接线与常数未被绕过；
4. 本地离线复核：`--gate` 得 `SHALE_OIL_PARAMS_RESULT: PASS`、`--selftest` 1 正例 + 33 变异全拒 + 3 不变性保持、`--md` 与路线稿表格一致。

### 2.4 科学结论（否定命题）

- 文献 rep 值：$k_{eff}=3.2\times10^{-18}$ m²（首版反推低 2–3 个量级）、$\mu=5\times10^{-3}$ Pa·s、$S=1\times10^{-9}$ Pa⁻¹、$H_f=42$ m、$t_f=7200$ s、$Q=12$ m³/min。
- 三档一致参数下 $\mathcal{K}=0.0794$–$0.1334$（无一进入韧性主导 $\mathcal{K}>3$）、$\mathcal{C}=0.00365$–$0.03146$（滤失约束不激活）、$T\ll1$（存储主导）——**文献推翻了首版数值假设，但没有推翻定性分区**。
- $T$ 的真实角点包络 $T\in[6.0\times10^{-10},1.8\times10^{-3}]$，跨 **6.48 个数量级**；M7 固化该跨度，禁止把区间偷窄成单点以伪装出「清晰」答案。
- **不能交付排量优选表的三条理由**：① 窗口错配（B0/B1/B2 为起裂前存储主导，$\mathcal{K}/\mathcal{C}$ 是传播段标度，B1 A4 已证 $G_c,k,\mu$ 对 $p_{bulk}(t)$ 不可辨识）；② 无可优化目标（$p_w(t)\approx p_0+(Q_0/U)t$，$\beta$ 对 $Q$ 严格线性，无阈值/拐点/极大值）；③ 不确定性过宽。
- 新增风险 **R8**：把离线文献标度当现场排量建议，或把传播段 $\mathcal{K}/\mathcal{C}$ 结论写进起裂前窗口。这是继 #124（`pmax` 局部过冲）、#127（标度阶段错配）之后第三次同型纪律约束。

---

## 3. #133：B3 起裂准入（AT2 韧性阶梯）

### 3.1 从 #132 的否定结论到 N5′ 更正

#132 实测末档 `gc2p7`：$\alpha=1.000000$、$d_{max}/d_{mean}=1.000004$，即**全域均匀软化而非局域化裂纹**。旧 N5 要求「每个 rung 都必须局域化」，等于向诊断实验索要期望的物理结果——判据前提错误，不是代码 bug。四条排除代码缺陷的证据：N1 的 EXO ↔ `Gate 3.2` 逐位一致（6 组 × 11 记录、356 节点 tri3）；N2 压力斜坡全程 $\beta=f/S$ 不变；N3 三种独立 $\alpha_{eq}$ 实现一致 $<10^{-10}$ 且 $C_1/C_3$ 源码契约在位；N6 表面能幂律与注入功不变。

**N5 → N5′ 不是放宽阈值**：B0 参考必须局域化（$3<d_{max}/d_{mean}<20$，实测 `6.996393`，证明判别器有效）；末档必须落入 `LOCALISED`（$>1.5$）或 `UNIFORM_SOFTENING`（$\le1.05$）两个互斥带之一，**中间模糊带 $(1.05,1.5]$ 一律 FAIL**；UNIFORM 分支强制状态纪律声明在位。$-5$ 残差平台显式报告条数（#132 为 26 条），不作硬判。沿用 #127→#128（L4）的更正范式：**只改判据语义，不改 deck / Fortran / 任何物理观测量**。

### 3.2 N1–N8 可证伪命题

| 判据 | 命题 |
|---|---|
| N1 | 6 组 EXO/log 各 11 条记录、EXO `pmax` 与 `Gate 3.2` 逐位一致、damage $\in[0,1]$、无 fatal |
| N2 | $\beta=f/S$ 全程不变：低损伤 3 档 $<10^{-6}$、高损伤档 $<5\times10^{-3}$（容纳 $k(v)$ 反馈） |
| N3 | 三种独立实现（闭式 / 代数变形 / $C_3$ 数值求根）差 $<10^{-10}$；$C_1$ 在残差与 Jacobian 两处都在位；#127 表面能比值校准；低损伤 $\alpha\sim G_c^{-1}$ |
| N4 | 至少一档 $\alpha\ge\alpha_{stiff}=1-1/\sqrt2$；damage 随 $1/G_c$ 单调；增长 $\ge10^3$；解析 $G_c^*$ 落在实测夹逼区间 1 个数量级邻域 |
| N5′ | 见 3.1：机制分类 + 模糊带硬拒 + 状态纪律强制 |
| N6 | 低损伤段 $\text{surface}\sim G_c^{-1}$、注入功 $W$ 不变、能量行有限 |
| N7 | 5 个 deck 单行差异且键为 `fracturetoughness`；AT2/C1/C3 源码契约；$\ge14$ 次 EXO 重建；`gb3_rc` 延迟硬裁决 |
| N8 | B0/B1/B2 verdict 仍 PASS、生产 51 步日志无 B3 deck 名、`B3_REGIME` / `B3_INITIATION_ACCESS` 状态声明在位 |

### 3.3 验收依据

1. `gb3_rc` 由 `PIPESTATUS[0]` 捕获并在 step 15 末尾延迟裁决 ⇒ step `success` 蕴含 `GATE415_B3_RESULT: PASS`（N1–N8 全绿）；
2. 其余 20 组 verdict 与全部 mpirun rc 同在该 step ⇒ 蕴含原有判据与生产路径零扰动；
3. 更正轮用 **#132 真实 artifact 全量回放**（截断 zip 以 local-header 流式解压恢复 80 个成员，仅 `gc_half.gen` 由同包 B1 `gc.gen` 补齐）：`--gate` 得 N1–N8 全 PASS、`B3_REGIME: UNIFORM_SOFTENING`、`B3_INITIATION_ACCESS: NOT_ESTABLISHED`，且 **20 个 `gb3_*` 指标与 #132 原始 job logs 逐位一致、零漂移**；
4. 本地 `--selftest` 1 正例 + 32 变异（含 5 不变性）全 PASS；`--predict` 三实现一致；静态 `git diff --check`、`py_compile`、workflow 关键串计数、`Run Terzaghi` 整块 `bash -n` 全通过。

### 3.4 实测否定结论（#132 真实读数，#133 复验未改变）

| 观测量 | 实测 | 解读 |
|---|---|---|
| `alpha_gc2p7` | `1.000000e+00` | AT2 越阈值后 snap 到完全损伤 |
| `local_gc2p7` ($d_{max}/d_{mean}$) | `1.000004` | **均匀软化，不是裂纹** |
| `local_b0_inject` | `6.996393` | 起裂前参考场强局域化 ⇒ 判别器有效 |
| `dmg_power_slope` / `surface_power_slope` | `-1.006199` / `-1.009178` | $\alpha\sim G_c^{-1}$、$\text{surface}\sim G_c^{-1}$ 被真实求解器证实 |
| `gc_bracket_lo/hi` | `2.7 / 27` N/m | 夹逼区间含解析 $G_c^*=9.650647$ N/m（$a(\alpha)=1/2$ 定义） |
| `beta_dev_gc2p7` | `0` | 降 $G_c$ 四个数量级仍对 $p(t)$ 逐位不可辨识 |
| `surface_ratio_gc_real/pred` | `2.000410 / 1.999408` | #127 独立校准仍在位 |
| `damage_sneg5_msgs` | `26` | 残差平台显式报告，未隐藏 |

解析锚 $2\Psi=1.998714\times10^{1}$ Pa 由 #125 真实 damage 场反演；高损伤档 $\alpha_{real}/\alpha_{eq}$ 偏离（末两档 `1.8538 / 1.6754`）是「$\alpha\sim O(1)$ 时刚度退化改变 $\Psi$ 空间分布」的预期结果，不是错误。

**能力边界**：B3 **能**给出刚度显著退化与 $G_c^*$ 夹逼，**不能**给出局域化断裂；不主张 Santillán 排量阈值、Detournay M–K 顶点、$\mathcal{K}/\mathcal{C}$ 分区、缝长/缝网面积/起裂压力绝对值/排量优选表。真正进入传播段需另行立项 B3′：空间非均匀驱动（预裂纹/缺口/非均匀体力）、AT1 损伤律评估、传播段观测量定义、`-5` 残差平台影响评估，且保持既有 20 组 verdict 与 $\beta=f/S$ 零扰动。

---

## 4. 证据等级：逻辑蕴含级为何成立，缺口在哪

两次验收的 `conclusion=success` 与「verdict 文本逐字节复核」是**两件不同的事**：

| 端点 / 途径 | 返回 | 影响 |
|---|---|---|
| `actions/runs/{id}` | `200`，`conclusion=success` | run 级事实可得 |
| `actions/runs/{id}/jobs` | 匿名限流后可用（`403 rate limit exceeded` → 重置后取得 step 状态） | 可确认 step 14/15/16/17 均 `success` |
| `actions/artifacts/{aid}/zip` | `401 Requires authentication` | 拿不到 `build-diagnostics`（#133 为 38,125,641 B），无法逐字节复核 verdict 文本 |
| `actions/jobs/{jid}/logs` | `403 Must have admin rights` | 无法读 job logs |
| 网页 `suites/.../artifacts/...`、`checks/.../logs/15`、run `/logs` | `404` / "Sign in to view logs" | 匿名网页路径不可用 |

**为什么仍然可以判定验收通过**：本 workflow 采用「延迟 `exit 1`」结构——所有判据的真实退出码都由 `PIPESTATUS[0]` 捕获，并在 step 15 末尾统一裁决。因此 step 15 `success` 在逻辑上**蕴含**该 step 内全部 verdict PASS；这是独立于日志可读性的结构事实。补充证据是本地离线回放与静态接线检查（见 2.3 / 3.3）。

**缺口与补证条件**：若日后取得 `actions:read` 凭据或登录态，应补做一次 verdict 文本逐字节比对，重点：
- #130：`SHALE_OIL_PARAMS_RESULT: PASS`、8 条 `=== M* ... PASS`、`so_T_span_decades=6.477121e+00`、`so_knum_rep_over_b0=1.667390e+01`；
- #133：`GATE415_B3_RESULT: PASS`、`B3_REGIME: UNIFORM_SOFTENING`、`B3_INITIATION_ACCESS: NOT_ESTABLISHED`，以及 20 组非 B3 verdict 与 #132 的逐字节一致性。

补证只提升证据等级，**不改变验收结论，也不改变任何物理结论**。

---

## 5. 复现与核查命令

```bash
# #130：真实文件判据 + 表格再生成 + 自检
python vDef/data/check_shale_oil_params.py
python vDef/data/check_shale_oil_params.py --md
python vDef/data/check_shale_oil_params.py --selftest

# #133：解析锚 + 自检（--gate 需 CI artifact 目录）
python vDef/data/check_gate415_b3_initiation.py --predict
python vDef/data/check_gate415_b3_initiation.py --selftest
python vDef/data/check_gate415_b3_initiation.py --gate <artifact_dir>

# 接线静态核查
grep -n 'gso_rc\|gb3_rc' .github/workflows/build.yml
```

## 6. 提交与文件索引

| 提交 | 内容 |
|---|---|
| `37bc756` | #130：页岩油文献参数重算表 A′–C′ + M1–M8 判据接入（20 组 verdict） |
| `abd5821` | #131：验收记录文档收尾 |
| `a4424dc` | #132：B3 起裂准入阶梯 N1–N8 首轮（唯一 N5 失败，物理前提更正依据） |
| `1d4307e` | #133：N5′ 机制分类更正 + 状态纪律接线（21 组 verdict） |

原始记录：`RATE_FRACABILITY_ROADMAP.md` 第 14 节（14.9 为 #130 验收）、第 15 节（15.9 为 #132 实测、15.11 为 #133 证据边界）；`FLOW_FIX_CHECKLIST.md` 研究项条目与「CI #132 / #133」专节。

## 7. 状态纪律声明

本文所有数字分为两类：**CI 真实读数**（2.4 / 3.4 节标注来源者）与**离线文献复算**（表 A′–C′、$G_c^*$ 解析锚）。后者不是求解器验证结论；文献参数是公开区间，不是区块实测；$k_{core}$ 与 $k_{eff}$、static 与 dynamic 模量不可混用。**不得用于现场配产或施工排量优选**。B3 的 `PASS` 只证明机制诊断与求解器一致，**不证明起裂已建立**。
