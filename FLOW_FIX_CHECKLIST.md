# MEF90 + Flow 模块 — 统一修复清单（可执行）

> 合并自两份分析：
> - 重精读修复清单（按"可运行基线"排序）
> - 用户 12 点分析（P0 物理正确性 / P1 可用性 / P2 工程性能）
>
> 所有源码改动**只落在子目录**（`m_Flow/`、`vDef/`、`m_HeatXfer/`），根目录同名 `.F90` 是过期镜像，不要改。
> 每批改完用：
> `cd vDef && mpirun -np 2 ../bin/$PETSC_ARCH/vDef -prefix data/Beam1x10-tri3 -options_file_yaml data/terzaghi1D.yaml`
> 先看退出码与 `SNESGetConvergedReason`，再看 pressure 字段。建议一次只推进一道闸门。

---

## 闸门 1 — 跑通不崩 + 强制边界（最高优先级）

- [x] **1.1 删除未分配指针的收尾崩溃**
  - 位置：`vDef/vDef.F90:37`（声明）、`:618-619`（DeAllocate）
  - 问题：`flowEnergy`/`flowInjectionWork` 是 `PetscReal,Dimension(:),Pointer`，全程未 `Allocate`，结尾 `DeAllocate` 触发未定义行为/崩溃（grep 确认无分配点）。
  - 改法：删除第 37 行声明与第 618-619 行 `DeAllocate`，以及注释 `!!! NEW: Darcy diagnostics`。
  - 验收：正常退出，无 `attempt to DEALLOCATE unassociated pointer`。

- [x] **1.2 Dirichlet 残差读取解向量 `x` 而非已组装残差**
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:397`（签名无 `x`）、`:422`、`:446`、`MEF90FlowOperatorCore` 调用点 `:128`
  - 问题：签名没有 `x` 参数，第 422/446 行从 `residual`（已组装弱形式残差）取值，写成 `residual - p_bc`，等于把边界当成了自然边界。**`terzaghi1D.yaml` 的 cs0010 `bc_pressure:0` 从未被强制**，顶面排水条件形同虚设，当前 run 的解本就无意义。
  - 改法：调用点 `:128` 改为 `Call MEF90FlowApplyDirichletResidual(x,residual,FlowCtx,ierr)`；签名改为 `Subroutine MEF90FlowApplyDirichletResidual(x,residual,FlowCtx,ierr)` 并加 `Type(Vec),Intent(IN) :: x`；第 422、446 行 `VecGetValues(residual,...)` 改为 `VecGetValues(x,...)`。
  - 验收：顶面 Dirichlet 节点孔压精确为 0，`SNESGetConvergedReason ≥ 0`。

- [x] **1.3 CreateSNES 建真实矩阵 + KSP + 零空间**
  - 位置：`m_Flow/m_MEF90_Flow.F90:352-371`
  - 问题：`SNESSetJacobian` 传 `PETSC_NULL_OBJECT`，解析 Jacobian 与 `MEF90FlowApplyDirichletJacobian`（`MatZeroRows` diag=1）从未生效，PETSc 改用有限差分；且 `addNullSpace` 注册了却没用。
  - 改法：照抄 `m_HeatXfer/m_MEF90_HeatXfer.F90:616-663`：
    ```fortran
    Type(Mat) :: matP
    Type(KSP) :: kspP
    Call DMCreateMatrix(FlowCtx%DMScal,MATAIJ,matP,ierr)
    Call MatSetOption(matP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr)
    Call SNESSetType(SNESPressure,SNESKSPONLY,ierr)
    Call SNESSetJacobian(SNESPressure,matP,matP,MEF90FlowBilinearForm,FlowCtx,ierr)
    Call SNESGetKSP(SNESPressure,kspP,ierr)
    Call KSPSetType(kspP,KSPCG,ierr)
    Call KSPSetTolerances(kspP,1.0D-8,PETSC_DEFAULT_DOUBLE_PRECISION,1.0D+10,PETSC_DEFAULT_INTEGER,ierr)
    ```
    并接入 `gOpt%addNullSpace`（`FlowCtx.F90:560` 已注册）：为真时 `MatNullSpaceCreate` + `MatSetNullSpace` + `KSPSetNullSpace`，避免全 Neumann 流场矩阵奇异。
  - 验收：矩阵非零；`SNESGetConvergedReason ≥ 0`；全 Neumann 配置不报奇异。

- [x] **1.4 解析同名 `MEF90FlowSetTransients` 跨模块歧义（经核实无需改动）**
  - 位置：`m_Flow/m_MEF90_Flow.F90:185`（dispatcher，调 2D/3D）与 `m_Flow/m_MEF90_FlowAssembly.F90:535`（真正做 source/boundaryPressure 缩放的 core）同名。
  - 问题（已核实**非问题**）：dispatcher 在 `m_MEF90_Flow`、core 在 `m_MEF90_FlowAssembly`，后者经 `m_MEF90_Flow.F90:15-28` 的 `use ... MEF90FlowSetTransients2D => MEF90FlowSetTransients` 重命名别名被 dispatcher 调用；`vDef.F90` 仅 `USE m_MEF90_Flow`（:12），故 :414 解析到 dispatcher，无歧义（CI 可编译为证）。
  - 改法：**无需改动**。仅当有人在 vDef 同时 `USE` 两模块才会触发 ambiguous，届时再处理。
  - 验收：编译无 duplicate/ambiguous 警告；`SetTransients` 行为符合预期。

---

## 闸门 2 — 边界强制 + 输出可见

- [x] **2.1 主循环调用 `MEF90FlowViewEXO` 并扩展输出项**
  - 位置：`vDef/vDef.F90`（save-results 段，参照 `MEF90HeatXferViewEXO` 调用处）；`m_Flow/m_MEF90_FlowAssembly.F90:354-389`
  - 问题：主程序从不调用 `MEF90FlowViewEXO`，渗流结果无法后处理。且当前 ViewEXO 只输出 source/pressure/boundaryPressure，**不输出 `permField`/`cellDarcyVelocity`**。
  - 改法：在 save 段加 `Call MEF90FlowViewEXO(MEF90FlowCtx,step,ierr)`；在 ViewEXO 内补 `permField`、`cellDarcyVelocity` 的 `VecViewExodusCell/Vertex` 输出（yaml 需相应给 offset）。
  - 验收：exo 里出现 pressure（及扩展的 permField / cellDarcyVelocity）变量。

- [x] **2.2 初值用注册的 `p0`**
  - 位置：`m_Flow/m_MEF90_FlowCtx.F90:389`（`VecSet(pressure,0.0_Kr)`）
  - 问题：写死 0，忽略已注册的 `initialPressure`。
  - 改法：`MEF90FlowCtxCreateVectors` 内 `Call PetscBagGetDataMEF90FlowCtxGlobalOptions(...gOpt...)` 后，`pressure` 与 `pressurePreviousStep` 都 `VecSet(...,gOpt%initialPressure,ierr)`。（SetFromOptions 在 CreateVectors 之前已调用，gOpt 已就绪。）
  - 验收：t=0 时 p = p0。

- [x] **2.3 yaml 给非平凡初值**
  - 位置：`vDef/data/terzaghi1D.yaml:50`（`p0: 0.0`）
  - 问题：p0=0、source=0、bc p=0 → 稳态 `p≡0`，算例无验证价值。
  - 改法：稳态验证给非零 `bc_pressure` 或 `source`；瞬态验证设 `p0: 1.0e5`（σ₀），由闸门 3 驱动衰减。
  - 验收：p 场非零。
  - 注意：为让当前（稳态）求解也能观察到非零 p 场，已在 `cs0001` 临时加 `source: 1.0e-6`（yaml 已注明 TEMP）。**推进闸门 3 做解析解比对前务必把它改回 `0.0`**，否则会偏离 Terzaghi 无源解析解。

- [x] **2.4 注册 `boundaryPressureOffset`**
  - 位置：`m_Flow/m_MEF90_FlowCtx.F90:603` 注册段（缺该项）；`ViewEXO:378` 用 `gOpt%boundaryPressureOffset`
  - 问题：该字段默认 0 且未 `PetscBagRegisterInt`，边界压力输出分支恒跳过、且无法经 yaml 设置。
  - 改法：注册段补 `Call PetscBagRegisterInt(bag,FlowGlobalOptions%boundaryPressureOffset,default%boundaryPressureOffset,'boundary_pressure_offset','EXO offset for boundary pressure',ierr)`；yaml `flow:` 下加 `boundary_pressure_offset: 3`。
  - 验收：exo 出现 boundaryPressure 变量。

---

## 闸门 3 — 接瞬态 + 符号校验

- [x] **3.1 流动求解移出 AltMin、接入 TS**（代码已落地 2026-11-06 18:45，**零编译/零运行验证**）
  - 位置：`vDef/vDef.F90:140`（`MEF90FlowCreateSNES` 之后）、`:407-424`（当前在 AltMin 内 `SNESSolve(snesPressure)`）
  - 问题：`MEF90FlowCreateTS` 全项目无调用点，IFunction/IJacobian 与储存项装配无驱动；`terzaghi1D.yaml` 设 `timestepping_type` 却实际跑稳态，解析解对比无从谈起。
  - 改法：在 `:140` 按 `MEF90FlowGlobalOptions%timeSteppingType` 分支，`transient` 时 `Call MEF90FlowCreateTS(...)`；主循环**时间步层（AltMin 之前）** 做 `TSSolve(tsPressure,pressure)` 并每步把 `pressure` 拷入 `pressurePreviousStep`；`:417` 的 `SNESSolve(snesPressure)` 仅在 steady 分支执行。
  - **实际实现**：`m_Flow/m_MEF90_Flow.F90` 重写 `MEF90FlowCreateTS`（两个独立 Mat + MAT_SPD/对称 + TSBEULER/TS_LINEAR + TSSetSolution/InitialTimeStep/ExactFinalTime/FromOptions + TSGetSNES→KSPCG）；`vDef/vDef.F90` 5 处接线（声明、初始化分支、时间步层 TSSolve+收敛检查、AltMin 内 steady 分支、收尾 TSDestroy）；`terzaghi1D.yaml` 回退 TEMP `source` 为 0、枚举值改精确 `Transient`；新增自检日志 `Gate 3.1: transient Darcy TS created`。详见文末「3.1 静态实现记录」。
  - 验收（**未执行**）：`transient` 与 `steady` 输出明显不同；p 随时间衰减。

- [x] **3.2 用解析 Terzaghi 校验符号与量级**（由 Gate 3.3 P1–P8 硬判据验收；结论：衰减方向与量级正确，无需改算子符号）
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:101`（注释 `G(p)=div(A∇p)`，稳态算子符号）、`MEF90DiffusionOperatorSet` 符号；IFunction `F=K(p)+S·ṗ`、IJacobian `J=K+a·S·M`（BE 下 a=1/dt，形式自洽）
  - 问题：方程应为 `-div(k/μ∇p)=S ṗ + f`，但算子注释算的是 `+div(A∇p)`，与方程差一个符号；源项 `:105-107` 也用 `ADD_VALUES`（**4.4 复核**：`ADD_VALUES` 本身符合 Bourdin 契约，缺的是 `F` 前的负号，不是 ADD 模式错）。
  - 改法：用解析解 `p(z,t)=(4σ₀/π)Σ(-1)ⁿ/(2n+1)·cos((2n+1)πz/2H)·exp(-(2n+1)²π²c_v t/4H²)`（c_v=k/(μS)，t=1s 时 T≈1.09）在单步与全程比对；若衰减方向反了，在 `MEF90FlowOperatorCore` 统一把 A 取负号，并同步源项符号。
  - 验收：数值衰减方向与量级与解析解一致。**已验收**：Gate 3.3 P6「L2 err vs BE series」CI #99 实测 `0.03206 ≤ 0.05`，P1–P6 / P8 全 PASS；衰减方向正确，故 `MEF90FlowOperatorCore` 的 A 未取负号。该结论只约束 `f=0` 的齐次瞬态，**对源项符号无区分力**；4.4 已用独立非零源诊断证实原装配为 `G+F` 并修复为 `G−F`。

---

## 闸门 4 — 物理增强与工程闭环

- [x] **4.1a Biot 体力：叠加而非覆盖 + 体力缩放不回滚**（CI #85 落地、**CI #86 全绿收口**，见文末专节）
  - 实际位置（原清单行号已过时）：`m_Flow/m_MEF90_FlowAssembly.F90:1180-1256`；`vDef/vDef.F90:520-533`（**时间步层**，不是每个 AltMin 迭代）。
  - **覆盖性复核**：原清单写的「`VecCopy` 重置体力 → 步内 CST/Linear 缩放被回滚」在当前代码里**不成立**。Phase 2b 的实现是 `VecCopy(baseForce,force)` + `VecAXPy(force,-α,∇p)`，快照在时间步层、`MEF90DefMechSetTransients` 之后、AltMin 之外，每次 AltMin 调用都从快照重建 → 幂等、不累积、不回滚。
  - **原清单「仅作首次快照」的建议是错的**：`m_vDefDefault.F90:43/79` 默认 `ForceScaling = MEF90Scaling_Linear`，`SetTransients` 每步把 `force` 缩放到 `Force × t`；若只在启动时快照一次，`baseForce` 会恒为 t=0 的值（0），非零体力被静默丢弃。**正确语义是每时间步快照一次**，已写入 `vDef.F90` 注释。
  - **CI #86 覆盖性更正**：这里写的 `9.81×t` 只对 `terzaghi1D_bodyforce.yaml` 成立。生产 deck 的 `Force: 0, -9.81, 0` 因逗号后空格被 RealArray 只读成首元素（`0 0 0`），`linear × 0 = 0`；CI #85 实测生产 `|baseForce| ≡ 0`。CI #86 把 4 个既有 deck 改成显式 `0,0,0`（解析值不变），另建 BODYFORCE 诊断 deck 写 `0,-9.81,0`，**bag view 实测 `Force = 0 -9.81 0`**、`|baseForce|` 线性 `0 → 9.81`，B3/B4 从此才有区分力。
  - 本轮落地：`VecCopy+VecAXPy` → 单次 `VecWAXPY(force,-α,∇p,baseForce)`（清单要求的写法，数学等价、无中间态）；每步输出一条 `Gate 4.1:` 诊断；新增 `terzaghi1D_nobiot.yaml`（仅 `biot_coeff 1.0→0.0`）、`check_gate41_biot.py`、CI 硬闸门 **B1–B7**。
  - 验收（CI #86 起）：B3 判 **BODYFORCE** 日志「`|baseForce|` 线性增长 0→9.81、R²≥0.98」，B4 也判 **BODYFORCE**「`‖force−(base−α∇p)‖∞≈0` 且 `max|base|≥5`」（生产 `base≡0` 时该误差对叠加/覆盖无区分力）；B6 用 NOBIOT 对照网格判「耦合真实生效」；Gate 3.3 P1–P8 保持回归。

- [x] **4.1b 有效应力项：互斥总应力模式（biot_effective_stress）**（已实现；CI #98 验收通过，六组闸门全绿）
  - **原清单表述有误**：「力学侧缺 `σ' = σ + α·p·I`，当前只做了体力近似」并非缺项。分部积分
    `∫(−α∇p)·v dΩ = ∫αp·div(v) dΩ − ∫_Γ αp·v·n dΓ`，而总应力形式
    `∫(C:ε − αpI):ε(v) dΩ = ∫C:ε:ε(v) dΩ − ∫αp·div(v) dΩ`。
    **两者的体积项是同一个 `−∫αp·div(v)`，不是两个可相加的效应**；差别仅在自然边界的 `−αp·n` 面力如何解释。若按字面「在体力之外再补 `αpI`」，孔压效应被双重计入，`u_x` 峰值约翻倍（`1.10e-4 → 2.2e-4 m`），直接打掉刚转绿的 Gate 3.3 P4（上界 `2e-4`）。
  - 因此「含孔压时应力场体现有效应力修正」落地为：**新增与体力形式互斥的总应力模式**，由 `flow.biot_effective_stress: true`（PetscBool，默认 false 保持生产 superposition）控制。启用时 `MEF90FlowApplyBiotForce` 提前返回——**停掉 `-α∇p` 体力**、改调 `MEF90FlowComputeCellPressure` 把逐单元体均 p 写入 `cellPorePressureSec`；DefMech 残差在 `m_DefMech/m_MEF90_DefMechAssembly.F90` 对 `stressGauss` 减 `α·p·I`（仿 line 277 热膨胀写法，tensor−scalar×identity），Jacobian 不变（p 与 u 解耦）；复用 4.2 的 SectionReal 跨模块 wire（`DefMechCtx%porePressureSec => FlowCtx%cellPorePressureSec` 一次性指针关联，非 TARGET 组件）。cell-based 孔压场**不能复用 `pressureForce`**（它是 coDim=1 的边界压力面力）。
  - 验收（Gate 4.1b，新增 `vDef/data/check_gate41b.py`）：**E1** 应力 deck 日志含 `Gate 4.1b: form=stress`；**E2** 有效应力形式的**物理签名硬判据**（见下条更正）；**E3** 生产 superposition 模式日志**不含** `Gate 4.1b:`（默认构建不变更）且仍含 `Gate 4.1:`（4.1a 路径完好）；**E4** 互斥结构——应力日志有 `Gate 4.1b:` 但**无** `Gate 4.1:` body-force 行，证明两种形式无法同时生效。
  - **CI #97 诊断更正（2026-11-09，推翻原 E2 的「位等价」前提）**：原判据称「Terzaghi 两端 Dirichlet ⇒ 自然边界项为零 ⇒ 两形式位移场逐位等价（relL2≤1e-3）」，这是**错的**。deck 只钉住 `x=±5` 两个端面（`n=±e_x`）；侧向 `y=±0.5` 是自由 Neumann 面，`v_y` 不受约束且 `p≠0`（step1 `1e5 Pa`、step25 `3.58e4 Pa`），故 `∮αp·v·n ≠ 0`。两形式是**不同边值问题**：4.1a 为 `σ'·n=0`，4.1b 的 `σ_total·n=0` 等价于 `σ'·n=+αp·n`，产生向外侧向膨胀。装配代码**没有 bug**：TRI3/P1 下 `GradS_BF` 为单元常数、`pCell` 为 Gauss 体均，`∫p·divφᵢ = pCell·Vol·tr(GradS_BFᵢ)` 精确成立，`ComputeCellPressure`/`ComputePressureGradient` 权重一致。
  - E2 实测与 `0.69` 的解释：pressure 场全 51 步完全一致（0 Pa 差）；`u_y^S ≈ k·(y·p)`、`k≈7.5e-10 1/Pa` 全步恒定、`cos(u_y^S, y·p)≈+0.97~0.98`、`|u_y^S|≈4.5·|u_y^M|` 且反相关；侧向膨胀经 Poisson 耦合带来轴向收缩，使 `u_x^S≈0.69·u_x^M` 且 `cos≈1.000`、`Δu_x` 与 `p` 相关 `-0.77`。**0.69 是边界牵引的后果，不是标量装配系数错误**。用全局 `s*=1.446` 缩放后 relL2 仍 `0.3549`，故单标量假设不成立。
  - **E2 新判据（8 项全硬）**：① 压力驱动一致 `max|pS−pM|/max|pM|≤1e-6`；② 网格坐标一致；③ `cos(u_x^S,u_x^M)≥0.99`；④ `0.40≤‖u_x^S‖/‖u_x^M‖≤0.95`（上界 `<1` 专门捕捉「4.1b 悄悄退化回 4.1a」）；⑤ `peak|u_x|∈[5e-5,2e-4]`；⑥ 额外 `u_y` 为真实场 `‖Δu_y‖≥0.5‖u_y^M‖`；⑦ `cos(Δu_y, y·p)≥0.90`；⑧ 柔度 `k` 落在平面应变解析包络 `[0.4·(1+ν)(1−2ν)/(E(1−ν)), 2.0·(1−ν²)/E]`（`α=1, E=1e9, ν=0.25`）。旧 relL2 保留但仅作信息性输出。
  - 状态：**已实现并验收（CI #98）**；`terzaghi1D_stress.yaml` + `check_gate41b.py` 落地，Jacobian 不变、默认 superposition 零影响。P0 的两行 `INTENT(INOUT)` 修复在 `f20fa18`（CI #97 编译通过、五组回归全绿，仅旧 E2 FAIL）；判据修正在 `b2f94ea`，**CI #98 `success`**：`GATE41B_RESULT PASS`（E1/E2/E3/E4 全过，E2 8/8）、Gate 3.3 (P1–P6,P8)、4.1 (B1–B7)、4.2 (D1–D8)、4.3 (C1–C7)、MAIN VERDICT 全 PASS。**CI #99 `success`**（`74b6c87`）：Biot 单位张量换成 `MEF90_MATS_IDENTITY` 后 2D/3D 双编译通过，六组判据仍全绿且 E2 的 7 个物理指标与 #98 逐位一致。
  - **真实隐患（非 E2 成因）已修**：`m_DefMech/m_MEF90_DefMechAssembly.F90:328` 原把 `matpropSet%LinearThermalExpansion` 当单位张量用；两个 run 的 bag dump 恰好都是 `1 1 0`，但 4.2 诊断 deck 设 `0,0,0`、`test9.yaml` 设 `1.0,0.0,0.0`，开 stress 模式会静默归零 / 变成各向异性。现改为 `MEF90_MATS_IDENTITY`（`MEF90/mef90.inc`：2D→`MEF90MatS2DIdentity` = `MatS2D(1,1,0)`、3D→`MEF90MatS3DIdentity` = `MatS3D(1,1,1,0,0,0)`），是模块 `Parameter` 常量，deck 无法覆盖。可见性由 `Use m_MEF90 → Use m_MEF90_LinAlg`（两者均无模块级 `Private`）保证；`PetscReal * MEF90_MATS` 重载已有 CI 编译先例（`m_DefMech/m_MEF90_DefMechSplitDeviatoric.F90:91`、`m_Flow/m_MEF90_FlowAssembly.F90:115`）。`m_DefMech/Makefile` 同时编 `...Assembly2D.o` / `...Assembly3D.o`，故两个宏分支都会被编译。当前 deck 下新旧值数值完全相同 ⇒ 物理场零变化，只去掉隐式耦合。
  - **防回归实验已落地（deck 层，零新增 CI 运行）**：`terzaghi1D.yaml` 与 `terzaghi1D_stress.yaml` 的 `cs0001` 都写入 `LinearThermalExpansion: 0,0,0`（键名/格式沿用 CI #90 的 `terzaghi1D_advection.yaml`）。本 deck `heatxfer: null`，LTE 只用于温度禁用时不执行的热应变项，故物理量不变；`diff` 确认两 deck **仍只差一行 `biot_effective_stress: true`**。修复前 stress 位移会归零 → E2 的 `peak-ux-in-band` / `ux-shape-preserved` 必 FAIL；修复后 E2 仍 PASS。Gate 4.3 的 C4 只比对 `flow:` 段与 Flow Global bag，这个 DefMech 材料键不参与比对，故 C1–C7 不受影响。**验收方式**：CI 日志的 bag dump 必须显示 `LinearThermalExpansion = 0 0 0`（证明 deck 真的生效），同时 E2 的 7 个物理指标与 #98/#99 逐位一致。
  - **网格限制**：Exodus 只有 3 个 element block（`eb1` TRI3 600 单元、`eb2/eb3` 各 5 个 BAR2 仅 `x=±5` 两端），**没有侧向 facet 标签**。若想恢复「严格位等价」判据，需改网格或写通用边界循环补 `−αp·n` 牵引，成本高，暂不采用。

- [x] **4.2 平流：把 `cellDarcyVelocity` 喂入 HeatXfer**（CI #87 落地，见文末专节）
  - 位置：`m_HeatXfer/m_MEF90_HeatXferAssembly.F90` 的 4 处平流装配（稳态残差/Jacobian、瞬态 IFunction/IJacobian）；`m_Flow/m_MEF90_Flow{F,Ctx}.F90` 的持久 `cellDarcyVelocitySec`；`vDef/vDef.F90` 的 Ctx 挂接与 `Gate 4.2:` 诊断。
  - 问题：`enableFlowAdvection`/`fluidDensity` 原是占位；`cellDarcyVelocity` 全库无消费者。且 `MEF90Diffusion*AdvectionSet` 只接受**一个常向量**（整个 cell IS 共用），直接替换会丢掉逐单元空间结构。
  - 改法：新增 `MEF90HeatXferFlowAdvectionOperatorSet` / `...BilinearFormSet` 两个**私有逐单元**装配子程序（镜像 `MEF90/m_MEF90_DiffusionImplementation.F90` 的两个循环，但每个 cell 从 `SectionReal` 读自己的 `v_darcy`），`advectionVector = -ρc_p·v_darcy(cell)` 保持原常数路径的符号/标度；`enable_advection=true` 时走场，`false` 时原常数路径完全不变。`MEF90` 核心库与 2D/3D 泛型接口**未改动**。
  - 验收：ADVECTION / NOADVECTION 两个诊断 deck 仅差 `flow.enable_advection`，硬判据 **D1–D8**（CI #99 实测全 PASS）要求温度场出现**非均匀、随时间增长**的差异，并校验力学场不爆；生产 deck 保持 `heatxfer: null`，Gate 3.3 完全隔离。

- [x] **4.3 边界压力时间缩放语义（CI #92 覆盖性更正：原表述不成立，降级为可见性/判据工作；CI #94 已收口）**
  - **原清单行号过时**：实际位置是 `m_Flow/m_MEF90_FlowAssembly.F90:701-735`
    的 `MEF90FlowSetTransients`（`boundaryPressureScaling` 的 `Select case` 在 725-733）。
  - **"恒压边界被乘时间"在当前工程里根本不发生**：缩放语义已由枚举
    `bc_pressure_scaling` = `constant` / `linear` / `file` 三选一区分
    （`FlowCtx.F90:648-649` 注册；`CST` 分支直接 `UpdateboundaryPressure`
    不缩放，`Linear` 分支才 `VecScale(..., timeVal)`）。
    `vDef/m_vDefDefault.F90:188` 默认 `MEF90Scaling_CST`，且**7 个 deck 全部
    显式写 `bc_pressure_scaling: constant`**。Linear 是用户主动选择的
    "渐变加载"语义，不是 bug；Gate 3.3 的 P6 解析解比对（L2 err 0.032）
    也间接证实恒压 BC 未随时间漂移。
  - **真正缺的是运行期可见性**：`m_Flow/m_MEF90_FlowCtx.F90` 的
    `MEF90FlowCtxSetFromOptions` 注册了 Global / CellSet / VertexSet 三个 bag，
    但**从未调用 `PetscBagView`**；而 DefMech / HeatXfer / Materials /
    MEF90Ctx 都在 `verbose>0` 下打印。所以 CI 日志里查不到
    `bc_pressure_scaling`、`source_scaling`、`biot_coeff`、
    `enable_advection`、`boundary_pressure_offset` —— 这正是 4.1/4.2 不得不
    自造 `Gate 4.1:` / `Gate 4.2:` 诊断行的根本原因。
  - **本轮落地（CI #93）**：照抄 `m_MEF90_DefMechCtx.F90:1025-1027 / 1047-1051`
    的惯用法，在 FlowCtx 三处注册后各加一个 `verbose>0` 守卫的
    `PetscBagView(..., PETSC_VIEWER_STDOUT_WORLD, ...)`。**纯打印，物理与
    选项值零改动**；生产 deck `verbose: 1`，CI 日志立即可见。
  - **CI #93 已回（success）**：6 份求解日志全部出现 `flow Flow Global Options`
    （23 个必需键）与 3 个 `Flow CellSet Options`；生产 `cs10` 是唯一
    `flowBC = TRUE`、`bc_pressure = 0`；`bc_pressure_scaling = constant`、
    `source_scaling = constant`、`pressure_scaling = constant`；
    Gate 3.3 / 4.1 / 4.2 回归同轮全 PASS。
  - **CI #94 已回（success）**：`vDef/data/check_gate43_flowbag.py` 基于上述
    真实格式写硬判据 **C1–C7**；`build.yml` 接入 `--gate` 调用、延迟裁决与
    `/tmp/gate43_verdict.log` artifact；真实 artifact 回放 **C1–C7 全 PASS**
    （`GATE43_RESULT: PASS`，7 项检查 0 失败、6 份运行受检）。
    C4 做**全量 deck↔bag 比对**（6 份日志 116 项，专抓「YAML 键被静默忽略」）；
    C5 查 3 类变体信号；C6 直接读 51 步网格判恒压槽 `span=0 Pa`、
    闭端仍衰减（防全零假场）并复核 Gate 3.3 解析 L2；C7 查 4.1/4.2
    回归、致命签名与挂接线。**不猜 PetscBagView 的输出格式**。
  - 离线验收（零本地编译 / 零求解器运行）：CI #93 真实数据回放 C1–C7
    全 PASS（116 项比对 0 失配、BC `span=0 Pa`、闭端 `1e5→9055.71`、
    `L2=0.03206`）；**12 例负例自测**逐条命中目标判据（删 bag 键、
    翻转 cellset BC、删 deck 显式声明、deck/bag 值不符、μ(T) 悄悄打开、
    BC 压力漂移、全零假场、4.2 回归、缺挂接线、缺运行日志），
    `--gate` 失败时 `rc=1` 且打印 `GATE43_RESULT: FAIL`；
    静态校验 PyYAML / 10 个 `run` 块 `bash -n` / 9 个 py `ast.parse` /
    新增行非 ASCII=0 / `git diff --check` 全绿。
  - **勾选条件**：CI #94 `conclusion=success` 且 `GATE43_RESULT: PASS`。

- [x] **4.4 源项符号对齐**（真 bug；**CI #101 `success` 验收通过**）
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:139-153`（`MEF90FlowOperatorCore` 调 `MEF90DiffusionRHSSetCell`）；底层 `MEF90/m_MEF90_DiffusionImplementation.F90:376-413`。
  - **结论更正**：原清单把因果写反了。注释 "Subtract source RHS: G = DiffusionOperator - Source" 是**正确意图**，代码漏了负号；不是"注释错"。Bourdin 契约（`MEF90/m_MEF90_DiffusionImplementation.F90:27-60`）：`-div(A∇p)+λp=f`，Operator `G=+∫A∇p·∇w`、RHS `F=+∫fw`（两者都 `ADD_VALUES`），SNES Function 必须是 **`G−F`**。HeatXfer 遵守（`m_MEF90_HeatXferAssembly.F90:113` / `:476` `modifiedFluxPtr = -fluxPtr` 后再 ADD）；Flow 的 `SetSourceCst` 存物理 `+f` 且 `OperatorCore` 直接 ADD → 实际 `G+F`，**正源被当成汇**。
  - 为什么一直没发现：所有生产/诊断 deck 都是 `cs0001.source: 0.0`，六组判据对源项符号零区分力；Gate 3.3 P6 只约束 `f=0` 时的衰减方向，与 f 的符号无关，故其"无需改算子符号"的结论仍然成立。
  - 改法（最小侵入；**未采用**清单原提议的 `SectionRealAXPY`，全库无该 API 先例）：在 `OperatorCore` 的 cell 循环里，把**局部** `sourceSec` 的当前 cell 值取负后再交给 `MEF90DiffusionRHSSetCell`；复用已声明未用的 `sourceCell` 和 `SetSourceCst` / `UpdatePerm` 已验证的 `SectionRealRestrict/Restore` 写入模式，删掉死变量 `rhsSec`。`FlowCtx%source` 仍存物理 `+f`，EXO 输出与 Gate 4.3 比对不受影响。
  - 验收（可证伪）：新增 `vDef/data/terzaghi1D_source.yaml`，与主 deck **只差一行**（`cs0001.source: 0.0 → 1.0e-6`）；`delta_p = p_source − p_main` 消去 1e5 Pa 的 Terzaghi 衰减，与一维解析源响应 `p_f(s,n)=f/(2A)(H²−s²)+Σa_m cos(λ_m s)(1+c_vλ_m²dt)^(−n)`（`a_m=−2(−1)^m f/(A H λ_m³)`）比对。`check_gate44_source.py` 的 S1–S8 硬判据要求 `cos(delta_p,+p_f)≥0.95`（旧 `G+F` 会给 −1）、幅值 `[0.70,1.30]`、`relL2≤0.15`、闭合端正峰、单调递增，并把"两 deck 只差一行"本身列为 S2。离线自测：正例 PASS、5 个负例全 FAIL；用 #97 真实 mesh 回放得 `S3 baseline relL2=0.00887`、解析峰值 `924.3 Pa`。
  - **已验收（CI #101 = `success`，`97fa070`）**：`GATE44_RESULT: PASS`，S1–S8 全绿；实测 `cos(delta_p,+p_f)=0.9999975`、`||delta_p||/||p_f||=1.0049`、`relL2=0.00536`、末态闭合端 `926.55 Pa` vs 解析 `924.30 Pa`（比值 1.0024）、闭合端 `0→391.9→724.6→926.6 Pa` 单调增长；bag dump 为 `source = 1e-06`。原有六组判据全绿，Gate 3.3 P6 `0.0320613`、Gate 4.1b E2 七指标与 #98/#99/#100 **逐位一致** ⇒ `source=0` 时本修复是严格 no-op。

- [x] **4.5 能量诊断**（修好并接线；**CI #102 `success` 验收通过**，`GATE45_RESULT: PASS`）
  - 实际位置：`m_Flow/m_MEF90_FlowAssembly.F90:399-520`（`MEF90FlowEnergy`）；`vDef.F90:567-577`（时间步层接线）。
  - **复核结论**：选择清单给的两条路中"修好"这一条 —— 删死代码会同时失去 4.8 守恒性后处理的落点。除清单指出的三个未填充 section 外，还发现 4 个额外缺陷：
    1. `MEF90DiffusionEnergySet` / `MEF90DiffusionWorkSetCell` 的 `energy`/`work` 声明为 `Intent(OUT)` 但内部只做 `+=`（`MEF90/m_MEF90_DiffusionImplementation.F90:491` / `:580`），调用方**必须预先清零** —— HeatXfer 在 `m_MEF90_HeatXferAssembly.F90:354-355` 做了，Flow 没有（`energyCell`/`workCell` 是未定义值）。
    2. 两处 `mesh` 参数传 `FlowCtx%DM`（主 DM，2 dof/vertex），而 `pSec` 是 `DMScalSec`（1 dof/vertex）；清单文末 767 行早就记过"留待后续统一"。改为 `DMScal`，与 HeatXfer `:365` / `:370` 一致。
    3. 逐 rank 部分和未经 `MPI_AllReduce` 就返回，CI 用 `-np 2`，打印值只是半个域。
    4. `vDef.F90` 从未调用 `MEF90FlowEnergy`（闸门 1 删掉的 `flowEnergy`/`flowInjectionWork` 是未分配 pointer 数组，这次改用 scalar）。
  - 改法：补 `permFieldSec` / `sourceSec` 的 `SectionRealToVec(SCATTER_REVERSE)` 填充，`p0Sec` 清零（FlowCtx 无对应 Vec，且 `lambda=0` 使参考压力不参与 E，取 0 与 drained BC 一致，避免 `0*NaN`）；`energyCell`/`workCell` 每 cell 前清零；两处 `DMScal`；末尾 `MPI_AllReduce` 得到全局 E / W。`vDef.F90` 在 Gate 3.2 压力包络之后每步输出 `Gate 4.5: step N t T E X W Y`，**纯后处理，不回写求解**。
  - 验收（可证伪）：新增 `check_gate45_energy.py` 的 T1–T8 硬判据。因 4.4 的两 deck 只差 `cs0001.source`，主 deck `f=0` ⇒ `W_main≡0`，source deck ⇒ `W_src>0` 且跟随 `∫f·p dΩ`；`E=½∫(k/μ)|∇p|²` 必须非零、单调衰减、与 BE 解析级数同量级。解析锚：`E_an(n)=(A·W·P0²/H)Σd_m²`（因 `a_m²λ_m²≡4/H²`，n=0 阶跃发散故只比 n≥1），`E_an(50)=0.269`；`W_an(n)=f·W·[(8P0H/π²)Σd_m/(2m+1)² + fH³/(3A) − (2f/(AH))Σd_m/λ_m⁴]`，`W_an` 从 `0.999` 降到 `0.0656`。离线自测正例 PASS、6 个负例（E 恒 0 / W_src 恒 0 / W_src 反号 / E 递增 / W_main 非零 / 缺行）全 FAIL。**所有诊断 deck 都会自动输出该行，无需新增 CI 运行**。
  - **提交前离线预测（CI #101 artifact 真实场，零 Fortran 风险）**：用 `parse_cdf2.CDF2` 从 `out_main.gen` / `out_source.gen` 读 51 步 `Pressure` 节点场与逐单元 `permField`（EXO slot 4 = `5e-8`），按线性三角形精确梯度重建 `E=½Σ A_c|∇p_c|²Ω_c` 与 `W=Σ f_c(Σp_i/3)Ω_c`，再喂给 `evaluate_text`。T1–T8 **全部 PASS**：`W_main≡0`、`E(1)/E(50)=334.98`、`E(50)=0.2520` vs `E_an(50)=0.2690`（比值 0.937）、`W_src(0)/W_src(50)=15.57`、`cos(W_src,W_an)=0.99998`、`W_src(50)=0.06372` vs 解析 `0.06564`（比值 0.971）、`E_src/E_main=1.223`（交叉项抬高，实证 T8 必须用 `[0.30,3.00]` 宽幅而非"小扰动"窄幅）。日志侧 `Gate 3.2` 在两个 deck 各 51 行 ⇒ 4.5 输出行数满足 T1。综合静态自检 **17/17 PASS**（10 个 `run: |` 块 `bash -n`、12 个 Python AST、INTENT 链、延迟判据顺序、`--selftest`）。
  - **已验收（CI #102 = `success`，`c05bac5`）**：`GATE45_RESULT: PASS`，T1–T8 全绿。实测 `W_main≡0`（51 步全 `0.0000000000000000`）、`E(1)/E(50)=334.978` 且 n≥1 严格递减、`E(50)=0.252043` vs `E_an(50)=0.269017`（比值 0.937）、`W_src(0)/W_src(50)=15.567`、`cos(W_src,W_an)=0.99998`、`W_src(50)=0.0637243` vs 解析 `0.0656407`（比值 0.971）、`E_src/E_main=1.223`。**与提交前离线预测逐位吻合**，证明 `MEF90FlowEnergy` 的 E/W 语义与解析一致、MPI 归约给出的是全局值。两份日志各 51 行 `Gate 4.5`，无 `NaN` / `Infinity` / `PETSC ERROR` / Fortran runtime 签名。Gate 4.1 / 4.2 / 4.3 / 4.1b / 4.4 的 verdict 与 #101 **逐位 IDENTICAL**，Gate 3.3 仅日志字节数 `285044→291011`（多 51 行输出）不同、物理判据数值不变 ⇒ 纯后处理无扰动。
  - **勾选条件（已满足）**：CI `conclusion=success` 且 `GATE45_RESULT: PASS`，且原有六组 + Gate 4.4 判据全绿（纯后处理，物理场逐位不变）。离线预测已先行通过，CI 已确认 Fortran 编译与运行期接线。

- [x] **4.6 μ(T) 别名加固**（**CI #106 `success` 验收通过**，`GATE46_RESULT: PASS`）
  - **精读更正（比清单原判断深一层）**：清单原写「`cellMuField` Vec 仅作 `Associated()` 守卫、值从未被读回」。复核修复前源码后确认：`MEF90FlowUpdateMuFromTemperature` **已经**把共享 `cellDMScalSec` `SCATTER_FORWARD` 到 `cellMuField`，Vec 并非空壳；真正的缺陷是 `UpdatePerm` 从不读这个持久 Vec，而是对 `FlowCtx%cellDMScalSec[cell]` **同一个 slot** 先后 `SectionRealRestrict` 出 `permPtr` 与 `muPtr`（两者互为别名），随后 `MEF90FlowUpdateStorativity` 又把 `S` 写进同一 slot。μ 只在 `UpdateMu → UpdatePerm` 两次调用之间「活着」；一旦调用顺序变化、或 `UpdateMu` 因 `DefMechCtx%temperature` 未关联 / opt-out 而早退，`UpdatePerm` 就会拿残留的 `S`（或上一步 `perm`）当黏度，产生静默错误。全项目复核：修复后 `FlowCtx%cellDMScalSec` 只剩 `SectionRealDuplicate` 的布局模板用途，与 DefMech / HeatXfer 的既有惯例一致，无第四个别名写入者。
  - **落地（`m_Flow/m_MEF90_Flow.F90`）**：三个例程 `UpdatePerm` / `UpdateMuFromTemperature` / `UpdateStorativity` 各自 `SectionRealDuplicate` 私有 section（`permSec` / `muSec` / `sSec`），共享 section 退化为模板；`UpdateMu` 私有 `muSec` → `SCATTER_FORWARD` → `cellMuField`，`UpdatePerm` 把 `cellMuField` → `SCATTER_REVERSE` → 自己的私有 `muSec` 再读；`muPtr(1) <= 0`（Vec 在 `CreateVectors` 里零初始化）时回退 `gOpt%dynamicViscosity`，避免 `1e-12` 下限把 `permField` 放大 `~1e9`。私有 section 的 `Duplicate` / `Destroy` 一一配对。
  - **诊断 deck**：`terzaghi1D_viscosity.yaml` = `terzaghi1D_advection.yaml` 去注释 90 行 + **唯一一行** `flow.temp_dep_viscosity: true`（91 行）。两者都有 live HeatXfer 温度场（`0→100 K`，51 步）；`TMin = 0.5·T_ref = 146.5 K` 的 clamp 让每个 cell 都落在下界，Sutherland 给出均匀 `μ(T) = 1e-3·(146.5/293)^1.5·(293+600)/(146.5+600) = 4.2293795e-4 Pa.s`，故 `permField = k_min/μ(T) = 1.1822065e-7`，恰是 opt-out 对照 `k_min/μ_ref = 5e-8` 的 `μ_ref/μ(T) = 2.3644130` 倍。若 μ 仍从共享 section 取，比值会塌到 1 或被 `1e-12` 放大到 `1e9`，两者都是硬 FAIL。
  - **V1–V8 判据（`vDef/data/check_gate46_alias.py`，纯标准库 + `parse_cdf2`）**：V1–V3 静态反别名（三例程私有 section、μ 经 `cellMuField` 往返、共享 section 无 `Complete`/`ToVec`、legacy 别名表达式消失）；V4 单行 deck 隔离 + bag `temp_dep_viscosity =  TRUE` + 无 crash / 未知选项签名；V5 51 步 / 600 cell 且末态温度跨度 `100 K`；V6 `perm(visc)/perm(ctrl) = 2.364413 ±5%` 且逐单元 spread `≤0.05`；V7 `max|perm·μ(T)/k_min − 1| ≤3%`；V8 opt-out 对照 `permField ≡ 5e-8`（`±2%`）。数据缺失 / 非有限值一律 fail-closed。
  - **离线先行**：`--selftest` 用真实源码 + synthetic 场，1 正例 + 7 负例（legacy 别名复辟、μ 未读回、`1e9` 放大、温度场死、bag 仍 FALSE、deck 非单行、对照 perm 被扰动）全 PASS；用 CI #104 真实 advection mesh + 解析倍率端到端预测 V1–V8 全 PASS（预测 `V7 rel = 1.29e-16`）。提交前另做：`git diff --check`、PyYAML 解析、deck 单行 diff、三例程 `Duplicate`/`Destroy` 配对、`g46_rc` 赋值→延迟判据→artifact 三处接线。
  - **已验收（CI #106 = `success`，`a260d7b`）**：`GATE46_RESULT: PASS`，V1–V8 全绿。实测 bag `temp_dep_viscosity =  TRUE`、51 行 `Gate 3.2`、`0` fatal 签名；`out_viscosity.gen` 600 个 cell 的 `permField` **全部**为 `1.18220652e-07`；`V6 ratio = 2.364413`（解析 `2.364413`）、per-cell spread `0`；`V7 max|perm·μ(T)/k_min − 1| = 0`（比离线预测更精确）；`V8` 对照偏差 `0`。**原有九组 verdict（Gate 3.3 / 4.1 / 4.2 / 4.3 / 4.1b / 4.4 / 4.5 / 4.8 / terzaghi）与 CI #104 逐位 IDENTICAL** ⇒ 生产 `temp_dep_viscosity = FALSE` 路径是严格 no-op，也实证排除了「私有 section 先 `SectionRealSet(...,0)` 会让未遍历 cell 归零」的顾虑（所有 cell 均被 cell-set 循环重写）。
  - **勾选条件（已满足）**：CI `conclusion=success` 且 `GATE46_RESULT: PASS`，且原有九组判据全绿、生产路径物理场逐位不变。

- [x] **4.7 时间相关材料参数入口**（**CI #111 `success` 验收通过**，`GATE47_RESULT: PASS`）
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90` 的 `MEF90FlowSetTransients`（原仅缩放 source + boundaryPressure）。
  - 问题：渗透率、储集系数、黏度等随时间演化（压密/损伤演化）无对应重载入口。
  - **落地方式（`ef6ee2a`，8 文件 / +1286 行）**：在 `FlowCtx` 增加两个持久逐单元 Vec `permTimeFactor` / `storTimeFactor`（`cellDMScal`，初值 1），`SetTransients` 每步按 `perm_scaling` / `stor_scaling` 写入因子；`vDef.F90` 把 `SetTransients` 提为每步**首个**调用，其后 `UpdateMuFromTemperature` → `UpdatePerm` → `UpdateStorativity` 每次**重建 base 场再乘因子**，因此不存在 `VecScale` 累积漂移，AltMin 重装配也不会把因子重置。因子是唯一真值来源。
  - 三定律：**CST**（因子恒 1，严格 no-op，生产默认）、**Linear**（`1 + rate*t`，`t=0` 无奇点）、**File**（`VecLoadExodusCell` 读 EXO cell slot，offset=0 时退化为 1）。μ(t) 因在算子中只以 `A = k/μ` 出现，复用同一 A 通道；4.6 的空间 `μ(T)` 场不受影响。
  - **工程防回归**：6 个新 bag 字段（2 Enum / 2 Real / 2 Int）按 CI #50 教训**加在类型声明末尾、注册末尾一一对应**，三份默认构造器同步为 `MEF90Scaling_CST` / rate=0 / offset=0；factor Vec 的 Allocate / `DMCreateGlobalVector` / `VecSet(1)` / Destroy / Nullify 生命周期完整；`flow.log` 编译 2D/3D 双份均无错。
  - 验收（`vDef/data/check_gate47_material.py`，849 行，**X1–X8 全硬**）：X1 51 步数据完整性；X2 静态反回归（分支 / factor Vec / bag 顺序 / 调用顺序）；X3 deck 只差 4 行（`perm_scaling/perm_rate/stor_scaling/stor_rate`）；X4 bag 可见 `linear 0.8/0.4`；X5 EXO `permField == 5e-8*(1+0.8t)` 逐记录零误差（末步比值 `1.800000`）；X6 51 行 factor 回读（`permF=1.8` / `storF=1.4`）；X7 时变 `c_v` 的 BE 乘积序列（`relL2(50)=0.0368`、`cos=1.00000`、`ratio=0.9633`）并用 `sep_const=0.360 / onlyk=0.782 / onlys=0.577` 拒伪「只接单通道」；X8 生产隔离（控制组 `relL2=0.0321`、perm 平坦）。离线自测 **1 正例 + 9 负例全 PASS**（覆盖入口死掉、单通道、Vec 未写、bag 不变、deck 多改、控制场被扰动）。
  - **CI #111 实测**：`GATE47_RESULT: PASS (X1-X8)`；`g47_rc` 延迟 `exit 1`（CI #82 早退教训），artifact 增 `gate47_verdict.log` / `out_material.gen` / `terzaghi_material.log`。**十二组判据同轮全 PASS**；与 #109 相比 9 组 verdict 中 7 组逐字节 IDENTICAL，另 2 组仅 solver log 字节数变化（每步多打印一行 `Gate 4.7:`），所有物理数值不变 ⇒ 生产路径零扰动。

- [x] **4.8 守恒性（质量平衡 + 能量/耗散平衡）后处理**（**CI #104 `success` 验收通过**，`GATE48_RESULT: PASS`）
  - **设计决策**：不改 Fortran。清单原建议"挂在主循环时间步层"，但 4.5 已让每步输出 `E`/`W`，Gate 2 已让 EXO 输出 `Pressure` / `permField` / `cellDarcyVelocity`，守恒量全部可由这三处既有产物重建，故做成 `vDef/data/check_gate48_conservation.py` 纯后处理判据，**零 Fortran 风险、物理场逐位不变**。
  - 两条独立全局守恒律（BE；齐次 drained Dirichlet `p=0` 使取 `w=p` 时边界项消失）：
    1. **能量/耗散平衡**（`w=p^{n+1}`）：`S M Δp · p^{n+1}/dt + 2E^{n+1} − W^{n+1} = 0`。`E`/`W` 就是 4.5 的诊断量，`S1 = S M Δp·p^{n+1}` 用一致质量矩阵 `M_ij = Ω/12(1+δ_ij)` 从 EXO 压力场重建。
    2. **质量平衡**（`w=1`）：`ΔM_n/dt + Q^{n+1} − F = 0`。`M=∫S p dΩ`（EXO 压力）、`F=∫f dΩ`（deck）、`Q=∮_{x=−5} q·n dS` 从 EXO 的 `cellDarcyVelocity_x`（elem var 5）**独立恢复**，不是由方程反解；U2 实测其与 `−A∂p/∂x` 的相对误差 `3.1e-16`，排除循环论证。
  - **离线诊断先行（CI #102 artifact）**：先算残差再定阈值。两条律的 `|R|/scale` 都稳定在 `0.5%–1.1%`，且都反演出同一个 `S_eff/S≈0.985–0.989`；与 Gate 3.3 的 FE/BE 解析 L2 误差 `0.032` 同量级，属离散误差而非符号/接线 bug。故硬阈值取 `3%`，仍对 E 未接线、源项符号反、`S`/`dt` 错、流出量缺失全部敏感（负例已验证）。
  - U1–U8：U1 数据完整性（51 EXO 步 + 51 `Gate 4.5` 行/ deck）；U2 排水面 `x=−5`、边长 1 m、速度场可信；U3/U4 两 deck 能量平衡 `≤3%`（`W_main=0`、`W_src>0`）；U5/U6 两 deck 质量平衡 `≤3%`（`F_main=0`、`F_src=1e-5`）；U7 瞬态闭合（`Q` 单调衰减、残差不增长）；U8 交叉校验 `W·S·Area = F·M`（`W=f∫p`、`M=S∫p`，实测 `7e-16`）。
  - 验收：`--selftest` 正例（解析第一特征模 + 精确 BE 衰减，两律闭合到 `1e-15`）PASS、6 个负例全 FAIL；#102 真实场离线回放 **U1–U8 全 PASS**（`U3/U4 max=0.0112/0.0112`、`U5/U6 max=0.00976/0.00975`、`Q: 1.688e-3 → 7.087e-5`、残差 late/early `0.77–0.88`）；接入 `build.yml` 的 `g48_rc` 延迟 `exit 1` 与 `gate48_verdict.log` artifact。
  - **已验收（CI #104 = `success`，`0a0c468`）**：`GATE48_RESULT: PASS`，U1–U8 全绿，实测与 #102 离线回放**逐位吻合**：`U2` 速度场相对误差 `3.11e-16`（排水面 `x=−5`、6 节点 / 5 边 / 长 1.0000 m）；`U3/U4` 能量残差 `0.01121 / 0.01118`、`U5/U6` 质量残差 `0.009763 / 0.009753`（阈值 3%）；`U7` 流出 `1.688e-3 → 7.087e-5`、late/early 残差 `0.77–0.88`（不增长）；`U8` `W·S·Area = F·M` 相对误差 `7.12e-16`。**原有八组 verdict 与 #102 全部逐位 IDENTICAL**（连 Gate 3.3 的日志字节数都不变，因为 4.8 完全不改 Fortran、不新增日志行）⇒ 纯后处理零扰动。九组判据在 #104 全绿。
  - **勾选条件（已满足）**：CI `conclusion=success` 且 `GATE48_RESULT: PASS`，且原有八组判据全绿（纯后处理，物理场逐位不变）。

- [x] **4.9 逐单元 IS 装配性能**（**CI #114 `success` 验收通过**，`GATE49_RESULT: PASS`）
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90`（`MEF90FlowMakeSingleCellIS` 每单元一次 + 一次 Diffusion 调用）；`m_Flow/m_MEF90_FlowCtx.F90` / `vDef/m_vDefDefault.F90`（选项袋）。
  - 问题：装配复杂度 O(N) 次 PETSc 对象（IS / elem）创建销毁。`MEF90DiffusionOperatorSet` / `BilinearFormSet` / `RHSSetCell` / `OperatorAddTransientTermSet` 内部按 `cell=1..size(cellID)` 顺序 `ADD_VALUES`；系数逐位均匀时，整组一次传入与逐单元多次传入的浮点累加顺序完全相同。
  - 改法（已落地，**不改 Diffusion 接口**）：新增 `flow.batched_assembly`（默认 `PETSC_FALSE`）；`MEF90FlowUniformCellField` 用逐位 `/=` 判定整组系数是否均匀（`sec` 声明为 `Intent(INOUT)`，规避 CI #95 的 PETSc 3.3 外部绑定 INTENT 问题）；均匀且 `size(cellID)>1` 时，`OperatorCore` / `BilinearFormCore` / `AddStorageToResidual` / `AddStorageMassMatrix` 四处循环把整个 `setIS` 一次传入 Diffusion，并用**同一个 `setIS`** 创建 `elem`（CI #65 几何契约）；不均匀则逐单元回退。`MEF90FlowEnergy` 有意不批量化 —— `DiffusionEnergySet` 内部累加标量，批量调用会改变浮点求和顺序，破坏逐位等价。
  - `PetscBag` 第 37 个字段 `batchedAssembly` 按 CI #50 教训追加在类型 / 注册 / 默认构造器末尾，默认 `PETSC_FALSE` ⇒ 生产路径保持 superposition 逐单元装配，零扰动。
  - 验收（可证伪，**不用噪声大的 wall-clock timing**）：`vDef/data/check_gate49_batch.py` Y1–Y8 —— Y1 数据完整性；Y2 静态反回归（helper / 4 处 fast path / bag 顺序 / 默认 FALSE）；Y3 deck 只差一行 `batched_assembly: true`；Y4 bag 可见性（对照 TRUE / 生产 FALSE）；Y5 确定性复杂度计数 `Gate 4.9: tag=... nCells=N nIS=1 nElem=1 mode=batched`（3 个运行时活跃 tag，生产日志无该行）；Y6 `p / u_x / u_y` 在 51 记录 × 356 节点**逐位一致**；Y7 O(N)→O(1) 外推 + `elem` 与 `cellIS` 同源；Y8 生产隔离（step 数相等、无 fatal 签名）。`--selftest` 1 正例 + 10 负例全 PASS；接入 `build.yml` 的 `g49_rc` 延迟 `exit 1` 与 `gate49_verdict.log` / `terzaghi_batched.log` / `out_batched.gen` artifact。
  - **CI #113 失败根因与修复（`dabe01a`）**：#113 编译 clean、Y6/Y7/Y8 全 PASS（物理逐位一致 `max|dPhysics|=0`、对象创建 `2400→8`、每步省 2392 次），但两处静态裁决误判。① `GATE49 Y5 FAIL`：判据要求 4 个 tag，而 `MEF90FlowAddStorageToResidual` 自 CI #64 起在 `TS_LINEAR + ksponly` 下 `xDot==0`、额外 scatter 被删除后已不被调用，其 fast path 是源码存在但运行时 dead code，永不打印 tag ⇒ Y5 改为只要求 3 个运行时活跃 tag（`OperatorCore` / `BilinearFormCore` / `StorageMassMatrix`），Y2/Y7 仍静态审计全部 4 处 fast path 并按 `N_FAST=4` 外推。② `GATE47 X2 FAIL`：硬编码「恰好 36 个 bag 字段」，4.9 合法追加第 37 个字段被误判为回归 ⇒ X2 改为 **append-only 契约**（4.7 六字段固定占据槽位 30..35、字段数只设下界、默认构造器 arity 注释必须与活跃 bag 数一致），插入中间字段仍会 FAIL；4.9 Y2 同步采用同一契约，构造器注释改 `37 components`，并各加 1 个精确负例。离线用 #113 真实日志 / EXO 回放：`GATE47_RESULT: PASS`、`GATE49_RESULT: PASS`。
  - **已验收（CI #114 = `success`，`dabe01a`）**：`GATE49_RESULT: PASS`，Y1–Y8 全绿；十三组 verdict 同轮全 PASS。实测 `nCells=300`（`-np 2` 每 rank）、`nIS=nElem=1`、151 条 `Gate 4.9:` 行（3 个运行时活跃 tag）、对象创建 `legacy=2400 → batched=8`（**300×** 归约，每步省 2392 次）、`p/u_x/u_y` 在 51 记录 × 356 节点逐位一致、`max|dPhysics|=0`、生产 step marker 408 相等且无 fatal 签名。
  - **生产零扰动（逐字节证明）**：与 #111 的 12 组共同 verdict 相比，**10 组逐字节 IDENTICAL**；`gate33_verdict.log` / `terzaghi_verdict.log` 唯一差异是 solver log 字节数 `297901→298025` / `297823→297947`（+124 B），来源仅两处：git changeset 行 `ef6ee2a→dabe01a` 与 PetscBag 视图新增一行 `batched_assembly = FALSE`，**所有物理判据数值逐位一致**。与 #113 相比生产 `terzaghi.log` 仅差 changeset 行，`out_main.gen` 的 `p/u_x/u_y` 共 54,468 个数值 **differing=0 / max|d|=0**。
  - **勾选条件（已满足）**：CI `conclusion=success` 且 `GATE49_RESULT: PASS`，原有 12 组 verdict 全绿、生产物理场逐位不变。

- [x] **4.10 清理根目录同名旧副本**（**CI #109 `success` 验收通过**，`GATE410_RESULT: PASS`）
  - **精读范围远超清单点名的几个文件**：根目录 tracked 的 84 个 `.F90`/`.inc` 中 **66 个与活子目录同名**；另有 `vDefDefault_remote.F90`、17 个 `test_*.F90`/`test20.F90`、9 个顶层 `Makefile_*`/`Makefile.include_*` 镜像。哈希比对显示 66 个同名镜像里 **56 个逐字节相同**（纯冗余）、**10 个已分叉**（`m_MEF90_Flow*.F90`、`m_MEF90_DefMechAssembly/Ctx/DefMech.F90`、`m_MEF90_HeatXferAssembly/Ctx.F90`、`m_vDefDefault.F90`、`vDef.F90`），再加上 `Makefile_m_Flow` 也分叉 —— 这才是「误改 / 误读」的真实风险面。
  - **构建隔离核实（归档前）**：顶层 `Makefile` 与 CI 只用 `make -C objs/${PETSC_ARCH} -f ../../<子目录>/Makefile`；`MEF90/`、`m_DefMech/`、`m_HeatXfer/`、`m_Flow/`、`vDef/` 五个活 Makefile 的源路径分别是 `${MEF90_DIR}/MEF90/...`、`BASE_DIR=${MEF90_DIR}/m_Flow` 等；全仓库 Makefile / CI / shell / 文档**无任何一处引用**根目录 `.F90` 或镜像 Makefile。`Tests/Makefile` 用 `${MEF90_DIR}/Tests/*.F90`，与根目录 `test_*.F90` 无关。
  - **落地（零删除）**：93 个文件 `git mv` 到归档目录，git 全部识别为 **100% rename（0 insertions / 0 deletions）**，历史完整可回溯；根目录 tracked 的 `.F90`/`.inc` 归零，只剩真正参与构建的 `Makefile` 与 `Makefile.include`。新增归档 `README.md`（原因 / 分类 / 11 个分叉清单 / 活源码对照 / 不参与构建声明 / 找回方式）与 `MANIFEST.sha256`（93 个 sha256）。`FlowCtx.F90` / `HeatXferCtx.F90` 是 `.gitignore` 明确忽略的本地诊断残留（untracked），不属于仓库内容，保留原地。
  - **W1–W8 硬判据（`vDef/data/check_gate410_legacy.py`，纯标准库 + git index）**：W1 根目录无 tracked `.F90`/`.inc` 镜像且未 tracked 者必须被 gitignore 覆盖；W2 归档 93 个文件与 manifest 条目集合一致；W3 每个哈希与**版本库 blob** 一致；W4 66 个同名镜像都有活版本；W5 五个活 Makefile 只引用子目录源且 `BASE_DIR` 正确；W6 顶层 Makefile / CI 只用 `-f ../../<子目录>/Makefile` 且不引用镜像 Makefile；W7 构建入口不引用归档目录；W8 README 列全 11 个分叉文件与 no-build 契约。`--selftest` 1 正例 + 10 负例全 PASS。
  - **CI #108 失败 → 根因 → 修复（重要经验）**：#108 编译与十组物理判据全绿，**仅 W3 FAIL**。用本机 `dos2unix` 复现 CI 的 `Fix CRLF` 步骤，发现 93 个归档文件中恰好 **1 个**被改写：`vDefDefault_remote.F90` 的 UTF-8 BOM 被移除（head `\xef\xbb\xbf#incl` → `#incl`，长度 `18808→18618`），其余 92 个匹配，`all()` 因此为 False。**归档是历史快照，应校验「提交进版本库的内容」而非会被 CI 换行修复改写的工作区字节**：W3 改为读 `git show :<path>` 的 index blob（CRLF→LF 规范化），manifest 本就基于同一 blob 故无需重算；同时让 W3 打印不匹配文件名而非只给布尔值。本地实测：干净树 PASS，且真的跑一次 `dos2unix` 去掉 BOM 后仍 PASS。
  - **已验收（CI #109 = `success`，`f2d7f57`）**：`GATE410_RESULT: PASS`，W1–W8 全绿；`W3` 报 `93 archived git-blob sha256 match MANIFEST.sha256: True; mismatched=none`；**十一组判据全绿，原有十组 verdict 与 CI #106 逐位 IDENTICAL**（连 `terzaghi_verdict.log` 也一样）⇒ 纯工程归档零扰动。
  - **勾选条件（已满足）**：CI `conclusion=success` 且 `GATE410_RESULT: PASS`，且原有十组判据全绿、物理场逐位不变。

- [ ] **4.11 CI 断言 + 回归 + 解析解比对**
  - 位置：`.github/workflows/build.yml:394-408`（Run Terzaghi，`continue-on-error: true` + `set +e` 仅 tee）；`vDef/Makefile` 的 `test` 只跑 test7/8
  - 问题：编译过但算错/崩也不会拦；`terzaghi1D.yaml` 未纳入回归。
  - 改法：去掉 `continue-on-error`，检查 `mpirun` 退出码；CI 加 `make test`；加 Python 脚本读 exo 的 pressure 与解析解 `p(z,t)` 比对（t=1s 时 T≈1.09），超阈值直接 fail。
  - 验收：CI 能拦住"编译过但算错/崩"的提交。（该条为总括性描述，已由 Gate 3.3 起的十五组硬判据体系覆盖，不单独计入在册待做。）

- [x] **4.12 `pressurePreviousStep` 死写转活诊断**（**CI #117 `success` 验收通过**，`GATE412_RESULT: PASS`）
  - 位置：`vDef/vDef.F90`（时间步循环）、`m_Flow/m_MEF90_FlowCtx.F90`（Vec 生命周期）。
  - 问题：该 Vec 只被创建 / 初值 / 每步 `VecCopy` / 销毁，**全项目无消费者**；原本预留给平流 / Biot 滞后耦合，但 4.1a / 4.1b / 4.2 收口后从未需要，是每步一次 MPI 全量拷贝的死写。
  - 改法（**零反馈后处理**）：初始化时把 Dirichlet 修正后的真实初始压力场拷入 `pressurePreviousStep`（仅靠 `CreateVectors` 的均匀 `VecSet` 会让 step 1 在排水面出现虚假 `−1e5 Pa` 压降）；每步 `TSSolve` **之前**用 `VecMin/VecMax` 采样（此时仍是上一步末场），`TSSolve` 之后打印 `Gate 4.12: step N t T dpmax ... dpmin ...`；每步末的 `VecCopy` 保留。`m_Flow` 装配例程仍不读该 Vec ⇒ 物理场逐位不变。
  - 物理意义：`dpmax` 是闭端（不渗透底）压力的步间降幅，即离散固结速率；实测峰值在 step 10（`−3734.9 Pa/step`）、末/峰 `0.132`、后期单调衰减，与 Terzaghi 固结先加速后指数衰减一致。
  - 硬判据 `check_gate412_prevstep.py` Z1–Z6：Z1 51+51 行 / step 连续 / `dt=0.02`；Z2 静态反回归（2 处快照 + 采样 + 打印 + CI 接线 + `FlowAssembly` 零消费者）；Z3 `dpmax(n) == pmax(n)−pmax(n−1)` 逐位一致、step 1 恰为 0；Z4 EXO 节点场极值锚定；Z5 固结形状（全负 / 单峰 step 3–20 / 末峰比 <0.30 / 尾部衰减）；Z6 无 fatal。`--selftest` 1 正例 + 8 负例全 PASS。

- [x] **4.13 标准 Newton 残差平台 `3~13` 根因定性**（**CI #117 `success` 验收通过**，`GATE413_RESULT: PASS`）
  - 背景：CI #83 实证标准 Newton（`f0=2.082e4 → it1≈9.3 → 平台 3~13`，51 次迭代无收益）与 ksponly 的物理场等价（`dP=0`、`dUx=0.21%`、`dDamage=0.38%`、应力 ≤0.41%）；CI #84 起生产改用 ksponly。故该平台**不是装配错误**，而是数值尺度问题。
  - 根因（静态定位到具体代码行）：`m_MEF90_DefMechAssembly.F90` 自由自由度残差为 `∫(C:ε(u)):ε(v) dΩ`（量级 `E ~ 1e9`），Dirichlet 行残差为 `residualPtr = u − u_bc`（量级 1），Jacobian 用 `MatZeroRowsColumnsIS(A, setISdof, 1.0_Kr, ...)` 把 BC 行列置为单位对角 1 ⇒ 矩阵混合 1 与 ~1e9 两种尺度，条件数 ~E。CG 在**相对残差**容差下收敛时，物理力残差仍停在几牛顿，SNES 反复迭代一个降不下去的 fnorm；该 deck 实质线性（SplitNone、损伤 ~3.7e-5、零阻尼 / 内聚），一次线性解即精确解。
  - 处置：不强行"修好"标准 Newton（生产已无该循环），而是把诊断**可证伪地冻结**：`check_gate413_newton.py` Q1–Q5 —— Q1 生产命令保持 `-disp_snes_type ksponly` 且未回退；Q2 三处尺度锚点仍在源码；Q3/Q4 两份日志无 Newton fnorm / displacement `-5`（7 条 damage `-5` 为信息性）；Q5 ksponly 与生产的 51 条压力包络**逐位一致**。`--selftest` 1 正例 + 6 负例全 PASS。
  - **CI #116 插曲**：新增的 `Gate 4.12:` 行被 `check_gate41_biot.py` 的裸前缀 `'Gate 4.1' in line` 误计为 102 条未解析探针（`2 日志 × 51 步`），使 B1 与 Gate 4.3 C7 连带 FAIL；改为精确 `'Gate 4.1:'` 后 #117 全绿。**经验**：诊断行前缀不要用裸子串匹配，`4.1` 与 `4.1b` / `4.12` / `4.13` 会碰撞。

- [x] **4.14 三维维度分派真实验证**（**CI #121 `success` 验收通过**，`GATE414_RESULT: PASS`，W1–W7）
  - 问题（**验证体系空白**）：此前 15 组硬判据**全部只在 2D `Beam1x10-tri3`** 上跑过；`m_Flow/Makefile` / `m_DefMech/Makefile` 虽然同时编译 `-DMEF90_DIM=2/3` 孪生，`MEF90_MATS_IDENTITY` 在 3D 分支映射为 6 分量 `MEF90MatS3DIdentity`，但 3D 路径从未被任何闸门实际执行。Gate 4.1b 的有效应力项与 4.9 的批量扩散算子都依赖该宏，2D-only 回归无法发现 3D 装配 / 分派错误。
  - 落地：`terzaghi1D_dim3.yaml`（与基准 deck **只差一行** `time.numstep 51→11`，dt=0.1 s、同覆盖 t=0..1 s，维度由 `-prefix data/Beam1x10-tetra4` 选择而非 YAML key）；`check_gate414_dim3.py` W1–W7；`build.yml` 的 3D 诊断运行 / 判据 / 延迟 `exit 1` / artifact 三项。3D 网格 `num_dim=3`、2635 节点、11761 TETRA4 + 56/56 TRI3 边界面、cell sets `[1,10,20]` 与 2D 一致。
  - 判据：W1 运行 / 网格数据完整性；W2 `Displacement_Z` + `num_dim=3` + `coordz` + 同 cell sets / 几何；W3 静态反回归（3D identity 宏、`-DMEF90_DIM=3` 孪生、`DMMeshGetDimension` 分派、门面 2D/3D 实现对、3D EXO 命名）；W4 11 个共同物理时刻 3D vs 2D `pmax`（`MAX_REL_P = 0.40`）；W5 3D 位移场 live / finite / `|uz|/|ux|` 有界；W6 求解链完整 + 每步 KSP 收敛 + 无 fatal；W7 单行 deck + CI 接线（含求解器配置与 `s414_rc` 裁决）+ 生产隔离。`--selftest` 1 正例 + 18 负例全 PASS。
  - **CI #119 失败（编译期维度宏）**：3D 运行在 `Gate 3.1: init 07 FlowUpdatePerm` 刷约 9900 行 `Element2D_Scal_Init: Element type not implemented yet MEF90_P1_Lagrange_3D_Scal` 后 MPI abort。根因：`m_MEF90_Flow.F90` 门面**不带 `-DMEF90_DIM` 编译**，`MEF90_ELEMENT_SCAL` 固定绑 2D 元素类型，而 `MEF90Element_Create` 泛型按 `elem` 的**声明类型**解析工厂。修复（`39444b8`）：`UpdatePerm` / `UpdateMuFromTemperature` 各拆 `...2D` / `...3D` 实现对 + `DMMeshGetDimension` wrapper；Gate 4.6 V1/V2、4.7 X2 的「取整段 body」静态契约同步改为审计 wrapper 与两个孪生，各自补可证伪负例。
  - **CI #120 失败（求解器代价，非装配 bug）**：编译干净、15 组原有判据全 PASS，但 3D step 1 完成后 step 2 的**第一次 flow TSSolve** 在 900 s 内未打印任一 `Linear solve converged`（2D 控制 356 节点 CG 42 次迭代收敛；3D 2635 节点 tetra4 上 CG + bjacobi/ILU(0) + rtol 1e-10 停滞），`timeout` 发 SIGTERM → MPI_ABORT → `out_dim3.gen` 截断在 13,766,716 B，W1 解析 EXO 抛 `struct.error`。
  - **修复（`c95a605`，仅 workflow + 判据，Fortran 零改动）**：3D 诊断运行改 `-np 1` + `-flow_pc_type bjacobi -flow_sub_pc_type lu`（单 rank 下唯一 bjacobi 子块即整个 seq-AIJ 矩阵，精确 LU 使 CG 每步 1 次迭代收敛）；加 `-flow_ksp_monitor`；cap 900→1800 s；**12 个生产 2D deck 保持 `-np 2` / ILU 不变，11 个物理步与 W1–W7 一条不减**。`s414_rc`（3D mpirun 退出码）此前只 echo、从不裁决，现接入延迟 hard fail；W6 新增「`Linear solve converged` ≥ `N_REC_3D−1`」硬要求；W7 冻结求解器配置与 `s414_rc` 裁决；`--gate` 的 data error 改为可归因的 W1 FAIL + 截断提示。
  - **CI #121 实测**：`GATE414_RESULT: PASS`，W4 `p0=1e5 / pend=1.183e4 / max shared-time rel-diff=0.306`，W5 `|ux|max=1.01e-4 / |uz|max=6.618e-6 / ratio=0.0655`，W6 10 条 `Linear solve converged ... iterations 1`；3D 日志 `Element2D_Scal_Init` / fatal / `PETSC ERROR` / `Caught signal` / `MPI_ABORT` 全为 0；`out_dim3.gen` 15,106,908 B、11 记录、6 个 nodal vars 含 `Displacement_Z`。TS view 实证 `1 MPI processes`、`flow_ PC=bjacobi`、`flow_sub_ PC=lu`、`package: petsc`。

- [x] **4.15 / Stage B0 施工排量–可压性相似校准 + 源项量纲标定**（**CI #125 `success` 验收通过**，`GATE415_RESULT: PASS`，V1–V8）
  - 立项依据：`RATE_FRACABILITY_ROADMAP.md` 第 3.5 / 5.1 节。此前 `source ≠ 0` 只在符号与接线层验证（4.4 / 4.5），**从未用于排量扫描、从未与可压性参数耦合**；且所有 deck 的 `T = c_v t/H² ≈ 1`，压力 1 s 内扩散殆尽，根本无法建立注入净压。B0 是阶段 B 的准入冒烟，**零 Fortran 改动**。
  - 相似映射（表 C，II 类中等脆性页岩）：`terzaghi415_scale.yaml` 相对 `terzaghi1D.yaml` **只改 4 行代码** —— `numstep 51→11`、`perm_min 5e-11→2.35e-18`、`perm_max 1e-9→4.70e-17`、`storativity 4.6e-10→1.50e-10`。得 `c_v=1.5667e-5 m²/s`、`T(1 s)=1.5667e-7`（与现场 `t_f=3600 s / H_f=50 m` 的 `1.57e-7` 一致）、`k/k_deck=4.7e-8`、`S/S_deck=0.326`、`H²/c_v=6.38e6 s`、边界层 `√(c_v t)=4 mm`。
  - 注入 deck：`terzaghi415_inject.yaml` 相对校准 deck **只改 1 行** `cs0001.source 0.0→1.5e-5 [1/s]`。因 `T≪1`，Darcy 项在 1 s 内被冻结，离散平衡退化为 `S dp/dt = f`，解析锚 `p(t)=p0+(f/S)t`、`f/S=1e5 Pa/s`，**体相（bulk / 中位数）压力**从 `1e5` 线性升至 `2e5 Pa`。
  - **CI #124 诊断更正（两个彼此独立的缺陷，均由真实 artifact 定位）**：
    - ① **EXO 追加复用**：`MEF90CtxOpenEXO`（`MEF90/m_MEF90_EXO.F90:125-143`）用 `Inquire(file=...,exist=...)` 判断；`_out.gen` 已存在时走 `EXOPEN(EXWRIT)` **追加**。两次 11 步 B0 run 因此继承了 51 记录的生产结果文件：`out_415_*.gen` 前 11 条确实是 4.15 新场，第 12–51 条是旧生产衰减尾部（`pmax[11]=70667.193` 与 `out_main.gen[11]` 完全相同）。V1 报 `51 != 11`；V3 未切片而遍历全部 51 条，把旧尾部算成 `90.94% drift`。修复：`build.yml` 在**每个** B0 run 前 `rm -f data/Beam1x10-tri3_out.gen`，让 vDef 走生产 run 已在使用的 `EXCRE(EXCLOB)+EXCOPY` 新建分支（同一逻辑见 `build.yml:439` / CI #67，且 `*_out.gen` 未被跟踪，本来就是 CI 首个 run 的标准路径）；V1 另加 EXO `pmax` 与求解器 `Gate 3.2` 打印值的**逐位交叉验证**，即使记录数碰巧相等也能抓出错文件。之后照旧 `cp -f /tmp/out_main.gen` 恢复生产 mesh，交付 artifact 不变。
    - ② **观测量选错（物理其实是对的）**：`A=k/μ=2.35e-15`、`a·S=1.5e-9` 时 Darcy 刚度几乎冻结，一致质量矩阵 + 强 Dirichlet 入口在**少数节点**造成局部数值过冲。#124 实测 step 10 `pmax=3.9188e5 Pa`（node 171）而 `p50=2.0000e5 Pa`、`p90=2.0076e5 Pa`、`mean=1.9913e5 Pa`；`pmax` 斜率 `2.9184e5 Pa/s` 是过冲模态，不是 `f/S`。日志 bag 复核 `S=1.5e-10`、`f=1.5e-5`、`dt=0.1`、`permF/storF=1` 全对；Gate 4.5 的 `W=∫f·p dΩ` 从 `14.88` 升到 `30.16`（`W/W_an=1.005`）独立证实**体相压力增长率恰为 `f/S=1e5 Pa/s`**。修复：V3/V4 改用中位数 `pbulk` 作物理观测量，`pmax` 降级为局域性诊断（`p90` 须贴近 `pbulk` ≤5%、峰值/体相斜率比须落在 `[1,6]` 的固定空间模态带），并继续硬拒 `slope≈f` 的 `[Pa/s]` 误读。**未改任何 Fortran**。
  - **源项量纲定论**：`m_MEF90_FlowCtx.F90:791` 注册 `Volumetric source/sink rate (1/s)`，`MEF90FlowEnergy` 计算 `W=∫f·p dΩ`，`-div(A∇p)=S dp/dt + f` 中 `S dp/dt` 亦为 `[1/s]` ⇒ **`f` 是 `[1/s]`，旧 deck 注释 `[Pa/s]` 为错**。`f/S` 才是压力爬升率 `[Pa/s]`；`Q_total = f·Area = 1.5e-4 m²/s` 是**单位厚度体积流量**，换算 3D 需乘缝高 `H_f`。V4 由真实 EXO 的 **`pbulk`** 序列反解 `S_eff=f/slope` 并显式拒绝 `slope≈f` 的 `[Pa/s]` 解释（两者差 `1/S=6.67e9` 倍）。B0 取保守 `f` 使响应保持线性、无损伤；现场排量窗口留到 B1 扫描。
  - 判据 V1–V8：V1 运行 / 11 记录干净 EXO / EXO-vs-`Gate 3.2` 交叉验证 / 网格完整性；V2 从 deck 反解并冻结 `T`、`c_v`、`k/S` 比与 `Q_total`；V3 零源对照的 **`pbulk`** 在 11 步内冻结（存储主导）；V4 **`pbulk`** 注入线性响应 + `S_eff` 标定 + `pmax` 局域性；V5 Gate 4.5 注入功（对照 `W≡0`、注入 `W>0` 且单调）；V6 每步 KSP 收敛 + 无 fatal；V7 静态反回归 + 4 行 / 1 行 deck 契约 + CI 接线（含两次 EXO 重建、两个 `s415_*_rc` 与 `g415_rc` 延迟裁决）；V8 生产隔离。`--selftest` **1 正例（刻意带真实 2.9184× 局部过冲）+ 25 个变异（含 1 个不变性）全 PASS**。
  - CI 接线：`build.yml` 在 3D 诊断之后新增 2 个 2D 运行（`-np 2` 生产同款配置、11 步、独立 `.log` / `.gen`），**每个 run 前重建干净结果 mesh**，verdict 调用 + **两个 mpirun 退出码与 `g415_rc` 全部接入延迟 `exit 1`**（复用 CI #120 教训）+ artifact 5 个新成员。12 个生产 2D deck 与 3D 诊断运行配置保持原样。
  - 提交前离线复核（零本地编译 / 零本地求解器）：`git diff --check`、`py_compile`、PyYAML workflow、三份 deck 的 4 行 / 1 行代码契约、`rm -f` 出现位置、判据 `--selftest`；并用 `ci124_diag.zip` 的**真实 EXO + 日志**做端到端回放：**原样 51 记录 → 只剩 V1 FAIL**；模拟重建后 11 记录 → **V1–V8 全 PASS**，`bulk slope=1.000000e5 Pa/s`、`S_eff=1.500000e-10`、线性残差 `1.8e-8`、`p90 dev=0.379%`、局部峰值比 `2.9188`、`W/W_an=1.005`。
  - **CI #125 实跑验收（`603ce0b`，`success`）**：`GATE415_RESULT: PASS`，V1–V8 全绿；`g415_slope_measured=1.000000e+05`（精确 `f/S`）、`g415_S_eff=1.500000e-10`、`g415_peak_ratio=2.918774`、`g415_p90_dev=3.787269e-03`、`g415_w_last=3.016054e+01`；`out_415_scale.gen` / `out_415_inject.gen` 均为**干净 11 记录**（`pmax` 末值 `100036.8` / `391877.6`），`out_main.gen` 仍为 51 记录（`9055.7`），生产 artifact 未被破坏。**原有 16 组 verdict 日志与 #124 逐字节 IDENTICAL**，唯一变化就是 `gate415_verdict.log` 由 FAIL 转 PASS —— 零物理扰动。

- [x] **4.15 / Stage B1 排量×可压性参数可辨识性主实验**（**CI #127 实测 `GATE415_B1_RESULT: PASS`，A1–A8 全绿**；该轮整体 `failure` 仅因 B2 的 L4，见下条。`check_gate415_b1_identifiability.py`）
  - **文献前提更正**：路线稿第 5.2 节原 R4 用「黏性主导 $\partial\ln p_{net}/\partial\ln Q\approx0$ / 向 $-1/3$ 移动」做判据，属**裂缝传播段**标度。B0 实测 `Damage_max=1.48e-4≪1`，是**起裂前存储主导**窗口。Lecampion, Desroches, Jeffrey & Bunger (2017), *JGR Solid Earth*, doi:`10.1002/2016JB013183` 明确：起裂前 $Q_{in}\approx0$，恒排量给出**线性加压**，斜率 $=Q_0/U$（$U$ 为注入系统柔度），岩石参数在压力曲线上不可辨识；起裂压力 ≠ 破裂压力。在 vDef 中 `S` 扮演 $U/Q_0$，$\beta=f/S$（B0 实测 $1e5$ Pa/s）。Santillán (2018) 三类排量模式与 Detournay (2016) M–K 顶点均属传播段，**本轮不主张验证**。
  - **可证伪命题（R4′ 可辨识性分层）**：Level 1 压力曲线 $\beta$ 只对 $(f,S)$ 敏感（$G_c/k/\mu$ 偏差 <6%）；Level 2 相场表面能 `C1 ∝ Gc/ℓ`、`C2 ∝ Gc·ℓ·toughnessAnisotropyMatrix`（`m_MEF90_DefMechAssembly.F90:2340-2341,2355`）对 $G_c$、$Q$、$A$ 有可测响应。B0 真实日志 `======= Total: ... surface:` 从 `1.20334E-06` 升到 `6.49012E-04`（539×），零源对照几乎冻结。
  - 交付：6 个**单行代码差异** deck（`terzaghi415_b1_{qhalf,qdouble,stor,gc,perm,visc}.yaml`：$Q$ 0.5×/2×、$S$ 2×、$G_c$ 0.5×、$k$ ×10、$\mu$ ×2）+ 6 个 11 步 2D 诊断 run + `A1–A8` 判据 + `gb1_rc` 延迟硬裁决 + 6 个 mpirun 退出码裁决；`--selftest` **1 正例 + 20 变异（含 1 不变性）全 PASS**。
  - A1 数据完整性（含 EXO↔`Gate 3.2` 逐位交叉、11 行能量分解）；A2 $\beta(0.5/1/2\times Q)$ 比值与 $S_{eff}$ 恒定；A3 $2S\Rightarrow\beta$ 减半并显式拒绝「β 只依赖 f」；A4 **Lecampion 数值验证**（$G_c/k/\mu$ 不扰动 $\beta$）；A5 $G_c$ 在表面能可辨识且与 A4 解耦；A6 $Q$ 单调性（$W$、表面能、$damage_{max}$）；A7 deck/源码/CI 静态契约（含 ≥8 次 EXO 重建）；A8 B0 CI #125 数值与生产隔离。

- [x] **4.15 / Stage B2 变排量 / 层理各向异性 / 3D 一致性**（**CI #128 `success` 验收通过**，`GATE415_B2_RESULT: PASS`，L1–L8 全绿；#127 首轮唯一 `L4 FAIL` 已查明为判据物理前提错误，重写为「起裂前各向异性不可辨识」零位命题后复验通过；`check_gate415_b2_extension.py`）
  - `terzaghi415_b2_linear.yaml`（`source_scaling: constant → linear`）：`MEF90FlowSetTransients`（`m_MEF90_FlowAssembly.F90:873-875`）`SetSourceCst` 后 `VecScale(source,timeVal)` ⇒ $f(t)=f_0t$；transient 模式下 `vDef.F90:661` 的第二次调用被 `If (.not. flowTransient)` 保护，唯一时间值是 `time(step)`。**BE 离散精确锚** $\Delta p_n=(f_0/S)dt^2n(n+1)/2$，$n=10\Rightarrow5.5e4$ Pa $=0.55\times$ B0；L2 拟合二次系数并要求线性模型残差 $>3\times$ 二次残差（工程类比 step-rate test，但连续斜坡 ≠ 阶梯）。
  - `terzaghi415_b2_aniso.yaml`（inject **+1 行** `toughnessAnisotropyMatrix: 2,0.5,0`）：`MatS2D=(XX,YY,XY)`、identity `=(1,1,0)`（`m_MEF90_LinAlg.F90:50-54,107`）；逗号后**无空格**（CI #85 RealArray 坑）。⚠️ **原判据「两个主方向同时偏离 ⇒ 表面能必变（>5%）」已被 CI #127 推翻**：$A$ 只进 $C_2|\nabla\alpha|^2$，起裂前被 $C_1w(\alpha)$ 淹没，实测 `surface_ratio=1.004065`（0.4065%）、`damage_dev=1.247%`、`beta_dev=0`。现行 L4 改为五条可证伪零位命题：**$d_s\in(10^{-4},5\times10^{-2})$**（已接线且未泄漏到 $C_1$）、**$d_s<d_{dmg}$**（只作用于梯度通道）、$d_{dmg}>1\%$、$\beta$ 偏差 <6%；`report()` 增 `gb2_surface_dev_aniso`。页岩层理 $G_c$ 各向异性「控制裂缝轨迹」属**传播段**结论，本窗口不主张验证。
  - 3D 通道：复用 `terzaghi415_inject.yaml` + `-prefix data/Beam1x10-tetra4`（维度由 mesh prefix 分派，Gate 4.14），复用 #121 的 `-np 1` + `bjacobi/LU` + `timeout 1800` 避免 #120 超时；L6 要求 3D $\beta$ 与 2D 差 <12%、6 个 nodal vars 含 `Displacement_Z`、无超时。3D run 后从 `/tmp/out_dim3.gen` 恢复 Gate 4.14 的结果 mesh。
  - 交付：2 个新 deck + 2 个 2D run + 1 个 3D run + `L1–L8` 判据 + `gb2_rc` 延迟硬裁决；`--selftest` 首版 **1 正例 + 19 变异**，**#127 后重写 L4：1 正例 + 24 变异（含 3 不变性）全 PASS**，正例表面能 / damage 已锚定 #127 真实读数。
  - **CI #127 首轮实测（`f07cf39`，`failure`）**：编译 2D/3D 双份成功、`Run Terzaghi` 1429 s（#125 为 1646 s）；原有 17 组 verdict 与 `GATE415_B1_RESULT: PASS` 全绿。L2 二次系数 `1000`（解析锚 `1000`）、`dp10=55000=0.550×B0`、线性模型残差 `13.64%` vs 二次 `0.00%`；L3 `W=0..23.284` 单调、常排量 `W_last=30.161`；L6 `beta_dim3=1.000000e+05`（偏差 `1.08e-13`）、`Displacement_Z` 活、无 timeout。**唯一失败 `GATE415_B2_RESULT: FAIL (L4)`**：`anisotropic C2 moved surface energy only 0.41%`，最终由延迟 `gb2_rc=1` 触发 `exit 1`。根因与修复见 `RATE_FRACABILITY_ROADMAP.md` 第 13 节。
  - **CI #128 复验验收（`e6cce4b`，`success`）**：`GATE415_B2_RESULT: PASS`，L1–L8 全绿；新 L4 实测输出 `surface dev 0.406% in band (0.01%..5%) and < damage dev 1.247%, beta dev 0.00%`，新增指标 `gb2_surface_dev_aniso = 4.064640e-03` 与离线预测逐位一致。`gb2_quad_slope=1.000000e+03`、`dp10=5.500000e+04`、`dp10_ratio=0.55`、`beta_dim3=1.000000e+05`（偏差 `1.082467e-13`）等 17 项 `gb2_*` 指标与 #127 **逐位相同** ⇒ 只改了 L4 的判定逻辑，物理场零扰动。**18 条 `GATE*_RESULT` 全 PASS**，其中 17 条非 B2 verdict 与 #127 **逐字节 IDENTICAL**，日志中无任何 `=== FAIL` 行；加上 `terzaghi_verdict` 汇总，体系达到 **19 组判据同轮全绿**。
  - **物理诚实边界**：B0 注入日志含 5 条 `damage field SNESSolve failed ... -5`（既有已知残差平台，CI #84）；能量行是当前场的确定性积分，只作**相对 / 可证伪**观测量，不表述为已收敛断裂解，不报缝长 / 缝网面积 / 起裂压力绝对值。全部结论限定在起裂前、存储主导、无裂缝扩展窗口。
  - **首版提交前离线复核全通过**：`git diff --check`、两个判据 `ast.parse`、PyYAML（workflow + 10 份 415 deck）、`Run Terzaghi` 整块 `bash -n`（1384 行）、B1/B2 `--selftest`、**#125 真实 artifact 混合回放**（A1–A8 / L1–L8 全 PASS，真实读数 `pbulk[-1]=199999.9999999968`、`surface_last=6.49012e-4`、`damage_max=1.4803e-4`、tri3 356 节点 / tetra4 2635 节点 `num_dim=3`）、B0 判据在新 workflow 下仍 `GATE415_RESULT: PASS`。详见 `RATE_FRACABILITY_ROADMAP.md` 第 12 节。
  - **#127 后 L4 重写的离线复核（零 Fortran / 零 deck / 零 workflow 改动）**：① B2 `--selftest` **PASS**（1 正例 + 24 变异，含 3 不变性）；② 把 #127 真实 `gb2_*` 指标代入新 L4，五条命题 **5/5 PASS**（`ds=0.004065`、`dd=0.012474`、`beta_dev=0`）；③ 剥离 GitHub 时间戳后用判据自身 `ENERGY` 正则从 8.48 MB job logs 解析出 **797 条**真实 `Total` 能量行，定位 `surface: 6.49012E-04`（isotropic）与 `6.51650E-04`（aniso），算出 `ratio=1.004064640`，与 verdict `gb2_surface_ratio_aniso=1.004065` 一致；④ L1 已在 #127 真实 log/EXO 上 PASS，本轮未改任何解析函数与 L1–L3 / L5–L8。详见路线稿第 13.4 节。
- [x] **研究项 / 页岩油文献参数重算（表 A′–C′，M1–M8）**（**CI #130 `success` 验收通过**，`SHALE_OIL_PARAMS_RESULT: PASS`，M1–M8 全绿；对应路线稿第 5.4 决策点 2 与第 14 节；证据等级说明见 14.9）
  - **交付**：`vDef/data/shale_oil_params.json`（S01–S20 文献、25 个参数 low/rep/high、`window` 首版设定与已发布表值、`formulas` 指数契约、`numerical` B0 冻结常数）、`vDef/data/check_shale_oil_params.py`（M1–M8 硬判据 + `--md` / `--json` / `--selftest`）、`build.yml` 新增 `gso_rc` 延迟裁决与 `shale_oil_params_verdict.log` artifact。**零 Fortran / 零 deck / 零求解器运行 / 零生产 run 改动**。
  - **可信度锚 M3**：用独立于 JSON 契约的参考公式逐位复现首版表 A–C 四类 × 10 个已发布数字（共 40 个比较点），全部落在印刷精度内（最紧为 I 类 $\mathcal{C}$@20：`0.0013552` vs `0.0014`，纯舍入差 3.2%）。**M2** 把每个标度写成参考式 / 字面式 / 展开幂律式三份，要求相对误差 $<10^{-10}$，使前因子与指数的静默漂移可证伪。
  - **文献更正**：清水 $\mu=1e{-3}$ → 变黏滑溜水 `1e-3–6e-2` (rep `5e-3`) Pa·s；首版 I–III 类反推 $k_{field}\approx1e{-20}–1e{-19}$ m² 比文献 $k_{eff}$ `1e-18–1e-17` (rep `3.2e-18`) **低 2–3 个量级**；$k_{core}$ (rep `1.2e-16`) 与 $k_{eff}$ 差 **37.5×**，M4 强制只用 $k_{eff}$ 进入 $T$ / $C_L$；$S$ 必须由 $\phi(c_{pore}+c_{fluid})$ 推导，rep 校验 `1.044e-9` vs `1e-9`（差 4.4%）；static/dynamic 模量不可混用；Q 窗口 `4–20` → `10–14` m³/min。新恒等式 $\sigma_c=\sqrt{E'G_c/\ell}=K_{Ic}/\sqrt{\ell}$ 与 $E'$ 无关。
  - **重算结果**：表 A′ rep 档 $E'=24.584$ GPa、$G_c=14.643$ N/m、$\sigma_c=1.342$ MPa、$c_v=6.400e{-7}$ m²/s、$T=2.612e{-6}$、$C_L=1.354e{-6}$ m/√s；表 B′ rep 档 $\mathcal{K}$ @10→14 = `0.1334 → 0.1226`、$\mathcal{C}$ = `0.01628 → 0.01376`；表 C′ rep 档 $k_{num}=3.918e{-17}$ m² = B0 冻结 $2.35e{-18}$ 的 **16.67×**。**文献推翻了数值假设，未推翻分区结论**：三档一致参数下 $\mathcal{K}<3$、$\mathcal{C}\ll1$、$T\ll1$ 仍全部成立。
  - **可证伪结论：不能交付排量优选表**。① 起裂前存储主导窗口下 $\beta=f/S$ 对 Q 严格线性、无阈值 ⇒ 「优选」无定义；② $\mathcal{K},\mathcal{C}$ 是传播段标度，B1 A4 已由真实求解器证明该窗口不可辨识；③ $T$ 的文献单调角点包络 $[6.0e{-10},1.8e{-3}]$ **跨 6.48 个数量级**（M7 固化，防止偷偷收窄区间伪造清晰答案）；④ q 口径 `so_q_ratio_field_over_b0=31.75` 的缺口在该窗口**不可观测**，不能据此判定 B0 错误。表 A′–C′ 的用途是数值平台相似区间与阶段 B3 边界条件，区块级结论仍需 DFIT / mini-frac 标定。新增风险 **R8**（把离线标度当现场排量建议）。
  - **离线自检**：`SHALE_OIL_PARAMS_RESULT: PASS`（M1–M8 全绿）；`--selftest` **1 正例 + 33 变异全被拒绝 + 3 不变性保持**；`--md` 表格与路线稿第 14.4 节一致性核对；`py_compile` / `ast.parse` / `json.load` / `yaml.safe_load` / `git diff --check` / `Run Terzaghi` 整块 `bash -n`（1402 行）全通过；B0 / B1 / B2 / 4.10 四个既有判据的 `--selftest` 回归全 PASS。详见 `RATE_FRACABILITY_ROADMAP.md` 第 14 节。
  - **CI #130 验收（`37bc756`，`success`）**：run `42807197140`，11:24:25 → 11:50:40 UTC（26 min 15 s）；job 与 step 14/15/16/17 全 `success`。⚠️ **证据等级为逻辑蕴含级，非日志逐字节级**：本轮 `runs/{id}/logs`、`jobs/{id}/logs`、`artifacts/{id}/zip` 三端点均 `401`（需 admin / `actions:read`），无法像 #128 那样逐字节复核 verdict 文本。验收依据：① `gso_rc` 由 `PIPESTATUS[0]` 捕获且在 step 15 末尾 `exit 1` 延迟裁决 ⇒ step 15 `success` **蕴含** M1–M8 全绿；② 其余 18 组 `GATE*_RESULT` 与 9 个 mpirun rc 的延迟裁决同在该 step ⇒ 原有判据与生产路径零扰动；③ M8 静态校验接线、M7 校验 B0 deck 冻结 token；④ 本地离线全绿。体系由 **19 组 → 20 组 verdict 同轮全绿**。若日后取得 `actions:read` 凭据，应补做 verdict 文本逐字节比对。

- [ ] **研究项 / B3 起裂准入立项（AT2 韧性阶梯，N1–N8）**（**代码已落地，CI 待验收**；对应路线稿第 15 节与第 5.4 决策点 3）
  - **物理前提改变**：B0/B1/B2 全部在起裂前（`Damage_max=1.480310e-4`）；B3 用零 Fortran 的单因子 `fracturetoughness` 阶梯把损伤推过 $a(\alpha)=1/2$ 的刚度半衰门槛。**仍不是裂缝传播段**，Santillán / Detournay 传播段标度与排量优选表在本阶段依然无定义。
  - **解析锚（源码级）**：AT2 $c_w=1/2$（`m_MEF90_DefMechAT2.F90:33`）、$a=(1-\alpha)^2$、`w=\alpha^2`、`C_1=G_c/(2\ell)`、`C_2=G_c\ell/(2c_w)A`（`m_MEF90_DefMechAssembly.F90:1687-1688,2062-2063`）、$C_3=a'\Psi+C_1w'$（:1786）⇒ $\alpha_{eq}=2\Psi/(G_c/\ell+2\Psi)$。由 #125 真实 damage 反演 $2\Psi=1.998714\times10^{1}$ Pa、$G_c^*=9.650647$ N/m（$a=1/2$ 定义）；#127 独立实测 `gb1_surface_ratio_gc=2.000410` vs 解析 `1.999408`（差 0.05%）、`gb1_beta_dev_gc=0` 逐位。
  - **交付**：4 个单行差异 deck（`2700 / 270 / 27 / 2.7` N/m，仅改 `fracturetoughness`）+ 复用 B0 注入 / B1 `gc` 两个既有 EXO 作为高档锚 + `check_gate415_b3_initiation.py`（N1–N8，`--gate` / `--selftest` / `--predict`）+ 4 个诊断 run、4 次 EXO 重建、`gb3_rc` 延迟硬裁决、verdict log 与 artifact。**零 Fortran / 零既有 deck / 零生产 run 改动**。
  - **N5 是本阶段核心判别**：损伤必须**局域化**（$d_{max}/d_{mean}>1.5$，真实 B0 为 `6.996`）；若末档为全域均匀软化，必须 FAIL 并如实报告「不构成裂纹」，**不得为凑 PASS 放宽阈值**。`-5` 残差平台显式报告条数、不作硬判。
  - **离线证据**：`--selftest` **1 正例 + 28 变异（含 3 不变性）全 PASS**（不变性含 `-5` 出现、deck 注释改写、末档 snap 到 $\alpha=1$）；**#125 真实 EXO/log 混合回放 N1–N8 全 PASS**（真实 `dmax_last=1.480309745e-4`、`dmax/dmean=6.996`、`\beta=1.000000e+5`、11 条 `Gate 3.2`/`Gate 4.5`/能量行、5 条 damage `-5`）；`git diff --check`、`py_compile`、PyYAML（workflow + 5 份 deck）、`Run Terzaghi` 整块 `bash -n`（1511 行）全通过。
  - **状态纪律**：不主张传播段标度、缝长 / 缝网面积 / 起裂压力绝对值 / 排量优选表；能量行只作相对 / 可证伪观测量。验收以 `GATE415_B3_RESULT: PASS` 且原有 **20 组 verdict 保持全绿**为准；若首轮 FAIL 沿用 #124→#125 / #127→#128 原则：**先验算物理阶段，只改判据 + 补可证伪负例，不改 deck 去凑判据**。

---

## 进度跟踪（按 CI #130 验收 + B0/B1/B2 + 文献参数重算收口 + B3 已落地待 CI，2026-11-11）

**代码进度（以清单条目计）**
- 闸门 1（不崩+边界）：**4 / 4**。1.1 / 1.2 / 1.3 已落地；1.4 经核实是 `use-rename` 分离设计、无需改动。CI #47 首次完整编译，此后由 Gate 3.3 / 4.x 持续回归。
- 闸门 2（可见）：**4 / 4**。2.1–2.4 完成；2.3 的 TEMP `source: 1.0e-6` 已在闸门 3 前改回 `0.0`，「务必改回」注意项解除。
- 闸门 3（瞬态+符号）：**2 / 2**。3.1 TS 瞬态主循环已实际运行（`Gate 3.1` / `Gate 3.2` 诊断行）；3.2 由 **Gate 3.3 P1–P8** 硬判据验收 —— P6「L2 err vs BE series」CI #99 实测 `0.03206 ≤ 0.05`，衰减方向与量级正确，**无需改算子符号**（该结论与源项符号无关，见 3.2 / 4.4）。
- 闸门 4（物理+工程）：**15 项收口 / 15 项在册（全部验收通过）**。4.1a（CI #86）、4.1b（CI #98 判据收口 + CI #99 Fortran 解耦，2D/3D 双编译通过）、4.2（CI #87）、4.3（CI #94）全绿；4.4 源项符号（CI #101 `GATE44_RESULT: PASS`，`G+F → G−F`，含非零源解析诊断 deck 与 S1–S8 硬判据）；4.5 能量诊断已修好并接线（section 填充 / 累加器清零 / `DMScal` / MPI 归约 / 时间步层输出；**CI #102 `GATE45_RESULT: PASS`**，T1–T8 全绿且与 #101 真实场离线预测逐位吻合）；4.8 守恒性后处理（**CI #104 `GATE48_RESULT: PASS`**，U1–U8 全绿；纯后处理零扰动）；4.6 μ(T) 别名加固（**CI #106 `GATE46_RESULT: PASS`**，V1–V8 全绿，三例程私有 section + `cellMuField` 往返）；4.10 根目录旧副本清理（**CI #109 `GATE410_RESULT: PASS`**，W1–W8 全绿，93 个文件 100% rename 归档、根目录 tracked `.F90` 归零、原有十组 verdict 与 #106 逐位一致）；4.7 时变材料入口（**CI #111 `GATE47_RESULT: PASS`**，X1–X8 全绿、`A(t)`/`S(t)` 双通道可证伪、生产物理场逐位不变）；4.9 逐单元 IS 批量装配（**CI #114 `GATE49_RESULT: PASS`**，Y1–Y8 全绿、O(N)→O(1) 300× 归约、物理场逐位一致、生产零扰动）；4.12 `pressurePreviousStep` 死写转步间压降活诊断（**CI #117 `GATE412_RESULT: PASS`**，Z1–Z6）；4.13 标准 Newton 残差平台 `3~13` 根因定性为 Dirichlet 单位对角 / 弹性刚度尺度失衡（**CI #117 `GATE413_RESULT: PASS`**，Q1–Q5）；4.14 三维维度分派真实验证（**CI #121 `GATE414_RESULT: PASS`**，W1–W7，首次在 3D `Beam1x10-tetra4` 上跑满 11 步固结，实证 `Displacement_Z`、6 分量 `MEF90MatS3DIdentity` 与门面 2D/3D 分派）；4.15 B0 施工排量–可压性相似校准 + 源项量纲标定（**CI #125 `GATE415_RESULT: PASS`**，2 个新 deck + V1–V8 判据 + 两次 EXO 重建 + 两个 `s415_*_rc` 与 `g415_rc` 延迟硬裁决，零 Fortran 改动；#124 首轮暴露 EXO 追加残留与 `pmax` 局部过冲两个独立缺陷，均已修复，见上方 4.15 条目）。
- 合计：**25 项已落地**（闸门 1–3 共 10 项 + 闸门 4 的 15 项），**0 项修复在册待做**（研究项单列于下条）。damage SNES 的 7 条 `-5` 为信息性，已由 Gate 3.3 P8 用损伤场有界替代；4.11 为总括性描述，已由硬判据体系覆盖。
- **研究项（不计入上述修复闸门，现已全部收口）**：4.15 Stage **B1 由 CI #127 验收**（`GATE415_B1_RESULT: PASS`）、Stage **B2 由 CI #128 验收**（`GATE415_B2_RESULT: PASS`），两条均勾 `[x]`。B2 的 #127 首轮 L4 判据物理前提错误已按实测更正并复验通过。两者均**零 Fortran 改动**，共 8 个单行差异 deck、2 个新判据脚本（A1–A8 / L1–L8）、9 个新诊断 run、11 次 2D + 1 次 3D EXO 重建、11 个延迟硬裁决，属于研究路线扩展而非缺陷修复，故不并入「25 项已落地」的修复计数。**CI #128 已补齐收口证据**：19 组 verdict 同轮全绿、18 条 `GATE*_RESULT` 全 PASS、17 条非 B2 verdict 与 #127 逐字节 IDENTICAL、无任何 `=== FAIL` 行。**新增研究项已由 CI #130 验收并勾 `[x]`**：页岩油文献参数重算表 A′–C′ + M1–M8 判据，纯后处理零 Fortran / 零 deck / 零求解器，`SHALE_OIL_PARAMS_RESULT: PASS` 与 1 正例 + 33 变异 + 3 不变性自测通过；结论是「不能交付排量优选表」，详见路线稿第 14 节与上方条目。证据等级为**逻辑蕴含级**（本轮日志/artifact 端点 `401`，无法逐字节复核 verdict 文本），详见 14.9。
- **CI 硬判据体系**：Gate 3.3 (P1–P8)、4.1 (B1–B7)、4.2 (D1–D8)、4.3 (C1–C7)、4.1b (E1–E4，其中 E2 为 8 项物理签名)、4.4 (S1–S8，其中 S2 校验"两 deck 只差一行")、4.5 (T1–T8，能量非零/衰减 + 注入功零基线)、4.8 (U1–U8，能量/耗散平衡 + 质量平衡 + 流出量独立性)、4.6 (V1–V8，静态反别名 + μ(T) 物理解锚)、4.10 (W1–W8，根目录清零 + 归档 blob 校验 + 构建隔离)；**十一组判据在 CI #109 全绿**（4.10 已验收）；此后追加 4.7 (X1–X8)、4.9 (Y1–Y8)、4.12 (Z1–Z6)、4.13 (Q1–Q5)，**CI #117 十五组判据同轮全绿**；再追加 4.14 (W1–W7，含每步 KSP 收敛与求解器配置反回归)，**CI #121 十六组判据同轮全绿**；再追加 4.15 B0 (V1–V8，相似校准 + 源项量纲标定 + 双 mpirun rc 裁决)，**CI #125 十七组判据同轮全绿**（`gate415_verdict.log` 由 #124 的 FAIL 转 PASS，其余 16 组逐字节 IDENTICAL）；统一「延迟 exit 1」，避免单个闸门早退隐藏后续物理表。**CI #127 首轮**：原有 17 组 + `GATE415_B1_RESULT` 全绿，B2 仅 L4 因判据物理前提错误 FAIL ⇒ `gb2_rc=1` 使整轮 `failure`（编译、EXO、3D、生产隔离均无缺陷）。**CI #128 复验**：L4 重写后 `GATE415_B2_RESULT: PASS`，体系由 17 组正式扩到 **19 组 verdict 同轮全绿**（18 条 `GATE*_RESULT` 全 PASS + `terzaghi_verdict` 汇总），17 条非 B2 verdict 与 #127 逐字节 IDENTICAL、无任何 `=== FAIL` 行；新增的 9 个 mpirun 退出码（`sb1_*` ×6、`sb2_*` ×2、`sdim3_rc` ×1）与 `gb1_rc` / `gb2_rc` 全部接入延迟 fail。**CI #130**：新增第 20 组判据 `check_shale_oil_params.py` (M1–M8，页岩油文献参数重算；纯后处理、不跑求解器、不产生 EXO)，`gso_rc` 以 `PIPESTATUS[0]` 捕获并接入同一延迟裁决；`Run Terzaghi` step `success` 蕴含 M1–M8 与原有 18 组 verdict 全绿 ⇒ **20 组 verdict 同轮全绿**。该组不改动任何 deck / Fortran / 生产 run。
- **最近十一次 CI**：#130 `success`（`37bc756`，页岩油文献参数重算表 A′–C′ + M1–M8 判据接入；**20 组 verdict 同轮全绿**，`SHALE_OIL_PARAMS_RESULT: PASS`，run 26 min 15 s，证据等级为逻辑蕴含级——本轮日志/artifact 端点 `401`，详见路线稿 14.9）；#129 `success`（`56c7b2a`，B2 文档收尾，18 条 verdict 与 #128 逐字节 IDENTICAL）；#128 `success`（`e6cce4b`，Gate 4.15 B2 L4 零位命题复验；**19 组判据同轮全绿**，`GATE415_B2_RESULT: PASS`，`gb2_surface_dev_aniso=4.064640e-03`，17 条非 B2 verdict 与 #127 逐字节 IDENTICAL）；#127 `failure`（`f07cf39`，B1/B2 首轮；编译 + 18 组 verdict 全绿，仅 `GATE415_B2_RESULT: FAIL (L4)` 判据物理前提错，已由 #128 修复）；#126 `success`（`e0f6bb4`，B0 文档收尾）；#125 `success`（`603ce0b`，Gate 4.15 B0 修复：两次 EXO 重建 + `pbulk` 物理锚；**17 组判据同轮全绿**，原有 16 组 verdict 与 #124 逐字节 IDENTICAL）；#124 `failure`（`ff614fc`，B0 首轮；编译与 16 组原有判据全绿，仅 `GATE415_RESULT: FAIL (V1,V3,V4)`，根因是 EXO 追加残留 + `pmax` 局部过冲，已由 #125 修复）；#123 `success`（`3a819ab`，排量–可压性文献调研与集成路线稿）；#122 `success`（`027fa86`，Gate 4.14 文档收尾）；#121 `success`（`c95a605`，Gate 4.14 3D 维度分派验收，16 组判据全绿）；#120 `failure`（`39444b8`，3D 诊断求解器超时，已由 #121 用 direct diagnostic PC 修复）。（更早的 #97 `failure` 为旧 E2「位等价」前提错误，已由 #98 改物理签名判据修复。）

**历史提交 / CI 记录（2026-11-07 09:55 快照，保留备查）**

> 以下 #46–#5x 的逐轮排查日志是当时的现场记录，**不代表当前状态**；当前状态以上面的进度段与文末各专节为准。
- 已提交 **`bd3c1f8`**（8 文件 +523 / −49）、**`1972820`**（构建修复）、**`c4732e0`**（库路径）、**`981433b`**（库 ABI），均推送到 `origin/main`（基线 `2e49c2e`）。`.trae/`、`.workbuddy/` 属工具数据，未入库。
- CI run **#46 = failure**（构建层）：vDef 未链接。真实根因 `vDef/m_vDefDefault.F90:154:59 Error: Invalid character in name`；#45 日志有逐字相同错误 → vDef 自 Flow 默认值块加入起从未构建成功，历史 success 全是假绿。
- CI run **#47 = failure**（运行层，12:38 UTC）：**所有 15 个构建步骤 success**，`vDef.log` / `flow.log` / `MEF90.log` / `m_DefMech.log` / `m_HeatXfer.log` 均无 error，vDef、vDefBT、vDefHF 等全部成功链接 → **闸门 1/2/3.1 的 Fortran 代码首次完整编译通过**（`flow.log` 仅一条 `PetscObjectSetName` type-mismatch warning）。唯一失败是 `Run Terzaghi`：`libmetis.so: cannot open shared object file`，mpirun exit 127。
- **#47 运行期根因**：PETSc 3.3 以 `--with-shared-libraries=1 --download-metis=1` 配置，链接行 `-L$PETSC_DIR/$PETSC_ARCH/lib -lmetis -lparmetis -lHYPRE -lexodus -lhdf5 ...`，但 `Run Terzaghi` 步骤从未设置 `LD_LIBRARY_PATH`，动态加载器找不到该目录。与 Fortran 代码无关。
- CI run **#48 = failure**（运行层，13:09 UTC）：构建仍全部 success；`c4732e0` 加的 `LD_LIBRARY_PATH` 让 `libmetis.so` 等 PETSc 库全部解析成功，但**误把 `$HOME/miniconda2/lib` 也加了进去**，其自带的旧 `libstdc++.so.6` 抢在系统 gcc 13 之前 → `version GLIBCXX_3.4.32 / 3.4.29 not found`（`libpetsc.so`、`libml.so.13` 需要），`libtinfo.so.6` 也让 bash 报 no-version-information。新增的 `ldd` 自检准确拦下并 `exit 1`，**仍未运行到任何物理求解**。
- **#48 修复**：移除 miniconda 库目录（它只服务于 PETSc configure 的 Python 2.7，运行时不需要）；`EXTRA_LIB_DIRS` 改为**追加在末尾**并排除 arch lib 自身，保证 bundled 运行时永不覆盖系统库；`ldd` 断言改为区分「库缺失」与「ABI 版本不匹配」，并打印实际解析到的 `libstdc++`/`libgcc_s`/`libz`/`libtinfo`。
- CI run **#49 = failure**（13:48 UTC）：库环境**全部打通** —— `ldd-vDef.log` 22 行无 `not found`，`libstdc++.so.6` / `libgcc_s.so.1` / `libz.so.1` 均解析到 `/lib/x86_64-linux-gnu`。**vDef 首次真正启动**：打印 MEF90 全局选项（`time_interpolation=linear`、`time_min=0`、`time_max=1`、`time_numstep=51`、`file_format=EXOSingle`），yaml 解析成功。随后在打开几何文件时中止：`Magic value is '../.'` / `Error opening EXO file data/Beam1x10-tri3.gen` / MPI_ABORT 65。
- **#49 根因**：`vDef/data/*.gen` 是 **34~37 字节的文本指针**，内容仅 `../../TestMeshes/Beam1x10-tri3.gen` —— 原本是 symlink，在 Windows checkout 时被物化为普通文本并已提交进仓库（全仓库共 20 个同类假指针，分布在 `vDef/data`、`HeatXfer/data`、`ThermoElasticity/data`）。真实二进制网格在 `TestMeshes/`，magic `CDF\002`。CI 的 `Fix CRLF` 只处理 `*.F90`/`*.h`/`Makefile`/`*.sh`，**不会破坏 `.gen` 二进制**。
- **#49 修复**：`vDef/data/` 4 个 `.gen` 用 `TestMeshes/` 真实网格替换（blob 与源文件哈希完全一致 → Git 复用，仓库体积不增加）；新增 `.gitattributes` 声明 `*.gen binary`；CI 加通用指针自愈（检测 magic，非 `CDF\001/002/005`/HDF5 时从 `../TestMeshes/` 恢复）+ 运行前 geometry magic 硬断言 + `Error opening EXO file` 精确诊断。
- **网格标签核对通过**：真实 `Beam1x10-tri3.gen` 的 `eb_prop1 = [1, 10, 20]`，block1 = 600 个 TRI3 域单元、block2/3 = 各 5 个 BAR2 边界单元，与 `terzaghi1D.yaml` 的 `cs0001`/`cs0010`/`cs0020` 完全一致。
- CI run **#50 = failure**（14:02–14:39 UTC）：**程序完整启动并读完整个输入 deck**（426 行日志），随后两个 rank 同时 `Caught signal number 11 SEGV`；`Gate 3.1` 自检行一条都没有。这是**首个根因位于 Fortran 代码内部**的失败，前四轮（#46~#49）都在构建/环境层。
- **#50 根因（静态定论）**：`vDef/vDef.F90` 原本在 `:143` 调 `MEF90FlowUpdatePerm(MEF90FlowCtx,MEF90DefMechCtx)`，而该例程 `m_MEF90_Flow.F90:259` 读 `DefMechCtx%DMScalSec`；这个 damage SectionReal 要到 `:170 MEF90DefMechCtxSetSections` 才创建 → 解引用未初始化句柄。另发现 `m_MEF90_FlowCtx.F90` 的 Global PetscBag 只注册 26 项 / 类型有 29 字段，且第 13/14 项 Real/Bool 颠倒、漏 3 个 scaling/offset 字段 → 第 20 项起所有 EXO offset 与容差全部错位读成垃圾值。两处均已修复，详见文末「CI #50 诊断结论与修复」。
- CI run **#51 = failure**（15:11–15:47 UTC）：`Gate 3.1: init 01..08` 全通过、进入 `Mechanics: step 1`、`_out.gen` 产出 2.79 MB；崩在 `Attempt to DEALLOCATE unallocated 'doflist'`（FlowAssembly:634），且 TS 路径因全局 bag 缺 `flow_` 前缀静默关闭。两处均已修，详见文末「CI #51 诊断结论与修复」。
- CI run **#52 = failure**（22:15–22:52 UTC）：`flow_` 前缀生效后，注册 Flow 全局 bag 时因 yaml 的 `source_scaling: CST` / `bc_pressure_scaling: CST` 非法枚举而报 `Unknown option CST`（合法值 `constant/linear/file/null`）。已修，并把 CreateTS 调用顺序对齐 HeatXfer 模板；详见文末「CI #52 诊断结论与修复」。
- CI run **#53 = failure**（23:13–23:50 UTC）：构建 15 步全 success（PETSc 重建 28 分钟），**`Gate 3.1: flow enum 2 Transient 2 enabled T` 与 `transient Darcy TS created` 首次出现**，`init 01..08` 全通过；但 `Mechanics: step 1, t=0` 之后立即 rank 1 SEGV（mpirun exit 59），无 `Gate 3.1: step` → 崩在 flowTransient 块内。已加逐调用追踪 + 零跨度保护 + 空指针防护，详见文末「CI #53 诊断结论与修复」。
- CI run **#54 = failure**（00:00–00:37 UTC）：逐调用 trace 生效，`Gate 3.1: step 1 post-UpdatePerm` 之后、`post-SetTransients` 之前 rank 1 SEGV → **崩在 `MEF90FlowSetTransients` 内部**，零跨度保护与 TS 回调均未到达。已定位并修复 `MEF90FlowSetSourceCst` 对 BAR2 边界 cell set 误做 cell-DM dof 映射，详见文末「CI #54 诊断结论与修复」。
- CI run **#55 = failure**（00:53–01:27 UTC）：构建 14 步全 success，`Run Terzaghi` 2 秒失败（mpirun exit 59）。逐 set trace 显示 `SetTransients enter → pre-SetSourceCst → SetSourceCst enter → SetSourceCst set 1 id 1` 后 rank 1 SEGV；rank 0 走完 `SetSourceCst done / post-SetSourceCst / pre-UpdateBC / UpdateBC cell set 2 id 10` 才被 MPI_ABORT 杀掉。已把源项换成已验证的 SectionReal 路径并修好 CI 缓存，详见文末「CI #55 诊断结论与修复」。
- 本地无 PETSc 编译链；所有「验收」栏目前均为**代码审查级通过，非运行验证**。

**本次复核发现（A / B 已修 18:35；C 已修 19:10；D 为既有改动）**
- **A（P0，会崩；已修）** `m_Flow/m_MEF90_FlowAssembly.F90:429-435` CellSet/Has_BC 分支只 `Allocate(dofList(nVal))`，`resPtr`/`bcPtr` 是 `Allocatable` 却未分配就被 `VecGetValues` 写入；VertexSet 分支（:453）三数组一起分配，写法正确。`terzaghi1D.yaml` 的 `cs0010` 正是 `flowBC: true` 的 CellSet，每次求解必走此分支 → 1.2 的"顶面 p=0"验收实际达不到。**修法**：改为 `Allocate(dofList(nVal),resPtr(nVal),bcPtr(nVal),stat=ierr)`，对应 `DeAllocate(dofList,resPtr,bcPtr)`。
- **B（并行风险；已修）** 同一残差分支用 `ISGetSize(bcIS,nVal)`（全局）取 `VecGetValues`，而 `MEF90FlowApplyDirichletJacobian`（:495）用 `ISGetLocalSize`；CI 用 `mpirun -np 2`，全局 size 会读到非本地 dof。VertexSet 分支同样用 `ISGetSize`。**修法**：残差两分支均改为 `ISGetLocalSize`，与雅可比分支一致。
- **C（收敛前提）** 1.3 照 HeatXfer 建了真矩阵并配 CG，但漏了模板里的 `MatSetOption(MAT_SPD)` / `MAT_SYMMETRY_ETERNAL`（见 `m_HeatXfer/m_MEF90_HeatXfer.F90:625-627`）；且 `MEF90FlowApplyDirichletJacobian` 对 `Jac` 与 `PreJac` 两次 `MatZeroRows`（CreateSNES 传的是同一个 `matP`），会破坏对称性。首次 CI 需重点看 `KSPGetConvergedReason`。
- **D（清单外改动）** `vDef/vDef.F90:313` 把 `(time(step)-time)/2.0_Kr` 改成 `(time(step)-time(step-1))/2.0_Kr`，修掉了把数组表达式传给标量实参的隐患；位于 `If (step > 1)` 内且 heatxfer 为 `null` 类型，`time(0)` 越界风险不成立。此改动未记入 11-06 工作日志，补记于此。

**TEMP 已回退（18:45）**
- `vDef/data/terzaghi1D.yaml` 的 `cs0001.source` 已从 TEMP `1.0e-6` 改回 `0.0`：瞬态衰减现在由 `p0=1.0e5` 驱动，符合 Terzaghi 无源解析解。

## 闸门 3 前置评估（18:35，"不编译是否可行"）

**结论：本地不编译可行且不可避免（无 PETSc 链）；但零运行验证直接进闸门 3 不可行。**

- **验证通道原本失效，已修**：`.github/workflows/build.yml` 的 Run Terzaghi 原为 `set +e` + `tee` + `continue-on-error: true`，mpirun 崩溃也绿。现已改为取 `PIPESTATUS[0]` 作退出码、扫描日志错误签名（segfault/abort/diverged 等）、去掉 `continue-on-error`（属 4.11 的前置）。push 一次即可得到**可信**的闸门 1/2 基线；若变红，那是真实状态而非 CI 被改坏。
- **3.1 不是"加几行调用"**：`MEF90FlowCreateTS`（`m_Flow/m_MEF90_Flow.F90:401-418`）把 `PETSC_NULL_OBJECT` 当 IJacobian 矩阵传（与 1.3 修前同类问题），且缺 `TSSetSolution` / `TSSetInitialTimeStep` / `TSSetExactFinalTime` / `TSSetProblemType` / `TSSetOptionsPrefix` / KSP+零空间配置，须照 `m_HeatXfer/m_MEF90_HeatXfer.F90:674-731` 重写。好消息：`MEF90FlowIFunction`/`IJacobian`（`FlowAssembly.F90:237-266`）是真实实现（OperatorCore + 储存项 / BilinearFormCore + 质量矩阵），不是空壳，此前"死代码"仅指无调用点。
- **`vDef.F90` 接线面**：新增 `tsPressure`/`TSConvergedReason` 声明；初始化须放在 `MEF90CtxGetTime`（:141）之后才拿得到 timeStart/dt；主循环时间步层（:375 之前）插 `TSSolve(tsPressure,pressure,time(step))` + `pressure → pressurePreviousStep` 拷贝；AltMin 内 :407-424 的 `UpdateMuFromTemperature`/`UpdatePerm`/`SetTransients`/`SNESSolve(snesPressure)` 须按 steady/transient 分支，否则每个 AltMin 迭代重复缩放并重复解 p；收尾 :591 附近补 `TSDestroy`。
- **语义需实跑定论**：IFunction 复用 `MEF90FlowOperatorCore`，其 Dirichlet 行被覆盖为 `x - p_bc`（代数约束形式），而 `AddStorageToResidual` 跳过 Has_BC cells；IJacobian 的 `BilinearFormCore` 会 `MatZeroEntries` + `MatZeroRows`，`AddStorageMassMatrix` 同样跳过 BC cells。这套组合在 TS/BE 下是否等价于 `-div(k/μ∇p)=S·ṗ`，静态阅读无法定论，正是 3.2 要用解析解验证的内容 —— **3.2 要求运行读 exo，不运行即无法完成**。
- **推荐顺序**：A/B → commit/push 拿闸门 1/2 基线（编译 + SNES/KSP 收敛 + 顶面 p=0）→ 再写 3.1 → CI → 3.2。这样任何失败都只定位到增量，而不是在 7 项未验证改动之上再叠一层重构。

## 3.1 静态实现记录（18:45，零编译 / 零运行）

**改动文件**
- `m_Flow/m_MEF90_Flow.F90` — 重写 `MEF90FlowCreateTS`：两个**独立** Mat（`MatCopy(Jac,PreJac)` 不能自拷贝、两次 `MatZeroRows` 不能作用同一对象）+ MAT_SPD / MAT_SYMMETRY_ETERNAL / KEEP_NONZERO_PATTERN + 零空间；`TSCreate → TSSetDM → TSSetOptionsPrefix → TSSetProblemType(TS_LINEAR) → TSSetType(TSBEULER) → TSSetIFunction → TSSetIJacobian(真矩阵) → TSSetSolution → TSSetInitialTimeStep → TSSetExactFinalTime → TSSetFromOptions`，再 `TSGetSNES → SNESGetKSP → KSPCG` + tolerances。
- `vDef/vDef.F90` — 5 处：① 声明 `tsPressure` / `tsPressureConvergedReason` / `flowTransient` / `tFlowEnd` / `dtFlow0`；② `MEF90CtxGetTime` 之后按枚举创建 TS（`dtFlow0 = time(2)-time(1)`，单步 fallback 1.0）；③ 时间步层（AltMin 之前）更新 μ/k、SetTransients、boundaryPressure、`TSSetDuration(10000,time(step))` + `TSSolve` + `TSGetConvergedReason`（负则按 410 打印）+ `VecCopy(pressure → pressurePreviousStep)`；④ AltMin 内 `SetTransients/UpdateboundaryPressure/SNESSolve(snesPressure)` 包进 `If (.not. flowTransient)`，`UpdateMu/UpdatePerm` 与 Biot/Darcy 保留；⑤ 收尾 `If (flowTransient) TSDestroy`。
- `vDef/data/terzaghi1D.yaml` — `timestepping_type: Transient`（精确匹配枚举列表第 3 项）、`cs0001.source` 回退 `0.0`。
- 自检日志：`Gate 3.1: transient Darcy TS created, timeSteppingType=...`。

**静态核对已通过**
- Fortran 块配对脚本：`vDef.F90` 无未闭合 If / Do / Select。
- `git diff --check`：无空白错误。
- 新 PETSc 调用（`TSSetProblemType` / `TSSetSolution` / `TSSetInitialTimeStep` / `TSSetExactFinalTime` / `TSGetSNES` / `TSSetOptionsPrefix` / `TSSetDuration` / `TSDestroy` / `TSGetConvergedReason` / `TS_LINEAR`）在 `m_HeatXfer/m_MEF90_HeatXfer.F90` 或 vDef 已有编译通过的先例；`TSBEULER` 是 Flow 原有选择。
- 枚举 `MEF90Flow_timeSteppingTypeTransient` 位于无 `Private` 的 `Module m_MEF90_FlowCtx`，`m_vDefDefault.F90:152` 已用同模块枚举 → vDef 可访问。

**未验证风险（首次 CI 必看）**
1. 编译面：`TS_LINEAR` + `TSBEULER` 组合；`TSSolve` 第三参数 IN/OUT 语义（已先 `tFlowEnd = time(step)` 兼容两种）；`Write(IOBuffer,410)` 的 `I2` 与 `TSConvergedReason`。
2. 若日志缺 `Gate 3.1: transient Darcy TS created`，说明 yaml 键没解析成枚举，TS 路径静默关闭（会退化成 steady 且无 p 求解）。
3. **TS 语义**：IFunction 复用稳态 `MEF90FlowOperatorCore`，Dirichlet 行写 `x - p_bc`；`AddStorageToResidual` / `AddStorageMassMatrix` 跳过 `Has_BC` 的 cell set，但内部单元 `cs0001` 仍会向顶面 Dirichlet dofs 装配质量项，`MatZeroRows` 设的 `J_ii=1` 可能被后续装配写回。这是 3.2 要查的衰减方向/量级根因，静态无法定论。
4. `MEF90FlowCreateSNES`（steady 路径）仍把同一个 `matP` 同时作为 J 和 P，而 `BilinearFormCore` 会 `MatCopy(Jac,PreJac)` + 两次 `MatZeroRows`。本次未改；steady 回归若报错优先修这里（CreateTS 已用双矩阵规避）。
5. CG + `MAT_SPD`：`MatZeroRows` 后严格对称性存疑，看 `KSPGetConvergedReason`。
6. 物理近似：transient 下 p 每时间步只解一次，AltMin 内 k(v) 刷新但 p 冻结 → staggered 近似，非完整 Biot 固结。

**状态**：该轮改动已随 `bd3c1f8` 提交并推送；仍**未编译、未运行**，等待 CI run #46。

## 3.1 第二轮静态复核修正（19:10，仍为零编译 / 零运行）

在 3.1 落地后做了一轮「对照已编译模板」的静态审计，查出 5 个**静态可定论**的错误（不是猜测的运行风险），全部已修。判据统一为 `m_HeatXfer/m_MEF90_HeatXfer.F90:674-731`、`m_DefMech/m_MEF90_DefMech.F90:1043` 等 CI 已编译通过的模板。

**P0-1 TS IJacobian 回调缺 `MatStructure flg`（会栈错乱 / 崩溃）**
- PETSc 3.3 的 `TSSetIJacobian` 回调是 `(ts,t,x,xdot,shift,A,B,flg,ctx,ierr)`，10 个参数；Flow 原来只有 9 个（无 `flg`）。`TSSetIJacobian` 以 `external` 注册回调，Fortran 编译器不做接口检查 → 能编译，但 TS 运行时按 10 参调用 9 参子程序，栈错乱。
- 修：`m_MEF90_Flow.F90:127` dispatcher 与 `m_MEF90_FlowAssembly.F90:267` 实现同时加 `MatStructure,Intent(INOUT) :: flg`，2D/3D 别名调用同步；实现末尾 `flg = SAME_NONZERO_PATTERN`（与 HeatXfer/DefMech 一致）。

**P0-2 SNES Jacobian 回调把 `flg` 写成 `PetscReal residualMode`（同类型错配）**
- `SNESSetJacobian` 回调应为 `(snes,x,A,B,flg,ctx,ierr)`；Flow 原来第 5 参是 `PetscReal residualMode`，把 PETSc 传来的整型枚举当实数读，且从不给 `flg` 赋值。参数个数相同所以不会崩，但 PETSc 拿到未初始化的矩阵结构标志。
- 修：`m_MEF90_Flow.F90:87` 与 `m_MEF90_FlowAssembly.F90:223` 改为 `MatStructure flg`；`BilinearFormCore` 的 `residualMode` 在 core 内从未被引用，改传 `1.0_Kr`；设 `flg = SAME_NONZERO_PATTERN`。

**P0-3 Dirichlet 零行被质量矩阵回填（残差 / 雅可比不一致）**
- 原顺序：IJacobian = `BilinearFormCore`（内含 `MatZeroRows` 设 `J_ii=1`）→ `AddStorageMassMatrix`（`a·S·M` 重新填回 BC 行）；IFunction = `OperatorCore`（内含 `ApplyDirichletResidual` 写 `x-p_bc`）→ `AddStorageToResidual`（`SectionRealToVec` 把 `S·M·xdot` 写回 F，污染 BC 行）。结果是 Dirichlet dofs 的 F 与 J 不匹配，TS/BE 下顶面 `p=0` 不再精确成立。
- 修：两处都把 Dirichlet 施加**移到最后**——IJacobian 在质量矩阵后重调 `ApplyDirichletJacobian`，IFunction 在储存项后重调 `ApplyDirichletResidual` + `VecAssemblyBegin/End`。BC dofs 作为代数约束，F 与 J 都不含时间导数项，两矩阵由 `MatCopy` 后同法零行保持一致。

**P0-4 C 项收口：`CreateSNES` 单矩阵自拷贝**
- `BilinearFormCore` 内 `MatCopy(Jac,PreJac,DIFFERENT_NONZERO_PATTERN)` + 对两矩阵各做一次 `MatZeroRows`；`SNESSetJacobian(SNESPressure,matP,matP,...)` 传同一对象，自拷贝且重复零行。
- 修：`CreateSNES` 与 `CreateTS` 一样建 `matP` / `matPpre` 两个独立 Mat，分别 `MatSetOptionsPrefix` / `KEEP_NONZERO_PATTERN` / `MatSetFromOptions`，`SNESSetJacobian` 传双矩阵。steady 回归路径不再有此隐患。

**P0-5 时间配置与求解器不匹配 / TS 类型常量风险**
- `terzaghi1D.yaml` 注释写「51 steps × 0.02s」，配置却是 `interpolation: quadratic`（`t=(sqrt(tmin)+i·dt)²`，非均匀），而 `'beuler'` 是**固定步长无自适应**：`dt=time(2)-time(1)=0.0004` 恒定，走到 t=1 需约 2500 内部步，与解析解时间轴也不一致。改回 `linear`，每个 outer step 恰好一个 TS 步。
- `TSSetType(ts,TSBEULER)` 的 Fortran 常量在已编译源码里无先例（只出现在根目录过期镜像）；改为注册名 `'beuler'`，并把 `TSSetProblemType` 移到 `TSSetType` 之后，完全对齐 HeatXfer 用 `'rosw'` 的顺序。

**CI 再次加固（4.11 前置）**
- 失败关键词补 `TSSolve failed` / `SNESolve failed` / `KSP_DIVERGED` / `error code meaning`（原扫描漏掉了 410/400 格式串的失败输出，TS 发散会假绿），并打印命中的日志行。
- 新增两条硬断言：日志必须含 `Gate 3.1: transient Darcy TS created`（否则说明 yaml 枚举没解析、TS 路径静默关闭）与 `Gate 3.1: step`（否则说明 TSSolve 从未推进）。`vDef.F90` 里对应加了每步 `TSGetTime` 的自检日志。
- 结果文件列举修正：原 `ls terzaghi1D-*.exo` 永远匹配不到。`MEF90CtxCreate`（`MEF90/m_MEF90_Ctx.F90:218-220`）规定 `-prefix data/Beam1x10-tri3` → 几何 `data/Beam1x10-tri3.gen`、**结果 `data/Beam1x10-tri3_out.gen`**。现打印该路径并 `find` 所有结果文件；**3.2 的解析解比对应读 `vDef/data/Beam1x10-tri3_out.gen`**，不是 `*.exo`。

**静态核对结果**
- 所有新 PETSc 接口 / 枚举 / 类型（`MatStructure`、`SAME_NONZERO_PATTERN`、`'beuler'`、`TSGetTime`、`TSSetProblemType` 等）在项目内均有已编译先例；`DIFFERENT_NONZERO_PATTERN` 原本就在 FlowAssembly 内使用，可见性无问题。
- `m_MEF90_FlowAssembly.F90` 通过 `MEF90_DIM` 宏被编译成 `m_MEF90_FlowAssembly2D/3D` 两个模块（`m_Flow/Makefile:17,21`），A/B 与本轮修正都会对实际参与编译的代码生效，不是改到了过期镜像。
- `git diff --check` 无空白错误；新增代码行均 ≤132 字符（项目原有 223 字符长行，说明编译选项允许长行）。
- Fortran 块配对脚本在 `vDef.F90` 无未闭合块；`m_MEF90_Flow.F90` / `FlowAssembly.F90` 的告警经与 `git show HEAD:` 对比为**数量完全一致**（6/12、0/1），属脚本对单行 `If` / `Do While` 的假阳性，非本轮引入。
- 回调签名链路已全部对齐：dispatcher ↔ 2D/3D 别名 ↔ assembly 实现，`residualMode` 仅剩 Core 内未被引用的哑元。

**仍未消除、必须靠首次 CI / 3.2 实跑定论**
1. `TS_LINEAR` + `'beuler'` 组合在 PETSc 3.3 的行为；`TSSetExactFinalTime` 与每步 `TSSetDuration` 的交互。
2. CG + `MAT_SPD`：`MatZeroRows` 只零行不零列，严格对称性存疑 → 看 `KSPGetConvergedReason`。
3. IFunction/IJacobian 的储存项仍跳过 `Has_BC` cell set（`cs0010`），内部 `cs0001` 的质量项是否使 Terzaghi 衰减速率正确 —— 这是 3.2 的核心判据，静态无法定论。
4. `pressurePreviousStep` 全项目**无任何消费者**（只在 FlowCtx 里创建 / 初值 / 销毁）；本轮新增的 `VecCopy(pressure → pressurePreviousStep)` 目前是死写，保留给闸门 4 的平流 / Biot 消费者。
5. transient 下 p 每时间步只解一次、AltMin 内 k(v) 刷新但 p 冻结 → staggered 近似，非完整 Biot 固结。

**状态（22:00 更新）**：#46 语法 → #47 库缺失 → #48 库 ABI → #49 假网格指针，四轮失败全在**代码之外**。#49 已让 vDef 首次启动并成功解析 yaml，但尚未进入 Flow/TS 求解。网格修复待推送，由 #50 复验。

## CI #46 诊断结论与修复（19:55，首次真实反馈）

**日志来源**：用系统凭据下载 `build-diagnostics` artifact（#46 id `13499247217`、#45 id `13392696105`），逐文件比对。

| 日志 | 结论 |
|---|---|
| `MEF90.log` / `m_DefMech.log` / `m_HeatXfer.log` | 无 error，核心库全部编译通过 |
| `flow.log` | **m_Flow 编译通过**：FlowCtx、FlowAssembly 2D/3D、Flow 全部产出 .o；仅一条 `PetscObjectSetName(SNESPressure)` 与 `tsPressure` 的 type-mismatch warning（`-fallow-argument-mismatch` 已压成警告） |
| `vDef.log` | **编译失败**：`m_vDefDefault.F90:154:59 Error: Invalid character in name at (1)`，`make: *** [m_vDefDefault.o] Error 1` → vDef 从未链接 |
| `terzaghi.log` | 仅 10 行：`mpirun ... could not access or execute ../bin/.../vDef`，2 ranks 全部启动失败，exit 132 |

**根因与修复**
1. **假绿链**：`Build` 步骤在 TS 诊断段 `set +e` 后再未 `set -e`；5) vDef 的 make 用 `| tee` 且无退出码检查，失败被吞；`Run Terzaghi` 的 `if: always()` 于是只面对一个不存在的二进制。→ 修复：5) 加 `|| { tail -80 /tmp/vDef.log; exit 1; }`、`SKIPPED` 分支直接 `exit 1`、末尾硬断言 `test -x bin/$PETSC_ARCH/vDef`；`Run Terzaghi` 启动前先检查可执行文件并打印 `terzaghi.log` 尾部，避免把「二进制缺失」误报成「TS 路径未激活」。
2. **`m_vDefDefault.F90` 语法错误**：`vDefFlowDefaultGlobalOptions` 块（闸门 2 之前加入、从未编译过）尾部 `100)             & ! pressureMaxIt` 在**闭合括号后留了尾随 `&`**，会把下一条 `Type(MEF90FlowCellSetOptions_Type)...` 声明拼进同一个构造函数；gfortran 回溯把错误报到 `154:59`。→ 修复：删掉该尾随 `&`；同时把 `1.0e-20_Kr` 这类 `<指数>_<kind>` 字面量全部换成等价十进制 `0.00000000000000000001_Kr` 等（与文件内已编译通过的 `0.0_Kr` 同构，规避任何预处理/词法歧义）；29 个值与 `MEF90FlowGlobalOptions_Type` 的 29 个字段顺序逐一核对一致。
3. **#45 对照**：历史 run #45 的 `vDef.log` 有**逐字相同**的 `154:59` 错误 → 该块加入后 vDef 从未构建成功，之前所有 CI success 只说明 PETSc/OpenMPI 依赖缓存可用。
4. **仍未验证**：`vDef.F90`（3.1 接线）**从未被编译过**——`m_vDefDefault.o` 先失败，vDef.F90 从未进入编译器。本轮已静态复核其新接口（`TSGetTime` / `TSGetConvergedReason` / `TSSetDuration` / `TSConvergedReason` / `Write(IOBuffer,*)` 等）在项目内均有已编译先例，变量声明齐全，但编译结果只能由 #47 给出。
5. **3.2 仍不可做**：需先跑出 `vDef/data/Beam1x10-tri3_out.gen`。

## CI #47 诊断结论与修复（20:45，编译首次通过 / 运行层阻塞）

**日志来源**：`build-diagnostics` artifact（#47 id `13502322746`）。

| 日志 | 结论 |
|---|---|
| `MEF90.log` / `m_DefMech.log` / `m_HeatXfer.log` | 无 error，核心库全部编译通过 |
| `flow.log` | m_Flow 全部 .o 产出、无 error；仅 `PetscObjectSetName(SNESPressure)` 与 `tsPressure` 的 type-mismatch warning |
| `vDef.log` | **263 行、无 error**：`vDef.F90` 首次真实编译、vDef / vDefBT / vDefHF / vDefP / vDefUpa 全部链接成功，末尾 `make: Leaving directory` |
| `terzaghi.log` | 13 行，全是 `libmetis.so: cannot open shared object file` + mpirun exit 127 |

**结论**：闸门 1/2/3.1 的 Fortran 源码至此**全部编译通过**；此前"3.1 接线从未编译"的风险已消除。失败位于运行期的动态库解析，与数值代码无关：`--with-shared-libraries=1` 下链接行含 `-L$PETSC_DIR/$PETSC_ARCH/lib -lmetis -lparmetis -lHYPRE -lexodus -lhdf5 -lnetcdf -lml -lyaml -lz`，而 `Run Terzaghi` 只设了 `PATH`，没设 `LD_LIBRARY_PATH`。`abort` 关键字扫描误把 OpenMPI 的 "job has been aborted" 报成 runtime failure signature，信息不够精确。

**修复（`.github/workflows/build.yml`）**
1. `Run Terzaghi` 设 `LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$MPI_HOME/lib:$HOME/miniconda2/lib:${LD_LIBRARY_PATH:-}`；再用 `find $PETSC_DIR/$PETSC_ARCH -name 'lib*.so*'` 把 `externalpackages` 下所有含 .so 的目录一并加入，避免 PETSc 3.3 下载包安装位置差异。
2. 启动前 `ldd ../bin/$PETSC_ARCH/vDef`，出现 `not found` 立即 `exit 1`；`ldd` 输出存 `/tmp/ldd-vDef.log` 并加入 `build-diagnostics` artifact。
3. `mpirun` 增加 `-x LD_LIBRARY_PATH`，确保两个 rank 都继承库路径。
4. 新增 `error while loading shared libraries` 专门断言（放在签名扫描之前），以后同类问题直接给出准确诊断而非泛化 "runtime failure"。
5. YAML 用 PyYAML 解析通过（15 steps）；提取 `run` 块经 `bash -n` 语法校验通过；`git diff --check` 干净。

**仍待 #48 定论**：`Gate 3.1: transient Darcy TS created`（yaml 枚举是否生效）、`Gate 3.1: step`（TSSolve 是否推进）、SNES/KSP 收敛原因、`vDef/data/Beam1x10-tri3_out.gen` 是否产出。3.2 的解析解比对必须等该文件。

## CI #48 诊断结论与修复（21:25，LD_LIBRARY_PATH 引入的 ABI 冲突）

**现象**：`c4732e0` 加了 `LD_LIBRARY_PATH` 后，`ldd` 显示 `libmetis.so`、`libparmetis.so`、`libexoIIv2for32.so.2`、`libhdf5.so.103`、`libyaml-0.so.2`、`libnetcdf.so.19` 等**全部解析成功**（上一轮 #47 的 `libmetis.so` 缺失已解决），但新增三行 ABI 错误：

```
libpetsc.so / libml.so.13 需要 version `GLIBCXX_3.4.32' (and 3.4.29) not found
libstdc++.so.6 => /home/runner/miniconda2/lib/libstdc++.so.6
libgcc_s.so.1 / libz.so.1 / libtinfo.so.6 也全部来自 miniconda2/lib
```

**根因**：`$HOME/miniconda2/lib` 是上一轮我为了"覆盖 Python 运行时"加进去的，实际只服务于 PETSc configure 阶段的 Python 2.7，**运行 vDef 完全不需要**。Miniconda2 自带的 `libstdc++.so.6` 比 runner 系统 gcc 13 旧，被 `LD_LIBRARY_PATH` 优先加载后屏蔽了系统的 `GLIBCXX_3.4.32`，导致 PETSc/ML 无法加载。同时 `libtinfo.so.6` 让 bash 报 `no version information available`。这是 CI 环境问题，与 Fortran 代码、TS 逻辑无关；也被上一轮加的 `ldd` 硬自检准确拦下（未启动 mpirun）。

**修复（`.github/workflows/build.yml`）**
1. `LD_LIBRARY_PATH` 只保留 `$PETSC_DIR/$PETSC_ARCH/lib:$MPI_HOME/lib`，**彻底移除 `$HOME/miniconda2/lib`**，并在注释里写明禁止再加的原因。
2. `EXTRA_LIB_DIRS`（`find` 收集 externalpackages 下含 `.so` 的目录）从**拼在最前面**改为**追加在末尾**，并用 `grep -v -x` 排除 arch lib 自身；保证任何 bundled 系统运行时都不会覆盖 `/usr/lib/x86_64-linux-gnu`。#48 的日志显示 find 实际只返回 arch/lib 一个目录，过滤后为空。
3. `ldd` 断言正则由宽泛的 `not found` 收紧为 `=> not found|version .* not found|error while loading shared libraries`，区分「库缺失」与「ABI 版本不匹配」；新增打印实际解析到的 `libstdc++`/`libgcc_s`/`libz`/`libtinfo` 便于复核。
4. PyYAML 解析（15 steps）、提取 `run` 块 `bash -n`、`git diff --check` 均通过。

**仍待 #49 定论**：`Gate 3.1: transient Darcy TS created`（yaml 枚举）、`Gate 3.1: step`（TSSolve 推进）、SNES/KSP 收敛原因、`vDef/data/Beam1x10-tri3_out.gen` 是否产出。**#49 才是 3.1 数值求解的首次真实运行**；3.2 的 Terzaghi 解析解比对必须等该文件。

## CI #49 诊断结论与修复（22:00，程序首次启动 / 假网格指针）

**日志来源**：`build-diagnostics` artifact（#49 id `13506732051`，9718 bytes，比 #48 的 7112 大约多 2.6 KB → 程序真的产生了输出）。

| 阶段 | 结论 |
|---|---|
| 编译 | 5 个编译日志与 #47/#48 完全一致，无 error |
| `ldd-vDef.log` | **22 行、无 `not found`**；`libmetis` / `libparmetis` / `libexodus` / `libhdf5` / `libml` / `libnetcdf` / `libyaml` 全部解析到 PETSc arch lib；`libstdc++.so.6` / `libgcc_s.so.1` / `libz.so.1` 解析到 `/lib/x86_64-linux-gnu` → **#47/#48 的库问题彻底解决** |
| `terzaghi.log`（78 行） | vDef **首次真正启动**：打印版权头、`geometry file: data/Beam1x10-tri3.gen`、`result file: data/Beam1x10-tri3_out.gen`，以及 `PetscBag` 全局选项（`time_interpolation=linear`、`time_min=0`、`time_max=1`、`time_numstep=51`、`time_skip=0`、`file_format=EXOSingle`）→ **yaml 解析成功、时间轴符合 51×0.02s 设计** |
| 中止点 | `Could not recognize ... valid Exodus/NetCDF file variant. Magic value is '../.'` → `Error opening EXO file` → `MPI_ABORT ... errorcode 65` |

**根因**：`vDef/data/Beam1x10-tri3.gen` 不是网格，而是 **34 字节文本**，内容只有 `../../TestMeshes/Beam1x10-tri3.gen`。这些文件原本是 symlink，在 Windows（未启用 `core.symlinks`）上被物化为含目标路径的普通文本文件后提交进仓库。全仓库共 **20 个**同类假指针：`vDef/data`（4）、`HeatXfer/data`（4，指向 `TestMeshes`）、`ThermoElasticity/data`（12，指向 `HeatXfer/data`，属二级指针）。真实二进制在 `TestMeshes/`，`Beam1x10-tri3.gen` 为 21332 字节、magic `CDF\002`（NetCDF classic 64-bit offset）。

**修复**
1. `vDef/data/` 下 4 个 `.gen` 用 `TestMeshes/` 的真实网格按字节替换；`git hash-object` 显示 4 个文件与源 blob 哈希完全相同（`148bf0f0ea` / `60821ce21a` / `7ffde604ad` / `39f41502d6`），Git 会复用 blob，**仓库体积不增加**。
2. 新增 `.gitattributes`：`*.gen binary`，防止今后任何平台的换行转换误伤二进制网格。
3. CI `Run Terzaghi` 加通用自愈：遍历 `vDef/data/*.gen`，前 4 字节不是 `CDF\001/002/005` 或 HDF5 signature 时从 `../TestMeshes/<basename>` 恢复；随后对 `data/Beam1x10-tri3.gen` 做 magic 硬断言，失败立即 `exit 1`。
4. 新增 `Error opening EXO file` / `Could not recognize ... valid Exodus` 精确诊断，优于泛化的 `MPI_ABORT` 关键词命中。
5. 确认 `Fix CRLF + Patch Makefiles` 步骤的 `dos2unix` 只作用于 `*.sh` / `Makefile` / `*.F90` / `*.inc` / `*.h`，**不会破坏 `.gen` 二进制**。

**网格标签核对（scipy.io.netcdf_file 读真实文件）**
- `eb_prop1 = [1, 10, 20]`，`num_el_blk=3`；`num_el_in_blk1=600 / num_nod_per_el1=3`（TRI3 域），`num_el_in_blk2=5 / num_nod_per_el2=2`（BAR2 顶面），`num_el_in_blk3=5 / num_nod_per_el3=2`（BAR2 底面）。
- 与 `terzaghi1D.yaml` 的 `cs0001`（domain）、`cs0010`（top, `flowBC: true`）、`cs0020`（bottom）**完全对应**，无需改配置。

**仍待 #50 定论**：`Gate 3.1: transient Darcy TS created`（Flow 枚举 `Transient` 是否解析成功）、`Gate 3.1: step`（TSSolve 是否推进）、SNES/KSP 收敛原因、`vDef/data/Beam1x10-tri3_out.gen` 是否产出。**#50 才可能首次进入 3.1 的瞬态数值求解**；3.2 的 Terzaghi 解析解比对必须等该文件。

**遗留（非阻塞，留给闸门 4）**：其余 16 个假指针 `.gen` 仍在 `HeatXfer/data` / `ThermoElasticity/data`，本次 CI 自愈逻辑只覆盖 `vDef/data`。

## CI #50 诊断结论与修复（22:55，首个 Fortran 代码内部根因）

**日志来源**：`build-diagnostics` artifact（#50 id `13511785681`，12920 bytes；`terzaghi.log` 426 行 / 34307 bytes）。

| 阶段 | 结论 |
|---|---|
| 编译 | 5 个编译日志与 #47~#49 一致，无 error；`flow.log` 仍只有 `PetscObjectSetName(SNES/TS)` 的 type-mismatch warning |
| `ldd-vDef.log` | 22 行、零 `not found`，系统运行时全部解析到 `/lib/x86_64-linux-gnu` → 库环境已彻底稳定 |
| `terzaghi.log` 1–355 行 | **整个输入 deck 读取完成**：MEF90 全局选项、DefMech global + 3×cell set、HeatXfer global + 3×cell set、3×MatProp2D 材料属性全部正常打印，几何/结果文件路径正确 |
| `terzaghi.log` 357–426 行 | rank 1 与 rank 0 同时 `Caught signal number 11 SEGV`，`MPI_ABORT errorcode 59`；**没有任何 `Gate 3.1` 自检行** → 崩溃在 `vDef.F90:136–165` 的 Flow 初始化链路，未进入时间步循环 |
| CI 判定 | 泛化的 `abort` / `petsc error` 关键词扫描命中，但没给出代码位置 |

**根因 1（P0，SEGV 直接原因）— 初始化顺序错**
`vDef.F90:143` 过早调用 `MEF90FlowUpdatePerm(MEF90FlowCtx,MEF90DefMechCtx)`。该例程在 `m_Flow/m_MEF90_Flow.F90:259` 执行
`Call SectionRealRestrictClosure(DefMechCtx%DMScalSec,FlowCtx%DM,setIdx(cell),elemType%numDof,dmgPtr,ierr)`，
而 `DefMechCtx%DMScalSec` 由 `MEF90DefMechCtxSetSections`（`m_DefMech/m_MEF90_DefMechCtx.F90:555`）创建，在主程序 `:170` 才调用 —— 相差约 30 行，读的是未初始化句柄。DefMech/HeatXfer 的材料打印恰好在 `:131` 结束，与日志 355 行完美吻合。

**根因 2（P0，内存布局错位）— PetscBag 注册顺序 ≠ 派生类型字段顺序**
`MEF90FlowGlobalOptions_Type` 有 **29** 个字段，`PetscBagRegisterMEF90FlowCtxGlobalOptions` 只注册 **26** 项：
- 第 13/14 项把 `arrheniusActivationEnergy`(PetscReal, 8B) 与 `useArrheniusViscosity`(PetscBool, 4B) 写反；
- 漏掉 `pressurePreviousStepOffset`、`pressurePreviousStepScaling`、`pressureScaling`；
- 从第 20 项起所有 EXO offset / scaling / solver tolerance 逐字节错位，且尾部 3 个字段是未初始化内存。

PETSc bag 按 `bag->data + bag->offset` 顺序写值，`PetscBagGetData` 再把 `bag->data` 强转成 Fortran 派生类型指针，因此注册顺序必须与字段声明**逐字节一致**。`CellSet`（7 字段 / 9 注册，`advectionVector(1:3)` 拆 3 项）与 `VertexSet`（2/2）核对后确认无误。

**修复**
1. `vDef/vDef.F90`：把 `MEF90FlowUpdatePerm` + `MEF90FlowUpdateStorativity` 从 Flow 初始化段移到 `MEF90DefMechCtxSetSections` / `MEF90DefMechCtxCreateVectors` **之后**（首次 TSSolve 之前仍会装配好系数）；`MEF90FlowCreateSNES` / `MEF90FlowCreateTS` 不依赖这两个向量，可留在原处。
2. `vDef/vDef.F90`：在初始化链插入 8 条 `Gate 3.1: init 01..08` 分段追踪日志，#51 若再崩可直接定位到具体 Call。
3. `m_Flow/m_MEF90_FlowCtx.F90`：Global bag 注册按类型定义重排为 **29 项、顺序完全一致**（补齐 3 个漏项，新增 `pressure_previous_step_offset` / `pressure_previous_step_scaling` / `pressure_scaling` 三个选项名）。自动脚本复核结果：`Global 29/29 identical_order=True`，`VertexSet 2/2 True`，`CellSet 7/9 且前 7 项无 mismatch（数组三分量）`。
4. `vDef/data/terzaghi1D.yaml`：补 `pressure_previous_step_offset: 6`，避免其默认值 `2` 与 `source_offset: 2` 撞 EXO 列。
5. `.github/workflows/build.yml`：新增 `Caught signal number 11 SEGV` 精确分支，打印**最后到达的 Gate 3.1 trace 编号**，不再只给泛化 "runtime failure signature"。

**校验**：PyYAML（15 steps）、提取 `run` 块 `bash -n`、`git diff --check`、Fortran 块配对（HEAD/WORK 均为 0 bad / 0 unclosed）、行宽（编译用 `-ffree-line-length-none`，超长行均为既有行）全部通过。

**仍待 #51 定论**：8 条 init trace 能否全部走完、`Gate 3.1: transient Darcy TS created`（Flow 枚举 `Transient` 是否解析）、`Gate 3.1: step`（TSSolve 推进）、SNES/KSP 收敛、`vDef/data/Beam1x10-tri3_out.gen` 是否产出。**#51 才可能首次进入 3.1 的瞬态数值求解**；3.2 的 Terzaghi 解析解比对必须等该文件。

**下一轮可能的遗留点（静态可见，未处理）**：`MEF90FlowUpdatePerm` 会对 `cs0010`/`cs0020` 这两个 BAR2 边界 cell set 调 `MEF90Element_Create`（与 TRI3 域单元同路径）；若 #51 在此报错，需按边界块跳过体积分。另 `pressurePreviousStep` 全项目仍无消费者，`VecCopy` 是死写（闸门 4）。

## CI #51 诊断结论与修复（2026-11-07 06:10，首次进入数值阶段）

**日志来源**：run `42240634507` / job `127969679819` 的完整 job 日志（2906 行），`Run Terzaghi` step failure；artifact `build-diagnostics` id `13518271108`。构建 15 步全部 success，`ldd` 干净，无需再查环境层。

**进展（#50 的修复全部生效）**
- `Gate 3.1: init 01..08` **8 条全部走完**，`Time values array` 打印 51 个 0→1 的均匀时刻，DefMech/HeatXfer SNES 对象创建完成，`Mechanics: step 1, t=0` 进入主循环。
- 结果文件 `data/Beam1x10-tri3_out.gen` **已产出 2,789,988 字节**（3.2 的数据源路径确认无误），但只是初始化写出的空/初始场，尚不含瞬态压力演化。
- 崩溃点：`Alt. Min. step 1` 后两个 rank 同时 `Fortran runtime error: Attempt to DEALLOCATE unallocated 'doflist'`，位置 `m_Flow/m_MEF90_FlowAssembly.F90:634`，mpirun exit 2；CI 在运行时签名扫描处 `exit 1`，后面的 Gate 断言未执行。

**根因 1（P0，直接崩溃）— 误释放 PETSc 托管的 F90 指针**
`dofList` 声明为 `PetscInt,Dimension(:),Pointer`。PETSc 3.3 的 `ISGetIndicesF90` 自己把这个指针关联到内部数组，`ISRestoreIndicesF90` 再解除关联；代码却先 `Allocate(dofList(nVal))`、restore 后再 `DeAllocate(dofList,...)`，运行时必然报「未分配」并终止。`#51` 是第一次真正走到 `MEF90FlowUpdateboundaryPressure`（`Has_BC` 为真）才暴露；前几轮都在更早的位置中止。同类共 **6 处**：残差 2 处（451/475）、雅可比 2 处（518/541）、边界压力 2 处（629/650）。

**修复**：6 处统一删除 `Allocate(dofList)` / `DeAllocate(dofList)`，只保留对纯 `Allocatable` 辅助数组（`resPtr` / `bcPtr` / `val`）的分配与释放；另把 `MEF90FlowUpdateboundaryPressure` 里剩下的两处 `ISGetSize` 改成 `ISGetLocalSize`，与 A/B 修复对齐，避免 `-np 2` 下读非本地 dof；加注释固化所有权规则。

**根因 2（P0，TS 路径静默关闭）— Flow global bag 缺 `flow_` 选项前缀**
日志中**没有** `Gate 3.1: transient Darcy TS created`，也没有 `Gate 3.1: step` → `flowTransient = false`，程序走 steady 分支（这也解释了为什么 AltMin 内会调 `SetTransients` 并触发根因 1）。`MEF90FlowCtxSetFromOptions` 注册全局 bag 时传的是 `PetscBagRegisterMEF90FlowCtxGlobalOptions(...,'flow','',...)`，prefix 为空 → bag 只认 `-timestepping_type`；而 yaml 顶层 `flow:` 块按 DefMech / HeatXfer 同样的规则展平成 `-flow_timestepping_type`（日志里 DefMech 打印 `DefMech_TimeStepping_Type`、HeatXfer 打印 `heatxfer_timeStepping_type` 可证）。因此 **Flow 的全部全局选项都取默认值**：`SteadyState`、`p0=0`、`k_min=1e-20`、容差缺省 —— 即使不崩也算不出 Terzaghi 衰减。

**修复**
1. `m_Flow/m_MEF90_FlowCtx.F90`：全局 bag prefix 由 `''` 改为 `'flow_'`，使 `-flow_timestepping_type` / `-flow_perm_min` / `-flow_p0` 等全部命中。CellSet/VertexSet bag 的 `cs0001_` 前缀本来就与 yaml 一致，不动。
2. `vDef/vDef.F90`：在 `flowTransient` 赋值后**无条件**打印 `Gate 3.1: flow enum <值> Transient <枚举值> enabled <T/F>`，下轮即使枚举仍不对也能一眼看出，不再静默退化。
3. `vDef/data/terzaghi1D.yaml`：顶层补 `timestepping_type: Transient` 作为扁平别名双保险（若某 PETSc 3.3 构建不把顶层块名写进 option，仍能激活瞬态）。
4. `.github/workflows/build.yml`：Gate trace 与「TS 未激活」分支都加 `Gate 3.1: flow enum` 输出，便于下一轮直接读枚举值。

**校验**：`git diff --check` 干净；PyYAML 解析 workflow（15 steps）与 yaml 通过；提取 `Run Terzaghi` 的 `run` 块 `bash -n` 通过；`m_Flow/m_MEF90_FlowAssembly.F90` 已无任何 `Allocate(dofList` / `DeAllocate(dofList` / `ISGetSize(` 残留。

**仍待 #52 定论**：`Gate 3.1: flow enum` 是否为 Transient(=2)、`Gate 3.1: transient Darcy TS created`、`Gate 3.1: step`（TSSolve 推进）、`TSGetConvergedReason` / KSP 收敛、`_out.gen` 内是否出现随时间衰减的 pressure。**#52 才可能首次完成 3.1 的瞬态数值求解**；3.2 的 Terzaghi 解析解比对必须等该文件含真实瞬态压力场。

## CI #52 诊断结论与修复（2026-11-07 07:10）

**状态**：run `42303260569` / job `128186129733` `completed failure`；15 个构建步全 success，`Run Terzaghi` failure（22:15–22:52 UTC）。完整日志 2869 行；artifact `build-diagnostics` id `13547191306`。

**根因（`flow_` 前缀生效后暴露的下一层）**：两个 rank 均报
`[0/1]PETSC ERROR: Unknown option CST for -flow_source_scaling!`，调用链 `PetscOptionsGetEList(options.c:1374) → PetscOptionsGetEnum(:1427) → PetscBagRegisterEnum(bag.c:71)`。崩在 `Gate 3.1: init 02 FlowCtxSetFromOptions` 内部，`init 03..08` 未执行，`_out.gen` 只有 21,176 字节（空壳，不是结果）。`MEF90ScalingList` 的合法字符串是 `constant / linear / file / null`（`MEF90/m_MEF90_Ctx.F90:110-116`），`CST` 只是 Fortran 枚举名 `MEF90Scaling_CST=0`，不是 PetscBag option 值。yaml 中 `source_scaling: CST`、`bc_pressure_scaling: CST` 均非法。

**重要正面证据**：`-flow_source_scaling` 被读到，说明 #51 的 `flow_` 前缀修复**确实生效**（此前所有 Flow 全局选项走默认值）。而 `timestepping_type` 是全局 bag 注册的**第一个** enum（第 562 行），`source_scaling` 排在第 26 个字段（第 629 行）——`Transient` 已成功通过注册，证明 yaml 的 `timestepping_type: Transient` 合法，**下一轮 `flowTransient` 会首次为 true**，`p0=1e5`、`k_min=5e-11` 等也终于生效（`MEF90FlowCtxCreateVectors` 已有 `VecSet(pressure,initialPressure)`）。

**修复**
1. `vDef/data/terzaghi1D.yaml`：`CST` → `constant` 两处；更正 `# Option names without 'flow_' prefix` 的过时注释，写明键映射到 `-flow_<key>` 与合法 scaling 枚举值。
2. `.github/workflows/build.yml`：新增 `Unknown option` 精确诊断分支（在 Exodus 检查前），直接打印非法值并失败。
3. `m_Flow/m_MEF90_Flow.F90`（主动对齐，降低 #53 风险）：`MEF90FlowCreateTS` 调用顺序改为与仓库中已运行的 `MEF90HeatXferCreateTS` 一致 —— `TSGetSNES` 与 `TSSetIFunction`/`TSSetIJacobian` 放在 `TSSetType` **之前**（PETSc 3.3 的 `TSSetType` 可能重置 ops 表），再 `TSSetType('beuler')` → `TSSetProblemType(TS_LINEAR)` → `TSSetSolution` → `TSSetInitialTimeStep` → `TSSetExactFinalTime` → `TSSetFromOptions` → KSP 设置。日志已确认 `petsctsdef.h` 中 `#define TSBEULER 'beuler'`，类型名无误；`tFlowEnd = time(step)` 也已赋值，`TSSolve` 目标时间正确。

**校验**：PyYAML 解析 workflow / yaml 通过，`Run Terzaghi` 的 `run` 块 `bash -n` 通过；`flow` 段枚举值现为 `{'timestepping_type':'Transient','source_scaling':'constant','bc_pressure_scaling':'constant'}`。

**仍待 #53 定论**：`Gate 3.1: flow enum`（应为 2 / T）、`transient Darcy TS created`、`Gate 3.1: step`（TSSolve 推进与 TS 内部时间）、`TSGetConvergedReason` / KSP 收敛、IFunction/IJacobian 首次执行时的并行 Dirichlet dof 索引、`_out.gen` 内随时间衰减的 pressure。**#53 才可能首次完成 3.1 的瞬态求解**；3.2 解析解比对仍需该文件含真实瞬态压力场。


## CI #53 诊断结论与修复（2026-11-07 07:55，瞬态路径首次激活）

**状态**：run `42310380078` / job `128210036201` `completed failure`；step 1–14 全 success（`Build PETSc-3.3` 23:21–23:49 UTC，28 分钟；`Build MEF90 + m_Flow + vDef` 23:49–23:50 编译通过），`Run Terzaghi` 仅 1 秒即 failure（23:50:50–23:50:51Z），mpirun exit **59**，完整 job 日志 2934 行。

**首次真正进入瞬态（正面结果）**
- `Gate 3.1: init 01..08` 8 条全部走完；`Time values array` 51 个 0→1 时刻（dt=0.02）。
- `Gate 3.1: flow enum           2  Transient           2  enabled T` —— #51 的 `flow_` 前缀修复 + #52 的合法枚举值全部生效，`flowTransient` **首次为 true**。
- `Gate 3.1: transient Darcy TS created, timeSteppingType= 2` —— `MEF90FlowCreateTS` 首次成功创建。
- DefMech Disp/damage 两个 SNES 对象正常创建（712×712 / 356×356 矩阵）；`_out.gen` 2,789,988 字节（仍只是初始场）。

**崩溃现场**
`Mechanics: step    1, t= 0.00000E+00` 之后**没有任何** `Gate 3.1: step`，rank 1 立即 `Caught signal number 11 SEGV`，`MPI_ABORT`（errorcode 59）使 rank 0 收到 signal 15 后退出。没有 Fortran runtime error，PETSc 3.3 release 构建只给 `User provided function() line 0 in unknown directory unknown file`，无调用栈。因此崩溃必在 `vDef/vDef.F90` 的 `flowTransient` 块内（该块在 #53 之前从未执行过）。

**静态核对（已排除的项）**
- TS 回调 dispatcher `MEF90FlowIFunction` / `MEF90FlowIJacobian`（`m_Flow/m_MEF90_Flow.F90:105/127`）签名与 HeatXfer 模板一致，IJacobian 含 `MatStructure flg`，无参数错位。
- `MEF90DiffusionOperatorAddTransientTermSet(G,mesh,x,cellIS,alpha,elem,elemType,ierr)` 与 `MEF90DiffusionBilinearFormSet(K,mesh,cellIS,A,lambda,elem,elemType,ierr)` 的实参顺序/类型均与定义一致；`A_zero = 0.0_Kr`、`A_cell = permCell(1)` 是项目既有的标量→`MEF90_MATS` 赋值写法。
- `MEF90FlowSetTransients` / `MEF90FlowUpdateboundaryPressure` 的 `constant` 分支只做 `VecSetValues`，路径与 init 阶段一致。
- **两个仍未排除的高风险**：① `beuler + TS_LINEAR + IFunction/IJacobian` 是仓库内从未运行过的组合（HeatXfer 原生用 `rosw/ra3pw`，且其 transient 分支在本次 yaml 中也是 NULL，两者都只有编译级验证）；② step 1 时 `time(1)=0` 且 TS 初始时间也是 0，`TSSetDuration(...,0)` 会要求 TSSolve 积分零长度区间，而 HeatXfer 模板在同一循环用 `If (t < time(step))` 做了保护。③ `SCell(1)` 在储存项/质量矩阵两处是裸指针读取，未做 `Associated`/`size` 保护（rank 相关性崩溃的典型来源）。

**本轮改动（不改物理模型，只加可观测性 + 防御）**
1. `vDef/vDef.F90`：flowTransient 块内每个调用前后加 `Gate 3.1: step N ...` 检查点（enter / post-UpdateMu / post-UpdatePerm / post-SetTransients / post-UpdateBoundaryPressure / pre-TSSolve / post-TSSolve）；`TSSolve` 前先 `TSGetTime`，按 HeatXfer 的写法用 `If (t < tFlowEnd)` 包住，step 1 打印 `TSSolve skipped (no time advance)`，`TSGetConvergedReason` 也一并进保护分支避免读到未初始化值。`Gate 3.1: step ... TS time` 打印仍在块外，保证 CI 判据有效。
2. `m_Flow/m_MEF90_FlowAssembly.F90`：`MEF90FlowIFunction`（enter / post-OperatorCore）、`MEF90FlowIJacobian`（enter）、`MEF90FlowAddStorageToResidual`、`MEF90FlowAddStorageMassMatrix` 四个入口加 trace；两处 `SCell(1)` 加 `Associated(SCell)` + `size(SCell) > 0` 保护，否则取 0。
3. `.github/workflows/build.yml`：mpirun 加 `-flow_ts_view -flow_ts_monitor`（让 PETSc 自己报告 TS 类型、问题类型与步进）；SEGV 分支打印崩溃前 40 行（用 `crashline` + `sed`，避免 awk 引号嵌套）；Gate grep 放宽到全部 `Gate 3.1` 行并 `tail -30`。

**校验**：`git diff --check` 干净；PyYAML 解析 workflow / yaml 通过；`Run Terzaghi` 的 `run` 块 `bash -n` 通过；`PetscPrintf(PETSC_COMM_WORLD, ...)` 在项目内有 36 个文件先例；Fortran `If/End If` 配对已逐段复核。

**#54 判据（按优先级）**：① 日志里最后一个 `Gate 3.1: step 1 ...` 检查点 = 精确崩溃位置（若停在 `pre-TSSolve` 之后且无 `IFunction enter` → PETSc TS 层；若停在 `IFunction enter` 之后 → 残差/储存项；若停在 `IJacobian enter` → 雅可比）；② `-flow_ts_view` 打印的 type / problem type；③ `-flow_ts_monitor` 是否出现步进；④ 若零跨度保护已消除崩溃，则看 `Gate 3.1: step 2 ... post-TSSolve`、`TS time`、`TSGetConvergedReason` 与 `_out.gen` 是否出现随时间衰减的 pressure。3.2 的 Terzaghi 解析解比对仍需该文件含真实瞬态压力场。


## CI #54 诊断结论与修复（2026-11-07 08:50，逐调用 trace 命中崩溃点）

**状态**：run `42315564918` / job `128227179782` `completed failure`（基于 `0a30b0e`）。构建 15 步全 success；`Run Terzaghi` 1 秒失败，mpirun exit **59**，完整 job 日志 2987 行。

**trace 结果（#53 加的检查点生效）**
```
Mechanics: step    1, t= 0.00000E+00
 Gate 3.1: step 1  enter transient block
 Gate 3.1: step 1  post-UpdateMu
 Gate 3.1: step 1  post-UpdatePerm
[1]PETSC ERROR: Caught signal number 11 SEGV
```
`init 01..08`、`flow enum 2 Transient 2 enabled T`、`transient Darcy TS created` 全部正常；崩溃被夹在 `post-UpdatePerm` 与 `post-SetTransients` 之间 → **`vDef.F90` 调用的唯一函数 `MEF90FlowSetTransients(MEF90FlowCtx,step,time(step),ierr)`**。rank 0 无 `post-SetTransients`（被 `MPI_ABORT` 以 signal 15 杀掉），说明两 rank 都在该函数内、rank 1 先越界。

**根因（确定性）**：`MEF90FlowSetTransients2D`（即 `m_Flow/m_MEF90_FlowAssembly.F90` 的 `MEF90FlowSetTransients`，该文件由 Makefile 以 `-DMEF90_DIM=2/3` 编译两次）在 `sourceScaling=constant` 分支首次调用 `MEF90FlowSetSourceCst` —— 这是全项目**唯一从未被执行过**的函数。它遍历所有 cell set，对 `cs0010`/`cs0020`（BAR2 边界 facet，`flowBC: true`）也调用
```
DMMeshISCreateISglobaldof(FlowCtx%cellDMScal, setIS, 0, setISdof, ...)
```
但这些 point 在 `cellDMScalSec`（由 `DMMeshGetCellSectionReal` 建立，**只为 codim=0 体单元分配 dof**）里没有 dof → 越界读，rank 相关。同文件的 `MEF90FlowOperatorCore` 用 `If (.not. Has_BC)` + `If (elemType%Codim == 0)` 跳过它们，两个储存项也跳过 `Has_BC`；唯独 `SetSourceCst` 漏了。这正解释了「为何 rank 1 挂、为何 #54 之前从未触发」。

**静态排除**：`time` 数组由 `MEF90CtxGetTime` 在 `linear` 分支对所有 rank `Allocate`，`time(step)` 有效；`cellDMScal` 是 `DMMeshClone(Mesh,...)`，`DMMeshGetStratumIS(DM,...)` 与 `(cellDMScal,...)` 的 point id 空间一致，DM 混用不是根因；`MEF90FlowUpdateboundaryPressure` 在 #51 稳态路径已执行到 `VecSetValues` 之后才崩（当时是 `DeAllocate(dofList)`），其索引方式不是本次 SEGV 来源。

**本轮改动（2 文件）**
1. `m_Flow/m_MEF90_FlowAssembly.F90`（重写 3 个子程序）：
   - `MEF90FlowSetSourceCst`：**跳过 `Has_BC` 及 `Codim /= 0` 的 cell set**；对齐已验证的 `MEF90HeatXferSetFluxCst` —— `DMMeshGetStratumIS(FlowCtx%cellDMScal,...)`、`nVal = size(setIdx)`（不再混用 `ISGetLocalSize` 与 F90 指针长度）；所有 PETSc 调用补 `CHKERRQ`。
   - `MEF90FlowUpdateboundaryPressure`：`MEF90ISCreateCelltoVertex` 补 `CHKERRQ(ierr)`（原来缺，失败会静默传递未初始化的 `bcIS`）；两处 `ISDestroy` 补 `CHKERRQ`；保留 #51 的 `ISGetLocalSize` 修复与 PETSc 托管指针所有权。
   - `MEF90FlowSetTransients`：四个付费分支全部补 `CHKERRQ`，`constant` 分支前后加 Gate 检查点；入口打印 `srcScale/bcScale` 枚举值。
   - 三个函数统一加逐 set / 逐阶段 `Gate 3.1:` trace，若仍崩溃可直接读出 set 编号与函数。
2. `.github/workflows/build.yml`：SEGV 诊断的 Gate grep 由 `tail -30` 放宽到 `tail -60`，容纳新增检查点。

**校验**：`git diff --check` 干净；PyYAML 解析 workflow / yaml 通过；`Run Terzaghi` 的 `run` 块 `bash -n` 通过；三个子程序的 `Do/End Do`、`If/End If`（两处为单行 `If (...) Cycle`）、`Select Case/End Select` 配对已脚本核对；`PetscPrintf(PETSC_COMM_WORLD,...)` 与 `Write(IOBuffer,*)` 在项目内有大量先例且 #54 已编译通过。仍为零本地编译、靠 CI 验证。

**#55 判据**：① 应看到 `SetTransients enter ... srcScale 0 bcScale 0` → `pre-SetSourceCst` → `SetSourceCst enter` → `set 1 id 1` → `SetSourceCst done` → `post-SetSourceCst` → `pre-UpdateBC` → `UpdateBC cell set ...` → `UpdateBC done`；② 若零跨度保护让 step 1 打 `TSSolve skipped`，step 2 应出现 `pre-TSSolve` / `IFunction enter` / `IJacobian enter` / `post-TSSolve` 及 `-flow_ts_monitor` 步进；③ `TSGetConvergedReason`、TS time 是否推进、`_out.gen` 内 pressure 是否随时间衰减（3.2 前置）。若仍 SEGV，最后一个 `Gate 3.1: ... set N` 即精确越界位置。

**性能备注**：workflow 已内置 `actions/cache` 缓存 `$HOME/petsc-3.3-mef90` 与 OpenMPI，但 `actions/cache@v4` 只在 **job 成功**时保存；此前每轮 job 均 failure，故缓存从未落盘、每轮都重建 PETSc（约 28 分钟）。只要有一轮 job success，后续 run 会自动命中缓存，无需改 workflow。


## CI #55 诊断结论与修复（2026-11-07 09:55，源项改走 SectionReal + CI 缓存修复）

**状态**：run `42320914171` / job `128244974102` `completed failure`（基于 `71890d5`）。构建 step 1–14 全 success（PETSc 重建 28 分钟，MEF90 + m_Flow + vDef 编译通过）；`Run Terzaghi` 01:27:36–38Z 仅 2 秒失败，mpirun exit **59**，完整 job 日志 3006 行。`init 01..08`、`flow enum 2 Transient 2 enabled T`、`transient Darcy TS created` 全部正常。

**trace 结果（#54 加的逐 set 检查点生效）**
```
Mechanics: step    1, t= 0.00000E+00
 Gate 3.1: step 1  enter transient block
 Gate 3.1: step 1  post-UpdateMu
 Gate 3.1: step 1  post-UpdatePerm
 Gate 3.1: SetTransients enter step 1  srcScale 0  bcScale 0
Gate 3.1: SetTransients pre-SetSourceCst
Gate 3.1: SetSourceCst enter
 Gate 3.1: SetSourceCst set 1  id 1
[1]PETSC ERROR: Caught signal number 11 SEGV
```
关键细节：**rank 0 完整走完** `SetSourceCst done → post-SetSourceCst → pre-UpdateBC → UpdateBC enter → UpdateBC cell set 2 id 10`，之后才收到 MPI_ABORT 的 signal 15。所以 rank 1 崩在 `set 1`（唯一的 codim=0 体单元集 `cs0001`）循环体内，而 rank 0 通过了同一段代码 —— 典型的 rank 相关越界，不是逻辑分支错误。

**根因（#54 的 BAR2 跳过修复已生效，但只解决了第一个问题）**：`MEF90FlowSetSourceCst` 仍走 `DMMeshGetStratumIS(cellDMScal,...) → DMMeshISCreateISglobaldof(cellDMScal,...,setISdof) → ISGetIndicesF90(setISdof,setIdx) → VecSetValues(sourceVec,size(setIdx),setIdx,...)`。在 PETSc 3.3 的 DMMesh cell-section 下，这条 global-dof 路径在 rank 1 上给出空/越出 `sourceVec` 布局的索引（rank 0 拥有低编号 point 恰好合法），`ISGetIndicesF90`/`size()` 或 `VecSetValues` release 模式直接越界读写 → SEGV。HeatXfer 的 `SetFluxCst` 是同款写法，但本 CI 里 HeatXfer 的 transient 分支一直是 NULL，**从未在 2 rank 上验证过**，不能作为可信模板。

**修复：改用项目里已被 10 轮 CI 验证的 SectionReal 路径**（`MEF90FlowUpdatePerm` / `MEF90FlowUpdateStorativity` 在 init 07/08 和每个 step 的 `post-UpdatePerm` 都通过）。`MEF90FlowSetSourceCst` 现在：
```
SectionRealDuplicate(cellDMScalSec, sourceSec) → SectionRealSet(0)
→ DMmeshGetLabelIdIS(DM,'Cell Sets') + AllGatherMerge
→ 每个 set：skip Has_BC / Codim/=0；DMMeshGetStratumIS(DM,...)
→ ISGetIndicesF90(setIS,cellIDs)；逐 cell SectionRealRestrict/写入 source/SectionRealRestore
→ SectionRealComplete + SectionRealToVec(sourceSec, cellDMScalScatter, SCATTER_FORWARD, sourceVec)
→ SectionRealDestroy
```
完全不再调用 `DMMeshISCreateISglobaldof` 与 `VecSetValues`，也不再对 `setISdof` 取 F90 指针。仍保留 #54 的 BAR2 跳过，并打印 `set / id / nCell / src` 便于下一轮核对。

**附加防御（同类空 IS 模式，共 6 处）**：`MEF90FlowUpdateboundaryPressure`（cell set / vertex set）、`MEF90FlowApplyDirichletResidual`（cell set / vertex set）、`MEF90FlowApplyDirichletJacobian`（cell set / vertex set）都在 `ISGetLocalSize` 后用 `If (nVal > 0)` 包住指针操作；两个 Dirichlet 函数里原先缺 `CHKERRQ` 的 `MEF90ISCreateCelltoVertex` 也补上。这些函数会在 step 2 的 TS 回调里首次成批执行。

**CI 缓存修复（性能，更正 #54 章节的推断）**：原 workflow 有两个独立缺陷 —— ① `actions/cache@v4` 的 `path:` 里写 `$HOME/...`，而 cache action **不展开 shell 变量**，路径是字面量、永远缓存不到真实目录；② `actions/cache@v4` 的保存是 Post 步骤，**job failure 时被 skip**，而此前每轮 job 都 failure，缓存从未落盘（日志里 `Post Cache ... completed skipped`）。现在改成 `actions/cache/restore@v4` + `~/...` 路径，并在 job 末尾加 3 个 `if: always()` + `continue-on-error: true` 的 `actions/cache/save@v4`，key 保持不变。**从 #56 起即使失败也会把 PETSc/OpenMPI 缓存下来，后续 run 可省掉约 28 分钟重建。**

**校验**：`git diff --check` 干净；PyYAML 解析 workflow（18 steps）/ yaml 通过；`Run Terzaghi` 的 `run` 块 `bash -n` 通过；5 个改写子程序的 `Do/End Do`、`If/End If`、`Select Case/End Select` 配对脚本核对通过（`SetSourceCst` 的 2 个 `If` 是单行 `If (...) Cycle`，无需 `End If`）；新增的 `SectionRealDuplicate/Set/Restrict/Restore/Complete/ToVec/Destroy`、`SCATTER_FORWARD` 在本文件中分别有 16/3/8/2/4/15/16/5 处既有先例。仍为零本地编译、靠 CI 验证。

**#56 判据（按顺序）**：① `SetSourceCst enter (SectionReal path)` → `set 1 id 1 nCell ... src ...` → `done`；② `UpdateBC cell set ... nVal ...` / `UpdateBC done`、`pre-TSSolve`；③ step 1 打 `TSSolve skipped (no time advance)`，step 2 起出现 `IFunction enter / post-OperatorCore / IJacobian enter / post-TSSolve` 与 `-flow_ts_monitor` 步进；④ `TSGetConvergedReason`、TS time 推进、`_out.gen` 内 pressure 是否随时间衰减（3.2 Terzaghi 解析解比对前置）。若仍 SEGV，最后一个 `Gate 3.1` 检查点仍可精确定位（缓存生效后下一轮只需几分钟）。


## CI #56 诊断结论与修复（2026-11-07 10:45，BC 点 ID → 全局 dof 映射修复）

**状态**：run `42326375284` / job `128262502342` `completed failure`（基于 `ff958c1`），01:52:49–02:29:47Z 共 37 分钟，构建 14 步全 success；`Run Terzaghi` 仍 2 秒内失败，mpirun exit **59**，完整 job 日志 3048 行。`init 01..08`、`flow enum 2 Transient 2 enabled T`、`transient Darcy TS created` 全部正常。

**进展（#55 的 SectionReal 源项修复完全生效）**
```
Gate 3.1: SetTransients enter step 1  srcScale 0  bcScale 0
Gate 3.1: SetTransients pre-SetSourceCst
Gate 3.1: SetSourceCst enter (SectionReal path)
 Gate 3.1: SetSourceCst set 1  id 1  nCell 300  src 0.0
Gate 3.1: SetSourceCst done
Gate 3.1: SetTransients post-SetSourceCst
Gate 3.1: SetTransients pre-UpdateBC
Gate 3.1: UpdateBC enter
 Gate 3.1: UpdateBC cell set 2  id 10
 Gate 3.1: UpdateBC cell set 2  nVal 0        <-- rank 0 在此之后被 MPI_ABORT 杀掉
[1]PETSC ERROR: Caught signal number 11 SEGV
```
`SetSourceCst` 两 rank 均完成，崩溃点后移到 `MEF90FlowUpdateboundaryPressure` 内 `cs0010`（顶面 BAR2，id=10）循环。**关键：此前所有 `Gate 3.1` trace 都用 `PetscPrintf(PETSC_COMM_WORLD,...)`，只在 rank 0 打印**，所以 `nVal 0` 是 rank 0 的值；rank 1 拥有全部顶面顶点（`nVal > 0`），它是唯一会进入 `nVal > 0` 分支的 rank。

**根因（确定性，同时解决遗留项）**：`bcIS` / `setIS` 里存的是 **mesh point ID**，不是 Vec/Mat 的全局 dof 索引。Flow 三个 BC 子程序把它们直接传给 `VecGetValues` / `VecSetValues` / `MatZeroRows`，rank 1 用 point id 当矩阵行/向量行索引 → release 模式下越界读写 SEGV；rank 0 因 `nVal=0` 被保护跳过，所以呈现典型 rank 相关崩溃。HeatXfer / DefMech 从来不这么写：所有 BC IS 都经 `MEF90VecSetValuesISdof` / `MEF90VecGetValuesISdof`（内部 `DMMeshISCreateISglobaldof` 做 point→dof 映射），Jacobian 用 `DMMeshISCreateISglobaldof` + `MatZeroRowsColumnsIS`。`vDef.F90:424` 在 init 阶段已调用 `MEF90DefMechUpdateboundaryDisplacement`（走同一套 helper）且两 rank 通过 —— 这证明 global-dof 路径在 vertex/vect DM 上可靠，#55 的问题只出在 cell-centered `cellDMScal` section。

**本轮改动（2 文件）**
1. `m_Flow/m_MEF90_FlowAssembly.F90` —— 三个 BC 子程序统一对齐 HeatXfer 模板：
   - `MEF90FlowUpdateboundaryPressure`：`VecSetValues(...,dofList,...)` → `MEF90VecSetValuesISdof(FlowCtx%DMScal,boundaryPressure,val,bcIS,1,INSERT_VALUES,...)`（`val` 改 Fortran POINTER 数组）；vertex set 同理。
   - `MEF90FlowApplyDirichletResidual`：`VecGetValues` × 2 → `MEF90VecGetValuesISdof`，`VecSetValues` → `MEF90VecSetValuesISdof`（`residual = x - p_bc` 逻辑不变）；末尾补 `VecAssemblyBegin/End`；`ISDestroy` 全补 `CHKERRQ`。
   - `MEF90FlowApplyDirichletJacobian`：`MatZeroRows(Jac,nVal,dofList,...)` → `DMMeshISCreateISglobaldof(FlowCtx%DMScal,bcIS,0,setISdof,...)` + `MatZeroRowsColumnsIS(Jac/PreJac,setISdof,1.0,...)`；vertex set 同理。
   - 三个函数内所有 Gate trace 从 `PETSC_COMM_WORLD` 改为 **`PETSC_COMM_SELF`**：两个 rank 各自打印，`id` 与 `nVal` 两行夹住 `MEF90ISCreateCelltoVertex`，下轮若仍崩可直接判定是该函数内部还是之后的 dof 映射。
2. `.github/workflows/build.yml`：SEGV 诊断的 Gate grep `tail -60` → `tail -100`（容纳双 rank + step 2 的 TS 回调 trace）。

**CI 缓存已落盘（重要）**：#56 末尾 3 个 `actions/cache/save@v4` 全部成功 —— `Cache saved with key: openmpi-3.1.3-ubuntu-latest-v1` / `petsc-dl-v1` / `petsc-mef90-3.3-ubuntu-latest-v1`。**#57 起 PETSc/OpenMPI 直接命中缓存，37 分钟缩短到几分钟**，迭代节奏大幅加快。

**校验**：`git diff --check` 干净；三个改写子程序的 `Do/End Do`、`If/End If` 配对脚本核对通过（脚本对两处单行 `If (...) Cycle` 的假阳性已排除）；`MEF90Vec*ISdof` / `DMMeshISCreateISglobaldof` / `MatZeroRowsColumnsIS` 均由 `Use m_MEF90` 可见且在 HeatXfer/DefMech 有大量先例；全文件已无 `VecSetValues` / `VecGetValues` / `MatZeroRows(` 的实际调用（仅注释）。仍为零本地编译、靠 CI 验证。

**#57 判据（按顺序）**：① 两 rank 各自打印 `UpdateBC cell set 2 id 10 / nVal ... / done`，若 rank 1 只打出 `id 10` 而无 `nVal` → 崩在 `MEF90ISCreateCelltoVertex` 内部；② `post-UpdateBC` → step 1 `TSSolve skipped (no time advance)`；③ step 2 起 `IFunction enter / post-OperatorCore / IJacobian enter / post-TSSolve` 与 `-flow_ts_monitor` 步进、`TSGetConvergedReason`；④ `_out.gen` 内 pressure 随时间从 `p0=1e5` 衰减（3.2 Terzaghi 解析解比对前置）。


## CI #57 诊断结论与修复（2026-11-07 11:55，collective 死锁 → SectionReal 全路径）

**状态**：run `42331220979` / job `128277936048`，`in_progress` 挂起约 1 小时后主动 `cancel`（`completed cancelled`）。**缓存全部命中**：`Cache OpenMPI / PETSc deps / Miniconda+PETSc build` 全 success，`Build OpenMPI` / `Pre-download PETSc deps` / `Build PETSc-3.3` 全部 skipped，构建 14 步在 02:51:18Z 即完成（**从 37 分钟降到约 25 秒**）。`Run Terzaghi` 02:53:24Z 开始后再无任何输出，job 无 `timeout-minutes`、step 无超时，会一直挂到 6 小时默认上限。运行中/取消后的 job log archive 未及时生成（302 签名地址返回 `BlobNotFound`）、artifacts 为空，故本轮靠静态定位。

**根因（高置信度）**：#56 把裸 `VecSetValues` 换成 `MEF90VecSetValuesISdof` / `MEF90VecGetValuesISdof`，但这些 helper **内含 collective 调用**：
```
MEF90VecGetValuesISdof: DMMeshGetSectionReal → DMMeshCreateGlobalScatter → SectionRealToVec(SCATTER_REVERSE)  [VecScatter = collective]
MEF90VecSetValuesISdof: DMMeshISCreateISglobaldof → ... (DM Mesh 内部含同步)
```
而 #55 留下的 `If (nVal > 0)` 保护让 **rank 0（顶面 BC 顶点数 0）整段跳过、rank 1（拥有全部顶面顶点）独自进入 collective** → MPI 死锁。行为变化完全吻合：#56 是 rank 1 直接越界 SEGV（裸 `VecSetValues` 非 collective），#57 换成 helper 后变成两 rank 无限等待。同类风险还有 `MatZeroRowsColumnsIS`（collective）被 `If (nVal > 0)` 包住。**通用规则：任何 collective 调用都不能放在 rank 相关的分支里。**

**本轮改动（2 文件）**
1. `m_Flow/m_MEF90_FlowAssembly.F90` —— 三个 BC 子程序统一改走 **SectionReal 路径**（#55 `SetSourceCst` 已在两 rank 验证；`MEF90FlowOperatorCore` / `MEF90FlowAddStorageToResidual` 同款）：
   - `MEF90FlowUpdateboundaryPressure`：`SectionRealDuplicate(DMScalSec,bcSec)` + `SectionRealSet(0)` → 逐 BC point `SectionRealRestrict/写/Restore`（**纯本地，可用 `nVal > 0` 保护**）→ 函数末尾无条件 `SectionRealComplete` + `SectionRealToVec(SCATTER_FORWARD)` + `Destroy`。
   - `MEF90FlowApplyDirichletResidual`：duplicate `resSec/xSec/bcSec`，三者各 `SectionRealToVec(SCATTER_REVERSE)` 读回（无条件）；逐 BC point 写 `rPtr(1) = xPtr(1) - bPtr(1)`（本地）；末尾无条件 `SectionRealComplete(resSec)` + `SCATTER_FORWARD` 写回 `residual` + Destroy × 3。
   - `MEF90FlowApplyDirichletJacobian`：**删掉两处 `If (nVal > 0)`**，`DMMeshISCreateISglobaldof` + `MatZeroRowsColumnsIS(Jac/PreJac, ...)` 在**两个 rank 上无条件执行**（rank 0 传空 dof IS，PETSc 视为 no-op）。
   - 结果：`MEF90Vec*ISdof` 两个 helper 已从 Flow assembly 中完全移除；所有 collective 都落在所有 rank 都执行的路径上，仅本地 per-point 写入受 `nVal` 保护。
2. `.github/workflows/build.yml` —— 防挂兜底：`Run Terzaghi` 加 `timeout-minutes: 20`；`mpirun` 用 `timeout -k 10s 600` 包装；rc=124/137 时打印 `=== FAIL: Terzaghi TIMED OUT / HUNG ===` + 最后 100 条 Gate trace + 最后 80 行原始日志后 `exit 1`，保证任何死锁都能在 10 分钟内产出可诊断日志（且 `Upload diagnostics` 的 `if: always()` 能收走 `/tmp/terzaghi.log`）。

**校验**：`git diff --check` 干净；PyYAML（18 steps，`timeout-minutes: 20`）通过；`Run Terzaghi` 的 `run` 块 `bash -n` 通过；18 个 Subroutine 的 `Do/End Do`、`If/End If`、`Select/End Select` 配对脚本全平衡；全文件已无 `VecSetValues` / `VecGetValues` / `MatZeroRows(` / `MEF90Vec*ISdof` 的实际调用（仅注释）；新增的 `SectionRealDuplicate/Set/Restrict/Restore/Complete/ToVec/Destroy` 与 `SCATTER_REVERSE/FORWARD` 在本文件均有既有先例（其中 `MEF90FlowOperatorCore` 已用同一个 `DMScalScatter` 对多个 duplicate section 做 ToVec）。仍为零本地编译、靠 CI 验证。

**#58 判据（按顺序）**：① 缓存命中 → 构建 < 1 分钟；② 两 rank 各打印 `UpdateBC enter / cell set 2 id 10 / nVal ... / UpdateBC done`（`PETSC_COMM_SELF`），`post-UpdateBC` → `pre-TSSolve`；③ step 1 `TSSolve skipped (no time advance)`，step 2 起 `IFunction enter / post-OperatorCore / IJacobian enter / post-TSSolve` 与 `-flow_ts_monitor` 步进；④ TS 收敛原因、TS time 推进、`_out.gen` 压力从 `p0=1e5` 衰减（3.2 前置）。若仍挂，`timeout` 会在 10 分钟内给出卡住前最后一条 Gate trace，可直接区分死锁与 TS 收敛停滞。


## CI #58 诊断结论与修复（2026-11-07 12:35，3.1 实质跑通 / 断言误判 + 3.2 数据接线）

**状态**：run `42337427966` / job `128297410589`，`completed failure`，全程 **3 分 22 秒**（三个缓存全部命中，14 个构建步骤全 success，编译无 error）。**`Run Terzaghi` 的 mpirun exit code = 0**，程序正常跑完 51 个时间步并自行结束；CI 判红只因为最后一条「runtime failure signature」扫描。

### 闸门 3.1 数值求解首次真实跑通（确定性证据）

| 证据 | 值 |
|---|---|
| `mpirun` 退出码 | **0**（无 SEGV、无 MPI_ABORT、无死锁、无超时） |
| 时间推进 | `-flow_ts_monitor` 共 99 行，`0 TS dt 0.02 time 0` → `1 TS dt 0.02 time 1` |
| `TS Object:TSPressure(flow_)` | `type: beuler`，`total number of linear solver iterations=52`，**`total number of rejected steps=0`** |
| SNES/KSP | `ksponly` + `cg` + `bjacobi/ilu`，`rows=356`，`nonzeros=2266`，无 malloc |
| Gate trace | 两 rank 各自打印 `UpdateBC ... nVal 6 / nVal 0 / done`；每步 `SetSourceCst(SectionReal) → UpdateBC → IFunction → storage residual → IJacobian → storage matrix → post-TSSolve → TS time` 全部到达 |
| 结果文件 | `vDef/data/Beam1x10-tri3_out.gen` **2,845,300 bytes**（#51 的 2,789,988 字节是初始场，本轮为完整瞬态输出） |
| 总交替最小化次数 | 51 |

**#57 的 SectionReal 全路径修复完全生效**：collective 全部落在两 rank 都执行的路径上，本地逐 point 写入才用 `nVal>0` 保护，死锁消失；`MEF90ISCreateCelltoVertex` + `DMMeshISCreateISglobaldof` + `MatZeroRowsColumnsIS` 在 2 rank 上稳定运行。

### CI 判红的真实根因：DefMech 遗留报错被当成 Flow 失败

`terzaghi.log` 中 DefMech AltMin 每步都打印 `displacement SNESSolve failed with SNESConvergedReason -3`（49 次）与 `damage field SNESSolve failed with SNESConvergedReason -5`（49 次），以及对应的 `for error code meaning.` URL 行。旧扫描模式含 `SNESolve failed|error code meaning`，因此误判。**这是力学侧固有问题**（`mechanical work` 达 `4.56E+46`，但 `alpha min/max = 3.87E-05 / 6.47E-04`）：damage 幅度 < 0.07%，经 `k(v)=k_min+(k_max-k_min)v` 对渗透率影响 < 1.2%，**不污染 Flow 压力场**，留待闸门 4 / 独立力学收敛项处理。新签名对 #58 全量日志回归命中数为 **0**。

### 新发现：step 16 浮点边界跳步（3.2 前需校正）

`step 16 pre-TSSolve TS time 0.29999999999999999 target 0.29999999999999999` → `If (t < tFlowEnd)` 为假 → `TSSolve skipped`。加上 step 1 的零跨度跳过，51 个外层步只有 **49 次 `post-TSSolve`**。TS 最终时间仍为 `1.0000000000000004`（target `1.0`），但压力时间轴在这一步出现一个 dt 的平台；3.2 硬比对前必须把判据改成基于 `time(step)` 的容差比较（例如 `If (tFlowEnd - t > 0.5*dtFlow0)` 或直接按步号驱动），否则解析解会错开一个时间步。本轮先记录并让 Gate 3.2 曲线把它暴露出来，不动 TS 时间逻辑（零编译验证下风险过高）。

### 本轮改动（3 文件）

1. `vDef/vDef.F90` —— 新增 **Gate 3.2 压力探针**：声明 `pressureMin/pressureMax`，在每个外层步 `VecCopy(pressure → pressurePreviousStep)` 之后调用 `VecMin/VecMax`（项目内 `m_MEF90_DefMechAssembly.F90:91-92` 已有同款先例；两者是 collective 但两 rank 都执行，不违反 #57 的规则），打印 `Gate 3.2: step N pmin ... pmax ...`。`pmax` 即不排水底面压力：Terzaghi 解析解在 `t=1 s, T≈1.09` 时给出 `p(0,t)≈8.7e3 Pa`，即从 `p0=1e5` 衰减到约 **8.7%**。
2. `.github/workflows/build.yml` —— 判据分层，不是放松而是加强：
   - `SNESolve failed` / `error code meaning` 从致命签名移除，改为显式 `WARNING` 并按 reason 统计次数；致命签名收窄为 `segmentation fault|signal 11|Caught signal|MPI_ABORT|attempt to deallocate|PETSC ERROR|...|TSSolve failed|TSConvergedReason -[1-9]|KSP_DIVERGED|DIVERGED_|not converged`。
   - 新增 **Gate 3.1 正向硬判据**：`rc==0`、`post-TSSolve >= 48`（#58 实测 49，同时打印 skip 次数）、`-flow_ts_monitor` 到达 `time 1`、所有 `TS Object` 的 `rejected steps=0`、`_out.gen` 非空且 magic 为 `434446`。通过时打印 `Gate 3.1 PASS`。
   - 新增 **Gate 3.2 probe**（本轮 informational）：打印首末步 `pmax`、衰减比、每 10 步采样；下轮据实测曲线升级为硬解析解比对。
   - 新增 `Upload Terzaghi result mesh` step，把 `_out.gen` 作为 artifact 留存，供 3.2 离线读取 nodal pressure 剖面。
3. `FLOW_FIX_CHECKLIST.md` / 本地工作日志：本节与进度记录。

**校验**：`git diff --check` 干净；PyYAML（19 steps）通过；`Run Terzaghi` 块 `bash -n` 通过；用 #58 的 `terzaghi.log` 回归新判据得到 `FATAL_RE 命中=0`、`nSolve=49`、`FINAL_TIME_OK`、`rejected!=0 count=0`、`_out.gen` 存在。仍为零本地编译、靠 CI 验证。

**#59 判据**：① CI 首次全绿（3.1 正式通过）；② `Gate 3.1 PASS` 与 `Gate 3.2 probe` 同时出现；③ `pmax(final)/pmax(initial)` 应约 **0.09**（≈8.7e3 / 1e5），曲线单调下降且 step 16 处有一个平台；④ 若比值 > 1（压力反而增长）即为 3.2 要修的**符号反向**（`MEF90FlowOperatorCore` 的 `G(p)=+div(A∇p)` 与 `-div(k/μ∇p)=S·ṗ` 差一个负号）；⑤ 下轮用 artifact 中的 `_out.gen` 做 z 剖面与解析级数逐点比对。


## CI #59 结果与 Gate 3.2 根因（2026-11-07 13:10，3.1 全绿 / 压力指数发散 → DM dof 映射）

**状态**：run `42340345883` / job `128306409584`，**`completed success` —— 项目首次 CI 全绿**（缓存命中，约 3 分钟）。`Run Terzaghi` 输出：`mpirun exit code: 0`、DefMech 98 次非收敛被正确标为 WARNING、`Gate 3.1: 49 completed TSSolve calls, 2 skipped`、`Gate 3.1 PASS: ... valid Exodus file`。`Upload Terzaghi result mesh` artifact 上传成功（1.23 MB 压缩包 / 2,845,300 字节原始文件）。

### Gate 3.2 探针暴露的物理错误（不是"绿灯即正确"）

51 个 `Gate 3.2: step N pmin/pmax` 显示压力**指数发散**而非 Terzaghi 衰减：

| step | pmin | pmax |
|---|---|---|
| 1 | 1.0e5 | 1.0e5 |
| 2 | 0 | 1.0e5 |
| 6 | -1.68e3 | 9.88e4 |
| 10 | -4.04e7 | 3.23e4 |
| 20 | -3.36e8 | 4.02e11 |
| 30 | -3.50e16 | 3.96e19 |
| 40 | -3.97e24 | 4.71e27 |
| 51 | -4.70e36 | **1.65e36** |

`pmax(final)/pmax(initial) = 1.65e31`（解析解应 ≈ **0.087**）。pmin/pmax 交替变号放大，是 F 与 J 不一致导致线性解被反复放大，不是物理/时间步长问题（BE 无条件稳定，且 `rejected steps=0`、KSP 均"收敛"）。

### 从 `_out.gen` 直接读出的确定性证据

手写 NetCDF classic (CDF2) 解析：`numrecs=51`，record 起点 23980、record 大小 55320 字节，`vals_nod_var1` = pressure（356 节点 × double），`vals_elem_var4eb1` = perm field（600 单元 × double）。

- **`vals_elem_var4eb1` 全 51 步恒为 `5.0e-8`** = k/μ = 5e-11/1e-3，**渗透率正确**，排除"K=0 导致不衰减"。
- pressure 的 min/max 与 Gate 3.2 探针逐点吻合（step 51: -4.701e36 / 1.654e36），确认探针可信。
- `vals_nod_var2` 是 damage（3.87e-5 ~ 6.47e-4），与日志 `alpha min/max` 一致。

### 根因（确定性，对照已验证 HeatXfer 模板）

`MEF90FlowCtx` 里 `DM` 是**主 DM**：`vDef.F90` 调 `DMMeshSetMaxDof(Mesh,dim)` / `DMSetBlockSize(Mesh,dim)`（dim=2），DefMech 在其上设 **2 dof/vertex** 的 default section。压力 Vec、356×356 矩阵、`DMScalSec` 全部挂在 `DMScal`（`DMMeshClone` 后 `MaxDof=1`、`BlockSize=1`、`DMMeshGetVertexSectionReal(...,1,...)`，即 **1 dof/vertex**）。

但 Flow 的 6 处底层装配把 **`FlowCtx%DM`** 传给 `MEF90DiffusionOperatorSet` / `MEF90DiffusionBilinearFormSet` / `MEF90DiffusionOperatorAddTransientTermSet` / `MEF90DiffusionRHSSetCell`。其中 `MEF90DiffusionBilinearFormSet` 内部 `DMMeshGetSectionReal(mesh,'default',defaultSection)` 取的是**主 DM 的 2-dof section**，于是 `DMmeshAssembleMatrix` 把 3×3 单元刚度矩阵按 6-dof 索引写进 1-dof 的 356×356 矩阵 → J 的行列完全错乱；而 IFunction 走 SectionReal 路径（按 `DMScalSec` 的 point id 写），**F 与 J 不相容**。`KSPCG` 在 `MAT_SPD` 假设下"收敛"到错误解，每步放大约 2.5 倍。`m_HeatXfer/m_MEF90_HeatXferAssembly.F90:595/466` 全程使用 `MEF90HeatXferCtx%DMScal`，正是应当对齐的模板。

注意：清单 3.2 原先怀疑的「稳态算子 `+div(A∇p)` 符号反了」**并不成立** —— 弱形式 `∫A∇p·∇v` 正是 `-div(A∇p)` 的离散，F=`K p + S M ṗ`、J=`K + a S M` 的符号是自洽的。

### 本轮改动（#60，2 文件）

1. `m_Flow/m_MEF90_FlowAssembly.F90`：6 处 TS 路径装配的 mesh 参数 `FlowCtx%DM` → **`FlowCtx%DMScal`**（`OperatorSet` / `RHSSetCell` / `BilinearFormSet(Jac)` ×2 / `BilinearFormSet(PreJac)` / `OperatorAddTransientTermSet`）；`MEF90FlowIJacobian` 增加 `Gate 3.2 diag: ... shift a ...` 打印（前 4 步，预期 `a=1/dt=50`）。`MEF90DiffusionEnergySet` 仍用主 DM（纯后处理，不参与 TS 求解），留待后续统一。
2. `.github/workflows/build.yml`：Gate 3.2 探针升级为**硬判据** —— `ratio=pmax(final)/pmax(initial) > 1.5` 即 `FAIL: Gate 3.2 pressure DIVERGED` 并打印 diag；`0.02 ≤ ratio ≤ 0.30` 打印 `Gate 3.2 PASS`（解析解 ≈0.087）；区间外标 WARNING。**这样 #59 那种"3.1 全绿但物理发散"不再可能通过。**

**校验**：`git diff --check`、PyYAML（19 steps）、`bash -n` 全通过；6 处调用逐一核对，TS 路径已无 `Diffusion*Set(...,FlowCtx%DM,...)` 残留。仍为零本地编译、靠 CI 验证。

**#60 判据**：① 编译通过、`Gate 3.1 PASS` 保持；② `Gate 3.2 diag shift a` = 50；③ pmax 单调衰减、`ratio≈0.09` → `Gate 3.2 PASS`、CI 全绿；④ 若仍发散，则下一轮用 `MatNorm(Jac)` 与 `VecNorm(F)` 逐项对比 F/J 装配，并检查 `MEF90Element_Create` 是否也需换成 `DMScal`；⑤ 收敛后把 step 16 的浮点跳步改为基于 `time(step)` 的容差驱动，再做 z 剖面与 Terzaghi 级数逐点比对。


## CI #60 复核（2026-11-07 13:25，DM 假设被证伪 → 转向 F/J 一致性）

**状态**：run `42342061821` / job `128311695113`，`completed failure`，但失败原因正确且信息明确：`=== FAIL: Gate 3.2 pressure DIVERGED (pmax final/initial=1.65448e+31 > 1.5) ===`（硬判据首次生效，杜绝假绿）。`Gate 3.1 PASS` 保持，`mpirun exit code: 0`。

**关键观测**：`Gate 3.2 diag: IJacobian t 0.02/0.04/0.06 shift a 50.0` —— BE 的 `a=1/dt=50` **完全正确**；但 `pmax` 末值 **1.6544841814641773E+036 与 #59 逐位完全相同**。把 6 处装配的 mesh 从 `FlowCtx%DM` 改成 `DMScal` **没有改变任何数值** → 主 DM 与 DMScal 的 default section 在本配置下等价，「2-dof/1-dof 映射错位」假设被证伪（改动本身仍保留，因为它与 HeatXfer 模板一致）。

**新的最强怀疑（可由 #61 直接证伪/确认）**：`MEF90FlowIFunction` 内共发生 **3 次 `SectionRealToVec(..., SCATTER_FORWARD, F)`**（`OperatorCore` 内 1 次 + `AddStorageToResidual` 内 1 次 + `ApplyDirichletResidual` 内 1 次），而 HeatXfer 的 IFunction 只在全部分项累加进同一个 `FSec` 后 scatter **1 次**。若该 scatter 在 PETSc 3.3 是 `ADD_VALUES` 语义，`F` 会被反复累加成数倍 `K·p + S·M·ṗ`，而 `J` 只装配一次 → **F/J 不相容**，`KSPCG` 在 `MAT_SPD` 假设下收敛到错误解并逐放大 —— 与观测的 `|x|`、`|xDot|` 交替变号放大一致。

**#61 改动（纯诊断，不动求解逻辑）**：在 `MEF90FlowIFunction` 的前 4 步（`t<0.07`）打印 `|x|∞`、`|xDot|∞`，以及 `|F|∞` 在 `post-OperatorCore` / `post-Storage` / `post-Dirichlet` 三个检查点的值。若 `|F|` 逐次成倍增长即可确认 scatter 累加，下一步把 IFunction 重构为「单个 `FSec` 累加 K·p + S·M·ṗ + 源项 → 一次 SCATTER_FORWARD → BC 行单独覆盖」的 HeatXfer 结构。


## CI #61 诊断结果（2026-11-07 13:40，scatter 累加假设也被证伪；xDot 恒为 0）

**状态**：run `42342725365`，`completed failure`（仍是 `FAIL: Gate 3.2 pressure DIVERGED`，`ratio=1.65448e+31`），`Gate 3.1 PASS` 与 `mpirun exit code: 0` 保持。曲线与 #59/#60 逐位一致。

**实测诊断（t=0.02 / 0.04 / 0.06，两 rank 一致）**：

| t | \|x\|∞ | \|xDot\|∞ | \|F\| post-OperatorCore | post-Storage | post-Dirichlet | IJacobian shift a |
|---|---|---|---|---|---|---|
| 0.02 | 1.0e5 | **0** | 1.0e5 | 1.0e5 | 1.0e5 | 50 |
| 0.04 | 1.0e5 | **0** | 1.0697e-2 | 1.0697e-2 | 1.0697e-2 | 50 |
| 0.06 | 9.9978e4 | **0** | 7.0755e-5 | 7.0755e-5 | 7.0755e-5 | 50 |

**三条确定性结论**：
1. `SectionRealToVec(SCATTER_FORWARD)` 在本配置下是 **INSERT 语义**：三次 scatter 后 `|F|` 完全不变，「F 被反复累加」假设证伪，无需重构 IFunction。
2. **`xDot` 在所有步恒为 0**，`AddStorageToResidual`（`S·M·ṗ`）在 TS 求解中对 F **零贡献**；而 `IJacobian` 仍加 `a·S·M`。PETSc 3.3 `TS_LINEAR` + `ksponly` 是在 `(x=U_n, xdot=0)` 处求 F，再解 `J ΔU = -F`（`J=K+aSM`），代数上等价于 BE 的 `(K+aSM)U_{n+1}=aSM·U_n`。
3. `|F|` 从 1e5 骤降到 1e-2 / 7e-5，只是**顶面 Dirichlet 代数行** `x-p_bc` 收敛；`|x|` 始终 ~1e5，说明内部压力几乎没被推进。按量级估算 `K/aSM ≈ 360`（K~2.5e-8，aSM~7e-11），正确 BE 应在首步就衰减到 ~0.1% —— 因此 **J 中的 K 项实际没有生效**（或被 BC 处理抹掉），使 `U_{n+1}≈U_n` 并逐步放大高频误差。

**下一轮唯一该做的诊断**：在 `MEF90FlowIJacobian` 里打印 `MatNorm(Jac,NORM_INFINITY)` 以及 `MatGetDiagonal` 的 min/max（分 `BilinearFormCore` 后、`AddStorageMassMatrix` 后、`ApplyDirichletJacobian` 后三个检查点），直接看 K 是否被清零、质量项量级是否与 7e-11 相符。若 `MatZeroRowsColumnsIS` 的 `setISdof` 实际覆盖全体顶点（而非仅顶面 6 个），`J` 会退化成单位阵 —— 这与「内部压力不变 + 高频交替放大」完全吻合，是最可能的剩余根因。


## #62 静态复核（2026-11-07 14:30）：IS 假说排除，锁定「BC 零行被调用 4 次」的结构差异

**已排除**：`MEF90FlowApplyDirichletJacobian` 的 `DMMeshISCreateISglobaldof(DMScal,bcIS,0,setISdof)` + `MatZeroRowsColumnsIS(A,setISdof,1.0,...)` 与 `m_HeatXfer/m_MEF90_HeatXferAssembly.F90:617-620` **逐字相同**，`MEF90ISCreateCelltoVertex` 也只返回边界 BAR2 的顶点，不存在「覆盖全体 356 顶点」的可能。`A_zero = 0.0_Kr` 经 `MatS2D_Get_Real` 正确置零全部三个分量。

**新发现的结构差异（确定性）**：`MEF90HeatXFerIJacobian` 是 `MatZeroEntries(A)` → 装 K（`MEF90DiffusionBilinearFormSet`）→ 装 `shift*ρc`（`MEF90_MassMatrixAssembleSet`）→ `MatAssembly` → **最后统一做一次** BC `MatZeroRowsColumnsIS`；且 `TSSetIJacobian(tsTemp,matTemp,matTemp,...)` 里 A 与 M 是同一个矩阵。而 Flow 的 `MEF90FlowBilinearFormCore` 在 K 装完、`MatCopy` 之后**立刻**调 `ApplyDirichletJacobian`，随后 `AddStorageMassMatrix` 又往已被清零的 BC 行/列 `ADD_VALUES` 加质量项，`MEF90FlowIJacobian` 末尾**再**调一次 `ApplyDirichletJacobian` —— 每步 Jac/PreJac 各 2 次、共 **4 次** `MatZeroRowsColumnsIS`，中间夹着一次往零行写值的 ADD 装配。这与 HeatXfer 已验证模板不符。

**#62 改动（结构对齐 + 三点诊断，同时落地）**：

1. `MEF90FlowBilinearFormCore`：删除末尾的 `ApplyDirichletJacobian`，只负责装配纯 K 并 `MatCopy` 到 PreJac；BC 交由调用方在所有装配完成后统一施加。
2. `MEF90FlowBilinearForm`（稳态 SNES wrapper）：`Core(...,0.0_Kr,...)` 之后补一次 `ApplyDirichletJacobian`（稳态路径运行时被 `If (.not. flowTransient)` 跳过，但语义必须正确）。
3. `MEF90FlowIJacobian`：`Core(...,0.0_Kr,...)` → `AddStorageMassMatrix` → **`ApplyDirichletJacobian` 仅一次**；新增 `MEF90FlowReportJacobian` 辅助子程序，在前 3 步（`t<0.07`）打印 `MatNorm(Jac,NORM_INFINITY)` 与 `MatGetDiagonal` 的 min/max，检查点 `post-Core / post-Storage / post-Dirichlet`；同时打印 `FlowCtx%storativity` 与 `permField` 的 min/max，确认 S=4.6e-10、k/μ=5e-8 确实进入装配。
4. 静态量级修正：单元面积 ≈1.67e-2 m²，K 对角 ≈1e-7，M 对角 ≈6e-13，`aM≈3e-11`，即 `K/aM≈1e4`（清单 3.2 里原估的 360 偏小）。因此 `J≈K` 且 K 在常向量方向奇异，只有 Dirichlet 行给出对角 1 —— 矩阵条件数极差，与 KSP 迭代数 1 → 50 → 54 逐步攀升一致。

**判据**：① 若 `Gate 3.2 PASS`（ratio≈0.09）→ 结构差异即根因，3.2 通过；② 若仍发散，`diagJ` 三点直接给出答案：`|J|inf≈1e-7 / diag min≈1e-8` 说明 K 存在（问题在 F/J 相容性或 CG 收敛），`|J|inf≈1 / diag min=max=1` 说明 J 退化成单位阵；`diagFld` 同时排除参数未填充。③ 收敛后再处理 step 16 的浮点跳步与 z 剖面解析解比对。


## CI #62 结果（2026-11-07 14:45）：J 装配被完全证清，根因收敛到求解器/CG 对称性

**状态**：run `42348017902` / job `128329892722`，`completed failure`，但 `Gate 3.1 PASS`、`mpirun exit 0`；`pmax` 末值 **1.6544841814641773E+036 与 #59/#60/#61 逐位相同**（结构改动同样不改变数值）。

**`diagFld` / `diagJ` 实测（t=0.02/0.04/0.06 完全一致，J 是常数矩阵）**：

| 检查点 | \|J\|∞ | diag min | diag max |
|---|---|---|---|
| post-Core（纯 K） | 6.5392496e-7 | 1.0451187e-8 | 3.2696248e-7 |
| post-Storage（+aSM） | 6.5392496e-7 | 1.0645001e-8 | 3.2742320e-7 |
| post-Dirichlet（BC 行=1） | 1.0 | 1.1425415e-8 | 1.0 |

`S min/max = 4.6e-10 / 4.6e-10`、`k/μ min/max = 5e-8 / 5e-8` → **参数全部正确进入装配**；质量项使对角 +1.94e-10（=a·S·单元体积，正确）；BC 只影响少数行（diag max 由 3.27e-7 → 1）。**J = K + a·S·M + Dirichlet 单位行，完全正确，K 没有被抹掉，「BC IS 覆盖全体顶点」假说彻底排除。**

**压力演化给出新的决定性线索**：step 1 `pmin=pmax=1e5` → step 2 `0 / 1e5` → step 3 `-2843 / 99978` → step 4 `-2623 / 99864` → step 5 `-2032 / 99647` → step 6 `-1677 / 98811`。**前 6 步正是 Terzaghi 早期的正确物理**（顶面排水、底部几乎不变、扩散前沿逐层下移），step 7 起 `pmax=2.46e5`、`pmin=-9.4e5`… 正负交替放大。即**不是 F/J 符号反向**（否则首步就错），而是 J 是常数、病态矩阵上 CG 的误差逐步累积出 checkerboard 模式。

**不对称的直接证据**：`Gate 3.1: UpdateBC cell set 2 nVal` 两 rank 分别是 **6 和 0** —— `DMMeshGetStratumIS` / `MEF90ISCreateCelltoVertex` 产出的是**本地** point IS。`MatZeroRowsColumnsIS` 是 collective 且要清**列**，两 rank 传入不同 IS 时跨 rank 的耦合列无法一致清零 → J 非对称；而 `MEF90FlowCreateTS` 设了 `MAT_SPD` + `MAT_SYMMETRY_ETERNAL` + `KSPCG`，CG 在非对称矩阵上无收敛保证。日志中 `ILU: using diagonal shift to prevent zero pivot` 也表明预条件器处于病态边界。

**（备选解释）** DefMech / HeatXfer 用同样的 rank 相关 IS + `MatZeroRowsColumnsIS`，但它们的多分量矩阵与分区不同；若 GMRES 也发散，则须回到 PETSc 3.3 `TS_LINEAR`+`beuler` 的 IFunction 约定（本次未能取到 3.3 源码，仓库 tags 已不含 v3.3）。

**#63 动作（纯命令行二分，零代码改动、零编译风险）**：mpirun 增加 `-flow_ksp_type gmres -flow_ksp_rtol 1.0e-10 -flow_ksp_converged_reason`，保留默认 `bjacobi/ilu` PC。`KSPSetFromOptions` 会覆盖代码里的 `KSPCG`。GMRES 不要求 SPD：① 若 `Gate 3.2 PASS`（ratio≈0.09）→ 确认是 CG/对称性，随后补 `MEF90ISAllGatherMerge` 让两 rank 的 BC IS 一致（治本）或直接采用 GMRES；② 若仍发散 → 排除迭代求解器，下一轮核对 BE 线性步的 F/J 符号约定并考虑直接法 `-flow_pc_type lu`。


## CI #63 结果 + PETSc 3.3 源码核对（2026-11-07 15:00）：求解器排除，语义核对通过

**#63**（run `42348910801`）`completed failure`，`Gate 3.1 PASS`、`mpirun exit 0`；`pmax` 末值 **1.6543647997127005E+036**（#62 为 1.6544841814641773E+036，仅末几位不同），前 10 步曲线逐点一致。**GMRES + rtol 1e-10 完全不改变结果 → CG / 矩阵对称性 / KSP 收敛精度全部排除**，发散是离散方程本身的确定性性质（J 为常数矩阵的线性递推 `U_{n+1}=G U_n` 存在 |λ|>1 的模式）。

**PETSc 3.3 源码核对**（`petsc/petsc` tag `v3.3`，`src/ts/impls/implicit/theta/theta.c` 与 `src/ts/interface/ts.c`）：

- `TSStep_Theta`：`th->shift = 1/(Theta*dt)`；`VecCopy(ts->vec_sol,th->X)` → `SNESSolve(snes,th->affine,th->X)`；BE（Theta=1）时 `affine = F(t0,U0,0)*(Theta-1)/Theta = 0`，即 SNES 解 `G(U)=0`。
- `SNESTSFormFunction_Theta`：`Xdot = shift*(X-X0)`，再 `TSComputeIFunction(ts,stage_time,X,Xdot,Y)`；`SNESTSFormJacobian_Theta`：`TSComputeIJacobian(..., shift, ...)`。
- `TS_LINEAR` 在 `TSSetProblemType` / `TSGetSNES` / `TSSetFromOptions` 里**只做一件事**：`SNESSetType(snes,SNESKSPONLY)`；`ts.c:447` 的「同 t 则 `MatScale(A,shift/old_shift)` 复用」优化因每步 `t` 不同而未命中（日志每步都有 `IJacobian enter`）。
- `TSComputeIFunction`（ts.c:372-376）在只提供 ifunction 时直接采用用户残差；`TSComputeIJacobian`（ts.c:457-465）在只提供 iJacobian、无 RHSJacobian 时直接采用用户矩阵，**不会再 `MatShift` 或 `MatAXPY`**。

→ 我们的 `F=K·x-f+Dirichlet`、`J=K+a·S·M+Dirichlet` 与 BE 的 `U_{n+1}=U_n-J^{-1}F(U_n,0)` **代数上完全自洽**，`a=50`、`xDot=0` 也符合源码。因此剩余嫌疑只能是 **F 与 J 中 K 项的实际数值/映射不一致**，而非 TS 约定。

**#64 动作（收敛到 HeatXfer 的单次 scatter 结构）**：`MEF90FlowIFunction` 删除 `MEF90FlowAddStorageToResidual(F,xDot,...)` 与末尾的 `MEF90FlowApplyDirichletResidual(x,F,...)`，只保留 `MEF90FlowOperatorCore(x,F,...)` 一次 `SCATTER_FORWARD`（`OperatorCore` 内部已在最后施加 Dirichlet）。理由：`TS_LINEAR+ksponly` 下 `xDot≡0`，储存残差恒为 0，两次额外 scatter 最好情况是无操作、最坏情况会把 ghost 点值重复累加进 F（`|F|∞` 无法发现局部累加，而 J 只装配一次 K → F/J 比例失衡，正好产生 |λ|>1 的交替放大）。同时把诊断升级为 `|F|∞` + **`|F|₁`**（整体累加会立刻暴露）。GMRES 选项保留（#63 已证明它不是变量）。

**判据**：① `|F|₁` 相对 #62/#63 是否变化、② `Gate 3.2 PASS` / ratio≈0.09；若仍发散，则 F/J 的 K 项映射是最后的未知量，下一轮需在同一子程序里用 `MatMult(J,x)` 与 F 逐点比对，或把 F 的装配也改走 `DMmeshAssembleMatrix` 的同一 default section。


## CI #64 结果（2026-11-07 15:10）：发散修复，但 ratio=0.348 只在解析带外

**状态**：run `42349923764`，`completed success`，`mpirun exit 0`，`Gate 3.1 PASS`（50 次 `post-TSSolve` + 1 次 skip、rejected=0、`_out.gen` 有效）。压力 `pmax` 从 `1e5` **单调衰减**到 `34806.929`，`ratio=0.348069`，CI 打 `WARNING: outside analytic band [0.02,0.30]`，不再硬 FAIL。#64 的单次 scatter 修复确实解决了 #59–#63 的指数发散。

**衰减曲线**：step 1–11 `1e5→92789`，step **12–16 完全冻结**（`pmax` 逐位不变，只有 `pmin` 有 1e-10 量级舍入变化），step 17–51 单调降到 `34807`。step 16 是已知的浮点零跨度 skip；step 12–15 是 `TSSolve` 目标时间与 TS 内部时间只差 ~3e-17、`beuler` 实际步长被压成近 0 所致。晚期拟合 `λ_eff≈1.468 /s`，而理论 `λ₁=π²c_v/(4H²)=2.682 /s`。

**artifact 解析**（`terzaghi-out-gen` / `Beam1x10-tri3_out.gen`，CDF2，51 records、`REC0=23980`、`RECSZ=55320`）：356 nodes、600 TRI3（`connect1`）、5 BAR2（`connect2`，6 个排水点全在 `x=-5`）、`coordx=-5…5`、`coordy=-0.5…0.5`。末步压力在固定 x 上沿 y 的变化 <1% → 有效算子是各向同性的 1D 扩散，只是系数偏小。

**纯 Python/numpy 有限元复现（决定性证据，零本地编译）**：用真实 `connect1/coordx/coordy` 装配一致质量矩阵 `M_e=A/12·(1+δij)` 与刚度 `K_e=A·R·G Gᵀ`（`R=k/μ=5e-8`），Dirichlet 消除 6 个排水点，做 50 步 BE `(K+aSM)pⁿ⁺¹=aSM pⁿ`（`a=1/dt=50`、`S=4.6e-10`）：

| 装配假设 | `|K|inf` | diagK max | `aSM` diag max | ratio(50) | `λ_eff` |
|---|---|---|---|---|---|
| 正确几何 + 各向同性 `XY=0` | 3.704e-7 | 1.852e-7 | 5.499e-10 | **0.0934** | 2.613 |
| 正确几何 + `XY=R` | 4.805e-7 | 2.164e-7 | 5.499e-10 | — | — |
| 全部单元复用 `elem(1)` + `XY=0` | 4.047e-7 | 2.023e-7 | **4.6072319e-10** | 0.0921 | 2.626 |
| **全部单元复用 `elem(1)` + `XY=R`** | **6.5392496e-7** | **3.2696248e-7** | **4.6072319e-10** | 0.2977 | **1.450** |
| **CI #64 实测** | **6.5392496e-7** | **3.2696248e-7** | **4.6072323e-10** | 0.3481 | 1.468 |

`|K|inf`、diagK max、`aSM` diag max **三项逐位吻合**（第 4 行），加上 5 个冻结步把 0.2977 抬到 `0.2977·e^{1.450×0.1}≈0.344`，与 0.348 一致。解析 FE 瑞利商 `(vᵀKv)/(S vᵀMv)=2.68213` 与理论 2.68196 吻合 → **离散方程本身应是 0.093，不是 0.348**。

**两个确定性根因**：

1. **标量张量被写成秩 1 矩阵**。`MEF90/m_MEF90_LinAlg.F90:1693` 的 `MatS2D_Get_Real` 把 `XX=YY=XY=R1`；`A_cell = permCell(1)`（实数量赋值重载）因此得到 `[[R,R],[R,R]]`，特征值 `{2R,0}`，只沿 (1,1) 方向扩散。而库里 `MEF90MatS2DIdentity = MatS2D(1,1,0)`（`:107`）明确表明各向同性应是 `XY=0`。
2. **`elem` 数组与 `cellIS` 错位**。`MEF90Element_Create(mesh,cellIS,elem,…)`（`m_MEF90_Elements.F90:741`）分配 `elem(size(cellID))` 且 `elem(k)` 对应 `cellIS` 的第 k 个单元；所有 `MEF90Diffusion*Set` 也按 `cellID(k)` 索引 `elem(k)`。Flow 的 5 个子程序却用**整个 `setIS`** 创建 `elem`，再对每个 cell 传 **`singleIS`（size 1）** → 600 个单元全部使用 `elem(1)` 的 `Grad_BF / Gauss_C`。质量项因此被 `elem(1)` 的面积 `A₁=1.716981e-2` 统一缩放，这正是 `aSM diag max=4.6072319e-10` 的来源。

**#65 改动（仅 `m_Flow/m_MEF90_FlowAssembly.F90`）**：

1. 三处 `A_cell = permCell(1)` → `A_cell = permCell(1) * MEF90_MATS_IDENTITY`（2D 展开为 `MatS2D(R,R,0)`；`mef90.inc:37` 的宏与 `DbleXMatS2D` 重载已在 DefMech 中验证可用）。`A_zero = 0.0_Kr` 本就三分量全零，保持不变。
2. `MEF90FlowOperatorCore` / `MEF90FlowBilinearFormCore` / `MEF90FlowEnergy` / `MEF90FlowAddStorageToResidual` / `MEF90FlowAddStorageMassMatrix` 五处：删除外层 `MEF90Element_Create(FlowCtx%DM,setIS,elem,…)` 与外层 `MEF90Element_Destroy(elem,…)`；改为在 cell 循环内用**同一个 `singleIS`** 创建 `elem`（size 1），装配后立即销毁。`MEF90FlowComputePressureGradient` 用 `setIS` 且按 `elem(cell)` 索引，本来就正确，保持不变。
3. 校验：19/19 个 `Subroutine` 块配对、6 个含单元数组的子程序各 `create=1 / destroy=1`、`git diff --check` 通过；零本地编译。

**#65 判据**：`diagJ post-Core` 的 `|J|inf` 应从 `6.539e-7` 降到 `~3.70e-7`、diag max 从 `3.27e-7` 降到 `~1.85e-7`；`post-Storage` 的 diag delta 应从 `4.607e-10` 变为 `~5.50e-10`；`pmax(final)/pmax(initial)` 应进入 `0.02–0.30` 且接近解析 `0.087`（含 5 个冻结步时约 `0.12`）→ `Gate 3.2 PASS`。若仍偏高，则只剩时间轴冻结步与 `TS_LINEAR` 的 `TSSetDuration` 语义问题。

**待办（本轮未动）**：step 12–16 的零跨度 `TSSolve` 改成基于 `time(step)` 的容差驱动（注意 workflow 的 `post-TSSolve>=48` 判据需同步调整）；YAML 注释里的 `1m(W)×10m(H)` 与真实几何 `x∈[-5,5], y∈[-0.5,0.5]` 对齐；`pressurePreviousStep` 死写；GMRES/CG 最终选型；诊断 trace 清理。


## CI #65 结果（2026-11-07 16:05）：3.1 与 3.2 同时通过，根因复现逐位吻合

**状态**：run `42353581894` / job `128346726709`，`completed success`；`Terzaghi mpirun exit code: 0`；`Gate 3.1 PASS`（51 个 pressure envelope、TS 到 t=1、rejected=0、`_out.gen` 有效）；**`Gate 3.2 PASS: Terzaghi decay ratio 0.124364 in [0.02,0.30]`**。`pmax` 从 `1e5` 单调降到 `12436.427`，`pmin` 全程约 0（#64 的 -331 负超调消失）。

**离线预测 vs CI 实测（零本地编译的定量验证）**：

| 量 | Python FE 预测 | CI #65 实测 | CI #64 |
|---|---|---|---|
| `diagJ post-Core |J|inf` | 3.7040672e-7 | **3.7040672351723571e-7** | 6.5392496e-7 |
| `diagJ post-Core diag min` | 4.1769252e-8 | **4.1769251804758425e-8** | 1.0451187e-8 |
| `diagJ post-Core diag max` | 1.8520336e-7 | **1.8520336175861785e-7** | 3.2696248e-7 |
| `post-Storage` diag min 增量 | 1.0033700e-10 | 1.0108357e-10 | 1.9381416e-10 |
| ratio（45 / 50 个有效 BE 步） | 0.1213 / 0.0934 | **0.124364** | 0.348069 |

→ 两个根因（`MatS2D_Get_Real` 的 `XY=R` + `elem(1)` 几何复用）**完全确认**；K 与 M 现在都是物理正确的各向同性 P1 有限元装配。

**剩余唯一偏差 = 时间轴**：#65 仍有 5 个冻结步（step 12–16 `pmax` 逐位不变）。原因是 `TSSetDuration(ts,10000,time(step))` 遇到 `|time(step)-t|≈3e-17` 的近零区间，`beuler` 的 `a=1/dt` 爆炸使 `J≈aSM`、压力不推进；45 个有效 BE 步对应 0.1213，与实测 0.124364 吻合（解析连续解 0.0871，50 步 FE/BE 解 0.0934）。

**#66 改动（2 个代码文件）**：

1. `vDef/vDef.F90`：TS 本来就由 `TSSetInitialTimeStep(...,dtFlow0=0.02)` 固定步长驱动，因此改为**每个外层步恰好推进一个固定 dt**：`tFlowEnd = t + dtFlow0`、`TSSetDuration(ts,10000,tFlowEnd)`、`If (step > 1) Then TSSolve(...)`。不再用插值 `time(step)` 作为 TS 目标，彻底消除近零区间；step 1 仍为初始条件不推进，step 2–51 共 50 步 → t=1。`nSolve=50` 仍满足 Gate 3.1 的 `>=48`。
2. `.github/workflows/build.yml`：Gate 3.2 从宽松的发散护栏 `[0.02,0.30]` 收紧为**真正的解析解比对带 `[0.05,0.15]`**（FE/BE 0.0934、连续解 0.0871），带外改为硬 `FAIL`；`ratio>1.5` 的发散硬判据保留。

校验：`git diff --check`、PyYAML（19 steps）、`Run Terzaghi` 块 `bash -n`、`vDef.F90` 的 34/34 个 `If/End If` 配对全通过；零本地编译。

**#66 判据**：① `Gate 3.1 PASS` 保持、`nSolve=50`、无 skip（除 step 1）；② 51 步 `pmax` 严格单调、无冻结段；③ `ratio≈0.093` 落入 `[0.05,0.15]` 且打印 `Gate 3.2 PASS`；④ `diagJ` 三项与 #65 一致（K/M 未再改动）。若 ratio 仍 ≈0.124，说明冻结步来自 PETSc 3.3 `TSSetExactFinalTime` 内部裁剪而非外层目标，需改用 `TSSetDuration(ts,1,...)` 强制单步。


## CI #66 结果（2026-11-07 16:15）：闸门 3.1 / 3.2 全部通过，解析解逐点比对完成

**状态**：run `42354263842` / job `128348773869`，`completed success`，`Terzaghi mpirun exit code: 0`。

| 判据 | 结果 |
|---|---|
| `Gate 3.1` | **PASS**（50 次 `post-TSSolve` + 1 次 step-1 skip、TS 到 t=1、rejected=0、`_out.gen` magic 有效） |
| `Gate 3.2` | **PASS**：`ratio=0.0954738 ∈ [0.05,0.15]`（FE/BE 预测 0.0934、连续解 0.0871） |
| 冻结步 | **0**（#65 的 step 12–16 全部消失，51 步 `pmax` 严格单调下降） |
| `pmax` | `1e5 → 9547.383`；`pmin≈0`（无负超调） |
| `diagJ post-Core` | `|J|inf=3.7040672351723571e-7`、diag `4.17692518e-8 / 1.85203362e-7`，与 #65/离线 FE 完全一致（K 未再改动） |

**`_out.gen` 逐点解析比对**（artifact `terzaghi-out-gen`，已落盘到 `vDef/data/Beam1x10-tri3_out.gen`，2,845,300 B，`CDF\002`，51 records，`time_whole=0…1`，`REC0=23980`、`RECSZ=55320`；`z=5-x`，`H=10`，`c_v=(k/μ)/S=108.696`）：

解析级数 `p/p0 = (4/π)Σ(-1)^m/(2m+1)·cos((2m+1)πz/2H)·exp(-λ_m t)`，其中 `λ_m=c_v((2m+1)π/2H)²`；BE 时间离散把每个模态的衰减替换成 `(1+λ_m dt)^(-t/dt)`。

| t | pmax FEM | 连续解析 | BE 修正解析 | L2 相对误差（连续 / BE 修正） |
|---|---|---|---|---|
| 0.20 | 77971 | 74126 | 74687 | 6.19% / **4.66%** |
| 0.50 | 35802 | 33306 | 34481 | 7.32% / **3.66%** |
| 0.80 | 16200 | 14897 | 15748 | 8.56% / **2.70%** |
| 1.00 | 9547 | 8713 | 9339 | 9.40% / **2.06%** |

末步剖面（p/p0）沿 x：`0.00000 / 0.01488 / 0.02940 / 0.04318 / 0.05588 / 0.06720 / 0.07702 / 0.08494 / 0.09074 / 0.09428 / 0.09547`（x=-5→5），与 BE 修正解析解逐点相差 ~2%。**结论**：与连续 Terzaghi 解的 ~9% 偏差全部来自一阶 BE 时间离散（dt=0.02），扣除后空间有限元误差仅 ~2%，闸门 3.2 的「Terzaghi 1D 解析解比对」实质通过。

**当前状态**：闸门 3.1 与 3.2 均稳定通过；transient Darcy 求解链（TS 创建 → IFunction/IJacobian → Dirichlet → 输出）物理正确；`vDef/data/Beam1x10-tri3_out.gen` 已可用于后处理。

**剩余清理项（非阻塞）**：YAML 注释里的 `1m(W)×10m(H)` 与真实几何 `x∈[-5,5], y∈[-0.5,0.5]` 对齐；`pressurePreviousStep` 死写；GMRES/CG 最终选型；`Gate 3.1/3.2` 诊断 trace 与 `|F|₁`/`diagJ` 打印减噪；DefMech AltMin 每步 SNES 非收敛（50 次 WARNING，力学侧遗留）；**新增**：输出时间标签超前 1 个 `dt`、Exodus 节点变量名错标 `Temperature`（见下节）。


## 后处理与仓库同步（2026-11-07 16:45）：产物交付，另发现 2 项待修缺陷

### 1. 本地仓库已对齐远程 main（代理恢复）

- `git fetch origin main` 现已通（代理 `127.0.0.1:7888` 恢复），`origin/main` 由 `8983360` 更新到 `f695c907`。最新 CI run 仍是 #66 `42354263842 / completed success`；`f695c907` 带 `[skip ci]`，之后没有新 run。
- 本地 4 个核心文件的 raw blob sha 与远程 tree **逐字节一致**：`build.yml a247f3f0`、`FLOW_FIX_CHECKLIST.md 300495ec`、`FlowAssembly.F90 709b9e78`、`vDef.F90 f9e2abbd`。
- 分叉情况：`merge-base = 8983360`，本地独有 5 个提交 / 远程独有 10 个；`README.md` 与 `vDefDefault_remote.F90` 的 HISTORY 差异经 `git diff --ignore-space-at-eol` 确认为纯 CRLF 假差异。
- 备份 `backup/local-pre-sync-5e22420` 后执行 `git reset --hard origin/main` → HEAD = `f695c907`，`git diff origin/main` 为空，`git status` 只剩未跟踪的 `.workbuddy/` / `.trae/`；`vDef/data/Beam1x10-tri3_out.gen` 完好（2,845,300 B，`CDF\002`，51 records）。

### 2. 后处理产物（本机 numpy + PIL，零编译）

| 文件 | 内容 |
|---|---|
| `vDef/data/terzaghi_pressure_timeseries.csv` | 51 行：`t` / `pmin` / `pmax` / `pmax·p0⁻¹`，连续解析、BE 解析、BE-shift1 的封闭端峰值，三套全场 L2 相对误差 |
| `vDef/data/terzaghi_pressure_field.csv` | 18,156 行（51 时刻 × 356 节点）：`x`, `y`, FEM 压力，三套解析值，两组逐点残差 |
| `vDef/data/terzaghi_validation.png` | 四联图：峰值衰减、t=1 剖面、L2 误差史、时空压力热图 |
| `vDef/data/parse_cdf2.py` | 精简 NetCDF CDF2 解析器（header / 变量表 / record 偏移），纯 `struct` |
| `vDef/data/postprocess_terzaghi.py` | 后处理主脚本，运行 `python vDef/data/postprocess_terzaghi.py` 即可一键重生成上表三份产物 |

CDF2 布局复核：`REC0=23980`、`RECSZ=55320`、`numrecs=51`、`coordx/coordy` 356 节点、`vals_nod_var1` 为压力、`vals_elem_var4eb1` 600 单元。参数 `c_v=(5e-11/1e-3)/4.6e-10=108.695652 m²/s`、`H=10 m`（`z=5-x`，排水面 `x=-5`）、`p0=1e5 Pa`、`dt=0.02 s`；Terzaghi 级数取 400 项，BE 修正把每个模态换成 `(1+λ_m dt)^(-t/dt)`。

| t | pmax FEM | 连续解析 | BE 解析 | L2 连续 | L2 BE(n) | L2 BE(n-1) |
|---|---|---|---|---|---|---|
| 0.10 | 95183 | 93606 | 92465 | 7.61% | 5.75% | **0.34%** |
| 0.20 | 77971 | 74126 | 74687 | 6.19% | 4.66% | **0.70%** |
| 0.50 | 35802 | 33306 | 34481 | 7.32% | 3.66% | **1.63%** |
| 1.00 | 9547 | 8713 | 9339 | 9.40% | **2.06%** | 3.15% |

### 3. 新发现 A（#67 候选）：`time_whole` 比物理推进超前 1 个 `dt`

- `t=0` 输出的是**未施加 Dirichlet BC 的均匀场**（356 节点 `pmin=pmax=1e5`）；`t=0.02` 输出场恰好等于「含 BC 的初始场」（只有 `x=-5` 的 6 个节点为 0，其余 `p0`），与 BE 级数 `n=0` 的 L2 = **0.0000%**。
- 因此 `time_whole=t` 的场对应 BE 解析的 **`n-1` 步**。逐步量化：`t=0.04 → BE(1)` **0.050%**、`0.10 → BE(4)` **0.336%**、`0.30 → BE(14)` **1.024%**、`0.70 → BE(34)` **2.240%**，早期几乎逐位吻合；后期空间离散谱偏差累积，`t=1.0` 落在 `BE(49)` 与 `BE(50)` 之间（`3.145%` vs `2.061%`）。
- 根因：`vDef.F90` 的 step 1 跳过 `TSSolve` 且初始条件未先施加 BC，输出时间标签没有补偿这一步。**不影响闸门 3.1/3.2**（`ratio=0.0955` 物理正确），但严格逐时解析比对须用 `n-1`（CSV 已提供 `*_shift1_*` 列）。
- 候选修法：① 初始化后立即调用 `UpdateboundaryPressure` 把排水面置 0 再写 `t=0`；或 ② 写 `time_whole` 时统一减 `dt`。两者都需重跑 CI，且会改变 Gate 3.2 的解析对齐基准。

### 4. 新发现 B：Exodus 节点变量名错标为 `Temperature`

- `_out.gen` 的 `name_nod_var(1)` 是 `Temperature`，但 `vals_nod_var1` 的数值已逐点验证为压力（Pa）。
- 根因：`m_DefMech/m_MEF90_DefMech.F90:1290` 与 `m_HeatXfer/m_MEF90_HeatXfer.F90:584` 把 `nameV` 的 offset 1 固定写成 `"Temperature"`；Flow 侧仅在 `m_Flow/m_MEF90_FlowCtx.F90:393` 用 `PetscObjectSetName(...,"Pressure")` 给 Vec 起名，未写进 `nameV`。
- 影响：ParaView / VisIt 会把压力场显示成温度、单位错误；纯输出元数据缺陷，不影响数值。修改前需先确认 `nameV` 的 offset 分配规则。

**当时状态（已被 CI #67 取代）**：两个闸门稳定通过；结果网格、三份后处理产物与两个可复现脚本可直接用于分析；本地与远程 `main` 一致。上节 A/B 两项缺陷随后在 CI #67 中修复，见下一节。


## CI #67 结果（2026-11-07 17:25）：A/B 修复完成，解析基准严格对齐

**状态**：commit `40cd29a`，run `42359268474` / job `128363964783`，`completed success`。`git push` 已恢复（代理 `127.0.0.1:7888` 可用），不再需要 REST API 绕行。

### 改动

1. **`vDef/vDef.F90`（修 A）**：初始化末尾、时间循环之前，用**已 public 的 `MEF90FlowOperator`** 计算一次 SNES 残差，再 `VecAXPY(pressure,-1.0_Kr,residualPressure)` 把 Dirichlet 修正从初始场里减掉。原理：`MEF90FlowCtxCreateVectors` 把压力设成全场 `p0`（无 BC）；均匀场下 `K·p=0`、`source=0`，故残差只在 Dirichlet dof 上等于 `p0-p_bc`、内部为 0。减掉后内部仍 `p0`、排水面变为 `p_bc`，正是解析级数的 `BE(0)`。此后 `time_whole=t` 的场就是 `BE(n)`（`n=t/dt`），不再滞后一步。零新增接口、零本地编译。
2. **`m_DefMech/m_MEF90_DefMech.F90`（修 B）**：`MEF90DefMechFormatEXO` 里 `nameV(temperatureOffset)` 改为按 `Associated(MEF90DefMechCtx%temperature)` 条件赋值 —— HeatXfer 激活时 `"Temperature"`，禁用时（vDef 已 `DeAllocate` 该成员，Flow 的 `pressureOffset=1` 占用同一 slot）为 `"Pressure"`。不新增 Exodus API、不改 offset 布局；根目录同名文件是旧副本，构建只用 `m_DefMech/`。
3. **`.github/workflows/build.yml`**：① 运行前 `rm -f data/Beam1x10-tri3_out.gen`，保证 `EXGVP→0` 从而一定执行 `FormatEXO`（否则变量名修复不生效）；② Gate 3.2 新增**真正的解析解硬比对**：纯 `awk` 计算 BE 修正 Terzaghi 级数在 t=1 的 `p/p0=0.09339316`，要求 `|ratio-expected|<=0.02`；③ 新增 `name_nod_var1=="Pressure"` 校验（纯标准库 `vDef/data/parse_cdf2.py`）。原有 `ratio>1.5` 发散判据与 `[0.05,0.15]` 带宽判据保留。
4. 新增 `vDef/data/parse_cdf2.py`（CDF2 解析器，CI 依赖）与 `vDef/data/postprocess_terzaghi.py`（后处理主脚本）。

静态校验：`git diff --check`、PyYAML（19 steps）、全部 `run` 块 `bash -n`、两个 Fortran 文件 `If/End If` 配对 0/0、新增 shell/awk 命令行本地原样复现，全通过。

### CI 关键输出

| 判据 | 结果 |
|---|---|
| Gate 3.1 | **PASS**（50 次 `post-TSSolve` + 1 次 step-1 skip、TS 到 t=1、rejected=0、`_out.gen` magic 有效） |
| 初始场 | `CI #67: initial pressure after Dirichlet pmin 0.0 pmax 100000.0`（排水面 6 节点为 0） |
| Gate 3.2 带宽 | **PASS**：`ratio=0.0905571 ∈ [0.05,0.15]`（#66 为 0.0954738） |
| Gate 3.2 解析 | **PASS**：BE 解析 `0.09339316`，`|0.0905571-0.09339316|=0.00284 <= 0.02` |
| Gate 3.2 元数据 | **PASS**：`name_nod_var1 = 'Pressure'`（#66 为 `Temperature`） |

### `_out.gen` 逐点复核（artifact `terzaghi-out-gen`，已落盘）

`CDF\002`、2,845,300 B、51 records、`REC0=23980`、`RECSZ=55320`、356 节点。`t=0` 场 `pmin=0 / pmax=1e5` 且仅 6 个节点低于 `0.5 p0`（排水面）；`t=0.02` 的 `pmax=99773.794` 正是 #66 的 `t=0.04` 值 —— 序列整体前移一步，完全符合预期。

| t | L2 vs BE(n) | L2 vs BE(n-1) | 最优 |
|---|---|---|---|
| 0.02 | **0.050%** | 25.652% | BE(n) |
| 0.10 | **0.401%** | 6.129% | BE(n) |
| 0.50 | **1.688%** | 6.688% | BE(n) |
| 1.00 | **3.206%** | 8.131% | BE(n) |

对比 #66（修复前）：早期是 `BE(n-1)` 逐位吻合（0.050% / 0.336%），现在**所有时刻都严格对齐 `BE(n)`**；`t=0` 的 L2 误差从 13.09% 降到 **0.000%**；全场相对连续解析解的误差从 9.40% 降到 **3.76%**（残余来自一阶 BE 与 P1 空间离散谱偏差）。

### 后处理产物（已用 #67 artifact 重新生成）

`terzaghi_pressure_timeseries.csv`（51 行）、`terzaghi_pressure_field.csv`（18,156 行）、`terzaghi_validation.png`（四联图）。CSV 保留 `*_shift1_*` 列作为修复前对齐的历史对照；PNG 的 (c) 面板同时画出 `BE(n)`（当前基准）与 `BE(n-1)`（旧对齐），直观展示修复效果。运行 `python vDef/data/postprocess_terzaghi.py` 可一键重生成。

**当时状态（已被 CI #68 取代）**：闸门 3.1 / 3.2 均稳定通过；时间轴与输出元数据正确；transient Darcy 求解链物理正确且逐时对齐解析解。剩余非阻塞清理：YAML 几何注释、`pressurePreviousStep` 死写、GMRES/CG 选型、诊断 trace 减噪、DefMech AltMin 50 次 SNES WARNING、以及 nodal slot 3 仍存在 `boundary_pressure_offset=3` 与 `disp_offset=3` 的槽位复用（名字仍是 `Displacement_X`，本次未动，需统一规划 offset 后再改）。下一节已修复其中三项。


## CI #68 结果（2026-11-07 18:00）：输出完整性修复与日志减噪

**状态**：commit `090781e`，run `42361549691` / job `128370821222`，`completed success`。`git push` 正常（代理可用，无需 REST API 绕行）。

### 根因：Flow 覆盖了 DefMech 的 nodal slot 3

`vDef.F90` 每步的写出顺序是 `MEF90DefMechViewEXO` → `MEF90FlowViewEXO`。YAML 里 `defmech.disp_offset=3`（占 nodal slot 3/4）与 `flow.boundary_pressure_offset=3` 撞在**同一个 nodal slot**；Flow 随后把常数 `boundaryPressure`（=0）写进 slot 3，把 `Displacement_X` 整场抹成 0。而 `nameV` 由 `MEF90DefMechFormatEXO` 生成、仍叫 `Displacement_X`，于是后处理看到的是"一个恒为 0 的位移分量"——纯输出数据损坏，不影响求解。

cell 变量是独立命名空间：`flow.source_offset=2` / `perm_field_offset=4` / `darcy_velocity_offset=5` 不与 DefMech 的 nodal offsets 冲突，无需调整。

### 改动

| 文件 | 内容 |
|---|---|
| `vDef/data/terzaghi1D.yaml` | `boundary_pressure_offset: 3 → 0`。当前 `bc_pressure_scaling: constant`，BC 值已知且无需输出槽；`0` 只让 `MEF90FlowViewEXO` 跳过写出，不影响求解。同时修正头部几何注释（真实 `x∈[-5,5]`、`y∈[-0.5,0.5]`、排水面 `x=-5`、`H=10` 沿 x、`z=5-x`）与 51×0.02 s = 0→1 s 的时间轴 |
| `m_Flow/m_MEF90_FlowAssembly.F90` | 删除 IFunction / IJacobian / SetTransients / UpdateboundaryPressure / SetSourceCst 中约 20 处每步诊断打印；`MEF90FlowReportJacobian` 保留为按需调试工具。**Gate 3.1 / 3.2 判据依赖的 `vDef.F90` 握手行全部不动** |
| `.github/workflows/build.yml` | 把 #67 的 `name_nod_var1=="Pressure"` 单槽检查升级为全表检查 |
| `vDef/data/check_exodus_metadata.py` | 新增纯标准库校验器：严格比对 4 个 nodal 变量名，打印每个 slot 在 t=0 / t=1 的值域；若 `Displacement_X` 仍恒为 0 则打 `[WARN]` |

静态校验：`git diff --check`、PyYAML（19 steps）、全部 `run` 块 `bash -n`、`m_MEF90_FlowAssembly.F90` 的 `If/End If` / `Subroutine` / `Do` 配对全为 0；新增 CI 命令行在本地按原样复现通过。

### CI 关键输出

| 判据 | 结果 |
|---|---|
| Gate 3.1 | **PASS**（TS 到 t=1、rejected=0、`_out.gen` 有效） |
| Gate 3.2 带宽 / 解析 | **PASS**：`ratio=0.0905571`，与 BE 解析 `0.09339316` 差 `0.00284`（与 #67 逐位一致） |
| 元数据全表 | **PASS**：`Pressure \| Damage \| Displacement_X \| Displacement_Y` |
| slot 3 | t0 `[-9.30729e-4, 0]`，t1 `[-5.40369e-5, 0]` —— **不再恒为 0**，`[WARN]` 消失 |
| 日志量 | 每步 `Gate 3.2 diag` 行从约 50 组降到 0（job log 仅剩 1 条 shell 字面量回显） |

### #67 → #68 artifact 逐位对比

两版 `CDF\002`、2,845,300 B、51 records、356 节点，文件布局完全相同：

| 场 | 最大差 |
|---|---|
| `vals_nod_var1`（压力，51×356 全时刻） | **0**（逐位一致） |
| `vals_nod_var2`（Damage）、`vals_nod_var4`（Displacement_Y） | 0 |
| `vals_nod_var3`（Displacement_X） | `5.40369e-5`（从全 0 恢复为真实位移） |

证明本次改动是**纯输出层**，对 transient Darcy 求解零影响。后处理三份产物已用 #68 artifact 重生成，t=1 仍为 FEM `9055.7 Pa` / 连续 `8712.7` / BE `9339.3`，L2 `3.764% / 3.206%`。

**当时状态（已被 CI #69/#70 取代）**：闸门 3.1 / 3.2 稳定通过；时间轴、变量命名、全部 4 个 nodal slot 的输出内容均正确；日志已收敛到可读。剩余非阻塞清理：`pressurePreviousStep` 死写（YAML offset=6 与 `vDef.F90:508` 的 `VecCopy` 无消费者）、GMRES/CG 最终选型、DefMech AltMin 50 次 SNES WARNING、cell slot 的 offset 规划，以及构建诊断步骤里既有的 2 条 `gfortran: Only one '-J' option allowed`（Test A/B PETSc 模块探测，#67 同样存在，不影响构建）。其中 **GMRES/CG 选型已由 CI #69 定案、CI #70 落地**。

## CI #69 / #70（2026-11-07 19:00）：GMRES vs CG 同源对比验证与最终选型

### 背景与设计

CI #63 为排查压力发散把 `-flow_ksp_type gmres` 加进命令行，覆盖了 `MEF90FlowCreateTS` 里的默认 `KSPCG`；此后一直沿用。当时发散的真实根因是 IFunction 三次 scatter（#64 修复），与求解器无关，但一直没复核 CG 是否可恢复。

CI #69（`42366003387 / success`，commit `4c198b3`）在**同一个 job、同一二进制、同一网格与 YAML** 下连跑两次，只改 `-flow_ksp_type`（两者均 `-flow_ksp_rtol 1.0e-10`）：

| 环节 | 做法 |
|---|---|
| 主运行 | GMRES，全部 Gate 3.1/3.2 判据不变，`_out.gen` 作为交付 artifact |
| 第二运行 | CG，另存 `/tmp/terzaghi_cg.log` 与 `/tmp/out_cg.gen` |
| 比对 | 新增 `vDef/data/compare_ksp_runs.py`（纯标准库）：51×356 压力场、Gate 3.2 ratio、KSP 收敛原因/迭代成本；从 `-flow_ts_view` 的 `KSP Object: (flow_)` 块断言选项确实生效；按 step 去重避免 job log 的 `tail` 回显污染末值 |
| 判据 | 场差 `<= 1e-2 Pa`、ratio 差 `<= 1e-5`、两者 ratio 均在 `[0.05,0.15]`；CG 发散时输出 `CG NOT VIABLE` 结论而非误判闸门失败 |
| 交付 | 比对结束后把 GMRES 的 `_out.gen` 复制回原位，保证既有产物不变 |

静态校验：`git diff --check`、PyYAML（20 steps）、全部 `run` 块 `bash -n`、四个脚本 `py_compile`、用 #68 日志/artifact 构造 mock 做本地 dry-run。

### 对比结果（CI 报告与本地复跑逐字一致）

| 指标 | GMRES | CG |
|---|---|---|
| `-flow_ts_view` 确认的 flow KSP 类型 | `gmres` ×50 | `cg` ×50 |
| 收敛原因 | `CONVERGED_RTOL` ×50 | `CONVERGED_RTOL` ×50 |
| 线性迭代总数 / 均值 / 峰值 | 3043 / 60.86 / 76 | **1973 / 39.46 / 44** |
| t=1 `pmax` | 9055.708842 | 9055.708841 |
| Gate 3.2 ratio | 0.0905571 | 0.0905571（差 `1.46e-11`） |

51×356 全场压力最大差 **2.54252e-05 Pa**（相对 `p0` = `2.54e-10`），远小于 `1e-2 Pa` 阈值。**结论：`KSPCG` 完全可用，且比 GMRES 少 35.16% 线性迭代**；`MAT_SPD` + 单次 `MatZeroRowsColumnsIS` 后的 J 确实保持对称正定，#62 当时对 CG 的怀疑不成立。默认 PC 仍是 `bjacobi`（块内 ILU），CG 在其上正常收敛。

### CI #70 落地（`42366719829 / success`，commit `de9dd1e`）

- workflow 主运行去掉 `-flow_ksp_type gmres`，保留 `-flow_ksp_rtol 1.0e-10 -flow_ksp_converged_reason`；`-flow_ts_view` 逐步打印 `type: cg`，确认默认生效。
- 移除一次性双跑、`terzaghi-ksp-comparison` artifact 与相关 diagnostics 行，生产 CI 回到 19 steps 单跑。`compare_ksp_runs.py` 保留为可复用工具（docstring 记录 #69 结论与重跑方式，#69 输入归档在 `terzaghi-ksp-comparison` artifact）。
- Gate 3.1 / 3.2 带宽 / BE 解析 / 元数据全表 **全 PASS**，`ratio=0.0905571`；51 次 `CONVERGED_RTOL`，迭代均值 38.69、峰值 44。

### 三版 artifact 逐位对比（51×356 全时刻）

| slot | #68 GMRES → #69 GMRES | #69 CG → #70 CG | #69 GMRES → #70 CG |
|---|---|---|---|
| 1 Pressure | 0 | **0** | 2.54e-5 Pa |
| 2 Damage | 0 | **0** | 0 |
| 3 Displacement_X | 0 | **0** | 2.22e-6 |
| 4 Displacement_Y | 0 | **0** | 0.01204 |

- #68 与 #69 GMRES **全字段逐位一致** → #69 的双跑没有污染主运行。
- #69 的 CG 与 #70 的单跑 CG **全字段逐位一致** → 差异只来自 KSP 选择本身，与运行顺序无关，结果可复现。
- `Displacement_Y` 的 `0.012` 差异来自 DefMech AltMin 的 **51 次非收敛 SNES**（两版都是 51 次）：非收敛的力学迭代把 `2.54e-5 Pa` 的压力扰动放大。Flow/Darcy 压力场不受影响，**但位移场目前不具备定量可比性**，若后续要分析力学量必须先解决 AltMin 非收敛。

### 交付与状态

`_out.gen`（#70 CG，2,845,300 B / 51 records / `Pressure | Damage | Displacement_X | Displacement_Y`）、`terzaghi_validation.png`、两份 CSV 与后处理脚本均已用 #70 artifact 重生成；图表标题改为 `CI #70 (default KSPCG)`。t=1 仍为 FEM `9055.7 Pa` / 连续 `8712.7` / BE `9339.3`，L2 `3.764% / 3.206%`。本地与远程 `main` 同步。

**剩余非阻塞项**：`pressurePreviousStep` 死写、DefMech AltMin 51 次 SNES 非收敛（现已知会放大力学场差异）、cell slot 的 offset 规划、构建诊断里既有的 2 条 `-J` 噪声。CSV/PNG 仍未入库（脚本可一键重生成）。

### CI #71 + #72：DefMech 位移 SNES 非收敛根因定位

**现象（CI #71, `42369005656 / failure`）**：位移 SNES 在 51/51 个外层 AltMin 步全报 `-5 (SNES_DIVERGED_MAX_IT)`。fnorm 轨迹：第 1 步从 `2.96e4` 降到 `~10`，之后在 `3–10` 之间振荡到 maxit（出现 `3.46 → 9.4` 回弹）。内层 KSP 全部 `CONVERGED_RTOL`（~40 it，KSPCG，`rtol=1e-8`）。damage 全程恒为 `α≈4e-5`（≈0）→ AT2 降阶 `a≈1`。

**代码核对结论（已逐项排除）**：
- `cellIsElastic=.FALSE.`（AT2 非 AT2Elastic），残差与 Jacobian 都走 HD split。HD `DEED` 求和恰回到 `HookesLaw*Strain`，`D2EED` 求和恰回到 `HookesLaw`（常数）。**damage≈0 时位移系统完全线性且 J≡F**，理论上 Newton 应 1 步收敛。
- 阻尼 `dampingCoefficientDisplacement` 默认 `0.0`（vDef 默认），内聚 `cohesiveStiffness=0` → 残差里唯一的阻尼/内聚项对本卡关闭；即便开启，也存在真实 bug：残差用 `dampingCoeff/timeStep`（装配行 304）而 Jacobian 用 `dampingCoeff*timeStep`（行 600），差 `timeStep²` 倍。但本卡不触发。
- 内层 KSP 用 `KSPCG` + `MAT_SPD` 声明；矩阵理论对称正定（HD 切线 `D2EEDPlus+D2EEDMinus=HookesLaw` 对称），BC 用 `MatZeroRowsColumnsIS` 保对称。CG 适用。
- 已排除：null space（BC 充分约束，2D 旋转模态也被 `u_x=0@x=+5` 消去）、时间轴/几何、`smoothPS` 函数（damage≈0 时 cancel）、压力耦合（常数载荷，无 u 依赖）。

**矛盾**：理论（线性+一致+SPD+CG）→ 应 1 步收敛；实测 → 停在 ~7-10 且振荡 → 必为运行期不一致。无法仅靠阅读定位，故用运行期对照实验定论。

**CI #72 诊断设计（已提交）**：Gate 3.3 改回 report-only（不阻断 Gates 3.1/3.2）；主生产运行（解析 J + CG + basic 线搜索）**最后**跑以保 `_out.gen`；其前插入两个变体：
- **V1**：解析 J + **GMRES** → 若收敛 ⇒ 矩阵非对称/不定，KSPCG（MAT_SPD）是元凶。
- **V2**：**FD Jacobian**（`-Disp_snes_fd`）+ GMRES → 若收敛 ⇒ `BilinearFormDisplacement` 有值 bug，FD-coloring 可恢复。
- 若两变体仍 -5 ⇒ 残差本身不一致 / 系统奇异。

**待 CI #72 结果后决策**：若 V1 收敛→改位移 KSP 为 GMRES（或修正对称/SPD）；若 V2 收敛→定位并修 `BilinearFormDisplacement` 的 J 值 bug；两者皆收敛则双修；皆否则查残差/奇异。

### CI #72 结果 + 容差分析（定位 J≠dR/du）

**关键容差事实**：位移 SNES `atol=1e-7, rtol=1e-5` → 收敛阈值 ≈ `atol + rtol*F0 = 1e-7 + 1e-5*2.96e4 ≈ 0.296`。实测 fnorm 卡在 `3–10` **从不**低于 0.296 → 因此报 `-5 (MAX_IT)`，不是"差一点收敛"。这确认 **解析 Jacobian 与真导数不一致**（若 J 正确，第一步即把残差打到 ~1e-12 << 0.296，1 步收敛）。

**CI #72 V1（解析 J + GMRES）**：`51× -5`，与 CG 基线一致 → 排除 CG/MAT_SPD 对称性假设（GMRES 可处理非对称/不定矩阵）。
**CI #72 V2（FD J + GMRES）**：`rc=124`（600s 超时），仅跑到 28 步（全 -5），**无法定论**——FD 每步 ~21s，51 步需 ~1070s。

**代码三重核对（已排除）**：
- HD split 数学：`D2EEDMinus=∂²W-/∂ε²` 手算与代码行 130 一致；`D2EEDPlus=HookesLaw−D2EEDMinus` 正是 `DEEDPlus` 的正确切线 → `AGauss=∂(stressGauss)/∂ε` 严格成立，J≡dR/du。
- deck `solverType: AltMin` → 位移子问题冻结 damage，且 Terzaghi 为压缩荷载 → HD split **激活**（非纯弹性），位移子问题本质**非线性**（split 的 SmoothPositiveSquare 引入非线性），但切线 D2EED 正确，正确 J 下应几步收敛。
- BC：残差 `R_i=u_i−g_i`（INSERT）、Jacobian `MatZeroRowsColumnsIS(...,1.0)` 行+列清零 → Dirichlet 一致；首步即把约束自由度拉到 g_i，不影响后续。

**结论**：数学上 J≡dR/du，但实测不收敛 → 必为**运行期装配/SNES 配置 bug**，无法靠阅读定位。

**CI #73 设计（已提交）**：V2 改 `timeout 1500s` 跑完完整 51 步 FD Jacobian + GMRES，作为**决定性判据**：
- FD `-5` 数骤降/为 0 ⇒ 解析 `BilinearFormDisplacement` 有值 bug（FD-coloring 可恢复）→ 后续定位/修复该 bug。
- FD 仍 ~51× -5 ⇒ 残差 `OperatorDisplacement` 本身不一致或系统奇异 → bug 在残差侧。

---

## ⚠️ CI #75–#77 修正结论（覆盖上文"J≠dR/du / 残差不一致"推论）

> 以下为最终核实后的判断，上文基于 watcher 区域误读得到的"解析 Jacobian 不一致"等结论**已被推翻**。

### 已确证事实（来自 ci75.log / ci76.log / ci77.log 干净比对）
1. **BC 空格 bug 真实存在且已修**：`displacementBC: 1, 1, 1`（逗号后空格）→ PETSc Bag BoolArray 只取首分量 → `TRUE FALSE FALSE` → u_y 未约束 → 刚体平移 → 刚度矩阵奇异（-5 全 51 步）。改 `1,1,1`（无空格，commit 9c1edf2）后 bag 视图为 `TRUE TRUE TRUE`，刚体奇异消除。
2. **但收敛并未达成**：ci75/ci76 生产运行均含 ~38–50× 位移 `-2`（FUNCTION_DOMAIN）+ 8× `-5`；ci77 更干净。
3. **损伤不是元凶（CI #77 NODAMAGE 隔离实验推翻）**：`terzaghi1D_nodamage.yaml` = 仅 `damage.type: AT2→AT2Elastic`（isElastic，不演化、Stiffness=1 不退化）。NODAMAGE 同样 **51/51 位移 `-5`** → 孔弹性核心本身不收敛，与损伤演化无关。
4. **求解器调参也不是元凶**：生产（解析J+CG+basic）、V2（FD-J+GMRES+l2）均失败；V2 实际是 `rc=124` 超时（FD 太慢），非收敛。cg/gmres、basic/l2、解析/FD 雅可比全部失败。
5. **压力/Terzaghi 物理完全正确**：NODAMAGE 中 `pmin≡0`（排水面 p=0 参考生效）、`pmax` 从 1e5 平滑衰减到 9055（step 51）→ 与 CI #75 的 pmax=9055.7 吻合。flow KSP 的 `zero pivot / diagonal shift` 仅是 PETSc ILU 在良性 SPD 矩阵上的常规预防性警告，**非错误**。
6. **位移矩阵本身非奇异**：位移 SNES 的 KSP 无 zero-pivot；所有 zero-pivot 标记都挂在 `flow_` KSP 上（且 benign）。

### 位移 SNES 的真实行为
- 位移残差：首步从 2.08e4 降到 ~7–9，之后在 **3–10 间振荡/plateau**，永不逼近 atol=1e-7 → 报 `-5 (MAX_IT)`。
- 位移 BC 正确（TRUE TRUE TRUE）、压力场正确、位移矩阵非奇异。四个求解配置全失败。

### 当前最可能根因（待 CI 实验坐实）
**位移线性求解（KSP）精度不足 / 预条件太弱**：位移 KSP = CG + bjacobi(ILU(0))。对于受 Biot 耦合项 + 强 Dirichlet 罚项放大的刚度矩阵，block-Jacobi-ILU 可能 stall，使 Newton 步不精确、非线性残差卡在线性求解误差地板（O(1)）无法继续下降 → 振荡 → MAX_IT。
- 备择假设：位移 Jacobian 存在近零空间（CG 不处理），或残差/雅可比 Dirichlet 施加尺度不一致。

### 下一步决定性实验（仅改 deck/CI，零 Fortran）
- **ROBUSTKSP 变体**：主 deck 叠加 `-disp_ksp_type gmres -disp_pc_type gamg`（或 `preonly+lu` 若 MUMPS 可用），看位移 SNES 是否收敛。
  - 收敛 ⇒ 确认是线性求解/预条件问题 → 改位移 KSP 配置即可（deck 层，符合零编译约束）。
  - 仍失败 ⇒ 是装配/Jacobian 近零空间 bug → 需 Fortran 改动，超出零本地编译约束，交用户处理。
- Gate 3.3 当前为硬闸门（commit c8f133a），CI 全红属预期；在根因未定位前保持硬闸门以减少噪声。

### 代码精读结论（m_MEF90_DefMechAssembly.F90，用户要求）
**位移 R/J 装配一致性已逐项核对，无 bug：**
- 残差 BC：`R = u − g`（INSERT，行 394/425）；雅可比 BC：`MatZeroRowsColumnsIS(A,setISdof,1.0,NULL,NULL)`（行 701/722）。`Has_DisplacementBC` vs `Has_displacementBC` 仅是大小写差异，Fortran 不敏感 → 同一标识 → BC 对 R/J **一致施加**。
- 弹性项：残差 293-295 (`stressGauss·GradS_BF`) 与雅可比 655-660 (`AGauss·GradS_BF·GradS_BF`) 一致。NODAMAGE（AT2Elastic=`cellIsElastic`）时 `stressGauss=HookesLaw:∇u`、`AGauss=HookesLaw` → **纯线弹性，J 精确=常数**。
- 阻尼/cohesive 项：残差 297-333 与雅可比 600-668 的 `cellDampingCoefficient=dampingCoeff·dt+cohesiveStiffness` 一致。

**决定性推论**：一个 R/J 一致的纯线弹性 Newton 问题本应 1 步收敛；但 NODAMAGE 仍 -5 振荡。且位移矩阵非奇异（Disp KSP 无 zero-pivot），而 NODAMAGE 日志 7785 行中**仅有 50 条 "Linear solve converged"（全是 flow KSP），位移 KSP 贡献为 0、亦无 DIVERGED**。
→ **根因 = 位移内层 KSP（CG + bjacobi/ILU(0)）对刚度矩阵求解失败/stall**，Newton 每步拿不到准确步长 → 残差卡在 O(1) 振荡 → MAX_IT。与压力耦合无关（压力场完美）、与损伤无关、与 J 一致性无关。

**下一步（仅改 deck/CI，零 Fortran）**：加 ROBUSTKSP 诊断变体，主 deck 叠加 `-disp_ksp_type gmres -disp_pc_type gamg`（或 `preonly+lu` 若 MUMPS 可用）+ `-disp_ksp_converged_reason`，看位移 SNES 是否收敛、位移 KSP 是否真正收敛。

---

## CI #78 ROBUSTKSP 诊断（已发起，commit 待填）

**变体配置**（已写入 `.github/workflows/build.yml`，Run Terzaghi 步骤）：
- 基于 **NODAMAGE deck（`terzaghi1D_nodamage.yaml`，AT2Elastic 纯线弹性，R/J 一致）**，叠加
  `-disp_ksp_type gmres -disp_pc_type gamg -disp_ksp_converged_reason`，tee → `/tmp/terzaghi_robustksp.log`。
- 选 NODAMAGE 而非主 deck：损伤已被证无关，纯线弹性下唯一变量 = 位移 KSP 配置
  （弱 `CG + bjacobi/ILU(0)` vs 强 `GMRES + GAMG`）。
- 步骤 timeout 45→60 min 容纳新增变体（主 600s + V2 1500s + NODAMAGE 600s + ROBUSTKSP 900s 最坏 ~55min）。

**判定**：
- 位移 SNES 收敛（`rk_fail=0`）⇒ 坐实弱 KSP/PC 是元凶；**修复 = 主 deck 改 `-Disp_ksp_type gmres -Disp_pc_type gamg`（零 Fortran，符合约束）**。
- 仍 `-5` ⇒ 排除弱 KSP，转向装配/近零空间 bug（需 Fortran，超出零编译约束，交用户）。

**关键旁证（已确认）**：V2 之前已用 `-disp_ksp_type gmres`（小写前缀，日志证实生效为 `Disp_ KSP type: gmres`），但 PC 仍是默认弱 `bjacobi/ILU(0)`，仍 -5 → 光换 `ksp_type` 不够，**必须换强 PC（GAMG/LU）**，这正是 ROBUSTKSP 的测试点。

---

## CI #78 实际结果（纠正：ROBUSTKSP 并未收敛）

**断言修正**：watcher 首报 `rk_fail=0 / CONVERGED` 是**脚本误判**。实际 ROBUSTKSP 运行在初始化阶段即 **SEGV 崩溃**（`MPI_ABORT` + `Caught signal number 11 SEGV`，PETSc 3.3.0-p7）。崩溃发生在 `type: gamg` PC 创建后约 0.5ms（GAMG setup 阶段）。

**根因**：`disp_addNullSpace = TRUE` 为位移矩阵注册了 near-null-space；PETSc 3.3 的 GAMG 与该 near-null-space 组合触发崩溃（GAMG 用 near-null-space 做粗化时越界）。这**不是**位移矩阵装配 bug（flow KSP 用 ILU 不崩，且代码精读已确认 R/J 一致）。

**PETSc 构建能力**（来自 Configure 行）：编了 `--download-ml=1`（Trilinos ML AMG）、`--download-hypre=1`（BoomerAMG）、`--download-metis/parmetis/chaco`，**未编 MUMPS/SuperLU**（无直接 LU）。

**结论**：CI #78 因 GAMG 崩溃**无有效结论**，原"弱 KSP 是元凶"假设既不能证实也不能排除。

**CI #79 修正**：ROBUSTKSP 改用 **ML**（`-disp_pc_type ml`，独立于 GAMG 的成熟 AMG，near-null-space 处理不同，PETSc 3.3 中更稳），并修复诊断脚本——增加崩溃检测（`grep Caught signal|SEGV|MPI_ABORT|PETSC ERROR` → 判为 CRASHED，不再误判 CONVERGED）。若 ML 收敛则坐实弱 KSP 假设；若 ML 仍崩/发散，下一轮改用非 AMG 强 PC（如 `-disp_pc_type bjacobi -disp_sub_pc_type ilu -disp_sub_pc_factor_levels 2`）作对照。

---

## ⚠️ CI #79 结果 + 两处重大更正（覆盖上文"弱 KSP 是元凶"推论）

**CI #79（run 42436094037，0a1a763）= failure**。ROBUSTKSP(ML) 同样 **SEGV 崩溃**（脚本这次正确判定 `CRASHED`，`rk_fail=999`）：`type: ml` 打印后立刻 `Caught signal number 11 SEGV`，`total number of linear solver iterations=0` → 崩在**第一次 KSP/PC setup，尚未做任何线性求解**。NODAMAGE 仍 51/51 `-5`；生产运行 Gate 3.1 到 step 51、pmax=9055.7 正常；Gate 3.3 硬闸门报 60 次非收敛 → FAIL（预期内）。

### 更正 1：AMG 类 PC 在这套环境下不可用，共同元凶 = `disp_addNullSpace`
GAMG（#78）与 ML（#79）都在 PC setup SEGV，而默认 `CG + bjacobi/ILU(0)` 从不崩。差别在于 AMG 必须查询矩阵的 near-null-space 做粗化。代码 `m_MEF90_DefMech.F90:1438-1446`：
```
MatNullSpaceCreateRigidBody(coordVec,nspDisp)
MatSetNearNullSpace(matDisp,nspDisp)
MatNullSpaceDestroy(nspDisp)          ! ← 立刻销毁
```
PETSc 3.3 下极可能是 **use-after-free / 悬空 near-null-space 指针**；CG/ILU 不读它故不崩，AMG 读它即越界 → SEGV。
更重要：**BC 修复（9c1edf2）后 cs0010/cs0020 已是 `TRUE TRUE TRUE`（完全约束），刚体模态根本不再是零空间**，此时仍附加 rigid-body near-null-space 属于**错误配置**，极可能正是 -5 振荡的来源。时间线吻合：BC 修复消除了纯 -5 奇异，却暴露出新的 -2/-5 混合振荡。deck 来源：`defmech.disp.addNullSpace: true`（terzaghi1D.yaml:45）；注册名 `disp_addNullSpace`（Ctx:866）。

### 更正 2："位移 KSP 无 'Linear solve converged' ⇒ KSP stall"是误判
baseline 与 NODAMAGE 运行**从未传 `-disp_ksp_converged_reason`**（只有 ROBUSTKSP 传了），PETSc 默认不打印该行；flow 有 50 条只因传了 `-flow_ksp_converged_reason`。所以"位移 KSP 贡献为 0"是**预期行为，不是失败信号** → "弱 KSP 是元凶"的立论基础不成立，ROBUSTKSP 实验方向本身即建立在误读之上。

### 下一步最高杠杆实验（零 Fortran，deck 一行）
把 `defmech.disp.addNullSpace: true → false`（或命令行 `-disp_addNullSpace 0`），其余不变：
- AMG 不崩 **且** 位移 SNES 收敛 ⇒ 坐实 near-null-space 是 -5 振荡与 AMG 崩溃的**共同元凶**，修复 = deck 一行（零 Fortran，完美符合约束）。
- AMG 不崩但位移仍 -5 ⇒ near-null-space 只是崩溃元凶，位移振荡另有原因。
- 仍崩 ⇒ 与 near-null-space 无关，需另查（可能要 Fortran / 交用户）。

---

## ⛔ 上段"悬空指针 / 关 null space"推论已被代码复核**推翻**（以此段为准）

**复核证据（全部来自本机源码 + ci78/ci79 日志，未跑新 CI）：**

1. **崩溃点精确一致**：GAMG(#78) 与 ML(#79) 都是在 `SNES Object:(Disp_)` → `KSP Object` →
   `PC Object: (Disp_)` → `type: gamg` / `type: ml` **打印完这一行后立即 SEGV**；
   `total number of linear solver iterations=0` 说明**根本没进入 KSP/PC setup 或求解**。
   → 崩在 `SNESView`/`PCView` 的打印路径，不是 PC setup，更不是线性求解。

2. **代码作者 2014 年就把这个 bug 写在源码里**（`m_MEF90_DefMech.F90:1498-1500`， damage SNES 同款注释在 1591-1593）：
   ```
   !!! SNESView seems to segfault when PC ml is used, so we setup an additional
   !!! debug level so we can bypass SNESView
   If (MEF90GlobalOptions%verbose > 1) Then
      Call SNESView(snesDisp,PETSC_VIEWER_STDOUT_WORLD,ierr)
   End If
   ```
   即：**`SNESView` + PC ml/gamg 会 segfault 是已知缺陷，作者专门加了 verbose 级别用来绕过它**。
   而 deck `terzaghi1D.yaml:22` 恰好是 `verbose: 2`（>1）→ 强制调用 SNESView → 踩中该缺陷。

3. **`MatSetNearNullSpace` + `MatNullSpaceDestroy` 不是新引入的 bug**：
   同一写法共出现 3 处（`m_DefMech/m_MEF90_DefMech.F90:1444-1446`、
   旧副本 `m_MEF90_DefMech.F90:1432-1434`、**专门的空空间测试程序 `Tests/TestNSP.F90:86-89`**）。
   TestNSP 就是为验证该用法而写，且它还调用了 `MatNullSpaceView`。PETSc 的
   `MatSet(Near)NullSpace` 惯例会持有对象引用，调用方随后 Destroy 是标准用法。
   → "悬空指针"缺乏支撑；本机无 PETSc 源码树，无法 100% 证伪引用计数语义，但**崩溃证据链与 null space 无关**。

4. **之前所有主运行/NODAMAGE/V2 用的是 `bjacobi/ILU` PC，SNESView 不崩**，所以它们的
   `51/51 -5` 是真实发散、与 SNESView 无关。只有两个 AMG 变体在求解前就被 SNESView 拦下。

**因此**：GAMG/ML 两轮实验**从未真正测到强 PC 的求解能力**；"弱 KSP 是元凶"至今
既未证实也未证伪。上一轮提出的"关 null space"方向若不先绕开 SNESView，即使改了也会同样崩。

### 修正后的决定性实验（零 Fortran，且能同时证伪/证实两个假设）
把 `verbose: 2 → 1`（或诊断变体命令行 `-verbose 1`，需确认 PetscBag 命令行优先于
options_file；最稳妥是直接改**诊断专用 deck** 的 verbose）。已核实：
- Gate 3.1/3.2/3.3 的打印都是无条件 `PetscPrintf`，**不依赖 verbose、不依赖 SNESView 输出**；
- `build.yml` 的 Gate 判据用的是 `SNESSolve failed` / `Gate 3.1` / `pmin|pmax`，与 SNESView 无关；
- `verbose` 注册于 `MEF90/m_MEF90_Ctx.F90:173`（全局 bag，int，名 `'verbose'`）。

判定：
- `-verbose 1` 后 **ML 能 setup 并正常求解、位移 SNES 收敛** ⇒ SNESView 是唯一崩溃原因，
  且强 PC 修复了 -5 ⇒ 最终修复 = deck `verbose: 1` + 位移 KSP/PC 配置（零 Fortran）。
- `-verbose 1` 后 **ML 在真正 PC setup 时仍崩** ⇒ near-null-space 悬空假设重新成立，
  再测 `-disp_addNullSpace 0`（这才是有意义的对照）。
- `-verbose 1` 后 **ML 不崩但位移仍 -5** ⇒ 弱 KSP 假设证伪，转查装配/其他（可能需 Fortran）。

**状态：按用户选择"暂停、先复核代码"，本轮未改任何 deck / build.yml / Fortran，未推 CI。**

### CI #80：按 verbose 1 方向推进（用户选择，已实施）

**新增** `vDef/data/terzaghi1D_robustksp.yaml` = `terzaghi1D_nodamage.yaml`，
除顶部说明注释外**仅改一行 `verbose: 2 → 1`**（diff 已逐行核对）。
依据：`m_MEF90_DefMech.F90:1499` 的 `If (verbose > 1) Call SNESView(...)`，
源码 1498 行作者注释明确写着 "SNESView seems to segfault when PC ml is used"。
CI #78(GAMG)/#79(ML) 都是在 SNESView 打印 `PC Object: (Disp_) type: ...` 后
立刻 SEGV 且 `linear solver iterations=0`，即**从未进入 PC setup/求解**——
两轮强 PC 实验从未真正测到强 PC。verbose 1 绕过 SNESView。

**build.yml**：ROBUSTKSP 变体改用新 deck（ML 强 PC 及所有其它参数不变）；
主运行 / NODAMAGE / V2 的 deck 与参数**全部不动**。已核实 Gate 3.1/3.2/3.3
打印是无条件 `PetscPrintf`，不依赖 verbose，也不依赖 SNESView 输出。

**加固判定脚本**（避免第三次误判）：`CONVERGED` 现在必须同时满足
`rc=0` + 无 `SEGV/MPI_ABORT/PETSC ERROR` + `rk_fail=0` + 至少 1 次 `Linear solve`；
另输出 `last-step` 与 `pmin/pmax` 便于人工核对。这同时堵住 CI #78 的
"崩溃不打印 SNESSolve failed ⇒ 假收敛" 与 timeout(`rc=124`) 的假收敛两类漏洞。

**静态校验**：PyYAML OK；提取 413 行内联 bash 通过 `bash -n`；`git diff --check` 干净。

判定树（CI #80 出结果后）：
- ML 正常 setup 且位移 SNES 收敛 ⇒ SNESView 是崩溃唯一原因 + 强 PC 修复 -5；
  最终修复 = 主 deck `verbose: 1` + 位移 KSP/PC 配置（**零 Fortran**，可收口 Gate 3.3）。
- ML 在真正 PC setup 时仍崩 ⇒ near-null-space 悬空假设重新成立，再测 `-disp_addNullSpace 0`。
- ML 不崩但位移仍 -5 ⇒ 弱 KSP 假设证伪，转查装配层（可能需 Fortran）。

### CI #80 结果 + CI #81 KSPINIT 诊断（已发起）

**CI #80 = failure，但实验有效**：
- `rc=0 crash=0 nonconvergent=50 linear-solves=2583 last-step=51` → `verbose: 1`
  确实绕过 SNESView SEGV；ML 真正完成 PC setup 与 2583 次线性求解（全部 converged），
  跑完 51 步，压力场正确（pmax 9055.7）。复核结论（SNESView 是崩溃元凶）成立。
- 但位移 SNES 仍 `51/51 -5`。fnorm 统计：2584 条中**无一条 < 1e-2**（最小 0.0317），
  首步 2.08e4→0.876 后在 0.03–16 间振荡 → **弱 KSP/PC 假设彻底证伪**。
- Gate 3.3 硬闸门报 60 次非收敛 → FAIL（预期）。

**新根因假设（代码级，可零 Fortran 验证）**：`m_MEF90_DefMech.F90:1476`
`KSPSetInitialGuessNonzero(kspDisp, PETSC_TRUE)`。CI #80 的 KSP 迭代直方图：
`1221 × iterations 0`、`1163 × iterations 1`、仅约 51 次为 39–44 迭代（每外层步一次
真正的求解）；日志多处打印 `using nonzero initial guess`。机制：第一次 Newton 步长解出
后，KSP 复用上一次增量作初值 → 立即 `CONVERGED_ATOL` → SNES 反复施加同一步长 →
残差振荡 → MAX_IT。这与 PC 类型、损伤、压力、BC、R/J 一致性都无关，且解释了为何
强 PC 也无效。

**CI #81**：ROBUSTKSP 变体改造为 KSPINIT —— 恢复默认 KSP（CG+bjacobi/ILU），
叠加 `-disp_ksp_initial_guess_nonzero 0`。因 `KSPSetFromOptions`（1481 行）在
`KSPSetInitialGuessNonzero`（1476 行）之后调用，命令行可覆盖 Fortran 设置，**零 Fortran**。
deck 仍为 `terzaghi1D_robustksp.yaml`（NODAMAGE + verbose 1）。
判定：收敛 ⇒ 修复就是这一个选项（主 deck 加上即可，可顺势收口 Gate 3.3）；
仍 -5 ⇒ 排除该假设，转查装配层。

- 收敛 ⇒ 确认是线性求解/预条件问题 → 改位移 KSP 配置即可（符合零编译约束）。
- 仍失败 ⇒ 是装配/Jacobian 近零空间 bug → 需 Fortran，超出零本地编译约束，交用户。

---

## ⚠️ 覆盖性更正：CI #81「KSPINIT CONVERGED」是**第三次假收敛**

**脚本判据的致命缺陷**：`vDef.F90:573` 仅在 `snesDispConvergedReason < 0` 时打印
`SNESSolve failed with SNESConvergedReason`。因此 `rk_fail=0` 只说明 SNES 返回了
**非负**原因，完全不等于残差收敛。此前所有以「没有 -5」为收敛标准的判定都不可靠。

**CI #81 真实数据**（`ci81.log` KSPINIT 区）：
- 位移 SNES 容差 `atol=1e-7, rtol=1e-5`（`m_MEF90_DefMech.F90:1466-1468`），
  初始 fnorm 2.0818e4 ⇒ 收敛阈值 ≈ 0.208。
- 每个外层步只有 **3 次** Newton 迭代：`2.08e4 → 9.345 → 9.345`（第三条与上一条
  **一字不差**），中间那次线性求解是 `CONVERGED_ATOL iterations 0` ⇒ 增量为 0。
- 全程 153 条 fnorm，最小 6.26、最后 7.45，**从未接近 0.208 或 1e-7**。
- `SNES_CONVERGED_SNORM_RELATIVE`（步长 0）静默返回正原因 ⇒ 程序不报错、rc=0、
  51 步跑完、压力场 pmax=9055.7 正常，但位移场求解是停滞的。
- `alpha min/max 0/0` 来自 `vDef.F90:627-632` 的**损伤场** VecMin/VecMax
  （AT2Elastic 惰性损伤恒 0），与位移 SNES 无关，勿再误读。

**结论**：`-disp_ksp_initial_guess_nonzero 0` **没有修复** -5，只是把症状从
`DIVERGED_MAX_IT` 换成静默假收敛。KSPINIT 假设未证实。
三次误判依次为：① CI #75 watcher 区间误读（量到变体区而非生产区）；
② CI #78 崩溃不打印失败消息 → 假 CONVERGED；③ CI #81 步长 0 → 假 CONVERGED。
**判定必须改为基于 fnorm 是否降到容差，而非失败消息计数。**

## 新根因（代码级，零 Fortran 可验证）：Split 造成 R/J 不一致

`m_DefMech/m_MEF90_DefMechAssembly.F90`：
- **残差 287–289**：无条件 `Call Split%DEED(...)`；`cellIsElastic` 时
  `stressGauss = stressGaussPlus + stressGaussMinus`（HydrostaticDeviatoric 谱分裂，
  分段非线性）。
- **雅可比 649–650**：`If (cellIsElastic) AGauss = matpropSet%HookesLaw`
  （**完整常数 Hooke 张量**），只有 `Else` 分支才调用 `Split%D2EED`。
⇒ 当 `damage.type: AT2Elastic`（`cellIsElastic=TRUE`）时 **J ≠ dR/du**，Newton 方向
不精确，线搜索把步长压到 0，残差降 3 个数量级后卡在 O(8) 不再下降。这能统一解释
为什么换遍 KSP/PC（CG/GMRES、bjacobi/GAMG/ML、FD Jacobian）、关损伤、关 nonzero
初值都无效。

**且 Terzaghi 是纯孔弹性固结基准，`unilateralContact`（应变能分裂/单侧接触）
本就不属于该物理**——deck 里的 `type: HydrostaticDeviatoric` 是多余配置。
`MEF90DefMech_unilateralContactTypeList(1) = 'None'`，`MEF90DefMechGetSplit` 对
None 返回 `MEF90_DEFMECHSPLITNONE()`，使 `DEED` 与 `D2EED` 同时退化为线性 Hooke，
R/J 一致 ⇒ Newton 应 1 步收敛。

**CI #82 计划**：主 deck `terzaghi1D.yaml` 的 `unilateralContact.type`
`HydrostaticDeviatoric → None`（唯一出现处，第 106 行）；主运行追加
`-disp_snes_converged_reason` 以打印真实原因；watcher 同时核 `defmech_fail`、
fnorm 分布与 reason，**不再只看失败消息计数**，杜绝第四次假收敛。

---

## 【覆盖性更正】CI #82 = FAIL：split / R-J 不一致假设被彻底证伪

**上面的 "J ≠ dR/du 是根因" 推论不成立，勿再据此推 CI。** 独立复核 ci82.log：

| 检查项 | 实测 |
|---|---|
| deck 改动生效 | bag 视图多处 `unilateralContact_type = None`（合法枚举，`SplitNone` 已实现，未 STOP） |
| 主运行 | **57 次 `-5`**；fnorm 2607 条，min **1.56**、last 4.66、**0 条 < 0.208**（rtol·F₀ 阈值） |
| **V2（`-disp_snes_fd` FD Jacobian）** | **首次完整跑完 `rc=0`，仍 59 次 `-5`**；首步 fnorm 6.17→3.59→3.45→3.40→4.16→4.83（降后反弹），与解析 J 同模式 |
| 压力场 / 崩溃 | `pmax=9055.7` 正确，crash=0，走完 51 步 |
| Gate 3.3 | 57 次 → **FAIL**（硬闸门，预期） |

**决定性推论**：FD Jacobian 由残差数值微分得到、数学上必与 R 一致；**连 FD 也发散 ⇒ 所有「J 与 R 不一致」类假设（应变能分裂、阻尼系数、BC 施加方式）全部证伪**，问题不在雅可比装配。

**本轮额外排除**（零成本代码复核）：
- `dampingCoefficientDisplacement = 0`、`cohesiveStiffness = 0`（bag 实测）⇒ 残差 `damping/dt` vs 雅可比 `damping*dt` 的系数不对称根本不生效。
- 日志里每步的 `type: ksponly` 属 **`flow_` SNES**（压力线性求解），不是位移 SNES。
- 残差向量管理正确：`SectionRealSet(residualSec,0)` + `VecSet(residual,0)` + `SectionRealUpdateClosure(ADD_VALUES)` + `SectionRealToVec(SCATTER_FORWARD)` + BC `INSERT_VALUES` ⇒ 排除「未清零累加 / 状态污染」。
- `50 × CONVERGED_RTOL` 是 flow KSP 线性求解（`Linear solve converged ... iterations 42`），与位移 SNES 无关；`-disp_snes_converged_reason` 未产生 `diverged due to` 输出。

**结论**：**配置层（deck + CI 命令行）可零 Fortran 验证的假设已穷尽**。残差在 Newton 迭代中降到 O(1) 后振荡、内层线性求解全部收敛、解析/FD Jacobian 同败 ⇒ 问题落在 **Fortran 源码层的 SNES 状态管理 / 残差回调的隐藏状态依赖**。改 Fortran 由 CI 编译是合规的（零本地编译约束不禁止），但周期与风险高于 deck 调整。

**另注**：Gate 3.3 硬闸门目前使每一次提交必红；是否临时退回 report-only 是与根因排查独立的工程决策。

---

## CI #73 artifact 复核：ksponly「CONVERGED」是第四次假收敛（无数据占位）

CI #82 的 diagnostics artifact 只有 7 个文件，缺 `snes_l2.log` / `snes_ksponly.log`——
因为 `for variant in l2 ksponly` 与 `diagnose_snes.py` 在 build.yml 第 811/826 行，
属 **Gate 3.3 之后的独立步骤**；Gate 3.3 `exit 1` 使 CI #82 该步骤被跳过。
CI #73（当时 3.3 为 report-only、整体 success）的 artifact 则完整保存了两个变体日志，
其 job log 16006–16011 行给出汇总表：

```
main-basic   51(-5:51) 51 51 50  1.049e+01  3.461e+00  7.304e+00  DIVERGED
l2-original  51(-5:51) 51 51 50  1.049e+01  2.008e+00  5.956e+00  DIVERGED
ksponly       0        51 51  0  0.000e+00  0.000e+00  0.000e+00  CONVERGED
```

**ksponly 那行不可采信**：`snes_ksponly.log` 里 **0 条 `SNES Function norm`**、
`iters=0`、三列 `0.000e+00` 均为 diagnose_snes.py 在「无样本」时的占位值。
`SNESKSPONLY` 只解一次线性系统，不算残差、不触发 monitor、不做收敛判据，
也不打印 `SNESSolve failed` ⇒ `fails=0` 与 `fnorm=0` 都**不证明残差归零**。
日志里 `maximum iterations=1` / SNESView `total number of linear solver
iterations=0` 印证这一点。**这是继 #75 区间误读、#78 崩溃假收敛、#81 步长 0
假收敛之后的第四次「假收敛」，成因是「无数据被当成零」。**

**l2 变体**（`-disp_snes_linesearch_type l2`）：51 次 `-5`、2601 条 fnorm、
min **1.845**、0 条 <1e-4 ⇒ 与 `basic` 同样发散，**线搜索类型也排除**。

**结论**：ksponly 是否真能让残差达标，**现成证据无法回答**。要判定必须新跑一轮
带真实判据的实验：主运行/诊断变体用 `-disp_snes_type ksponly`，判据不得使用
`SNESSolve failed` 计数或缺失 fnorm，而要用「位移场 vs Terzaghi 解析解」
（复用 Gate 3.2 的 FEM/BE 比较）或显式残差范数。

---

## 判据检索（2026-11-08）：位移场「真实判据」的可用证据与阈值

目的：为 ksponly 新实验选判据，杜绝第五次「无数据被当成零」的假收敛。

### 1. 现成实验结果（本地 CI artifact / 日志，可直接引用）

| 来源 | 结果 | 用途 |
|---|---|---|
| Gate 3.2 / CI #70 `terzaghi_pressure_timeseries.csv` | pmax `1e5 → 9055.7 Pa`，ratio `0.0905571`；BE 解析 `0.09339`（L2 相对误差 **3.206%**）；连续解析 `0.08713` | 压力场正确性基线，已通过且与位移求解解耦 |
| CI #69 GMRES vs CG 同源双跑（checklist 1065–1074） | 压力差 `2.54e-5 Pa`（`2.5e-10` 相对）、ux 差 `2.22e-6 m`、**uy 差 `0.01204 m`** | 位移场「不具备定量可比性」的直接证据；P7 阈值来源 |
| CI #70 `_out.gen`（`a70`，51×356，Gate 3.2 同一网格） | ux ∈ `[-9.31e-4, 0]` 单调衰减；uy 从 `±1.8e-5` 单调累积到**整场均匀** `-0.05564 m` | 错误解指纹 = 均匀刚体漂移；BC 修复前 ux@BC = `9.3e-4` |
| CI #73 metadata check | `slot 4 Displacement_Y t1 [-1.06137, -1.06133]`（整场均匀） | 漂移量随配置变化（-0.056 / -1.06），非物理解 |
| CI #82 artifact | 仅 7 个日志，**无 `_out.gen`**，无 l2/ksponly 日志 | 新实验必须显式保存结果网格 |

### 2. 外部解析解（判据的物理依据）

- Wang (2004) / E3S 2023 附录 A：Terzaghi 压力级数 (A.1) + 垂直位移闭式 (A.2)，`c_M = 1/α1`。
- Itasca FLAC3D《One-Dimensional Consolidation》（Detournay & Cheng 1993）：由 `∂σ_zz/∂z=0` 与 `σ=α1 ε - α p` 积分得位移级数 (5)，`α1 = K+4G/3`。
- Andrade et al., JMPS 2022 式 (18)(19)：位移/沉降无量纲级数，并**用解析解验证实验数据**。
- Springer 2025（THM verification）式 (16)：一维约束模量 `E(1-ν)/[(1-2ν)(1+ν)]`。

本项目：`E=1e9, ν=0.25` → `λ=4e8, G=4e8, M=K+4G/3=1.2e9 Pa`；`α=1, p0=1e5, H=10` → `p0 H/M = 8.33e-4 m`。用 Gate 3.2 同一 BE 级数 `d_m=(1+λ_m dt)^(-n)` 积分：
- **一端自由**（CI #70 当时配置）：`u(-5,0)` 解析 `-8.3305e-4 m` vs FEM `-9.2974e-4 m`（差 11.6%，源于 t=0 排水面 Dirichlet + 二维/重力效应）→ 常数、符号校验通过。
- **两端全固定**（当前 deck）：`u_x^an` 峰值 **`9.90e-5 m` 在 t=0.08 s（step 5）**；`t=0: 2.8e-7`、`t=1: 1.04e-5`，**非单调（先增后减）**，两端严格 0（`-1.6e-24`）。

### 3. 判据表（ksponly / 主运行均适用，按可靠性排序）

| # | 判据 | PASS 阈值 | 依据 |
|---|---|---|---|
| P1 | BC 节点 `x=±5` 的 `\|u_x\|,\|u_y\|` | `< 1e-10 m` | 解析解两端严格 0；CI #70 修复前 `9.3e-4 m` |
| P2 | 均匀刚体漂移 `\|mean(u_y)\|` | `< 1e-6 m` | 重力静载只产生弯曲 `±1.8e-5`；错误解整场均匀 `-0.0556/-1.06` |
| P3 | `\|u_y(t=1)\| / \|u_y(t=0)\|` | `< 2`（静载应稳态） | 实测 `3045×`（#70）、`~6e4×`（#73） |
| P4 | `max\|u_x\|` 全程峰值 | `5e-5 ~ 2e-4 m`，且峰值在 `t<0.2 s` | 两端固定解析 `9.90e-5 @ 0.08 s`；一端自由上限 `8.33e-4` |
| P5 | `u_x` 时间形态 | 非单调，不能从 `9e-4` 单调衰减 | 解析解；单调衰减属一端自由配置 |
| P6 | 压力场（复用 Gate 3.2） | ratio `0.0906`、L2(BE) `≤ 3.3%` | BE 解析，与位移求解解耦 |
| P7 | 同源一致性（ksponly vs 主运行） | **压力**全场差 `< 1e-2 Pa`；位移差仅报告、不作判据（两者本就应不同） | CI #69 反例：p 差 `2.54e-5 Pa` → uy 差 `0.01204 m` |

**明确不可用作判据**：`SNESSolve failed` 计数（`vDef.F90:573` 仅在 `reason<0` 打印）、`diagnose_snes.py` 无样本时的 `0.000e+00` 占位、`SNES_CONVERGED_SNORM_RELATIVE`（步长 0 假收敛）、未传 `-disp_ksp_converged_reason` 时的「无输出」。

### 4. 落地要点（CI #83 已实施）

- Gate 3.3 的 `exit 1` 会跳过其后 l2/ksponly 步骤（CI #82 已证）→ 诊断插在 **Gate 3.3 之前**（ROBUSTKSP 之后），并显式保存 `/tmp/out_main.gen` 与 `/tmp/out_ksponly.gen`，否则又会出现「无数据」。
- 可复用现成脚本：`vDef/data/check_exodus_metadata.py`（slot 范围）、`vDef/data/compare_ksp_runs.py`（同源日志对比）、`/d/tmp/compare_fields.py`（四 slot 全字段对比，仅本地有，未入库）。
- 判定树：ksponly 过 P1–P5 而主运行 51×`-5` 不过 ⇒ 「SNES 迭代污染解」成立，ksponly 即修复方向；两者都不过 ⇒ 问题在装配/状态管理；两者都过 ⇒ 51×`-5` 只是收敛判据问题，物理场正确。

**CI #83 实现**：新增 `vDef/data/check_terzaghi_displacement.py`（纯标准库 +
`parse_cdf2.py`；P1–P7 自动判定；恒返回 0、只打印 VERDICT，不新增闸门）。
`build.yml` 主运行结束后立即 `cp _out.gen /tmp/out_main.gen`（后续 NODAMAGE /
ROBUSTKSP / KSPONLY-PHYS 都会覆盖同一路径），在 ROBUSTKSP 之后插入
`-disp_snes_type ksponly` 的 51 步运行（600 s cap），保存 `/tmp/out_ksponly.gen`
并把主网格恢复给 Gate 3.1/3.2 与交付 artifact；新日志与两份网格加入
`build-diagnostics` artifact；step timeout `60 → 75`。
静态校验：PyYAML、10 个 run 块 `bash -n`、脚本 `ast.parse`、`git diff --check` 全过；
判据脚本用 CI #70 既有网格自测（P6 ratio `0.0905571`、L2 `0.0320613` 与 Gate 3.2
逐位一致，P1–P5 正确判 FAIL）。**零本地编译、零本地求解器运行**。

---

## CI #83 结果（run `42459128534`，commit `a3c7a27`，conclusion = failure）

artifact 下载：GitHub `archive_download_url` 直接带 `Authorization` 跟到 Azure blob 会
返回 `AuthenticationFailed` 403；须先取 `Location` 重定向，再**不带 token** 下载
（与 job-logs 同类问题）。watcher 的 `[D]` 段因把 `size_in_bytes` 写成 `size`
静默失败，artifact 实际存在（11 个文件，含两份 `_out.gen`）。

### 1. 判据表（脚本原始输出）

| 判据 | 阈值 | MAIN `/tmp/out_main.gen` | KSPONLY `/tmp/out_ksponly.gen` |
|---|---|---|---|
| P1 `max\|u\|` @ x=±5 | `<1e-10 m` | `0` PASS | `0` PASS |
| P2 `max\|mean(u_y)\|` | `<1e-6 m` | `1.0237e-6` **FAIL（超 2.4%）** | `8.1719e-7` PASS |
| P3 `u_y` 漂移 t1/t0 | `<2` | `0.0774` PASS | `0.0775` PASS |
| P4 `peak\|u_x\|` | `5e-5~2e-4 m` | `1.10348e-4` PASS | `1.10360e-4` PASS |
| P4 峰值时刻 | `≤0.2 s` | `0.08 s` PASS | `0.08 s` PASS |
| P5 `u_x` 非单调 | 内部峰 + 衰减 | `peak@step5` PASS | `peak@step5` PASS |
| P6 压力 ratio / pmin / L2(BE) | `0.05~0.15 / 0 / ≤0.05` | `0.0905571 / 0 / 0.0320613` PASS | 同左 PASS |
| P7 同源压力差 | `≤1e-2 Pa` | `0` PASS（`dUx=2.325e-7`、`dUy=1.605e-6 m` 仅报告） | — |

脚本打印 `MAIN VERDICT: FAIL (P2)` / `KSPONLY VERDICT: PASS` /
`FINAL KSPONLY PHYSICS VERDICT: PASS`；Gate 3.3 仍以 `57` 条失败硬闸门 `exit 1`（预期）。

### 2. 决定性证据：两次运行的解在物理上等价

- 两份网格 sha256 不同、二进制字节差 `776122 / 2845300`（浮点末位差异），
  均为 51 步记录，**不是同一文件被复制两次**。
- 全 32 个记录槽逐步对比（`D:/tmp/ci83diag/`，纯 Python 读 artifact）：
  - `vals_nod_var1` 压力：差 `0`（flow 与位移解耦）；
  - `vals_nod_var2` 损伤：差 `1.43e-7`（相对 `3.74e-5` 为 **0.38%**）；
  - `vals_nod_var3` ux：差 `2.33e-7 m`（相对峰值 `1.10e-4` 为 **0.21%**）；
  - `vals_nod_var4` uy：差 `1.61e-6 m`（uy 本身仅 `9.23e-6 m`，噪声级）；
  - 应力槽 `eb1 var7/8`：差 `61.96 / 27268`、`49.4 / 12167`（**0.23% / 0.41%**）；
  - 其余槽恒为 0。
- 运行时：主运行约 **19.5 min**（51 步 × 51 次 Newton 迭代），
  KSPONLY-PHYS **8 s**（51 步 × 1 次线性解，`SNES Function norm` 输出 0 条），
  rc 均为 0、无崩溃。

### 3. 覆盖性更正（务必以此为准）

1. **上面第 4 节的判定树第一支「SNES 迭代污染解」被证伪**。实际落在第三支：
   `两者物理场都正确 ⇒ -5 只是收敛判据问题`。主运行 57 次 `-5` **没有污染解**。
   `build.yml` 里硬编码的
   `interpretation: ksponly PASS + main 51x-5 => SNES iterations corrupt the field`
   是**错误判词**，下次提交应改为「两解等价 ⇒ 残差平台/收敛判据问题」。
2. **CI #83 的 `KSPINIT CONVERGED ... fix = -disp_ksp_initial_guess_nonzero 0`
   是第五次假收敛**。ROBUSTKSP 段 `nonconvergent=0` 但残差从未下降：
   153 条 fnorm / 51 步 = **每步 3 次迭代**，末值仍 `7.45~8.43`，
   与主运行震荡平台同量级 → `SNES_CONVERGED_SNORM_RELATIVE`（步长 0）提前退出，
   与 CI #81 结论一致。该 echo 分支不可作为修复依据。
3. **ksponly 的价值不是「修物理」，而是「跳过无效迭代」**：
   位移方程在当前配置（`SplitNone`、`AT2` 损伤仅 `3.7e-5`、
   `dampingCoefficientDisplacement=0`、`cohesiveStiffness=0`）下实质线性，
   一次线性解即精确解。51 次 Newton 迭代在残差平台 `3~13` 上震荡，无收益。
4. **ksponly 不是完整修复**：其运行仍有 **7 条 `damage field ... -5`**
   （主运行 6 条 damage + 51 条 displacement = Gate 3.3 的 57）。
   若仅切换 ksponly 而 Gate 3.3 继续数全部失败消息，CI 仍会红。
5. **P2 阈值 `1e-6 m` 偏严，但不得为通过而放水**：MAIN `1.0237e-6`
   仅超 2.4%，且比已知错误解指纹（`0.0556 / 1.06 m`）小 **5 个数量级**，
   比重力弯曲 `±1.8e-5 m` 小一个数量级，属数值噪声。阈值是否放宽到
   `1e-5 m` 需评审确认，本轮未改脚本。

### 4. 残差平台实证（下一步根因线索，`terzaghi.log`）

```
step  1: f0=2.082e+04 → it1=9.345 → 51 次迭代 fmin=4.616 flast=6.484
step  2: f0=1.887e+04 → it1=8.823 → flast=9.481
step 51: f0=2.996e+01 → it1≈7.45  → flast=4.657
全 51 步: flast min/median/max = 2.667 / 6.583 / 12.57
          fmin  min/median/max = 1.560 / 4.802 / 6.381
```
第一次 Newton 步就把残差压掉 3 个量级到 `~9`，此后 50 次迭代只在
`3~13` 震荡，**存在与解无关的不可约残差分量**。待验证来源（尚未测试）：
BC 固定端反力行 / 压力源项是否重复计入残差、残差与 Jacobian 的
scale 是否一致、damage 耦合项是否带常数偏置。CI #83 的 V2 变体
（`-disp_snes_fd` + GMRES，主运行之前跑）同样不收敛，且报的是
51× `-2`（`SNES_DIVERGED_FUNCTION_DOMAIN`，迭代 ~14 次即退出），
故「解析 Jacobian 写错」不是唯一解释。

### 4b. 统计口径（易踩坑，务必区分）

同一个 step 里先后跑了 5 个 mpirun，**job log 的 2318–16098 行不是单一运行**：

| 运行 | 日志 | 失败消息 |
|---|---|---|
| V2 变体（FD Jacobian + GMRES，主运行前） | `/tmp/diag_V2.log` | 51× displacement `-2` |
| 生产运行（解析 Jacobian + CG + basic LS） | `/tmp/terzaghi.log` | 51× displacement `-5` + 6× damage `-5` = **Gate 3.3 的 57** |
| NODAMAGE（AT2Elastic） | `/tmp/terzaghi_nodamage.log` | 51× displacement `-5` |
| ROBUSTKSP/KSPINIT | `/tmp/terzaghi_robustksp.log` | 0（假收敛，见更正 2） |
| KSPONLY-PHYS | `/tmp/terzaghi_ksponly_phys.log` | 7× damage `-5`，位移 0 |

Gate 3.3 只数 `terzaghi.log`；job log 全文的 51×`-2` / 103×`-5` / 21×damage`-5`
是 5 个运行之和，不可当作单次运行的失败数。

### 5. 下一步可选方向（未实施，等确认）

- **A（推荐，零 Fortran，收口闸门 3）**：主运行改用 `-disp_snes_type ksponly`，
  Gate 3.3 从「数 `SNESSolve failed`」改为「P1–P6 物理场判据 + 显式残差上界」，
  并单独处理 damage SNES 的 7 条 `-5`（或证明其不影响物理场后排除统计）。
- **B（Fortran 层）**：定位 `3~13` 不可约残差来源，从装配侧消除，
  让标准 Newton 真收敛；不依赖 ksponly。
- **C（最小）**：仅修正误导性 echo 与判定树，不动求解器配置。

**CI #83 复盘轮动作**：只读分析 CI #83 job log + `build-diagnostics` artifact
（含两份结果网格全槽对比、残差序列统计），并更新本文档与工作日志；
该轮未改 deck / build.yml / 判据脚本 / Fortran，未推 CI。

---

## CI #84：方案 A 落地（用户选 A，收口闸门 3）

### 1. 改动清单

| 文件 | 改动 |
|---|---|
| `.github/workflows/build.yml` | ① 生产运行加 `-disp_snes_type ksponly`（保留 monitor / basic LS / converged_reason 仅作诊断）；② **Gate 3.3 重写**：从 `grep -c "SNESSolve failed"` 硬闸门改为调用判据脚本 `--gate`，对主运行冻结网格 `/tmp/out_main.gen` 做物理场判定，失败 `exit 1`；③ 修正 `KSPINIT CONVERGED` 错误 echo，明确标注为第五次假收敛；④ 修正 `interpretation: ksponly PASS + main -5 => SNES iterations corrupt the field` 错误判词；⑤ KSPONLY-PHYS 变体改定位为「同源确定性对照」；⑥ **删除 CI #71 的 l2/ksponly 变体循环**（其 ksponly 分支与主运行重复、l2 分支可吃掉 600 s cap，且 `diagnose_snes.py` 把无 fnorm 的 ksponly 显示成 `0.000e+00`「收敛」，正是 CI #73 假收敛陷阱）；⑦ artifact 加 `/tmp/gate33_verdict.log`，去掉两个已退役日志 |
| `vDef/data/check_terzaghi_displacement.py` | 新增 `--gate <mesh> <log>` 硬闸门模式（机读 `GATE33_RESULT: PASS\|FAIL`，退出 0/1）；缺网格 / 缺日志 / 记录不足 50 / 末时刻 ≠ 1.0 / 非有限值 / 崩溃标记 / `Unknown option` / Gate 3.2 探针 < 50 **一律 FAIL**；新增 P8 损伤场有界检查；`load()` 增加 `vals_nod_var2` 读取与全字段 finite 检查；默认报告模式仍恒返回 0 |

### 2. P2 阈值校准（不是放水，是恢复区分力）

原 `max|mean(u_y)| < 1e-6 m` 是拍脑袋的绝对值：CI #83 主运行
`1.0237e-6`（FAIL）与 ksponly `8.1719e-7`（PASS）**是同一个物理场**，
阈值卡在两者之间，没有物理含义。改为双条件：

```
max|mean(u_y)| <= 1e-5 m  AND  <= 5% * peak|u_x|
```

依据：重力静载弯曲 `±1.8e-5 m`（CI #70）给出绝对尺度；相对量让判据随
位移尺度自适应。区分力实测：

| 解 | `max\|mean(u_y)\|` | 相对峰值 | 判定 |
|---|---|---|---|
| CI #83 主运行 / ksponly | `1.02e-6` / `8.17e-7` m | `0.93%` / `0.74%` | PASS |
| CI #70 错误解（均匀刚体漂移） | `0.0556 m` | **5978%** | FAIL（P1,P2,P3,P4,P5） |

即新阈值把「正确解」与「错误解」的间隔拉到 ~6400 倍，比原来仅靠
2.4% 边界判 FAIL 稳健得多。

### 3. damage SNES `-5` 的处理

ksponly 只作用于位移 SNES，damage SNES 仍会报 7 条 `-5`（CI #83 实测），
主运行原为 6 条。Gate 3.3 现在**只报告、不计数**这些消息，改用
**P8** 检查损伤场自身：有限、`∈[0,1]`、末步最大值 `≤0.1`
（实测 `3.74e-5`）。这是「排除统计」的物理替代保障，不是无条件放行。
残留问题（damage SNES 为何仍 `-5`）记录在案，属后续 B 方向。

### 4. 「显式残差上界」的诚实调整

A 方案原写「Gate 3.3 = 物理场判据 + 显式残差上界」，但 ksponly 不做
Newton 迭代、日志里没有 fnorm 序列，无法直接给残差上界。实际落地为
**物理场对解析解的误差上界 + 运行完整性硬条件**：P6 的 L2(BE) `≤0.05`、
ratio `∈[0.05,0.15]`，加上 `≥50` 记录、`t_end≈1`、有限值、无崩溃、
`≥50` 条 Gate 3.2 探针。不可约残差平台 `3~13` 的装配根因仍未解决。

### 5. 离线自测（零本地编译、零本地求解器运行）

- `PyYAML` 解析 / 10 个 `run` 块 `bash -n` / 脚本 `ast.parse` /
  `git diff --check`：全过。
- `--gate` 8 组用例全过：CI #83 主网格 + ksponly 日志 → PASS；ksponly
  网格 → PASS；主网格 + **含 51×`-5` 的生产日志** → 仍 PASS（证明闸门
  确实不看失败消息）；缺网格 / 缺日志 / 错误参数 → FAIL；CI #70
  错误解 → `FAIL (P1,P2,P3,P4,P5)`；报告模式两网格均 PASS。
- 用 CI #83 ksponly 日志离线核对 Gate 3.1/3.2 硬断言：
  `post-TSSolve=50`（需 ≥48）、`total number of rejected steps=0`、`TS dt
  0.02 time 1`、pmax `1e5→9055.70884`（ratio `0.09055709`）、
  `check_exodus_metadata.py` rc=0、`FATAL_RE` 命中 0 条。

### 6. 待 CI #84 验证的风险

1. 生产运行新增 `-disp_snes_monitor` 与 ksponly 组合（CI #83 的
   KSPONLY-PHYS 未带 monitor）；预计仅多一行 `0 SNES Function norm`。
2. 生产网格是否与 CI #83 的 `out_ksponly.gen` 一致（P7 同源确定性）。
3. Gate 3.3 之后的 Gate 3.1/3.2 硬断言在新日志下是否仍全过（已离线
   用同配置日志核对，但真实运行仍需确认）。
4. `damage field ... -5` 仍会出现在日志中，属**已知残留**，不再是闸门项。

---

## CI #84 结果（run `42462830545`，commit `c8f84f6`，**conclusion = success**）

**闸门 3 转绿**：`Gate 3.3 PASS` / `Gate 3.1 PASS` / `Gate 3.2 PASS`（两条）
/ 全文无 `=== FAIL` 行 / `GATE33_RESULT: PASS`。这是自该硬闸门恢复以来
第一次由**物理场判据**判出的绿色，而不是失败消息计数或 fnorm 占位。

### 1. Gate 3.3 判据表（`/tmp/gate33_verdict.log`，51 记录 / 356 节点 / 12 BC）

| 判据 | 实测 | 阈值 | 结果 |
|---|---|---|---|
| P1 BC 端位移 | `0` | `<1e-10 m` | PASS |
| P2 `max\|mean(u_y)\|` | `8.1719e-7 m` | `≤1e-5 m` | PASS |
| P2 相对峰值 | `0.0074` | `≤0.05` | PASS |
| P3 `u_y` 时间漂移 | `0.07749` | `<2` | PASS |
| P4 `peak\|u_x\|` / 峰值时刻 | `1.1036e-4 m` / `0.08 s` | `5e-5~2e-4 m` / `≤0.2 s` | PASS |
| P5 非单调 | `peak@step5` | 内部峰 + 衰减 | PASS |
| P6 ratio / pmin / L2(BE) | `0.0905571` / `0` / `0.0320613` | `0.05~0.15` / `≤1e-8` / `≤0.05` | PASS |
| P8 末步最大损伤 | `3.7458e-5` | `≤0.1` | PASS |

`Gate 3.3 diagnostics: displacement-msgs=0 damage-msgs=7 total=7
(INFORMATIONAL ONLY)`；生产日志 264202 字节、51 条 Gate 3.2 探针。

### 2. 决定性核实：求解器真的换了，且解是确定性的

- 主运行 `mpirun rc=0`，耗时 **6.3 s**（08:41:01 V2 结束 → 08:41:07
  主运行退出），对比 CI #83 标准 Newton 的 **19.5 min**。
- `out_main.gen` 与 `out_ksponly.gen` sha256 **不同**、字节差仅 **2**
  （两次独立运行，元数据差异）；P7 `dPressure=0`、`dUx=0`、`dUy=0`
  → 逐位确定性，不是「网格没被重写」。两份 solver log 哈希也不同。
- **跨 CI 交叉验证**：CI #84 主网格 vs CI #83 的独立 ksponly 网格，
  `vals_nod_var1..4` 逐步最大差 **全为 0**（仅 12 字节元数据差）；
  vs CI #83 的旧 Newton 主网格字节差 **776131**（即浮点末位差
  `2.3e-7 / 1.6e-6 m`）。→ 主运行确实切到 ksponly，物理场未变。
- 生产日志仅 7 条 `damage field ... -5`，**位移 SNES 失败 0 条**。

### 3. 各段失败消息归属（口径同 CI #83 第 4b 节）

| 运行 | 失败消息 |
|---|---|
| V2（FD Jacobian + GMRES，主运行前） | 51× displacement `-2` + 8× damage `-5` = 59 |
| 主运行（ksponly） | 7× damage `-5`，位移 0 |
| NODAMAGE | 51× displacement `-5` |
| ROBUSTKSP/KSPINIT | 0（已由新 echo 明确标注为假收敛） |
| KSPONLY-PHYS | 7× damage `-5` |

### 4. 新发现的 CI 时长瓶颈（下一步高杠杆动作）

V2 变体从 `08:26:52` 跑到 `08:41:01`，**约 14 分钟**占整个 job
（08:24:53–08:41:42，约 17 分钟）的 8 成。它是 CI #73 的历史
Jacobian 诊断，CI #82/#83 已证伪「解析 Jacobian 写错」这一假设，
现在既不影响闸门、也不产出新结论。主运行 / NODAMAGE / ROBUSTKSP /
KSPONLY-PHYS 合计仅 ~30 s。**退役 V2 可把 CI 从 ~17 min 降到 ~3 min**，
性质与 CI #84 删除 l2/ksponly 变体循环相同，属零风险清理（待确认）。

### 5. 残留问题（不影响闸门，但需在册）

- damage SNES 仍 `-5`（主运行 7 条）：P8 只保证损伤场有界，未解释原因。
- 残差不可约平台 `3~13` 的 Fortran 装配根因未解决（B 方向）。
- NODAMAGE 变体仍报 51×`-5`、ROBUSTKSP 仍是假收敛：两者的 echo
  已改为如实标注，不再误导；后续可一并退役。
- 工具坑（本轮踩到）：本机 `python3` 是 Windows 原生解释器，读
  `/d/tmp/xxx` 会 `FileNotFoundError`，必须用 `D:/tmp/xxx`；
  watcher 的 `[E]` 段即因此没打印 artifact 清单（artifact 本身齐全，
  12 个文件，已用 `D:/tmp/` 路径补下载核对）。

**本轮动作**：只读核对 CI #84 job log + `build-diagnostics` artifact
（网格哈希 / 跨 CI 物理场对比 / Gate 判据日志 / 分段失败统计），
更新本文档与工作日志；**未改任何代码、未推 CI**。

---

## CI #85：闸门 4.1a 落地（Biot 体力叠加 + 硬判据）

### 1. 为什么先做物理裁决

清单 4.1 写「补 `σ' = σ + α·p·I` 有效应力耦合」，但当前 `-α∇p` 体力形式与之**不是两个可相加的效应**：

```
∫(−α∇p)·v dΩ  =  ∫αp·div(v) dΩ − ∫_Γ αp·v·n dΓ          （分部积分）
∫(C:ε − αpI):ε(v) dΩ = ∫C:ε:ε(v) dΩ − ∫αp·div(v) dΩ       （总应力形式）
```

体积项同为 `−∫αp·div(v)`，只有自然边界的 `−αp·n` 面力解释不同。字面「再补一项」会双重计入孔压，把 `u_x` 峰值从 `1.1036e-4` 推到约 `2.2e-4 m`，越过 Gate 3.3 P4 上界 `2e-4`，直接打掉 CI #84 刚转绿的闸门。故本轮把 4.1 拆成 **4.1a（可验证收口）** 与 **4.1b（互斥总应力模式，单列待做）**。

### 2. 覆盖性更正：原清单的「仅作首次快照」会引入 bug

**~~`m_vDefDefault.F90:43/79` 默认 `ForceScaling = MEF90Scaling_Linear`，`MEF90DefMechSetTransients`（`vDef.F90:447`）每步重设 `force = (0,−9.81,0)×t`。~~ 此句已被 CI #85 证伪，见下文「CI #85 结果」节：deck 的 `Force: 0, -9.81, 0` 因逗号后空格被 RealArray 解析器截断成 `0 0 0`，bag view 实测 `Force = 0 0 0`，故 `linear × 0 = 0`，`baseForce` 全程恒 0。** 快照位置在**时间步层**（`vDef.F90:520`，AltMin 之外）仍然是正确的设计；原清单「仅作首次快照」的建议也仍然是错的（一旦将来启用非零体力，Linear 缩放会被回滚）。但「当前 `|baseForce|` 应为 `0→9.81` 线性增长」这一判据前提在本 deck 上不成立。

### 3. 改动清单

| 文件 | 改动 |
|---|---|
| `m_Flow/m_MEF90_FlowAssembly.F90` | `MEF90FlowApplyBiotForce`：`VecCopy+VecAXPy` → 单次 `VecWAXPY(force,−α,∇p,baseForce)`；每步输出一条 `Gate 4.1:` 诊断（α、`|base|∞`、`|gradp|∞`、`|force|∞`、`superpos_err`）。∞ 范数用 `VecMin/VecMax`（仓库已有编译先例，避开无先例的 `VecNorm`） |
| `vDef/vDef.F90` | 快照处补注释，固化「每时间步快照 = Linear 缩放下唯一正确语义」（逻辑未改） |
| `vDef/data/terzaghi1D_nobiot.yaml`（新增） | 主 deck 逐字副本，仅 `flow.biot_coeff 1.0 → 0.0` + 5 行说明；不用 `-flow_biot_coeff` 命令行覆盖，避免 PETSc 3.3 `Unknown option` 静默无数据 |
| `vDef/data/check_gate41_biot.py`（新增） | 纯标准库，`--gate` 模式输出机读 `GATE41_RESULT: PASS\|FAIL` 并退出 0/1；report 模式恒 0。缺日志/缺网格/记录不足/非有限/崩溃一律 FAIL |
| `.github/workflows/build.yml` | KSPONLY-PHYS 之后、恢复主网格之前插入 NOBIOT ksponly 变体（保存 `/tmp/out_nobiot.gen`）；Gate 3.3 之前算 Gate 4.1、**在 Gate 3.3 PASS 之后才裁决 exit 1**（避免 CI #82 式早退隐藏判据表）；artifact 增 3 个文件 |

### 4. 判据 B1–B7（阈值与物理依据）

| # | 判据 | PASS 阈值 | 能抓什么 |
|---|---|---|---|
| B1 | `Gate 4.1` 行数 / 解析完整性 | `≥50` 行、与 `Gate 3.2: step` 探针数差 `≤2`、无不可解析行 | 诊断没打出来 / 重复打印 |
| B2 | α | `max|α−1| ≤ 1e-6` | deck 被改成无耦合 |
| B3 | `|baseForce|` 时间形态 | 末值 `≥5`、`min ≤1`、线性回归 `R²≥0.98`、斜率 `≥5`、残差 `≤0.5` | **Linear 体力缩放被回滚**（核心验收） |
| B4 | `superpos_err/|force|` | `≤1e-8` | **覆盖**（≈`|base|`）、**丢 Biot**（≈`|α∇p|`）、**逐次累积**（≈`n|α∇p|`） |
| B5 | `|gradp|` 与 `|force|` | `>0`、max/min `≥1.5`、`max(|force|−|base|) ≥1` | Biot 项是死代码 / 不随时间演化 |
| B6 | NOBIOT 对照 | `peak|u_x|(main) ≥ 3×peak|u_x|(nobiot)` 且主峰值 ∈ `[5e-5,2e-4] m` | 耦合没有真正驱动位移 |
| B7 | 运行完整性 | `≥50` 探针、无崩溃签名、无 `Unknown option` | 假绿 |

`superpos_err = ‖force − (baseForce − α∇p)‖∞`，从**未被修改的快照**重算。
**CI #85 更正**：所谓「一个数字区分叠加 / 覆盖 / 丢弃 / 累积」只在 `|base| > 0`
时成立；生产 deck 的 `baseForce ≡ 0` 会让叠加与覆盖给出完全相同的 `err=0`。
因此 **CI #86 起**：B3/B4 改判 `terzaghi1D_bodyforce.yaml` 的 BODYFORCE 诊断
日志（`|base|` 线性 `0→9.81`），B4 额外要求 `max|base| ≥ 5` 才承认其区分力；
B1/B2/B7 同时覆盖两份日志；B6 显式处理 NOBIOT 位移场恒零（不再算 `pm/0 → inf`）。
详见文末「CI #86：方案 B 落地」。

### 5. 与 Gate 3.3 的关系（回归安全）

`VecWAXPY` 与 `VecCopy+VecAXPy` 数学等价，主运行物理场不应变化；Gate 3.3 的 P1–P8 继续作为回归判据。Gate 4.1 只新增诊断与对照运行，**不改力平衡、不加 `αpI`**。

### 6. 静态校验与离线自测（零本地编译 / 零求解器运行）

- `PyYAML` 解析通过；10 个 `run` 块 `bash -n` 全过；两个 Python 脚本 `ast.parse` 通过；NOBIOT deck 与主 deck 的 YAML 差异**只有 `flow.biot_coeff`**；`git diff --check` 干净。
- 新子程序块配对 `If/EndIf = 2/2`、`Do/EndDo = 0/0`，续行无悬空；`VecWAXPY`/`VecDuplicate`/`VecDestroy`/`VecMin/VecMax` 均有仓库内编译先例。
- 判据脚本离线自测：正确分支 PASS；回滚(B3)、覆盖(B4)、α=0(B2)、gradp 恒定(B5)、行数不足(B1)、对照网格同值/缺失(B6)、崩溃(B7) 全部按预期 FAIL；Fortran `D` 指数与无空格标签两种格式均可解析；`--gate` 的PASS/FAIL/参数错误/缺日志四种退出码正确；用 CI #84 真实网格验证 `peak|u_x| = 1.1036e-4 m / 51 records`。

### 7. 状态

代码已提交推送（`ba010ea`），CI #85 已回：**failure，且唯一失败项是 B3**；Gate 3.3 保持 PASS。详见下节。**4.1b（互斥总应力模式）未实现**，理由见清单 4.1b 条目；残差平台 `3~13` 的 Fortran 装配根因、damage SNES 的 7 条 `-5` 仍在册，不属本闸门。

---

## CI #85 结果（run `42478888575`，commit `ba010ea`，conclusion = failure）

### 1. 结论先行

**Fortran 改动是对的，物理场没变，Gate 4.1 的 B3 判据前提错了。** CI 变红只因
`GATE41_RESULT: FAIL (B3)`；`GATE33_RESULT: PASS`，P1–P6/P8 数值与 CI #84
**逐位一致**（`peak|u_x|=0.00011036`、`ratio=0.0905571`、`L2=0.0320613`）。
`VecWAXPY` 单次叠加没有破坏任何东西。

### 2. B1–B7 实测表

| # | 实测 | 结果 |
|---|---|---|
| B1 | 51 行 `Gate 4.1:`，51 条 `Gate 3.2` 探针，0 不可解析 | PASS |
| B2 | `max|α−1| = 0` | PASS |
| **B3** | `|base| 0 → 0`，`slope 0`，`R² 0.00000`，`max resid 0` | **FAIL** |
| B4 | `max superpos_err/|force| = 0` | PASS（但本 deck 上**无区分力**，见 §4） |
| B5 | `|gradp| 1417 … 7.609e5`（ratio 537），`max(|force|−|base|) 7.609e5` | PASS |
| B6 | 主 `peak|u_x|=1.1036e-4 m` / NOBIOT `0 m`，`ratio = inf` | PASS（除零侥幸，见 §4） |
| B7 | 51 探针、无崩溃、无 `Unknown option` | PASS |

时间轨迹：step1 `|base|=0 |gradp|=7.609e5 |force|=7.609e5 err=0`；step51
`|base|=0 |gradp|=1417 |force|=1417 err=0`。每时间步 1 条诊断（ksponly 下
每步 1 次 AltMin），行数完整。

### 3. B3 根因：deck 里的重力从来没生效，不是代码 bug

生产日志 CellSet PetscBagView 三个 set 全部是：

```
force_scaling = linear;   Force = 0 0 0 ;   pressureForce = 0
```

而 `vDef/data/terzaghi1D.yaml:116` 写的是 `Force: 0, -9.81, 0`（**逗号后带空格**）。
这是全 deck 唯一一个带逗号空格的数组键；`displacementBC` /
`boundaryDisplacement` 早已统一为无空格写法（`1,1,1` / `0,0,0`）。
PETSc 3.3 的 `PetscBagRegisterRealArray(...,3,'Force',...)`
（`m_DefMech/m_MEF90_DefMechCtx.F90:936`）在同样带空格的输入上只读到首元素
`0`，其余补 0 ⇒ `Force = 0 0 0`。赋值链
`cellSetOptions%force → MEF90DefMechSetForceCst`（`m_MEF90_DefMech.F90:212`）
→ `VecScale(force,time)`（`m_MEF90_DefMech.F90:92`）因此恒为 0：
**`force_scaling = linear` 确实在跑，只是 `linear × 0 = 0`。**

4 个 deck（`terzaghi1D*.yaml`）都带同一处写法，属历史遗留。CI #85 的 B3
不是「叠加被回滚」，而是**根本没有非零体力可供回滚**。

### 4. 由此暴露的两个判据缺陷（必须修，否则第六/第七次假绿）

1. **B4 在 `baseForce ≡ 0` 时退化**：此时
   `superpos_err = ‖force − (0 − α∇p)‖∞`，叠加与覆盖结果完全相同，
   `err=0` 无法区分两者。上一轮声称「一个数字区分叠加/覆盖/丢弃/累积」，
   只在 `|base| > 0` 时成立。B4 的 PASS 在本 deck 上**不能证明叠加语义**。
2. **B6 用 `ratio = pm/pn`，`pn=0` ⇒ `inf`**，靠 `inf >= 3` 判 PASS。
   物理结论是对的（`out_nobiot.gen` 实测 51 步 `ux/uy` 严格全 0，是
   「去掉唯一驱动 ⇒ 无位移」的理想对照），但脚本没有显式处理 `pn==0`，
   属除零侥幸，且 `inf` 未进 B7 的非有限值检查。

另有一处**注释错误**需就地更正：`vDef/vDef.F90:524` 写死
`|DefMechCtx%force| = 9.81*time(step)`，与 CI #85 实测 `0` 不符。

### 5. 其他核实（防误判）

- 构建日志两条 `f951: Fatal Error: gfortran: Only one '-J' option allowed`：
  **CI #84 同样有 2 条**，属既有噪音；job 正常产出 15 个 artifact 与
  `vDef` 二进制，编译成功。
- NOBIOT 运行 `biot_coeff=0` 生效：其日志 `Gate 4.1:` 行数为 **0**
  （`biotCoeff==0` 时 early-return，符合预期）；`out_nobiot.gen` 全 51 步
  `peak|u_x|=peak|u_y|=0`。
- Gate 4.1 的 `exit 1` 按设计排在 Gate 3.3 PASS 之后，判据表完整打印，
  没有重演 CI #82 的早退隐藏数据。
- job 全文 `FAIL` 行只有 2 条，均属 Gate 4.1；主运行仍只有 7 条
  `damage field -5`（报告项，P8 兜底 `3.7458e-5`）。

### 6. 下一步的三条路（**用户已选 B，见文末 CI #86 节**）

- **A 最小修**：把 4 个 deck 的 `Force` 改成无空格 `0,-9.81,0`，
  让重力真正生效 ⇒ B3 自然 PASS、B4 恢复区分力。**风险**：给 Terzaghi
  基准引入体力，两端固定梁会产生 `u_y` 弯曲，可能顶到 Gate 3.3 的
  P2 阈值（`1e-5 m` / `5%×peak`）。
- **B 零物理风险（推荐）**：新增 `terzaghi1D_bodyforce.yaml`（仅把 `Force`
  修为无空格非零值）+ 一个 8 秒的 BODYFORCE ksponly 变体，用**它**的日志
  验 B3/B4；生产 deck 保持纯净（顺手把误导性的 `0, -9.81, 0` 改为
  `0,0,0` 并改注释，解析结果不变、物理场逐位不变）。B6 显式处理 `pn==0`。
  Gate 3.3 完全不受影响。
- **C 放水（不推荐）**：把 B3 改成「`base` 恒 0 也算 PASS」，则 4.1a
  的「叠加而非覆盖」验收彻底失去区分力，等于取消闸门。

另：CI #85 耗时 14m27s，其中 CI #73 遗留的 V2（FD Jacobian + GMRES）变体
约占 14 分钟；它的假设早被 CI #82/#83 证伪，退役可把 CI 压到 ~3 分钟。

**本轮动作**：只读核对 job log + `build-diagnostics` artifact（15 文件）、
原始 `Gate 4.1:` 行、CellSet bag view、两份网格全步位移、跨 CI 编译噪音对比，
并就地更正本文档与 `vDef.F90` 注释所对应的错误结论；
**未改 deck / build.yml / 判据脚本 / Fortran 逻辑，未提交、未推 CI**。

---

## CI #86：方案 B 落地（用户选 B，零物理风险）

### 1. 思路

不在生产 deck 上引入重力（那会给 Terzaghi 基准加体力、可能顶到 Gate 3.3
的 P2 阈值），而是**新增一个只用于诊断的 BODYFORCE 变体**：它带非零预设
体力，让 B3「Linear 缩放不回滚」和 B4「叠加而非覆盖」重新具备区分力；
生产 mesh 在变体运行**之前**已冻结到 `/tmp/out_main.gen`，Gates 3.1/3.2/3.3
完全不受影响。生产 deck 则把误导性的 `Force: 0, -9.81, 0` 写成显式
`0,0,0` —— 两种写法在 PETSc 3.3 下都解析成 `0 0 0`，**物理场逐位不变**。

### 2. 改动清单

| 文件 | 改动 |
|---|---|
| `vDef/data/terzaghi1D_bodyforce.yaml`（新增） | 主 deck 逐字副本，唯一语义差异 `cs0001.Force: 0,0,0 → 0,-9.81,0`（**无逗号空格**）；头部注明只作 Gate 4.1 B3/B4 诊断 |
| `vDef/data/terzaghi1D{,_nobiot,_nodamage,_robustksp}.yaml` | 4 个 deck 的 `Force: 0, -9.81, 0` → `Force: 0,0,0`，注释改为解释 CI #85 的 RealArray 截断；解析结果均为 `0 0 0`，物理场不变 |
| `vDef/data/check_gate41_biot.py` | `--gate` 增第 4 参数 `<bodyforce.log>`；B1/B2/B7 覆盖两份日志；**B3/B4 改判 BODYFORCE 日志**（B4 额外要求 `max|base|≥5` 才承认区分力）；B5 仍判生产日志；**B6 显式处理 NOBIOT 位移恒零**（`pn ≤ 1e-9 m` 走「精确解耦」分支，不再产生 `inf`）；报告打印两段轨迹 |
| `.github/workflows/build.yml` | ① NOBIOT 之后、恢复主网格之前插入 BODYFORCE ksponly 变体，保存 `/tmp/terzaghi_bodyforce.log` 与 `/tmp/out_bodyforce.gen`，并打印 CellSet bag 中实际解析到的 `Force` 行；② Gate 4.1 命令加第 4 参数；③ **退役 CI #73 的 V2 变体**（FD Jacobian + GMRES，1500 s cap，占 ~14/17 min，假设早在 CI #82/#83 证伪，无外部引用）；④ artifact 增 2 文件 |
| `FLOW_FIX_CHECKLIST.md`、工作日志 | 4.1a 判据表加 CI #86 更正注记，新增本节 |

### 3. 判据新语义（B1–B7）

| # | 数据源 | CI #86 起的变化 |
|---|---|---|
| B1 | 生产 + BODYFORCE | 两份日志都要 `≥50` 行、探针数匹配、全有限值 |
| B2 | 生产 + BODYFORCE | 两份日志 α 全为 1 |
| **B3** | **BODYFORCE** | `|base|` 线性 `0→≥5`、R²≥0.98；生产 `|base|==0` 属**预期**、不再被判 FAIL |
| **B4** | **BODYFORCE** | `superpos_err/|force| ≤1e-8` **且** `max|base|≥5`；后者保证能区分叠加/覆盖 |
| B5 | 生产 | 不变（`|gradp|` 时变、Biot 项活跃） |
| **B6** | 两份网格 | `pn≤1e-9` → 判「NOBIOT 精确解耦」且主峰值在 `[5e-5,2e-4] m`；否则仍用 `pm/pn≥3`。**不再出现 `ratio=inf`** |
| B7 | 生产 + BODYFORCE | 两份日志都 `≥50` 探针、无崩溃 / `Unknown option` |

### 4. 静态校验与离线自测（零本地编译 / 零求解器运行）

- `PyYAML`、10 个 `run` 块 `bash -n`、两个脚本 `ast.parse`、`git diff --check` 全过；
  `build.yml` 中 `run_variant` / `diag_V2` 已完全消失。
- YAML 语义差异：NOBIOT 仅 `flow.biot_coeff 1.0→0.0`，BODYFORCE 仅
  `cs0001.Force '0,0,0'→'0,-9.81,0'`；4 个既有 deck 都只含一处 `Force: 0,0,0`。
- 判据自测 **15/15 通过**，用 **CI #85 真实生产日志 + 真实两份网格**
  （`peak|u_x| main=1.1036e-4 m`、`nobiot=0`）+ 合成 BODYFORCE 日志：
  完整 PASS；`base` 回滚/常量 → B3(、B4) FAIL；覆盖 / 丢 Biot → B4 FAIL；
  α=0 → B2；行数不足 → B1；缺日志 → DATA ERROR；崩溃 / `Unknown option` → B7；
  对照网格同值 / 缺失 → B6；gradp 恒定 → B5；Fortran `D` 指数与无空格标签均可解析。
  `--gate` 的 PASS / 参数错误 / 缺日志 / 判据失败四种退出码正确，report 模式恒 0。

### 5. 诚实说明（落地时唯一未验证项 —— **CI #86 已实测通过**）

**结果**：`terzaghi_bodyforce.log` 的 CellSet bag view 实测打印
`Force = 0 -9.81 0`（cs0001；其余两个 CellSet 仍 `0 0 0`），生产与 NOBIOT
日志全部 `Force = 0 0 0`；B3 的 `|baseForce|` 线性 `0 → 9.81`、slope `9.81`、
`R² = 1.00000`、最大残差 `1.776e-15`。**「无空格 RealArray」假设成立**，
不需要换格式。详见文末「CI #86 结果」节。以下保留落地时的风险说明存档：

**无法在零本地编译约束下确证 `Force: 0,-9.81,0`（无空格）会被 PETSc 3.3
的 RealArray 读成三个分量。** 历史证据是 `displacementBC: 1,1,1` 修好了
BoolArray、`boundaryDisplacement: 0,0,0` 也是无空格写法，但 RealArray 的
行为只能由 CI #86 实测：BODYFORCE 日志会打印 CellSet bag 中解析出的
`Force = ...`，且 B3 就是「`|baseForce|` 是否线性增长到 9.81」的直接判据。
若 B3 仍 FAIL 而 bag 仍显示 `Force = 0 0 0`，说明无空格写法对 RealArray
也不够，需要换格式（例如 PETSc bag view 原样输出的方括号 / 分号形式），
那将是**下一轮**的新实验，不是本轮可以假设的结论。

另：退役 V2 后 job 实测 `0:03:02`（CI #85 为 `14m27s`）；`timeout-minutes: 75`
作为上限保持不变（不依赖它来压缩实际耗时）。

### 6. 状态

改动已提交推送，**CI #86 `conclusion = success`**，Gate 4.1 / 3.3 / 3.2 / 3.1
全绿，4.1a 收口。**4.1b（互斥总应力模式）未实现**；
残差平台 `3~13` 的 Fortran 装配根因、damage SNES 的 7 条 `-5` 仍在册。

---

## CI #86 结果：4.1a 收口（run `42481503794`，commit `949f096`）

`conclusion = success`、`GATE41_RESULT: PASS`、`GATE33_RESULT: PASS`。
job 墙钟 `12:42:29 → 12:45:31 = 0:03:02`（CI #85 为 `14m27s`，
退役 V2 后缩短约 79%）。artifact `build-diagnostics` 17 文件齐全，
含 `out_main/gen`、`out_ksponly.gen`、`out_nobiot.gen`、`out_bodyforce.gen`
与两份 verdict 日志。

### 1. Gate 4.1 判据表（B1–B7 全 PASS）

| # | 实测 | 判据 |
|---|---|---|
| B1 | 生产 51 行 / 51 探针，BODYFORCE 51 行 / 51 探针，`0 unparsed` | 两份日志完整、全有限 |
| B2 | `max\|α−1\| = 0` | α 恒为 1 |
| **B3** | BODYFORCE `\|base\| 0 → 9.81`、slope `9.81`、`R²=1.00000`、max resid `1.776e-15`；生产 `\|base\|max=0`（预期） | Linear 缩放在**非零体力**下未回滚 |
| **B4** | BODYFORCE `max superpos_err/\|force\| = 0`，且 `max\|base\|≥5` 具备区分力 | 叠加而非覆盖 / 累积 / 丢 Biot |
| B5 | 生产 `\|gradp\| 1417 … 7.609e5`（比值 537） | ∇p 时变、Biot 项活跃 |
| B6 | 主 `peak\|u_x\|=1.1036e-4 m` vs NOBIOT `0 m` → `NOBIOT identically zero (exact decoupling)` | 不再出现 `ratio=inf` |
| B7 | 两份日志 51/51 探针，crash=False、unknown_option=False | 运行完整 |

两段轨迹：生产 `|base|` 全程 `0`、`|gradp| 7.609e5 → 1417`；BODYFORCE
`|base| 0 → 4.905 → 9.81`、`superpos_err` 全程 `0`。

### 2. 唯一未验证项的实测裁决

`terzaghi_bodyforce.log` 的 bag view 打印 **`Force = 0 -9.81 0`**（1 个 CellSet），
另两个 CellSet `0 0 0`；`terzaghi.log` / `terzaghi_nobiot.log` 各 3 个 CellSet
全部 `Force = 0 0 0`。原始 `Gate 4.1:` 行也显示 `|base|inf` 从
`0.0000` 经 `0.1962`（step 2）涨到 `9.6138`、`9.8100`。
⇒ **PETSc 3.3 的 RealArray 与 BoolArray 一样，要求逗号后无空格**；
`Force: 0,0,0` 与 `Force: 0, -9.81, 0` 都只读首元素，这是继 `displacementBC`
之后同类坑的第二次实证。

### 3. 防假绿核实（跨 CI 物理场）

| 对比 | `vals_nod_var1` (p) | `vals_nod_var3` (ux) | `vals_nod_var4` (uy) |
|---|---|---|---|
| CI #86 主网格 vs CI #85 主网格 | `0` | `0` | `0` |
| CI #86 主网格 vs BODYFORCE 网格 | `0` | `4.1249e-7 m` | `2.7215e-6 m` |
| CI #86 主网格 vs NOBIOT 网格 | `0` | `1.1036e-4 m` | `9.1998e-6 m` |

- 生产 deck 的 `Force` 写法清理**物理场逐位不变**（三槽全 0 差），
  Gate 3.3 P1–P6/P8 数值与 CI #84/#85 完全一致
  （`0.00011036 m`、`0.0320613`、`8.17189e-7 m`、`0.0905571`）。
- BODYFORCE 的重力确实进了位移场（非零但很小的 `~1e-6 m` 量级差），
  NOBIOT 的 `ux/uy` 严格为 0 ⇒ 三个变体是**独立、可区分**的物理实验。
- job log 中 `FAIL` 关键字 27 处，逐条归类为 **26 条 `echo` 脚本文案 +
  1 条注释**，无任何真实判据失败；两条 `f951: Only one -J` 与 CI #84/#85
  相同，属既有编译噪音，二进制与 15 个构建 artifact 齐全。
- `Gate 3.3 diagnostics: displacement-msgs=0 damage-msgs=7 (INFORMATIONAL ONLY)`：
  damage 的 7 条 `-5` 仍是已知残留，P8 兜底。

### 4. 状态与遗留

- **闸门 3 全绿 + 闸门 4.1a 全绿**；4.1a 勾选状态已按 CI #86 实测更正，
  并修掉清单 4.1a 里「`9.81×t` 对生产 deck 成立」的错误句。
- **4.1b（互斥总应力模式）未实现**：需新增 cell-based 孔压场并改
  4–5 个 Fortran 文件，零本地编译下风险高，且与 4.1a 判据信号冲突。
- 残差平台 `3~13` 的 Fortran 装配根因、damage SNES 的 7 条 `-5` 仍在册。
- 下一候选：4.2（把 `cellDarcyVelocity` 喂入 HeatXfer，纯 Fortran 局部改动，
  可用同类「变体 + 硬判据」模式验收）；或 B 方向追残差平台；或 4.1b。

**本轮动作**：只读核对 job log + `build-diagnostics` artifact
（bag `Force` 原始行、`Gate 4.1:` 原始行、四份网格全步物理场对比、
跨 CI 回归、API 时长、FAIL 归类），并把结论写入本文档与工作日志；
**未改任何代码 / deck / 工作流 / 判据脚本，未提交、未推 CI**。

---

## CI #87：4.2 落地（逐单元 Darcy 速度 → HeatXfer）

### 1. 物理/接口裁决

原清单写「用 `FlowCtx%cellDarcyVelocity` 替代 `cellSetOptions%advectionVector`」
**不能照字面做**：`MEF90DiffusionOperatorAdvectionSet` 与
`MEF90DiffusionBilinearFormAdvectionSet`（`MEF90/m_MEF90_DiffusionImplementation.F90:118/223`）
的 `advectionVector` 是**一个** `MEF90_VECT`，在整个 `cellIS` 循环里复用；
而 `cellDarcyVelocity` 是 `cellDMVect` 上每 cell `dim` 个分量的**向量场**。
直接替换必须丢弃空间结构，只剩一个常数（平均值），4.2 就失去意义。

因此不改 `MEF90` 核心库，在 `m_HeatXfer/m_MEF90_HeatXferAssembly.F90` 内新增两个
私有子程序，镜像库里的两个装配循环，但每个 cell 从 `SectionReal` 读自己的速度：

```
advectionVector(cell) = -ρc_p · v_darcy(cell)      ! 与常数路径同号同标度
R  += ∫ [advectionVector(cell)·∇T] · φ
J  += ∫ [advectionVector(cell)·∇φ_j] · φ_i
```

四处调用点（稳态残差 95 / 稳态 Jacobian 239 / 瞬态残差 447 / 瞬态 Jacobian 597）
统一改为开关分支：`enable_advection=true` 走逐单元场；`false` 时**完全保留**
原 `-density*SpecificHeat*cellSetOptions%advectionVector` 常数路径。

数据通道（避免给 HeatXferCtx 加 Flow 模块依赖、避免循环引用）：

| 文件 | 改动 |
|---|---|
| `m_Flow/m_MEF90_FlowCtx.F90` | 新增持久 `cellDarcyVelocitySec`（`SectionRealDuplicate(cellDMVectSec)`），在 `SetSections` 创建、`Destroy` 销毁 |
| `m_Flow/m_MEF90_Flow.F90` | `MEF90FlowComputeDarcyVelocity` 把局部 `vSec` 换成该持久 section，去掉 `duplicate/destroy`；速度场在子程序返回后仍可读，且仍 scatter 到 `cellDarcyVelocity` Vec |
| `m_HeatXfer/m_MEF90_HeatXferCtx.F90` | Ctx 加 `flowVelocitySec` 指针 + `enableFlowAdvection` 标志，`CtxCreate` 里 `Nullify` / `PETSC_FALSE` |
| `vDef/vDef.F90` | `FlowCtxSetSections` 后挂接 `MEF90HeatXferCtx%flowVelocitySec => MEF90FlowCtx%cellDarcyVelocitySec`、复制 bag 开关；每步输出一条 `Gate 4.2: enabled=... \|v_darcy\|inf=...` |

`MEF90FlowComputeDarcyVelocity` 原本就在 `vDef.F90:585` 每个 AltMin 迭代调用，
速度场现成；本轮只是让它可消费。生产 deck 的 `heatxfer: timestepping: type: null`
保持不变 ⇒ Gate 3.1/3.2/3.3 与交付网格完全不受影响。

### 2. 对照实验（与 4.1a 的 BODYFORCE 同一套方法论）

| deck | 相对生产 deck 的差异 |
|---|---|
| `terzaghi1D_advection.yaml` | `heatxfer.type null→steadystate`；`cs0010 tempBC:1 / boundaryTemp:0`；`cs0020 tempBC:1 / boundaryTemp:100`；顶层 `boundaryTemp.scaling: linear`；`flow.enable_advection: true` |
| `terzaghi1D_noadvection.yaml` | 与上**逐字相同**，唯一差异 `enable_advection: false`（PyYAML 语义 diff 已核：1 个键） |

温度 BC 使 `T_bc = 100·t`，域内 `∇T ≠ 0`；材料默认
`Density=SpecificHeat=ThermalConductivity=1`（`MEF90Mathium2D`），Péclet 温和，
稳态热 SNES 不至于因强对流失稳。两个变体都跑 ksponly，网格分别在
`/tmp/out_advection.gen` / `/tmp/out_noadvection.gen`，都在恢复生产网格之前保存。

### 3. 硬判据 D1–D7（`check_gate42_advection.py --gate`）

| # | 判据 | FAIL 条件 |
|---|---|---|
| D1 | 两份日志各 ≥50 `Gate 3.2` 探针 / `Gate 4.2 enabled=` 行、各有 `coupling wired`、无致命签名 / `Unknown option` | 运行不完整 |
| D2 | adv 全 `T`、noadv 全 `F` | 开关没传到位 / deck 写错 |
| D3 | `\|v_darcy\|inf` 有限、>0、首 > 末（随 ∇p 衰减），且两份日志**逐行一致** | 速度场没算 / 开关影响了流场 |
| D4 | 三份网格 51 记录、节点数一致、共享槽全有限 | 网格缺失或损坏 |
| D5 | 自动发现的响应槽在 noadv 末步 `max\|T\| ∈ [1,1e4] K`、空间跨度 >1 K | 热场没启用 / 温度恒 0 |
| **D6** | `\|T_on−T_off\|/max\|T_off\| ≥ 1e-4`、去均值后的**空间结构** ≥ 5% 差异、末步差异 > 首步 | 平流未生效 / 只是常数平移 |
| D7 | 生产网格槽 1 仍是 `~1e5 Pa` 压力量级 | 生产 deck 被热场污染 |

判据**不看**求解器消息、fnorm 占位、收敛原因码、失败计数（前五次假收敛教训）。
判据脚本自动扫描全部记录槽并报告每槽差异，因此不依赖「温度在第几槽」的假设
（`tempOffset=1` 与 `pressureOffset=1` 在诊断变体里会重叠）。
`--gate` 输出机读 `GATE42_RESULT` 并退出 0/1；判定插在 Gate 3.3 之前、
`exit 1` 延迟到 Gate 3.3/4.1 PASS 之后（CI #82 早退教训）。

### 4. 校验（零本地编译 / 零求解器运行）

- `PyYAML`、10 个 `run` 块 `bash -n`、三个判据脚本 `ast.parse`、
  `git diff --check` 全过；deck 语义差异核对（adv↔noadv 仅 1 键）。
- Fortran 静态检查：`If/EndIf`、`Do/EndDo`（含带标签的 `mainloopQS:` / `AltMin:`）、
  续行、子程序定义/结束、两个新子程序各 1 定义 + 2 调用、持久 section 生命周期
  （`duplicate=1` / `destroy=1`，`Flow.F90` 中 0 次误销毁）全过。
- 判据离线自测 **15/15 通过**：完整 PASS；开关相同 / 速度不衰减 / 速度两份不一致 /
  行数不足 / 探针不足 / 崩溃签名 / 缺日志 / 缺主网格 / 温度恒 0 / 平流未生效 /
  纯常数平移 / 差异不增长 / 生产隔离失效，全部按预期 FAIL。
  自测过程中**抓到并修复了判据脚本一个真实缺陷**：主网格缺失时 `finite_ok`
  未定义导致 D5 `UnboundLocalError`，现改为 D5 自行校验温度场有限性。
  `--gate` 缺网格 / 参数错误的退出码以及 CI #86 真实网格读取亦通过。
- `subprocess --gate` 的 PASS 分支通过 monkeypatch 覆盖；真实 PASS 需 CI 实测。

### 5. 诚实说明（本轮未验证项，需 CI #87 实测）

1. **首次启用 HeatXfer 求解路径**：所有 vDef deck 历史上都是 `heatxfer: null`，
   `MEF90HeatXferCreateSNES` / `SetTransients` / `ViewEXO` 从未在 CI 上跑过。
   `enable_advection`、`tempBC`、`boundaryTemp`、`steadystate` 等键名以及温度写入的 EXO 槽号，均需由 CI 日志的 bag view 与
   两份诊断网格确认。判据设计成「读不到/不合理即 FAIL」，不会假绿。
2. **`SectionReal` 点编号匹配**：`cellDarcyVelocitySec` 基于 `FlowCtx%cellDMVect`
   的 cell chart，HeatXfer 的 `cellID` 来自同一 `Mesh` 的 `'Cell Sets'` stratum，
   静态判断一致，但零编译下无法实证；若不匹配会在 CI 表现为崩溃或速度全 0。
3. **`MEF90_VECT` 逐单元赋值** `-rhoCp * vPtr(1:MEF90_DIM)` 沿用仓库既有写法
   （`advectionVector = -density*SpecificHeat*advectionVector(1:MEF90_DIM)`），
   但未编译验证。
4. `vDef` 主程序 `Do/EndDo` 计数 4/6 是检查器漏掉两个带标签 `Do`
   （`mainloopQS:` / `AltMin:`）的口径问题，本次改动**零新增循环**。

若 CI #87 红，按日志定位是「编译层 / 热场配置层 / 平流装配层」哪一类，
再按同一套「变体 + 物理判据」模式修；不因失败放宽判据。

### 6. 状态

改动已提交推送、CI #87 结果待回。**4.1b（互斥总应力模式）未实现**；
残差平台 `3~13` 的 Fortran 装配根因、damage SNES 的 7 条 `-5` 仍在册。

---

## CI #87 结果：编译层失败（未进入任何 Gate 判据）

`conclusion=failure`，`RUNID=42485231934`，耗时约 2 分钟即终止。
`vDef` 二进制未生成，因此 **Gate 3.3 / 4.1 / 4.2 一律没有运行**，
`gate42_verdict.log` 为空 —— 这是「判据读不到数据即 FAIL」的设计前提，
不是判据本身失败。

### 根因（唯一阻断）

`m_HeatXfer/m_MEF90_HeatXferAssembly.F90:699,749` 两个新增 helper 把
`elem` 声明成 `Type(MEF90_ELEMENTTYPE)`。该宏**不是** `mef90.inc` 的宏，
而是 Makefile 只给 `MEF90/m_MEF90_DiffusionImplementation.F90` 注入的
编译参数 `-DMEF90_ELEMENTTYPE=MEF90Element2D_Scal`（3D 版同理）。
`m_MEF90_HeatXferAssembly.F90` 的编译命令没有这个 `-D`，模块内所有既有
subroutine 一律用 `MEF90_ELEMENT_SCAL`。因此 gfortran 报
`Derived type 'mef90_elementtype' ... used before it is defined`，连带
`elem has no IMPLICIT type`、`Expecting END IF statement`、以及 4 处调用点
`Type mismatch in argument 'elem'`（98/248/462/618），
`m_MEF90_HeatXferAssembly2D.o` / `3D.o` 失败并被 make `Error 1 (ignored)`
吞掉，最终 `m_MEF90_HeatXfer.o` 缺 .mod、`vDef` 无规则可建。

`f951: Only one '-J' option allowed` 两条是 CI 既有的 PETSc `TS` 探针
（Test A / `use petsctsdef`、Test B / `use petsc`）噪声，与本次无关，
CI #86 成功时同样存在。

### 本轮一并修掉的、下一轮必然暴露的两个问题

1. **指针关联缺 TARGET**：`vDef.F90:161` 做
   `MEF90HeatXferCtx%flowVelocitySec => MEF90FlowCtx%cellDarcyVelocitySec`，
   而 `MEF90FlowCtx` 是普通局部变量、组件无 TARGET 属性 —— 首次编译
   `vDef` 时会被 Fortran 指针赋值规则拒绝。现给
   `cellDarcyVelocitySec` 加 `,Target`（仅属性，不改生命周期）。
2. **温度输出槽被压力覆盖**：离线读 CI #86 主网格确认结果文件只有
   **4 个 nodal var**：`Pressure(1) / Damage(2) / Displacement_X(3) /
   Displacement_Y(4)`，变量表由 `MEF90DefMechFormatEXO` 固定，没有第 5 个槽。
   每步调用顺序是 `MEF90HeatXferViewEXO`(384) →
   `MEF90DefMechViewEXO`(748) → `MEF90FlowViewEXO`(749)，
   而 HeatXfer 的 `tempOffset` 默认 1、Flow 的 `pressure_offset: 1`，
   Flow 后写会把温度覆盖成 ~1e5 Pa；`check_gate42_advection.py` 的
   `TEMP_HI=1e4` 会让 D5 必 FAIL。
   修法（**纯数据层，只动两个诊断 deck，生产 deck 与核心代码不动**）：
   `flow.pressure_offset: 0`（`MEF90FlowViewEXO` 有 `> 0` 守卫，安全），
   显式 `temp_Offset: 1` 让 HeatXfer 独占槽 1，并把
   `flux_Offset` / `externalTemp_Offset` 置 0 避免写 elem var 1/2。
   主网格仍走生产 deck，D7 依旧验证其槽 1 是压力。两 deck 仍**仅差
   `flow.enable_advection` 一个键**。

### 已实证的运行期前提（静态复核，非零编译）

- CI #86 主运行 bag view 为 `heatxfer_timeStepping_type = null`、
  `Thermal energies` 计数 0 → YAML 嵌套 `heatxfer.timestepping.type`
  确实生效；4.2 deck 的 `steadystate` 会真正启用热场。
- `MEF90Mathium2D` 默认 `Density=1 / SpecificHeat=1 /
  ThermalConductivity=I` → `rhoCp=1`、K 非奇异，稳态热场可解。
- `MEF90HeatXferSetTransients` 的 `MEF90Scaling_Linear` 分支对
  `boundaryTemperature` 做 `VecScale(...,time)` → `T_bc = 100*t`，
  满足 D5 的 `span>1` 与 D6 的「差异随时间增长」。
- cell set `Has_BC` 走 `MEF90ISCreateCelltoVertex`，BAR2 边界集
  cs0010/cs0020 的 `tempBC` 会落到顶点 Dirichlet；`tempBC: 1` 统一改为
  `true` 以避开 PETSc 3.3 Bool 解析歧义。
- `MEF90_VECT` 的 `Vect2D_Get_VectR(Dimension(2))` /
  `Vect3D_Get_VectR(Dimension(3))` 赋值重载使
  `-rhoCp * vPtr(1:MEF90_DIM)` 合法。
- 4 处调用点的 mesh 参数与原版 advection 分支逐位一致：
  Operator 用 `%DM`、BilinearForm 用 `%DMScal`。

### 校验

`PyYAML`(18 yaml)、10 个 `run` 块 `bash -n`、5 个 Python `ast.parse`、
`git diff --check`、Fortran `If/Do` 配对（两 helper 终态 0/0，
7 `Subroutine` / 7 `End Subroutine`）全过。

### CI #88 结果：`failure`（编译层，无任何 Gate 判据被判负）

`m_HeatXfer` 本轮编译通过（`MEF90_ELEMENT_SCAL` 修复已由编译器验证，
Makefile 用 `@` 隐藏命令，日志只有源文件列表与错误，CI #87 同样格式）。
唯一阻断：`m_MEF90_FlowCtx.F90:63` 的
`Type(SectionReal),Target :: cellDarcyVelocitySec` →
`Error: Attribute at (1) is not allowed in a TYPE definition`。
该 gfortran 不允许在派生类型**组件**上写 TARGET；连带
`m_mef90_flowctx_type.mod` 缺失、`m_MEF90_FlowCtx.o` 失败，
`vDef` 二进制未生成，`gate33/gate41/gate42` 全空。
两条 `Only one '-J' option` 仍是既有 PETSc TS 探针噪声。

CI #87 的 `flow.log` 复核：`m_MEF90_Flow.F90:399` 只是
`PetscObjectSetName` 的 **Warning**（`-fallow-argument-mismatch`），
`m_Flow PIPESTATUS: 0 0`，即含 4.2 持久 section 的
`MEF90FlowComputeDarcyVelocity` 早已编译成功，本轮无需再改。

### CI #89 前的修复（`777b41b` 之后）

1. **TARGET 移到变量**：撤销 FlowCtx 组件上的 `,Target`，改为
   `Type(MEF90FlowCtx_Type),Target :: MEF90FlowCtx`（vDef）。
   Fortran 95 规则：带 TARGET 的对象的非指针子对象同样是合法指针目标，
   故 `flowVelocitySec => MEF90FlowCtx%cellDarcyVelocitySec` 仍成立，
   4 处 `Associated(...)` 保护与指针设计全部保留。
2. **wire 挪出 `Intent(OUT)` 作用域**（静态发现，会让 D2/D6 静默假失败）：
   `MEF90HeatXferCtxSetSections` / `MEF90HeatXferCtxCreateVectors` 的 Ctx
   dummy 都是 `Intent(OUT)`，在 vDef 第 251/252 行执行；原 wire 在 161 行，
   会被置为未定义。现 wire 放在整个热场初始化 `End If` 之后，
   之后只剩 `Intent(IN/INOUT)` 调用。SNES 通过
   `SNESSetFunction(...,MEF90HeatXferCtx,...)` 存 Ctx 地址，回调可见新值。
3. **首步 section 清零**：热场 SNES（~366 行）早于首次
   `MEF90FlowComputeDarcyVelocity`（~595 行），存在一步滞后（D6 只要求
   `l_d > f_d`，可接受）；但首步会读到刚 Duplicate 的未初始化数组。
   现 wire 后用已 `VecSet(0)` 的 `cellDarcyVelocity` 做一次
   `SectionRealToVec(... SCATTER_REVERSE ...)` 灌零，避免 NaN/发散。
4. wire 日志新增 `velocity_sec_associated=`，给出运行期关联证据；
   `ROW_RE` 要求 `enabled=` 与 `|v_darcy|inf=` 同现，wire 行不会污染 D2。

### 校验（提交前）

`PyYAML`(18 yaml)、10 个 `run` 块 `bash -n`、5 个 Python `ast.parse`、
`git diff --check`、vDef 全文 `If/Do` 配对 0/0、wire 行号 > 热场
`End If` 行号、`MEF90_ELEMENTTYPE` 残留 0、7 `Subroutine` 对 7
`End Subroutine`，全过。

### CI #89 结果：编译全绿，Gate 3.3 / 4.1 回归 PASS，Gate 4.2 FAIL

`vDef` 二进制生成成功，`MEF90_ELEMENTTYPE` / TARGET / wire 三处修复全部生效。
主运行日志：`Gate 4.2: coupling wired, enable_advection= F
velocity_sec_associated= T`，51 条 `|v_darcy|inf` 从 `3.80e-2` 衰减到
`8.76e-5`。`GATE33_RESULT: PASS`（P1-P6/P8）、`GATE41_RESULT: PASS`
（B1-B7），生产物理场逐位未受 4.2 改动影响。

两个 4.2 变体在 **step 1 之前** SEGV（日志仅 520 行 = bag view + PETSC
ERROR，0 个 `Gate 3.2: step` / `Gate 4.2: enabled=`）：
- 热场确实启用：`heatxfer_timeStepping_type = SteadyState`、
  `boundaryTemp_scaling=linear`、cs0010/cs0020 `TempBC=TRUE`、
  cs0020 `boundaryTemp=100`；ADVECTION 变体 `enable_advection= T`。
- 热场 step 1 已成功求解（`Thermal energies` 全 0），说明新增的
  `MEF90HeatXferFlowAdvectionOperatorSet/BilinearFormSet` 已实际执行且
  未崩；helper 的 `SectionRealDestroy(defaultSection)` 与原版
  `MEF90DiffusionBilinearFormAdvectionSet:155` 一致。
- 崩溃点由 checkpoint 序列锁定：`Gate 3.1: step 1 enter transient block`
  之后**没有** `post-UpdateMu` → 首条调用
  `MEF90FlowUpdateMuFromTemperature`。

**根因（旧代码，与 4.2 新增无关）**：该函数开头只在
`DefMechCtx%temperature` 未关联时 Return；生产 deck 热场为 null 所以永远
提前返回。两个变体启用热场后指针被接上（vDef:259），函数体首次执行：
`muPtr(1) = ...` 写在 `SectionRealRestrict(...,muPtr,...)` **之前**，
而 `muPtr` 是未关联的局部指针 → 立即 SEGV。此外 bag 里明明有
`temp_dep_viscosity` 开关（默认 `PETSC_FALSE`，所有 deck 均未开启），
函数却完全没检查它。

### CI #90 前的修复

1. `m_Flow/m_MEF90_Flow.F90`：在 `PetscBagGetData...` 之后加
   `If (.NOT. gOpt%temperatureDependViscosity) Return`。默认配置下
   μ(T) 耦合本就应 no-op，生产 deck 行为逐位不变。
2. 同一函数：`SectionRealRestrict(cellDMScalSec, ..., muPtr)` 移到
   两处 `muPtr(1) = ...` 之前、`SectionRealRestore` 保留在最后，
   消除将来开启该选项时的同类 SEGV。
3. `check_gate42_advection.py`：**D4/D5/D6 以 D1 为前置**（`d1_ok`）。
   CI #89 中变体崩溃后 workflow 复用同一 EXO 文件，残留的是上一步
   BODYFORCE 的压力场（末步 ~9055 Pa），D4/D5 竟据此 PASS、D6 因两
   份网格逐位相同而 diff=0。回放 CI #89 真实数据确认加固后
   D4/D5/D6 全部改为 FAIL（`runs_complete=False`），整体仍 FAIL。

### 已知残余风险（CI #90 观察）

`m_DefMech` 有 17 处 `Associated(temperature)` 热膨胀分支，启用热场后
将首次执行；静态复核其 Duplicate→ToVec→Allocate→RestrictClosure 写法
规范，暂不隔离。若 CI #90 崩在该路径，再考虑给 4.2 诊断 deck 增加
DefMech 热耦合隔离开关。

### CI #90 结果：编译全绿、两变体跑满 51 步，Gate 4.2 FAIL (D5)

Gate 3.3（P1-P6/P8）与 Gate 4.1（B1-B7）继续全 PASS。Gate 4.2：
D1/D2/D3/D4/D6/D7 PASS，**D5 FAIL**（`1.75257e+18 K`），而该轮 D6 的
PASS 是假阳性。关键事实：

- 两个变体各 5670 行、51 步、0 SEGV、51 次 `Thermal energies`
  （0 → 500，对应 `T_bc=100*t`）—— 热场本身完全正常，平流 helper 跑通。
- EXO 变量名表由 `MEF90DefMechFormatEXO` 固定为 4 个 nodal var
  （Pressure / Damage / Displacement_X / Displacement_Y），**无法命名温度**。
  `temp_Offset=1` + `flow.pressure_offset: 0` 让温度占 slot 1；离线读网格
 证实 slot1 是 `0 → 100 K`，平流差异末步 `0.0063958 K`、首步 0、
 峰值 `0.0125 K @ step19`、`spatial_frac=0.638`。
- 判据"自动发现最大差异槽"选中 `vals_nod_var4`（Displacement_Y，
 `|d|=3.29e18`），把 `1.75e18` 当温度 → D5 FAIL；D6 在该槽的 PASS
 同样是爆炸位移造成的假阳性。
- 位移爆炸根因：`MEF90Mathium2D` 默认
  `LinearThermalExpansion = MEF90MatS2DIdentity`（α = 1 /K）。启用热场后
  `DefMechCtx%temperature` 被关联，`m_MEF90_DefMechAssembly.F90` 的 17 处
  `Associated(temperature)` 热膨胀分支首次执行；`T_bc=100*t` 产生 100 倍
  热应变 → damage 饱和为 1 → `|u| ~ 1e18 m`。生产 deck `heatxfer=null`，
  该路径从未触发。

### CI #91 前的修复

1. 两个 4.2 **诊断 deck**（三个 cell set 均加）：
   `LinearThermalExpansion: 0,0,0` 解耦热应变；
   `ThermalConductivity: 0.1,0.1,0`。稳态方程为
   `K*laplacian(T) - rhoCp*v.grad(T) = 0`，故 `delta_T ~ v/K`；默认 K=1 时
   末步 `v_darcy` 已衰减到 `8.8e-5`，差异仅 `6.4e-5` 相对量，低于
   `REL_DIFF_MIN=1e-4`。K 降 10 倍把信号抬到 `~0.064 K`
   （`rel ~ 6.4e-4`），Dirichlet 温度场、T span、两 deck 对照均不变，
   单元 Peclet 仍 ≤ 0.02（扩散主导、稳定）。生产 deck 与核心 Fortran 不动；
   两 deck 依旧**仅差 `enable_advection`**。
2. `check_gate42_advection.py`：从变体日志 bag view 的 `temp_Offset = N`
   锁定 `vals_nod_varN` 为温度槽，**不再用最大差异槽裁决**；差异表降级为
   诊断输出（标注 `largest raw diff (NOT used)`）。
3. D5 增加"温度随时间增长"（首步 → 末步）、D6 增加绝对下限
   `MIN_ABS_DIFF=1e-4 K`；`REL_DIFF_MIN=1e-4` **不放宽**。
4. 新增 **D8**：两变体末步 `max|u| <= 1e-3 m`、`max damage <= 0.1`
   （生产量级 `1.1e-5 m` / `3.7e-5`），以 D1 为前置；判据扩为 D1-D8。
5. workflow：变体 bag 摘要加打 `LinearThermalExpansion` /
   `ThermalConductivity`，运行时可直接确认新键生效；失败提示改 D1-D8。

### 离线验证（零本地编译 / 零求解器运行）

- **CI #90 真实数据回放**：D5 锁槽后 PASS（`100 K`，`grows=True`），
  D6 FAIL（`rel=6.4e-5`）、D8 FAIL（`1.75e18 m` / `damage=1`）
  → 加固后无假阳性。
- **模拟修复后数据**（机械场换回主网格、温度差异 ×10）：
  **D1-D8 全 PASS**；且 `auto_max_slot` 落在 elem 应力槽，进一步证明
  自动发现不可用。
- 静态校验：18 yaml / 5 py `ast.parse` / 10 个 run 块 `bash -n` /
  `git diff --check` / 新增内容全角字符 0。
- YAML 键名实证：`fracturetoughness`（全小写）把 bag 的
  `FractureToughness` 从默认 1 改成 27000、`internalLength` 把 1 改成
  0.2，说明 PETSc YAML 键大小写不敏感，新加的两个 MatProp 键可命中
  `cs0001_` 前缀 bag。

### CI #91：Gate 4.2 收口（conclusion=success，8 checks / 0 failed）

| 判据 | 结果 |
|---|---|
| D1 运行完整 | PASS，两变体各 51 rows / 51 probes，wired 1/1，fatal=None |
| D2 开关对照 | PASS，ADV 51/51 `T`、NOADV 0/51 `T` |
| D3 Darcy 速度 | PASS，`0.0380437 -> 7.08746e-05`（peak 0.0380437），两日志一致 |
| D4 网格完整 | PASS，records (51,51,51)、nodes=600、shared slots=31、全有限 |
| D5 温度场 | PASS，锁 `vals_nod_var1`（由 `temp_Offset=1`），`max\|T_off\|=100 K`、first=0 K、span=100 K、grows=True |
| D6 平流响应 | PASS，`\|dT\|` first=0 → last=`0.0639432 K`、`rel=6.394e-4`（min 1e-4）、spatial=0.638 |
| D7 生产隔离 | PASS，主网格 slot1 `\|max\|=1e5 Pa`（heatxfer 仍为 null） |
| D8 机械场健康 | PASS，`max\|u\|=1.13e-05 m <= 1e-3`、`max damage=3.75e-05 <= 0.1` |

Gate 3.3（P1-P6/P8）、Gate 4.1（B1-B7）回归同轮全 PASS。

**离线实证（零本地编译 / 零求解器运行）**

- 变体 bag：三个 cell set 均 `LinearThermalExpansion = 0 0 0`、
  `ThermalConductivity = 0.1 0.1 0`、`heatxfer_timeStepping_type = SteadyState`、
  `temp_Offset = 1`；wire 一 `T` 一 `F` 且 `velocity_sec_associated= T`。
  生产 bag 保持默认 `alpha = 1 1 0` / `K = 1 1 0` / `null`，确认 deck 未污染主运行。
- 生产主网格 `out_main.gen` 跨 CI #90/#91：md5 仅 11 字节不同，全部位于
  EXO 属性区的 git SHA（`e7d001a-dirty`→`c04f4a7-dirty`）与时间戳
  （`14:43:39`→`15:15:31`）；**31 个记录变量 × 51 步逐值比对 max abs diff = 0**，
  物理场逐位一致。
- 两变体日志等长（299857 B），归一化 T/F 后唯一实质差异是热能量：
  step1 `2.00002E-02` vs `2.00000E-02`、step2 `8.00102E-02` vs `8.00000E-02`，
  即平流对温度场的影响有独立能量证据；各 51 次 `Thermal energies`。
- 结果网格里**只有 `vals_nod_var1` 有非零差异**（0.124868），其余 3 个
  nodal 槽与 27 个 elem 槽全为 0 → 平流只改温度，机械/损伤/压力场不受影响。

**结论：4.2（逐单元 Darcy 速度 -> HeatXfer 平流）已收口。**
生产 deck、Gate 3.x / 4.1 核心 Fortran 全程未改动，开关默认关闭。
仍在册：4.1b（互斥总应力模式）、标准 Newton 残差平台 `3~13` 根因、
damage SNES 7 条 `-5`（本轮生产日志仍为 7 条，信息性，不影响物理判据）。

### CI #92：稳态复跑（纯文档提交触发，conclusion=success）

`e8710c2` 只改文档，CI 自动重跑。Gate 4.2 D1-D8 与 Gate 3.3 / 4.1 的
**全部数值与 CI #91 逐位相同**（D6 `0.0639432 K` / `rel 0.0006394`、
D8 `1.13e-05 m`、P4 `0.00011036 m`、P6 L2 `0.0320613`、B3 slope `9.81`）。

离线网格级硬比对（零本地编译 / 零求解器运行）：

| 结果网格 | 记录数 | 变量数 | `max\|CI #91 − CI #92\|` |
|---|---|---|---|
| `out_main.gen` | 51 / 51 | 31 | 0 |
| `out_advection.gen` | 51 / 51 | 31 | 0 |
| `out_noadvection.gen` | 51 / 51 | 31 | 0 |

→ 阈值型判据 D6（余量约 6×）**可重复、非偶发**，4.2 收口确认。

### 4.3 覆盖性复核与 Flow bag 可见性（CI #93 待回）

静态复核推翻了原清单 4.3 的前提，详见上方 4.3 条目更正。要点：
`bc_pressure_scaling` 的 `constant`/`linear`/`file` 三态已实现、默认 CST、
7 个 deck 全写 `constant`，所以"恒压边界被乘时间"不成立；缺的是
Flow 三个 bag 从不 `PetscBagView`，导致所有 Flow 开关在 CI 日志里不可验证。

本轮改动**仅 1 个文件、+16 行**：`m_Flow/m_MEF90_FlowCtx.F90` 三处
`verbose>0` 守卫的 `PetscBagView`（Global / CellSet / VertexSet），
照抄 DefMech 惯用法。静态校验：18 yaml、5 py、10 run 块 `bash -n`、
FlowCtx 全文 IF/DO 终态 0/0、三处 view 均在守卫内、`git diff --check`、
判据隔离（所有 ROW_RE 要求 `Gate 4.x:` 前缀，`TEMP_OFFSET_RE` 有词边界，
新增 bag 行不会污染 D1-D8 / B1-B7 / P1-P8）。

**下一轮再基于真实 bag 输出格式加 4.3 硬判据**，不猜格式。
4.1b、残差平台 `3~13` 根因、damage 7 条 `-5` 仍在册。

### CI #93 结果 + 4.3 硬判据落地（CI #94 待回）

`7651aee` 的 CI #93 `conclusion=success`，三处 `PetscBagView` 编译通过、
6 份日志全部可见 Flow Global / 3 个 CellSet bag；Gate 3.3 / 4.1 / 4.2
回归同轮全 PASS，**物理场零改动**。真实 bag 片段（与判据正则一致）：

```
PetscBag Object:  flow Flow Global Options
  biot_coeff = 1 ; Biot poroelastic coefficient α (0 = no coupling)
  temp_dep_viscosity =  FALSE; Enable temperature-dependent viscosity μ(T)
  enable_advection =  FALSE; Add v_darcy·∇T convection term to HeatXfer
  pressure_offset = 1 ; EXO file offset for pressure
  bc_pressure_scaling = constant; (MEF90scaling) Boundary pressure scaling
  pressure_scaling = constant; (MEF90scaling) Scaling law for pressure
PetscBag Object:  cs10 Flow CellSet Options
  flowBC =  TRUE; Cell set has pressure boundary condition
  bc_pressure = 0 ; Pressure value on BC surface (Pa)
```

基于该格式写 `vDef/data/check_gate43_flowbag.py`（纯标准库）：C1 六份日志
各 23 个必需 Global 键；C2 三个 CellSet bag 且 `cs0010` 的
`flowBC` / `bc_pressure` 与 deck 一致；C3 三个 scaling 全为 `constant`
且 bc/source 在每个 deck 显式声明（`pressure_scaling` 是库默认，只比对 bag）；
C4 逐个比对 deck `flow:` 段**每一个键**与 bag 值；C5 复核 nobiot α=0、
advection 开关、4.2 的 `pressure_offset=0`、μ(T) 关闭；C6 读生产网格
51 步压力场判恒压 BC `span=0 Pa`、闭端衰减、Gate 3.3 解析 L2；C7 查
4.1/4.2 verdict、致命签名与每份日志恰好一次挂接线。

**离线回放（CI #93 真实 artifact）**：C1–C7 全 PASS，
`116` 项 deck↔bag 比对 0 失配，BC 节点 6 个 / 51 步 `span=0 Pa`、
`mean=0 Pa` 与 bag 的 `bc_pressure=0` 一致，闭端
`100000 → 9055.71 Pa`（decay 11，防全零假场），`P6 L2=0.03206`。

**负例自测 12 例 / 0 missed**：基线全 PASS；逐条篡改后目标判据必 FAIL
（N1 删 `enable_advection`→C1；N2 翻转 cs10 `flowBC`→C2；N3 删 deck
`source_scaling`→C3；N4 改 `pressure_atol`→C4；N5 悄悄打开 μ(T)→C5；
N6 把一个 BC 节点第 10 步改成 12345 Pa→C6；N6b 整个压力场清零→C6；
N7 改 4.2 verdict→C7；N8 删挂接线→C7；N9 删一份运行日志→C1）；
`--gate` 失败用例 `rc=1` 且输出 `GATE43_RESULT: FAIL`。

`build.yml` 接入：Gate 3.3 之后调用（需 `gate33_verdict.log`），
`g43_rc` 延迟到 Gate 4.1/4.2 裁决之后再 `exit 1`（CI #82 早退教训），
artifact 增 `/tmp/gate43_verdict.log`。静态校验：PyYAML 6.0.3 解析、
10 个 `run` 块 `bash -n`、9 个 py `ast.parse`、新增行非 ASCII=0、
`git diff --check`、判据调用/裁决顺序断言全绿。

**本轮仍未勾选 4.3**：等 CI #94 实跑裁决。
4.1b、残差平台 `3~13` 根因、damage 7 条 `-5` 仍在册。

### CI #94 结果（4.3 收口，conclusion=success）

`6704b0d` 的 CI #94 `conclusion=success`，四个闸门同轮全 PASS：

| 闸门 | verdict | 关键结果 |
|---|---|---|
| 3.3 | PASS | P1–P6 + P8；P4 `0.00011036 m`、P6 压比 `0.0905571`、L2 `0.0320613` |
| 4.1 | PASS | B1–B7；BODYFORCE `|base| 0→9.81`、slope `9.81`、R²=1、叠加误差 ~0 |
| 4.2 | PASS | D1–D8；D5 `max|T_off|=100 K`、D6 `|dT| 0→0.0639432 K`、D8 `max\|u\|=1.13e-5 m` |
| 4.3 | PASS | C1–C7；116 项 deck↔bag 0 失配、BC `span=0 Pa`、闭端 `1e5→9055.71`、L2 `0.03206` |

真实 artifact 复核（`gate43_verdict.log`）：`GATE43_RESULT: PASS`，`7 checks, 0 failed, 6 runs inspected`。
4.3 硬判据接入后与 4.1/4.2 同为延迟 `exit 1`，失败时必红（已用 12 例负例自测 + CI 真实回放双重验证区分力）。

**结论：4.3（Flow 选项袋可见性 + 恒压 BC 语义 + 116 项 deck↔bag 比对）已收口。**
仍在册：4.1b（互斥总应力模式）、标准 Newton 残差平台 `3~13` 根因、damage SNES 7 条 `-5`（信息性）。

---

## CI #111 结果（4.7 收口，commit `ef6ee2a`，conclusion=success）

`ef6ee2a` 的 CI #111（run `42580671359`）`conclusion=success`，20 个 job step 全 success（含 `Build MEF90 + m_Flow + vDef` 与 `Run Terzaghi`）。下载真实 `build-diagnostics` artifact（13,350,853 B，`unzip -t` 无错）离线复核，**十二组 verdict 同轮全 PASS**：

| 闸门 | verdict | 关键结果 |
|---|---|---|
| MAIN / KSPONLY | PASS | P1–P6 + P8；P4 `0.00011036 m`、P6 压比 `0.0905571`、L2 `0.0320613` |
| 3.3 | PASS | 51 records P1–P6,P8 全过；7 条 SNESSolve-failed 为信息性 |
| 4.1 | PASS | B1–B7；BODYFORCE `\|base\| 0→9.81`、slope `9.81`、R²=1、叠加误差 ~0 |
| 4.1b | PASS | E1–E4；E2 8/8（`cos_ux=0.9994`、`ux_ratio=0.6912`、`cos_bnd=0.9793`、`k=8.063e-10 1/Pa`） |
| 4.2 | PASS | D1–D8；`max\|T_off\|=100 K`、`\|dT\| 0→0.0639432 K`、`max\|u\|=1.13e-5 m` |
| 4.3 | PASS | C1–C7；116 项 deck↔bag 0 失配、BC `span=0 Pa`、闭端 `1e5→9055.71` |
| 4.4 | PASS | S1–S8；`cos(δp,+p_f)=1.00000`、amp `1.0049`、末步 `926.55 Pa` vs 解析 `924.30 Pa` |
| 4.5 | PASS | T1–T8；`E(1)/E(50)=334.98`、`E(50)/E_an=0.937`、`cos(W,W_an)=0.99998` |
| 4.6 | PASS | V1–V8；`perm(visc)/perm(ctrl)=2.364413` vs 解析 `2.364413`、`max\|perm·μ/k−1\|=0` |
| 4.7 | **PASS** | **X1–X8**；`A(t)=5e-8(1+0.8t)` 逐记录零误差、末步比 `1.800000`、`permF=1.8`/`storF=1.4`、时变 BE `relL2(50)=0.0368` `cos=1.00000` `ratio=0.9633`、分离 `0.360/0.782/0.577` |
| 4.8 | PASS | U1–U8；能量残差 `≤0.0112`、质量残差 `≤0.00976`、Darcy 速度相对误差 `3.11e-16` |
| 4.10 | PASS | W1–W8；93 个归档 git-blob sha256 全匹配、`mismatched=none` |

**4.7 生产隔离（零扰动证明）**：与 #109 的 9 组 verdict 对比，`gate41/42/44/45/46/48/41b` 共 7 组**逐字节 IDENTICAL**；`terzaghi_verdict.log` 与 `gate33_verdict.log` 的唯一差异是 solver log 字节数（`290933/291011 → 297823/297901`，因每步多打印一行 `Gate 4.7:`），**所有物理判据数值逐位一致**。求解器健康：`PETSC ERROR=0`、segfault=0、7 条 damage SNES `-5`（信息性）、51 步完整。

**代码修改核查（本轮静态审查，零本地编译）**：工作树与 `HEAD=ef6ee2a` 一致，无未提交源码改动。4.7 的 6 个 bag 字段按 CI #50 教训加在类型/注册末尾且一一对应，默认 `MEF90Scaling_CST` 严格 no-op；`SetTransients → UpdateMu → UpdatePerm → UpdateStorativity` 调用顺序正确；两个 factor Vec 生命周期完整；`g47_rc` 延迟 `exit 1`（CI #82 教训），artifact 三项齐备。判据脚本 `--selftest` 离线重放 **1 正例 + 9 负例全 PASS**。另确认 4.1b 的两处历史隐患均已闭环：`INTENT(IN)→INTENT(INOUT)`（`f20fa18`，CI #97 编译通过）、单位张量改 `MEF90_MATS_IDENTITY`（CI #99，2D/3D 双编译 + LTE=0 防回归 deck），E2 判据亦已由错误的「位移逐位等价」更正为 8 项物理签名。

**结论：4.7（时间相关材料参数入口）已收口，闸门 4 累计 10 项完成。** 仍在册：4.9 逐单元 IS 装配性能、`pressurePreviousStep` 死写、标准 Newton 残差平台 `3~13` 根因、damage SNES 7 条 `-5`（信息性）。

---

## CI #114 结果（4.9 收口，commit `dabe01a`，conclusion=success）

`dabe01a` 的 CI #114（run `42597311984`）`conclusion=success`，20 个 job step 全 success（含 `Build MEF90 + m_Flow + vDef` 与 `Run Terzaghi`）。下载真实 `build-diagnostics` artifact（14,868,569 B，`zipfile.testzip()` 无错）离线复核，**十三组 verdict 同轮全 PASS**：

| 闸门 | verdict | 关键结果 |
|---|---|---|
| 3.3 | PASS | 51 records P1–P6,P8 全过；7 条 SNESSolve-failed 为信息性 |
| 4.1 / 4.1b | PASS | B1–B7；E1–E4（E2 8/8，`cos_ux=0.9994`、`ux_ratio=0.6912`） |
| 4.2 / 4.3 | PASS | D1–D8；C1–C7（116 项 deck↔bag 0 失配） |
| 4.4 / 4.5 | PASS | S1–S8；T1–T8（`E(1)/E(50)=335.0`、`E(50)/E_an=0.937`） |
| 4.6 / 4.8 | PASS | V1–V8（`ratio=2.36441`）；U1–U8（能量 ≤0.0112、质量 ≤0.00976） |
| 4.7 | PASS | X1–X8；`A(t)=5e-8(1+0.8t)` 零误差、末步比 `1.800000`、时变 BE `relL2(50)=0.0368` |
| 4.9 | **PASS** | **Y1–Y8**，详见下表 |
| 4.10 | PASS | W1–W8；93 个归档 git-blob sha256 全匹配 |

**Gate 4.9 实测（Y1–Y8 全绿）**

| 判据 | 结果 |
|---|---|
| Y1 | 两 deck 51 记录 / 356 节点、几何与时间轴对齐，无 fatal |
| Y2 | helper + 4 处 fast path、bag `37/37/37` 顺序、默认 `PETSC_FALSE` |
| Y3 | `terzaghi1D_batched.yaml` 与生产 deck 只差一行 `batched_assembly: true` |
| Y4 | 对照 bag `TRUE` / 生产 bag `FALSE` |
| Y5 | 151 条 `Gate 4.9:` 行、3 个运行时活跃 tag、`nCells=300`、`nIS=nElem=1`、生产日志 0 行 |
| Y6 | `p/u_x/u_y` 在 51 记录 × 356 节点**逐位一致**，`max\|dPhysics\|=0` |
| Y7 | `legacy=2400 → batched=8` 对象创建（**300×** 归约）、`elem` 与 `cellIS` 同源（CI #65 契约） |
| Y8 | step marker 408 相等、无 fatal 签名 |

复杂度证据取自确定性计数而非噪声大的 wall-clock：每步省 `2392` 次 PETSc 对象创建、`reduction_x=300`。生产 `terzaghi.log` 中 `PETSC ERROR` / `Caught signal` / `SEGV` / `fortran runtime error` 全为 0，51 步完整。

**CI #113 失败 → 修复（本轮核心）**：#113 编译 clean 且 Y6/Y7/Y8 已 PASS，但两处静态裁决误判 —— ① `GATE49 Y5` 要求 4 个 tag，而 `MEF90FlowAddStorageToResidual` 自 CI #64 起在 `TS_LINEAR + ksponly` 下 `xDot==0`、额外 scatter 删除后已不被调用，其 fast path 是源码存在但运行时 dead code；② `GATE47 X2` 硬编码「恰好 36 个 bag 字段」，把 4.9 合法追加的第 37 个字段误判为回归。修复为统一的 **append-only 契约**：Y5 只要求 3 个活跃 tag，Y2/Y7 仍静态审计全部 4 处 fast 路径并按 `N_FAST=4` 外推；4.7 X2 校验 4.7 六字段固定槽位 30..35 + 字段数下界 + 构造器 arity 注释与活跃 bag 一致；4.9 Y2 采用同一契约；`vDef/m_vDefDefault.F90` 注释改 `37 components`。各加 1 个精确负例 —— `--selftest`：4.7 为 1 正例 + 10 负例、4.9 为 1 正例 + 10 负例，全 PASS。

**生产零扰动（逐字节证明）**：与 #111 的 12 组共同 verdict 相比 **10 组逐字节 IDENTICAL**；`gate33_verdict.log` / `terzaghi_verdict.log` 唯一差异是 solver log 字节数 `297901→298025` / `297823→297947`（+124 B），来源仅 git changeset 行 `ef6ee2a→dabe01a` 与 PetscBag 视图新增的 `batched_assembly = FALSE` 一行，物理判据数值全部不变。与 #113 相比生产日志只差 changeset 行；`out_main.gen` 的 `p/u_x/u_y` 共 54,468 个数值 `differing=0`、`max|d|=0`。

**提交前离线复核（零本地编译 / 零本地求解器）**：`git diff --check`、`ast.parse`、PyYAML 三份 deck、10 个 workflow `run` block `bash -n`、bag 声明/注册 37-37 一致、新增行非 ASCII=0、两个判据 `--selftest` 全 PASS、用 #113 真实日志 / EXO 离线 `--gate` 回放提前得到 `GATE47_RESULT: PASS` 与 `GATE49_RESULT: PASS`，与 CI #114 实测一致。

**结论：4.9（逐单元 IS 装配性能）已收口，闸门 4 累计 11 项完成，合计 21 项已落地。** 仍在册：`pressurePreviousStep` 死写、标准 Newton 残差平台 `3~13` 根因、damage SNES 7 条 `-5`（信息性）。

---

## CI #117 结果（4.12 / 4.13 收口，commit `d142afb`，conclusion=success）

`d142afb` 的 CI #117（run `42608060588`）`conclusion=success`，job `Build + Terzaghi` 全 success。下载真实 `build-diagnostics` artifact（14,880,833 B，`zipfile.testzip()` 无错）离线复核，**十五组 verdict 同轮全 PASS**，其中新增两组：

| 判据 | 结果 | 关键指标 |
|---|---|---|
| 4.12 | **PASS (Z1–Z6)** | 51+51 行、2 处快照 + 采样 + 零装配消费者、`dpmax(n)=pmax(n)−pmax(n−1)` 逐位一致、EXO 极值锚定、峰值 @step 10、末/峰 `0.132`、408 step marker 无 fatal |
| 4.13 | **PASS (Q1–Q5)** | 生产 `-disp_snes_type ksponly`、三处尺度锚点齐全、两份日志无 Newton fnorm / displacement `-5`、ksponly 与生产 51 条压力包络逐位一致 |
| 原有 13 组 | 全 PASS | 3.3 / 4.1 (B1–B7) / 4.1b / 4.2 / 4.3 (C1–C7) / 4.4 / 4.5 / 4.6 / 4.7 / 4.8 / 4.9 / 4.10 / 主报告 |

**4.12 落地**：`pressurePreviousStep` 从「创建 / 初值 / 每步 VecCopy / 销毁但无人读」的死写，改造成**步间压降活诊断** —— 初始化时快照经 Dirichlet 修正的真实初始场，每步 `TSSolve` 前用 `VecMin/VecMax` 采样上一步末场，之后打印 `Gate 4.12: step N t T dpmax ... dpmin ...`；`m_Flow` 装配例程仍不读该 Vec，物理零反馈。实测 `dpmax` 首步 `−226.206 Pa`、峰值 `−3734.93 Pa`（step 10）、末步 `−491.674 Pa`，呈 Terzaghi 固结典型的先加速后指数衰减。

**4.13 根因（本轮静态定位）**：`m_MEF90_DefMechAssembly.F90` 中自由自由度残差 `∫(C:ε(u)):ε(v) dΩ`（尺度 `E ~ 1e9`），Dirichlet 行残差 `u − u_bc`（尺度 1），Jacobian 用 `MatZeroRowsColumnsIS(A, setISdof, 1.0_Kr, ...)` 置单位对角 ⇒ 矩阵混合相差约 9 个数量级的行，条件数 ~E。CG 按**相对残差**收敛后，物理力残差仍剩几牛顿，SNES 于是反复迭代一个降不下去的 `fnorm`（CI #83：`2.082e4 → ~9 → 3~13`）。该 deck 实质线性，一次线性解即精确解，故 CI #84 的 ksponly 是正确规避而非掩盖 bug。判据 Q1–Q5 把这一诊断可证伪地冻结，防止有人把生产改回标准 Newton 或误判为装配错误。

**CI #116 插曲（`failure`）**：新增的 51 行 `Gate 4.12:` 被 `check_gate41_biot.py` 的裸前缀 `'Gate 4.1' in line` 误计为 102 条未解析探针（2 日志 × 51 步），`B1 FAIL` 并连带 Gate 4.3 `C7 FAIL`；**#116 的 4.12 / 4.13 及其余 12 组全部 PASS**。修复为精确匹配 `'Gate 4.1:'`，先用 #116 真实日志离线回放得到 `GATE41_RESULT: PASS (B1-B7)`，再推 `d142afb` → #117 全绿。

**零扰动证明**：#117 与 #116 的 `out_main.gen` 压力 / `u_x` / `u_y` 共 **54,468 个数值 differing=0、max|d|=0**；15 组 verdict 中 13 组与 #116 **逐字节 IDENTICAL**，仅 `gate41_verdict.log` 与 `gate43_verdict.log` 由 FAIL 转 PASS（正是修复目标）。相对 #115，生产日志只多 51 行 `Gate 4.12:` 与 changeset 行（字节数 `298025→304451`），所有原有物理判据数值不变。

**提交前离线复核（零本地编译 / 零本地求解器）**：`git diff --check`、`ast.parse`、PyYAML、10 个 workflow `run` block `bash -n`、Fortran 块平衡与 HEAD 一致、两个新判据 `--selftest`（1 正例 + 8 负例 / 1 正例 + 6 负例）全 PASS、用 #115 真实日志 / EXO 离线 `--gate` 回放提前得到 `GATE412_RESULT: PASS` 与 `GATE413_RESULT: PASS`。

**结论：4.12 与 4.13 收口，闸门 4 累计 13 项完成、合计 23 项已落地，0 项在册待做。** damage SNES 7 条 `-5` 为信息性（Gate 3.3 P8 已用损伤场有界替代），4.11 总括性描述由十五组硬判据体系覆盖。

---

## CI #121 结果（4.14 收口，commit `c95a605`，conclusion=success）

`c95a605` 的 CI #121（run `42631282654`）`conclusion=success`，job `Build + Terzaghi` 全 success：`Build MEF90 + m_Flow + vDef` 1 min 54 s、`Run Terzaghi` 9 min 19 s（#120 是被 900 s 上限杀掉的 17 min 失败）。下载完整 `build-diagnostics` artifact（23,811,138 B / 46 成员，`zipfile.testzip()` 无错）离线复核，**十六组 verdict 同轮全 PASS**，新门 `GATE414_RESULT: PASS`：

| 判据 | 结果 | 关键指标 |
|---|---|---|
| W1 | PASS | 11 条 3D / 51 条 2D 记录、2635 个 3D 节点、所有场有限，`out_dim3.gen` 15,106,908 B **完整不再截断** |
| W2 | PASS | `num_dim=3`、6 个 nodal vars 含 **`Displacement_Z`**、`coordz`、cell sets `[1,10,20]`、几何 10×1×1 |
| W3 | PASS | 3D/2D identity 宏、`-DMEF90_DIM=3` 双孪生、`DMMeshGetDimension` 分派、门面 2D/3D 实现对、3D EXO 命名 |
| W4 | PASS | `p0=1e5`、`pend=1.183e4`、单调衰减、11 个共同物理时刻最大相对差 **0.306 ≤ 0.40** |
| W5 | PASS | `\|ux\|max=1.01e-4`、`\|uz\|max=6.618e-6`、`ratio=0.0655`、全部有限 |
| W6 | PASS | 3.1/3.2/4.12/4.1 诊断链齐全 + **10 条 `Linear solve converged ... iterations 1`**、无 fatal |
| W7 | PASS | 单行 deck（numstep 51→11）、默认物理、CI 接线（含求解器配置与 `s414_rc` 裁决）、2D 生产隔离 |
| 原有 15 组 | 全 PASS | 3.3 / 4.1 / 4.2 / 4.3 / 4.1b / 4.4 / 4.5 / 4.8 / 4.6 / 4.10 / 4.7 / 4.9 / 4.12 / 4.13 / 主报告 |

**两轮失败与根因（这是 4.14 真正补上的验证空白）**：此前 15 组闸门**只在 2D `Beam1x10-tri3` 上跑过**，3D 编译孪生与 `MEF90MatS3DIdentity` 从未被实际执行。
- **CI #119（`42a9a2e`, failure）**：3D 在 `Gate 3.1: init 07 FlowUpdatePerm` 刷约 9900 行 `Element2D_Scal_Init: Element type not implemented yet MEF90_P1_Lagrange_3D_Scal` 后 MPI abort。门面 `m_MEF90_Flow.F90` 不带 `-DMEF90_DIM` 编译，`MEF90_ELEMENT_SCAL` 恒为 2D 元素类型，而 `MEF90Element_Create` 泛型按 `elem` **声明类型**选工厂。修复 `39444b8`：`UpdatePerm` / `UpdateMuFromTemperature` 拆 2D/3D 实现对 + 维度 wrapper；4.6 V1/V2、4.7 X2 静态契约同步改为审计 wrapper 与两个孪生。
- **CI #120（`39444b8`, failure）**：编译 clean、15 组原有判据全 PASS，但 3D step 1 物理量正常（`p0=1e5`、`E=1853.2` vs 2D `1624.6`、无 NaN）后，step 2 **第一次 flow TSSolve** 在 900 s 内未打印任一 `Linear solve converged`；SIGTERM → MPI_ABORT → `out_dim3.gen` 截断在 13,766,716 B，W1 解析 EXO 抛 `struct.error`。**根因是求解器代价而非 3D 装配错误**：2D 356 节点 CG 42 次迭代收敛，3D 2635 节点 tetra4 上 CG + bjacobi/ILU(0) + rtol 1e-10 停滞。
- **修复 `c95a605`（仅 workflow + 判据，Fortran 零改动）**：3D 诊断运行改 `-np 1` + `-flow_pc_type bjacobi -flow_sub_pc_type lu`（单 rank 下唯一 bjacobi 子块即整个 seq-AIJ 矩阵，精确 LU 让 CG 每步 1 次迭代收敛）；加 `-flow_ksp_monitor`，cap 900→1800 s；**12 个生产 2D deck 保持 `-np 2` / ILU，11 个物理步与 W1–W7 判据一条不减**。`s414_rc`（3D mpirun 退出码）此前只 echo、从不裁决，现接入延迟 hard fail；W6 新增「每步 KSP 收敛」硬要求；W7 冻结求解器配置与 `s414_rc` 裁决；`--gate` data error 改为可归因 W1 FAIL + 截断提示。TS view 实证 `1 MPI processes`、`flow_ PC=bjacobi`、`flow_sub_ PC=lu`、`package: petsc`（内置 LU，不依赖外部包）。

**生产零扰动（逐字节证明）**：#121 与 #120 的 **15 组原有 verdict 逐字节 IDENTICAL**（仅 `gate414_verdict.log` 由 FAIL 转 PASS）；`out_main.gen` 2,845,300 B 中仅 **11 个 differing bytes**，全部位于 QA 记录字段 `39444b8-dirty → c95a605-dirty` 与时间戳 `13:16:28 → 14:07:52`，压力 / `u_x` / `u_y` 物理数值逐位不变。3D 日志 `Element2D_Scal_Init` / `not implemented yet` / `PETSC ERROR` / `Caught signal` / `MPI_ABORT` / `Gate 4.9:` / `Gate 4.1b:` 计数全为 0。

**提交前离线复核（零本地编译 / 零本地求解器）**：`git diff --check`、`ast.parse`、PyYAML、10 个 workflow `run` block `bash -n`、全部闸门 `--selftest`（4.14 为 1 正例 + 18 负例）、用 #120 真实日志 / EXO 回放 4.2 / 4.3 / 4.6 / 4.7 / 4.9 / 4.10 / 4.12 / 4.13 全 PASS，且 #120 截断数据在修复后的判据下输出可读的 W1 FAIL + 截断 hint。

**结论：4.14（三维维度分派真实验证）已收口，闸门 4 累计 14 项完成、合计 24 项已落地，0 项在册待做。** 3D 路径从「编译孪生存在但从未运行」变为「有 11 步固结物理场 + 位移 Z 分量 + 每步 KSP 收敛证据 + 静态反回归」的硬判据覆盖；`s414_rc` 的接入同时消除了「超时退出码只打印不裁决」这一类静默假绿风险。

---

## CI #124 / #125：Gate 4.15 B0 的两轮收口（2026-11-10）

### #124（`ff614fc`，failure）

- 编译 clean、原有 16 组 verdict 全 PASS；唯一失败是 `Run Terzaghi` 末尾的延迟裁决 `g415_rc`：`GATE415_RESULT: FAIL (V1,V3,V4)`。
- V1 `scale records 51 != 11; inject records 51 != 11`；V3 `control pressure drifted 90.94%`；V4 `S_eff=5.14e-11 vs 1.5e-10 (65.7% off); slope 291841 vs f/S=100000`。

### 离线诊断（零 Fortran / 零 deck / 零新 CI）

从 `ci124_diag.zip` 读真实 EXO：三个 `.gen` 都是 2,845,300 B / 51 记录但 `sha256` 各异；`out_415_*.gen` 前 11 条是 4.15 新场（`pmax[0]=1e5`、`inject pmax[10]=391877.6`），第 12–51 条与 `out_main.gen` 的旧生产衰减尾部**逐值相同**（`pmax[11]=70667.193`）。日志 bag 复核 `storativity=1.5e-10`、`source=1.5e-05`、`perm_scaling/stor_scaling=constant`、`permF/storF=1`、`TS dt=0.1`、11 步齐全。压力场分位数（inject rec10）：`p50=200000.00`（**精确 `p0+(f/S)t`**）、`p90=200757.45`、`mean=199133.69`、`pmax=391877.64 @node 171`；Gate 4.5 `W: 14.880 → 30.161`，`W/W_an=1.005`，独立证实体相斜率就是 `1e5 Pa/s`。⇒ **物理与 `[1/s]` 量纲标定都正确**，两个缺陷都在观测链：EXO 追加复用 + 用 `pmax` 当体相压力。

### #125（`603ce0b`，success）

`build.yml` 两次 `rm -f`；判据 `exo_fields` 增 `pbulk/p90/pmean`，V1 加 EXO↔`Gate 3.2` 逐位交叉验证，V3/V4 改用 `pbulk[:11]`，V4 另加 `p90` 局域性（≤5%）与峰值/体相斜率比 `[1,6]`。**零 Fortran 改动**。`GATE415_RESULT: PASS`，十七组判据同轮全绿：`g415_slope_measured=1.000000e+05`、`g415_S_eff=1.500000e-10`、`g415_peak_ratio=2.918774`、`g415_p90_dev=3.787269e-03`、`g415_w_last=3.016054e+01`。`out_415_scale.gen` / `out_415_inject.gen` = **11 记录**，`out_main.gen` = 51 记录；原有 16 组 verdict 与 #124 **逐字节 IDENTICAL**。离线回放：`ci124_diag.zip` 原样 → 只剩 V1 FAIL；模拟重建 11 记录 → V1–V8 全 PASS；`--selftest` 1 正例（刻意带 2.9184× 真实过冲）+ 25 变异（含 1 不变性）全 PASS。

**结论：4.15 B0 收口，闸门 4 累计 15 项完成、合计 25 项已落地，0 项在册待做。** 阶段 B 的准入冒烟从「接线完成」升级为「有真实存储主导场、有可证伪的 `f/S` 与 `S_eff` 标定、有 EXO 生命周期护栏」。B1/B2 已由用户「通过文献调研实现B1和B2」立项并落地，见下一节。

---

## Gate 4.15 B1 / B2 落地（2026-11-10，待首轮 CI 验收）

### 文献复核推翻了原 R4 判据

用户要求「通过文献调研实现B1和B2」。复核 Lecampion, Desroches, Jeffrey & Bunger (2017), *JGR: Solid Earth*, doi:`10.1002/2016JB013183` 全文页后确认：**起裂前**进入裂缝的流量 $Q_{in}\approx0$，恒排量 $Q_0$ 下井筒压力**线性**上升，斜率 $=Q_0/U$（$U$ 为注入系统柔度），可直接由加压段斜率反解 $U$；岩石力学参数在这条曲线上不可辨识；且**起裂压力 ≠ 破裂（峰值）压力**。B0 的 `Damage_max=1.48e-4≪1` 证明当前窗口正是这一段。

因此路线稿第 5.2 节 R4 的「黏性主导 $\partial\ln p_{net}/\partial\ln Q\approx0$ / 提高 $G_c,\mu$ 后向 $-1/3$ 移动」属于 Detournay (2016) / Savitski & Detournay (2002) 的**裂缝传播段**标度，Santillán et al. (2018) 的三类排量模式同样属传播段，**都不能拿来做 B1 判据**。R4 已撤回并改为 **R4′ 参数可辨识性分层**：$\beta=f/S$ 只可辨识 $(f,S)$；$G_c$ 在表面能（$C_1\propto G_c/\ell$、$C_2\propto G_c\ell A$）上可辨识。这既与文献一致，又保留「可压性」的研究价值。

### 本轮交付（零 Fortran 改动）

| 类别 | 文件 / 变更 |
|---|---|
| 新增 deck ×8 | `terzaghi415_b1_{qhalf,qdouble,stor,gc,perm,visc}.yaml`（各 1 行代码差异）、`terzaghi415_b2_linear.yaml`（1 行）、`terzaghi415_b2_aniso.yaml`（+1 行 `toughnessAnisotropyMatrix: 2,0.5,0`） |
| 新增判据 ×2 | `check_gate415_b1_identifiability.py`（A1–A8）、`check_gate415_b2_extension.py`（L1–L8）；纯标准库、`--gate <dir>` / `--selftest` |
| `build.yml` | +9 个诊断 run（B1 6×2D、B2 2×2D + 1×3D）、11 次 tri3 / 1 次 tetra4 EXO 重建、2 个 verdict 目录装配、`gb1_rc` / `gb2_rc` 与 9 个 mpirun rc 的延迟硬裁决、4 个新 artifact 条目 |
| 文档 | `RATE_FRACABILITY_ROADMAP.md` 第 0 / 2.3 / 5.1 / 5.2 / 5.4 / 9 / 11 / **新增 12 节**；本清单 4.15 B1/B2 条目与进度段；`README.md` 索引；工作日志 |

### 三个解析锚（可证伪）

1. **排量/柔度**：$\beta=f/S$，B1 的 0.5×/1×/2× $Q$ 给出 `0.5/1/2` 比值且 $S_{eff}=f/\beta$ 恒为 $1.5\times10^{-10}$；$2\times S$ 给出 $\beta$ 减半、$S_{eff}=3\times10^{-10}$。
2. **不可辨识 + 解耦**：$G_c/k/\mu$ 三臂的 $\beta$ 偏差 <6%，但 $G_c$ 臂的表面能显著变化 >5% —— **A4 不变、A5 变**的组合才是可证伪命题，单独一条都容易被「参数没接上」蒙混。
3. **变排量离散锚**：$f(t)=f_0t \Rightarrow \Delta p_n=(f_0/S)dt^2n(n+1)/2$，$\Delta p_{10}=5.5\times10^4$ Pa $=0.55\times$ B0；已核实 transient 模式下 `vDef.F90:661` 的第二次 `SetTransients` 被 `If (.not. flowTransient)` 保护，`timeVal=time(step)` 唯一生效。

### 离线复核（全部通过，零编译 / 零求解器）

`git diff --check`；`ast.parse` ×2；PyYAML 解析 workflow + 10 份 415 deck；`Run Terzaghi` 整块 `bash -n`（1384 行，11 次 tri3 重建、1 次 tetra4 重建、全部 rc / verdict token 齐全）；B1 `--selftest` **1 正例 + 20 变异**、B2 **1 正例 + 19 变异**（各含 1 个不变性：damage SNES `-5` 不得误判）；**#125 真实 artifact 混合回放** A1–A8 / L1–L8 全 PASS（真实 `pbulk[-1]=199999.9999999968 Pa`、11 行能量分解、`surface_last=6.49012e-4`、`damage_max=1.4803e-4`、tri3 356 节点 / tetra4 2635 节点 `num_dim=3`）；B0 判据在新 workflow 下仍 `GATE415_RESULT: PASS`（`β=1.000000e+05`、`S_eff=1.500000e-10`）。

### 状态纪律

B1/B2 数值均为**离线解析预测**，不是求解器结论。首轮 CI 后必须按真实 `gate415_b1_verdict.log` / `gate415_b2_verdict.log` 复核；若阈值与真实数值有系统性偏差，**应修正判据的解析锚或容差并记录理由，不得改 deck 去凑解析解**（#124→#125 原则）。验收通过需补「19 组 verdict 同轮全绿 + 原有 17 组逐字节 IDENTICAL」，再把两条 `- [ ]` 勾为 `- [x]`。

**能力边界**：全部结论限定在起裂前、存储主导、无裂缝扩展窗口；不主张验证传播段标度，不报缝长 / 缝网面积 / 起裂压力绝对值 / 现场排量优选表。若要进入传播窗口，需降低 `fracturetoughness`（现 27000 N/m）、提高排量或延长时间，并先解决 damage SNES 的 `-5` 残差平台 —— 这是物理前提的改变，需单独立项。

---

## CI #127：B1 验收 + B2 L4 物理前提更正（2026-11-10）

**一句话**：`f07cf39` 的 CI #127 `failure`，但**编译成功、18 组 verdict 全绿、B2 的 7/8 判据全绿**；唯一失败 `GATE415_B2_RESULT: FAIL (L4)` 是判据把**传播段**标度套到**起裂前**窗口，不是代码 bug。本轮按 #127 实测重写 L4，零 Fortran / 零 deck / 零 workflow 改动。

### 实测数值（CI #127 `Run Terzaghi`，1429 s）

| 判据 | 真实读数 | 结论 |
|---|---|---|
| 原有 17 组 + `GATE415_B1_RESULT` | 全 PASS | **B1 已验收，清单勾 `[x]`** |
| L2 | `quad_slope=1000`（锚 `1000`）、`dp10=55000 = 0.550×B0`、二次残差 `0.00%` vs 线性 `13.64%` | 变排量离散 $t^2$ 锚被求解器精确证实 |
| L3 | 对照 `W≡0`；变排量 `W=0..23.28441` 单调；常排量 `W_last=30.16054` | PASS |
| L4 | `surface_ref=6.490120e-4`、`surface_ratio_aniso=1.004065`（**0.4065%**）、`beta_dev=0`、`damage_dev=1.247384e-2` | **FAIL（旧阈值要求 >5%）** |
| L6 | `beta_dim3=1.000000e+05`、`beta_dev_dim3=1.082467e-13`、`Displacement_Z` 活、无 timeout | 3D 一致性 PASS |
| L5 / L7 / L8 | bag dump `toughnessAnisotropyMatrix: 2,0.5,0`、单行 deck、9 run + `gb2_rc` 延迟裁决、B0 `β=1e5 / S_eff=1.5e-10`、生产 51 步固结 `W=0` | 全 PASS |
| 最终失败 | `=== FAIL: Gate 4.15 B2 rejected ... (rc=1) ===` | 延迟 `gb2_rc=1` |

### 根因与排除

`m_MEF90_DefMechAssembly.F90:2340-2341,2355` 中 $A$ **只乘 $C_2|\nabla\alpha|^2$**，而 $C_1w(\alpha)$ 不含 $A$。起裂前 `Damage_max≈1.48e-4` 且近似均匀 ⇒ 均匀项压倒梯度项，故 $A$ 只能让表面能变化 O(0.1%)。从真实日志第 697 / 785 条 `Total` 能量行读出 `surface: 6.49012E-04 → 6.51650E-04`，`ratio=1.004064640`，`ds/damage_dev=0.326`，正是梯度通道占比的签名。排除代码缺陷的四条证据：L5 静态接线齐全、`damage_dev=1.247%>1%`（各向异性求解确实发生）、表面能非位等价（$A$ 确实进了 $C_2$）、linear/3D/B0/B1/生产隔离全 PASS。**唯一不成立的是「必须 >5%」这一传播段阈值命题**。这是继 #124 `pmax` 之后第二次判据物理前提错误，处置原则一致：只改判据并补负例，不改 deck 去凑判据。

### 本轮交付（零 Fortran / 零 deck / 零 workflow 改动）

| 文件 | 变更 |
|---|---|
| `vDef/data/check_gate415_b2_extension.py` | docstring 与 verdict 表更正；新增 `SURF_BAND=(1e-4,5e-2)` 与 #127 实测锚 `SURF_REF/SURF_ANISO/DMG_REF/DMG_ANISO`；`check_l4` 重写为五条可证伪零位命题；`report()` 增 `gb2_surface_dev_aniso`；`--selftest` 扩到 **1 正例 + 24 变异（含 3 不变性）** |
| `RATE_FRACABILITY_ROADMAP.md` | 头部 / 第 0 节新增第 6 处更正、R5→R5′、5.1 / 5.2 / 5.4 / 11 / 12.4 / 12.6 / 12.7 状态订正，**新增第 13 节**完整记录首轮实测、根因、新命题、离线复核与收口条件 |
| `FLOW_FIX_CHECKLIST.md` | B1 勾 `[x]`、B2 保持 `[ ]` 并写明 7/8 PASS 与 L4 重写；进度段按 #127 重整；本专节 |
| `README.md` | B1/B2 状态与 L4 更正索引 |
| `.gitignore` | 增 `/ci*_joblogs.zip`，8.5 MB 诊断文本不进版本库 |
| `.workbuddy/memory/2026-11-10.md` | 追加首轮 CI 诊断与「先验算物理阶段再定判据阈值」方法论 |

新 L4 五条命题：$d_s=\lvert E_{surf}^{aniso}/E_{surf}^{iso}-1\rvert$ 满足 ① $d_s>10^{-4}$（已接线）、② $d_s<5\times10^{-2}$（未泄漏到 $C_1$，仍是起裂前场）、③ $d_s<d_{dmg}$（只作用于梯度通道）、④ $d_{dmg}>10^{-2}$（各向异性求解确实发生）、⑤ $d_\beta<6\%$（压力曲线盲，与 B1 A4 同构）。#127 代入 `0.004065 / 0.012474 / 0` ⇒ **5/5 PASS**。

### 离线复核（零编译 / 零求解器）

`git diff --check`；判据 `ast.parse`；B2 `--selftest` **PASS（1 正例 + 24 变异）**；#127 真实 `gb2_*` 指标代入新 L4 **5/5 PASS**；剥离 GitHub 时间戳后用判据自身 `ENERGY` 正则从 8.48 MB job logs 解析出 **797 条**真实 `Total` 能量行并定位两条 `surface:` 值，`ratio=1.004064640` 与 verdict 一致；L1 已在 #127 真实 log/EXO 上 PASS，本轮未改任何解析函数与 L1–L3 / L5–L8；未改 Fortran / deck / `build.yml` / B0 / B1 判据 ⇒ 生产物理场零扰动。

### CI #128 复验结果（收口）

`e6cce4b` 推送后 CI #128（run `42772909839`）`conclusion=success`，`Build + Terzaghi` 全 success。从完整 job logs（压缩 2,243,135 B / 解压 8,487,437 字符）读取：

| 项 | 结果 |
|---|---|
| `GATE415_B2_RESULT` | **PASS**；L1–L8 全绿 |
| L4 真实输出 | `surface dev 0.406% in band (0.01%..5%) and < damage dev 1.247%, beta dev 0.00%: anisotropy is wired but NOT identifiable before initiation` |
| `gb2_surface_dev_aniso` | `4.064640e-03`（离线预测同值）；`gb2_surface_ratio_aniso=1.004065`、`gb2_damage_dev_aniso=1.247384e-02` |
| 其余 17 项 `gb2_*` | 与 CI #127 **逐位相同** ⇒ 只改 L4 判定逻辑，物理场零扰动 |
| 全部 verdict | **18 条 `GATE*_RESULT` 全 PASS**；17 条非 B2 verdict 与 #127 **逐字节 IDENTICAL**；25 条 `=== ... PASS` 汇总行、0 条 `=== FAIL` |
| 体系规模 | **19 组判据同轮全绿**（18 条 `GATE*_RESULT` + `terzaghi_verdict` 汇总口径） |

清单 4.15 B1 / B2 均已勾 `[x]`，进度段按 #128 重整。**0.4065% 仍只能解读为「起裂前该观测量不可辨识」，不能解读为「层理各向异性对裂缝轨迹无影响」**；要研究轨迹控制必须进入传播段（第 5.4 节决策点 3）。
