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

- [ ] **3.2 用解析 Terzaghi 校验符号与量级**
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:101`（注释 `G(p)=div(A∇p)`，稳态算子符号）、`MEF90DiffusionOperatorSet` 符号；IFunction `F=K(p)+S·ṗ`、IJacobian `J=K+a·S·M`（BE 下 a=1/dt，形式自洽）
  - 问题：方程应为 `-div(k/μ∇p)=S ṗ + f`，但算子注释算的是 `+div(A∇p)`，与方程差一个符号；源项 `:105-107` 也用 `ADD_VALUES`（见 4.4）。
  - 改法：用解析解 `p(z,t)=(4σ₀/π)Σ(-1)ⁿ/(2n+1)·cos((2n+1)πz/2H)·exp(-(2n+1)²π²c_v t/4H²)`（c_v=k/(μS)，t=1s 时 T≈1.09）在单步与全程比对；若衰减方向反了，在 `MEF90FlowOperatorCore` 统一把 A 取负号，并同步源项符号。
  - 验收：数值衰减方向与量级与解析解一致。

---

## 闸门 4 — 物理增强与工程闭环

- [ ] **4.1 Biot 体力：叠加而非覆盖 + 补有效应力项**
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:948-953`（`VecCopy(baseForce,force)` 再叠加 `-α∇p`）；`vDef/vDef.F90:382-384`（每个 AltMin 迭代开始拍 `baseForce` 快照）
  - 问题：`baseForce` 是步初快照，`VecCopy` 重置体力 → DefMech 步内按 CST/Linear 缩放的体力被回滚成步初值。且力学侧缺 `σ' = σ + α·p·I` 的有效应力耦合（当前只做了体力近似）。
  - 改法：改为 `VecCopy(DefMechCtx%force, baseForce)` 仅作首次快照，AltMin 内用 `VecWAXPY(force, -biotCoeff, cellPressureGradient, baseForce)` 之类实现**叠加**而非覆盖；在 DefMech 应力装配处补 `α·p·I` 项。
  - 验收：开启体力缩放后结果不回滚；含孔压时应力场体现有效应力修正。

- [ ] **4.2 平流：把 `cellDarcyVelocity` 喂入 HeatXfer**
  - 位置：`m_HeatXfer/m_MEF90_HeatXferAssembly.F90:95-97 / 239-241 / 447-449 / 597-599`（已用常数 `advectionVector` 调 `MEF90DiffusionOperatorAdvectionSet`）；`m_Flow/...cellDarcyVelocity` 全库无消费者（grep 确认）
  - 问题：`enableFlowAdvection`/`fluidDensity` 是占位：HeatXfer 已有常数平流能力，但 Flow 算出的随空间变化的 `cellDarcyVelocity` 从未送入。
  - 改法：在 `MEF90HeatXferAssembly` 的平流分支，当 `enableFlowAdvection` 为真时用 `FlowCtx%cellDarcyVelocity` 替代 `cellSetOptions%advectionVector` 作为平流速度传入 `MEF90DiffusionOperatorAdvectionSet`。
  - 验收：开启 `enable_advection` 后热场受 v_darcy 平流影响；或若暂不支持，明确标注该选项。

- [ ] **4.3 边界压力时间缩放语义**
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:560-562`（`MEF90Scaling_Linear` 下 `VecScale(boundaryPressure, timeVal)`）
  - 问题：恒压边界（如大气压 0）物理上通常应保持常数，Linear 模式对其乘时间值语义可疑。
  - 改法：区分"恒压边界"与"渐变加载边界"两类选项；恒压走无缩放分支，渐变才乘时间/文件。
  - 验收：恒压 BC 在瞬态中保持常数；渐变 BC 按设定 ramp。

- [ ] **4.4 源项符号对齐**
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:105-107`（调 `MEF90DiffusionRHSSetCell`，底层 `ADD_VALUES`）
  - 问题：注释写"Subtract source RHS: G = DiffusionOperator - Source"，但底层是 `ADD_VALUES`，实际 `G = K + f`。
  - 改法：定死物理约定（`-div(k/μ∇p)=f` 对应 `K - f`），用独立 `rhsSec` 组装后 `SectionRealAXPY(resSec,-1,rhsSec)`；同步方程文档。
  - 验收：加均匀源时 p(z) 符号/量级与一维解析解一致。

- [ ] **4.5 能量诊断**
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:275-345`（`permFieldSec`/`sourceSec`/`p0Sec` 从未从 Vec 填入 → 能量恒 0）；主程序未调用 `MEF90FlowEnergy`
  - 改法：若做能量诊断，补三行反向 `SectionRealToVec` 散射填充；否则连同 `vDef.F90:37` 的变量一起删干净。
  - 验收：能量项非零（或确认已移除相关死代码）。

- [ ] **4.6 μ(T) 别名加固（脆弱但当前不崩）**
  - 位置：`m_Flow/m_MEF90_Flow.F90:270-278`（`UpdatePerm` 从 `cellDMScalSec` 同时取 `permPtr`/`muPtr`，二者同内存；`cellMuField` Vec 仅作 `Associated()` 守卫、值从未被读回）
  - 问题：`μ` 与 `perm` 靠"写入顺序"区分，依赖 `UpdateMuFromTemperature` 每次先行调用；一旦调用顺序变化即错。
  - 改法：为 μ 单独开 Section，或从持久 `FlowCtx%cellMuField` scatter 出独立 `muSec` 再读，消除顺序依赖。
  - 验收：开启 `temp_dep_viscosity` 后 `permField` 逐单元 = k(v)/μ(T)，与调用顺序无关。

- [ ] **4.7 时间相关材料参数入口**
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90:535-566`（`SetTransients` 仅缩放 source + boundaryPressure）
  - 问题：渗透率、储集系数、黏度等随时间演化（压密/损伤演化）无对应重载入口。
  - 改法：在 `SetTransients` 增加 permeability/storativity/viscosity 的 CST/Linear/File 缩放分支（或复用 `UpdatePerm`/`UpdateStorativity` 的时变重载）。
  - 验收：材料参数可随时间场演化。

- [ ] **4.8 守恒性（质量平衡）后处理**
  - 位置：新增诊断，建议挂在主循环时间步层
  - 问题：缺直接验证装配符号与边界条件正确性的手段。
  - 改法：每步计算 注入量 − 流出量 − 储存量变化，输出到日志/监控。
  - 验收：稳态下守恒残差趋零；瞬态下逐步闭合。

- [ ] **4.9 逐单元 IS 装配性能**
  - 位置：`m_Flow/m_MEF90_FlowAssembly.F90`（`MEF90FlowMakeSingleCellIS` 每单元一次 + 一次 Diffusion 调用）
  - 问题：装配复杂度 O(N) 次 PETSc 对象创建/销毁。
  - 改法：把 `A` 从"单元常数"放宽到"按 Gauss 点取值"（需改 Diffusion 接口或自写装配），或批量复用 IS。
  - 验收：大网格装配耗时显著下降。

- [ ] **4.10 清理根目录同名旧副本**
  - 位置：根目录 `m_MEF90_Flow.F90`(637行)、`m_MEF90_FlowAssembly.F90`(937行)、`FlowCtx.F90`、`m_MEF90_FlowCtx.F90`、`m_vDefDefault.F90` 等（与 `m_Flow/`、`vDef/` 分叉，首行 `Use` 列表不同）
  - 问题：同名副本极易误改/误读。
  - 改法：删除或移入归档目录；确认顶层 Makefile 只引用子目录。
  - 验收：根目录无参与编译的过期 `.F90` 副本。

- [ ] **4.11 CI 断言 + 回归 + 解析解比对**
  - 位置：`.github/workflows/build.yml:394-408`（Run Terzaghi，`continue-on-error: true` + `set +e` 仅 tee）；`vDef/Makefile` 的 `test` 只跑 test7/8
  - 问题：编译过但算错/崩也不会拦；`terzaghi1D.yaml` 未纳入回归。
  - 改法：去掉 `continue-on-error`，检查 `mpirun` 退出码；CI 加 `make test`；加 Python 脚本读 exo 的 pressure 与解析解 `p(z,t)` 比对（t=1s 时 T≈1.09），超阈值直接 fail。
  - 验收：CI 能拦住"编译过但算错/崩"的提交。

---

## 进度跟踪（更新 2026-11-06 19:55，CI #46 首个可信失败后）

**代码进度**
- 闸门 1（不崩+边界）：3 项改动完成 + C 项（单矩阵自拷贝 / 回调错配）已补，1.4 免改 — 3 / 4
- 闸门 2（可见）：4 / 4
- 闸门 3（瞬态+符号）：1 / 2（3.1 静态实现 + 二轮正确性修正；**3.1 代码至今从未被编译或运行**；3.2 必须运行）
- 闸门 4（物理+工程）：0 / 11（4.11 的「CI 可信化」前置已做，主体未动）
- 合计：8 项功能已落地，12 项待做

**提交 / CI 状态**
- 已提交 **`bd3c1f8`**（8 个文件，+523 / −49，含清单文档），2026-11-06 19:18 推送到 `origin/main`（基线 `2e49c2e`）。`.trae/`、`.workbuddy/` 属工具数据，未入库。
- CI run **#46 = failure**（11:44 UTC）。`Build MEF90 + m_Flow + vDef` 显示 success 是**假绿**：该步骤中途 `set +e` 后再未恢复，vDef 编译失败被静默吞掉；`Run Terzaghi` 因二进制不存在报 `mpirun ... unable to launch`（exit 132），随后的 Gate 3.1 自检断言只是连带误报 `transient TS path was NOT activated`。
- **真实根因**：`vDef/m_vDefDefault.F90:154:59` `Error: Invalid character in name` → `m_vDefDefault.o` 编译失败、vDef 未链接。对比历史 **#45** 的 `build-diagnostics`，有**逐字相同**的错误 → vDef 自 Flow 默认值块加入起从未成功构建，历史 success 全部只证明依赖链可用。
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

**状态（19:55 更新）**：CI #46 给出首个可信结论——**失败发生在构建层，3.1 的瞬态代码根本没进入运行**。修复已写入工作区，待提交推送后由 #47 复验。

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
