# legacy/root-mirror — 根目录过期源码镜像归档（Gate 4.10）

归档日期：2026-11-09（来源状态：CI #107 / `627fe7d` 之后）
归档方式：`git mv`，**未删除任何文件**，git 历史完整可回溯。

## 为什么归档

本仓库是 Bourdin `mef90/vDef` 的 fork。新增的 `m_Flow/` 达西渗流模块以及
闸门 1–4 的所有修复都只落在子目录；历史上根目录仍保留了一整套同名
`.F90` / `.inc` 源码镜像和顶层 `Makefile_*` 镜像。精读确认：

- 顶层 `Makefile` 与 CI 只用
  `make -C objs/${PETSC_ARCH} -f ../../<子目录>/Makefile`；
  `MEF90/`、`m_DefMech/`、`m_HeatXfer/`、`m_Flow/`、`vDef/` 五个 Makefile
  的源路径分别是 `${MEF90_DIR}/MEF90/...`、
  `BASE_DIR=${MEF90_DIR}/m_Flow` 等，**没有任何规则引用根目录 `.F90`**。
- 66 个同名 Fortran 镜像中 **56 个与活源码逐字节相同**（纯冗余），
  **10 个已分叉**；另有 1 个 Makefile 镜像 `Makefile_m_Flow` 也已分叉。
  这些分叉副本既不会被编译、又长得像真源码，是最容易误改 / 误读的陷阱。
- 根目录 `test_*.F90` / `test20.F90` 是旧独立测试脚本；`Tests/Makefile`
  引用的是 `${MEF90_DIR}/Tests/*.F90`，与它们无关。
- 全仓库（Makefile / CI / shell / 文档）**没有任何一处引用**这些镜像
  Makefile 或旧测试脚本。

## 归档内容（93 个文件）

| 类别 | 数量 | 说明 |
| --- | --- | --- |
| 与活子目录同名的 `.F90` / `.inc` 镜像 | 66 | 56 个逐字节相同 + 10 个分叉 |
| `vDefDefault_remote.F90` | 1 | 根目录独有的 `vDefDefault` 远端镜像 |
| `test_*.F90` + `test20.F90` | 17 | 旧独立测试脚本，无人引用 |
| `Makefile_*` / `Makefile.include_*` 镜像 | 9 | 顶层 Makefile 镜像，无人引用 |

归档后根目录 tracked 的 `.F90` / `.inc` 为 **0**，只剩真正参与构建的
`Makefile` 与 `Makefile.include`。`MANIFEST.sha256` 记录全部 93 个文件的
sha256，由 W3 判据逐项校验。哈希基于 **git index 中的 blob 内容**（再做
CRLF→LF 规范化），即校验「提交进版本库的内容」而非工作区字节。这一点很
关键：CI 的 `Fix CRLF` 步骤会对全仓库 `.F90` 跑 `dos2unix`，它会移除
`vDefDefault_remote.F90` 的 UTF-8 BOM，从而改变该文件的工作区字节，但
版本库内容不变（CI #108 的 W3 误判由此而来）。因此不要用
`sha256sum -c` 直接校验工作区文件。

## 已分叉文件（最需要警惕，切勿当活源码使用）

| 归档文件 | 活源码 |
| --- | --- |
| `m_MEF90_DefMech.F90` | `m_DefMech/m_MEF90_DefMech.F90` |
| `m_MEF90_DefMechAssembly.F90` | `m_DefMech/m_MEF90_DefMechAssembly.F90` |
| `m_MEF90_DefMechCtx.F90` | `m_DefMech/m_MEF90_DefMechCtx.F90` |
| `m_MEF90_Flow.F90` | `m_Flow/m_MEF90_Flow.F90` |
| `m_MEF90_FlowAssembly.F90` | `m_Flow/m_MEF90_FlowAssembly.F90` |
| `m_MEF90_FlowCtx.F90` | `m_Flow/m_MEF90_FlowCtx.F90` |
| `m_MEF90_HeatXferAssembly.F90` | `m_HeatXfer/m_MEF90_HeatXferAssembly.F90` |
| `m_MEF90_HeatXferCtx.F90` | `m_HeatXfer/m_MEF90_HeatXferCtx.F90` |
| `m_vDefDefault.F90` | `vDef/m_vDefDefault.F90` |
| `vDef.F90` | `vDef/vDef.F90` |
| `Makefile_m_Flow` | `m_Flow/Makefile` |

## 活源码位置

| 归档前缀 | 活源码目录 |
| --- | --- |
| `m_MEF90_Flow*`, `m_MEF90_FlowCtx` | `m_Flow/` |
| `m_MEF90_DefMech*`, `mef90DefMech.inc` | `m_DefMech/` |
| `m_MEF90_HeatXfer*` | `m_HeatXfer/` |
| `m_MEF90*` 基础库, `mef90.inc` | `MEF90/` |
| `vDef*`, `m_vDefDefault.F90` | `vDef/` |
| `ThermoElasticity.F90`, `ThermoElastoPlasticity.F90` | 同名子目录 |

## 不参与构建的证明

`vDef/data/check_gate410_legacy.py` 的 W1–W8 硬判据已接入
`.github/workflows/build.yml`（统一「延迟 `exit 1`」）：

- W1 根目录不再有任何 tracked 的 `.F90` / `.inc` 镜像
- W2 归档目录文件数 / manifest 条目数完整
- W3 `MANIFEST.sha256` 每个哈希与磁盘文件一致
- W4 每个同名镜像都能在活子目录找到对应文件
- W5 五个活 Makefile 的所有 `.F90` 依赖都指向子目录
- W6 顶层 Makefile 与 CI 只使用子目录 Makefile，不引用镜像
- W7 任何构建入口都不引用 `legacy/`
- W8 本 README 存在且列明全部分叉文件

## 如何找回 / 使用

- 查看历史内容：`git show <归档前提交>:<原文件名>`，或直接在本目录打开。
- **不要把这些文件复制回根目录**：任何修复都应改子目录里的活源码，
  否则既不会进入编译，又会重新引入误改风险。
- 根目录的 `FlowCtx.F90` / `HeatXferCtx.F90` 不在本归档中：它们是
  `.gitignore` 明确忽略的本地诊断残留（untracked），不是仓库内容。
