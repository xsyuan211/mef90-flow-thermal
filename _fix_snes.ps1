# 修复 m_MEF90_FlowAssembly.F90 的 SNES INTEGER(4) vs INTEGER(8) 类型错配
# 策略：提取 OperatorCore/BilinearFormCore，TS 回调直接调 Core（不带 SNES 参数）

$file = "d:\SoftBuild\Mef90+Flow+Dy\m_Flow\m_MEF90_FlowAssembly.F90"
$content = Get-Content $file -Raw

# 1. 把 MEF90FlowOperator 重命名为 MEF90FlowOperatorCore（去掉 SNES 参数，作为私有）
#    并把原来的 MEF90FlowOperator 变成薄包装
$oldOperatorHeader = @"
   Subroutine MEF90FlowOperator(snesPressure,x,residual,FlowCtx,ierr)
      PetscInt,Intent(IN)                                :: snesPressure   ! SNES opaque handle; unused inside
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
"@

$newOperatorHeader = @"
   Subroutine MEF90FlowOperatorCore(x,residual,FlowCtx,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
"@

$content = $content.Replace($oldOperatorHeader, $newOperatorHeader)

# 在 MEF90FlowOperatorCore 的 End Subroutine 之后，插入一个包装版本
$oldOperatorEnd = @"
      Call SectionRealDestroy(pSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowOperator
"@

$newOperatorEnd = @"
      Call SectionRealDestroy(pSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowOperatorCore


   !!!! SNES callback wrapper — matches PETSc SNESSetFunction signature
   Subroutine MEF90FlowOperator(snes,x,residual,FlowCtx,ierr)
      PetscInt,Intent(IN)                                :: snes
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Call MEF90FlowOperatorCore(x,residual,FlowCtx,ierr)
   End Subroutine MEF90FlowOperator
"@

$content = $content.Replace($oldOperatorEnd, $newOperatorEnd)

# 2. 同理处理 MEF90FlowBilinearForm
$oldBilinearHeader = @"
   Subroutine MEF90FlowBilinearForm(snesPressure,x,Jac,PreJac,residualMode,FlowCtx,ierr)
      PetscInt,Intent(IN)                                :: snesPressure   ! SNES opaque handle; unused inside
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      PetscReal,Intent(IN)                               :: residualMode
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
"@

$newBilinearHeader = @"
   Subroutine MEF90FlowBilinearFormCore(x,Jac,PreJac,residualMode,FlowCtx,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      PetscReal,Intent(IN)                               :: residualMode
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
"@

$content = $content.Replace($oldBilinearHeader, $newBilinearHeader)

$oldBilinearEnd = @"
      Call SectionRealDestroy(permFieldSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowBilinearForm
"@

$newBilinearEnd = @"
      Call SectionRealDestroy(permFieldSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowBilinearFormCore


   !!!! SNES callback wrapper — matches PETSc SNESSetJacobian signature
   Subroutine MEF90FlowBilinearForm(snes,x,Jac,PreJac,residualMode,FlowCtx,ierr)
      PetscInt,Intent(IN)                                :: snes
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      PetscReal,Intent(IN)                               :: residualMode
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Call MEF90FlowBilinearFormCore(x,Jac,PreJac,residualMode,FlowCtx,ierr)
   End Subroutine MEF90FlowBilinearForm
"@

$content = $content.Replace($oldBilinearEnd, $newBilinearEnd)

# 3. 修改 TS 回调 — 调用 Core 版本，不再传 PETSC_NULL_OBJECT
$content = $content.Replace(
    '      Call MEF90FlowOperator(PETSC_NULL_OBJECT,x,F,FlowCtx,ierr)',
    '      Call MEF90FlowOperatorCore(x,F,FlowCtx,ierr)'
)
$content = $content.Replace(
    '      Call MEF90FlowBilinearForm(PETSC_NULL_OBJECT,x,Jac,PreJac,1.0_Kr,FlowCtx,ierr)',
    '      Call MEF90FlowBilinearFormCore(x,Jac,PreJac,1.0_Kr,FlowCtx,ierr)'
)

Set-Content $file $content -NoNewline
Write-Host "DONE: Core + wrapper pattern applied"
