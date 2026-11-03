#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_FlowAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_FlowCtx
   Use m_MEF90_DiffusionInterface   !!! Reuse Bourdin's well-tested diffusion assembly
   Use m_MEF90_DefMechCtx_Type      !!! Biot coupling: write f_p into DefMechCtx%force
   Implicit none

   Private
   Public :: MEF90FlowOperator,           &   ! SNES residual callback
             MEF90FlowBilinearForm,       &   ! SNES Jacobian callback
             MEF90FlowIFunction,          &   ! TS residual callback
             MEF90FlowIJacobian,          &   ! TS Jacobian callback
             MEF90FlowEnergy,             &   ! Darcy dissipation + injection work
             MEF90FlowViewEXO,            &   ! Save pressure/perm/source to Exodus
             MEF90FlowUpdateboundaryPressure, & ! Set Dirichlet p from PetscBag
             MEF90FlowSetTransients,      &   ! Apply time-scaling to source + BC
             MEF90FlowSetSourceCst,       &   ! Fill source from cell-set PetscBag
             MEF90FlowComputePressureGradient, & ! Cell-centered ∇p volume average
             MEF90FlowApplyBiotForce          ! Biot f_p = -α∇p → DefMechCtx%force

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  INNER LOOP HELPERS — per-cell assembly of Diffusion terms
!!!  MEF90DiffusionOperatorSet/BilinearFormSet take CONSTANT A (MEF90_MATS)
!!!  For per-cell-varying A = permField(cell)/mu, we build a single-cell IS each loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Build an IS containing exactly one cell (cellID)
   Subroutine MEF90FlowMakeSingleCellIS(mesh,cellID,singleCellIS,ierr)
      Type(DM),Intent(IN)                                 :: mesh
      PetscInt,Intent(IN)                                  :: cellID
      Type(IS),Intent(OUT)                                :: singleCellIS
      PetscErrorCode,Intent(OUT)                           :: ierr
      PetscInt,Dimension(1)                               :: idArr
      idArr(1) = cellID
      Call ISCreateGeneral(PETSC_COMM_SELF,1,idArr,PETSC_COPY_VALUES,singleCellIS,ierr)
   End Subroutine MEF90FlowMakeSingleCellIS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SNES RESIDUAL:  G(p) = DiffusionOperator(p) - SourceRHS  +  Dirichlet BC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowOperator"

   Subroutine MEF90FlowOperator(snesPressure,x,residual,FlowCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesPressure
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(SectionReal)                                  :: pSec, resSec
      Type(SectionReal)                                  :: permFieldSec, sourceSec, rhsSec
      Type(IS)                                           :: CellSetGlobalIS, setIS, singleIS
      PetscInt,dimension(:),Pointer                      :: setID, cellID
      PetscInt                                           :: set, cell, QuadratureOrder, nVal
      Type(MEF90FlowCellSetOptions_Type),pointer        :: cellSetOptions
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
      Type(MEF90Element_Type)                            :: elemType
      Type(MEF90_MATS)                                   :: A_cell
      PetscReal,Dimension(1)                              :: permCell, sourceCell

      !!! Create sections
      Call SectionRealDuplicate(FlowCtx%DMScalSec,pSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%DMScalSec,resSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,permFieldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,sourceSec,ierr);CHKERRQ(ierr)

      !!! Scatter global Vecs → local Sections
      Call SectionRealToVec(pSec,FlowCtx%DMScalScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(permFieldSec,FlowCtx%cellDMScalScatter,SCATTER_REVERSE,FlowCtx%permField,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(sourceSec,FlowCtx%cellDMScalScatter,SCATTER_REVERSE,FlowCtx%source,ierr);CHKERRQ(ierr)
      Call SectionRealSet(resSec,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (.not. cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
            QuadratureOrder = elemType%order * 2
            Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
                 cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)

            If (elemType%Codim == 0) Then
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr)

                  !!! A = permField(cell) (k/μ scalar)
                  Call SectionRealRestrict(permFieldSec,cellID(cell),permCell,ierr);CHKERRQ(ierr)
                  A_cell = permCell(1)

                  !!! G(p) = div(A∇p) with λ=0 (steady storage handled separately)
                  Call MEF90DiffusionOperatorSet(resSec,FlowCtx%DM,pSec,singleIS, &
                       A_cell,0.0_Kr,elem,elemType,ierr)

                  !!! Subtract source RHS: G = DiffusionOperator - Source
                  Call MEF90DiffusionRHSSetCell(resSec,FlowCtx%DM,sourceSec,singleIS, &
                       elem,elemType,ierr)

                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If

            Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr)
         End If
      End Do

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      Call SectionRealComplete(resSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(resSec,FlowCtx%DMScalScatter,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)

      !!!
      !!! DIRICHLET BC: G(p_bc) = p - p_bc = 0
      !!!
      Call MEF90FlowApplyDirichletResidual(residual,FlowCtx,ierr)

      Call VecAssemblyBegin(residual,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(residual,ierr);CHKERRQ(ierr)

      Call SectionRealDestroy(permFieldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(sourceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(resSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(pSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowOperator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SNES JACOBIAN:  K = dG/dp = Diffusion Bilinear Form (per-cell varying A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowBilinearForm"

   Subroutine MEF90FlowBilinearForm(snesPressure,x,Jac,PreJac,residualMode,FlowCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesPressure
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      PetscReal,Intent(IN)                               :: residualMode
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(IS)                                           :: CellSetGlobalIS, setIS, singleIS
      PetscInt,dimension(:),Pointer                      :: setID, cellID
      PetscInt                                           :: set, cell, QuadratureOrder
      Type(MEF90FlowCellSetOptions_Type),pointer        :: cellSetOptions
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      Type(MEF90Element_Type)                            :: elemType
      Type(MEF90_MATS)                                   :: A_cell
      Type(SectionReal)                                  :: permFieldSec
      PetscReal,Dimension(1)                              :: permCell

      Call MatZeroEntries(Jac,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,permFieldSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(permFieldSec,FlowCtx%cellDMScalScatter,SCATTER_REVERSE,FlowCtx%permField,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (.not. cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
            QuadratureOrder = elemType%order * 2
            Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
                 cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
            If (elemType%Codim == 0) Then
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr)
                  Call SectionRealRestrict(permFieldSec,cellID(cell),permCell,ierr);CHKERRQ(ierr)
                  A_cell = permCell(1)
                  Call MEF90DiffusionBilinearFormSet(Jac,FlowCtx%DM,singleIS, &
                       A_cell,0.0_Kr,elem,elemType,ierr)
                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If
            Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr)
         End If
      End Do

      Call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatCopy(Jac,PreJac,DIFFERENT_NONZERO_PATTERN,ierr);CHKERRQ(ierr)

      !!! DIRICHLET BC: K_ii = 1 for all BC dofs
      Call MEF90FlowApplyDirichletJacobian(Jac,PreJac,FlowCtx,ierr)

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)
      Call SectionRealDestroy(permFieldSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowBilinearForm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  TS RESIDUAL + JACOBIAN — reuse Operator/BilinearForm, add storage term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowIFunction"

   Subroutine MEF90FlowIFunction(ts,t,x,xDot,F,FlowCtx,ierr)
      Type(TS),Intent(IN)                                :: ts
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x, xDot
      Type(Vec),Intent(INOUT)                            :: F
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Call MEF90FlowOperator(PETSC_NULL_OBJECT,x,F,FlowCtx,ierr)
      Call MEF90FlowAddStorageToResidual(F,xDot,FlowCtx,ierr)
   End Subroutine MEF90FlowIFunction


#undef __FUNCT__
#define __FUNCT__ "MEF90FlowIJacobian"

   Subroutine MEF90FlowIJacobian(ts,t,x,xDot,a,Jac,PreJac,FlowCtx,ierr)
      Type(TS),Intent(IN)                                :: ts
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x, xDot
      PetscReal,Intent(IN)                               :: a
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      !!! First, assemble steady Darcy Jacobian K = dG/dp
      Call MEF90FlowBilinearForm(PETSC_NULL_OBJECT,x,Jac,PreJac,1.0_Kr,FlowCtx,ierr)
      !!! Add a·S·M (mass matrix):  K += a·S·BFB
      Call MEF90FlowAddStorageMassMatrix(Jac,PreJac,a,FlowCtx,ierr)
   End Subroutine MEF90FlowIJacobian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  ENERGY:  E = ½∫ (k/μ)|∇p|² dΩ   W = ∫ f·p dΩ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowEnergy"

   Subroutine MEF90FlowEnergy(pressure,timeVal,FlowCtx,Energy,InjectionWork,ierr)
      Type(Vec),Intent(IN)                               :: pressure
      PetscReal,Intent(IN)                               :: timeVal
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscReal,Intent(OUT)                              :: Energy, InjectionWork
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(SectionReal)                                  :: pSec, p0Sec, permFieldSec, sourceSec
      Type(IS)                                           :: CellSetGlobalIS, setIS, singleIS
      PetscInt,dimension(:),Pointer                      :: setID, cellID
      PetscInt                                           :: set, cell, QuadratureOrder
      Type(MEF90FlowCellSetOptions_Type),pointer        :: cellSetOptions
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      Type(MEF90Element_Type)                            :: elemType
      Type(MEF90_MATS)                                   :: A_cell
      PetscReal,Dimension(1)                              :: permCell
      PetscReal                                          :: energyCell, workCell

      Energy = 0.0_Kr
      InjectionWork = 0.0_Kr

      Call SectionRealDuplicate(FlowCtx%DMScalSec,pSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,p0Sec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,permFieldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,sourceSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(pSec,FlowCtx%DMScalScatter,SCATTER_REVERSE,pressure,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (.not. cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
            QuadratureOrder = elemType%order * 2
            Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
                 cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)

            If (elemType%Codim == 0) Then
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr)
                  Call SectionRealRestrict(permFieldSec,cellID(cell),permCell,ierr);CHKERRQ(ierr)
                  A_cell = permCell(1)

                  Call MEF90DiffusionEnergySet(energyCell,pSec,p0Sec,FlowCtx%DM, &
                       A_cell,0.0_Kr,singleIS,elem,elemType,ierr)
                  Energy = Energy + energyCell

                  Call MEF90DiffusionWorkSetCell(workCell,pSec,FlowCtx%DM, &
                       sourceSec,singleIS,elem,elemType,ierr)
                  InjectionWork = InjectionWork + workCell

                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If
            Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      Call SectionRealDestroy(pSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(p0Sec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(permFieldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(sourceSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowEnergy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  EXODUS OUTPUT — pressure (vertex), source (cell), permField (cell)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowViewEXO"

   Subroutine MEF90FlowViewEXO(FlowCtx,step,ierr)
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(Vec)                                          :: localVec

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)

      !!! Cell-centered fields
      If (gOpt%sourceOffset > 0 .AND. Associated(FlowCtx%source)) Then
         Call VecViewExodusCell(FlowCtx%cellDMScal,FlowCtx%source,FlowCtx%MEF90Ctx%IOcomm, &
                                FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%sourceOffset,ierr);CHKERRQ(ierr)
      End If

      !!! Vertex-centered: pressure + boundaryPressure
      Call DMGetLocalVector(FlowCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      If (gOpt%pressureOffset > 0 .AND. Associated(FlowCtx%pressure)) Then
         Call DMGlobalToLocalBegin(FlowCtx%DMScal,FlowCtx%pressure,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd  (FlowCtx%DMScal,FlowCtx%pressure,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(FlowCtx%DMScal,localVec,FlowCtx%MEF90Ctx%IOcomm, &
                                  FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%pressureOffset,ierr);CHKERRQ(ierr)
      End If
      If (gOpt%boundaryPressureOffset > 0 .AND. Associated(FlowCtx%boundaryPressure)) Then
         Call DMGlobalToLocalBegin(FlowCtx%DMScal,FlowCtx%boundaryPressure,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd  (FlowCtx%DMScal,FlowCtx%boundaryPressure,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(FlowCtx%DMScal,localVec,FlowCtx%MEF90Ctx%IOcomm, &
                                  FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%boundaryPressureOffset,ierr);CHKERRQ(ierr)
      End If
      Call DMRestoreLocalVector(FlowCtx%DMScal,localVec,ierr);CHKERRQ(ierr)

      If (FlowCtx%MEF90Ctx%rank == 0) Then
         Call EXUPDA(FlowCtx%MEF90Ctx%fileExoUnit,ierr)
      End If
   End Subroutine MEF90FlowViewEXO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  DIRICHLET BC — for SNES residual and Jacobian
!!!  Residual:  G(dof) = x(dof) - bc(dof)
!!!  Jacobian:  K(dof,:) = 0,  K(dof,dof) = 1.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Subroutine MEF90FlowApplyDirichletResidual(residual,FlowCtx,ierr)
      Type(Vec),Intent(INOUT)                             :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                  :: FlowCtx
      PetscErrorCode,Intent(OUT)                           :: ierr

      Type(IS)                                            :: CellSetGlobalIS, bcIS
      Type(IS)                                            :: VertexSetGlobalIS, setIS, cellSetIS
      PetscInt,dimension(:),Pointer                        :: setID, cellIDs, dofList
      PetscInt                                            :: set, nVal
      PetscReal,Dimension(:)                              :: resPtr, bcPtr
      Type(MEF90FlowCellSetOptions_Type),pointer          :: cellSetOptions
      Type(MEF90FlowVertexSetOptions_Type),pointer        :: vertexSetOptions

      !!! Cell sets with Has_BC → Cell→Vertex IS
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),cellSetIS,ierr);CHKERRQ(ierr)
            Call MEF90ISCreateCelltoVertex(FlowCtx%DMScal,PETSC_COMM_WORLD,cellSetIS,bcIS,ierr)
            Call ISGetSize(bcIS,nVal,ierr);CHKERRQ(ierr)
            Allocate(dofList(nVal),stat=ierr)
            Call ISGetIndicesF90(bcIS,dofList,ierr);CHKERRQ(ierr)
            Call VecGetValues(residual,nVal,dofList,resPtr,ierr);CHKERRQ(ierr)
            Call VecGetValues(FlowCtx%boundaryPressure,nVal,dofList,bcPtr,ierr);CHKERRQ(ierr)
            resPtr = resPtr - bcPtr    !!! G = x - p_bc on Dirichlet dofs
            Call VecSetValues(residual,nVal,dofList,resPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(bcIS,dofList,ierr);CHKERRQ(ierr)
            DeAllocate(dofList)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(cellSetIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      !!! Vertex sets with Has_BC → override
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxVertexSetOptions(FlowCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DMScal,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            Call ISGetSize(setIS,nVal,ierr);CHKERRQ(ierr)
            Allocate(dofList(nVal),resPtr(nVal),bcPtr(nVal),stat=ierr)
            Call ISGetIndicesF90(setIS,dofList,ierr);CHKERRQ(ierr)
            Call VecGetValues(residual,nVal,dofList,resPtr,ierr);CHKERRQ(ierr)
            Call VecGetValues(FlowCtx%boundaryPressure,nVal,dofList,bcPtr,ierr);CHKERRQ(ierr)
            resPtr = resPtr - bcPtr
            Call VecSetValues(residual,nVal,dofList,resPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(setIS,dofList,ierr);CHKERRQ(ierr)
            DeAllocate(dofList,resPtr,bcPtr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr)
   End Subroutine MEF90FlowApplyDirichletResidual


   Subroutine MEF90FlowApplyDirichletJacobian(Jac,PreJac,FlowCtx,ierr)
      Type(Mat),Intent(INOUT)                             :: Jac, PreJac
      Type(MEF90FlowCtx_Type),Intent(IN)                  :: FlowCtx
      PetscErrorCode,Intent(OUT)                           :: ierr

      Type(IS)                                            :: CellSetGlobalIS, bcIS
      Type(IS)                                            :: VertexSetGlobalIS, setIS, cellIS
      PetscInt,dimension(:),Pointer                        :: setID
      PetscInt                                            :: set, nVal, i
      PetscReal                                           :: diagVal
      PetscInt,Dimension(:),Pointer                       :: dofList
      Type(MEF90FlowCellSetOptions_Type),pointer          :: cellSetOptions
      Type(MEF90FlowVertexSetOptions_Type),pointer        :: vertexSetOptions

      !!! DIRICHLET JACOBIAN: K(i,i) = 1 for all BC dofs (via MatSetValues below)

      !!! Cell sets with Has_BC
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_BC) Then
            !!! cellIS = cell-set StratumIS, bcIS = Cell→Vertex dof IS, dofList = vertex dof IDs
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),cellIS,ierr);CHKERRQ(ierr)
            Call MEF90ISCreateCelltoVertex(FlowCtx%DMScal,PETSC_COMM_WORLD,cellIS,bcIS,ierr)
            Call ISGetLocalSize(bcIS,nVal,ierr);CHKERRQ(ierr)
            Allocate(dofList(nVal),stat=ierr)
            Call ISGetIndicesF90(bcIS,dofList,ierr);CHKERRQ(ierr)
            !!! MatZeroRows zeroes rows dofList(1:nVal), sets diagonal = 1.0 (standard Dirichlet Jacobian pattern)
            Call MatZeroRows(Jac,nVal,dofList,1.0_Kr,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call MatZeroRows(PreJac,nVal,dofList,1.0_Kr,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(bcIS,dofList,ierr);CHKERRQ(ierr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
            DeAllocate(dofList)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      !!! Vertex sets with Has_BC
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxVertexSetOptions(FlowCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DMScal,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            Call ISGetLocalSize(setIS,nVal,ierr);CHKERRQ(ierr)
            Allocate(dofList(nVal),stat=ierr)
            Call ISGetIndicesF90(setIS,dofList,ierr);CHKERRQ(ierr)
            Call MatZeroRows(Jac,nVal,dofList,1.0_Kr,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call MatZeroRows(PreJac,nVal,dofList,1.0_Kr,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(setIS,dofList,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            DeAllocate(dofList)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr)

      Call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyBegin(PreJac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd(PreJac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowApplyDirichletJacobian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SET TRANSIENTS — apply time-scaling law to source + pressure BC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowSetTransients"

   Subroutine MEF90FlowSetTransients(FlowCtx,step,timeVal,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      PetscInt,Intent(IN)                                :: step
      PetscReal,Intent(IN)                               :: timeVal
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)

      Select case (gOpt%sourceScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(FlowCtx%cellDMScal,FlowCtx%source,FlowCtx%MEF90Ctx%IOcomm, &
                                    FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%sourceOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90FlowSetSourceCst(FlowCtx%source,FlowCtx,ierr)
            Call VecScale(FlowCtx%source,timeVal,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90FlowSetSourceCst(FlowCtx%source,FlowCtx,ierr)
      End Select

      Select case (gOpt%boundaryPressureScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusVertex(FlowCtx%DMScal,FlowCtx%boundaryPressure,FlowCtx%MEF90Ctx%IOcomm, &
                                     FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%boundaryPressureOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90FlowUpdateboundaryPressure(FlowCtx%boundaryPressure,FlowCtx,ierr)
            Call VecScale(FlowCtx%boundaryPressure,timeVal,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90FlowUpdateboundaryPressure(FlowCtx%boundaryPressure,FlowCtx,ierr)
      End Select
   End Subroutine MEF90FlowSetTransients


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  UPDATE BOUNDARY PRESSURE — fill boundaryPressure Vec from PetscBag options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowUpdateboundaryPressure"

   Subroutine MEF90FlowUpdateboundaryPressure(boundaryPressure,FlowCtx,ierr)
      Type(Vec),Intent(INOUT)                            :: boundaryPressure
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90FlowVertexSetOptions_Type),pointer       :: vertexSetOptions
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS, cellSetIS, bcIS, VertexSetGlobalIS, setIS
      PetscInt,dimension(:),Pointer                      :: setID, dofList
      PetscInt                                           :: set, nVal
      PetscReal,Dimension(:)                             :: val

      Call VecSet(boundaryPressure,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),cellSetIS,ierr);CHKERRQ(ierr)
            Call MEF90ISCreateCelltoVertex(FlowCtx%DMScal,PETSC_COMM_WORLD,cellSetIS,bcIS,ierr)
            Call ISGetSize(bcIS,nVal,ierr);CHKERRQ(ierr)
            Allocate(dofList(nVal),val(nVal),stat=ierr)
            Call ISGetIndicesF90(bcIS,dofList,ierr);CHKERRQ(ierr)
            val = cellSetOptions%boundaryPressure
            Call VecSetValues(boundaryPressure,nVal,dofList,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(bcIS,dofList,ierr);CHKERRQ(ierr)
            DeAllocate(dofList,val)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(cellSetIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxVertexSetOptions(FlowCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DMScal,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            Call ISGetSize(setIS,nVal,ierr);CHKERRQ(ierr)
            Allocate(dofList(nVal),val(nVal),stat=ierr)
            Call ISGetIndicesF90(setIS,dofList,ierr);CHKERRQ(ierr)
            val = vertexSetOptions%boundaryPressure
            Call VecSetValues(boundaryPressure,nVal,dofList,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(setIS,dofList,ierr);CHKERRQ(ierr)
            DeAllocate(dofList,val)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr)

      Call VecAssemblyBegin(boundaryPressure,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(boundaryPressure,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowUpdateboundaryPressure


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SET SOURCE CST — fill FlowCtx%source from cellSetOptions%source
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowSetSourceCst"

   Subroutine MEF90FlowSetSourceCst(sourceVec,FlowCtx,ierr)
      Type(Vec),Intent(INOUT)                            :: sourceVec
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS, setIS, setISdof
      PetscInt,dimension(:),Pointer                      :: setID, setIdx
      PetscInt                                           :: set, nVal
      PetscReal,Dimension(:),Pointer                     :: val

      Call VecSet(sourceVec,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call DMMeshISCreateISglobaldof(FlowCtx%cellDMScal,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)   !!! separate ptr
         Call ISGetLocalSize(setISdof,nVal,ierr);CHKERRQ(ierr)
         Allocate(val(nVal),stat=ierr)
         val = cellSetOptions%source
         Call VecSetValues(sourceVec,nVal,setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(val)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      Call VecAssemblyBegin(sourceVec,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(sourceVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowSetSourceCst


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  ADD STORAGE MASS TERM — for transient Darcy
!!!  Residual += S·ṗ    (using DiffusionOperatorAddTransientTermSet with alpha = cell S)
!!!  Jacobian += a·S·BFB (using DiffusionBilinearFormSet with A=0, lambda = a·S)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowAddStorageToResidual"

   Subroutine MEF90FlowAddStorageToResidual(F,pDot,FlowCtx,ierr)
      Type(Vec),Intent(INOUT)                             :: F
      Type(Vec),Intent(IN)                                :: pDot
      Type(MEF90FlowCtx_Type),Intent(IN)                  :: FlowCtx
      PetscErrorCode,Intent(OUT)                           :: ierr

      Type(SectionReal)                                   :: FSec, pDotSec, storativitySec
      Type(IS)                                           :: CellSetGlobalIS, setIS, singleIS
      PetscInt,dimension(:),Pointer                       :: setID, cellID
      PetscInt                                            :: set, cell, QuadratureOrder
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellSetOptions
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer       :: elem
      Type(MEF90Element_Type)                             :: elemType
      PetscReal,Dimension(1)                              :: SCell
      PetscReal                                           :: S_avg
      Type(MEF90_MATS)                                    :: A_zero

      Call SectionRealDuplicate(FlowCtx%DMScalSec,FSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%DMScalSec,pDotSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,storativitySec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(FSec,FlowCtx%DMScalScatter,SCATTER_REVERSE,F,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(pDotSec,FlowCtx%DMScalScatter,SCATTER_REVERSE,pDot,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(storativitySec,FlowCtx%cellDMScalScatter,SCATTER_REVERSE,FlowCtx%storativity,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (.not. cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
            QuadratureOrder = elemType%order * 2
            Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
                 cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
            If (elemType%Codim == 0) Then
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr)
                  !!! alpha = S(cell), pre-multiply and pass to DiffusionOperatorAddTransientTermSet
                  Call SectionRealRestrict(storativitySec,cellID(cell),SCell,ierr);CHKERRQ(ierr)
                  S_avg = SCell(1)
                  Call MEF90DiffusionOperatorAddTransientTermSet(FSec,FlowCtx%DM,pDotSec, &
                       singleIS,S_avg,elem,elemType,ierr)
                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If
            Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr)
         End If
      End Do

      Call SectionRealComplete(FSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(FSec,FlowCtx%DMScalScatter,SCATTER_FORWARD,F,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)
      Call SectionRealDestroy(FSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(pDotSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(storativitySec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowAddStorageToResidual


#undef __FUNCT__
#define __FUNCT__ "MEF90FlowAddStorageMassMatrix"

   Subroutine MEF90FlowAddStorageMassMatrix(Jac,PreJac,a,FlowCtx,ierr)
      Type(Mat),Intent(INOUT)                             :: Jac, PreJac
      PetscReal,Intent(IN)                                :: a
      Type(MEF90FlowCtx_Type),Intent(IN)                  :: FlowCtx
      PetscErrorCode,Intent(OUT)                           :: ierr

      Type(IS)                                           :: CellSetGlobalIS, setIS, singleIS
      PetscInt,dimension(:),Pointer                       :: setID, cellID
      PetscInt                                           :: set, cell, QuadratureOrder
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellSetOptions
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer       :: elem
      Type(MEF90Element_Type)                             :: elemType
      PetscReal,Dimension(1)                              :: SCell
      PetscReal                                           :: lambda_val
      Type(MEF90_MATS)                                    :: A_zero

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (.not. cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
            QuadratureOrder = elemType%order * 2
            Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
                 cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
            If (elemType%Codim == 0) Then
               A_zero = 0.0_Kr
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr)
                  Call SectionRealRestrict(FlowCtx%storativity,cellID(cell),SCell,ierr);CHKERRQ(ierr)
                  lambda_val = a * SCell(1)
                  !!! K += 0*∇ + λ·BFB   (lambda term only → mass matrix)
                  Call MEF90DiffusionBilinearFormSet(Jac,FlowCtx%DM,singleIS, &
                       A_zero,lambda_val,elem,elemType,ierr)
                  Call MEF90DiffusionBilinearFormSet(PreJac,FlowCtx%DM,singleIS, &
                       A_zero,lambda_val,elem,elemType,ierr)
                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If
            Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr)
         End If
      End Do

      Call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyBegin(PreJac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd(PreJac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)
   End Subroutine MEF90FlowAddStorageMassMatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  CELL-CENTERED PRESSURE GRADIENT — volume-weighted average per cell
!!!  ∇p_cell = Σ_gauss (∇p_gauss * Gauss_C) / Σ_gauss Gauss_C
!!!  Biot body force  f_p = -α · ∇p_cell   (applied as volume force in DefMech)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowComputePressureGradient"

   Subroutine MEF90FlowComputePressureGradient(pressure,cellPressureGradient,FlowCtx,ierr)
      Type(Vec),Intent(IN)                               :: pressure
      Type(Vec),Intent(INOUT)                            :: cellPressureGradient
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(SectionReal)                                  :: pSec, gradpSec
      Type(IS)                                           :: CellSetGlobalIS, setIS
      PetscInt,dimension(:),Pointer                      :: setID, cellID
      PetscInt                                           :: set, cell, QuadratureOrder
      PetscInt                                           :: iDoF, iGauss, numDof, numGauss
      PetscInt                                           :: dim
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellSetOptions
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer     :: elem
      Type(MEF90Element_Type)                             :: elemType
      PetscReal,Dimension(:),Pointer                     :: pDofLoc
      PetscReal,Dimension(:),Pointer                     :: secPtr          !!! 独立 Section 指针
      PetscReal,Dimension(MEF90_DIM)                     :: gradpGauss, gradpCell
      PetscReal                                          :: volumeWeight, totalVolume

      dim = MEF90_DIM

      Call SectionRealDuplicate(FlowCtx%DMScalSec,pSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMVectSec,gradpSec,ierr);CHKERRQ(ierr)
      Call SectionRealSet(gradpSec,0.0_Kr,ierr);CHKERRQ(ierr)

      Call SectionRealToVec(pSec,FlowCtx%DMScalScatter,SCATTER_REVERSE,pressure,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         !!! 计算 ∇p 不需要跳过 BC set——BC set 的 cell 内压力也有梯度
         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
         QuadratureOrder = elemType%order * 2
         Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
              cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)

         If (elemType%Codim == 0) Then
            numDof  = elemType%numDof
            Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            Allocate(pDofLoc(numDof))
            Do cell = 1, size(cellID)
               !!! 读 cell 内所有 pressure Dof 值
               Call SectionRealRestrictClosure(pSec,FlowCtx%DM,cellID(cell),numDof,pDofLoc,ierr);CHKERRQ(ierr)
               numGauss = size(elem(cell)%BF,2)

               !!! volume-weighted average of ∇p over Gauss points
               gradpCell  = 0.0_Kr
               totalVolume = 0.0_Kr
                  Do iGauss = 1, numGauss
                  gradpGauss = 0.0_Kr
                  Do iDoF = 1, numDof
                     !!! elem(cell)%Grad_BF(iDoF,iGauss) is Vect2D/Vect3D
                     gradpGauss(1) = gradpGauss(1) + pDofLoc(iDoF) * elem(cell)%Grad_BF(iDoF,iGauss)%X
                     gradpGauss(2) = gradpGauss(2) + pDofLoc(iDoF) * elem(cell)%Grad_BF(iDoF,iGauss)%Y
#if MEF90_DIM == 3
                     gradpGauss(3) = gradpGauss(3) + pDofLoc(iDoF) * elem(cell)%Grad_BF(iDoF,iGauss)%Z
#endif
                  End Do
                  volumeWeight = elem(cell)%Gauss_C(iGauss)
                  gradpCell    = gradpCell + gradpGauss * volumeWeight
                  totalVolume  = totalVolume + volumeWeight
               End Do ! iGauss

               If (totalVolume > 0.0_Kr) Then
                  gradpCell = gradpCell / totalVolume
               End If

               !!! 写入 cell-centered Section — Restrict 获取指针，写值，再 Restore
               Call SectionRealRestrict(gradpSec,cellID(cell),secPtr,ierr);CHKERRQ(ierr)
               secPtr(1:dim) = gradpCell(1:dim)
               Call SectionRealRestore(gradpSec,cellID(cell),secPtr,ierr);CHKERRQ(ierr)
            End Do ! cell
            DeAllocate(pDofLoc)
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         End If

         Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr)
      End Do ! set

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      Call SectionRealComplete(gradpSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(gradpSec,FlowCtx%cellDMVectScatter,SCATTER_FORWARD,cellPressureGradient,ierr);CHKERRQ(ierr)

      Call SectionRealDestroy(gradpSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(pSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowComputePressureGradient


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  APPLY BIOT BODY FORCE  f_p = -α · ∇p  →  DefMechCtx%force
!!!  Uses VecAXPy so that any preset body forces (gravity, etc.) are preserved.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowApplyBiotForce"

   Subroutine MEF90FlowApplyBiotForce(FlowCtx,DefMechCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(IN)                :: FlowCtx
      Type(MEF90DefMechCtx_Type),Intent(INOUT)          :: DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      PetscReal                                          :: biotCoeff

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)
      biotCoeff = gOpt%biotCoefficient

      !!! α = 0 表示不启用 Biot 耦合，跳过
      If (biotCoeff == 0.0_Kr) Return

      !!! 计算 cell-centered ∇p
      Call MEF90FlowComputePressureGradient(FlowCtx%pressure, &
           FlowCtx%cellPressureGradient,FlowCtx,ierr);CHKERRQ(ierr)

      !!! 从预设 baseForce 重建，再叠加 Biot 体积力
      !!! f_p = -α · ∇p  →  force = baseForce + (-biotCoeff) * cellPressureGradient
      If (Associated(DefMechCtx%force)) Then
         !!! 先用 baseForce 重建 DefMechCtx%force（保留预设重力等）
         If (Associated(FlowCtx%baseForce)) Then
            Call VecCopy(FlowCtx%baseForce,DefMechCtx%force,ierr);CHKERRQ(ierr)
         Else
            Call VecSet(DefMechCtx%force,0.0_Kr,ierr);CHKERRQ(ierr)
         End If
         !!! 叠加 Biot 项：y = a*x + y
         Call VecAXPy(DefMechCtx%force,-biotCoeff,FlowCtx%cellPressureGradient,ierr);CHKERRQ(ierr)
      End If
   End Subroutine MEF90FlowApplyBiotForce


End Module MEF90_APPEND(m_MEF90_FlowAssembly,MEF90_DIM)D
