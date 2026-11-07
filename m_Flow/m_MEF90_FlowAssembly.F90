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
#define __FUNCT__ "MEF90FlowOperatorCore"

   Subroutine MEF90FlowOperatorCore(x,residual,FlowCtx,ierr)
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
      PetscReal,Dimension(:),Pointer                      :: permCell, sourceCell

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
            !!! CI #65: elem(k) is indexed by the k-th entry of the cellIS
            !!! passed to MEF90Element_Create. Every assembly call below uses
            !!! singleIS (size 1), so elem MUST be created from that same IS.
            !!! Creating it once from the whole setIS made all 600 cells reuse
            !!! elem(1)'s Grad_BF / Gauss_C (reproduced exactly by |K|inf in
            !!! CI #64). elem is now created/destroyed inside the cell loop.
            If (elemType%Codim == 0) Then
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr);CHKERRQ(ierr)
                  !!! CI #65: build a size-1 elem array from the SAME singleIS
                  !!! that the MEF90Diffusion*Set routines receive. Those
                  !!! routines index elem(k) by the position inside cellIS, so
                  !!! elem(1) now carries THIS cell's coordinates/geometry.
                  Call MEF90Element_Create(FlowCtx%DM,singleIS,elem,QuadratureOrder, &
                       cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)

                  !!! A = permField(cell) (k/μ scalar)
                  Call SectionRealRestrict(permFieldSec,cellID(cell),permCell,ierr);CHKERRQ(ierr)
                  !!! CI #65: an isotropic scalar permeability is R * Identity,
                  !!! i.e. MatS2D(R,R,0). MatS2D_Get_Real (used by the plain
                  !!! "A = real" assignment) sets XX=YY=XY=R, which is the
                  !!! singular rank-1 tensor [[R,R],[R,R]] and diffuses only
                  !!! along the (1,1) diagonal. That is the second half of the
                  !!! Gate 3.2 ratio error (CI #64 ratio=0.348 vs analytic
                  !!! 0.087). MEF90MatS2DIdentity is MatS2D(1,1,0).
                  A_cell = permCell(1) * MEF90_MATS_IDENTITY

                  !!! Gate 3.2 (CI #59): the pressure Vec/matrix/section live on
                  !!! DMScal (1 dof / vertex). The MAIN DM has a 2-dof default
                  !!! section (DefMech displacement), so passing FlowCtx%DM to
                  !!! these helpers makes DMmeshAssembleMatrix map a 3x3 element
                  !!! matrix through a 6-dof section -> wrong rows/cols, i.e.
                  !!! IFunction and IJacobian disagree -> exponential blow-up.
                  !!! HeatXfer uses MEF90HeatXferCtx%DMScal throughout; mirror it.
                  Call MEF90DiffusionOperatorSet(resSec,FlowCtx%DMScal,pSec,singleIS, &
                       A_cell,0.0_Kr,elem,elemType,ierr)

                  !!! Subtract source RHS: G = DiffusionOperator - Source
                  Call MEF90DiffusionRHSSetCell(resSec,FlowCtx%DMScal,sourceSec,singleIS, &
                       elem,elemType,ierr)

                  Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If

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
      Call MEF90FlowApplyDirichletResidual(x,residual,FlowCtx,ierr)

      Call VecAssemblyBegin(residual,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(residual,ierr);CHKERRQ(ierr)

      Call SectionRealDestroy(permFieldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(sourceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(resSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(pSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowOperatorCore


   !!!! SNES callback wrapper — matches PETSc SNESSetFunction signature
   Subroutine MEF90FlowOperator(snesHandle,x,residual,FlowCtx,ierr)
      Type(SNES),Intent(IN)                                :: snesHandle
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Call MEF90FlowOperatorCore(x,residual,FlowCtx,ierr)
   End Subroutine MEF90FlowOperator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SNES JACOBIAN:  K = dG/dp = Diffusion Bilinear Form (per-cell varying A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowBilinearFormCore"

   Subroutine MEF90FlowBilinearFormCore(x,Jac,PreJac,residualMode,FlowCtx,ierr)
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
      PetscReal,Dimension(:),Pointer                      :: permCell

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
            !!! CI #65: elem is built from singleIS inside the cell loop below,
            !!! NOT from the whole setIS (see MEF90FlowOperatorCore).
            If (elemType%Codim == 0) Then
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr);CHKERRQ(ierr)
                  !!! CI #65: build a size-1 elem array from the SAME singleIS
                  !!! that the MEF90Diffusion*Set routines receive. Those
                  !!! routines index elem(k) by the position inside cellIS, so
                  !!! elem(1) now carries THIS cell's coordinates/geometry.
                  Call MEF90Element_Create(FlowCtx%DM,singleIS,elem,QuadratureOrder, &
                       cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(permFieldSec,cellID(cell),permCell,ierr);CHKERRQ(ierr)
                  !!! CI #65: an isotropic scalar permeability is R * Identity,
                  !!! i.e. MatS2D(R,R,0). MatS2D_Get_Real (used by the plain
                  !!! "A = real" assignment) sets XX=YY=XY=R, which is the
                  !!! singular rank-1 tensor [[R,R],[R,R]] and diffuses only
                  !!! along the (1,1) diagonal. That is the second half of the
                  !!! Gate 3.2 ratio error (CI #64 ratio=0.348 vs analytic
                  !!! 0.087). MEF90MatS2DIdentity is MatS2D(1,1,0).
                  A_cell = permCell(1) * MEF90_MATS_IDENTITY
                  Call MEF90DiffusionBilinearFormSet(Jac,FlowCtx%DMScal,singleIS, &
                       A_cell,0.0_Kr,elem,elemType,ierr)
                  Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If
            Call ISDestroy(setIS,ierr)
         End If
      End Do

      Call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatCopy(Jac,PreJac,DIFFERENT_NONZERO_PATTERN,ierr);CHKERRQ(ierr)

      !!! CI #62: DO NOT zero-row here. In the transient path the storage
      !!! matrix (a*S*M) is added AFTER this call, which refills the BC rows
      !!! with mass entries; MEF90FlowIJacobian then had to zero-row AGAIN.
      !!! That gave 4 MatZeroRowsColumnsIS calls per step (Jac+PreJac, twice).
      !!! HeatXfer assembles K + shift*M FIRST and applies Dirichlet ONCE at
      !!! the very end. Mirror that: the callers now own BC enforcement.
      !!! (SNES steady path -> MEF90FlowBilinearForm, TS path -> IJacobian.)

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)
      Call SectionRealDestroy(permFieldSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowBilinearFormCore


   !!! SNES callback wrapper — matches PETSc 3.3 SNESSetJacobian signature:
   !!! (snes, x, A, B, MatStructure flg, ctx, ierr)
   Subroutine MEF90FlowBilinearForm(snesHandle,x,Jac,PreJac,flg,FlowCtx,ierr)
      Type(SNES),Intent(IN)                                :: snesHandle
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      !!! residualMode was a PetscReal dummy and is unused inside Core; the
      !!! assembly always uses the full Jacobian (scale 1.0).
      Call MEF90FlowBilinearFormCore(x,Jac,PreJac,0.0_Kr,FlowCtx,ierr);CHKERRQ(ierr)
      !!! CI #62: steady SNES path must enforce Dirichlet once, at the end.
      Call MEF90FlowApplyDirichletJacobian(Jac,PreJac,FlowCtx,ierr);CHKERRQ(ierr)
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90FlowBilinearForm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  TS RESIDUAL + JACOBIAN — reuse Operator/BilinearForm, add storage term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowIFunction"

   Subroutine MEF90FlowIFunction(tsHandle,t,x,xDot,F,FlowCtx,ierr)
      Type(TS),Intent(IN)                                :: tsHandle
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x, xDot
      Type(Vec),Intent(INOUT)                            :: F
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Character(len=256)                                 :: IOBuffer
      PetscReal                                          :: nX, nXDot, nF, nF1
      !!! CI #53 diagnostic: mark entry into the TS residual callback.
      Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: IFunction enter\n',ierr);CHKERRQ(ierr)
      !!! CI #64: SINGLE SCATTER, exactly like MEF90HeatXFerIFunction. The old
      !!! version performed THREE SectionRealToVec(SCATTER_FORWARD) into F
      !!! (OperatorCore, AddStorageToResidual, ApplyDirichletResidual). With
      !!! PETSc 3.3 TS_LINEAR + ksponly, xDot is identically 0 at the only
      !!! function evaluation (SNESTSFormFunction_Theta evaluates F at
      !!! X == X0, so Xdot = shift*(X-X0) = 0), hence the storage residual
      !!! is exactly zero and MEF90FlowOperatorCore already applies the
      !!! Dirichlet rows LAST. Both extra scatters were pure no-ops in the
      !!! best case and potential ghost-point accumulations in the worst case
      !!! (|F|inf cannot detect a partial accumulation). Reduce to one pass.
      Call MEF90FlowOperatorCore(x,F,FlowCtx,ierr);CHKERRQ(ierr)
      Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: IFunction post-OperatorCore\n',ierr);CHKERRQ(ierr)
      If (t < 0.07_Kr) Then
         Call VecNorm(x,NORM_INFINITY,nX,ierr);CHKERRQ(ierr)
         Call VecNorm(xDot,NORM_INFINITY,nXDot,ierr);CHKERRQ(ierr)
         Call VecNorm(F,NORM_INFINITY,nF,ierr);CHKERRQ(ierr)
         Call VecNorm(F,NORM_1,nF1,ierr);CHKERRQ(ierr)
         Write(IOBuffer,*) 'Gate 3.2 diag: t',t,' |x|',nX,' |xDot|',nXDot,' |F|inf',nF,' |F|1',nF1,'\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      End If
      Call VecAssemblyBegin(F,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(F,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowIFunction


#undef __FUNCT__
#define __FUNCT__ "MEF90FlowIJacobian"

   !!! Gate 3.1 fix: PETSc 3.3 requires the MatStructure output argument in
   !!! every TS IJacobian callback (cf. MEF90HeatXferIJacobian2D).
   Subroutine MEF90FlowIJacobian(tsHandle,t,x,xDot,a,Jac,PreJac,flg,FlowCtx,ierr)
      Type(TS),Intent(IN)                                :: tsHandle
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x, xDot
      PetscReal,Intent(IN)                               :: a
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Character(len=256)                                 :: IOBuffer
      PetscReal                                          :: sMin, sMax, kMin, kMax
      !!! CI #53 diagnostic: mark entry into the TS Jacobian callback.
      Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: IJacobian enter\n',ierr);CHKERRQ(ierr)
      !!! Gate 3.2 diagnostic: report the BE shift a = d(xdot)/d(x) for the
      !!! first few steps. Expected 1/dt = 50 for dt=0.02.
      If (t < 0.07_Kr) Then
         Write(IOBuffer,*) 'Gate 3.2 diag: IJacobian t',t,' shift a',a,'\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Call VecMin(FlowCtx%storativity,PETSC_NULL_INTEGER,sMin,ierr);CHKERRQ(ierr)
         Call VecMax(FlowCtx%storativity,PETSC_NULL_INTEGER,sMax,ierr);CHKERRQ(ierr)
         Call VecMin(FlowCtx%permField,PETSC_NULL_INTEGER,kMin,ierr);CHKERRQ(ierr)
         Call VecMax(FlowCtx%permField,PETSC_NULL_INTEGER,kMax,ierr);CHKERRQ(ierr)
         Write(IOBuffer,*) 'Gate 3.2 diagFld: t',t,' S min/max',sMin,sMax,' k/mu min/max',kMin,kMax,'\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      End If
      !!! First, assemble steady Darcy Jacobian K = dG/dp (no BC yet).
      Call MEF90FlowBilinearFormCore(x,Jac,PreJac,0.0_Kr,FlowCtx,ierr);CHKERRQ(ierr)
      If (t < 0.07_Kr) Then
         Call MEF90FlowReportJacobian('post-Core',t,Jac,x,ierr);CHKERRQ(ierr)
      End If
      !!! Add a·S·M (mass matrix):  K += a·S·BFB
      Call MEF90FlowAddStorageMassMatrix(Jac,PreJac,a,FlowCtx,ierr);CHKERRQ(ierr)
      If (t < 0.07_Kr) Then
         Call MEF90FlowReportJacobian('post-Storage',t,Jac,x,ierr);CHKERRQ(ierr)
      End If
      !!! CI #62: Dirichlet is applied ONCE, after every matrix contribution,
      !!! exactly like MEF90HeatXFerIJacobian. J(dof,:) = 0, J(dof,dof) = 1,
      !!! consistent with the residual F(dof) = x(dof) - p_bc.
      Call MEF90FlowApplyDirichletJacobian(Jac,PreJac,FlowCtx,ierr);CHKERRQ(ierr)
      If (t < 0.07_Kr) Then
         Call MEF90FlowReportJacobian('post-Dirichlet',t,Jac,x,ierr);CHKERRQ(ierr)
      End If
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90FlowIJacobian


#undef __FUNCT__
#define __FUNCT__ "MEF90FlowReportJacobian"
!!! CI #62 diagnostic: report |J|inf and the diagonal range at each assembly
!!! checkpoint. If K is present, |J|inf ~ 1e-7 and diag min ~ 1e-8, max = 1
!!! (the Dirichlet rows). If the BC IS wiped the whole matrix, |J|inf ~ 1 and
!!! every diagonal entry = 1 -> J is the identity and CG cannot diffuse.
   Subroutine MEF90FlowReportJacobian(tag,t,Jac,x,ierr)
      Character(len=*),Intent(IN)                        :: tag
      PetscReal,Intent(IN)                               :: t
      Type(Mat),Intent(IN)                               :: Jac
      Type(Vec),Intent(IN)                               :: x
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(Vec)                                          :: diagVec
      PetscReal                                          :: nJ, dMin, dMax
      Character(len=256)                                 :: IOBuffer

      Call VecDuplicate(x,diagVec,ierr);CHKERRQ(ierr)
      Call MatNorm(Jac,NORM_INFINITY,nJ,ierr);CHKERRQ(ierr)
      Call MatGetDiagonal(Jac,diagVec,ierr);CHKERRQ(ierr)
      Call VecMin(diagVec,PETSC_NULL_INTEGER,dMin,ierr);CHKERRQ(ierr)
      Call VecMax(diagVec,PETSC_NULL_INTEGER,dMax,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) 'Gate 3.2 diagJ: t',t,' ',Trim(tag),' |J|inf',nJ,' diag min/max',dMin,dMax,'\n'
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call VecDestroy(diagVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowReportJacobian


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
      PetscReal,Dimension(:),Pointer                      :: permCell
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
            !!! CI #65: elem(k) is indexed by the k-th entry of the cellIS
            !!! passed to MEF90Element_Create. Every assembly call below uses
            !!! singleIS (size 1), so elem MUST be created from that same IS.
            !!! Creating it once from the whole setIS made all 600 cells reuse
            !!! elem(1)'s Grad_BF / Gauss_C (reproduced exactly by |K|inf in
            !!! CI #64). elem is now created/destroyed inside the cell loop.
            If (elemType%Codim == 0) Then
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr);CHKERRQ(ierr)
                  !!! CI #65: build a size-1 elem array from the SAME singleIS
                  !!! that the MEF90Diffusion*Set routines receive. Those
                  !!! routines index elem(k) by the position inside cellIS, so
                  !!! elem(1) now carries THIS cell's coordinates/geometry.
                  Call MEF90Element_Create(FlowCtx%DM,singleIS,elem,QuadratureOrder, &
                       cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(permFieldSec,cellID(cell),permCell,ierr);CHKERRQ(ierr)
                  !!! CI #65: an isotropic scalar permeability is R * Identity,
                  !!! i.e. MatS2D(R,R,0). MatS2D_Get_Real (used by the plain
                  !!! "A = real" assignment) sets XX=YY=XY=R, which is the
                  !!! singular rank-1 tensor [[R,R],[R,R]] and diffuses only
                  !!! along the (1,1) diagonal. That is the second half of the
                  !!! Gate 3.2 ratio error (CI #64 ratio=0.348 vs analytic
                  !!! 0.087). MEF90MatS2DIdentity is MatS2D(1,1,0).
                  A_cell = permCell(1) * MEF90_MATS_IDENTITY

                  Call MEF90DiffusionEnergySet(energyCell,pSec,p0Sec,FlowCtx%DM, &
                       A_cell,0.0_Kr,singleIS,elem,elemType,ierr)
                  Energy = Energy + energyCell

                  Call MEF90DiffusionWorkSetCell(workCell,pSec,FlowCtx%DM, &
                       sourceSec,singleIS,elem,elemType,ierr)
                  InjectionWork = InjectionWork + workCell

                  Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If
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

      !!! Cell-centered diagnostics: permeability field k(v)/mu and Darcy velocity vector
      If (gOpt%permFieldOffset > 0 .AND. Associated(FlowCtx%permField)) Then
         Call VecViewExodusCell(FlowCtx%cellDMScal,FlowCtx%permField,FlowCtx%MEF90Ctx%IOcomm, &
                                FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%permFieldOffset,ierr);CHKERRQ(ierr)
      End If
      If (gOpt%darcyVelocityOffset > 0 .AND. Associated(FlowCtx%cellDarcyVelocity)) Then
         Call VecViewExodusCell(FlowCtx%cellDMVect,FlowCtx%cellDarcyVelocity,FlowCtx%MEF90Ctx%IOcomm, &
                                FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%darcyVelocityOffset,ierr);CHKERRQ(ierr)
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
   Subroutine MEF90FlowApplyDirichletResidual(x,residual,FlowCtx,ierr)
      Type(Vec),Intent(INOUT)                             :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                  :: FlowCtx
      PetscErrorCode,Intent(OUT)                           :: ierr
      Type(Vec),Intent(IN)                                :: x

      Type(IS)                                            :: CellSetGlobalIS, bcIS
      Type(IS)                                            :: VertexSetGlobalIS, setIS, cellSetIS
      PetscInt,dimension(:),Pointer                        :: setID, pointIDs
      PetscInt                                            :: set, nVal, v
      Type(SectionReal)                                   :: resSec, xSec, bcSec
      PetscReal,Dimension(:),Pointer                      :: rPtr, xPtr, bPtr
      Type(MEF90FlowCellSetOptions_Type),pointer          :: cellSetOptions
      Type(MEF90FlowVertexSetOptions_Type),pointer        :: vertexSetOptions

      !!! CI #57 ROOT CAUSE: MEF90VecGetValuesISdof / MEF90VecSetValuesISdof
      !!! contain COLLECTIVE calls (DMMeshCreateGlobalScatter +
      !!! SectionRealToVec). Guarding them with "If (nVal > 0)" made rank 1
      !!! (owns every top-surface BC vertex) enter a collective that rank 0
      !!! (nVal = 0) skipped -> MPI deadlock, i.e. the #56 SEGV turned into a
      !!! 1h hang. Follow the SectionReal pattern already proven on both ranks
      !!! in MEF90FlowOperatorCore / MEF90FlowAddStorageToResidual: every
      !!! collective runs UNCONDITIONALLY on both ranks, only the purely local
      !!! per-point writes are guarded by nVal.

      !!! Cell sets with Has_BC → Cell→Vertex IS
      Call SectionRealDuplicate(FlowCtx%DMScalSec,resSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%DMScalSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%DMScalSec,bcSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(resSec,FlowCtx%DMScalScatter,SCATTER_REVERSE,residual,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,FlowCtx%DMScalScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(bcSec,FlowCtx%DMScalScatter,SCATTER_REVERSE,FlowCtx%boundaryPressure,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),cellSetIS,ierr);CHKERRQ(ierr)
            Call MEF90ISCreateCelltoVertex(FlowCtx%DMScal,PETSC_COMM_WORLD,cellSetIS,bcIS,ierr);CHKERRQ(ierr)
            Call ISGetLocalSize(bcIS,nVal,ierr);CHKERRQ(ierr)
            !!! Purely local: a rank may own zero BC vertices.
            If (nVal > 0) Then
               Call ISGetIndicesF90(bcIS,pointIDs,ierr);CHKERRQ(ierr)
               Do v = 1, nVal
                  Call SectionRealRestrict(resSec,pointIDs(v),rPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(xSec,pointIDs(v),xPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(bcSec,pointIDs(v),bPtr,ierr);CHKERRQ(ierr)
                  rPtr(1) = xPtr(1) - bPtr(1)   !!! G = x - p_bc on Dirichlet dofs
                  Call SectionRealRestore(resSec,pointIDs(v),rPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(xSec,pointIDs(v),xPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(bcSec,pointIDs(v),bPtr,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(bcIS,pointIDs,ierr);CHKERRQ(ierr)
            End If
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(cellSetIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      !!! Vertex sets with Has_BC → override
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxVertexSetOptions(FlowCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DMScal,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            Call ISGetLocalSize(setIS,nVal,ierr);CHKERRQ(ierr)
            If (nVal > 0) Then
               Call ISGetIndicesF90(setIS,pointIDs,ierr);CHKERRQ(ierr)
               Do v = 1, nVal
                  Call SectionRealRestrict(resSec,pointIDs(v),rPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(xSec,pointIDs(v),xPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(bcSec,pointIDs(v),bPtr,ierr);CHKERRQ(ierr)
                  rPtr(1) = xPtr(1) - bPtr(1)
                  Call SectionRealRestore(resSec,pointIDs(v),rPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(xSec,pointIDs(v),xPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(bcSec,pointIDs(v),bPtr,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,pointIDs,ierr);CHKERRQ(ierr)
            End If
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

      !!! Collective scatter back — executed on EVERY rank, no nVal guard.
      Call SectionRealComplete(resSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(resSec,FlowCtx%DMScalScatter,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(resSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(bcSec,ierr);CHKERRQ(ierr)

      Call VecAssemblyBegin(residual,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(residual,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowApplyDirichletResidual


   Subroutine MEF90FlowApplyDirichletJacobian(Jac,PreJac,FlowCtx,ierr)
      Type(Mat),Intent(INOUT)                             :: Jac, PreJac
      Type(MEF90FlowCtx_Type),Intent(IN)                  :: FlowCtx
      PetscErrorCode,Intent(OUT)                           :: ierr

      Type(IS)                                            :: CellSetGlobalIS, bcIS, setISdof
      Type(IS)                                            :: VertexSetGlobalIS, setIS, cellIS
      PetscInt,dimension(:),Pointer                        :: setID
      PetscInt                                            :: set, nVal
      Type(MEF90FlowCellSetOptions_Type),pointer          :: cellSetOptions
      Type(MEF90FlowVertexSetOptions_Type),pointer        :: vertexSetOptions

      !!! CI #56: point-ID IS -> global dof IS, exactly as in
      !!! MEF90HeatXferBilinearForm / MEF90HeatXferIJacobian. MatZeroRows
      !!! used to receive raw mesh point IDs as matrix row indices.
      !!! DIRICHLET JACOBIAN: K(i,i) = 1 for all BC dofs.

      !!! Cell sets with Has_BC
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_BC) Then
            !!! cellIS = cell-set StratumIS, bcIS = Cell→Vertex POINT IS,
            !!! setISdof = global dof IS (the only valid MatZeroRows input).
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),cellIS,ierr);CHKERRQ(ierr)
            Call MEF90ISCreateCelltoVertex(FlowCtx%DMScal,PETSC_COMM_WORLD,cellIS,bcIS,ierr);CHKERRQ(ierr)
            Call ISGetLocalSize(bcIS,nVal,ierr);CHKERRQ(ierr)
            !!! CI #57: NO "If (nVal > 0)" guard here. MatZeroRowsColumnsIS is
            !!! collective; rank 0 owns zero BC vertices but MUST still call it
            !!! with an empty dof IS, otherwise rank 1 blocks forever.
            Call DMMeshISCreateISglobaldof(FlowCtx%DMScal,bcIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(Jac,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(PreJac,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      !!! Vertex sets with Has_BC
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxVertexSetOptions(FlowCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(FlowCtx%DMScal,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            Call ISGetLocalSize(setIS,nVal,ierr);CHKERRQ(ierr)
            !!! CI #57: both ranks must enter these collective calls (see above).
            Call DMMeshISCreateISglobaldof(FlowCtx%DMScal,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(Jac,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(PreJac,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

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
      Character(len=256)                                 :: IOBuffer

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) 'Gate 3.1: SetTransients enter step',step,' srcScale',gOpt%sourceScaling,' bcScale',gOpt%boundaryPressureScaling,'\n'
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

      Select case (gOpt%sourceScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(FlowCtx%cellDMScal,FlowCtx%source,FlowCtx%MEF90Ctx%IOcomm, &
                                    FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%sourceOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90FlowSetSourceCst(FlowCtx%source,FlowCtx,ierr);CHKERRQ(ierr)
            Call VecScale(FlowCtx%source,timeVal,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: SetTransients pre-SetSourceCst\n',ierr);CHKERRQ(ierr)
            Call MEF90FlowSetSourceCst(FlowCtx%source,FlowCtx,ierr);CHKERRQ(ierr)
            Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: SetTransients post-SetSourceCst\n',ierr);CHKERRQ(ierr)
      End Select

      Select case (gOpt%boundaryPressureScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusVertex(FlowCtx%DMScal,FlowCtx%boundaryPressure,FlowCtx%MEF90Ctx%IOcomm, &
                                     FlowCtx%MEF90Ctx%fileExoUnit,step,gOpt%boundaryPressureOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90FlowUpdateboundaryPressure(FlowCtx%boundaryPressure,FlowCtx,ierr);CHKERRQ(ierr)
            Call VecScale(FlowCtx%boundaryPressure,timeVal,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: SetTransients pre-UpdateBC\n',ierr);CHKERRQ(ierr)
            Call MEF90FlowUpdateboundaryPressure(FlowCtx%boundaryPressure,FlowCtx,ierr);CHKERRQ(ierr)
            Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: SetTransients post-UpdateBC\n',ierr);CHKERRQ(ierr)
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
      PetscInt,dimension(:),Pointer                      :: setID, pointIDs
      PetscInt                                           :: set, nVal, v
      Type(SectionReal)                                  :: bcSec
      PetscReal,Dimension(:),Pointer                     :: Ptr
      Character(len=256)                                 :: IOBuffer

      !!! CI #57: SectionReal path (same as the proven SetSourceCst). All
      !!! collective operations (Complete / ToVec) run unconditionally on
      !!! BOTH ranks at the end of this routine; only the local per-point
      !!! writes are guarded by nVal. No DMMeshISCreateISglobaldof, no
      !!! collective call inside an nVal-dependent branch -> no deadlock.
      Call VecSet(boundaryPressure,0.0_Kr,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%DMScalSec,bcSec,ierr);CHKERRQ(ierr)
      Call SectionRealSet(bcSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call PetscPrintf(PETSC_COMM_SELF,'Gate 3.1: UpdateBC enter\n',ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_BC) Then
            Write(IOBuffer,*) 'Gate 3.1: UpdateBC cell set',set,' id',setID(set),'\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
            Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),cellSetIS,ierr);CHKERRQ(ierr)
            !!! CI #54: this call used to have NO CHKERRQ, so a bad/unset bcIS
            !!! was silently passed to ISGetLocalSize and segfaulted on rank 1.
            Call MEF90ISCreateCelltoVertex(FlowCtx%DMScal,PETSC_COMM_WORLD,cellSetIS,bcIS,ierr);CHKERRQ(ierr)
            Call ISGetLocalSize(bcIS,nVal,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.1: UpdateBC cell set',set,' nVal',nVal,'\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
            !!! A rank may own zero vertices of a boundary set. Purely local
            !!! per-point writes, safe to guard (CI #55 rank-1 pattern).
            If (nVal > 0) Then
               Call ISGetIndicesF90(bcIS,pointIDs,ierr);CHKERRQ(ierr)
               Do v = 1, nVal
                  Call SectionRealRestrict(bcSec,pointIDs(v),Ptr,ierr);CHKERRQ(ierr)
                  Ptr(1) = cellSetOptions%boundaryPressure
                  Call SectionRealRestore(bcSec,pointIDs(v),Ptr,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(bcIS,pointIDs,ierr);CHKERRQ(ierr)
            End If
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(cellSetIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxVertexSetOptions(FlowCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Write(IOBuffer,*) 'Gate 3.1: UpdateBC vertex set',set,' id',setID(set),'\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
            Call DMMeshGetStratumIS(FlowCtx%DMScal,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            Call ISGetLocalSize(setIS,nVal,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.1: UpdateBC vertex set',set,' nVal',nVal,'\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
            If (nVal > 0) Then
               Call ISGetIndicesF90(setIS,pointIDs,ierr);CHKERRQ(ierr)
               Do v = 1, nVal
                  Call SectionRealRestrict(bcSec,pointIDs(v),Ptr,ierr);CHKERRQ(ierr)
                  Ptr(1) = vertexSetOptions%boundaryPressure
                  Call SectionRealRestore(bcSec,pointIDs(v),Ptr,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,pointIDs,ierr);CHKERRQ(ierr)
            End If
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

      !!! Collective scatter back — executed on EVERY rank, no nVal guard.
      Call SectionRealComplete(bcSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(bcSec,FlowCtx%DMScalScatter,SCATTER_FORWARD,boundaryPressure,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(bcSec,ierr);CHKERRQ(ierr)

      Call VecAssemblyBegin(boundaryPressure,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(boundaryPressure,ierr);CHKERRQ(ierr)
      Call PetscPrintf(PETSC_COMM_SELF,'Gate 3.1: UpdateBC done\n',ierr);CHKERRQ(ierr)
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
      Type(MEF90Element_Type)                            :: elemType
      Type(SectionReal)                                  :: sourceSec
      Type(IS)                                           :: CellSetGlobalIS, setIS
      PetscInt,dimension(:),Pointer                      :: setID, cellIDs
      PetscInt                                           :: set, cell, nCell
      PetscReal,Dimension(:),Pointer                     :: Sptr
      Character(len=256)                                 :: IOBuffer

      !!! CI #55 ROOT CAUSE: the previous implementation used
      !!! DMMeshISCreateISglobaldof(cellDMScal,...) + VecSetValues. On rank 1
      !!! the global-dof IS is empty / out of the cell Vec layout in PETSc 3.3
      !!! and ISGetIndicesF90 on it segfaults right after the "set 1 id 1"
      !!! trace. MEF90FlowUpdatePerm / UpdateStorativity (init 07/08 and every
      !!! step's post-UpdatePerm) already prove the SectionReal path is safe
      !!! on both ranks: write into a duplicate of cellDMScalSec, then
      !!! SectionRealToVec(SCATTER_FORWARD) into the cell-centered Vec.
      Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: SetSourceCst enter (SectionReal path)\n',ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,sourceSec,ierr);CHKERRQ(ierr)
      Call SectionRealSet(sourceSec,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         !!! cs0010/cs0020 are BAR2 boundary facets: Has_BC=T and no dof in the
         !!! cell section. OperatorCore and the storage terms skip them too.
         If (cellSetOptions%Has_BC) Cycle
         elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
         If (elemType%Codim /= 0) Cycle
         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setIS,cellIDs,ierr);CHKERRQ(ierr)
         nCell = size(cellIDs)
         Write(IOBuffer,*) 'Gate 3.1: SetSourceCst set',set,' id',setID(set),' nCell',nCell,' src',cellSetOptions%source,'\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Do cell = 1, nCell
            Call SectionRealRestrict(sourceSec,cellIDs(cell),Sptr,ierr);CHKERRQ(ierr)
            Sptr(1) = cellSetOptions%source
            Call SectionRealRestore(sourceSec,cellIDs(cell),Sptr,ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(setIS,cellIDs,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      Call SectionRealComplete(sourceSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(sourceSec,FlowCtx%cellDMScalScatter,SCATTER_FORWARD,sourceVec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(sourceSec,ierr);CHKERRQ(ierr)
      Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: SetSourceCst done\n',ierr);CHKERRQ(ierr)
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
      PetscReal,Dimension(:),Pointer                      :: SCell
      PetscReal                                           :: S_avg
      Type(MEF90_MATS)                                    :: A_zero

      Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: storage residual enter\n',ierr);CHKERRQ(ierr)
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
            !!! CI #65: elem is built from singleIS inside the cell loop below,
            !!! NOT from the whole setIS (see MEF90FlowOperatorCore).
            If (elemType%Codim == 0) Then
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr);CHKERRQ(ierr)
                  !!! CI #65: build a size-1 elem array from the SAME singleIS
                  !!! that the MEF90Diffusion*Set routines receive. Those
                  !!! routines index elem(k) by the position inside cellIS, so
                  !!! elem(1) now carries THIS cell's coordinates/geometry.
                  Call MEF90Element_Create(FlowCtx%DM,singleIS,elem,QuadratureOrder, &
                       cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
                  !!! alpha = S(cell), pre-multiply and pass to DiffusionOperatorAddTransientTermSet
                  Call SectionRealRestrict(storativitySec,cellID(cell),SCell,ierr);CHKERRQ(ierr)
                  !!! CI #53: a rank may own zero cells in a set; guard the
                  !!! restricted pointer before reading SCell(1).
                  If (Associated(SCell)) Then
                     If (size(SCell) > 0) Then
                        S_avg = SCell(1)
                     Else
                        S_avg = 0.0_Kr
                     End If
                  Else
                     S_avg = 0.0_Kr
                  End If
                  Call MEF90DiffusionOperatorAddTransientTermSet(FSec,FlowCtx%DMScal,pDotSec, &
                       singleIS,S_avg,elem,elemType,ierr)
                  Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If
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
      PetscReal,Dimension(:),Pointer                      :: SCell
      PetscReal                                           :: lambda_val
      Type(MEF90_MATS)                                    :: A_zero
      Type(SectionReal)                                   :: storativitySec

      Call PetscPrintf(PETSC_COMM_WORLD,'Gate 3.1: storage matrix enter\n',ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,storativitySec,ierr);CHKERRQ(ierr)
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
            !!! CI #65: elem is built from singleIS inside the cell loop below,
            !!! NOT from the whole setIS (see MEF90FlowOperatorCore).
            If (elemType%Codim == 0) Then
               A_zero = 0.0_Kr
               Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
               Do cell = 1, size(cellID)
                  Call MEF90FlowMakeSingleCellIS(FlowCtx%DM,cellID(cell),singleIS,ierr);CHKERRQ(ierr)
                  !!! CI #65: build a size-1 elem array from the SAME singleIS
                  !!! that the MEF90Diffusion*Set routines receive. Those
                  !!! routines index elem(k) by the position inside cellIS, so
                  !!! elem(1) now carries THIS cell's coordinates/geometry.
                  Call MEF90Element_Create(FlowCtx%DM,singleIS,elem,QuadratureOrder, &
                       cellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(storativitySec,cellID(cell),SCell,ierr);CHKERRQ(ierr)
                  If (Associated(SCell)) Then
                     If (size(SCell) > 0) Then
                        lambda_val = a * SCell(1)
                     Else
                        lambda_val = 0.0_Kr
                     End If
                  Else
                     lambda_val = 0.0_Kr
                  End If
                  !!! K += 0*∇ + λ·BFB   (lambda term only → mass matrix)
                  Call MEF90DiffusionBilinearFormSet(Jac,FlowCtx%DMScal,singleIS, &
                       A_zero,lambda_val,elem,elemType,ierr)
                  Call MEF90DiffusionBilinearFormSet(PreJac,FlowCtx%DMScal,singleIS, &
                       A_zero,lambda_val,elem,elemType,ierr)
                  Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
                  Call ISDestroy(singleIS,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            End If
            Call ISDestroy(setIS,ierr)
         End If
      End Do

      Call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyBegin(PreJac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd(PreJac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)
      Call SectionRealDestroy(storativitySec,ierr);CHKERRQ(ierr)
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
