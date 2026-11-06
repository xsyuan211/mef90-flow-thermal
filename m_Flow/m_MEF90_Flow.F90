#include "../MEF90/mef90.inc"
Module m_MEF90_Flow
#include "finclude/petscdef.h"
   Use petsc
   Use m_MEF90
   Use m_MEF90_FlowCtx
   Use m_MEF90_DefMechCtx     !!! NEW: for MEF90DefMechCtx_Type reference

   Use m_MEF90_FlowAssembly2D, &   !!! NEW: bridge to 2D dimmed assembly template
      MEF90FlowOperator2D               => MEF90FlowOperator, &
      MEF90FlowBilinearForm2D           => MEF90FlowBilinearForm, &
      MEF90FlowIFunction2D              => MEF90FlowIFunction, &
      MEF90FlowIJacobian2D              => MEF90FlowIJacobian, &
      MEF90FlowEnergy2D                 => MEF90FlowEnergy, &
      MEF90FlowViewEXO2D                => MEF90FlowViewEXO, &
      MEF90FlowUpdateboundaryPressure2D => MEF90FlowUpdateboundaryPressure, &
      MEF90FlowSetTransients2D          => MEF90FlowSetTransients, &
      MEF90FlowComputePressureGradient2D => MEF90FlowComputePressureGradient, &
      MEF90FlowApplyBiotForce2D         => MEF90FlowApplyBiotForce
   Use m_MEF90_FlowAssembly3D, &   !!! NEW: bridge to 3D dimmed assembly template
      MEF90FlowOperator3D               => MEF90FlowOperator, &
      MEF90FlowBilinearForm3D           => MEF90FlowBilinearForm, &
      MEF90FlowIFunction3D              => MEF90FlowIFunction, &
      MEF90FlowIJacobian3D              => MEF90FlowIJacobian, &
      MEF90FlowEnergy3D                 => MEF90FlowEnergy, &
      MEF90FlowViewEXO3D                => MEF90FlowViewEXO, &
      MEF90FlowUpdateboundaryPressure3D => MEF90FlowUpdateboundaryPressure, &
      MEF90FlowSetTransients3D          => MEF90FlowSetTransients, &
      MEF90FlowComputePressureGradient3D => MEF90FlowComputePressureGradient, &
      MEF90FlowApplyBiotForce3D         => MEF90FlowApplyBiotForce

   Implicit none

   Private
   !!! Public — Context lifecycle (defined in FlowCtx.F90 or here)
   Public :: MEF90FlowCtxSetFromOptions,     &
             MEF90FlowCtxCreateVectors,      &
             MEF90FlowCtxSetSections,        &
             MEF90FlowCtxDestroyVectors,     &
             MEF90FlowCtxDestroy,            &
             MEF90FlowUpdatePerm,            &   !!! k(v) = k_min + (k_max - k_min)*damage
             MEF90FlowUpdateStorativity,     &   !!! S per cell from options
             MEF90FlowUpdateMuFromTemperature,&   !!! NEW: μ(T) temperature-dependent viscosity
             MEF90FlowComputeDarcyVelocity,  &   !!! NEW: v_darcy = -(k/μ)·∇p + body force
             MEF90FlowCreateSNES,            &   !!! SNES nonlinear pressure solver
             MEF90FlowCreateTS,              &   !!! TS transient pressure solver
             MEF90FlowSetTransients,         &   !!! time-dependent source/BC update
             MEF90FlowUpdateboundaryPressure,&
             MEF90FlowOperator,              &   !!! Assembled via 2D/3D dispatch below
             MEF90FlowBilinearForm,          &
             MEF90FlowIFunction,             &
             MEF90FlowIJacobian,             &
             MEF90FlowEnergy,                &
             MEF90FlowViewEXO,               &
             MEF90FlowComputePressureGradient,& !!! Cell-centered ∇p volume average
             MEF90FlowApplyBiotForce           !!! Biot f_p = -α·∇p → DefMechCtx%force

Contains

!!! --------------------------------------------------------------------------------
!!!  Wrapper functions — dispatch 2D vs 3D based on mesh dimension (HeatXfer pattern)
!!! --------------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowOperator"
   Subroutine MEF90FlowOperator(snesPressure,x,residual,FlowCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesPressure
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowOperator2D(snesPressure,x,residual,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowOperator3D(snesPressure,x,residual,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowOperator

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowBilinearForm"
   !!! Gate 3.1 fix: PETSc 3.3 SNESSetJacobian callbacks take a MatStructure
   !!! output argument, not a PetscReal (cf. MEF90DefMechBilinearFormDamage
   !!! at m_DefMech/m_MEF90_DefMech.F90:1043). The former residualMode dummy
   !!! was reading the enum as a real and never setting flg.
   Subroutine MEF90FlowBilinearForm(snesPressure,x,Jac,PreJac,flg,FlowCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesPressure
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowBilinearForm2D(snesPressure,x,Jac,PreJac,flg,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowBilinearForm3D(snesPressure,x,Jac,PreJac,flg,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowIFunction"
   Subroutine MEF90FlowIFunction(tsHandle,t,x,xDot,F,FlowCtx,ierr)
      Type(TS),Intent(IN)                                :: tsHandle
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x, xDot
      Type(Vec),Intent(INOUT)                            :: F
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowIFunction2D(tsHandle,t,x,xDot,F,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowIFunction3D(tsHandle,t,x,xDot,F,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowIFunction

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowIJacobian"
   !!! Gate 3.1 fix: PETSc 3.3 TS IJacobian callbacks carry a MatStructure
   !!! output argument (see MEF90HeatXferIJacobian at
   !!! m_HeatXfer/m_MEF90_HeatXfer.F90:464). Omitting it makes TS call this
   !!! routine with one extra argument -> stack corruption at run time.
   Subroutine MEF90FlowIJacobian(tsHandle,t,x,xDot,a,Jac,PreJac,flg,FlowCtx,ierr)
      Type(TS),Intent(IN)                                :: tsHandle
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x, xDot
      PetscReal,Intent(IN)                               :: a
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowIJacobian2D(tsHandle,t,x,xDot,a,Jac,PreJac,flg,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowIJacobian3D(tsHandle,t,x,xDot,a,Jac,PreJac,flg,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowIJacobian

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowEnergy"
   Subroutine MEF90FlowEnergy(pressure,timeVal,FlowCtx,Energy,InjectionWork,ierr)
      Type(Vec),Intent(IN)                               :: pressure
      PetscReal,Intent(IN)                               :: timeVal
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscReal,Intent(OUT)                              :: Energy, InjectionWork
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowEnergy2D(pressure,timeVal,FlowCtx,Energy,InjectionWork,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowEnergy3D(pressure,timeVal,FlowCtx,Energy,InjectionWork,ierr)
      End If
   End Subroutine MEF90FlowEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowViewEXO"
   Subroutine MEF90FlowViewEXO(FlowCtx,step,ierr)
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowViewEXO2D(FlowCtx,step,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowViewEXO3D(FlowCtx,step,ierr)
      End If
   End Subroutine MEF90FlowViewEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowUpdateboundaryPressure"
   Subroutine MEF90FlowUpdateboundaryPressure(boundaryPressure,FlowCtx,ierr)
      Type(Vec),Intent(INOUT)                            :: boundaryPressure
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowUpdateboundaryPressure2D(boundaryPressure,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowUpdateboundaryPressure3D(boundaryPressure,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowUpdateboundaryPressure

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowSetTransients"
   Subroutine MEF90FlowSetTransients(FlowCtx,step,timeVal,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      PetscInt,Intent(IN)                                :: step
      PetscReal,Intent(IN)                               :: timeVal
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowSetTransients2D(FlowCtx,step,timeVal,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowSetTransients3D(FlowCtx,step,timeVal,ierr)
      End If
   End Subroutine MEF90FlowSetTransients

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Update permeability field from damage field and global options
!!!  k(v) = k_min + (k_max - k_min) * max(v, residual_damage)
!!!  permField(cell) = k(v) / mu   (cell-centered, averaged from damage)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowUpdatePerm"

   Subroutine MEF90FlowUpdatePerm(FlowCtx,DefMechCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: DefMechCtx
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellOpt
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer     :: elem
      Type(MEF90Element_Type)                             :: elemType
      Type(IS)                                            :: CellSetGlobalIS, setIS
      PetscInt,dimension(:),Pointer                       :: setIdx
      PetscInt,dimension(:),Pointer                       :: setID
      PetscInt                                            :: set, cell
      PetscReal                                           :: damageGauss, vGauss, k_v, k_eff
      PetscReal,Dimension(:),Pointer                      :: dmgPtr
      PetscReal,Dimension(:),Pointer                     :: permPtr, muPtr
      PetscInt                                            :: iDoF, iGauss, QuadratureOrder

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellOpt,ierr);CHKERRQ(ierr)

         !!! Base k for this cell set — use override if given, else k_min
         k_eff = gOpt%permeabilityMin
         If (cellOpt%permeabilityOverride > 0.0_Kr) Then
            k_eff = cellOpt%permeabilityOverride
         End If

         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)

         elemType = MEF90knownElements(cellOpt%ElemTypeShortID)
         QuadratureOrder = elemType%order * 2
         Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
              cellOpt%ElemTypeShortID,ierr);CHKERRQ(ierr)

         Allocate(dmgPtr(elemType%numDof))
         Do cell = 1, size(setIdx)
            !!! Restrict damage to this cell's DOFs
            Call SectionRealRestrictClosure(DefMechCtx%DMScalSec,FlowCtx%DM,setIdx(cell),elemType%numDof,dmgPtr,ierr);CHKERRQ(ierr)

            !!! Average damage at cell (Gauss-point averaging of v-field shape)
            vGauss = 0.0_Kr
            Do iGauss = 1, size(elem(cell)%Gauss_C)
               damageGauss = 0.0_Kr
               Do iDoF = 1, elemType%numDof
                  damageGauss = damageGauss + dmgPtr(iDoF) * elem(cell)%BF(iDoF,iGauss)
               End Do
               vGauss = vGauss + elem(cell)%Gauss_C(iGauss) * damageGauss
            End Do
            vGauss = vGauss / Sum(elem(cell)%Gauss_C)

            !!! Apply residual damage floor
            If (vGauss < gOpt%residualDamage) vGauss = gOpt%residualDamage

            !!! k(v) = k_min + (k_max - k_min)*v
            k_v = k_eff + (gOpt%permeabilityMax - k_eff) * vGauss

            !!! SectionRealRestrict → get writeable ptr, then write perm = k/μ
            Call SectionRealRestrict(FlowCtx%cellDMScalSec,setIdx(cell),permPtr,ierr);CHKERRQ(ierr)

            !!! Use temperature-dependent viscosity μ(T) if enabled, else constant μ_ref
            If (gOpt%temperatureDependViscosity .AND. Associated(FlowCtx%cellMuField)) Then
               Call SectionRealRestrict(FlowCtx%cellDMScalSec,setIdx(cell),muPtr,ierr);CHKERRQ(ierr)
               permPtr(1) = k_v / max(muPtr(1),1e-12_Kr)
            Else
               permPtr(1) = k_v / max(gOpt%dynamicViscosity,1e-12_Kr)
            End If

            Call SectionRealRestore(FlowCtx%cellDMScalSec,setIdx(cell),permPtr,ierr);CHKERRQ(ierr)
         End Do
         DeAllocate(dmgPtr)
         Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do

      Call SectionRealComplete(FlowCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(FlowCtx%cellDMScalSec,FlowCtx%cellDMScalScatter, &
           SCATTER_FORWARD,FlowCtx%permField,ierr);CHKERRQ(ierr)

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowUpdatePerm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Storativity per cell — store into FlowCtx%storativity Vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowUpdateStorativity"

   Subroutine MEF90FlowUpdateStorativity(FlowCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellOpt
      Type(IS)                                            :: CellSetGlobalIS, setIS
      PetscInt,dimension(:),Pointer                       :: setIdx
      PetscInt,dimension(:),Pointer                       :: setID
      PetscInt                                            :: set, cell
      PetscReal                                           :: S_eff
      PetscReal,Dimension(:),Pointer                      :: Sptr

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellOpt,ierr);CHKERRQ(ierr)
         S_eff = gOpt%storativity
         If (cellOpt%storativityOverride >= 0.0_Kr) S_eff = cellOpt%storativityOverride

         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Do cell = 1, size(setIdx)
            Call SectionRealRestrict(FlowCtx%cellDMScalSec,setIdx(cell),Sptr,ierr);CHKERRQ(ierr)
            Sptr(1) = S_eff
            Call SectionRealRestore(FlowCtx%cellDMScalSec,setIdx(cell),Sptr,ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do

      Call SectionRealComplete(FlowCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(FlowCtx%cellDMScalSec,FlowCtx%cellDMScalScatter, &
           SCATTER_FORWARD,FlowCtx%storativity,ierr);CHKERRQ(ierr)

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowUpdateStorativity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  CREATE SNES NONLINEAR SOLVER for pressure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCreateSNES"

   Subroutine MEF90FlowCreateSNES(FlowCtx,SNESPressure,residual,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      Type(SNES),Intent(OUT)                              :: SNESPressure
      Type(Vec),Intent(INOUT)                             :: residual
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(Mat)                                          :: matP, matPpre
      Type(KSP)                                          :: kspP
      Type(MatNullSpace)                                 :: nspP

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)

      !!! Gate 3.1 fix (checklist item C): MEF90FlowBilinearFormCore does
      !!! MatCopy(Jac,PreJac) and zero-rows BOTH matrices. Passing the same
      !!! Mat as A and B is therefore invalid — use two independent objects,
      !!! exactly as MEF90FlowCreateTS now does.
      Call DMCreateMatrix(FlowCtx%DMScal,MATAIJ,matP,ierr);CHKERRQ(ierr)
      Call DMCreateMatrix(FlowCtx%DMScal,MATAIJ,matPpre,ierr);CHKERRQ(ierr)
      Call MatSetOptionsPrefix(matP,'flow_',ierr);CHKERRQ(ierr)
      Call MatSetOptionsPrefix(matPpre,'flow_',ierr);CHKERRQ(ierr)
      Call MatSetOption(matP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matPpre,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (gOpt%addNullSpace) Then
         Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspP,ierr);CHKERRQ(ierr)
         Call MatSetNullSpace(matP,nspP,ierr);CHKERRQ(ierr)
      End If
      Call MatSetFromOptions(matP,ierr);CHKERRQ(ierr)
      Call MatSetFromOptions(matPpre,ierr);CHKERRQ(ierr)

      Call SNESCreate(FlowCtx%MEF90Ctx%comm,SNESPressure,ierr);CHKERRQ(ierr)
      Call SNESSetOptionsPrefix(SNESPressure,'flow_',ierr);CHKERRQ(ierr)
      Call SNESSetDM(SNESPressure,FlowCtx%DMScal,ierr);CHKERRQ(ierr)
      Call SNESSetType(SNESPressure,SNESKSPONLY,ierr);CHKERRQ(ierr)
      Call SNESSetFunction(SNESPressure,residual,MEF90FlowOperator,FlowCtx,ierr);CHKERRQ(ierr)
      Call SNESSetJacobian(SNESPressure,matP,matPpre,MEF90FlowBilinearForm,FlowCtx,ierr);CHKERRQ(ierr)
      Call SNESSetTolerances(SNESPressure,gOpt%pressureATol,gOpt%pressureRTol, &
           PETSC_DEFAULT_DOUBLE_PRECISION,gOpt%pressureMaxIt,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(SNESPressure,"SNESPressure",ierr);CHKERRQ(ierr)
      Call SNESSetFromOptions(SNESPressure,ierr);CHKERRQ(ierr)
      Call SNESGetKSP(SNESPressure,kspP,ierr);CHKERRQ(ierr)
      Call KSPSetType(kspP,KSPCG,ierr);CHKERRQ(ierr)
      Call KSPSetInitialGuessNonzero(kspP,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (gOpt%addNullSpace) Then
         Call KSPSetNullSpace(kspP,nspP,ierr);CHKERRQ(ierr)
      End If
      Call KSPSetTolerances(kspP,1.0D-8,PETSC_DEFAULT_DOUBLE_PRECISION,1.0D+10,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      Call KSPSetFromOptions(kspP,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowCreateSNES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  CREATE TS TRANSIENT SOLVER for pressure (theta-method/BDF)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCreateTS"

   Subroutine MEF90FlowCreateTS(FlowCtx,tsPressure,residual,timeStart,dt,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      Type(TS),Intent(OUT)                                :: tsPressure
      Type(Vec),Intent(INOUT)                             :: residual
      PetscReal,Intent(IN)                                :: timeStart, dt
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(Mat)                                          :: matP, matPpre
      Type(MatNullSpace)                                 :: nspP
      Type(SNES)                                         :: snesP
      Type(KSP)                                          :: kspP

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)

      !!! Gate 3.1: IJacobian needs REAL matrices. MEF90FlowBilinearFormCore does
      !!! MatCopy(Jac,PreJac,...) and MEF90FlowApplyDirichletJacobian zero-rows
      !!! both of them, so Jac and PreJac MUST be two distinct Mat objects.
      Call DMCreateMatrix(FlowCtx%DMScal,MATAIJ,matP,ierr);CHKERRQ(ierr)
      Call DMCreateMatrix(FlowCtx%DMScal,MATAIJ,matPpre,ierr);CHKERRQ(ierr)
      Call MatSetOptionsPrefix(matP,'flow_',ierr);CHKERRQ(ierr)
      Call MatSetOptionsPrefix(matPpre,'flow_',ierr);CHKERRQ(ierr)
      Call MatSetOption(matP,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matPpre,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matP,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matPpre,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matPpre,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (gOpt%addNullSpace) Then
         Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspP,ierr);CHKERRQ(ierr)
         Call MatSetNullSpace(matP,nspP,ierr);CHKERRQ(ierr)
      End If
      Call MatSetFromOptions(matP,ierr);CHKERRQ(ierr)
      Call MatSetFromOptions(matPpre,ierr);CHKERRQ(ierr)

      Call TSCreate(FlowCtx%MEF90Ctx%comm,tsPressure,ierr);CHKERRQ(ierr)
      Call TSSetDM(tsPressure,FlowCtx%DMScal,ierr);CHKERRQ(ierr)
      Call TSSetOptionsPrefix(tsPressure,'flow_',ierr);CHKERRQ(ierr)
      !!! Gate 3.1 fix: pass the registered type NAME, exactly as the already
      !!! compiling HeatXfer template does with 'rosw'. The Fortran constant
      !!! TSBEULER has no precedent in the compiled sources. Also follow the
      !!! HeatXfer order: TSSetType BEFORE TSSetProblemType.
      Call TSSetType(tsPressure,'beuler',ierr);CHKERRQ(ierr)
      Call TSSetProblemType(tsPressure,TS_LINEAR,ierr);CHKERRQ(ierr)
      Call TSSetIFunction(tsPressure,residual,MEF90FlowIFunction,FlowCtx,ierr);CHKERRQ(ierr)
      Call TSSetIJacobian(tsPressure,matP,matPpre,MEF90FlowIJacobian,FlowCtx,ierr);CHKERRQ(ierr)
      Call TSSetSolution(tsPressure,FlowCtx%pressure,ierr);CHKERRQ(ierr)
      Call TSSetInitialTimeStep(tsPressure,timeStart,dt,ierr);CHKERRQ(ierr)
      Call TSSetExactFinalTime(tsPressure,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(tsPressure,"TSPressure",ierr);CHKERRQ(ierr)
      Call TSSetFromOptions(tsPressure,ierr);CHKERRQ(ierr)

      Call TSGetSNES(tsPressure,snesP,ierr);CHKERRQ(ierr)
      Call SNESGetKSP(snesP,kspP,ierr);CHKERRQ(ierr)
      Call KSPSetType(kspP,KSPCG,ierr);CHKERRQ(ierr)
      Call KSPSetInitialGuessNonzero(kspP,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (gOpt%addNullSpace) Then
         Call KSPSetNullSpace(kspP,nspP,ierr);CHKERRQ(ierr)
      End If
      Call KSPSetTolerances(kspP,1.0D-8,PETSC_DEFAULT_DOUBLE_PRECISION,1.0D+10, &
           PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      Call KSPSetFromOptions(kspP,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowCreateTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  COMPUTE CELL-CENTERED ∇p — volume-weighted average per cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowComputePressureGradient"

   Subroutine MEF90FlowComputePressureGradient(pressure,cellPressureGradient,FlowCtx,ierr)
      Type(Vec),Intent(IN)                               :: pressure
      Type(Vec),Intent(INOUT)                            :: cellPressureGradient
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowComputePressureGradient2D(pressure,cellPressureGradient,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowComputePressureGradient3D(pressure,cellPressureGradient,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowComputePressureGradient


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  APPLY BIOT BODY FORCE  f_p = -α · ∇p  →  DefMechCtx%force
!!!  Internal dim dispatch — ComputePressureGradient depends on dim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowApplyBiotForce"

   Subroutine MEF90FlowApplyBiotForce(FlowCtx,DefMechCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      Type(MEF90DefMechCtx_Type),Intent(INOUT)           :: DefMechCtx
      PetscErrorCode,Intent(OUT)                          :: ierr
      PetscInt                                            :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowApplyBiotForce2D(FlowCtx,DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowApplyBiotForce3D(FlowCtx,DefMechCtx,ierr)
      End If
   End Subroutine MEF90FlowApplyBiotForce


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  COMPUTE CELL-CENTERED μ(T) — temperature-dependent viscosity
!!!  μ(T) formula options (same pattern as widely used in porous-media literature):
!!!    Sutherland (liquid):  μ(T) = μ_ref * (T/T_ref)^1.5 * (T_ref+C) / (T+C)
!!!    Arrhenius-like (gas): μ(T) = μ_ref * (T/T_ref)^0.5 * exp( Ea/R * (1/T - 1/T_ref) )
!!!  Uses Gauss-point averaged temperature from DefMechCtx%temperature (shared with HeatXfer)
!!!  Cell scalar written to FlowCtx%cellMuField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowUpdateMuFromTemperature"

   Subroutine MEF90FlowUpdateMuFromTemperature(FlowCtx,DefMechCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: DefMechCtx
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellOpt
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer       :: elem
      Type(MEF90Element_Type)                             :: elemType
      Type(IS)                                            :: CellSetGlobalIS, setIS
      PetscInt,dimension(:),Pointer                       :: setIdx
      PetscInt,dimension(:),Pointer                       :: setID
      PetscInt                                            :: set, cell, iDoF, iGauss, QuadratureOrder
      PetscReal                                           :: TCell, TGauss, TRef, TMin, TMax
      PetscReal                                           :: mu_ref, C, Ea_over_R, T_dimless
      PetscReal,Dimension(:),Pointer                      :: TPtr
      PetscReal,Dimension(:),Pointer                      :: muPtr

      !!! Must have temperature-enabled (HeatXfer active + shared pointer)
      If (.NOT. Associated(DefMechCtx%temperature)) Return

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      TRef   = max(gOpt%referenceTemperature,1.0_Kr)   !!! avoid T_ref=0
      mu_ref = gOpt%dynamicViscosity
      C      = gOpt%sutherlandConstant
      Ea_over_R = gOpt%arrheniusActivationEnergy
      TMin = TRef * 0.5_Kr   !!! clamp to prevent μ→∞ at T→0
      TMax = TRef * 3.0_Kr   !!! upper clamp for numerical safety

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellOpt,ierr);CHKERRQ(ierr)

         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)

         elemType = MEF90knownElements(cellOpt%ElemTypeShortID)
         QuadratureOrder = elemType%order * 2
         Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
              cellOpt%ElemTypeShortID,ierr);CHKERRQ(ierr)

         Allocate(TPtr(elemType%numDof))
         Do cell = 1, size(setIdx)
            !!! Restrict temperature field to this cell's DOFs
            Call SectionRealRestrictClosure(DefMechCtx%DMScalSec,FlowCtx%DM,setIdx(cell),elemType%numDof,TPtr,ierr);CHKERRQ(ierr)

            !!! Gauss-point average → cell mean temperature
            TCell = 0.0_Kr
            Do iGauss = 1, size(elem(cell)%Gauss_C)
               TGauss = 0.0_Kr
               Do iDoF = 1, elemType%numDof
                  TGauss = TGauss + TPtr(iDoF) * elem(cell)%BF(iDoF,iGauss)
               End Do
               TCell = TCell + elem(cell)%Gauss_C(iGauss) * TGauss
            End Do
            TCell = TCell / Sum(elem(cell)%Gauss_C)

            !!! Clamp
            If (TCell < TMin) TCell = TMin
            If (TCell > TMax) TCell = TMax

            !!! Compute μ(T)
            T_dimless = TCell / TRef
            If (gOpt%useArrheniusViscosity) Then
               !!! Arrhenius-like (gas): μ ∝ T^0.5 * exp(Ea/(R*T))
               muPtr(1) = mu_ref * sqrt(T_dimless) * exp(Ea_over_R * (1.0_Kr/TCell - 1.0_Kr/TRef))
            Else
               !!! Sutherland (liquid): μ ∝ T^1.5 * (T_ref+C)/(T+C)
               muPtr(1) = mu_ref * T_dimless * sqrt(T_dimless) * (TRef + C) / (TCell + C)
            End If

            !!! Write to cellMuField Section
            Call SectionRealRestrict(FlowCtx%cellDMScalSec,setIdx(cell),muPtr,ierr);CHKERRQ(ierr)
            Call SectionRealRestore(FlowCtx%cellDMScalSec,setIdx(cell),muPtr,ierr);CHKERRQ(ierr)
         End Do
         DeAllocate(TPtr)
         Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr)
      End Do

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      !!! Scatter Section → global Vec
      Call SectionRealComplete(FlowCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(FlowCtx%cellDMScalSec,FlowCtx%cellDMScalScatter, &
           SCATTER_FORWARD,FlowCtx%cellMuField,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowUpdateMuFromTemperature


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  COMPUTE CELL-CENTERED DARCY VELOCITY
!!!  v_darcy = -(k(x)/μ)·∇p + bodyForce_contribution
!!!  where:  k(x) comes from permField(cell) = k/μ, so v = -permField·∇p  (simpler!)
!!!  Uses cellPressureGradient (already computed) and permField (k/μ already there)
!!!  Result written to FlowCtx%cellDarcyVelocity (cellDMVect, dim components per cell)
!!!  Also adds gravity body force from advectionVector per cell set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowComputeDarcyVelocity"

   Subroutine MEF90FlowComputeDarcyVelocity(FlowCtx,DefMechCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: DefMechCtx
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellOpt
      Type(SectionReal)                                    :: permFieldSec, gradpSec, vSec
      Type(IS)                                            :: CellSetGlobalIS, setIS
      PetscInt,dimension(:),Pointer                       :: setIdx
      PetscInt,dimension(:),Pointer                       :: setID
      PetscInt                                            :: set, cell, dim
      PetscReal,Dimension(:),Pointer                      :: permCell
      PetscReal,Dimension(:),Pointer                      :: gradpPtr, vPtr
      PetscReal,Dimension(3)                              :: gravityContrib
      PetscBool                                           :: hasGravity

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)

      !!! Zero out velocity field
      Call VecSet(FlowCtx%cellDarcyVelocity,0.0_Kr,ierr)

      !!! Sections for read/write
      Call SectionRealDuplicate(FlowCtx%cellDMScalSec,permFieldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMVectSec,gradpSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(FlowCtx%cellDMVectSec,vSec,ierr);CHKERRQ(ierr)

      !!! Scatter global Vecs → local Sections
      Call SectionRealToVec(permFieldSec,FlowCtx%cellDMScalScatter,SCATTER_REVERSE,FlowCtx%permField,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(gradpSec,FlowCtx%cellDMVectScatter,SCATTER_REVERSE,FlowCtx%cellPressureGradient,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellOpt,ierr);CHKERRQ(ierr)

         !!! Check if this cell set has gravity/body force
         hasGravity = .FALSE.
         If (any(abs(cellOpt%advectionVector(1:dim)) > 0.0_Kr)) hasGravity = .TRUE.

         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)

         Do cell = 1, size(setIdx)
            !!! Read cell permeability (already = k/μ)
            Call SectionRealRestrict(permFieldSec,setIdx(cell),permCell,ierr);CHKERRQ(ierr)

            !!! Read cell-centered ∇p
            Call SectionRealRestrict(gradpSec,setIdx(cell),gradpPtr,ierr);CHKERRQ(ierr)

            !!! v_darcy = -permField·∇p   (permField already is k/μ)
            !!! Then add gravity contribution if enabled:  v_grav = permField·(ρ_f · g_vec)  [signed differently]
            Call SectionRealRestrict(vSec,setIdx(cell),vPtr,ierr);CHKERRQ(ierr)
            vPtr(1:dim) = -permCell(1) * gradpPtr(1:dim)

            If (hasGravity) Then
               gravityContrib(1:dim) = cellOpt%advectionVector(1:dim)
               vPtr(1:dim) = vPtr(1:dim) + permCell(1) * gravityContrib(1:dim) * gOpt%fluidDensity
            End If

            Call SectionRealRestore(vSec,setIdx(cell),vPtr,ierr);CHKERRQ(ierr)
         End Do

         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr)
      End Do

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr)

      !!! Scatter local Section → global Vec
      Call SectionRealComplete(vSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(vSec,FlowCtx%cellDMVectScatter,SCATTER_FORWARD,FlowCtx%cellDarcyVelocity,ierr);CHKERRQ(ierr)

      !!! Clean up Sections
      Call SectionRealDestroy(permFieldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(gradpSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(vSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowComputeDarcyVelocity


End Module m_MEF90_Flow
