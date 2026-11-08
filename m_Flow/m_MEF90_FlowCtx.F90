#include "../MEF90/mef90.inc"
Module m_MEF90_FlowCtx_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   Private
   Public :: MEF90FlowCtx_Type
   Public :: MEF90FlowGlobalOptions_Type
   Public :: MEF90FlowCellSetOptions_Type
   Public :: MEF90FlowVertexSetOptions_Type

   !!!
   !!!  Context type for Darcy flow in porous media
   !!!  Pressure-only formulation using standard FEM (same mathematical structure as HeatXfer)
   !!!  Equation:  -div( k(p,v)/mu  grad(p) ) = S dp/dt + f_source
   !!!  where:  k(v) = k_min + (k_max - k_min) * damage   (phase-field coupling)
   !!!
   Type MEF90FlowCtx_Type
      !!! vertex-based vecs
      Type(Vec),pointer                :: pressure
      Type(Vec),pointer                :: boundaryPressure
      Type(Vec),pointer                :: pressurePreviousStep

      !!! cell-based vecs
      Type(Vec),pointer                :: source              ! well/volume source term f_source
      Type(Vec),pointer                :: permField           ! k(v)/mu  assembled diffusivity per cell
      Type(Vec),pointer                :: storativity         ! storativity S per cell
      Type(Vec),pointer                :: cellPressureGradient! cell-centered ∇p for Biot coupling
      Type(Vec),pointer                :: baseForce           ! snapshot of DefMechCtx%force (gravity etc.) before Biot superposition
      Type(Vec),pointer                :: cellDarcyVelocity   ! v_darcy = -(k(x)/μ)·∇p + (k(x)/μ)·ρg  (cell-centered vector, for heat advection)
      Type(Vec),pointer                :: cellMuField         ! μ(T)  assembled viscosity per cell (for thermal→flow coupling)

      !!! PetscBag options
      PetscBag                         :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer    :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: VertexSetOptionsBag

      !!! Material properties shared with DefMech
      PetscBag,Dimension(:),Pointer    :: MaterialPropertiesBag

      !!! Mesh & Section infrastructure
      Type(MEF90Ctx_Type),pointer      :: MEF90Ctx
      Type(DM),pointer                 :: DM
      Type(DM)                         :: DMScal
      Type(DM)                         :: cellDMScal
      Type(DM)                         :: cellDMVect           ! clone for cell-centered vectors (∇p, force)
      Type(VecScatter)                 :: DMScalScatter
      Type(VecScatter)                 :: cellDMScalScatter
      Type(VecScatter)                 :: cellDMVectScatter
      Type(SectionReal)                :: DMScalSec
      Type(SectionReal)                :: cellDMScalSec
      Type(SectionReal)                :: cellDMVectSec
      !!! Gate 4.2: persistent cell-centered v_darcy section, kept alive
      !!! between MEF90FlowComputeDarcyVelocity calls so HeatXferAssembly can
      !!! read it (via MEF90HeatXferCtx%flowVelocitySec) during SNES/TS
      !!! residual + Jacobian assembly. SectionRealDuplicate of
      !!! cellDMVectSec creates an independent section with the same layout.
      !!! TARGET lives on the *variable* in vDef
      !!!   Type(MEF90FlowCtx_Type),Target :: MEF90FlowCtx
      !!! not on this component: gfortran rejects TARGET inside a TYPE
      !!! definition ("Attribute ... is not allowed in a TYPE definition",
      !!! CI #88). Fortran 95 rules make every non-pointer subobject of a
      !!! TARGET variable a valid pointer-assignment target, so
      !!!   MEF90HeatXferCtx%flowVelocitySec => MEF90FlowCtx%cellDarcyVelocitySec
      !!! remains valid and Associated(...) guards still work.
      Type(SectionReal)                :: cellDarcyVelocitySec

      !!! Solvers
      Type(SNES)                       :: SNESPressure
      Type(KSP)                        :: KSPPressure
      Type(TS)                         :: TSPressure

      !!! Energy viewer for diagnostics
      Type(PetscViewer)                :: globalEnergyViewer
   End Type MEF90FlowCtx_Type

   !!!
   !!!  Global options — registered via PetscBag, set from YAML
   !!!
   Type MEF90FlowGlobalOptions_Type
      PetscEnum                        :: timeSteppingType     ! NULL, SteadyState, Transient
      PetscBool                        :: addNullSpace

      !!! Hydraulic parameters (global defaults)
      PetscReal                        :: permeabilityMin      ! k_min  (m^2)
      PetscReal                        :: permeabilityMax      ! k_max  (m^2)
      PetscReal                        :: dynamicViscosity     ! mu(T0) reference viscosity at T_ref  (Pa.s)
      PetscReal                        :: storativity          ! S      (1/Pa)
      PetscReal                        :: initialPressure      ! p0     (Pa)
      PetscReal                        :: residualDamage       ! damage cutoff below which k = k_min only
      PetscReal                        :: biotCoefficient      ! Biot α — poroelastic coupling (0 = disabled, 1 = full)

      !!! Thermal→Flow coupling: μ(T) Sutherland-like / Arrhenius temperature-dependent viscosity
      PetscBool                        :: temperatureDependViscosity ! enable μ(T) coupling (requires HeatXfer active)
      PetscReal                        :: referenceTemperature       ! T_ref (K) where dynamicViscosity is defined
      PetscReal                        :: sutherlandConstant         ! Sutherland constant C (K) for liquid μ(T); air C~120K
      PetscReal                        :: arrheniusActivationEnergy  ! Ea/R (K) for gas-like μ(T) ∝ T^0.5 exp(Ea/(R*T))
      PetscBool                        :: useArrheniusViscosity      ! .true. use Arrhenius, .false. use Sutherland

      !!! Darcy velocity → HeatXfer advection coupling
      PetscBool                        :: enableFlowAdvection       ! add v_darcy·∇T convection term to HeatXfer (requires HeatXfer)
      PetscReal                        :: fluidDensity              ! ρ_f (kg/m^3) needed for buoyancy/gravity body force

      !!! I/O offsets & scalings
      PetscInt                         :: pressureOffset
      PetscInt                         :: sourceOffset
      PetscInt                         :: boundaryPressureOffset
      PetscInt                         :: pressurePreviousStepOffset
      PetscInt                         :: permFieldOffset
      PetscInt                         :: darcyVelocityOffset
      PetscInt                         :: sourceScaling
      PetscInt                         :: boundaryPressureScaling
      PetscInt                         :: pressurePreviousStepScaling
      PetscInt                         :: pressureScaling

      !!! Solver parameters
      PetscReal                        :: pressureATol
      PetscReal                        :: pressureRTol
      PetscInt                         :: pressureMaxIt
   End Type MEF90FlowGlobalOptions_Type

   !!!
   !!!  Per-cell-set options (material-specific overrides)
   !!!
   Type MEF90FlowCellSetOptions_Type
      PetscInt                         :: elemTypeShortID
      PetscReal                        :: source                  ! volumetric source/sink (1/s)
      PetscReal                        :: permeabilityOverride    ! if >0, overrides k_min for this cell set
      PetscReal                        :: storativityOverride     ! if >=0, overrides global S
      PetscBool                        :: Has_BC                  ! cell-face pressure BC
      PetscReal                        :: boundaryPressure        ! pressure value (Pa)
      PetscReal,Dimension(3)           :: advectionVector         ! gravity body force (m/s^2)
   End Type MEF90FlowCellSetOptions_Type

   !!!
   !!!  Per-vertex-set options
   !!!
   Type MEF90FlowVertexSetOptions_Type
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryPressure
   End Type MEF90FlowVertexSetOptions_Type

End Module m_MEF90_FlowCtx_Type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private wrappers to bind PetscBag to the custom Option types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module m_MEF90FlowGlobalOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_FlowCtx_Type
   Implicit None
   Private
   Public :: PetscBagGetDataMEF90FlowCtxGlobalOptions

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_FlowCtx_Type
         PetscBag                                           :: bag
         Type(MEF90FlowGlobalOptions_Type),pointer          :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataMEF90FlowCtxGlobalOptions(bag,data,ierr)
      PetscBag                                           :: bag
      Type(MEF90FlowGlobalOptions_Type),pointer          :: data
      PetscErrorCode                                    :: ierr
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90FlowCtxGlobalOptions
End Module m_MEF90FlowGlobalOptions_Private

Module m_MEF90FlowCellSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_FlowCtx_Type
   Implicit None
   Private
   Public :: PetscBagGetDataMEF90FlowCtxCellSetOptions

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_FlowCtx_Type
         PetscBag                                           :: bag
         Type(MEF90FlowCellSetOptions_Type),pointer         :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataMEF90FlowCtxCellSetOptions(bag,data,ierr)
      PetscBag                                           :: bag
      Type(MEF90FlowCellSetOptions_Type),pointer         :: data
      PetscErrorCode                                    :: ierr
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90FlowCtxCellSetOptions
End Module m_MEF90FlowCellSetOptions_Private

Module m_MEF90FlowVertexSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_FlowCtx_Type
   Implicit None
   Private
   Public :: PetscBagGetDataMEF90FlowCtxVertexSetOptions

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_FlowCtx_Type
         PetscBag                                           :: bag
         Type(MEF90FlowVertexSetOptions_Type),pointer       :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataMEF90FlowCtxVertexSetOptions(bag,data,ierr)
      PetscBag                                           :: bag
      Type(MEF90FlowVertexSetOptions_Type),pointer       :: data
      PetscErrorCode                                    :: ierr
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90FlowCtxVertexSetOptions
End Module m_MEF90FlowVertexSetOptions_Private

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Main context module — Enum lists + Create / SetSections / CreateVectors / Destroy
!!!                        + SetFromOptions (PetscBag registration)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module m_MEF90_FlowCtx
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_FlowCtx_Type
   Use m_MEF90FlowGlobalOptions_Private
   Use m_MEF90FlowCellSetOptions_Private
   Use m_MEF90FlowVertexSetOptions_Private
   Implicit none

   PetscSizeT,protected   :: sizeofMEF90FlowGlobalOptions
   PetscSizeT,protected   :: sizeofMEF90FlowCellSetOptions
   PetscSizeT,protected   :: sizeofMEF90FlowVertexSetOptions

   Enum,bind(c)
      enumerator  :: MEF90Flow_timeSteppingTypeNULL = 0,      &
                     MEF90Flow_timeSteppingTypeSteadyState,   &
                     MEF90Flow_timeSteppingTypeTransient
   End Enum
   Character(len = MEF90_MXSTRLEN),dimension(6),protected :: MEF90Flow_timeSteppingTypeList

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  PRIVATE INIT — compute sizes, populate enum list strings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCtxInitialize_Private"

   Subroutine MEF90FlowCtxInitialize_Private(ierr)
      PetscErrorCode,Intent(OUT)                         :: ierr
      Type(MEF90FlowGlobalOptions_Type)                  :: FlowGlobalOptions
      Type(MEF90FlowCellSetOptions_Type)                 :: FlowCellSetOptions
      Type(MEF90FlowVertexSetOptions_Type)               :: FlowVertexSetOptions
      character(len=1),pointer                           :: dummychar(:)
      PetscSizeT                                         :: sizeofchar

      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofMEF90FlowGlobalOptions  = size(transfer(FlowGlobalOptions,dummychar))  * sizeofchar
      sizeofMEF90FlowCellSetOptions  = size(transfer(FlowCellSetOptions,dummychar)) * sizeofchar
      sizeofMEF90FlowVertexSetOptions = size(transfer(FlowVertexSetOptions,dummychar)) * sizeofchar

      MEF90Flow_timeSteppingTypeList(1) = 'Null'
      MEF90Flow_timeSteppingTypeList(2) = 'SteadyState'
      MEF90Flow_timeSteppingTypeList(3) = 'Transient'
      MEF90Flow_timeSteppingTypeList(4) = 'MEF90Flow_timeSteppingType'
      MEF90Flow_timeSteppingTypeList(5) = '_MEF90Flow_timeSteppingType'
      MEF90Flow_timeSteppingTypeList(6) = ''
   End Subroutine MEF90FlowCtxInitialize_Private


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  CONTEXT CREATE — instantiate DM clones, bags, sections, scatter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCtxCreate"

   Subroutine MEF90FlowCtxCreate(FlowCtx,Mesh,MEF90Ctx,ierr)
      Type(MEF90FlowCtx_Type),Intent(OUT)              :: FlowCtx
      Type(DM),target,Intent(IN)                       :: Mesh
      Type(MEF90Ctx_Type),target,Intent(IN)            :: MEF90Ctx
      PetscErrorCode,Intent(OUT)                       :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer        :: FlowGlobalOptions
      Type(MEF90FlowCellSetOptions_Type),pointer       :: FlowCellSetOptions
      Type(MEF90FlowVertexSetOptions_Type),pointer     :: FlowVertexSetOptions
      Type(IS)                                         :: setIS
      PetscInt                                         :: numSet, set, dim
      Character(len=MEF90_MXSTRLEN)                    :: filename, IOBuffer

      Call MEF90FlowCtxInitialize_Private(ierr)
      Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)

      FlowCtx%DM      => Mesh
      FlowCtx%MEF90Ctx => MEF90Ctx

      !!! Clone DM for vertex scalar DOF (pressure) and cell scalar DOF (source/perm)
      Call DMMeshClone(Mesh,FlowCtx%DMScal,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(FlowCtx%DMScal,1,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(FlowCtx%DMScal,1,ierr);CHKERRQ(ierr)

      Call DMMeshClone(Mesh,FlowCtx%cellDMScal,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(FlowCtx%cellDMScal,1,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(FlowCtx%cellDMScal,1,ierr);CHKERRQ(ierr)

      Call DMMeshClone(Mesh,FlowCtx%cellDMVect,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(FlowCtx%cellDMVect,dim,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(FlowCtx%cellDMVect,dim,ierr);CHKERRQ(ierr)

      !!! Global options bag
      Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90FlowGlobalOptions,FlowCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)

      !!! Cell set option bags — per material/zone
      Call DMmeshGetLabelIdIS(Mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr)
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(FlowCtx%CellSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90FlowCellSetOptions,FlowCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      !!! Vertex set option bags — per BC surface
      Call DMmeshGetLabelIdIS(Mesh,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr)
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(FlowCtx%VertexSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90FlowVertexSetOptions,FlowCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      !!! Energy viewer (optional diagnostics)
      filename = trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.flowener'
      Call PetscViewerASCIIOpen(MEF90Ctx%comm,filename,FlowCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(FlowCtx%globalEnergyViewer, &
         "# step  time  pressure  injectionWork  storage  total\n",ierr);CHKERRQ(ierr)
      Call PetscViewerFlush(FlowCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)

      !!! Nullify Vec pointers until CreateVectors is called
      Nullify(FlowCtx%pressure)
      Nullify(FlowCtx%boundaryPressure)
      Nullify(FlowCtx%pressurePreviousStep)
      Nullify(FlowCtx%source)
      Nullify(FlowCtx%permField)
      Nullify(FlowCtx%storativity)
      Nullify(FlowCtx%cellPressureGradient)
      Nullify(FlowCtx%baseForce)
      Nullify(FlowCtx%cellDarcyVelocity)
      Nullify(FlowCtx%cellMuField)
   End Subroutine MEF90FlowCtxCreate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SET SECTIONS — link DM → Sections → Scatter for vertex & cell scalar data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCtxSetSections"

   Subroutine MEF90FlowCtxSetSections(FlowCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)            :: FlowCtx
      PetscErrorCode,Intent(OUT)                        :: ierr
      PetscInt                                          :: dim

      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      Call DMMeshGetVertexSectionReal(FlowCtx%DMScal,"default",1,FlowCtx%DMScalSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(FlowCtx%DMScal,"default",FlowCtx%DMScalSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(FlowCtx%DMScal,1,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(FlowCtx%DMScal,FlowCtx%DMScalSec,FlowCtx%DMScalScatter,ierr);CHKERRQ(ierr)

      Call DMMeshGetCellSectionReal(FlowCtx%cellDMScal,"default",1,FlowCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(FlowCtx%cellDMScal,"default",FlowCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(FlowCtx%cellDMScal,1,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(FlowCtx%cellDMScal,FlowCtx%cellDMScalSec,FlowCtx%cellDMScalScatter,ierr);CHKERRQ(ierr)

      Call DMMeshGetCellSectionReal(FlowCtx%cellDMVect,"default",dim,FlowCtx%cellDMVectSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(FlowCtx%cellDMVect,"default",FlowCtx%cellDMVectSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(FlowCtx%cellDMVect,dim,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(FlowCtx%cellDMVect,FlowCtx%cellDMVectSec,FlowCtx%cellDMVectScatter,ierr);CHKERRQ(ierr)
      !!! Gate 4.2: independent section that will hold v_darcy per cell.
      Call SectionRealDuplicate(FlowCtx%cellDMVectSec,FlowCtx%cellDarcyVelocitySec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowCtxSetSections


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  CREATE VECTORS — allocate all global Vecs and name them for output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCtxCreateVectors"

   Subroutine MEF90FlowCtxCreateVectors(FlowCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)            :: FlowCtx
      PetscErrorCode,Intent(OUT)                        :: ierr
      Type(MEF90FlowGlobalOptions_Type),pointer         :: gOpt

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)

      Allocate(FlowCtx%pressure,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%DMScal,FlowCtx%pressure,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%pressure,"Pressure",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%pressure,gOpt%initialPressure,ierr);CHKERRQ(ierr)

      Allocate(FlowCtx%pressurePreviousStep,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%DMScal,FlowCtx%pressurePreviousStep,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%pressurePreviousStep,"PressurePreviousStep",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%pressurePreviousStep,gOpt%initialPressure,ierr);CHKERRQ(ierr)

      Allocate(FlowCtx%boundaryPressure,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%DMScal,FlowCtx%boundaryPressure,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%boundaryPressure,"BoundaryPressure",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%boundaryPressure,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(FlowCtx%source,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%cellDMScal,FlowCtx%source,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%source,"Source",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%source,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(FlowCtx%permField,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%cellDMScal,FlowCtx%permField,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%permField,"PermField",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%permField,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(FlowCtx%storativity,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%cellDMScal,FlowCtx%storativity,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%storativity,"Storativity",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%storativity,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(FlowCtx%cellPressureGradient,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%cellDMVect,FlowCtx%cellPressureGradient,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%cellPressureGradient,"PressureGradient",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%cellPressureGradient,0.0_Kr,ierr);CHKERRQ(ierr)

      !!! baseForce — snapshot of user-preset body forces (gravity etc.)
      !!! DefMechCtx%force is on cellDMVect (dim-components per cell)
      Allocate(FlowCtx%baseForce,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%cellDMVect,FlowCtx%baseForce,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%baseForce,"BaseForce",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%baseForce,0.0_Kr,ierr);CHKERRQ(ierr)

      !!! cellDarcyVelocity — v_darcy = -(k(x)/μ)·∇p + body_force_term
      !!! cell-centered vector (dim-components), used by HeatXfer advection coupling
      Allocate(FlowCtx%cellDarcyVelocity,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%cellDMVect,FlowCtx%cellDarcyVelocity,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%cellDarcyVelocity,"DarcyVelocity",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%cellDarcyVelocity,0.0_Kr,ierr);CHKERRQ(ierr)

      !!! cellMuField — μ(T) assembled viscosity per cell (temperature-dependent flow coupling)
      !!! same DM as permField (cellDMScal, scalar per cell)
      Allocate(FlowCtx%cellMuField,stat=ierr)
      Call DMCreateGlobalVector(FlowCtx%cellDMScal,FlowCtx%cellMuField,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(FlowCtx%cellMuField,"ViscosityField",ierr);CHKERRQ(ierr)
      Call VecSet(FlowCtx%cellMuField,0.0_Kr,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowCtxCreateVectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  DESTROY VECTORS — reverse of CreateVectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCtxDestroyVectors"

   Subroutine MEF90FlowCtxDestroyVectors(FlowCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)            :: FlowCtx
      PetscErrorCode,Intent(OUT)                        :: ierr

      If (Associated(FlowCtx%pressure)) Then
         Call VecDestroy(FlowCtx%pressure,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%pressure); Nullify(FlowCtx%pressure)
      End If
      If (Associated(FlowCtx%pressurePreviousStep)) Then
         Call VecDestroy(FlowCtx%pressurePreviousStep,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%pressurePreviousStep); Nullify(FlowCtx%pressurePreviousStep)
      End If
      If (Associated(FlowCtx%boundaryPressure)) Then
         Call VecDestroy(FlowCtx%boundaryPressure,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%boundaryPressure); Nullify(FlowCtx%boundaryPressure)
      End If
      If (Associated(FlowCtx%source)) Then
         Call VecDestroy(FlowCtx%source,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%source); Nullify(FlowCtx%source)
      End If
      If (Associated(FlowCtx%permField)) Then
         Call VecDestroy(FlowCtx%permField,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%permField); Nullify(FlowCtx%permField)
      End If
      If (Associated(FlowCtx%storativity)) Then
         Call VecDestroy(FlowCtx%storativity,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%storativity); Nullify(FlowCtx%storativity)
      End If
      If (Associated(FlowCtx%cellPressureGradient)) Then
         Call VecDestroy(FlowCtx%cellPressureGradient,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%cellPressureGradient); Nullify(FlowCtx%cellPressureGradient)
      End If
      If (Associated(FlowCtx%baseForce)) Then
         Call VecDestroy(FlowCtx%baseForce,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%baseForce); Nullify(FlowCtx%baseForce)
      End If
      If (Associated(FlowCtx%cellDarcyVelocity)) Then
         Call VecDestroy(FlowCtx%cellDarcyVelocity,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%cellDarcyVelocity); Nullify(FlowCtx%cellDarcyVelocity)
      End If
      If (Associated(FlowCtx%cellMuField)) Then
         Call VecDestroy(FlowCtx%cellMuField,ierr);CHKERRQ(ierr)
         DeAllocate(FlowCtx%cellMuField); Nullify(FlowCtx%cellMuField)
      End If
   End Subroutine MEF90FlowCtxDestroyVectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  CONTEXT DESTROY — release bags, sections, scatter, DM clones, viewers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCtxDestroy"

   Subroutine MEF90FlowCtxDestroy(FlowCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(OUT)              :: FlowCtx
      PetscErrorCode,Intent(OUT)                       :: ierr
      PetscInt                                         :: set

      Call PetscBagDestroy(FlowCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
      Do set = 1, size(FlowCtx%CellSetOptionsBag)
         Call PetscBagDestroy(FlowCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(FlowCtx%CellSetOptionsBag)
      Do set = 1, size(FlowCtx%VertexSetOptionsBag)
         Call PetscBagDestroy(FlowCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(FlowCtx%VertexSetOptionsBag)

      Call PetscViewerDestroy(FlowCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)

      Nullify(FlowCtx%pressure)
      Nullify(FlowCtx%boundaryPressure)
      Nullify(FlowCtx%pressurePreviousStep)
      Nullify(FlowCtx%source)
      Nullify(FlowCtx%permField)
      Nullify(FlowCtx%storativity)
      Nullify(FlowCtx%cellPressureGradient)

      Call VecScatterDestroy(FlowCtx%DMScalScatter,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(FlowCtx%cellDMScalScatter,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(FlowCtx%cellDMVectScatter,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(FlowCtx%DMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(FlowCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(FlowCtx%cellDMVectSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(FlowCtx%cellDarcyVelocitySec,ierr);CHKERRQ(ierr)
      Call DMDestroy(FlowCtx%DMScal,ierr);CHKERRQ(ierr)
      Call DMDestroy(FlowCtx%cellDMScal,ierr);CHKERRQ(ierr)
      Call DMDestroy(FlowCtx%cellDMVect,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowCtxDestroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  PETSCBAG REGISTRATION — bind Fortran option struct fields to Petsc CLI options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90FlowCtxGlobalOptions"

   Subroutine PetscBagRegisterMEF90FlowCtxGlobalOptions(bag,name,prefix,default,ierr)
      PetscBag                                             :: bag
      Character(len=*),Intent(IN)                          :: prefix,name
      Type(MEF90FlowGlobalOptions_Type),Intent(IN)         :: default
      PetscErrorCode,Intent(OUT)                           :: ierr
      Type(MEF90FlowGlobalOptions_Type),pointer            :: FlowGlobalOptions

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(bag,FlowGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"Flow Global Options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,FlowGlobalOptions%timeSteppingType, &
           MEF90Flow_timeSteppingTypeList,default%timeSteppingType, &
           'timestepping_type','Darcy flow time-stepping type',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,FlowGlobalOptions%addNullSpace,default%addNullSpace, &
           'addNullSpace','Add null-space to flow SNES',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterReal(bag,FlowGlobalOptions%permeabilityMin,default%permeabilityMin, &
           'perm_min','Residual permeability (m^2)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%permeabilityMax,default%permeabilityMax, &
           'perm_max','Fully-broken permeability (m^2)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%dynamicViscosity,default%dynamicViscosity, &
           'viscosity','Dynamic viscosity (Pa.s)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%storativity,default%storativity, &
           'storativity','Storativity S (1/Pa)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%initialPressure,default%initialPressure, &
           'p0','Initial pore pressure (Pa)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%residualDamage,default%residualDamage, &
           'residual_damage','Damage threshold below which k = k_min only',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%biotCoefficient,default%biotCoefficient, &
           'biot_coeff','Biot poroelastic coefficient α (0 = no coupling)',ierr);CHKERRQ(ierr)

      !!! Thermal→Flow coupling: μ(T)
      !!! CRITICAL (CI #50): registration order MUST match the field order of
      !!! MEF90FlowGlobalOptions_Type byte-for-byte. PetscBag writes each value
      !!! at bag->data + bag->offset, and PetscBagGetData casts bag->data back
      !!! to the Fortran derived type. A Bool/Real swap or a missing field
      !!! shifts every later field (EXO offsets and solver tolerances) by
      !!! several bytes, so gOpt reads garbage. The old code swapped
      !!! arrheniusActivationEnergy/useArrheniusViscosity and omitted 3 fields.
      Call PetscBagRegisterBool (bag,FlowGlobalOptions%temperatureDependViscosity, &
           default%temperatureDependViscosity,'temp_dep_viscosity', &
           'Enable temperature-dependent viscosity μ(T)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%referenceTemperature, &
           default%referenceTemperature,'T_ref', &
           'Reference temperature for μ(T) = μ_ref at T_ref (K)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%sutherlandConstant, &
           default%sutherlandConstant,'sutherland_C', &
           'Sutherland constant C (K); air~120, water~600',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%arrheniusActivationEnergy, &
           default%arrheniusActivationEnergy,'Ea_over_R', &
           'Arrhenius Ea/R (K); air~110, water~5100',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool (bag,FlowGlobalOptions%useArrheniusViscosity, &
           default%useArrheniusViscosity,'use_arrhenius', &
           'Use Arrhenius instead of Sutherland for μ(T)',ierr);CHKERRQ(ierr)

      !!! Flow→Thermal advection coupling
      Call PetscBagRegisterBool (bag,FlowGlobalOptions%enableFlowAdvection, &
           default%enableFlowAdvection,'enable_advection', &
           'Add v_darcy·∇T convection term to HeatXfer',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%fluidDensity, &
           default%fluidDensity,'fluid_rho', &
           'Fluid density ρ_f (kg/m^3) for buoyancy term',ierr);CHKERRQ(ierr)

      !!! I/O offsets & scalings — all 10 PetscInt fields, in declaration order
      Call PetscBagRegisterInt (bag,FlowGlobalOptions%pressureOffset,default%pressureOffset, &
           'pressure_offset','EXO file offset for pressure',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,FlowGlobalOptions%sourceOffset,default%sourceOffset, &
           'source_offset','EXO file offset for source term',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,FlowGlobalOptions%boundaryPressureOffset,default%boundaryPressureOffset, &
           'boundary_pressure_offset','EXO offset for boundary pressure',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,FlowGlobalOptions%pressurePreviousStepOffset, &
           default%pressurePreviousStepOffset,'pressure_previous_step_offset', &
           'EXO offset for pressure at previous step',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,FlowGlobalOptions%permFieldOffset,default%permFieldOffset, &
           'perm_field_offset','EXO offset for permeability field k(v)/mu',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,FlowGlobalOptions%darcyVelocityOffset,default%darcyVelocityOffset, &
           'darcy_velocity_offset','EXO offset for Darcy velocity vector',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,FlowGlobalOptions%sourceScaling,MEF90ScalingList, &
           default%sourceScaling,'source_scaling','Source term scaling law',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,FlowGlobalOptions%boundaryPressureScaling,MEF90ScalingList, &
           default%boundaryPressureScaling,'bc_pressure_scaling','Boundary pressure scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,FlowGlobalOptions%pressurePreviousStepScaling,MEF90ScalingList, &
           default%pressurePreviousStepScaling,'pressure_previous_step_scaling', &
           'Scaling law for pressure at previous step',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,FlowGlobalOptions%pressureScaling,MEF90ScalingList, &
           default%pressureScaling,'pressure_scaling','Scaling law for pressure',ierr);CHKERRQ(ierr)

      !!! Solver params
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%pressureATol,default%pressureATol, &
           'pressure_atol','Pressure absolute tolerance',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowGlobalOptions%pressureRTol,default%pressureRTol, &
           'pressure_rtol','Pressure relative tolerance',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,FlowGlobalOptions%pressureMaxIt,default%pressureMaxIt, &
           'pressure_maxit','Pressure solver max iterations',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90FlowCtxGlobalOptions


#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90FlowCtxCellSetOptions"

   Subroutine PetscBagRegisterMEF90FlowCtxCellSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                              :: bag
      Character(len=*),Intent(IN)                           :: prefix,name
      Type(MEF90FlowCellSetOptions_Type),Intent(IN)         :: default
      PetscErrorCode,Intent(OUT)                            :: ierr
      Type(MEF90FlowCellSetOptions_Type),pointer            :: FlowCellSetOptions

      Call PetscBagGetDataMEF90FlowCtxCellSetOptions(bag,FlowCellSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"Flow CellSet Options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterInt (bag,FlowCellSetOptions%elemTypeShortID,default%elemTypeShortID, &
           'elem_short_id','Short ID for flow element type',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowCellSetOptions%source,default%source, &
           'source','Volumetric source/sink rate (1/s)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowCellSetOptions%permeabilityOverride, &
           default%permeabilityOverride, &
           'perm_override','Override k_min (m^2), 0 to use global',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowCellSetOptions%storativityOverride, &
           default%storativityOverride, &
           'storativity_override','Override S (1/Pa), -1 to use global',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,FlowCellSetOptions%Has_BC,default%Has_BC, &
           'flowBC','Cell set has pressure boundary condition',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowCellSetOptions%boundaryPressure,default%boundaryPressure, &
           'bc_pressure','Pressure value on BC surface (Pa)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowCellSetOptions%advectionVector(1),default%advectionVector(1), &
           'gravity_x','x-component of gravity/body force (m/s^2)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowCellSetOptions%advectionVector(2),default%advectionVector(2), &
           'gravity_y','y-component of gravity (m/s^2)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowCellSetOptions%advectionVector(3),default%advectionVector(3), &
           'gravity_z','z-component of gravity (m/s^2)',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90FlowCtxCellSetOptions


#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90FlowCtxVertexSetOptions"

   Subroutine PetscBagRegisterMEF90FlowCtxVertexSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                              :: bag
      Character(len=*),Intent(IN)                           :: prefix,name
      Type(MEF90FlowVertexSetOptions_Type),Intent(IN)       :: default
      PetscErrorCode,Intent(OUT)                            :: ierr
      Type(MEF90FlowVertexSetOptions_Type),pointer          :: FlowVertexSetOptions

      Call PetscBagGetDataMEF90FlowCtxVertexSetOptions(bag,FlowVertexSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"Flow VertexSet Options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterBool(bag,FlowVertexSetOptions%Has_BC,default%Has_BC, &
           'flow_has_bc','Vertex set has pressure Dirichlet BC',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,FlowVertexSetOptions%boundaryPressure, &
           default%boundaryPressure,'Pressure value (Pa)','Dirichlet BC pressure (Pa)',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90FlowCtxVertexSetOptions


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SET FROM OPTIONS — populate all bags with YAML/default values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCtxSetFromOptions"

   Subroutine MEF90FlowCtxSetFromOptions(FlowCtx,prefix,defaultGlobalOptions, &
                                          defaultCellSetOptions,defaultVertexSetOptions,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)                 :: FlowCtx
      Character(len=*),Intent(IN)                           :: prefix
      Type(MEF90FlowGlobalOptions_Type),Intent(IN)          :: defaultGlobalOptions
      Type(MEF90FlowCellSetOptions_Type),Intent(IN)         :: defaultCellSetOptions
      Type(MEF90FlowVertexSetOptions_Type),Intent(IN)       :: defaultVertexSetOptions
      PetscErrorCode,Intent(OUT)                            :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer             :: gOpt
      Type(MEF90CtxGlobalOptions_Type),pointer              :: MEF90CtxGlobalOptions
      Type(IS)                                              :: setIS
      PetscInt                                               :: set, numSet
      PetscInt,dimension(:),Pointer                          :: setID
      Character(len=MEF90_MXSTRLEN)                          :: setName, setprefix
      Type(MEF90FlowCellSetOptions_Type)                     :: myDefaultCellSetOptions
      Type(MEF90FlowVertexSetOptions_Type)                   :: myDefaultVertexSetOptions
      Type(MEF90Element_Type),Dimension(:),Pointer          :: ElemTypeScal

      Call PetscBagGetDataMEF90CtxGlobalOptions(FlowCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr);CHKERRQ(ierr)

      !!! Global bag — YAML top-level key is 'flow:' → option names are
      !!! '-flow_timestepping_type', '-flow_perm_min', ... exactly like
      !!! DefMech ('DefMech_*') and HeatXfer ('heatxfer_*'). The empty prefix
      !!! used before CI #51 made EVERY Flow global option fall back to its
      !!! default (SteadyState, p0=0, k_min=1e-20), so the TS path was disabled.
      Call PetscBagRegisterMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag, &
           'flow','flow_',defaultGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)

      !!! CI #93 / Gate 4.3 visibility: DefMech, HeatXfer, Materials and
      !!! MEF90Ctx all dump their bags under verbose>0, but Flow never did.
      !!! That is why bc_pressure_scaling / source_scaling / biot_coeff /
      !!! enable_advection / boundary_pressure_offset were absent from every
      !!! CI log and Gates 4.1/4.2 had to invent their own trace lines.
      !!! Pure print-out: NO physics, NO option values are changed.
      If (MEF90CtxGlobalOptions%verbose > 0) Then
         Call PetscBagView(FlowCtx%GlobalOptionsBag,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      End If

      !!! Cell set bags — per-set prefix + EXO elem type auto-detection (DefMech pattern)
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)

      !!! Override default elemTypeShortID from EXO file (DefMech line 1037 pattern)
      Call EXOGetCellSetElementType_Scal(FlowCtx%MEF90Ctx,ElemTypeScal,ierr)

      Do set = 1, numSet
         Write(setName,'(A,I0)') 'cs', setID(set)
         Write(setprefix,'(A,I4.4,A)') 'cs', setID(set), '_'
         myDefaultCellSetOptions = defaultCellSetOptions
         myDefaultCellSetOptions%elemTypeShortID = ElemTypeScal(set)%ShortID
         Call PetscBagRegisterMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set), &
              setName,setprefix,myDefaultCellSetOptions,ierr);CHKERRQ(ierr)
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Call PetscBagView(FlowCtx%CellSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      DeAllocate(ElemTypeScal)

      !!! Vertex set bags — per-set prefix (vsXXXX_ pattern)
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Do set = 1, numSet
         Write(setName,'(A,I0)') 'vs', setID(set)
         Write(setprefix,'(A,I4.4,A)') 'vs', setID(set), '_'
         myDefaultVertexSetOptions = defaultVertexSetOptions
         Call PetscBagRegisterMEF90FlowCtxVertexSetOptions(FlowCtx%VertexSetOptionsBag(set), &
              setName,setprefix,myDefaultVertexSetOptions,ierr);CHKERRQ(ierr)
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Call PetscBagView(FlowCtx%VertexSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowCtxSetFromOptions

End Module m_MEF90_FlowCtx
