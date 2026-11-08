#include "../MEF90/mef90.inc"
Program vDef
#include <finclude/petscdef.h>
   Use petsc
   Use m_MEF90
   !Use m_vDef
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_HeatXfer
   Use m_MEF90_FlowCtx
   Use m_MEF90_Flow
   use m_vDefDefault
   Implicit NONE   

   PetscErrorCode                                     :: ierr
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions

   !!! Defect mechanics contexts
   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions

   !!! HeatXfer contexts
   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Pointer      :: MEF90HeatXferGlobalOptions

   
   Type(DM),target                                    :: Mesh
   Type(IS)                                           :: setIS,CellSetGlobalIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   PetscReal,Dimension(:),Pointer                     :: time
   PetscReal,Dimension(:),Pointer                     :: thermalEnergySet,heatFluxWorkSet
   PetscReal,Dimension(:),Pointer                     :: elasticEnergySet,surfaceEnergySet,forceWorkSet,cohesiveEnergySet
   PetscReal,Dimension(:),Pointer                     :: elasticEnergy,surfaceEnergy,forceWork,cohesiveEnergy,totalMechanicalEnergy

   !!!! NEW: Flow (Darcy pressure) context for Phase 1 coupled HM simulation
   Type(MEF90FlowCtx_Type)                           :: MEF90FlowCtx
   Type(MEF90FlowGlobalOptions_Type),Pointer          :: MEF90FlowGlobalOptions
   Type(SNES)                                         :: snesPressure
   Type(Vec)                                          :: residualPressure
   !!! Gate 3.1 — transient Darcy pressure solver (only created when
   !!! MEF90FlowGlobalOptions%timeSteppingType == Transient)
   Type(TS)                                           :: tsPressure
   TSConvergedReason                                  :: tsPressureConvergedReason
   Logical                                            :: flowTransient
   PetscReal                                          :: tFlowEnd, dtFlow0
   !!! Gate 3.2: global pressure extrema per outer step. VecMax is the
   !!! undrained bottom pressure -> Terzaghi decay probe (zero local compile,
   !!! VecMin/VecMax already used in m_MEF90_DefMechAssembly.F90:91-92).
   PetscReal                                          :: pressureMin, pressureMax
   !!! Gate 4.2: per-step infinity norm of the cell-centered Darcy velocity
   PetscReal                                          :: darcyVmin, darcyVmax

   SNESConvergedReason                                :: snesDispConvergedReason
   Type(Vec)                                          :: residualDisp
   Type(Vec)                                          :: damageAltMinOld,displacementAltMinOld
   Type(Vec)                                          :: damageLB,damageUB   
   PetscReal,Dimension(:),Pointer                     :: damageArray,damageAltMinOldArray,damageLBArray,damageUBArray
   PetscInt                                           :: iDof
   PetscReal                                          :: SOROmega,mySOROmega
   SNESConvergedReason                                :: snesDamageConvergedReason
   Type(Vec)                                          :: residualDamage,localVec
   PetscInt                                           :: AltMinIter
   PetscReal                                          :: damageMaxChange

   Type(SNES)                                         :: snesTemp
   SNESConvergedReason                                :: snesTempConvergedReason
   Type(TS)                                           :: tsTemp
   TSConvergedReason                                  :: tsTempConvergedReason
   Type(TSAdapt)                                      :: tsAdaptTemp
   Type(Vec)                                          :: residualTemp

   PetscReal                                          :: tsTempInitialStep,tsTempInitialTime
   PetscInt                                           :: tsTempmaxIter
   PetscReal                                          :: t
          
   PetscBool                                          :: flg
   Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
   Type(PetscViewer)                                  :: logViewer
   Integer                                            :: numfield
   
   Integer                                            :: step
   PetscInt                                           :: dim
   PetscInt                                           :: AltMinStep = 0
   
   PetscReal                                          :: alphaMaxChange,alphaMin,alphaMax
   PetscLogStage                                      :: logStageHeat,logStageDamage,logStageDisplacement,logStageEnergy,logStageIO
      


   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call PetscPrintf(PETSC_COMM_WORLD," # vDef: numerical implementation of variational models of Defect Mechanics\n",ierr);CHKERRQ(ierr)
   Call PetscLogStageRegister('HeatXfer    ',logStageHeat,ierr);CHKERRQ(ierr)
   Call PetscLogStageRegister('Damage      ',logStageDamage,ierr);CHKERRQ(ierr)
   Call PetscLogStageRegister('Displacement',logStageDisplacement,ierr);CHKERRQ(ierr)
   Call PetscLogStageRegister('Energy      ',logStageEnergy,ierr);CHKERRQ(ierr)
   Call PetscLogStageRegister('IO          ',logStageIO,ierr);CHKERRQ(ierr)

   !!! Get all MEF90-wide options
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,vDefDefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90CtxGetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(Mesh,dim,ierr);CHKERRQ(ierr) 
   Call DMSetBlockSize(Mesh,dim,ierr);CHKERRQ(ierr)
   
   !!! Open output file
   Call MEF90CtxOpenEXO(MEF90Ctx,Mesh,ierr)
   
   !!! Create DefMech context, get all DefMech options
   Call MEF90DefMechCtxCreate(MEF90DefMechCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   If (dim == 2) Then
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,vDefDefMechDefaultGlobalOptions2D, &
                                         vDefDefMechDefaultCellSetOptions,vDefDefMechDefaultVertexSetOptions,ierr)
   Else
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,vDefDefMechDefaultGlobalOptions3D, &
                                         vDefDefMechDefaultCellSetOptions,vDefDefMechDefaultVertexSetOptions,ierr)
   End If
   Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
   
   !!! Create HeatXfer context, get all HeatXfer options
   Call MEF90HeatXferCtxCreate(MEF90HeatXferCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,vDefHeatXferDefaultGlobalOptions, &
                                       vDefHeatXferDefaultCellSetOptions,vDefHeatXferDefaultVertexSetOptions,ierr)
   Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get material properties bags
   If (dim == 2) Then
      Call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag,MEF90DefMechCtx%DMVect,MEF90Mathium2D,MEF90Ctx,ierr)
   Else
      Call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag,MEF90DefMechCtx%DMVect,MEF90Mathium3D,MEF90Ctx,ierr)
   End If   
   MEF90HeatXferCtx%MaterialPropertiesBag => MEF90DefMechCtx%MaterialPropertiesBag

   !!! NEW: Create Flow (Darcy) context after DefMechCtx
   !!! IMPORTANT: MEF90FlowUpdatePerm reads DefMechCtx%DMScalSec (the damage
   !!! SectionReal), which is not created until MEF90DefMechCtxSetSections
   !!! below. Calling it here dereferences an uninitialised section handle and
   !!! segfaults on every rank (CI #50). Both coefficient updates are therefore
   !!! deferred until right after the DefMech vectors exist.
   Call PetscPrintf(MEF90Ctx%Comm,"Gate 3.1: init 01 FlowCtxCreate\n",ierr);CHKERRQ(ierr)
   Call MEF90FlowCtxCreate(MEF90FlowCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call PetscPrintf(MEF90Ctx%Comm,"Gate 3.1: init 02 FlowCtxSetFromOptions\n",ierr);CHKERRQ(ierr)
   Call MEF90FlowCtxSetFromOptions(MEF90FlowCtx,PETSC_NULL_CHARACTER,vDefFlowDefaultGlobalOptions, &
                                     vDefFlowDefaultCellSetOptions,vDefFlowDefaultVertexSetOptions,ierr)
   Call PetscPrintf(MEF90Ctx%Comm,"Gate 3.1: init 03 FlowBagGetData\n",ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90FlowCtxGlobalOptions(MEF90FlowCtx%GlobalOptionsBag,MEF90FlowGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscPrintf(MEF90Ctx%Comm,"Gate 3.1: init 04 FlowSetSections\n",ierr);CHKERRQ(ierr)
   Call MEF90FlowCtxSetSections(MEF90FlowCtx,ierr);CHKERRQ(ierr)
   !!! Gate 4.2: wire the persistent cell-wise v_darcy section into the
   !!! HeatXfer context. MUST come after FlowCtxSetSections (which creates
   !!! cellDarcyVelocitySec) and after the Flow global bag is read. When
   !!! enable_advection is false HeatXferAssembly keeps using the legacy
   !!! constant cellSetOptions%advectionVector path, so this is a no-op.
   MEF90HeatXferCtx%flowVelocitySec => MEF90FlowCtx%cellDarcyVelocitySec
   MEF90HeatXferCtx%enableFlowAdvection = MEF90FlowGlobalOptions%enableFlowAdvection
   Write(IOBuffer,*) 'Gate 4.2: coupling wired, enable_advection=', &
                     MEF90HeatXferCtx%enableFlowAdvection,'\n'
   Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscPrintf(MEF90Ctx%Comm,"Gate 3.1: init 05 FlowCreateVectors\n",ierr);CHKERRQ(ierr)
   Call MEF90FlowCtxCreateVectors(MEF90FlowCtx,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90FlowCtx%pressure,residualPressure,ierr);CHKERRQ(ierr)
   Call PetscPrintf(MEF90Ctx%Comm,"Gate 3.1: init 06 FlowCreateSNES\n",ierr);CHKERRQ(ierr)
   Call MEF90FlowCreateSNES(MEF90FlowCtx,snesPressure,residualPressure,ierr);CHKERRQ(ierr)

   !!! Create time array from global options
   Call MEF90CtxGetTime(MEF90Ctx,time,ierr)

   !!! Gate 3.1 — build the transient Darcy solver ONCE, after the time array
   !!! exists. It is advanced on the time-step level (before AltMin), not inside
   !!! every alternate-minimization iteration. Steady runs keep using snesPressure.
   flowTransient = (MEF90FlowGlobalOptions%timeSteppingType == MEF90Flow_timeSteppingTypeTransient)
   !!! CI #51: print the enum unconditionally. If the yaml key still fails to
   !!! resolve, the log shows the actual value instead of silently going steady.
   Write(IOBuffer,*) 'Gate 3.1: flow enum',MEF90FlowGlobalOptions%timeSteppingType, &
                     ' Transient',MEF90Flow_timeSteppingTypeTransient,' enabled',flowTransient,'\n'
   Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
   If (flowTransient) Then
      If (size(time) > 1) Then
         dtFlow0 = time(2) - time(1)
      Else
         dtFlow0 = 1.0_Kr
      End If
      Call MEF90FlowCreateTS(MEF90FlowCtx,tsPressure,residualPressure,time(1),dtFlow0,ierr);CHKERRQ(ierr)
      !!! Gate 3.1 self-check: if this line is absent from the CI log, the yaml key
      !!! 'timestepping_type' did not resolve to the Transient enum and the TS path
      !!! is silently disabled.
      Write(IOBuffer,*) 'Gate 3.1: transient Darcy TS created, timeSteppingType=', &
                        MEF90FlowGlobalOptions%timeSteppingType,'\n'
      Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
   End If


   !!! Create sections, vectors, and solvers for DefMech Context
   Call MEF90DefMechCtxSetSections(MEF90DefMechCtx,ierr)
   Call MEF90DefMechCtxCreateVectors(MEF90DefMechCtx,ierr)

   !!! Gate 3.1 / CI #50: the Flow coefficients may only be assembled now that
   !!! DefMechCtx%DMScalSec exists. MEF90FlowUpdatePerm restricts the damage
   !!! field from that section; MEF90FlowUpdateStorativity fills S once before
   !!! the first TSSolve. Both used to be called 30 lines too early and crashed.
   Call PetscPrintf(MEF90Ctx%Comm,"Gate 3.1: init 07 FlowUpdatePerm (after DefMech sections)\n",ierr);CHKERRQ(ierr)
   Call MEF90FlowUpdatePerm(MEF90FlowCtx,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90FlowUpdateStorativity(MEF90FlowCtx,ierr);CHKERRQ(ierr)
   Call PetscPrintf(MEF90Ctx%Comm,"Gate 3.1: init 08 Flow coefficients ready\n",ierr);CHKERRQ(ierr)

   !!! CI #67 — apply the Dirichlet BC to the INITIAL pressure field.
   !!! MEF90FlowCtxCreateVectors sets p = p0 everywhere (no BC), so the
   !!! first TSSolve spent its whole step merely zeroing the drained face:
   !!! for a uniform field K*p = 0, source = 0, and F(U0,0) equals p0-p_bc
   !!! on the Dirichlet dofs and 0 in the interior. Every output record was
   !!! therefore labelled one dt too late: time_whole=t held the BE(n-1)
   !!! field (verified offline: t=0.04 matched BE(1) to 0.05%, t=0.10
   !!! matched BE(4) to 0.34%). Reuse the already-public SNES residual to
   !!! extract exactly that Dirichlet correction and subtract it ONCE before
   !!! the time loop: the interior stays p0, the drained face becomes p_bc.
   !!! From now on the field written at time_whole=t is the BE(n) field with
   !!! n = t/dt, and step 2..51 still take exactly one fixed-dt BE step each.
   Call MEF90FlowUpdateboundaryPressure(MEF90FlowCtx%boundaryPressure,MEF90FlowCtx,ierr);CHKERRQ(ierr)
   Call MEF90FlowOperator(snesPressure,MEF90FlowCtx%pressure,residualPressure,MEF90FlowCtx,ierr);CHKERRQ(ierr)
   Call VecAXPY(MEF90FlowCtx%pressure,-1.0_Kr,residualPressure,ierr);CHKERRQ(ierr)
   Call VecMin(MEF90FlowCtx%pressure,PETSC_NULL_INTEGER,pressureMin,ierr);CHKERRQ(ierr)
   Call VecMax(MEF90FlowCtx%pressure,PETSC_NULL_INTEGER,pressureMax,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) 'CI #67: initial pressure after Dirichlet pmin',pressureMin,' pmax',pressureMax,'\n'
   Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

   Call VecDuplicate(MEF90DefMechCtx%displacement,displacementAltMinOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%displacement,residualDisp,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualDisp,"residualDisp",ierr);CHKERRQ(ierr)
   Call MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx,MEF90DefMechCtx%snesDisp,residualDisp,ierr)
   Call VecDuplicate(MEF90DefMechCtx%damage,damageAltMinOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%damage,residualDamage,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualDamage,"residualDamage",ierr);CHKERRQ(ierr)
   Call MEF90DefMechCreateSNESDamage(MEF90DefMechCtx,MEF90DefMechCtx%snesDamage,residualDamage,ierr)
   DeAllocate(MEF90DefMechCtx%temperature)
      
   !!! Create sections, vectors, and solvers for HeatXfer Context
   If (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90HeatXfer_timeSteppingTypeNULL) Then
      Call MEF90HeatXferCtxSetSections(MEF90HeatXferCtx,ierr)
      Call MEF90HeatXferCtxCreateVectors(MEF90HeatXferCtx,ierr)
      Call VecDuplicate(MEF90HeatXferCtx%temperature,residualTemp,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(residualTemp,"residualTemp",ierr);CHKERRQ(ierr)
      Select Case(MEF90HeatXferGlobalOptions%timeSteppingType)
         Case (MEF90HeatXFer_timeSteppingTypeSteadyState)
            Call MEF90HeatXferCreateSNES(MEF90HeatXferCtx,snesTemp,residualTemp,ierr)
         Case (MEF90HeatXFer_timeSteppingTypeTransient)
            tsTempInitialStep = (time(size(time))-time(1)) / (size(time) - 1.0_Kr) / 10.0_Kr
            tsTempInitialTime = time(1)
            Call MEF90HeatXferCreateTS(MEF90HeatXferCtx,tsTemp,residualTemp,tsTempInitialTime,tsTempInitialStep,ierr)
            Call TSGetAdapt(tsTemp,tsAdaptTemp,ierr);CHKERRQ(ierr)
            Call TSAdaptSetFromOptions(tsAdaptTemp,ierr);CHKERRQ(ierr)
      End Select

      !!! Link the temperature field in the DefMechContext with that of the HeatXfer
      MEF90DefMechCtx%temperature => MEF90HeatXferCtx%temperature
   End If   
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(elasticEnergySet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   elasticEnergySet = 0.0_Kr
   Allocate(surfaceEnergySet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   surfaceEnergySet = 0.0_Kr
   Allocate(forceWorkSet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   forceWorkSet = 0.0_Kr
   Allocate(cohesiveEnergySet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   cohesiveEnergySet = 0.0_Kr
   Allocate(thermalEnergySet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   thermalEnergySet = 0.0_Kr
   Allocate(heatFluxWorkSet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   heatFluxWorkSet = 0.0_Kr
   
   Allocate(elasticEnergy(MEF90GlobalOptions%timeNumStep))
   elasticEnergy = 0.0_Kr
   Allocate(surfaceEnergy(MEF90GlobalOptions%timeNumStep))
   surfaceEnergy = 0.0_Kr
   Allocate(forceWork(MEF90GlobalOptions%timeNumStep))
   forceWork = 0.0_Kr
   Allocate(cohesiveEnergy(MEF90GlobalOptions%timeNumStep))
   cohesiveEnergy = 0.0_Kr
   Allocate(totalMechanicalEnergy(MEF90GlobalOptions%timeNumStep))
   totalMechanicalEnergy = 0.0_Kr
   
   !!!
   !!! Try to figure out if the file was formatted
   !!!
   If (MEF90Ctx%rank == 0) Then
      Call EXGVP(MEF90Ctx%fileExoUnit,"N",numfield,ierr)
   End If
   Call MPI_Bcast(numfield,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)   
   If (numfield == 0) Then
      Call MEF90DefMechFormatEXO(MEF90DefMechCtx,time,ierr)
   End If
   
   !!!
   !!! Actual computations / time stepping
   !!!
   If (.NOT. MEF90GlobalOptions%dryrun) Then
      !!! 
      !!! Reload state if time_skip > 0
      !!!
      If (MEF90GlobalOptions%timeSkip > 0) Then
         If ((MEF90HeatXferGlobalOptions%tempOffset > 0) .AND. (Associated(MEF90HeatXferCtx%temperature))) Then
            Call DMGetLocalVector(MEF90HeatXferCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
            Call VecLoadExodusVertex(MEF90HeatXferCtx%DMScal,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                     MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,MEF90GlobalOptions%timeSkip,MEF90HeatXferGlobalOptions%TempOffset,ierr);CHKERRQ(ierr)
            Call DMLocalToGlobalBegin(MEF90HeatXferCtx%DMScal,localVec,INSERT_VALUES,MEF90HeatXferCtx%Temperature,ierr);CHKERRQ(ierr)
            Call DMLocalToGlobalEnd(MEF90HeatXferCtx%DMScal,localVec,INSERT_VALUES,MEF90HeatXferCtx%Temperature,ierr);CHKERRQ(ierr)
            Call DMRestoreLocalVector(MEF90HeatXferCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
            If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeTransient) Then
               Call TSSetTime(tsTemp,time(MEF90GlobalOptions%timeSkip),ierr);CHKERRQ(ierr)
            End If
         End If
   
         Call DMGetLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
         Call VecLoadExodusVertex(MEF90DefMechCtx%DMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,MEF90GlobalOptions%timeSkip,MEF90DefMechGlobalOptions%DisplacementOffset,ierr);CHKERRQ(ierr)
         Call DMLocalToGlobalBegin(MEF90DefMechCtx%DMVect,localVec,INSERT_VALUES,MEF90DefMechCtx%Displacement,ierr);CHKERRQ(ierr)
         Call DMLocalToGlobalEnd(MEF90DefMechCtx%DMVect,localVec,INSERT_VALUES,MEF90DefMechCtx%Displacement,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
   
         Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
         Call VecLoadExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,MEF90GlobalOptions%timeSkip,MEF90DefMechGlobalOptions%DamageOffset,ierr);CHKERRQ(ierr)
         Call DMLocalToGlobalBegin(MEF90DefMechCtx%DMScal,localVec,INSERT_VALUES,MEF90DefMechCtx%Damage,ierr);CHKERRQ(ierr)
         Call DMLocalToGlobalEnd(MEF90DefMechCtx%DMScal,localVec,INSERT_VALUES,MEF90DefMechCtx%Damage,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      End If

      !!! Start iterations
      step = MEF90GlobalOptions%timeSkip+1

      mainloopQS: Do
         MEF90DefMechCtx%analysisTime = time(step)
         If (step > 1) Then
            MEF90DefMechCtx%timeStep = time(step) - time(step-1)
         Else 
            MEF90DefMechCtx%timeStep = 0.0_Kr
         End If

         Call PetscLogStagePush(logStageHeat,ierr);CHKERRQ(ierr)
         !!! Solve for temperature
         Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
            Case (MEF90HeatXFer_timeSteppingTypeSteadyState) 
               Write(IOBuffer,100) step,time(step)
               Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)

               !!! Update fields
               Call PetscLogStagePush(logStageIO,ierr);CHKERRQ(ierr)
               Call MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr)
               Call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr);
               Call PetscLogStagePop(ierr);CHKERRQ(ierr)
               !!! Solve SNES
               Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,MEF90HeatXferCtx%temperature,ierr);CHKERRQ(ierr)
               Call SNESGetConvergedReason(snesTemp,snesTempConvergedReason,ierr);CHKERRQ(ierr)
               If (snesTempConvergedReason < 0) Then  
                  Write(IOBuffer,400) "temperature",snesTempConvergedReason
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               End If

               !!! Compute thermal energy
               Call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,thermalEnergySet,heatFluxWorkSet,ierr);CHKERRQ(ierr)
               Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
               Call MEF90ISAllGatherMerge(MEF90Ctx%Comm,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
               Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
               Call PetscPrintf(MEF90Ctx%Comm,"\nThermal energies: \n",ierr);CHKERRQ(ierr)
               Do set = 1, size(setID)
                  Write(IOBuffer,101) setID(set),thermalEnergySet(set),heatFluxWorkSet(set),thermalEnergySet(set)-heatFluxWorkSet(set)
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
               Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
               Write(IOBuffer,102) sum(thermalEnergySet),sum(heatFluxWorkSet),sum(thermalEnergySet)-sum(heatFluxWorkSet)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

               !!! Save results
               Call PetscLogStagePush(logStageIO,ierr);CHKERRQ(ierr)
               Call MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
               Call PetscLogStagePop(ierr);CHKERRQ(ierr)
            Case (MEF90HeatXFer_timeSteppingTypeTransient)
               If (step > 1) Then
                  Write(IOBuffer,110) step,time(step)
                  Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
                  !!! Update fields
                  Call PetscLogStagePush(logStageIO,ierr);CHKERRQ(ierr)
                  Call MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr)
                  Call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr);
                  Call PetscLogStagePop(ierr);CHKERRQ(ierr)
                  !!! Make sure TS does not overstep
                  Call TSGetTime(tsTemp,t,ierr);CHKERRQ(ierr)
                  If (t < time(step)) Then
                     Call TSAdaptSetStepLimits(tsAdaptTemp,PETSC_DECIDE,(time(step)-time(step-1))/2.0_Kr,ierr);CHKERRQ(ierr)
                     !!! Something is up here. 
                     !!! replacing the constant 10000 with a variable leads to divergence of TSAdapt
                     !!! when using gcc
                     Call TSSetDuration(tsTemp,10000,time(step),ierr);CHKERRQ(ierr)
                     Call TSSolve(tsTemp,MEF90HeatXferCtx%temperature,time(step),ierr);CHKERRQ(ierr)
                     Call TSGetConvergedReason(tsTemp,tsTempConvergedReason,ierr);CHKERRQ(ierr)
                     If (tsTempConvergedReason < 0) Then  
                        Write(IOBuffer,410) "temperature",tsTempConvergedReason
                        Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                     End If
                     Call TSGetTime(tsTemp,t,ierr);CHKERRQ(ierr)
                     time(step) = t
                  Else
                     Write(IOBuffer,*) 'TS exceeded analysis time. Skipping step\n'
                     Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                  End If
               End If

               !!! Compute thermal energy
               Call PetscLogStagePush(logStageEnergy,ierr);CHKERRQ(ierr)
               Call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,thermalEnergySet,heatFluxWorkSet,ierr);CHKERRQ(ierr)
               Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
               Call MEF90ISAllGatherMerge(MEF90Ctx%Comm,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
               Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
               Call PetscPrintf(MEF90Ctx%Comm,"\nThermal energies: \n",ierr);CHKERRQ(ierr)
               Do set = 1, size(setID)
                  Write(IOBuffer,101) setID(set),thermalEnergySet(set),heatFluxWorkSet(set),thermalEnergySet(set)-heatFluxWorkSet(set)
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
               Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
               Write(IOBuffer,102) sum(thermalEnergySet),sum(heatFluxWorkSet),sum(thermalEnergySet)-sum(heatFluxWorkSet)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               Call PetscLogStagePop(ierr);CHKERRQ(ierr)
               !!! Save results
               Call PetscLogStagePush(logStageIO,ierr);CHKERRQ(ierr)
               Call MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
               Call PetscLogStagePop(ierr);CHKERRQ(ierr)
            Case (MEF90HeatXfer_timeSteppingTypeNULL)
               Continue
            Case default
               Write(IOBuffer,*) "Implemented HeatXfer timeSteppingType: ", MEF90HeatXferGlobalOptions%timeSteppingType, "\n"
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               STOP
         End Select
         Call PetscLogStagePop(ierr);CHKERRQ(ierr)

         Call PetscLogStagePush(logStageDisplacement,ierr);CHKERRQ(ierr)
         !!! Update previous time step fields:
         Call VecCopy(MEF90DefMechCtx%displacement,MEF90DefMechCtx%displacementPreviousStep,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90DefMechCtx%damage,MEF90DefMechCtx%damagePreviousStep,ierr);CHKERRQ(ierr)

         !!! Solve for displacement and damage
         Write(IOBuffer,200) step,time(step) 
         Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
         damageMaxChange = 1.0D+20
         Call MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx,MEF90DefMechCtx%snesDamage,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)

         !!! Update fields
         Call PetscLogStagePush(logStageIO,ierr);CHKERRQ(ierr)
         Call MEF90DefMechSetTransients(MEF90DefMechCtx,step,time(step),ierr)
         Call MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx%displacement,MEF90DefMechCtx,ierr)
         Call MEF90DefMechUpdateboundaryDamage(MEF90DefMechCtx%damage,MEF90DefMechCtx,ierr)
         Call PetscLogStagePop(ierr);CHKERRQ(ierr)

         !!! Gate 3.1 — TRANSIENT DARCY: advance p once per time step, BEFORE
         !!! alternate minimization. p is then frozen during AltMin; only the
         !!! material field k(v) is refreshed there. This is a staggered
         !!! approximation (no re-solve of p for every damage update).
         If (flowTransient) Then
            !!! CI #53: the first transient step SEGVs right after
            !!! "Mechanics: step 1". Print a checkpoint before/after EVERY
            !!! call in this block so the next CI log pinpoints the culprit.
            Write(IOBuffer,*) 'Gate 3.1: step',step,' enter transient block\n'
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call MEF90FlowUpdateMuFromTemperature(MEF90FlowCtx,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.1: step',step,' post-UpdateMu\n'
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call MEF90FlowUpdatePerm(MEF90FlowCtx,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.1: step',step,' post-UpdatePerm\n'
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call MEF90FlowSetTransients(MEF90FlowCtx,step,time(step),ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.1: step',step,' post-SetTransients\n'
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call MEF90FlowUpdateboundaryPressure(MEF90FlowCtx%boundaryPressure,MEF90FlowCtx,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.1: step',step,' post-UpdateBoundaryPressure\n'
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call TSGetTime(tsPressure,t,ierr);CHKERRQ(ierr)
            !!! CI #66: this TS has a FIXED dt from TSSetInitialTimeStep
            !!! (dtFlow0 = time(2)-time(1) = 0.02 s). Driving TSSetDuration
            !!! with the interpolated time(step) produced near-zero intervals
            !!! at steps 12-16 (|time(step)-t| ~ 3e-17): beuler then used
            !!! shift a=1/dt ~ 1e16, J ~ a*S*M, and the pressure froze for
            !!! FIVE outer steps. That is the sole remaining difference
            !!! between CI #65 (ratio 0.1244) and the FE/BE value 0.0934
            !!! (continuum Terzaghi 0.0871). Advance exactly ONE fixed-dt
            !!! backward-Euler step per outer step instead.
            tFlowEnd = t + dtFlow0
            Call TSSetDuration(tsPressure,10000,tFlowEnd,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.1: step',step,' pre-TSSolve TS time',t,' target',tFlowEnd,'\n'
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            !!! CI #53/#66: step 1 is the initial condition (t = time(1) = 0)
            !!! and must not advance the TS; steps 2..51 each take one step.
            If (step > 1) Then
               Call TSSolve(tsPressure,MEF90FlowCtx%pressure,tFlowEnd,ierr);CHKERRQ(ierr)
               Write(IOBuffer,*) 'Gate 3.1: step',step,' post-TSSolve\n'
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               Call TSGetConvergedReason(tsPressure,tsPressureConvergedReason,ierr);CHKERRQ(ierr)
               If (tsPressureConvergedReason < 0) Then
                  Write(IOBuffer,410) "pressure",tsPressureConvergedReason
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               End If
            Else
               Write(IOBuffer,*) 'Gate 3.1: step',step,' TSSolve skipped (no time advance)\n'
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            End If
            !!! Gate 3.1 self-check: report the TS internal time so the CI log
            !!! shows whether p was actually advanced on this step.
            Call TSGetTime(tsPressure,t,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.1: step',step,' TS time',t,' target',tFlowEnd,'\n'
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call VecCopy(MEF90FlowCtx%pressure,MEF90FlowCtx%pressurePreviousStep,ierr);CHKERRQ(ierr)
            !!! Gate 3.2: emit the pressure envelope of THIS outer step.
            !!! VecMin/VecMax are collective but every rank executes this
            !!! block, so no rank-dependent deadlock (CI #57 lesson).
            !!! pmax is the impermeable-bottom pressure: Terzaghi predicts it
            !!! decays from p0=1e5 towards ~8.7e3 at t=1 s (T~1.09).
            Call VecMin(MEF90FlowCtx%pressure,PETSC_NULL_INTEGER,pressureMin,ierr);CHKERRQ(ierr)
            Call VecMax(MEF90FlowCtx%pressure,PETSC_NULL_INTEGER,pressureMax,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) 'Gate 3.2: step',step,' pmin',pressureMin,' pmax',pressureMax,'\n'
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
         End If

         !!! Gate 4.1: snapshot the TIME-SCALED preset force (gravity etc.)
         !!! into baseForce, ONCE PER TIME STEP (not once per AltMin iter).
         !!! MEF90DefMechSetTransients above re-applies ForceScaling on every
         !!! outer step; the default is MEF90Scaling_Linear. NOTE (CI #85):
         !!! terzaghi1D.yaml writes `Force: 0, -9.81, 0` WITH spaces after the
         !!! commas, so PETSc 3.3's RealArray reader keeps only the first
         !!! element and the CellSet bag view reads `Force = 0 0 0`. Hence
         !!! |DefMechCtx%force| is currently 0 (= linear x 0), NOT 9.81*time.
         !!! A non-zero body force is exercised only by the Gate 4.1
         !!! BODYFORCE diagnostic deck. Snapshotting only once at startup
         !!! (as the original checklist 4.1 suggested) would STILL be wrong
         !!! once a non-zero force is enabled: Linear scaling would be pinned
         !!! at its t=0 value.
         !!! MEF90FlowApplyBiotForce rebuilds force = baseForce - alpha*grad(p)
         !!! from this snapshot on every AltMin call, so the operation is
         !!! idempotent (no accumulation) and the scaling is preserved.
         If (Associated(MEF90FlowCtx%baseForce) .AND. Associated(MEF90DefMechCtx%force)) Then
            Call VecCopy(MEF90DefMechCtx%force,MEF90FlowCtx%baseForce,ierr);CHKERRQ(ierr)
         End If

         Select case(MEF90DefMechGlobalOptions%timeSteppingType)
            Case (MEF90DefMech_TimeSteppingTypeQuasiStatic)
               Select case(MEF90DefMechGlobalOptions%SolverType)
                  Case(MEF90DefMech_SolverTypeAltMin)
                     !Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDamage,1,ierr);CHKERRQ(ierr)
                     Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDisp,1,ierr);CHKERRQ(ierr)

                     AltMin: Do AltMinIter = 1, MEF90DefMechGlobalOptions%maxit
                        AltMinStep = altminstep + 1
                        Write(IObuffer,208) AltMinIter
                        Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

                        !!! Solve SNES
                        If (mod(AltMinIter-1,MEF90DefMechGlobalOptions%PCLag) == 0) Then
                           !Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDisp,-2,ierr);CHKERRQ(ierr)
                           Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDamage,-2,ierr);CHKERRQ(ierr)
                        End If 

                        Call PetscLogStagePush(logStageDisplacement,ierr);CHKERRQ(ierr)
                        Call VecCopy(MEF90DefMechCtx%displacement,displacementAltMinOld,ierr);CHKERRQ(ierr)

                        !!! NEW — Phase 1 StepμT→k(v)→p: 热→渗流→力学 完整耦合链
                        !!! 完整顺序：T → μ(T) → k(v)/μ(T) → p → ∇p → Biot f_p → u + v
                        !!!
                        !!! Stepμ: temperature-dependent viscosity μ(T)
                        Call MEF90FlowUpdateMuFromTemperature(MEF90FlowCtx,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
                        !!! StepP: update k(v) (now uses cellMuField if temp_dep_viscosity enabled)
                        Call MEF90FlowUpdatePerm(MEF90FlowCtx,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
                        !!! Gate 3.1 — in transient mode p was already advanced on the
                        !!! time-step level; do NOT re-solve it here, and do not re-apply
                        !!! time scaling / boundary pressure every AltMin iteration.
                        If (.not. flowTransient) Then
                           Call MEF90FlowSetTransients(MEF90FlowCtx,step,MEF90DefMechCtx%analysisTime,ierr)
                           Call MEF90FlowUpdateboundaryPressure(MEF90FlowCtx%boundaryPressure,MEF90FlowCtx,ierr)
                           !!! StepP_solve: steady Darcy pressure
                           Call SNESSolve(snesPressure,PETSC_NULL_OBJECT,MEF90FlowCtx%pressure,ierr);CHKERRQ(ierr)
                        End If

                        !!! StepBiot: f_p = -α·∇p → DefMechCtx%force (内部先算 ∇p)
                        Call MEF90FlowApplyBiotForce(MEF90FlowCtx,MEF90DefMechCtx,ierr);CHKERRQ(ierr)

                        !!! Stepv_darcy: 为热对流项 v·∇T 预计算 Darcy 速度场
                        !!! 启用 enableFlowAdvection 时 HeatXferAssembly 将读取 cellDarcyVelocity
                        Call MEF90FlowComputeDarcyVelocity(MEF90FlowCtx,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
                        !!! Gate 4.2 diagnostic: one row per outer step. The
                        !!! velocity field is computed in BOTH variants; only
                        !!! the enable flag differs, so |v|inf must track the
                        !!! pressure gradient in both logs.
                        Call VecMin(MEF90FlowCtx%cellDarcyVelocity,PETSC_NULL_INTEGER,darcyVmin,ierr);CHKERRQ(ierr)
                        Call VecMax(MEF90FlowCtx%cellDarcyVelocity,PETSC_NULL_INTEGER,darcyVmax,ierr);CHKERRQ(ierr)
                        Write(IOBuffer,*) 'Gate 4.2: enabled=',MEF90HeatXferCtx%enableFlowAdvection, &
                             ' |v_darcy|inf=',max(abs(darcyVmin),abs(darcyVmax)), &
                             ' vmin=',darcyVmin,' vmax=',darcyVmax,'\n'
                        Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

                        Call SNESSolve(MEF90DefMechCtx%snesDisp,PETSC_NULL_OBJECT,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)
                        Call SNESGetConvergedReason(MEF90DefMechCtx%snesDisp,snesDispConvergedReason,ierr);CHKERRQ(ierr)
                        If (snesDispConvergedReason < 0) Then  
                           Write(IOBuffer,400) "displacement",snesDispConvergedReason
                           Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                        End If

                        !!! Over relaxation on u, skipping the first alternate minimization step
                        If ((AltMinIter > 1) .AND. (SOROmega /= 1.0_Kr)) Then
                           SOROmega = abs(MEF90DefMechGlobalOptions%SOROmega)
                           Call VecAXPBY(MEF90DefMechCtx%displacement,1.0_Kr - SOROmega,SOROmega,displacementAltMinOld,ierr);CHKERRQ(ierr)
                        EndIf

                        Call PetscLogStagePop(ierr);CHKERRQ(ierr)

                        Call PetscLogStagePush(logStageDamage,ierr);CHKERRQ(ierr)
                        Call VecCopy(MEF90DefMechCtx%damage,damageAltMinOld,ierr);CHKERRQ(ierr)
                        Call SNESSolve(MEF90DefMechCtx%snesDamage,PETSC_NULL_OBJECT,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
                        Call SNESGetConvergedReason(MEF90DefMechCtx%snesDamage,snesDamageConvergedReason,ierr);CHKERRQ(ierr)
                        If (snesDamageConvergedReason < 0) Then
                           Write(IOBuffer,400) "damage field",snesDamageConvergedReason
                           Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                        End If
                        !!! Over relaxation of the damage variable
                        If (AltMinIter > 1) Then
                           If ((MEF90DefMechGlobalOptions%SOROmega > 0.0_Kr) .AND. (MEF90DefMechGlobalOptions%SOROmega /= 1.0)) Then
                              mySOROmega = MEF90DefMechGlobalOptions%SOROmega
                              !!! LIMITED SOR
                              Call SNESVIGetVariableBounds(MEF90DefMechCtx%snesDamage,damageLB,damageUB,ierr);CHKERRQ(ierr)
                              Call VecGetArrayF90(damageLB,damageLBArray,ierr);CHKERRQ(ierr)
                              Call VecGetArrayF90(damageUB,damageUBArray,ierr);CHKERRQ(ierr)
                              Call VecGetArrayF90(damageAltMinOld,damageAltMinOldArray,ierr);CHKERRQ(ierr)
                              Call VecGetArrayF90(MEF90DefMechCtx%damage,damageArray,ierr);CHKERRQ(ierr)
                              Do iDof = 1, size(damageArray)
                                 If (damageArray(iDof) > damageAltMinOldArray(iDof)) Then
                                    mySOROmega = min(mySOROmega,(damageUBArray(iDof)-damageAltMinOldArray(iDof)) / (damageArray(iDof) - damageAltMinOldArray(iDof)))
                                 Else If (damageArray(iDof) < damageAltMinOldArray(iDof)) Then
                                    mySOROmega = min(mySOROmega,(damageLBArray(iDof)-damageAltMinOldArray(iDof)) / (damageArray(iDof) - damageAltMinOldArray(iDof)))
                                 End If
                              End Do
                              Call MPI_AllReduce(mySOROmega,SOROmega,1,MPIU_SCALAR,MPI_MIN,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
                              Call VecRestoreArrayF90(MEF90DefMechCtx%damage,damageArray,ierr);CHKERRQ(ierr)
                              Call VecRestoreArrayF90(damageAltMinOld,damageAltMinOldArray,ierr);CHKERRQ(ierr)
                              Call VecRestoreArrayF90(damageUB,damageUBArray,ierr);CHKERRQ(ierr)
                              Call VecRestoreArrayF90(damageLB,damageLBArray,ierr);CHKERRQ(ierr)
                              Call VecAXPBY(MEF90DefMechCtx%damage,1.0_Kr - SOROmega,SOROmega,damageAltMinOld,ierr);CHKERRQ(ierr)
                           Else If (MEF90DefMechGlobalOptions%SOROmega < 0.0_Kr) Then
                              !!! PROJECTED SOR
                              SOROmega = -MEF90DefMechGlobalOptions%SOROmega
                              Call VecAXPBY(MEF90DefMechCtx%damage,1.0_Kr - SOROmega,SOROmega,damageAltMinOld,ierr);CHKERRQ(ierr)
                              Call SNESVIGetVariableBounds(MEF90DefMechCtx%snesDamage,damageLB,damageUB,ierr);CHKERRQ(ierr)
                              Call VecPointwiseMax(MEF90DefMechCtx%damage,MEF90DefMechCtx%damage,damageLB,ierr);CHKERRQ(ierr)
                              Call VecPointwiseMin(MEF90DefMechCtx%damage,MEF90DefMechCtx%damage,damageUB,ierr);CHKERRQ(ierr)
                           EndIf
                        End If

                        Call VecMin(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,alphaMin,ierr);CHKERRQ(ierr)
                        Call VecMax(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,alphaMax,ierr);CHKERRQ(ierr)
                        Call VecAxPy(damageAltMinOld,-1.0_Kr,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
                        ! damageAltMinOld = damageAltMinOld - damage
                        Call VecNorm(damageAltMinOld,NORM_INFINITY,damageMaxChange,ierr);CHKERRQ(ierr)
                        Write(IOBuffer,209) alphamin,alphamax,damageMaxChange
                        Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                        Call PetscLogStagePop(ierr);CHKERRQ(ierr)

                        If (damageMaxChange <= MEF90DefMechGlobalOptions%damageATol) Then
                           EXIT
                        End If

                        If (mod(AltMinIter,25) == 0) Then
                           Call PetscLogStagePush(logStageIO,ierr);CHKERRQ(ierr)
                           Call MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
                           Call MEF90FlowViewEXO(MEF90FlowCtx,step,ierr)
                           Call PetscLogStagePop(ierr);CHKERRQ(ierr)
                        End If
                     End Do AltMin
                  Case(MEF90DefMech_SolverTypeQuasiNewton1,MEF90DefMech_SolverTypeQuasiNewton2)
                     Call PetscLogStagePush(logStageDisplacement,ierr);CHKERRQ(ierr)
                     Call SNESSolve(MEF90DefMechCtx%snesDisp,PETSC_NULL_OBJECT,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)
                     Call PetscLogStagePop(ierr);CHKERRQ(ierr)
                  Case default
                     Write(IOBuffer,*) "Unimplemented DefMech solver type: ", MEF90DefMechGlobalOptions%SolverType, "\n"
                     Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                     STOP
               End Select
            Case (MEF90DefMech_TimeSteppingTypeNULL)
               Continue
            Case default
               Write(IOBuffer,*) "Unimplemented DefMech time stepping: ", MEF90DefMechGlobalOptions%timeSteppingType, "\n"
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               STOP
         End Select
         Call PetscLogStagePop(ierr);CHKERRQ(ierr)

         !!! Compute energies
         elasticEnergySet  = 0.0_Kr
         forceWorkSet      = 0.0_Kr
         surfaceEnergySet  = 0.0_Kr
         cohesiveEnergySet = 0.0_Kr
         Call PetscLogStagePush(logStageEnergy,ierr);CHKERRQ(ierr)
         Call MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,elasticEnergySet,ierr);CHKERRQ(ierr)
         Call MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,forceWorkSet,ierr);CHKERRQ(ierr)
         Call MEF90DefMechSurfaceEnergy(MEF90DefMechCtx%damage,MEF90DefMechCtx,surfaceEnergySet,ierr);CHKERRQ(ierr)
         Call MEF90DefMechCohesiveEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,cohesiveEnergySet,ierr);CHKERRQ(ierr)
         elasticEnergy(step)  = sum(elasticEnergySet)
         forceWork(step)      = sum(forceWorkSet)
         surfaceEnergy(step)  = sum(surfaceEnergySet)
         cohesiveEnergy(step) = sum(cohesiveEnergySet)
         totalMechanicalEnergy(step) = elasticEnergy(step) - forceWork(step) + cohesiveEnergy(step) + surfaceEnergy(step)

         !!!
         !!! Print and save energies
         !!!
         Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
         Call MEF90ISAllGatherMerge(MEF90Ctx%Comm,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
         Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
         Call PetscPrintf(MEF90Ctx%Comm,"\nMechanical energies: \n",ierr);CHKERRQ(ierr)
         Do set = 1, size(setID)
            Write(IOBuffer,201) setID(set),elasticEnergySet(set),forceWorkSet(set),cohesiveEnergySet(set),surfaceEnergySet(set),elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            Write(IOBuffer,500) step,time(step),elasticEnergySet(set),forceWorkSet(set),cohesiveEnergySet(set),surfaceEnergySet(set),elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set)
            Call PetscViewerASCIIPrintf(MEF90DefMechCtx%setEnergyViewer(set),IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscViewerFlush(MEF90DefMechCtx%setEnergyViewer(set),ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
         Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
         Write(IOBuffer,202) elasticEnergy(step),forceWork(step),cohesiveEnergy(step),surfaceEnergy(step),totalMechanicalEnergy(step)
         Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

         Write(IOBuffer,500) step,time(step),elasticEnergy(step),forceWork(step),cohesiveEnergy(step),surfaceEnergy(step),totalMechanicalEnergy(step)
         Call PetscViewerASCIIPrintf(MEF90DefMechCtx%globalEnergyViewer,IOBuffer,ierr);CHKERRQ(ierr)
         Call PetscViewerFlush(MEF90DefMechCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)
         Call PetscLogStagePop(ierr);CHKERRQ(ierr)

         Call PetscLogStagePush(logStageIO,ierr);CHKERRQ(ierr)
         !!!
         !!! Save results and boundary Values
         !!!
         If (MEF90DefMechGlobalOptions%stressOffset > 0) Then
            Call MEF90DefMechStress(MEF90DefMechCtx%displacement,MEF90DefMechCtx,MEF90DefMechCtx%stress,ierr)
         End If
         Call MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
         Call MEF90FlowViewEXO(MEF90FlowCtx,step,ierr)
         Call PetscLogStagePop(ierr);CHKERRQ(ierr)
         !!!
         !!! Save performance log file
         !!!
         Call PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr);CHKERRQ(ierr)
         Call PetscLogView(logViewer,ierr);CHKERRQ(ierr)
         Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)

         If (step == MEF90GlobalOptions%timeNumStep) Then
            EXIT
         Else
            step = step + 1
         End If
      End Do MainloopQS
   End If
   !!! Clean up and exit nicely
   Select case(MEF90DefMechGlobalOptions%timeSteppingType)
   Case (MEF90DefMech_timeSteppingTypeQuasiStatic)
      Call SNESDestroy(MEF90DefMechCtx%snesDisp,ierr);CHKERRQ(ierr)
      Call SNESDestroy(MEF90DefMechCtx%snesDamage,ierr);CHKERRQ(ierr)
      Call VecDestroy(residualDisp,ierr);CHKERRQ(ierr)
      Call VecDestroy(residualDamage,ierr);CHKERRQ(ierr)
      !!! NEW: Flow solver cleanup
      If (flowTransient) Then
         Call TSDestroy(tsPressure,ierr);CHKERRQ(ierr)
      End If
      Call SNESDestroy(snesPressure,ierr);CHKERRQ(ierr)
      Call VecDestroy(residualPressure,ierr);CHKERRQ(ierr)
   End Select
   
   Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
   Case (MEF90HeatXFer_timeSteppingTypeSteadyState) 
      Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
      call VecDestroy(residualTemp,ierr);CHKERRQ(ierr)
   Case (MEF90HeatXFer_timeSteppingTypeTransient) 
      Call TSDestroy(tsTemp,ierr);CHKERRQ(ierr)
   End Select

   Write(IOBuffer,*) 'Total number of alternate minimizations:',AltMinStep,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call MEF90DefMechCtxDestroyVectors(MEF90DefMechCtx,ierr)
   Call MEF90FlowCtxDestroyVectors(MEF90FlowCtx,ierr) !!! NEW
   Nullify(MEF90HeatXferCtx%temperature)
   Call MEF90HeatXferCtxDestroyVectors(MEF90HeatXferCtx,ierr)
   Call VecDestroy(damageAltMinOld,ierr);CHKERRQ(ierr)
   Call VecDestroy(displacementAltMinOld,ierr);CHKERRQ(ierr)

   DeAllocate(elasticEnergySet)
   DeAllocate(surfaceEnergySet)
   DeAllocate(forceWorkSet)
   DeAllocate(cohesiveEnergySet)
   DeAllocate(thermalEnergySet)
   DeAllocate(heatFluxWorkSet)
   
   DeAllocate(elasticEnergy)
   DeAllocate(forceWork)
   DeAllocate(cohesiveEnergy)
   DeAllocate(surfaceEnergy)
   DeAllocate(totalMechanicalEnergy)

   DeAllocate(time)
   Call MEF90DefMechCtxDestroy(MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90FlowCtxDestroy(MEF90FlowCtx,ierr)         !!! NEW
   Call MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
   Call MEF90CtxCloseEXO(MEF90Ctx,ierr)
   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)

   If (.NOT. MEF90GlobalOptions%dryrun) Then
      Call PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr);CHKERRQ(ierr)
      Call PetscLogView(logViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)
   End If
   Call MEF90CtxDestroy(MEF90Ctx,ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize(ierr)
100 Format("\nHeat transfer: solving steady state step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
110 Format("\nHeat transfer: step ",I4,", until t=",ES12.5,"\n")
200 Format("\nMechanics: step ",I4,", t=",ES12.5,"\n")
201 Format("cell set ",I4,"  elastic energy: ",ES12.5," work: ",ES12.5," cohesive: ",ES12.5," surface: ",ES12.5," total: ",ES12.5,"\n")
202 Format("======= Total: elastic energy: ",ES12.5," work: ",ES12.5," cohesive: ",ES12.5," surface: ",ES12.5," total: ",ES12.5,"\n")
208 Format("   Alt. Min. step ",I5," ")
209 Format(" alpha min / max", ES12.5, " / ", ES12.5, ", max change ", ES12.5,"\n")
400 Format(" [ERROR]: ",A," SNESSolve failed with SNESConvergedReason ",I2,". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
410 Format(" [ERROR]: ",A," TSSolve failed with TSConvergedReason ",I2,". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
500 Format(I6, 6(ES16.5),"\n")

End Program vDef
