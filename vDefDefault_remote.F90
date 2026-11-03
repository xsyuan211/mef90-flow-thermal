#include "../MEF90/mef90.inc"
module m_vDefDefault
#include <finclude/petscdef.h>
   Use petsc
   Use m_MEF90
   Use m_MEF90_Parameters, Only: Kr   !!! explicit kind parameter for real literals
   Use m_MEF90_DefMechCtx
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_FlowCtx
   Implicit NONE   


!!! Default values of the contexts
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: vDefDefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         PETSC_FALSE,                   & ! validate
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0e0_Kr,                        & ! timeMin
                                                         1.0e0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         0,                             & ! timeSkip
                                                         1,                             & ! numCycle
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat

   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: vDefDefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_TimeSteppingTypeQuasiStatic, & ! solverType
                                                         MEF90DefMech_SolverTypeAltMin,            & ! timeSteppingType
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! CrackPressureOffset
                                                         0,                       & ! plasticStrainOffset
                                                         6,                       & ! StressOffset
                                                         9,                       & ! cumulatedPlasticDissipationOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! displacementLowerBoundScaling
                                                         MEF90Scaling_CST,        & ! displacementUpperBoundScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         MEF90Scaling_Linear,     & ! CrackPressureScaling
                                                         1e-3,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         10,                      & ! PCLag
                                                         1.0e0_Kr,                  & ! SOROmega
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         1.0e-3,                  & ! InjectedVolumeAtol
                                                         0.0e0_Kr,                  & ! dampingCoefficientDisplacement
                                                         0.0e0_Kr)                    ! dampingCoefficientDamage

   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: vDefDefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_TimeSteppingTypeQuasiStatic, & ! solverType
                                                         MEF90DefMech_SolverTypeAltMin,            & ! timeSteppingType
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! CrackPressureOffset
                                                         0,                       & ! plasticStrainOffset
                                                         7,                       & ! StressOffset
                                                         13,                      & ! cumulatedPlasticDissipationOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! displacementLowerBoundScaling
                                                         MEF90Scaling_CST,        & ! displacementUpperBoundScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         MEF90Scaling_Linear,     & ! CrackPressureScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         10,                      & ! PCLag
                                                         1.0e0_Kr,                  & ! SOROmega
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         1.0e-3,                  & ! InjectedVolumeAtol
                                                         0.0e0_Kr,                  & ! dampingCoefficientDisplacement
                                                         0.0e0_Kr)                    ! dampingCoefficientDamage

   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: vDefDefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                                & ! elemTypeShortIDDispl will be overriden
                                                         -1,                                                & ! elemTypeShortIDDamage will be overriden
                                                         [0.0e0_Kr,0.0e0_Kr,0.0e0_Kr],                            & ! force
                                                         0.0e0_Kr,                                            & ! pressureForce
                                                         0.0e0_Kr,                                            & ! CrackPressure
                                                         MEF90DefMech_damageTypeAT1,                        & ! damageType
                                                         MEF90DefMech_plasticityTypeNone,                   & ! plasticityType
                                                         MEF90DefMech_unilateralContactTypeNone,            & ! unilateralContactType
                                                         1.0e-5_Kr,                                            & ! unilateralContactHydrostatocDeviatoricGamma
                                                         PETSC_FALSE,                                       & ! unilateralContactHybrid
                                                         1.0e0_Kr,                                            & ! DamageATLinSoftk
                                                         1.25e0_Kr,                                           & ! DamageAT1expb
                                                         MEF90DefMech_drivingForceTypeNone,                 & ! drivingForceType
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],             & ! Has Displacement BC
                                                         [0.0e0_Kr,0.0e0_Kr,0.0e0_Kr],                            & ! boundary Displacement
                                                         [MEF90_NINFINITY,MEF90_NINFINITY,MEF90_NINFINITY], & ! displacement Lower Bound
                                                         [MEF90_INFINITY,MEF90_INFINITY,MEF90_INFINITY],    & ! displacement Upper Bound
                                                         PETSC_FALSE,                                       & ! Has Damage BC
                                                         PETSC_FALSE,                                       & ! CrackVolumeControlled
                                                         PETSC_FALSE,                                       & ! WorkControlled
                                                         0.e0_Kr)                                               ! Boundary Damage
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: vDefDefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],             & ! Has Displacement BC
                                                         [0.0e0_Kr,0.0e0_Kr,0.0e0_Kr],                            & ! boundary Displacement
                                                         [MEF90_NINFINITY,MEF90_NINFINITY,MEF90_NINFINITY], & ! displacement Lower Bound
                                                         [MEF90_INFINITY,MEF90_INFINITY,MEF90_INFINITY],    & ! displacement Upper Bound
                                                         PETSC_FALSE,                                       & ! Has Damage BC
                                                         0.0e0_Kr)                                              ! boundary Damage

   Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: vDefHeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                         MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
                                                         PETSC_FALSE,         & ! addNullSpace
                                                         1,                   & ! tempOffset
                                                         0.,                  & ! initialTemperature
                                                         MEF90Scaling_Linear, & ! boundaryTempScaling
                                                         0,                   & ! boundaryTempOffset
                                                         MEF90Scaling_Linear, & ! externalTempScaling
                                                         2,                   & ! externalTempOffset
                                                         MEF90Scaling_Linear, & ! fluxScaling
                                                         1)                     ! fluxOffset
   Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: vDefHeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                         -1,                   & ! elemTypeShortID will be overriden
                                                         0.0e0_Kr,               & ! flux
                                                         0.0e0_Kr,               & ! surfaceThermalConductivity
                                                         0.0e0_Kr,               & ! externalTemp
                                                         PETSC_FALSE,          & ! Has BC
                                                         0.0e0_Kr,               & ! boundaryTemp
                                                         [0.0e0_Kr,0.0e0_Kr,0.0e0_Kr]) ! advectionVelocity
                                                         
   Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: vDefHeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0e0_Kr)          ! boundaryTemp

!!! NEW: Flow (Darcy pressure) default options -- Phase 1
   Type(MEF90FlowGlobalOptions_Type),Parameter    :: vDefFlowDefaultGlobalOptions = MEF90FlowGlobalOptions_Type( &
                                                         MEF90Flow_timeSteppingTypeSteadyState, & ! timeSteppingType (steady Darcy)
                                                          PETSC_FALSE,          & ! addNullSpace
                                                         1.0e-20_Kr,     & ! permeabilityMin  (k_min, 1e-20 m^2 ~ 10 mD)
                                                         1.0e-14_Kr,     & ! permeabilityMax  (k_max, 1e-14 m^2 ~ 1000 D)
                                                         1.0e-3_Kr,      & ! dynamicViscosity (mu, water at 20degC)
                                                         4.6e-10_Kr,     & ! storativity S = (1/Kf + mu/Ks) ~ 4.6e-10 Pa^-1
                                                         0.0e0_Kr,         & ! initialPressure p0
                                                         0.0e0_Kr,         & ! residualDamage (damage < this -> k = k_min only)
                                                         0.0e0_Kr,         & ! biotCoefficient alpha (0 = no Biot coupling, set 1.0 for full poroelasticity)
                                                         1,              & ! pressureOffset (EXO)
                                                         1,              & ! sourceOffset
                                                         0,              & ! boundaryPressureOffset
                                                         2,              & ! pressurePreviousStepOffset
                                                         MEF90Scaling_CST, & ! sourceScaling
                                                         MEF90Scaling_CST, & ! boundaryPressureScaling
                                                         MEF90Scaling_Linear, & ! pressurePreviousStepScaling
                                                         MEF90Scaling_CST,    & ! pressureScaling
                                                         1e-8_Kr,        & ! pressureATol
                                                         1e-4_Kr,        & ! pressureRTol
                                                         100)             & ! pressureMaxIt
   Type(MEF90FlowCellSetOptions_Type),Parameter   :: vDefFlowDefaultCellSetOptions = MEF90FlowCellSetOptions_Type( &
                                                         -1,              & ! elemTypeShortID (overridden)
                                                         0.0e0_Kr,          & ! source (volumetric injection/sink)
                                                         0.0e0_Kr,          & ! permeabilityOverride (0 = use global)
                                                         -1.0e0_Kr,         & ! storativityOverride (-1 = use global)
                                                         PETSC_FALSE,     & ! Has BC
                                                         0.0e0_Kr,          & ! boundaryPressure
                                                         [0.0e0_Kr,0.0e0_Kr,0.0e0_Kr]) ! advectionVector (gravity)
   Type(MEF90FlowVertexSetOptions_Type),Parameter :: vDefFlowDefaultVertexSetOptions = MEF90FlowVertexSetOptions_Type( &
                                                         PETSC_FALSE,     & ! Has BC
                                                         0.0e0_Kr)          & ! boundaryPressure

end module m_vDefDefault

