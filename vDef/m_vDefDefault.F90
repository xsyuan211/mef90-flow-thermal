#include "../MEF90/mef90.inc"
module m_vDefDefault
#include <finclude/petscdef.h>
   Use petsc
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_FlowCtx   !!! NEW
   Implicit NONE   


!!! Default values of the contexts
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: vDefDefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         PETSC_FALSE,                   & ! validate
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
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
                                                         1.0_Kr,                  & ! SOROmega
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         1.0e-3,                  & ! InjectedVolumeAtol
                                                         0.0_Kr,                  & ! dampingCoefficientDisplacement
                                                         0.0_Kr)                    ! dampingCoefficientDamage

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
                                                         1.0_Kr,                  & ! SOROmega
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         1.0e-3,                  & ! InjectedVolumeAtol
                                                         0.0_Kr,                  & ! dampingCoefficientDisplacement
                                                         0.0_Kr)                    ! dampingCoefficientDamage

   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: vDefDefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                                & ! elemTypeShortIDDispl will be overriden
                                                         -1,                                                & ! elemTypeShortIDDamage will be overriden
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                            & ! force
                                                         0.0_Kr,                                            & ! pressureForce
                                                         0.0_Kr,                                            & ! CrackPressure
                                                         MEF90DefMech_damageTypeAT1,                        & ! damageType
                                                         MEF90DefMech_plasticityTypeNone,                   & ! plasticityType
                                                         MEF90DefMech_unilateralContactTypeNone,            & ! unilateralContactType
                                                         1.0D-5,                                            & ! unilateralContactHydrostatocDeviatoricGamma
                                                         PETSC_FALSE,                                       & ! unilateralContactHybrid
                                                         1.0_Kr,                                            & ! DamageATLinSoftk
                                                         1.25_Kr,                                           & ! DamageAT1expb
                                                         MEF90DefMech_drivingForceTypeNone,                 & ! drivingForceType
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],             & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                            & ! boundary Displacement
                                                         [MEF90_NINFINITY,MEF90_NINFINITY,MEF90_NINFINITY], & ! displacement Lower Bound
                                                         [MEF90_INFINITY,MEF90_INFINITY,MEF90_INFINITY],    & ! displacement Upper Bound
                                                         PETSC_FALSE,                                       & ! Has Damage BC
                                                         PETSC_FALSE,                                       & ! CrackVolumeControlled
                                                         PETSC_FALSE,                                       & ! WorkControlled
                                                         0._Kr)                                               ! Boundary Damage
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: vDefDefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],             & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                            & ! boundary Displacement
                                                         [MEF90_NINFINITY,MEF90_NINFINITY,MEF90_NINFINITY], & ! displacement Lower Bound
                                                         [MEF90_INFINITY,MEF90_INFINITY,MEF90_INFINITY],    & ! displacement Upper Bound
                                                         PETSC_FALSE,                                       & ! Has Damage BC
                                                         0.0_Kr)                                              ! boundary Damage

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
                                                         0.0_Kr,               & ! flux
                                                         0.0_Kr,               & ! surfaceThermalConductivity
                                                         0.0_Kr,               & ! externalTemp
                                                         PETSC_FALSE,          & ! Has BC
                                                         0.0_Kr,               & ! boundaryTemp
                                                         [0.0_Kr,0.0_Kr,0.0_Kr]) ! advectionVelocity
                                                         
   Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: vDefHeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp

!!! NEW: Flow (Darcy pressure) default options - Phase 1
!!! NOTE: there must be NO ampersand after the closing parenthesis (a trailing
!!!       & would splice the next Type declaration into this constructor).
!!!       Decimal literals are used instead of 1.0e-20_Kr style constants.
!!!       Value order must match MEF90FlowGlobalOptions_Type (29 components):
!!!       timeSteppingType, addNullSpace, permeabilityMin, permeabilityMax,
!!!       dynamicViscosity, storativity, initialPressure, residualDamage,
!!!       biotCoefficient, temperatureDependViscosity, referenceTemperature,
!!!       sutherlandConstant, arrheniusActivationEnergy, useArrheniusViscosity,
!!!       enableFlowAdvection, fluidDensity, pressureOffset, sourceOffset,
!!!       boundaryPressureOffset, pressurePreviousStepOffset, permFieldOffset,
!!!       darcyVelocityOffset, sourceScaling, boundaryPressureScaling,
!!!       pressurePreviousStepScaling, pressureScaling, pressureATol,
!!!       pressureRTol, pressureMaxIt.
   Type(MEF90FlowGlobalOptions_Type),Parameter    :: vDefFlowDefaultGlobalOptions = MEF90FlowGlobalOptions_Type( &
                                                         MEF90Flow_timeSteppingTypeSteadyState, &
                                                         PETSC_FALSE,    &
                                                         0.00000000000000000001_Kr, &  ! k_min = 1e-20 m^2
                                                         0.00000000000001_Kr, &        ! k_max = 1e-14 m^2
                                                         0.001_Kr,      &                ! mu_ref = 1e-3 Pa.s
                                                         0.00000000046_Kr, &            ! S = 4.6e-10 1/Pa
                                                         0.0_Kr,         &                ! p0
                                                         0.0_Kr,         &                ! residualDamage
                                                         0.0_Kr,         &                ! biotCoefficient
                                                         PETSC_FALSE,    &                ! temperatureDependViscosity
                                                         293.0_Kr,       &                ! T_ref
                                                         600.0_Kr,       &                ! Sutherland C
                                                         5100.0_Kr,      &                ! Ea/R
                                                         PETSC_FALSE,    &                ! useArrheniusViscosity
                                                         PETSC_FALSE,    &                ! enableFlowAdvection
                                                         1000.0_Kr,      &                ! rho_f
                                                         1,              &                ! pressureOffset
                                                         1,              &                ! sourceOffset
                                                         0,              &                ! boundaryPressureOffset
                                                         2,              &                ! pressurePreviousStepOffset
                                                         0,              &                ! permFieldOffset
                                                         0,              &                ! darcyVelocityOffset
                                                         MEF90Scaling_CST, &               ! sourceScaling
                                                         MEF90Scaling_CST, &               ! boundaryPressureScaling
                                                         MEF90Scaling_Linear, &            ! pressurePreviousStepScaling
                                                         MEF90Scaling_CST,    &            ! pressureScaling
                                                         0.00000001_Kr,  &                ! pressureATol = 1e-8
                                                         0.0001_Kr,      &                ! pressureRTol = 1e-4
                                                         100)              ! pressureMaxIt
   Type(MEF90FlowCellSetOptions_Type),Parameter   :: vDefFlowDefaultCellSetOptions = MEF90FlowCellSetOptions_Type( &
                                                         -1,              & ! elemTypeShortID (overridden)
                                                         0.0_Kr,          & ! source (volumetric injection/sink)
                                                         0.0_Kr,          & ! permeabilityOverride (0 = use global)
                                                         -1.0_Kr,         & ! storativityOverride (-1 = use global)
                                                         PETSC_FALSE,     & ! Has BC
                                                         0.0_Kr,          & ! boundaryPressure
                                                         [0.0_Kr,0.0_Kr,0.0_Kr]) ! advectionVector (gravity)
   Type(MEF90FlowVertexSetOptions_Type),Parameter :: vDefFlowDefaultVertexSetOptions = MEF90FlowVertexSetOptions_Type( &
                                                         PETSC_FALSE,     & ! Has BC
                                                         0.0_Kr)            ! boundaryPressure

end module m_vDefDefault