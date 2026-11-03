module m_test
   implicit none
   integer, parameter :: Kr = selected_real_kind(15)

   integer, parameter :: MEF90Flow_timeSteppingTypeSteadyState = 1
   integer, parameter :: PETSC_FALSE = 0
   integer, parameter :: MEF90Scaling_CST = 1
   integer, parameter :: MEF90Scaling_Linear = 2

   type :: MEF90FlowGlobalOptions_Type
      integer :: timeSteppingType
      integer :: addNullSpace
      real(Kr) :: permeabilityMin
      real(Kr) :: permeabilityMax
      real(Kr) :: dynamicViscosity
      real(Kr) :: storativity
      real(Kr) :: initialPressure
      real(Kr) :: residualDamage
      real(Kr) :: biotCoefficient
      integer  :: pressureOffset
      integer  :: sourceOffset
      integer  :: boundaryPressureOffset
      integer  :: pressurePreviousStepOffset
      integer  :: sourceScaling
      integer  :: boundaryPressureScaling
      integer  :: pressurePreviousStepScaling
      integer  :: pressureScaling
      real(Kr) :: pressureATol
      real(Kr) :: pressureRTol
      integer  :: pressureMaxIt
   end type

   type(MEF90FlowGlobalOptions_Type), parameter :: v = MEF90FlowGlobalOptions_Type( &
      MEF90Flow_timeSteppingTypeSteadyState, & ! timeSteppingType
      PETSC_FALSE,          & ! addNullSpace
      1.0e-20_Kr,     & ! permeabilityMin
      1.0e-14_Kr,     & ! permeabilityMax
      1.0e-3_Kr,      & ! dynamicViscosity
      4.6e-10_Kr,     & ! storativity
      0.0e0_Kr,         & ! initialPressure
      0.0e0_Kr,         & ! residualDamage
      0.0e0_Kr,         & ! biotCoefficient
      1,              & ! pressureOffset
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

end module
