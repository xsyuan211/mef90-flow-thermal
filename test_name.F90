module m_test
   implicit none
   integer, parameter :: Kr = selected_real_kind(15)
   integer, parameter :: PETSC_FALSE = 0
   integer, parameter :: MEF90Flow_timeSteppingTypeSteadyState = 1

   type :: T
      integer :: a1, a2
      real(Kr) :: r1, r2, r3, r4, r5, r6, r7
      integer :: i1, i2, i3, i4, i5, i6, i7, i8
      real(Kr) :: r8, r9
      integer :: i9
   end type

   type(T), parameter :: v = T( &
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
      1, &
      1, &
      2, &
      1, &
      1e-8_Kr,        & ! pressureATol
      1e-4_Kr,        & ! pressureRTol
      100)

end module
