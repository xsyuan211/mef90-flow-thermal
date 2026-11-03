module m_test
   implicit none
   integer, parameter :: Kr = selected_real_kind(15)
   integer, parameter :: PETSC_FALSE = 0

   type :: T2
      integer :: a
      real(Kr) :: b
   end type

   type(T2), parameter :: v = T2( &
      PETSC_FALSE,     & ! Has BC
      0.0e0_Kr)          & ! boundaryPressure

end module
