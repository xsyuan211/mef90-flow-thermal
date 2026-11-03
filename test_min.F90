module m_test
   implicit none
   integer, parameter :: Kr = selected_real_kind(15)
   integer, parameter :: PETSC_FALSE = 0

   type :: T
      integer :: a
      integer :: b
      real(Kr) :: c
   end type

   type(T), parameter :: v = T( &
      1, &
      PETSC_FALSE, &
      1.0e-20_Kr)

end module
