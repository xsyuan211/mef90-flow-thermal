module m_test
   implicit none
   integer, parameter :: Kr = selected_real_kind(15)
   integer, parameter :: PETSC_FALSE = 0
   integer, parameter :: CST = 1
   integer, parameter :: LIN = 2

   type :: T
      integer :: a1, a2
      real(Kr) :: r1, r2, r3, r4, r5, r6, r7
      integer :: i1, i2, i3, i4, i5, i6, i7, i8
      real(Kr) :: r8, r9
      integer :: i9
   end type

   type(T), parameter :: v = T( &
      1, &
      PETSC_FALSE, &
      1.0e-20_Kr, &
      1.0e-14_Kr, &
      1.0e-3_Kr, &
      4.6e-10_Kr, &
      0.0e0_Kr, &
      0.0e0_Kr, &
      0.0e0_Kr, &
      1, &
      1, &
      0, &
      2, &
      CST, &
      CST, &
      LIN, &
      CST, &
      1e-8_Kr, &
      1e-4_Kr, &
      100)

end module
