module m_test
   implicit none
   integer, parameter :: Kr = selected_real_kind(15)

   type :: T
      integer :: a1, a2
      real(Kr) :: r1, r2, r3, r4, r5, r6, r7
      integer :: i1, i2, i3, i4, i5, i6, i7, i8
      real(Kr) :: r8, r9
      integer :: i9
   end type

   type(T), parameter :: v = T( &
      1, & ! 1
      0, & ! 2
      1.0e-20_Kr, & ! 3
      1.0e-14_Kr, & ! 4
      1.0e-3_Kr, & ! 5
      4.6e-10_Kr, & ! 6
      0.0e0_Kr, & ! 7
      0.0e0_Kr, & ! 8
      0.0e0_Kr, & ! 9
      1, & ! 10
      1, & ! 11
      0, & ! 12
      2, & ! 13
      0, & ! 14
      0, & ! 15
      1, & ! 16
      0, & ! 17
      1e-8_Kr, & ! 18
      1e-4_Kr, & ! 19
      100) ! 20

end module
