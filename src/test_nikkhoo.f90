!==============================================================================
! test_nikkhoo.f90
! 
! Test program for the Nikkhoo & Walter (2015) triangular dislocation method
!==============================================================================

program test_nikkhoo
  use nikkhoo_walter
  implicit none
  
  ! Test parameters
  integer, parameter :: n_points = 12  ! 4 original + 8 new test points
  real(DP), dimension(n_points) :: x, y, z
  real(DP), dimension(3) :: p1, p2, p3
  real(DP) :: ss, ds, ts, mu, lambda
  real(DP), dimension(n_points, 6) :: stress, strain
  
  ! Initialize test data
  ! Original test points
  x = [-1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
       3.0_DP, -3.0_DP, -1.0_DP, -1.0_DP, 1.0_DP, -1.0_DP, -1.0_DP, 1.0_DP]
  y = [-1.0_DP, -1.0_DP, 1.0_DP, -1.0_DP, &
       -3.0_DP, 3.0_DP, -1.0_DP, 1.0_DP, -1.0_DP, -1.0_DP, 1.0_DP, -1.0_DP]
  z = [-1.0_DP, -1.0_DP, -1.0_DP, -2.0_DP, &
       -6.0_DP, -3.0_DP, -1.0_DP, -1.0_DP, -1.0_DP, -8.0_DP, -8.0_DP, -8.0_DP]
  
  p1 = [-1.0_DP, 0.0_DP, -1.0_DP]
  p2 = [1.0_DP, 0.0_DP, -1.0_DP]
  p3 = [0.0_DP, 1.0_DP, -1.0_DP]
  
  ss = 1.0_DP  ! Strike-slip
  ds = 0.5_DP  ! Dip-slip
  ts = 0.0_DP  ! Tensile-slip
  
  mu = 3.0e10_DP      ! Shear modulus
  lambda = 3.0e10_DP  ! Lame's first parameter
  
  ! Calculate stresses and strains
  call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain, n_points)
  

end program test_nikkhoo


