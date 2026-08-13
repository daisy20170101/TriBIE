program quick_test
  use nikkhoo_walter
  use, intrinsic :: ieee_arithmetic
  implicit none

  real(DP), dimension(3) :: p1, p2, p3
  real(DP) :: x, y, z, ss, ds, ts, mu, lambda
  real(DP), dimension(6) :: stress, strain

  ! Triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Slip and elastic parameters
  ss = 1.0_DP; ds = -1.0_DP; ts = 2.0_DP
  mu = 3.0e10_DP; lambda = 3.0e10_DP

  print *, 'Quick test - Points 8 and 9'
  print *, ''

  ! Point 8
  x = 3.0_DP; y = -3.0_DP; z = -6.0_DP
  call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Point 8: Exx =', strain(1)
  print *, 'Expected:     7.064e-4'
  print *, 'Difference:', abs(strain(1) - 7.064e-4_DP)
  print *, ''

  ! Point 9
  x = -3.0_DP; y = 3.0_DP; z = -3.0_DP
  call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Point 9: Exx =', strain(1)
  print *, 'Expected:     2.113e-4'
  print *, 'Difference:', abs(strain(1) - 2.113e-4_DP)

end program quick_test
