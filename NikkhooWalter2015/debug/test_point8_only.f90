program test_point8_only
  use nikkhoo_walter
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Triangle vertices
  real(DP), dimension(3) :: p1, p2, p3

  ! Test point 8
  real(DP) :: x, y, z
  real(DP) :: ss, ds, ts, mu, lambda
  real(DP), dimension(6) :: stress, strain

  ! Set up triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Set slip components
  ss = 1.0_DP
  ds = -1.0_DP
  ts = 2.0_DP

  ! Set elastic parameters
  mu = 3.0e10_DP
  lambda = 3.0e10_DP

  print *, '============================================================'
  print *, 'Testing POINT 8 ONLY with full debug output'
  print *, '============================================================'
  print *, ''

  ! Test point 8: (3.0, -3.0, -6.0)
  x = 3.0_DP
  y = -3.0_DP
  z = -6.0_DP

  print *, 'Point 8 coords: x=', x, ', y=', y, ', z=', z
  print *, ''
  print *, 'Calling tdstress_hs...'
  print *, ''

  call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)

  print *, ''
  print *, '============================================================'
  print *, 'RESULTS for Point 8'
  print *, '============================================================'
  print *, 'Strain(1) Exx = ', strain(1)
  if (ieee_is_nan(strain(1))) then
    print *, '  ✗ FAIL: Got NaN'
  else
    print *, '  ✓ PASS: Got valid number'
  end if
  print *, ''

end program test_point8_only
