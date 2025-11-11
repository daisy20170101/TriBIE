program debug_trimode
  use nikkhoo_walter
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Triangle vertices
  real(DP), dimension(3) :: p1, p2, p3

  ! Test points 8 and 9
  real(DP) :: x8, y8, z8, x9, y9, z9
  real(DP), dimension(6) :: stress, strain

  ! Slip and elastic parameters
  real(DP) :: ss, ds, ts, mu, lambda

  ! Set up triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Set up test points 8 and 9
  x8 = 3.0_DP
  y8 = -3.0_DP
  z8 = -6.0_DP

  x9 = -3.0_DP
  y9 = 3.0_DP
  z9 = -3.0_DP

  ! Set slip components
  ss = 1.0_DP
  ds = -1.0_DP
  ts = 2.0_DP

  ! Set elastic parameters
  mu = 3.0e10_DP
  lambda = 3.0e10_DP

  print *, '======================================='
  print *, 'Debug trimode for points 8 and 9'
  print *, '======================================='
  print *, ''

  ! Test point 8
  print *, 'Point 8: x=', x8, ', y=', y8, ', z=', z8
  call tdstress_hs(x8, y8, z8, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)

  print *, '  Strain(1) [Exx]:', strain(1)
  if (ieee_is_nan(strain(1))) then
    print *, '  --> NaN detected! Point classified as on-edge (trimode=0)'
  end if
  print *, ''

  ! Test point 9
  print *, 'Point 9: x=', x9, ', y=', y9, ', z=', z9
  call tdstress_hs(x9, y9, z9, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)

  print *, '  Strain(1) [Exx]:', strain(1)
  if (ieee_is_nan(strain(1))) then
    print *, '  --> NaN detected! Point classified as on-edge (trimode=0)'
  end if
  print *, ''

  print *, '======================================='
  print *, 'Analysis:'
  print *, '  If NaN appears, the points are being'
  print *, '  incorrectly classified as trimode=0'
  print *, '  (on triangle edge) due to floating-'
  print *, '  point precision in barycentric coords.'
  print *, '======================================='

end program debug_trimode
