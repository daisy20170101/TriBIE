program test_center
  use nikkhoo_walter
  implicit none

  integer, parameter :: DP = selected_real_kind(15, 307)

  ! Triangle vertices
  real(DP) :: p1(3), p2(3), p3(3)

  ! Test point - CENTER of triangle
  real(DP) :: x, y, z

  ! Slip components
  real(DP) :: ss, ds, ts

  ! Elastic parameters
  real(DP) :: mu, lambda

  ! Results
  real(DP) :: stress(6), strain(6)

  ! Initialize triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Center of triangle
  x = (p1(1) + p2(1) + p3(1)) / 3.0_DP
  y = (p1(2) + p2(2) + p3(2)) / 3.0_DP
  z = (p1(3) + p2(3) + p3(3)) / 3.0_DP

  print *, '=========================================='
  print *, 'Testing CENTER of Triangle'
  print *, '=========================================='
  print *, 'Triangle vertices:'
  print *, 'P1 = ', p1
  print *, 'P2 = ', p2
  print *, 'P3 = ', p3
  print *, ''
  print *, 'Center point:'
  print *, 'x = ', x
  print *, 'y = ', y
  print *, 'z = ', z
  print *, ''

  ! Slip components
  ss = 1.0_DP    ! Strike-slip
  ds = -1.0_DP   ! Dip-slip
  ts = 2.0_DP    ! Tensile-slip

  print *, 'Slip: ss=', ss, ' ds=', ds, ' ts=', ts

  ! Elastic parameters
  mu = 3.0e10_DP
  lambda = 3.0e10_DP

  print *, 'Elastic: mu=', mu, ' lambda=', lambda
  print *, ''

  ! Call tdstress_hs
  call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)

  print *, '=========================================='
  print *, 'RESULTS for CENTER:'
  print *, '=========================================='
  print *, 'Strain:'
  print *, '  Exx = ', strain(1)
  print *, '  Eyy = ', strain(2)
  print *, '  Ezz = ', strain(3)
  print *, '  Exy = ', strain(4)
  print *, '  Exz = ', strain(5)
  print *, '  Eyz = ', strain(6)
  print *, ''
  print *, 'Stress:'
  print *, '  Sxx = ', stress(1)
  print *, '  Syy = ', stress(2)
  print *, '  Szz = ', stress(3)
  print *, '  Sxy = ', stress(4)
  print *, '  Sxz = ', stress(5)
  print *, '  Syz = ', stress(6)
  print *, ''

  ! Check for NaN
  if (strain(1) /= strain(1)) then
    print *, 'ERROR: Exx is NaN for CENTER of triangle!'
    print *, 'This should NEVER happen - the center should be well-defined.'
  else
    print *, 'SUCCESS: Exx is finite for center'
  end if

end program test_center
