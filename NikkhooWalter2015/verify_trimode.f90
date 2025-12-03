program verify_trimode_classification
  use nikkhoo_walter
  implicit none

  ! Triangle vertices
  real(DP) :: p1(3), p2(3), p3(3)
  real(DP) :: ss, ds, ts, mu, lambda
  real(DP) :: stress(6), strain(6)

  ! Initialize
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ss = 1.0_DP
  ds = -1.0_DP
  ts = 2.0_DP
  mu = 3.0e10_DP
  lambda = 3.0e10_DP

  print *, '============================================================'
  print *, 'Testing Points 4, 5, 12, 15 - Expected FINITE values'
  print *, '============================================================'
  print *, ''
  print *, 'Expected results from reference:'
  print *, '  Point 4:  Exx = 0.000829157341339727'
  print *, '  Point 5:  Exx = 0.00114439668841158'
  print *, '  Point 12: Exx = 0.00441202690885827'
  print *, '  Point 15: Exx = -0.000914111766849476'
  print *, ''

  ! Point 4: (7.0, -1.0, -5.0)
  print *, '------------------------------------------------------------'
  print *, 'Point 4: (7.0, -1.0, -5.0)'
  print *, 'Expected: Exx = 0.000829157341339727'
  call tdstress_hs(7.0_DP, -1.0_DP, -5.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Got:      Exx =', strain(1)
  print *, ''

  ! Point 5: (-7.0, -1.0, -5.0)
  print *, '------------------------------------------------------------'
  print *, 'Point 5: (-7.0, -1.0, -5.0)'
  print *, 'Expected: Exx = 0.00114439668841158'
  call tdstress_hs(-7.0_DP, -1.0_DP, -5.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Got:      Exx =', strain(1)
  print *, ''

  ! Point 12: (1.0, -1.0, -1.0)
  print *, '------------------------------------------------------------'
  print *, 'Point 12: (1.0, -1.0, -1.0)'
  print *, 'Expected: Exx = 0.00441202690885827'
  call tdstress_hs(1.0_DP, -1.0_DP, -1.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Got:      Exx =', strain(1)
  print *, ''

  ! Point 15: (1.0, -1.0, -8.0)
  print *, '------------------------------------------------------------'
  print *, 'Point 15: (1.0, -1.0, -8.0)'
  print *, 'Expected: Exx = -0.000914111766849476'
  call tdstress_hs(1.0_DP, -1.0_DP, -8.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Got:      Exx =', strain(1)
  print *, ''

  print *, '============================================================'
  print *, 'Point 1 (Center): (-0.333, -0.333, -4.667)'
  print *, 'Expected: Exx = 0.0481047005255181'
  call tdstress_hs(-1.0_DP/3.0_DP, -1.0_DP/3.0_DP, -14.0_DP/3.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Got:      Exx =', strain(1)
  print *, ''

end program verify_trimode_classification
