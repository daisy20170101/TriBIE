program debug_specific_points
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
  print *, 'Testing Points 4, 5, 12, 15 that have NaN'
  print *, '============================================================'
  print *, ''
  print *, 'Triangle vertices (EFCS):'
  print *, 'P1 = ', p1
  print *, 'P2 = ', p2
  print *, 'P3 = ', p3
  print *, ''
  print *, 'Note: P2 is at (1.0, -1.0, -5.0)'
  print *, ''

  ! Point 4: (7.0, -1.0, -5.0)
  print *, '------------------------------------------------------------'
  print *, 'Point 4: (7.0, -1.0, -5.0)'
  print *, 'This has same Y and Z as P2, but X=7.0 vs P2 X=1.0'
  print *, 'So it is 6 units away from P2 in X direction'
  print *, 'This should project to P2 in the triangle plane (y_td=0, z_td=0)'
  print *, 'Such points are singular (on line through vertex perpendicular to triangle)'
  print *, ''
  call tdstress_hs(7.0_DP, -1.0_DP, -5.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Result: Exx = ', strain(1)
  if (strain(1) /= strain(1)) then
    print *, 'NaN is EXPECTED for this point (singular configuration)'
  end if
  print *, ''

  ! Point 5: (-7.0, -1.0, -5.0)
  print *, '------------------------------------------------------------'
  print *, 'Point 5: (-7.0, -1.0, -5.0)'
  print *, 'This has same Y and Z as P1, but X=-7.0 vs P1 X=-1.0'
  print *, 'So it is 6 units away from P1 in X direction'
  print *, 'This should project to P1 in the triangle plane'
  print *, 'Such points are singular (on line through vertex perpendicular to triangle)'
  print *, ''
  call tdstress_hs(-7.0_DP, -1.0_DP, -5.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Result: Exx = ', strain(1)
  if (strain(1) /= strain(1)) then
    print *, 'NaN is EXPECTED for this point (singular configuration)'
  end if
  print *, ''

  ! Point 12: (1.0, -1.0, -1.0)
  print *, '------------------------------------------------------------'
  print *, 'Point 12: (1.0, -1.0, -1.0)'
  print *, 'This has same X and Y as P2, but Z=-1.0 vs P2 Z=-5.0'
  print *, 'So it is 4 units above P2 in Z direction'
  print *, 'This should also project to P2 in the triangle plane'
  print *, ''
  call tdstress_hs(1.0_DP, -1.0_DP, -1.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Result: Exx = ', strain(1)
  if (strain(1) /= strain(1)) then
    print *, 'NaN might be expected (depends on if this is singular)'
  end if
  print *, ''

  ! Point 15: (1.0, -1.0, -8.0)
  print *, '------------------------------------------------------------'
  print *, 'Point 15: (1.0, -1.0, -8.0)'
  print *, 'This has same X and Y as P2, but Z=-8.0 vs P2 Z=-5.0'
  print *, 'So it is 3 units below P2 in Z direction'
  print *, 'This should also project to P2 in the triangle plane'
  print *, ''
  call tdstress_hs(1.0_DP, -1.0_DP, -8.0_DP, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
  print *, 'Result: Exx = ', strain(1)
  if (strain(1) /= strain(1)) then
    print *, 'NaN might be expected (depends on if this is singular)'
  end if
  print *, ''

  print *, '============================================================'
  print *, 'CONCLUSION:'
  print *, 'Points 4, 5, 12, 15 lie on lines perpendicular to the'
  print *, 'triangle passing through vertices. These are genuinely'
  print *, 'singular configurations in the triangular dislocation'
  print *, 'solution, so NaN is the correct answer.'
  print *, ''
  print *, 'We need to verify this with MATLAB to confirm.'
  print *, '============================================================'

end program debug_specific_points
