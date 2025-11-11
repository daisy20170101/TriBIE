program test_multiple_points
  use nikkhoo_walter
  implicit none

  integer, parameter :: DP = selected_real_kind(15, 307)
  integer :: i, n_pass, n_fail

  ! Triangle vertices
  real(DP) :: p1(3), p2(3), p3(3)

  ! Test points
  integer, parameter :: N_POINTS = 7
  real(DP) :: test_x(N_POINTS), test_y(N_POINTS), test_z(N_POINTS)
  character(len=50) :: test_names(N_POINTS)

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

  ! Define test points
  ! Point 1: Center
  test_x(1) = (p1(1) + p2(1) + p3(1)) / 3.0_DP
  test_y(1) = (p1(2) + p2(2) + p3(2)) / 3.0_DP
  test_z(1) = (p1(3) + p2(3) + p3(3)) / 3.0_DP
  test_names(1) = 'Center of triangle'

  ! Point 2: Well inside, offset from center
  test_x(2) = -0.5_DP
  test_y(2) = -0.5_DP
  test_z(2) = -4.8_DP
  test_names(2) = 'Inside triangle (offset from center)'

  ! Point 3: Far below triangle
  test_x(3) = 0.0_DP
  test_y(3) = 0.0_DP
  test_z(3) = -10.0_DP
  test_names(3) = 'Far below triangle'

  ! Point 4: Far to the side
  test_x(4) = 5.0_DP
  test_y(4) = 0.0_DP
  test_z(4) = -5.0_DP
  test_names(4) = 'Far to the side'

  ! Point 5: Above the triangle but still in half-space
  test_x(5) = 0.0_DP
  test_y(5) = 0.0_DP
  test_z(5) = -3.0_DP
  test_names(5) = 'Above triangle'

  ! Point 6: Near vertex P1 but not exactly on it
  test_x(6) = -0.9_DP
  test_y(6) = -0.9_DP
  test_z(6) = -4.9_DP
  test_names(6) = 'Near vertex P1'

  ! Point 7: Near edge P1-P2 but not on it
  test_x(7) = 0.0_DP
  test_y(7) = -0.95_DP
  test_z(7) = -4.9_DP
  test_names(7) = 'Near edge P1-P2'

  ! Slip components
  ss = 1.0_DP
  ds = -1.0_DP
  ts = 2.0_DP

  ! Elastic parameters
  mu = 3.0e10_DP
  lambda = 3.0e10_DP

  print *, '============================================================'
  print *, 'Testing Multiple Points'
  print *, '============================================================'
  print *, 'Triangle vertices:'
  print *, 'P1 = ', p1
  print *, 'P2 = ', p2
  print *, 'P3 = ', p3
  print *, ''
  print *, 'Parameters: ss=', ss, ' ds=', ds, ' ts=', ts
  print *, '            mu=', mu, ' lambda=', lambda
  print *, ''

  n_pass = 0
  n_fail = 0

  do i = 1, N_POINTS
    print *, '------------------------------------------------------------'
    print *, 'Point', i, ': ', trim(test_names(i))
    print *, 'Position: (', test_x(i), ',', test_y(i), ',', test_z(i), ')'

    call tdstress_hs(test_x(i), test_y(i), test_z(i), &
                     p1, p2, p3, ss, ds, ts, mu, lambda, &
                     stress, strain)

    print *, 'Result: Exx = ', strain(1)

    if (strain(1) /= strain(1)) then
      print *, '*** FAIL: NaN ***'
      n_fail = n_fail + 1
    else
      print *, '*** PASS: Finite value ***'
      n_pass = n_pass + 1
    end if
    print *, ''
  end do

  print *, '============================================================'
  print *, 'SUMMARY'
  print *, '============================================================'
  print *, 'Passed:', n_pass, ' out of', N_POINTS
  print *, 'Failed:', n_fail, ' out of', N_POINTS
  print *, ''

  if (n_fail > 0) then
    print *, 'WARNING: Some points returned NaN!'
    print *, 'This suggests a problem with the implementation.'
  else
    print *, 'SUCCESS: All points returned finite values.'
  end if

end program test_multiple_points
