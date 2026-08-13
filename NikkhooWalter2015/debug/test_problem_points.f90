program test_problem_points
  use nikkhoo_walter
  implicit none

  ! Triangle vertices
  real(DP) :: p1(3), p2(3), p3(3)

  ! Test points
  integer, parameter :: N_POINTS = 4
  real(DP) :: test_x(N_POINTS), test_y(N_POINTS), test_z(N_POINTS)
  integer :: point_nums(N_POINTS)
  character(len=50) :: descriptions(N_POINTS)

  ! Slip components
  real(DP) :: ss, ds, ts

  ! Elastic parameters
  real(DP) :: mu, lambda

  ! Results
  real(DP) :: stress(6), strain(6)
  integer :: i

  ! Initialize triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Define the problematic points from original test_casep
  ! Point 4: x=7.0, y=-1.0, z=-5.0
  test_x(1) = 7.0_DP
  test_y(1) = -1.0_DP
  test_z(1) = -5.0_DP
  point_nums(1) = 4
  descriptions(1) = 'Point 4: Far to the side'

  ! Point 5: x=-7.0, y=-1.0, z=-5.0
  test_x(2) = -7.0_DP
  test_y(2) = -1.0_DP
  test_z(2) = -5.0_DP
  point_nums(2) = 5
  descriptions(2) = 'Point 5: Far to the other side'

  ! Point 12: x=1.0, y=-1.0, z=-1.0
  test_x(3) = 1.0_DP
  test_y(3) = -1.0_DP
  test_z(3) = -1.0_DP
  point_nums(3) = 12
  descriptions(3) = 'Point 12: At P2 x-y, different z'

  ! Point 15: x=1.0, y=-1.0, z=-8.0
  test_x(4) = 1.0_DP
  test_y(4) = -1.0_DP
  test_z(4) = -8.0_DP
  point_nums(4) = 15
  descriptions(4) = 'Point 15: At P2 x-y, different z'

  ! Slip components
  ss = 1.0_DP
  ds = -1.0_DP
  ts = 2.0_DP

  ! Elastic parameters
  mu = 3.0e10_DP
  lambda = 3.0e10_DP

  print *, '============================================================'
  print *, 'Testing Problematic Points 4, 5, 12, 15'
  print *, '============================================================'
  print *, 'Triangle vertices:'
  print *, 'P1 = ', p1
  print *, 'P2 = ', p2
  print *, 'P3 = ', p3
  print *, ''
  print *, 'Observations:'
  print *, '  P2 = (1.0, -1.0, -5.0)'
  print *, '  Points 12 and 15 have x=1.0, y=-1.0 (same as P2 in x-y)'
  print *, '  Points 4 and 5 have y=-1.0 (same as P1 and P2 y-coord)'
  print *, ''

  do i = 1, N_POINTS
    print *, '------------------------------------------------------------'
    print *, trim(descriptions(i))
    print *, 'Position: (', test_x(i), ',', test_y(i), ',', test_z(i), ')'
    print *, ''

    call tdstress_hs(test_x(i), test_y(i), test_z(i), &
                     p1, p2, p3, ss, ds, ts, mu, lambda, &
                     stress, strain)

    print *, 'Result: Exx = ', strain(1)

    if (strain(1) /= strain(1)) then
      print *, '*** FAIL: NaN ***'
    else
      print *, '*** PASS: Finite value ***'
    end if
    print *, ''
  end do

  print *, '============================================================'
  print *, 'Analysis:'
  print *, '  If Points 12 and 15 fail: Issue with points aligned with'
  print *, '    P2 in x-y but different z (extended line from P2)'
  print *, '  If Points 4 and 5 fail: Issue with y=-1.0 (edge P1-P2 line)'
  print *, '============================================================'

end program test_problem_points
