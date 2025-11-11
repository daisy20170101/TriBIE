program test_casez
  use nikkhoo_walter
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Triangle vertices
  real(DP), dimension(3) :: p1, p2, p3

  ! Test points
  real(DP), dimension(15) :: x, y, z
  real(DP), dimension(15) :: exx_results

  ! Slip and elastic parameters
  real(DP) :: ss, ds, ts, mu, lambda
  real(DP), dimension(6) :: stress, strain

  integer :: i

  ! Set up triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Set up test points
  x = [-1.0_DP/3.0_DP, -1.0_DP/3.0_DP, -1.0_DP/3.0_DP, 7.0_DP, -7.0_DP, -1.0_DP, -1.0_DP, &
       3.0_DP, -3.0_DP, -1.0_DP, -1.0_DP, 1.0_DP, -1.0_DP, -1.0_DP, 1.0_DP]
  y = [-1.0_DP/3.0_DP, -1.0_DP/3.0_DP, -1.0_DP/3.0_DP, -1.0_DP, -1.0_DP, -3.0_DP, 3.0_DP, &
       -3.0_DP, 3.0_DP, -1.0_DP, 1.0_DP, -1.0_DP, -1.0_DP, 1.0_DP, -1.0_DP]
  z = [-3.0_DP, -14.0_DP/3.0_DP, -6.0_DP, -5.0_DP, -5.0_DP, -6.0_DP, -3.0_DP, &
       -6.0_DP, -3.0_DP, -1.0_DP, -1.0_DP, -1.0_DP, -8.0_DP, -8.0_DP, -8.0_DP]

  ! Set slip components
  ss = 1.0_DP   ! Strike-slip
  ds = -1.0_DP  ! Dip-slip
  ts = 2.0_DP   ! Tensile-slip

  ! Set elastic parameters
  mu = 3.0e10_DP
  lambda = 3.0e10_DP

  ! Print header
  print *, '=============================================='
  print *, 'Testing casez_log implementation'
  print *, '=============================================='
  print *, 'Triangle vertices:'
  print *, '  p1 = ', p1
  print *, '  p2 = ', p2
  print *, '  p3 = ', p3
  print *, ''
  print *, 'Slip components: ss=', ss, ' ds=', ds, ' ts=', ts
  print *, 'Elastic params: mu=', mu, ' lambda=', lambda
  print *, ''
  print *, '=============================================='
  print *, 'Point#      x          y          z         e_xx'
  print *, '=============================================='

  ! Calculate strain for each test point
  do i = 1, 15
    call tdstress_hs(x(i), y(i), z(i), p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
    exx_results(i) = strain(1)

    ! Print results
    if (ieee_is_nan(strain(1))) then
      write(*, '(I4, 3F11.4, A15)') i, x(i), y(i), z(i), '        NaN'
    else
      write(*, '(I4, 3F11.4, ES15.6)') i, x(i), y(i), z(i), strain(1)
    end if
  end do

  print *, '=============================================='

end program test_casez
