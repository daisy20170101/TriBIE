program test_casep
  use nikkhoo_walter
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Triangle vertices
  real(DP), dimension(3) :: p1, p2, p3

  ! Test points
  integer, parameter :: n_points = 15
  real(DP), dimension(n_points) :: x, y, z
  real(DP), dimension(n_points, 6) :: stress, strain

  ! Slip and elastic parameters
  real(DP) :: ss, ds, ts, mu, lambda

  integer :: i

  ! Set up triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Set up test points (same as test_casez)
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
  print *, '============================================================================='
  print *, 'Testing Nikkhoo-Walter Triangular Dislocation Method'
  print *, '============================================================================='
  print *, 'Triangle vertices:'
  print *, '  p1 = ', p1
  print *, '  p2 = ', p2
  print *, '  p3 = ', p3
  print *, ''
  print *, 'Slip components: ss=', ss, ' ds=', ds, ' ts=', ts
  print *, 'Elastic params: mu=', mu, ' lambda=', lambda
  print *, ''

  ! Calculate stress and strain for each test point
  do i = 1, n_points
    call tdstress_hs(x(i), y(i), z(i), p1, p2, p3, ss, ds, ts, mu, lambda, &
                     stress(i,:), strain(i,:))
  end do

  ! Output stress components table
  write(*,*) 'Stress components:'
  write(*,*) 'Point    X        Y        Z        Sxx       Syy       Szz       Sxy       Sxz       Syz'
  write(*,*) '-------------------------------------------------------------------------------------------'

  do i = 1, n_points
    write(*,'(I3,3F9.3,6F10.3)') i, x(i), y(i), z(i), stress(i,1), stress(i,2), stress(i,3), &
                                 stress(i,4), stress(i,5), stress(i,6)
  end do

  write(*,*) ''
  write(*,*) 'Strain components:'
  write(*,*) 'Point    Exx      Eyy      Ezz      Exy      Exz      Eyz'
  write(*,*) '--------------------------------------------------------'

  do i = 1, n_points
    write(*,'(I3,6F9.6)') i, strain(i,1), strain(i,2), strain(i,3), &
                         strain(i,4), strain(i,5), strain(i,6)
  end do

  print *, ''
  print *, '============================================================================='

end program test_casep
