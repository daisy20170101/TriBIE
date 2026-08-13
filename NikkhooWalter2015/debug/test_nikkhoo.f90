!==============================================================================
! test_nikkhoo.f90
! 
! Test program for the Nikkhoo & Walter (2015) triangular dislocation method
!==============================================================================

program test_nikkhoo
  use nikkhoo_walter
  implicit none
  
  ! Test parameters
  integer, parameter :: n_points = 7
  integer :: i
  real(DP), dimension(n_points) :: x, y, z
  real(DP), dimension(3) :: p1, p2, p3
  real(DP) :: ss, ds, ts, mu, lambda
  real(DP), dimension(n_points, 6) :: stress, strain
  
  ! Initialize test data
  x = [-1.0_DP/3.0_DP, -1.0_DP/3.0_DP, -1.0_DP/3.0_DP, 7.0_DP, -7.0_DP, -1.0_DP, -1.0_DP]
  y = [-1.0_DP/3.0_DP, -1.0_DP/3.0_DP, -1.0_DP/3.0_DP, -1.0_DP, -1.0_DP, -3.0_DP, 3.0_DP]
  z = [-3.0_DP, -14.0_DP/3.0_DP, -6.0_DP, -5.0_DP, -5.0_DP, -6.0_DP, -3.0_DP]
  
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]
  
  ss = 1.0_DP  ! Strike-slip
  ds = -1.0_DP  ! Dip-slip
  ts = 2.0_DP  ! Tensile-slip
  
  mu = 3.0e10_DP      ! Shear modulus
  lambda = 3.0e10_DP  ! Lame's first parameter
  
  ! Calculate stresses and strains (loop over calculation points)
  do i = 1, n_points
    call tdstress_hs(x(i), y(i), z(i), p1, p2, p3, ss, ds, ts, mu, lambda, &
                     stress(i, :), strain(i, :))
  end do
  
  ! Output results
  write(*,*) 'Nikkhoo & Walter (2015) Triangular Dislocation Test'
  write(*,*) '=================================================='
  write(*,*) 'Number of calculation points:', n_points
  write(*,*) 'Triangular dislocation vertices:'
  write(*,*) '  P1 =', p1
  write(*,*) '  P2 =', p2
  write(*,*) '  P3 =', p3
  write(*,*) 'Slip components: SS =', ss, ', DS =', ds, ', TS =', ts
  write(*,*) 'Elastic parameters: mu =', mu, ', lambda =', lambda
  write(*,*) ''
  write(*,*) 'Results:'
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

end program test_nikkhoo
