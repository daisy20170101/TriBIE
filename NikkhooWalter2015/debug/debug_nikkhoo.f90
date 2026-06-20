!==============================================================================
! debug_nikkhoo.f90
! 
! Debug program for the Nikkhoo & Walter (2015) triangular dislocation method
! This version includes detailed debug output to identify inconsistencies
!==============================================================================

program debug_nikkhoo
  use nikkhoo_walter
  implicit none
  
  ! Test parameters
  integer, parameter :: n_points = 1  ! Single point for detailed debugging
  integer :: i
  real(DP), dimension(n_points) :: x, y, z
  real(DP), dimension(3) :: p1, p2, p3
  real(DP) :: ss, ds, ts, mu, lambda
  real(DP), dimension(n_points, 6) :: stress, strain
  
  ! Debug variables
  real(DP) :: nu, bx, by, bz
  real(DP), dimension(3) :: vnorm, vstrike, vdip, ey, ez
  real(DP), dimension(3, 3) :: A
  real(DP), dimension(3) :: p1_td, p2_td, p3_td
  real(DP) :: x_td, y_td, z_td
  real(DP), dimension(3) :: e12, e13, e23
  real(DP) :: A_angle, B_angle, C_angle
  integer :: trimode
  logical :: casep_log, casen_log, casez_log
  
  ! Initialize test data (single point for debugging)
  x = [-1.0_DP/3.0_DP]
  y = [-1.0_DP/3.0_DP]
  z = [-14.0_DP/3.0_DP]
  
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]
  
  ss = 1.0_DP  ! Strike-slip
  ds = -1.0_DP  ! Dip-slip
  ts = 2.0_DP  ! Tensile-slip
  
  mu = 3.0e10_DP      ! Shear modulus
  lambda = 3.0e10_DP  ! Lame's first parameter
  
  ! Calculate Poisson's ratio
  nu = 1.0_DP / (1.0_DP + lambda / mu) / 2.0_DP
  
  ! Slip vector components
  bx = ts  ! Tensile-slip
  by = ss  ! Strike-slip
  bz = ds  ! Dip-slip
  
  write(*,*) '=== DEBUG: Nikkhoo & Walter (2015) Triangular Dislocation ==='
  write(*,*) 'Input Parameters:'
  write(*,*) '  Calculation point: (', x(1), ',', y(1), ',', z(1), ')'
  write(*,*) '  Triangle vertices:'
  write(*,*) '    P1 = (', p1(1), ',', p1(2), ',', p1(3), ')'
  write(*,*) '    P2 = (', p2(1), ',', p2(2), ',', p2(3), ')'
  write(*,*) '    P3 = (', p3(1), ',', p3(2), ',', p3(3), ')'
  write(*,*) '  Slip components: SS =', ss, ', DS =', ds, ', TS =', ts
  write(*,*) '  Elastic parameters: mu =', mu, ', lambda =', lambda
  write(*,*) '  Poisson ratio: nu =', nu
  write(*,*) '  Slip vector: (', bx, ',', by, ',', bz, ')'
  write(*,*) ''
  
  ! Calculate unit vectors
  ey = [0.0_DP, 1.0_DP, 0.0_DP]
  ez = [0.0_DP, 0.0_DP, 1.0_DP]
  
  ! Normal vector
  call cross_product(p2 - p1, p3 - p1, vnorm)
  vnorm = vnorm / norm2(vnorm)
  
  ! Strike vector
  call cross_product(ez, vnorm, vstrike)
  if (norm2(vstrike) < EPS) then
    vstrike = ey * vnorm(3)
    if (p1(3) > 0.0_DP) then
      vstrike = -vstrike
    end if
  end if
  vstrike = vstrike / norm2(vstrike)
  
  ! Dip vector
  call cross_product(vnorm, vstrike, vdip)
  
  write(*,*) '=== Coordinate System Vectors ==='
  write(*,*) '  vnorm = (', vnorm(1), ',', vnorm(2), ',', vnorm(3), ')'
  write(*,*) '  vstrike = (', vstrike(1), ',', vstrike(2), ',', vstrike(3), ')'
  write(*,*) '  vdip = (', vdip(1), ',', vdip(2), ',', vdip(3), ')'
  write(*,*) ''
  
  ! Transformation matrix (columns are unit vectors, matching sub_nikkhoo.f90)
  ! coord_trans uses transpose(A), so we need columns here
  A(:, 1) = vnorm
  A(:, 2) = vstrike
  A(:, 3) = vdip
  
  write(*,*) '=== Transformation Matrix A (EFCS to TDCS) ==='
  write(*,*) '  A(1,:) = (', A(1,1), ',', A(1,2), ',', A(1,3), ')'
  write(*,*) '  A(2,:) = (', A(2,1), ',', A(2,2), ',', A(2,3), ')'
  write(*,*) '  A(3,:) = (', A(3,1), ',', A(3,2), ',', A(3,3), ')'
  write(*,*) ''
  
  ! Transform coordinates to TDCS
  p1_td = 0.0_DP
  p2_td = 0.0_DP
  p3_td = 0.0_DP
  
  call coord_trans(x(1) - p2(1), y(1) - p2(2), z(1) - p2(3), A, x_td, y_td, z_td)
  call coord_trans(p1(1) - p2(1), p1(2) - p2(2), p1(3) - p2(3), A, p1_td(1), p1_td(2), p1_td(3))
  call coord_trans(p3(1) - p2(1), p3(2) - p2(2), p3(3) - p2(3), A, p3_td(1), p3_td(2), p3_td(3))
  
  write(*,*) '=== TDCS Coordinates ==='
  write(*,*) '  Calculation point: (', x_td, ',', y_td, ',', z_td, ')'
  write(*,*) '  Triangle vertices:'
  write(*,*) '    p1_td = (', p1_td(1), ',', p1_td(2), ',', p1_td(3), ')'
  write(*,*) '    p2_td = (', p2_td(1), ',', p2_td(2), ',', p2_td(3), ')'
  write(*,*) '    p3_td = (', p3_td(1), ',', p3_td(2), ',', p3_td(3), ')'
  write(*,*) ''
  
  ! Calculate unit vectors along TD sides
  e12 = (p2_td - p1_td) / norm2(p2_td - p1_td)
  e13 = (p3_td - p1_td) / norm2(p3_td - p1_td)
  e23 = (p3_td - p2_td) / norm2(p3_td - p2_td)
  
  write(*,*) '=== Unit Vectors Along TD Sides ==='
  write(*,*) '  e12 = (', e12(1), ',', e12(2), ',', e12(3), ')'
  write(*,*) '  e13 = (', e13(1), ',', e13(2), ',', e13(3), ')'
  write(*,*) '  e23 = (', e23(1), ',', e23(2), ',', e23(3), ')'
  write(*,*) ''
  
  ! Calculate angles
  A_angle = acos(dot_product(e12, e13))
  B_angle = acos(-dot_product(e12, e23))
  C_angle = acos(dot_product(e23, e13))
  
  write(*,*) '=== Triangle Angles ==='
  write(*,*) '  A_angle =', A_angle, 'rad =', A_angle * 180.0_DP / PI, 'deg'
  write(*,*) '  B_angle =', B_angle, 'rad =', B_angle * 180.0_DP / PI, 'deg'
  write(*,*) '  C_angle =', C_angle, 'rad =', C_angle * 180.0_DP / PI, 'deg'
  write(*,*) ''
  
  ! Determine configuration
  call trimode_finder(y_td, z_td, x_td, p1_td, p2_td, p3_td, trimode)
  
  casep_log = (trimode == 1)
  casen_log = (trimode == -1)
  casez_log = (trimode == 0)
  
  write(*,*) '=== Configuration ==='
  write(*,*) '  trimode =', trimode
  write(*,*) '  casep_log =', casep_log
  write(*,*) '  casen_log =', casen_log
  write(*,*) '  casez_log =', casez_log
  write(*,*) ''
  
  ! Calculate stresses and strains with detailed debugging
  write(*,*) '=== Calling tdstress_hs ==='
  call tdstress_hs(x(1), y(1), z(1), p1, p2, p3, ss, ds, ts, mu, lambda, &
                   stress(1, :), strain(1, :))
  
  write(*,*) '=== Final Results ==='
  write(*,*) 'Stress tensor:'
  write(*,*) '  Sxx =', stress(1, 1)
  write(*,*) '  Syy =', stress(1, 2)
  write(*,*) '  Szz =', stress(1, 3)
  write(*,*) '  Sxy =', stress(1, 4)
  write(*,*) '  Sxz =', stress(1, 5)
  write(*,*) '  Syz =', stress(1, 6)
  write(*,*) ''
  write(*,*) 'Strain tensor:'
  write(*,*) '  Exx =', strain(1, 1)
  write(*,*) '  Eyy =', strain(1, 2)
  write(*,*) '  Ezz =', strain(1, 3)
  write(*,*) '  Exy =', strain(1, 4)
  write(*,*) '  Exz =', strain(1, 5)
  write(*,*) '  Eyz =', strain(1, 6)
  write(*,*) ''

end program debug_nikkhoo
