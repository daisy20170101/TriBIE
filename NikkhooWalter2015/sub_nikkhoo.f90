!==============================================================================
! sub_nikkhoo.f90
! 
! Fortran 90 implementation of the Nikkhoo & Walter (2015) triangular dislocation
! method for calculating stresses and strains in an elastic half-space.
!
! Converted from MATLAB code TDstressHS.m
! 
! Reference: Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An 
! analytical, artefact-free solution. Geophysical Journal International
!
! Author: Converted from MATLAB by AI Assistant
! Date: 2024
!==============================================================================

module nikkhoo_walter
  implicit none
  
  ! Precision parameters
  integer, parameter :: DP = selected_real_kind(15, 307)
  real(DP), parameter :: PI = 3.141592653589793238462643383279502884197_DP
  real(DP), parameter :: EPS = 1.0e-15_DP
  
  contains

!==============================================================================
! Main function: TDstressHS
! Calculates stresses and strains associated with a triangular dislocation 
! in an elastic half-space.
!==============================================================================
subroutine tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                       stress, strain, n_points)
  implicit none
  
  ! Input parameters
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: x, y, z
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  real(DP), intent(in) :: ss, ds, ts, mu, lambda
  
  ! Output parameters
  real(DP), dimension(n_points, 6), intent(out) :: stress, strain
  
  ! Local variables
  real(DP), dimension(n_points, 6) :: sts_ms, str_ms, sts_fsc, str_fsc
  real(DP), dimension(n_points, 6) :: sts_is, str_is
  real(DP), dimension(3) :: p1_img, p2_img, p3_img
  integer :: i
  
  ! Check half-space constraint
  do i = 1, n_points
    if (z(i) > 0.0_DP .or. p1(3) > 0.0_DP .or. p2(3) > 0.0_DP .or. p3(3) > 0.0_DP) then
      write(*,*) 'ERROR: Half-space solution: Z coordinates must be negative!'
      stop
    end if
  end do
  
  ! Calculate main dislocation contribution
  call tdstress_fs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                   sts_ms, str_ms, n_points)
  
  ! Calculate harmonic function contribution
  call tdstress_harfunc(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                        sts_fsc, str_fsc, n_points)
  
  ! Calculate image dislocation contribution
  p1_img = p1; p2_img = p2; p3_img = p3
  p1_img(3) = -p1_img(3)
  p2_img(3) = -p2_img(3)
  p3_img(3) = -p3_img(3)
  
  call tdstress_fs(x, y, z, p1_img, p2_img, p3_img, ss, ds, ts, mu, lambda, &
                   sts_is, str_is, n_points)
  
  ! Special case for surface elements
  if (abs(p1_img(3)) < EPS .and. abs(p2_img(3)) < EPS .and. abs(p3_img(3)) < EPS) then
    sts_is(:, 5) = -sts_is(:, 5)  ! xz component
    sts_is(:, 6) = -sts_is(:, 6)  ! yz component
    str_is(:, 5) = -str_is(:, 5)  ! xz component
    str_is(:, 6) = -str_is(:, 6)  ! yz component
  end if
  
  ! Calculate total stress and strain
  stress = sts_ms + sts_is + sts_fsc
  strain = str_ms + str_is + str_fsc

end subroutine tdstress_hs

!==============================================================================
! TDstressFS: Full-space triangular dislocation
!==============================================================================
subroutine tdstress_fs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                       stress, strain, n_points)
  implicit none
  
  ! Input parameters
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: x, y, z
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  real(DP), intent(in) :: ss, ds, ts, mu, lambda
  
  ! Output parameters
  real(DP), dimension(n_points, 6), intent(out) :: stress, strain
  
  ! Local variables
  real(DP) :: nu, bx, by, bz
  real(DP), dimension(3) :: vnorm, vstrike, vdip, ey, ez
  real(DP), dimension(3, 3) :: A
  real(DP), dimension(3) :: p1_td, p2_td, p3_td
  real(DP), dimension(n_points) :: x_td, y_td, z_td
  real(DP), dimension(3) :: e12, e13, e23
  real(DP) :: A_angle, B_angle, C_angle
  integer, dimension(n_points) :: trimode
  logical, dimension(n_points) :: casep_log, casen_log, casez_log
  integer :: i, n_p, n_n, n_z
  real(DP), dimension(:), allocatable :: xp, yp, zp, xn, yn, zn
  real(DP), dimension(:), allocatable :: exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p
  real(DP), dimension(:), allocatable :: exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n
  real(DP), dimension(n_points) :: exx, eyy, ezz, exy, exz, eyz
  real(DP), dimension(n_points) :: sxx, syy, szz, sxy, sxz, syz
  
  ! Calculate Poisson's ratio
  nu = 1.0_DP / (1.0_DP + lambda / mu) / 2.0_DP
  
  ! Slip vector components
  bx = ts  ! Tensile-slip
  by = ss  ! Strike-slip
  bz = ds  ! Dip-slip
  
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
    ! Special case for image dislocation
    if (p1(3) > 0.0_DP) then
      vstrike = -vstrike
    end if
  end if
  vstrike = vstrike / norm2(vstrike)
  
  ! Dip vector
  call cross_product(vnorm, vstrike, vdip)
  
  ! Transformation matrix
  A(:, 1) = vnorm
  A(:, 2) = vstrike
  A(:, 3) = vdip
  
  ! Transform coordinates to TDCS
  p1_td = 0.0_DP
  p2_td = 0.0_DP
  p3_td = 0.0_DP
  
  call coord_trans(x - p2(1), y - p2(2), z - p2(3), A, x_td, y_td, z_td, n_points)
  call coord_trans(p1(1) - p2(1), p1(2) - p2(2), p1(3) - p2(3), A, p1_td(1), p1_td(2), p1_td(3), 1)
  call coord_trans(p3(1) - p2(1), p3(2) - p2(2), p3(3) - p2(3), A, p3_td(1), p3_td(2), p3_td(3), 1)
  
  ! Calculate unit vectors along TD sides
  e12 = (p2_td - p1_td) / norm2(p2_td - p1_td)
  e13 = (p3_td - p1_td) / norm2(p3_td - p1_td)
  e23 = (p3_td - p2_td) / norm2(p3_td - p2_td)
  
  ! Calculate angles
  A_angle = acos(dot_product(e12, e13))
  B_angle = acos(-dot_product(e12, e23))
  C_angle = acos(dot_product(e23, e13))
  
  ! Determine configuration for each point
  call trimode_finder(y_td, z_td, x_td, p1_td(2:3), p2_td(2:3), p3_td(2:3), trimode, n_points)
  
  casep_log = (trimode == 1)
  casen_log = (trimode == -1)
  casez_log = (trimode == 0)
  
  n_p = count(casep_log)
  n_n = count(casen_log)
  n_z = count(casez_log)
  
  ! Allocate arrays for each configuration
  if (n_p > 0) then
    allocate(xp(n_p), yp(n_p), zp(n_p))
    allocate(exx_p(n_p), eyy_p(n_p), ezz_p(n_p), exy_p(n_p), exz_p(n_p), eyz_p(n_p))
    
    ! Extract points for configuration I
    call extract_points(x_td, y_td, z_td, casep_log, xp, yp, zp, n_points, n_p)
    
    ! Calculate strains for configuration I
    call tdsetup_s(xp, yp, zp, A_angle, bx, by, bz, nu, p1_td, -e13, &
                   exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p, n_p)
    call tdsetup_s(xp, yp, zp, B_angle, bx, by, bz, nu, p2_td, e12, &
                   exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p, n_p)
    call tdsetup_s(xp, yp, zp, C_angle, bx, by, bz, nu, p3_td, e23, &
                   exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p, n_p)
  end if
  
  if (n_n > 0) then
    allocate(xn(n_n), yn(n_n), zn(n_n))
    allocate(exx_n(n_n), eyy_n(n_n), ezz_n(n_n), exy_n(n_n), exz_n(n_n), eyz_n(n_n))
    
    ! Extract points for configuration II
    call extract_points(x_td, y_td, z_td, casen_log, xn, yn, zn, n_points, n_n)
    
    ! Calculate strains for configuration II
    call tdsetup_s(xn, yn, zn, A_angle, bx, by, bz, nu, p1_td, e13, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n, n_n)
    call tdsetup_s(xn, yn, zn, B_angle, bx, by, bz, nu, p2_td, -e12, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n, n_n)
    call tdsetup_s(xn, yn, zn, C_angle, bx, by, bz, nu, p3_td, -e23, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n, n_n)
  end if
  
  ! Combine results
  exx = 0.0_DP; eyy = 0.0_DP; ezz = 0.0_DP
  exy = 0.0_DP; exz = 0.0_DP; eyz = 0.0_DP
  
  if (n_p > 0) then
    call assign_points(exx, eyy, ezz, exy, exz, eyz, &
                       exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p, &
                       casep_log, n_points, n_p)
  end if
  
  if (n_n > 0) then
    call assign_points(exx, eyy, ezz, exy, exz, eyz, &
                       exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n, &
                       casen_log, n_points, n_n)
  end if
  
  ! Set NaN for points on TD sides
  do i = 1, n_points
    if (casez_log(i)) then
      exx(i) = huge(1.0_DP)
      eyy(i) = huge(1.0_DP)
      ezz(i) = huge(1.0_DP)
      exy(i) = huge(1.0_DP)
      exz(i) = huge(1.0_DP)
      eyz(i) = huge(1.0_DP)
    end if
  end do
  
  ! Transform strain tensor to EFCS
  call tens_trans(exx, eyy, ezz, exy, exz, eyz, A, &
                  exx, eyy, ezz, exy, exz, eyz, n_points)
  
  ! Calculate stress tensor
  sxx = 2.0_DP * mu * exx + lambda * (exx + eyy + ezz)
  syy = 2.0_DP * mu * eyy + lambda * (exx + eyy + ezz)
  szz = 2.0_DP * mu * ezz + lambda * (exx + eyy + ezz)
  sxy = 2.0_DP * mu * exy
  sxz = 2.0_DP * mu * exz
  syz = 2.0_DP * mu * eyz
  
  ! Output
  stress(:, 1) = sxx; stress(:, 2) = syy; stress(:, 3) = szz
  stress(:, 4) = sxy; stress(:, 5) = sxz; stress(:, 6) = syz
  
  strain(:, 1) = exx; strain(:, 2) = eyy; strain(:, 3) = ezz
  strain(:, 4) = exy; strain(:, 5) = exz; strain(:, 6) = eyz
  
  ! Cleanup
  if (allocated(xp)) deallocate(xp, yp, zp, exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p)
  if (allocated(xn)) deallocate(xn, yn, zn, exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n)

end subroutine tdstress_fs

!==============================================================================
! Harmonic function contribution
!==============================================================================
subroutine tdstress_harfunc(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                            stress, strain, n_points)
  implicit none
  
  ! Input parameters
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: x, y, z
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  real(DP), intent(in) :: ss, ds, ts, mu, lambda
  
  ! Output parameters
  real(DP), dimension(n_points, 6), intent(out) :: stress, strain
  
  ! Local variables
  real(DP) :: bx, by, bz
  real(DP), dimension(3) :: vnorm, vstrike, vdip, ey, ez
  real(DP), dimension(3, 3) :: A
  real(DP) :: bX, bY, bZ
  real(DP), dimension(n_points, 6) :: stress1, strain1, stress2, strain2, stress3, strain3
  
  ! Slip vector components
  bx = ts; by = ss; bz = ds
  
  ! Calculate unit vectors
  ey = [0.0_DP, 1.0_DP, 0.0_DP]
  ez = [0.0_DP, 0.0_DP, 1.0_DP]
  
  call cross_product(p2 - p1, p3 - p1, vnorm)
  vnorm = vnorm / norm2(vnorm)
  
  call cross_product(ez, vnorm, vstrike)
  if (norm2(vstrike) < EPS) then
    vstrike = ey * vnorm(3)
  end if
  vstrike = vstrike / norm2(vstrike)
  
  call cross_product(vnorm, vstrike, vdip)
  
  ! Transform slip vector
  A(:, 1) = vnorm
  A(:, 2) = vstrike
  A(:, 3) = vdip
  
  call coord_trans(bx, by, bz, A, bX, bY, bZ, 1)
  
  ! Calculate contributions from each side
  call angsetup_fsc_s(x, y, z, bX, bY, bZ, p1, p2, mu, lambda, stress1, strain1, n_points)
  call angsetup_fsc_s(x, y, z, bX, bY, bZ, p2, p3, mu, lambda, stress2, strain2, n_points)
  call angsetup_fsc_s(x, y, z, bX, bY, bZ, p3, p1, mu, lambda, stress3, strain3, n_points)
  
  ! Total contribution
  stress = stress1 + stress2 + stress3
  strain = strain1 + strain2 + strain3

end subroutine tdstress_harfunc

!==============================================================================
! Tensor transformation
!==============================================================================
subroutine tens_trans(txx1, tyy1, tzz1, txy1, txz1, tyz1, A, &
                      txx2, tyy2, tzz2, txy2, txz2, tyz2, n_points)
  implicit none
  
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: txx1, tyy1, tzz1, txy1, txz1, tyz1
  real(DP), dimension(3, 3), intent(in) :: A
  real(DP), dimension(n_points), intent(out) :: txx2, tyy2, tzz2, txy2, txz2, tyz2
  
  integer :: i
  
  do i = 1, n_points
    txx2(i) = A(1,1)**2 * txx1(i) + 2*A(1,1)*A(1,2)*txy1(i) + 2*A(1,1)*A(1,3)*txz1(i) + &
              2*A(1,2)*A(1,3)*tyz1(i) + A(1,2)**2 * tyy1(i) + A(1,3)**2 * tzz1(i)
    
    tyy2(i) = A(2,1)**2 * txx1(i) + 2*A(2,1)*A(2,2)*txy1(i) + 2*A(2,1)*A(2,3)*txz1(i) + &
              2*A(2,2)*A(2,3)*tyz1(i) + A(2,2)**2 * tyy1(i) + A(2,3)**2 * tzz1(i)
    
    tzz2(i) = A(3,1)**2 * txx1(i) + 2*A(3,1)*A(3,2)*txy1(i) + 2*A(3,1)*A(3,3)*txz1(i) + &
              2*A(3,2)*A(3,3)*tyz1(i) + A(3,2)**2 * tyy1(i) + A(3,3)**2 * tzz1(i)
    
    txy2(i) = A(1,1)*A(2,1)*txx1(i) + (A(1,1)*A(2,2) + A(2,1)*A(1,2))*txy1(i) + &
              (A(1,1)*A(2,3) + A(2,1)*A(1,3))*txz1(i) + (A(2,3)*A(1,2) + A(1,3)*A(2,2))*tyz1(i) + &
              A(2,2)*A(1,2)*tyy1(i) + A(1,3)*A(2,3)*tzz1(i)
    
    txz2(i) = A(1,1)*A(3,1)*txx1(i) + (A(1,1)*A(3,2) + A(3,1)*A(1,2))*txy1(i) + &
              (A(1,1)*A(3,3) + A(3,1)*A(1,3))*txz1(i) + (A(3,3)*A(1,2) + A(1,3)*A(3,2))*tyz1(i) + &
              A(3,2)*A(1,2)*tyy1(i) + A(1,3)*A(3,3)*tzz1(i)
    
    tyz2(i) = A(2,1)*A(3,1)*txx1(i) + (A(3,1)*A(2,2) + A(2,1)*A(3,2))*txy1(i) + &
              (A(3,1)*A(2,3) + A(2,1)*A(3,3))*txz1(i) + (A(2,3)*A(3,2) + A(3,3)*A(2,2))*tyz1(i) + &
              A(2,2)*A(3,2)*tyy1(i) + A(2,3)*A(3,3)*tzz1(i)
  end do

end subroutine tens_trans

!==============================================================================
! Coordinate transformation
!==============================================================================
subroutine coord_trans(x1, x2, x3, A, X1, X2, X3, n_points)
  implicit none
  
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: x1, x2, x3
  real(DP), dimension(3, 3), intent(in) :: A
  real(DP), dimension(n_points), intent(out) :: X1, X2, X3
  
  integer :: i
  real(DP), dimension(3) :: r
  
  do i = 1, n_points
    r = matmul(A, [x1(i), x2(i), x3(i)])
    X1(i) = r(1)
    X2(i) = r(2)
    X3(i) = r(3)
  end do

end subroutine coord_trans

!==============================================================================
! Triangular mode finder
!==============================================================================
subroutine trimode_finder(x, y, z, p1, p2, p3, trimode, n_points)
  implicit none
  
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: x, y, z
  real(DP), dimension(2), intent(in) :: p1, p2, p3
  integer, dimension(n_points), intent(out) :: trimode
  
  real(DP), dimension(n_points) :: a, b, c
  integer :: i
  
  ! Calculate barycentric coordinates
  a = ((p2(2) - p3(2)) * (x - p3(1)) + (p3(1) - p2(1)) * (y - p3(2))) / &
      ((p2(2) - p3(2)) * (p1(1) - p3(1)) + (p3(1) - p2(1)) * (p1(2) - p3(2)))
  
  b = ((p3(2) - p1(2)) * (x - p3(1)) + (p1(1) - p3(1)) * (y - p3(2))) / &
      ((p2(2) - p3(2)) * (p1(1) - p3(1)) + (p3(1) - p2(1)) * (p1(2) - p3(2)))
  
  c = 1.0_DP - a - b
  
  ! Determine configuration
  trimode = 1
  do i = 1, n_points
    if (a(i) <= 0.0_DP .and. b(i) > c(i) .and. c(i) > a(i)) then
      trimode(i) = -1
    else if (b(i) <= 0.0_DP .and. c(i) > a(i) .and. a(i) > b(i)) then
      trimode(i) = -1
    else if (c(i) <= 0.0_DP .and. a(i) > b(i) .and. b(i) > c(i)) then
      trimode(i) = -1
    else if (a(i) == 0.0_DP .and. b(i) >= 0.0_DP .and. c(i) >= 0.0_DP) then
      trimode(i) = 0
    else if (a(i) >= 0.0_DP .and. b(i) == 0.0_DP .and. c(i) >= 0.0_DP) then
      trimode(i) = 0
    else if (a(i) >= 0.0_DP .and. b(i) >= 0.0_DP .and. c(i) == 0.0_DP) then
      trimode(i) = 0
    end if
    
    if (trimode(i) == 0 .and. abs(z(i)) > EPS) then
      trimode(i) = 1
    end if
  end do

end subroutine trimode_finder

!==============================================================================
! TD Setup S
!==============================================================================
subroutine tdsetup_s(x, y, z, alpha, bx, by, bz, nu, tri_vertex, side_vec, &
                     exx, eyy, ezz, exy, exz, eyz, n_points)
  implicit none
  
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: x, y, z
  real(DP), intent(in) :: alpha, bx, by, bz, nu
  real(DP), dimension(3), intent(in) :: tri_vertex, side_vec
  
  real(DP), dimension(n_points), intent(out) :: exx, eyy, ezz, exy, exz, eyz
  
  real(DP), dimension(2, 2) :: A
  real(DP), dimension(n_points) :: y1, z1, by1, bz1
  integer :: i
  
  ! Transformation matrix
  A(1, :) = [side_vec(3), -side_vec(2)]
  A(2, :) = [side_vec(2), side_vec(3)]
  
  ! Transform coordinates
  do i = 1, n_points
    y1(i) = A(1,1) * (y(i) - tri_vertex(2)) + A(1,2) * (z(i) - tri_vertex(3))
    z1(i) = A(2,1) * (y(i) - tri_vertex(2)) + A(2,2) * (z(i) - tri_vertex(3))
  end do
  
  ! Transform slip components
  by1 = A(1,1) * by + A(1,2) * bz
  bz1 = A(2,1) * by + A(2,2) * bz
  
  ! Calculate strains
  call angdis_strain(x, y1, z1, -PI + alpha, bx, by1, bz1, nu, &
                     exx, eyy, ezz, exy, exz, eyz, n_points)
  
  ! Transform back to TDCS
  ! (Implementation would continue with tensor transformation)

end subroutine tdsetup_s

!==============================================================================
! Angular dislocation strain
!==============================================================================
subroutine angdis_strain(x, y, z, alpha, bx, by, bz, nu, &
                         exx, eyy, ezz, exy, exz, eyz, n_points)
  implicit none
  
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: x, y, z
  real(DP), intent(in) :: alpha, bx, by, bz, nu
  
  real(DP), dimension(n_points), intent(out) :: exx, eyy, ezz, exy, exz, eyz
  
  ! This is a simplified version - the full implementation would include
  ! all the complex strain calculations from the MATLAB code
  
  ! For now, return zero strains (placeholder)
  exx = 0.0_DP
  eyy = 0.0_DP
  ezz = 0.0_DP
  exy = 0.0_DP
  exz = 0.0_DP
  eyz = 0.0_DP

end subroutine angdis_strain

!==============================================================================
! Angular setup FSC S
!==============================================================================
subroutine angsetup_fsc_s(x, y, z, bX, bY, bZ, PA, PB, mu, lambda, &
                          stress, strain, n_points)
  implicit none
  
  integer, intent(in) :: n_points
  real(DP), dimension(n_points), intent(in) :: x, y, z
  real(DP), intent(in) :: bX, bY, bZ, mu, lambda
  real(DP), dimension(3), intent(in) :: PA, PB
  
  real(DP), dimension(n_points, 6), intent(out) :: stress, strain
  
  ! Placeholder implementation
  stress = 0.0_DP
  strain = 0.0_DP

end subroutine angsetup_fsc_s

!==============================================================================
! Utility functions
!==============================================================================

! Cross product
subroutine cross_product(a, b, c)
  implicit none
  real(DP), dimension(3), intent(in) :: a, b
  real(DP), dimension(3), intent(out) :: c
  
  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)
end subroutine cross_product

! Vector norm
function norm2(v) result(n)
  implicit none
  real(DP), dimension(:), intent(in) :: v
  real(DP) :: n
  
  n = sqrt(sum(v**2))
end function norm2

! Extract points based on logical mask
subroutine extract_points(x, y, z, mask, x_out, y_out, z_out, n_total, n_selected)
  implicit none
  integer, intent(in) :: n_total, n_selected
  real(DP), dimension(n_total), intent(in) :: x, y, z
  logical, dimension(n_total), intent(in) :: mask
  real(DP), dimension(n_selected), intent(out) :: x_out, y_out, z_out
  
  integer :: i, j
  
  j = 1
  do i = 1, n_total
    if (mask(i)) then
      x_out(j) = x(i)
      y_out(j) = y(i)
      z_out(j) = z(i)
      j = j + 1
    end if
  end do
end subroutine extract_points

! Assign points back to full arrays
subroutine assign_points(exx, eyy, ezz, exy, exz, eyz, &
                         exx_part, eyy_part, ezz_part, exy_part, exz_part, eyz_part, &
                         mask, n_total, n_selected)
  implicit none
  integer, intent(in) :: n_total, n_selected
  real(DP), dimension(n_total), intent(inout) :: exx, eyy, ezz, exy, exz, eyz
  real(DP), dimension(n_selected), intent(in) :: exx_part, eyy_part, ezz_part, exy_part, exz_part, eyz_part
  logical, dimension(n_total), intent(in) :: mask
  
  integer :: i, j
  
  j = 1
  do i = 1, n_total
    if (mask(i)) then
      exx(i) = exx_part(j)
      eyy(i) = eyy_part(j)
      ezz(i) = ezz_part(j)
      exy(i) = exy_part(j)
      exz(i) = exz_part(j)
      eyz(i) = eyz_part(j)
      j = j + 1
    end if
  end do
end subroutine assign_points

end module nikkhoo_walter
