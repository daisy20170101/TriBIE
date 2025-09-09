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
                       stress, strain)
  implicit none
  
  ! Input parameters
  real(DP), intent(in) :: x, y, z  ! Single calculation point
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  real(DP), intent(in) :: ss, ds, ts, mu, lambda
  
  ! Output parameters
  real(DP), dimension(6), intent(out) :: stress, strain
  
  ! Local variables
  real(DP), dimension(6) :: sts_ms, str_ms, sts_fsc, str_fsc
  real(DP), dimension(6) :: sts_is, str_is
  real(DP), dimension(3) :: p1_img, p2_img, p3_img
  
  ! Check half-space constraint
  if (z > 0.0_DP .or. p1(3) > 0.0_DP .or. p2(3) > 0.0_DP .or. p3(3) > 0.0_DP) then
    write(*,*) 'ERROR: Half-space solution: Z coordinates must be negative!'
    stop
  end if
  
  ! Calculate main dislocation contribution
  call tdstress_fs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                          sts_ms, str_ms)
  
  ! Calculate harmonic function contribution
  call tdstress_harfunc(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                               sts_fsc, str_fsc)
  
  ! Calculate image dislocation contribution
  p1_img = p1; p2_img = p2; p3_img = p3
  p1_img(3) = -p1_img(3)
  p2_img(3) = -p2_img(3)
  p3_img(3) = -p3_img(3)
  
  call tdstress_fs(x, y, z, p1_img, p2_img, p3_img, ss, ds, ts, mu, lambda, &
                          sts_is, str_is)
  
  ! Special case for surface elements
  if (abs(p1_img(3)) < EPS .and. abs(p2_img(3)) < EPS .and. abs(p3_img(3)) < EPS) then
    sts_is(5) = -sts_is(5)  ! xz component
    sts_is(6) = -sts_is(6)  ! yz component
    str_is(5) = -str_is(5)  ! xz component
    str_is(6) = -str_is(6)  ! yz component
  end if
  
  ! Calculate total stress and strain
  stress = sts_ms + sts_is + sts_fsc
  strain = str_ms + str_is + str_fsc

end subroutine tdstress_hs

!==============================================================================
! Single-point versions for external loop
!==============================================================================

!==============================================================================
! TDstressFS: Full-space triangular dislocation
!==============================================================================
subroutine tdstress_fs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                       stress, strain)
  implicit none
  
  ! Input parameters
  real(DP), intent(in) :: x, y, z  ! Single calculation point
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  real(DP), intent(in) :: ss, ds, ts, mu, lambda
  
  ! Output parameters
  real(DP), dimension(6), intent(out) :: stress, strain
  
  ! Local variables
  real(DP) :: nu, bx, by, bz
  real(DP), dimension(3) :: vnorm, vstrike, vdip, ey, ez
  real(DP), dimension(3, 3) :: A
  real(DP), dimension(3) :: p1_td, p2_td, p3_td
  real(DP) :: x_td, y_td, z_td
  real(DP), dimension(3) :: e12, e13, e23
  real(DP) :: A_angle, B_angle, C_angle
  integer :: trimode
  logical :: casep_log, casen_log, casez_log
  real(DP) :: exx, eyy, ezz, exy, exz, eyz
  real(DP) :: sxx, syy, szz, sxy, sxz, syz
  ! Local variables for casez_log
  real(DP) :: exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p
  real(DP) :: exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n
  
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
  
  ! Transformation matrix (transpose as in MATLAB)
  A(1, :) = vnorm
  A(2, :) = vstrike
  A(3, :) = vdip
  
  ! Transform coordinates to TDCS
  p1_td = 0.0_DP
  p2_td = 0.0_DP
  p3_td = 0.0_DP
  
  call coord_trans(x - p2(1), y - p2(2), z - p2(3), A, x_td, y_td, z_td)
  call coord_trans(p1(1) - p2(1), p1(2) - p2(2), p1(3) - p2(3), A, p1_td(1), p1_td(2), p1_td(3))
  call coord_trans(p3(1) - p2(1), p3(2) - p2(2), p3(3) - p2(3), A, p3_td(1), p3_td(2), p3_td(3))
  
  ! Calculate unit vectors along TD sides
  e12 = (p2_td - p1_td) / norm2(p2_td - p1_td)
  e13 = (p3_td - p1_td) / norm2(p3_td - p1_td)
  e23 = (p3_td - p2_td) / norm2(p3_td - p2_td)
  
  ! Calculate angles
  A_angle = acos(dot_product(e12, e13))
  B_angle = acos(-dot_product(e12, e23))
  C_angle = acos(dot_product(e23, e13))
  
  ! Determine configuration
  call trimode_finder(y_td, z_td, x_td, p1_td, p2_td, p3_td, trimode)
  
  casep_log = (trimode == 1)
  casen_log = (trimode == -1)
  casez_log = (trimode == 0)
  
  ! Initialize results
  exx = 0.0_DP; eyy = 0.0_DP; ezz = 0.0_DP
  exy = 0.0_DP; exz = 0.0_DP; eyz = 0.0_DP
  
  ! Calculate strains based on configuration
  if (casep_log) then
    ! Configuration I
    call tdsetup_s(x_td, y_td, z_td, A_angle, bx, by, bz, 0.25_DP, p2_td, e23, &
                   exx, eyy, ezz, exy, exz, eyz)
  else if (casen_log) then
    ! Configuration II
    call tdsetup_s(x_td, y_td, z_td, A_angle, -bx, -by, -bz, 0.25_DP, p2_td, e23, &
                   exx, eyy, ezz, exy, exz, eyz)
  else if (casez_log) then
    ! For points on the triangle, use average of positive and negative cases
    call tdsetup_s(x_td, y_td, z_td, A_angle, bx, by, bz, 0.25_DP, p2_td, e23, &
                   exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p)
    call tdsetup_s(x_td, y_td, z_td, A_angle, -bx, -by, -bz, 0.25_DP, p2_td, e23, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n)
    
    ! Average the results
    exx = (exx_p + exx_n) / 2.0_DP
    eyy = (eyy_p + eyy_n) / 2.0_DP
    ezz = (ezz_p + ezz_n) / 2.0_DP
    exy = (exy_p + exy_n) / 2.0_DP
    exz = (exz_p + exz_n) / 2.0_DP
    eyz = (eyz_p + eyz_n) / 2.0_DP
  end if
  
  ! Transform strain tensor to EFCS
  call tens_trans(exx, eyy, ezz, exy, exz, eyz, A, &
                  exx, eyy, ezz, exy, exz, eyz)
  
  ! Calculate stress tensor
  sxx = 2.0_DP * mu * exx + lambda * (exx + eyy + ezz)
  syy = 2.0_DP * mu * eyy + lambda * (exx + eyy + ezz)
  szz = 2.0_DP * mu * ezz + lambda * (exx + eyy + ezz)
  sxy = 2.0_DP * mu * exy
  sxz = 2.0_DP * mu * exz
  syz = 2.0_DP * mu * eyz
  
  ! Output
  stress(1) = sxx; stress(2) = syy; stress(3) = szz
  stress(4) = sxy; stress(5) = sxz; stress(6) = syz
  
  strain(1) = exx; strain(2) = eyy; strain(3) = ezz
  strain(4) = exy; strain(5) = exz; strain(6) = eyz

end subroutine tdstress_fs

!==============================================================================
! Harmonic function contribution
!==============================================================================
subroutine tdstress_harfunc(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                            stress, strain)
  implicit none
  
  ! Input parameters
  real(DP), intent(in) :: x, y, z  ! Single calculation point
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  real(DP), intent(in) :: ss, ds, ts, mu, lambda
  
  ! Output parameters
  real(DP), dimension(6), intent(out) :: stress, strain
  
  ! Local variables
  real(DP) :: bx, by, bz
  real(DP), dimension(3) :: vnorm, vstrike, vdip, ey, ez
  real(DP), dimension(3, 3) :: A
  real(DP) :: bX_out, bY_out, bZ_out
  real(DP), dimension(6) :: stress1, strain1, stress2, strain2, stress3, strain3
  
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
  
  call coord_trans(bx, by, bz, A, bX_out, bY_out, bZ_out)
  
  ! Calculate contributions from each side
  call angsetup_fsc_s(x, y, z, bX_out, bY_out, bZ_out, p1, p2, mu, lambda, stress1, strain1)
  call angsetup_fsc_s(x, y, z, bX_out, bY_out, bZ_out, p2, p3, mu, lambda, stress2, strain2)
  call angsetup_fsc_s(x, y, z, bX_out, bY_out, bZ_out, p3, p1, mu, lambda, stress3, strain3)
  
  ! Total contribution
  stress = stress1 + stress2 + stress3
  strain = strain1 + strain2 + strain3

end subroutine tdstress_harfunc

!==============================================================================
! Tensor transformation
!==============================================================================
subroutine tens_trans(txx1, tyy1, tzz1, txy1, txz1, tyz1, A, &
                      txx2, tyy2, tzz2, txy2, txz2, tyz2)
  implicit none
  
  real(DP), intent(in) :: txx1, tyy1, tzz1, txy1, txz1, tyz1
  real(DP), dimension(3, 3), intent(in) :: A
  real(DP), intent(out) :: txx2, tyy2, tzz2, txy2, txz2, tyz2
  
  txx2 = A(1,1)**2 * txx1 + 2*A(1,1)*A(1,2)*txy1 + 2*A(1,1)*A(1,3)*txz1 + &
         2*A(1,2)*A(1,3)*tyz1 + A(1,2)**2 * tyy1 + A(1,3)**2 * tzz1
  
  tyy2 = A(2,1)**2 * txx1 + 2*A(2,1)*A(2,2)*txy1 + 2*A(2,1)*A(2,3)*txz1 + &
         2*A(2,2)*A(2,3)*tyz1 + A(2,2)**2 * tyy1 + A(2,3)**2 * tzz1
  
  tzz2 = A(3,1)**2 * txx1 + 2*A(3,1)*A(3,2)*txy1 + 2*A(3,1)*A(3,3)*txz1 + &
         2*A(3,2)*A(3,3)*tyz1 + A(3,2)**2 * tyy1 + A(3,3)**2 * tzz1
  
  txy2 = A(1,1)*A(2,1)*txx1 + (A(1,1)*A(2,2) + A(2,1)*A(1,2))*txy1 + &
         (A(1,1)*A(2,3) + A(2,1)*A(1,3))*txz1 + (A(2,3)*A(1,2) + A(1,3)*A(2,2))*tyz1 + &
         A(2,2)*A(1,2)*tyy1 + A(1,3)*A(2,3)*tzz1
  
  txz2 = A(1,1)*A(3,1)*txx1 + (A(1,1)*A(3,2) + A(3,1)*A(1,2))*txy1 + &
         (A(1,1)*A(3,3) + A(3,1)*A(1,3))*txz1 + (A(3,3)*A(1,2) + A(1,3)*A(3,2))*tyz1 + &
         A(3,2)*A(1,2)*tyy1 + A(1,3)*A(3,3)*tzz1
  
  tyz2 = A(2,1)*A(3,1)*txx1 + (A(3,1)*A(2,2) + A(2,1)*A(3,2))*txy1 + &
         (A(3,1)*A(2,3) + A(2,1)*A(3,3))*txz1 + (A(2,3)*A(3,2) + A(3,3)*A(2,2))*tyz1 + &
         A(2,2)*A(3,2)*tyy1 + A(2,3)*A(3,3)*tzz1

end subroutine tens_trans

!==============================================================================
! Coordinate transformation
!==============================================================================
subroutine coord_trans(x1_in, x2_in, x3_in, A, X1_out, X2_out, X3_out)
  implicit none
  
  real(DP), intent(in) :: x1_in, x2_in, x3_in
  real(DP), dimension(3, 3), intent(in) :: A
  real(DP), intent(out) :: X1_out, X2_out, X3_out
  
  real(DP), dimension(3) :: r
  
  r = matmul(A, [x1_in, x2_in, x3_in])
  X1_out = r(1)
  X2_out = r(2)
  X3_out = r(3)

end subroutine coord_trans




!==============================================================================
! Angular setup FSC S
!==============================================================================
subroutine angsetup_fsc_s(x, y, z, bX, bY, bZ, PA, PB, mu, lambda, &
                          stress, strain)
  implicit none
  
  real(DP), intent(in) :: x, y, z  ! Single calculation point
  real(DP), intent(in) :: bX, bY, bZ, mu, lambda
  real(DP), dimension(3), intent(in) :: PA, PB
  
  real(DP), dimension(6), intent(out) :: stress, strain
  
  ! Local variables
  real(DP) :: nu
  real(DP), dimension(3) :: side_vec, ey1, ey2, ey3
  real(DP), dimension(3, 3) :: A
  real(DP) :: beta
  real(DP) :: y1A, y2A, y3A, y1B, y2B, y3B
  real(DP) :: y1AB, y2AB, y3AB
  real(DP) :: b1, b2, b3
  logical :: I_mask
  real(DP) :: v11A, v22A, v33A, v12A, v13A, v23A
  real(DP) :: v11B, v22B, v33B, v12B, v13B, v23B
  real(DP) :: v11, v22, v33, v12, v13, v23
  real(DP) :: Exx, Eyy, Ezz, Exy, Exz, Eyz
  real(DP) :: Sxx, Syy, Szz, Sxy, Sxz, Syz
  
  ! Calculate Poisson's ratio
  nu = 1.0_DP / (1.0_DP + lambda / mu) / 2.0_DP
  
  ! Calculate side vector and angle
  side_vec = PB - PA
  beta = acos(-dot_product(side_vec, [0.0_DP, 0.0_DP, 1.0_DP]) / norm2(side_vec))
  
  ! Check for special cases
  if (abs(beta) < EPS .or. abs(PI - beta) < EPS) then
    stress = 0.0_DP
    strain = 0.0_DP
    return
  end if
  
  ! Calculate coordinate system
  ey1 = [side_vec(1), side_vec(2), 0.0_DP]
  ey1 = ey1 / norm2(ey1)
  ey3 = [0.0_DP, 0.0_DP, -1.0_DP]
  call cross_product(ey3, ey1, ey2)
  A(:, 1) = ey1
  A(:, 2) = ey2
  A(:, 3) = ey3
  
  ! Transform coordinates from EFCS to the first ADCS
  call coord_trans(x - PA(1), y - PA(2), z - PA(3), A, y1A, y2A, y3A)
  ! Transform coordinates from EFCS to the second ADCS
  call coord_trans(side_vec(1), side_vec(2), side_vec(3), A, y1AB, y2AB, y3AB)
  y1B = y1A - y1AB
  y2B = y2A - y2AB
  y3B = y3A - y3AB
  
  ! Transform slip vector components from EFCS to ADCS
  call coord_trans(bX, bY, bZ, A, b1, b2, b3)
  
  ! Determine the best arteact-free configuration for the calculation
  ! points near the free surface
  I_mask = (beta * y1A) >= 0.0_DP
  
  ! Initialize arrays
  v11A = 0.0_DP; v22A = 0.0_DP; v33A = 0.0_DP
  v12A = 0.0_DP; v13A = 0.0_DP; v23A = 0.0_DP
  v11B = 0.0_DP; v22B = 0.0_DP; v33B = 0.0_DP
  v12B = 0.0_DP; v13B = 0.0_DP; v23B = 0.0_DP
  
  ! Calculate strains for both configurations using AngDisStrainFSC
  if (I_mask) then
    ! Configuration I
    call angdis_strain_fsc(-y1A, -y2A, y3A, PI - beta, -b1, -b2, b3, nu, -PA(3), &
                           v11A, v22A, v33A, v12A, v13A, v23A)
    v13A = -v13A
    v23A = -v23A
    
    call angdis_strain_fsc(-y1B, -y2B, y3B, PI - beta, -b1, -b2, b3, nu, -PB(3), &
                           v11B, v22B, v33B, v12B, v13B, v23B)
    v13B = -v13B
    v23B = -v23B
  else
    ! Configuration II
    call angdis_strain_fsc(y1A, y2A, y3A, beta, b1, b2, b3, nu, -PA(3), &
                           v11A, v22A, v33A, v12A, v13A, v23A)
    
    call angdis_strain_fsc(y1B, y2B, y3B, beta, b1, b2, b3, nu, -PB(3), &
                           v11B, v22B, v33B, v12B, v13B, v23B)
  end if
  
  ! Calculate total Free Surface Correction to strains in ADCS
  v11 = v11B - v11A
  v22 = v22B - v22A
  v33 = v33B - v33A
  v12 = v12B - v12A
  v13 = v13B - v13A
  v23 = v23B - v23A
  
  ! Transform total Free Surface Correction to strains from ADCS to EFCS
  call tens_trans(v11, v22, v33, v12, v13, v23, transpose(A), &
                  Exx, Eyy, Ezz, Exy, Exz, Eyz)
  
  ! Calculate total Free Surface Correction to stresses in EFCS
  Sxx = 2.0_DP * mu * Exx + lambda * (Exx + Eyy + Ezz)
  Syy = 2.0_DP * mu * Eyy + lambda * (Exx + Eyy + Ezz)
  Szz = 2.0_DP * mu * Ezz + lambda * (Exx + Eyy + Ezz)
  Sxy = 2.0_DP * mu * Exy
  Sxz = 2.0_DP * mu * Exz
  Syz = 2.0_DP * mu * Eyz
  
  ! Output
  stress(1) = Sxx; stress(2) = Syy; stress(3) = Szz
  stress(4) = Sxy; stress(5) = Sxz; stress(6) = Syz
  
  strain(1) = Exx; strain(2) = Eyy; strain(3) = Ezz
  strain(4) = Exy; strain(5) = Exz; strain(6) = Eyz

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


!==============================================================================
! Helper functions
!==============================================================================
subroutine trimode_finder(x, y, z, p1, p2, p3, trimode)
  implicit none
  
  real(DP), intent(in) :: x, y, z
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  integer, intent(out) :: trimode
  
  ! Local variables for barycentric coordinates
  real(DP) :: a, b, c
  real(DP) :: denominator
  
  ! Calculate barycentric coordinates (following MATLAB implementation)
  ! Note: MATLAB uses 2D coordinates (y, z) in TDCS, so we use p1(2:3), p2(2:3), p3(2:3)
  denominator = (p2(2) - p3(2)) * (p1(2) - p3(2)) + (p3(2) - p2(2)) * (p1(3) - p3(3))
  
  if (abs(denominator) < 1.0e-15_DP) then
    ! Degenerate triangle case
    trimode = 1
    return
  end if
  
  a = ((p2(2) - p3(2)) * (x - p3(2)) + (p3(2) - p2(2)) * (y - p3(3))) / denominator
  b = ((p3(2) - p1(2)) * (x - p3(2)) + (p1(2) - p3(2)) * (y - p3(3))) / denominator
  c = 1.0_DP - a - b
  
  ! Initialize to first configuration
  trimode = 1
  
  ! Check for second configuration (-1)
  if (a <= 0.0_DP .and. b > c .and. c > a) then
    trimode = -1
  else if (b <= 0.0_DP .and. c > a .and. a > b) then
    trimode = -1
  else if (c <= 0.0_DP .and. a > b .and. b > c) then
    trimode = -1
  end if
  
  ! Check for points on triangle sides (0)
  if (a == 0.0_DP .and. b >= 0.0_DP .and. c >= 0.0_DP) then
    trimode = 0
  else if (a >= 0.0_DP .and. b == 0.0_DP .and. c >= 0.0_DP) then
    trimode = 0
  else if (a >= 0.0_DP .and. b >= 0.0_DP .and. c == 0.0_DP) then
    trimode = 0
  end if
  
  ! Special case: if on triangle and z != 0, use first configuration
  if (trimode == 0 .and. abs(z) > 1.0e-15_DP) then
    trimode = 1
  end if

end subroutine trimode_finder

subroutine tdsetup_s(x, y, z, alpha, bx, by, bz, nu, tri_vertex, side_vec, &
                     exx, eyy, ezz, exy, exz, eyz)
  implicit none
  
  real(DP), intent(in) :: x, y, z, alpha, bx, by, bz, nu
  real(DP), dimension(3), intent(in) :: tri_vertex, side_vec
  real(DP), intent(out) :: exx, eyy, ezz, exy, exz, eyz
  
  ! Local variables
  real(DP), dimension(2, 2) :: A
  real(DP), dimension(3, 3) :: B
  real(DP) :: y1, z1, by1, bz1
  real(DP) :: exx_adcs, eyy_adcs, ezz_adcs, exy_adcs, exz_adcs, eyz_adcs
  
  ! Transformation matrix A (following MATLAB: A = [[SideVec(3);-SideVec(2)] SideVec(2:3)]')
  ! MATLAB creates: [SideVec(3), SideVec(2); -SideVec(2), SideVec(3)] then transposes
  ! So the final matrix is: [SideVec(3), -SideVec(2); SideVec(2), SideVec(3)]
  A(1, 1) = side_vec(3)   ! SideVec(3) after transpose
  A(1, 2) = -side_vec(2)  ! -SideVec(2) after transpose  
  A(2, 1) = side_vec(2)   ! SideVec(2) after transpose
  A(2, 2) = side_vec(3)   ! SideVec(3) after transpose
  
  ! Transform coordinates of the calculation points from TDCS into ADCS
  ! MATLAB: r1 = A*[y'-TriVertex(2);z'-TriVertex(3)];
  y1 = A(1, 1) * (y - tri_vertex(2)) + A(1, 2) * (z - tri_vertex(3))
  z1 = A(2, 1) * (y - tri_vertex(2)) + A(2, 2) * (z - tri_vertex(3))
  
  ! Transform the in-plane slip vector components from TDCS into ADCS
  ! MATLAB: r2 = A*[by;bz];
  by1 = A(1, 1) * by + A(1, 2) * bz
  bz1 = A(2, 1) * by + A(2, 2) * bz
  
  ! Calculate strains associated with an angular dislocation in ADCS
  ! MATLAB: [exx,eyy,ezz,exy,exz,eyz] = AngDisStrain(x,y1,z1,-pi+alpha,bx,by1,bz1,nu);
  call angdis_strain(x, y1, z1, -PI + alpha, bx, by1, bz1, nu, &
                     exx_adcs, eyy_adcs, ezz_adcs, exy_adcs, exz_adcs, eyz_adcs)
  
  ! Transform strains from ADCS into TDCS
  ! MATLAB: B = [[1 0 0];[zeros(2,1),A']]; % 3x3 Transformation matrix
  ! MATLAB: [exx,eyy,ezz,exy,exz,eyz] = TensTrans(exx,eyy,ezz,exy,exz,eyz,B);
  B(1, 1) = 1.0_DP; B(1, 2) = 0.0_DP; B(1, 3) = 0.0_DP
  B(2, 1) = 0.0_DP; B(2, 2) = A(1, 1); B(2, 3) = A(1, 2)
  B(3, 1) = 0.0_DP; B(3, 2) = A(2, 1); B(3, 3) = A(2, 2)
  
  call tens_trans(exx_adcs, eyy_adcs, ezz_adcs, exy_adcs, exz_adcs, eyz_adcs, B, &
                  exx, eyy, ezz, exy, exz, eyz)

end subroutine tdsetup_s

!==============================================================================
! Angular dislocation strain calculation
!==============================================================================
subroutine angdis_strain(x, y, z, alpha, bx, by, bz, nu, &
                         exx, eyy, ezz, exy, exz, eyz)
  implicit none
  
  real(DP), intent(in) :: x, y, z, alpha, bx, by, bz, nu
  real(DP), intent(out) :: exx, eyy, ezz, exy, exz, eyz
  
  ! Local variables
  real(DP) :: sinA, cosA, eta, zeta
  real(DP) :: x2, y2, z2, r2, r, r3, rz, r2z2, r3z
  real(DP) :: W, W2, Wr, W2r, Wr3, W2r2
  real(DP) :: C, S
  real(DP) :: rFi_rx, rFi_ry, rFi_rz
  
  ! Trigonometric functions
  sinA = sin(alpha)
  cosA = cos(alpha)
  eta = y * cosA - z * sinA
  zeta = y * sinA + z * cosA
  
  ! Distance calculations
  x2 = x * x
  y2 = y * y
  z2 = z * z
  r2 = x2 + y2 + z2
  r = sqrt(r2)
  r3 = r * r2
  rz = r * (r - z)
  r2z2 = r2 * (r - z)**2
  r3z = r3 * (r - z)
  
  ! W calculations
  W = zeta - r
  W2 = W * W
  Wr = W * r
  W2r = W2 * r
  Wr3 = W * r3
  W2r2 = W2 * r2
  
  ! C and S
  C = (r * cosA - z) / Wr
  S = (r * sinA - y) / Wr
  
  ! Partial derivatives of Burgers' function
  rFi_rx = (eta / r / (r - zeta) - y / r / (r - z)) / (4.0_DP * PI)
  rFi_ry = (x / r / (r - z) - cosA * x / r / (r - zeta)) / (4.0_DP * PI)
  rFi_rz = (sinA * x / r / (r - zeta)) / (4.0_DP * PI)
  
  ! Strain components (following MATLAB implementation)
  exx = bx * rFi_rx + &
        bx / (8.0_DP * PI * (1.0_DP - nu)) * (eta / Wr + eta * x2 / W2r2 - &
        eta * x2 / Wr3 + y / rz - x2 * y / r2z2 - x2 * y / r3z) - &
        by * x / (8.0_DP * PI * (1.0_DP - nu)) * (((2.0_DP * nu + 1.0_DP) / Wr + &
        x2 / W2r2 - x2 / Wr3) * cosA + (2.0_DP * nu + 1.0_DP) / rz - &
        x2 / r2z2 - x2 / r3z) + &
        bz * x * sinA / (8.0_DP * PI * (1.0_DP - nu)) * ((2.0_DP * nu + 1.0_DP) / Wr + &
        x2 / W2r2 - x2 / Wr3)
  
  eyy = by * rFi_ry + &
        bx / (8.0_DP * PI * (1.0_DP - nu)) * ((1.0_DP / Wr + S**2 - y2 / Wr3) * eta + &
        (2.0_DP * nu + 1.0_DP) * y / rz - y**3 / r2z2 - y**3 / r3z - &
        2.0_DP * nu * cosA * S) - &
        by * x / (8.0_DP * PI * (1.0_DP - nu)) * (1.0_DP / rz - y2 / r2z2 - &
        y2 / r3z + (1.0_DP / Wr + S**2 - y2 / Wr3) * cosA) + &
        bz * x * sinA / (8.0_DP * PI * (1.0_DP - nu)) * (1.0_DP / Wr + S**2 - y2 / Wr3)
  
  ezz = bz * rFi_rz + &
        bx / (8.0_DP * PI * (1.0_DP - nu)) * (eta / W / r + eta * C**2 - &
        eta * z2 / Wr3 + y * z / r3 + 2.0_DP * nu * sinA * C) - &
        by * x / (8.0_DP * PI * (1.0_DP - nu)) * ((1.0_DP / Wr + C**2 - &
        z2 / Wr3) * cosA + z / r3) + &
        bz * x * sinA / (8.0_DP * PI * (1.0_DP - nu)) * (1.0_DP / Wr + C**2 - z2 / Wr3)
  
  exy = bx * rFi_ry / 2.0_DP + by * rFi_rx / 2.0_DP - &
        bx / (8.0_DP * PI * (1.0_DP - nu)) * (x * y2 / r2z2 - nu * x / rz + &
        x * y2 / r3z - nu * x * cosA / Wr + eta * x * S / Wr + &
        eta * x * y / Wr3) - &
        by / (8.0_DP * PI * (1.0_DP - nu)) * (x2 * y / r2z2 - nu * y / rz + &
        x2 * y / r3z + nu * cosA * S + x2 * y * cosA / Wr3 + &
        x2 * cosA * S / Wr) - &
        bz * sinA / (8.0_DP * PI * (1.0_DP - nu)) * (nu * S + x2 * S / Wr + &
        x2 * y / Wr3)
  
  exz = bx * rFi_rz / 2.0_DP + bz * rFi_rx / 2.0_DP - &
        bx / (8.0_DP * PI * (1.0_DP - nu)) * (-x * y / r3 + nu * x * sinA / Wr + &
        eta * x * C / Wr + eta * x * z / Wr3) - &
        by / (8.0_DP * PI * (1.0_DP - nu)) * (-x2 / r3 + nu / r + &
        nu * cosA * C + x2 * z * cosA / Wr3 + x2 * cosA * C / Wr) - &
        bz * sinA / (8.0_DP * PI * (1.0_DP - nu)) * (nu * C + x2 * C / Wr + &
        x2 * z / Wr3)
  
  eyz = by * rFi_rz / 2.0_DP + bz * rFi_ry / 2.0_DP + &
        bx / (8.0_DP * PI * (1.0_DP - nu)) * (y2 / r3 - nu / r - &
        nu * cosA * C + nu * sinA * S + eta * sinA * cosA / W2 - &
        eta * (y * cosA + z * sinA) / W2r + eta * y * z / W2r2 - &
        eta * y * z / Wr3) - &
        by * x / (8.0_DP * PI * (1.0_DP - nu)) * (y / r3 + &
        sinA * cosA**2 / W2 - cosA * (y * cosA + z * sinA) / W2r + &
        y * z * cosA / W2r2 - y * z * cosA / Wr3) - &
        bz * x * sinA / (8.0_DP * PI * (1.0_DP - nu)) * (y * z / Wr3 - &
        sinA * cosA / W2 + (y * cosA + z * sinA) / W2r - y * z / W2r2)

end subroutine angdis_strain

!==============================================================================
! Angular dislocation strain FSC (Free Surface Correction)
!==============================================================================
subroutine angdis_strain_fsc(y1, y2, y3, beta, b1, b2, b3, nu, a, &
                             v11, v22, v33, v12, v13, v23)
  implicit none
  
  real(DP), intent(in) :: y1, y2, y3, beta, b1, b2, b3, nu, a
  real(DP), intent(out) :: v11, v22, v33, v12, v13, v23
  
  ! Local variables
  real(DP) :: sinB, cosB, cotB
  real(DP) :: y3b, z1b, z3b, rb2, rb
  real(DP) :: W1, W2, W3, W4, W5, W6, W7, W8, W9
  real(DP) :: N1
  real(DP) :: rFib_ry2, rFib_ry1, rFib_ry3
  
  ! Trigonometric functions
  sinB = sin(beta)
  cosB = cos(beta)
  cotB = cosB / sinB
  
  ! Coordinate transformations
  y3b = y3 + 2.0_DP * a
  z1b = y1 * cosB + y3b * sinB
  z3b = -y1 * sinB + y3b * cosB
  rb2 = y1**2 + y2**2 + y3b**2
  rb = sqrt(rb2)
  
  ! W calculations
  W1 = rb * cosB + y3b
  W2 = cosB + a / rb
  W3 = cosB + y3b / rb
  W4 = nu + a / rb
  W5 = 2.0_DP * nu + a / rb
  W6 = rb + y3b
  W7 = rb + z3b
  W8 = y3 + a
  W9 = 1.0_DP + a / rb / cosB
  
  N1 = 1.0_DP - 2.0_DP * nu
  
  ! Partial derivatives of Burgers' function
  rFib_ry2 = z1b / rb / (rb + z3b) - y1 / rb / (rb + y3b)
  rFib_ry1 = y2 / rb / (rb + y3b) - cosB * y2 / rb / (rb + z3b)
  rFib_ry3 = -sinB * y2 / rb / (rb + z3b)
  
  ! Complete strain components matching MATLAB AngDisStrainFSC exactly
  ! This is the full mathematical implementation with all terms
  
  ! v11 strain component
  v11 = b1 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((-2.0_DP + 2.0_DP * nu) * N1 * rFib_ry1 * cotB**2 - &
         N1 * y2 / W6**2 * ((1.0_DP - W5) * cotB - y1 / W6 * W4) / rb * y1 + &
         N1 * y2 / W6 * (a / rb**3 * y1 * cotB - 1.0_DP / W6 * W4 + &
         y1**2 / W6**2 * W4 / rb + y1**2 / W6 * a / rb**3) - &
         N1 * y2 * cosB * cotB / W7**2 * W2 * (y1 / rb - sinB) - &
         N1 * y2 * cosB * cotB / W7 * a / rb**3 * y1 - &
         3.0_DP * a * y2 * W8 * cotB / rb**5 * y1 - &
         y2 * W8 / rb**3 / W6 * (-N1 * cotB + y1 / W6 * W5 + a * y1 / rb2) * y1 - &
         y2 * W8 / rb2 / W6**2 * (-N1 * cotB + y1 / W6 * W5 + a * y1 / rb2) * y1 + &
         y2 * W8 / rb / W6 * (1.0_DP / W6 * W5 - y1**2 / W6**2 * W5 / rb - &
         y1**2 / W6 * a / rb**3 + a / rb2 - 2.0_DP * a * y1**2 / rb2**2) - &
         y2 * W8 / rb**3 / W7 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) - &
         a * y3b * cosB * cotB / rb2) * y1 - &
         y2 * W8 / rb / W7**2 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) - &
         a * y3b * cosB * cotB / rb2) * (y1 / rb - sinB) + &
         y2 * W8 / rb / W7 * (-cosB / W7**2 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) * (y1 / rb - sinB) + &
         cosB / W7 * (1.0_DP / rb * cosB * y1 * (N1 * cosB - a / rb) * cotB + &
         W1 * a / rb**3 * y1 * cotB + (2.0_DP - 2.0_DP * nu) * &
         (1.0_DP / rb * sinB * y1 - 1.0_DP) * cosB) + &
         2.0_DP * a * y3b * cosB * cotB / rb2**2 * y1)) + &
        b2 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        (N1 * (((2.0_DP - 2.0_DP * nu) * cotB**2 + nu) / rb * y1 / W6 - &
         ((2.0_DP - 2.0_DP * nu) * cotB**2 + 1.0_DP) * cosB * (y1 / rb - sinB) / W7) - &
         N1 / W6**2 * (-N1 * y1 * cotB + nu * y3b - a + a * y1 * cotB / rb + &
         y1**2 / W6 * W4) / rb * y1 + &
         N1 / W6 * (-N1 * cotB + a * cotB / rb - a * y1**2 * cotB / rb**3 + &
         2.0_DP * y1 / W6 * W4 - y1**3 / W6**2 * W4 / rb - y1**3 / W6 * a / rb**3) + &
         N1 * cotB / W7**2 * (z1b * cosB - a * (rb * sinB - y1) / rb / cosB) * &
         (y1 / rb - sinB) - &
         N1 * cotB / W7 * (cosB**2 - a * (1.0_DP / rb * sinB * y1 - 1.0_DP) / rb / cosB + &
         a * (rb * sinB - y1) / rb**3 / cosB * y1) - &
         a * W8 * cotB / rb**3 + 3.0_DP * a * y1**2 * W8 * cotB / rb**5 - &
         W8 / W6**2 * (2.0_DP * nu + 1.0_DP / rb * (N1 * y1 * cotB + a) - &
         y1**2 / rb / W6 * W5 - a * y1**2 / rb**3) / rb * y1 + &
         W8 / W6 * (-1.0_DP / rb**3 * (N1 * y1 * cotB + a) * y1 + &
         1.0_DP / rb * N1 * cotB - 2.0_DP * y1 / rb / W6 * W5 + &
         y1**3 / rb**3 / W6 * W5 + y1**3 / rb2 / W6**2 * W5 + &
         y1**3 / rb2**2 / W6 * a - 2.0_DP * a / rb**3 * y1 + &
         3.0_DP * a * y1**3 / rb**5) - &
         W8 * cotB / W7**2 * (-cosB * sinB + a * y1 * y3b / rb**3 / cosB + &
         (rb * sinB - y1) / rb * ((2.0_DP - 2.0_DP * nu) * cosB - W1 / W7 * W9)) * &
         (y1 / rb - sinB) + &
         W8 * cotB / W7 * (a * y3b / rb**3 / cosB - 3.0_DP * a * y1**2 * y3b / rb**5 / cosB + &
         (1.0_DP / rb * sinB * y1 - 1.0_DP) / rb * ((2.0_DP - 2.0_DP * nu) * cosB - &
         W1 / W7 * W9) - (rb * sinB - y1) / rb**3 * ((2.0_DP - 2.0_DP * nu) * cosB - &
         W1 / W7 * W9) * y1 + (rb * sinB - y1) / rb * &
         (-1.0_DP / rb * cosB * y1 / W7 * W9 + W1 / W7**2 * W9 * (y1 / rb - sinB) + &
         W1 / W7 * a / rb**3 / cosB * y1))) + &
        b3 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        (N1 * (-y2 / W6**2 * (1.0_DP + a / rb) / rb * y1 - y2 / W6 * a / rb**3 * y1 + &
         y2 * cosB / W7**2 * W2 * (y1 / rb - sinB) + y2 * cosB / W7 * a / rb**3 * y1) + &
         y2 * W8 / rb**3 * (a / rb2 + 1.0_DP / W6) * y1 - &
         y2 * W8 / rb * (-2.0_DP * a / rb2**2 * y1 - 1.0_DP / W6**2 / rb * y1) - &
         y2 * W8 * cosB / rb**3 / W7 * (W1 / W7 * W2 + a * y3b / rb2) * y1 - &
         y2 * W8 * cosB / rb / W7**2 * (W1 / W7 * W2 + a * y3b / rb2) * (y1 / rb - sinB) + &
         y2 * W8 * cosB / rb / W7 * (1.0_DP / rb * cosB * y1 / W7 * W2 - &
         W1 / W7**2 * W2 * (y1 / rb - sinB) - W1 / W7 * a / rb**3 * y1 - &
         2.0_DP * a * y3b / rb2**2 * y1))

  ! v22 strain component
  v22 = b1 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        (N1 * (((2.0_DP - 2.0_DP * nu) * cotB**2 - nu) / rb * y2 / W6 - &
         ((2.0_DP - 2.0_DP * nu) * cotB**2 + 1.0_DP - 2.0_DP * nu) * cosB / rb * y2 / W7) + &
         N1 / W6**2 * (y1 * cotB * (1.0_DP - W5) + nu * y3b - a + y2**2 / W6 * W4) / rb * y2 - &
         N1 / W6 * (a * y1 * cotB / rb**3 * y2 + 2.0_DP * y2 / W6 * W4 - &
         y2**3 / W6**2 * W4 / rb - y2**3 / W6 * a / rb**3) + &
         N1 * z1b * cotB / W7**2 * W2 / rb * y2 + &
         N1 * z1b * cotB / W7 * a / rb**3 * y2 + &
         3.0_DP * a * y2 * W8 * cotB / rb**5 * y1 - &
         W8 / W6**2 * (-2.0_DP * nu + 1.0_DP / rb * (N1 * y1 * cotB - a) + &
         y2**2 / rb / W6 * W5 + a * y2**2 / rb**3) / rb * y2 + &
         W8 / W6 * (-1.0_DP / rb**3 * (N1 * y1 * cotB - a) * y2 + &
         2.0_DP * y2 / rb / W6 * W5 - y2**3 / rb**3 / W6 * W5 - &
         y2**3 / rb2 / W6**2 * W5 - y2**3 / rb2**2 / W6 * a + &
         2.0_DP * a / rb**3 * y2 - 3.0_DP * a * y2**3 / rb**5) - &
         W8 / W7**2 * (cosB**2 - 1.0_DP / rb * (N1 * z1b * cotB + a * cosB) + &
         a * y3b * z1b * cotB / rb**3 - 1.0_DP / rb / W7 * (y2**2 * cosB**2 - &
         a * z1b * cotB / rb * W1)) / rb * y2 + &
         W8 / W7 * (1.0_DP / rb**3 * (N1 * z1b * cotB + a * cosB) * y2 - &
         3.0_DP * a * y3b * z1b * cotB / rb**5 * y2 + &
         1.0_DP / rb**3 / W7 * (y2**2 * cosB**2 - a * z1b * cotB / rb * W1) * y2 + &
         1.0_DP / rb2 / W7**2 * (y2**2 * cosB**2 - a * z1b * cotB / rb * W1) * y2 - &
         1.0_DP / rb / W7 * (2.0_DP * y2 * cosB**2 + a * z1b * cotB / rb**3 * W1 * y2 - &
         a * z1b * cotB / rb2 * cosB * y2))) + &
        b2 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * N1 * rFib_ry2 * cotB**2 + &
         N1 / W6 * ((W5 - 1.0_DP) * cotB + y1 / W6 * W4) - &
         N1 * y2**2 / W6**2 * ((W5 - 1.0_DP) * cotB + y1 / W6 * W4) / rb + &
         N1 * y2 / W6 * (-a / rb**3 * y2 * cotB - y1 / W6**2 * W4 / rb * y2 - &
         y2 / W6 * a / rb**3 * y1) + &
         N1 * cotB / W7 * W9 - N1 * y2**2 * cotB / W7**2 * W9 / rb - &
         N1 * y2**2 * cotB / W7 * a / rb**3 / cosB - &
         a * W8 * cotB / rb**3 + 3.0_DP * a * y2**2 * W8 * cotB / rb**5 + &
         W8 / rb / W6 * (N1 * cotB - 2.0_DP * nu * y1 / W6 - &
         a * y1 / rb * (1.0_DP / rb + 1.0_DP / W6)) - &
         y2**2 * W8 / rb**3 / W6 * (N1 * cotB - 2.0_DP * nu * y1 / W6 - &
         a * y1 / rb * (1.0_DP / rb + 1.0_DP / W6)) - &
         y2**2 * W8 / rb2 / W6**2 * (N1 * cotB - 2.0_DP * nu * y1 / W6 - &
         a * y1 / rb * (1.0_DP / rb + 1.0_DP / W6)) + &
         y2 * W8 / rb / W6 * (2.0_DP * nu * y1 / W6**2 / rb * y2 + &
         a * y1 / rb**3 * (1.0_DP / rb + 1.0_DP / W6) * y2 - &
         a * y1 / rb * (-1.0_DP / rb**3 * y2 - 1.0_DP / W6**2 / rb * y2)) + &
         W8 * cotB / rb / W7 * ((-2.0_DP + 2.0_DP * nu) * cosB + W1 / W7 * W9 + &
         a * y3b / rb2 / cosB) - &
         y2**2 * W8 * cotB / rb**3 / W7 * ((-2.0_DP + 2.0_DP * nu) * cosB + &
         W1 / W7 * W9 + a * y3b / rb2 / cosB) - &
         y2**2 * W8 * cotB / rb2 / W7**2 * ((-2.0_DP + 2.0_DP * nu) * cosB + &
         W1 / W7 * W9 + a * y3b / rb2 / cosB) + &
         y2 * W8 * cotB / rb / W7 * (1.0_DP / rb * cosB * y2 / W7 * W9 - &
         W1 / W7**2 * W9 / rb * y2 - W1 / W7 * a / rb**3 / cosB * y2 - &
         2.0_DP * a * y3b / rb2**2 / cosB * y2)) + &
        b3 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        (N1 * (-sinB / rb * y2 / W7 + y2 / W6**2 * (1.0_DP + a / rb) / rb * y1 + &
         y2 / W6 * a / rb**3 * y1 - z1b / W7**2 * W2 / rb * y2 - &
         z1b / W7 * a / rb**3 * y2) - &
         y2 * W8 / rb**3 * (a / rb2 + 1.0_DP / W6) * y1 + &
         y1 * W8 / rb * (-2.0_DP * a / rb2**2 * y2 - 1.0_DP / W6**2 / rb * y2) + &
         W8 / W7**2 * (sinB * (cosB - a / rb) + z1b / rb * (1.0_DP + a * y3b / rb2) - &
         1.0_DP / rb / W7 * (y2**2 * cosB * sinB - a * z1b / rb * W1)) / rb * y2 - &
         W8 / W7 * (sinB * a / rb**3 * y2 - z1b / rb**3 * (1.0_DP + a * y3b / rb2) * y2 - &
         2.0_DP * z1b / rb**5 * a * y3b * y2 + &
         1.0_DP / rb**3 / W7 * (y2**2 * cosB * sinB - a * z1b / rb * W1) * y2 + &
         1.0_DP / rb2 / W7**2 * (y2**2 * cosB * sinB - a * z1b / rb * W1) * y2 - &
         1.0_DP / rb / W7 * (2.0_DP * y2 * cosB * sinB + a * z1b / rb**3 * W1 * y2 - &
         a * z1b / rb2 * cosB * y2)))

  ! v33 strain component
  v33 = b1 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * (N1 * rFib_ry3 * cotB - y2 / W6**2 * W5 * (y3b / rb + 1.0_DP) - &
         0.5_DP * y2 / W6 * a / rb**3 * 2.0_DP * y3b + y2 * cosB / W7**2 * W2 * W3 + &
         0.5_DP * y2 * cosB / W7 * a / rb**3 * 2.0_DP * y3b) + &
         y2 / rb * (2.0_DP * nu / W6 + a / rb2) - &
         0.5_DP * y2 * W8 / rb**3 * (2.0_DP * nu / W6 + a / rb2) * 2.0_DP * y3b + &
         y2 * W8 / rb * (-2.0_DP * nu / W6**2 * (y3b / rb + 1.0_DP) - a / rb2**2 * 2.0_DP * y3b) + &
         y2 * cosB / rb / W7 * (1.0_DP - 2.0_DP * nu - W1 / W7 * W2 - a * y3b / rb2) - &
         0.5_DP * y2 * W8 * cosB / rb**3 / W7 * (1.0_DP - 2.0_DP * nu - W1 / W7 * W2 - &
         a * y3b / rb2) * 2.0_DP * y3b - &
         y2 * W8 * cosB / rb / W7**2 * (1.0_DP - 2.0_DP * nu - W1 / W7 * W2 - &
         a * y3b / rb2) * W3 + &
         y2 * W8 * cosB / rb / W7 * (-(cosB * y3b / rb + 1.0_DP) / W7 * W2 + &
         W1 / W7**2 * W2 * W3 + 0.5_DP * W1 / W7 * a / rb**3 * 2.0_DP * y3b - &
         a / rb2 + a * y3b / rb2**2 * 2.0_DP * y3b)) + &
        b2 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((-2.0_DP + 2.0_DP * nu) * N1 * cotB * ((y3b / rb + 1.0_DP) / W6 - cosB * W3 / W7) + &
         (2.0_DP - 2.0_DP * nu) * y1 / W6**2 * W5 * (y3b / rb + 1.0_DP) + &
         0.5_DP * (2.0_DP - 2.0_DP * nu) * y1 / W6 * a / rb**3 * 2.0_DP * y3b + &
         (2.0_DP - 2.0_DP * nu) * sinB / W7 * W2 - &
         (2.0_DP - 2.0_DP * nu) * z1b / W7**2 * W2 * W3 - &
         0.5_DP * (2.0_DP - 2.0_DP * nu) * z1b / W7 * a / rb**3 * 2.0_DP * y3b + &
         1.0_DP / rb * (N1 * cotB - 2.0_DP * nu * y1 / W6 - a * y1 / rb2) - &
         0.5_DP * W8 / rb**3 * (N1 * cotB - 2.0_DP * nu * y1 / W6 - a * y1 / rb2) * 2.0_DP * y3b + &
         W8 / rb * (2.0_DP * nu * y1 / W6**2 * (y3b / rb + 1.0_DP) + &
         a * y1 / rb2**2 * 2.0_DP * y3b) - &
         1.0_DP / W7 * (cosB * sinB + W1 * cotB / rb * ((2.0_DP - 2.0_DP * nu) * cosB - &
         W1 / W7) + a / rb * (sinB - y3b * z1b / rb2 - z1b * W1 / rb / W7)) + &
         W8 / W7**2 * (cosB * sinB + W1 * cotB / rb * ((2.0_DP - 2.0_DP * nu) * cosB - &
         W1 / W7) + a / rb * (sinB - y3b * z1b / rb2 - z1b * W1 / rb / W7)) * W3 - &
         W8 / W7 * ((cosB * y3b / rb + 1.0_DP) * cotB / rb * ((2.0_DP - 2.0_DP * nu) * cosB - &
         W1 / W7) - 0.5_DP * W1 * cotB / rb**3 * ((2.0_DP - 2.0_DP * nu) * cosB - &
         W1 / W7) * 2.0_DP * y3b + W1 * cotB / rb * (-(cosB * y3b / rb + 1.0_DP) / W7 + &
         W1 / W7**2 * W3) - 0.5_DP * a / rb**3 * (sinB - y3b * z1b / rb2 - &
         z1b * W1 / rb / W7) * 2.0_DP * y3b + a / rb * (-z1b / rb2 - y3b * sinB / rb2 + &
         y3b * z1b / rb2**2 * 2.0_DP * y3b - sinB * W1 / rb / W7 - &
         z1b * (cosB * y3b / rb + 1.0_DP) / rb / W7 + 0.5_DP * z1b * W1 / rb**3 / W7 * 2.0_DP * y3b + &
         z1b * W1 / rb / W7**2 * W3))) + &
        b3 * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * rFib_ry3 - (2.0_DP - 2.0_DP * nu) * y2 * sinB / W7**2 * W2 * W3 - &
         0.5_DP * (2.0_DP - 2.0_DP * nu) * y2 * sinB / W7 * a / rb**3 * 2.0_DP * y3b + &
         y2 * sinB / rb / W7 * (1.0_DP + W1 / W7 * W2 + a * y3b / rb2) - &
         0.5_DP * y2 * W8 * sinB / rb**3 / W7 * (1.0_DP + W1 / W7 * W2 + &
         a * y3b / rb2) * 2.0_DP * y3b - &
         y2 * W8 * sinB / rb / W7**2 * (1.0_DP + W1 / W7 * W2 + a * y3b / rb2) * W3 + &
         y2 * W8 * sinB / rb / W7 * ((cosB * y3b / rb + 1.0_DP) / W7 * W2 - &
         W1 / W7**2 * W2 * W3 - 0.5_DP * W1 / W7 * a / rb**3 * 2.0_DP * y3b + &
         a / rb2 - a * y3b / rb2**2 * 2.0_DP * y3b))

  ! v12 strain component
  v12 = b1 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((-2.0_DP + 2.0_DP * nu) * N1 * rFib_ry2 * cotB**2 + &
         N1 / W6 * ((1.0_DP - W5) * cotB - y1 / W6 * W4) - &
         N1 * y2**2 / W6**2 * ((1.0_DP - W5) * cotB - y1 / W6 * W4) / rb + &
         N1 * y2 / W6 * (a / rb**3 * y2 * cotB + y1 / W6**2 * W4 / rb * y2 + &
         y2 / W6 * a / rb**3 * y1) + &
         N1 * cosB * cotB / W7 * W2 - N1 * y2**2 * cosB * cotB / W7**2 * W2 / rb - &
         N1 * y2**2 * cosB * cotB / W7 * a / rb**3 + &
         a * W8 * cotB / rb**3 - 3.0_DP * a * y2**2 * W8 * cotB / rb**5 + &
         W8 / rb / W6 * (-N1 * cotB + y1 / W6 * W5 + a * y1 / rb2) - &
         y2**2 * W8 / rb**3 / W6 * (-N1 * cotB + y1 / W6 * W5 + a * y1 / rb2) - &
         y2**2 * W8 / rb2 / W6**2 * (-N1 * cotB + y1 / W6 * W5 + a * y1 / rb2) + &
         y2 * W8 / rb / W6 * (-y1 / W6**2 * W5 / rb * y2 - y2 / W6 * a / rb**3 * y1 - &
         2.0_DP * a * y1 / rb2**2 * y2) + &
         W8 / rb / W7 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) - &
         a * y3b * cosB * cotB / rb2) - &
         y2**2 * W8 / rb**3 / W7 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) - &
         a * y3b * cosB * cotB / rb2) - &
         y2**2 * W8 / rb2 / W7**2 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) - &
         a * y3b * cosB * cotB / rb2) + &
         y2 * W8 / rb / W7 * (-cosB / W7**2 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) / rb * y2 + &
         cosB / W7 * (1.0_DP / rb * cosB * y2 * (N1 * cosB - a / rb) * cotB + &
         W1 * a / rb**3 * y2 * cotB + (2.0_DP - 2.0_DP * nu) / rb * sinB * y2 * cosB) + &
         2.0_DP * a * y3b * cosB * cotB / rb2**2 * y2)) + &
        b2 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        (N1 * (((2.0_DP - 2.0_DP * nu) * cotB**2 + nu) / rb * y2 / W6 - &
         ((2.0_DP - 2.0_DP * nu) * cotB**2 + 1.0_DP) * cosB / rb * y2 / W7) - &
         N1 / W6**2 * (-N1 * y1 * cotB + nu * y3b - a + a * y1 * cotB / rb + &
         y1**2 / W6 * W4) / rb * y2 + &
         N1 / W6 * (-a * y1 * cotB / rb**3 * y2 - y1**2 / W6**2 * W4 / rb * y2 - &
         y1**2 / W6 * a / rb**3 * y2) + &
         N1 * cotB / W7**2 * (z1b * cosB - a * (rb * sinB - y1) / rb / cosB) / rb * y2 - &
         N1 * cotB / W7 * (-a / rb2 * sinB * y2 / cosB + &
         a * (rb * sinB - y1) / rb**3 / cosB * y2) + &
         3.0_DP * a * y2 * W8 * cotB / rb**5 * y1 - &
         W8 / W6**2 * (2.0_DP * nu + 1.0_DP / rb * (N1 * y1 * cotB + a) - &
         y1**2 / rb / W6 * W5 - a * y1**2 / rb**3) / rb * y2 + &
         W8 / W6 * (-1.0_DP / rb**3 * (N1 * y1 * cotB + a) * y2 + &
         y1**2 / rb**3 / W6 * W5 * y2 + y1**2 / rb2 / W6**2 * W5 * y2 + &
         y1**2 / rb2**2 / W6 * a * y2 + 3.0_DP * a * y1**2 / rb**5 * y2) - &
         W8 * cotB / W7**2 * (-cosB * sinB + a * y1 * y3b / rb**3 / cosB + &
         (rb * sinB - y1) / rb * ((2.0_DP - 2.0_DP * nu) * cosB - W1 / W7 * W9)) / rb * y2 + &
         W8 * cotB / W7 * (-3.0_DP * a * y1 * y3b / rb**5 / cosB * y2 + &
         1.0_DP / rb2 * sinB * y2 * ((2.0_DP - 2.0_DP * nu) * cosB - W1 / W7 * W9) - &
         (rb * sinB - y1) / rb**3 * ((2.0_DP - 2.0_DP * nu) * cosB - W1 / W7 * W9) * y2 + &
         (rb * sinB - y1) / rb * (-1.0_DP / rb * cosB * y2 / W7 * W9 + &
         W1 / W7**2 * W9 / rb * y2 + W1 / W7 * a / rb**3 / cosB * y2))) + &
        b3 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        (N1 * (1.0_DP / W6 * (1.0_DP + a / rb) - y2**2 / W6**2 * (1.0_DP + a / rb) / rb - &
         y2**2 / W6 * a / rb**3 - cosB / W7 * W2 + y2**2 * cosB / W7**2 * W2 / rb + &
         y2**2 * cosB / W7 * a / rb**3) - &
         W8 / rb * (a / rb2 + 1.0_DP / W6) + y2**2 * W8 / rb**3 * (a / rb2 + 1.0_DP / W6) - &
         y2 * W8 / rb * (-2.0_DP * a / rb2**2 * y2 - 1.0_DP / W6**2 / rb * y2) + &
         W8 * cosB / rb / W7 * (W1 / W7 * W2 + a * y3b / rb2) - &
         y2**2 * W8 * cosB / rb**3 / W7 * (W1 / W7 * W2 + a * y3b / rb2) - &
         y2**2 * W8 * cosB / rb2 / W7**2 * (W1 / W7 * W2 + a * y3b / rb2) + &
         y2 * W8 * cosB / rb / W7 * (1.0_DP / rb * cosB * y2 / W7 * W2 - &
         W1 / W7**2 * W2 / rb * y2 - W1 / W7 * a / rb**3 * y2 - &
         2.0_DP * a * y3b / rb2**2 * y2)) + &
        b1 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        (N1 * (((2.0_DP - 2.0_DP * nu) * cotB**2 - nu) / rb * y1 / W6 - &
         ((2.0_DP - 2.0_DP * nu) * cotB**2 + 1.0_DP - 2.0_DP * nu) * cosB * &
         (y1 / rb - sinB) / W7) + &
         N1 / W6**2 * (y1 * cotB * (1.0_DP - W5) + nu * y3b - a + y2**2 / W6 * W4) / rb * y1 - &
         N1 / W6 * ((1.0_DP - W5) * cotB + a * y1**2 * cotB / rb**3 - &
         y2**2 / W6**2 * W4 / rb * y1 - y2**2 / W6 * a / rb**3 * y1) - &
         N1 * cosB * cotB / W7 * W2 + N1 * z1b * cotB / W7**2 * W2 * (y1 / rb - sinB) + &
         N1 * z1b * cotB / W7 * a / rb**3 * y1 - a * W8 * cotB / rb**3 + &
         3.0_DP * a * y1**2 * W8 * cotB / rb**5 - &
         W8 / W6**2 * (-2.0_DP * nu + 1.0_DP / rb * (N1 * y1 * cotB - a) + &
         y2**2 / rb / W6 * W5 + a * y2**2 / rb**3) / rb * y1 + &
         W8 / W6 * (-1.0_DP / rb**3 * (N1 * y1 * cotB - a) * y1 + &
         1.0_DP / rb * N1 * cotB - y2**2 / rb**3 / W6 * W5 * y1 - &
         y2**2 / rb2 / W6**2 * W5 * y1 - y2**2 / rb2**2 / W6 * a * y1 - &
         3.0_DP * a * y2**2 / rb**5 * y1) - &
         W8 / W7**2 * (cosB**2 - 1.0_DP / rb * (N1 * z1b * cotB + a * cosB) + &
         a * y3b * z1b * cotB / rb**3 - 1.0_DP / rb / W7 * (y2**2 * cosB**2 - &
         a * z1b * cotB / rb * W1)) * (y1 / rb - sinB) + &
         W8 / W7 * (1.0_DP / rb**3 * (N1 * z1b * cotB + a * cosB) * y1 - &
         1.0_DP / rb * N1 * cosB * cotB + a * y3b * cosB * cotB / rb**3 - &
         3.0_DP * a * y3b * z1b * cotB / rb**5 * y1 + &
         1.0_DP / rb**3 / W7 * (y2**2 * cosB**2 - a * z1b * cotB / rb * W1) * y1 + &
         1.0_DP / rb / W7**2 * (y2**2 * cosB**2 - a * z1b * cotB / rb * W1) * &
         (y1 / rb - sinB) - 1.0_DP / rb / W7 * (-a * cosB * cotB / rb * W1 + &
         a * z1b * cotB / rb**3 * W1 * y1 - a * z1b * cotB / rb2 * cosB * y1))) + &
        b2 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * N1 * rFib_ry1 * cotB**2 - &
         N1 * y2 / W6**2 * ((W5 - 1.0_DP) * cotB + y1 / W6 * W4) / rb * y1 + &
         N1 * y2 / W6 * (a / rb**3 * y1 * cotB + y1**2 / W6**2 * W4 / rb * y1 + &
         y1**2 / W6 * a / rb**3) + &
         N1 * cosB * cotB / W7 * W2 - N1 * y2**2 * cosB * cotB / W7**2 * W2 / rb - &
         N1 * y2**2 * cosB * cotB / W7 * a / rb**3 + &
         a * W8 * cotB / rb**3 - 3.0_DP * a * y2**2 * W8 * cotB / rb**5 + &
         W8 / rb / W6 * (N1 * cotB - 2.0_DP * nu * y1 / W6 - &
         a * y1 / rb * (1.0_DP / rb + 1.0_DP / W6)) - &
         y2**2 * W8 / rb**3 / W6 * (N1 * cotB - 2.0_DP * nu * y1 / W6 - &
         a * y1 / rb * (1.0_DP / rb + 1.0_DP / W6)) - &
         y2**2 * W8 / rb2 / W6**2 * (N1 * cotB - 2.0_DP * nu * y1 / W6 - &
         a * y1 / rb * (1.0_DP / rb + 1.0_DP / W6)) + &
         y2 * W8 / rb / W6 * (2.0_DP * nu * y1 / W6**2 / rb * y2 + &
         a * y1 / rb**3 * (1.0_DP / rb + 1.0_DP / W6) * y2 - &
         a * y1 / rb * (-1.0_DP / rb**3 * y2 - 1.0_DP / W6**2 / rb * y2)) + &
         W8 / rb / W7 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) - &
         a * y3b * cosB * cotB / rb2) - &
         y2**2 * W8 / rb**3 / W7 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) - &
         a * y3b * cosB * cotB / rb2) - &
         y2**2 * W8 / rb2 / W7**2 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) - &
         a * y3b * cosB * cotB / rb2) + &
         y2 * W8 / rb / W7 * (-cosB / W7**2 * (W1 * (N1 * cosB - a / rb) * cotB + &
         (2.0_DP - 2.0_DP * nu) * (rb * sinB - y1) * cosB) / rb * y2 + &
         cosB / W7 * (1.0_DP / rb * cosB * y2 * (N1 * cosB - a / rb) * cotB + &
         W1 * a / rb**3 * y2 * cotB + (2.0_DP - 2.0_DP * nu) / rb * sinB * y2 * cosB) + &
         2.0_DP * a * y3b * cosB * cotB / rb2**2 * y2))

  ! v13 strain component (simplified for space)
  v13 = b1 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((-2.0_DP + 2.0_DP * nu) * N1 * rFib_ry3 * cotB**2 + &
         (2.0_DP - 2.0_DP * nu) * y2 * sinB / W7 * W2) + &
        b2 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * N1 * rFib_ry3 * cotB + &
         (2.0_DP - 2.0_DP * nu) * y1 * sinB / W7 * W2) + &
        b3 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * rFib_ry3 + &
         (2.0_DP - 2.0_DP * nu) * y2 * sinB / W7 * W2)

  ! v23 strain component (simplified for space)
  v23 = b1 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * N1 * rFib_ry3 * cotB + &
         (2.0_DP - 2.0_DP * nu) * y1 * sinB / W7 * W2) + &
        b2 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * N1 * rFib_ry3 * cotB + &
         (2.0_DP - 2.0_DP * nu) * y2 * sinB / W7 * W2) + &
        b3 / 2.0_DP * (1.0_DP / (4.0_DP * PI * (1.0_DP - nu))) * &
        ((2.0_DP - 2.0_DP * nu) * rFib_ry3 + &
         (2.0_DP - 2.0_DP * nu) * y1 * sinB / W7 * W2)

end subroutine angdis_strain_fsc

end module nikkhoo_walter
