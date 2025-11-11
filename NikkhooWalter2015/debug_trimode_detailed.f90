program debug_trimode_detailed
  implicit none
  integer, parameter :: DP = selected_real_kind(15, 307)

  ! Triangle vertices in TDCS
  real(DP), dimension(3) :: p1_td, p2_td, p3_td

  ! Test points 8 and 9 coordinates
  real(DP) :: x8, y8, z8, x9, y9, z9

  ! Barycentric coordinates
  real(DP) :: a, b, c, denominator
  integer :: trimode

  ! Tolerances
  real(DP), parameter :: BARY_TOL = 1.0e-12_DP
  real(DP), parameter :: Z_TOL = 1.0e-10_DP

  ! For coordinate transformation
  real(DP), dimension(3) :: p1, p2, p3, vnorm, vstrike, vdip
  real(DP), dimension(3,3) :: A
  real(DP) :: x_td, y_td, z_td
  real(DP) :: bx, by, bz, nu

  print *, '============================================================'
  print *, 'Detailed trimode debug for points 8 and 9'
  print *, '============================================================'
  print *, ''

  ! Original triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Calculate normal and transformation (simplified from sub_nikkhoo.f90)
  vnorm = cross_product(p2 - p1, p3 - p1)
  vnorm = vnorm / norm_vec(vnorm)

  vstrike = p2 - p1
  vstrike = vstrike / norm_vec(vstrike)

  vdip = cross_product(vnorm, vstrike)

  ! Build transformation matrix
  A(1,:) = vnorm
  A(2,:) = vstrike
  A(3,:) = vdip

  ! Transform vertices to TDCS
  p1_td = matmul(A, p1)
  p2_td = matmul(A, p2)
  p3_td = matmul(A, p3)

  print *, 'Triangle in TDCS:'
  print *, '  p1_td = ', p1_td
  print *, '  p2_td = ', p2_td
  print *, '  p3_td = ', p3_td
  print *, ''

  ! Test point 8: (3.0, -3.0, -6.0)
  print *, '========== POINT 8 =========='
  x8 = 3.0_DP
  y8 = -3.0_DP
  z8 = -6.0_DP
  print *, 'Original coords: x=', x8, ', y=', y8, ', z=', z8

  ! Transform to TDCS
  x_td = A(1,1)*x8 + A(1,2)*y8 + A(1,3)*z8
  y_td = A(2,1)*x8 + A(2,2)*y8 + A(2,3)*z8
  z_td = A(3,1)*x8 + A(3,2)*y8 + A(3,3)*z8

  print *, 'TDCS coords: x_td=', x_td, ', y_td=', y_td, ', z_td=', z_td

  ! Calculate barycentric coordinates (note: called with y_td, z_td, x_td)
  call calc_trimode(y_td, z_td, x_td, p1_td, p2_td, p3_td, a, b, c, trimode)

  print *, 'Barycentric coords: a=', a, ', b=', b, ', c=', c
  print *, 'Sum (should be 1): ', a+b+c
  print *, 'trimode = ', trimode

  if (trimode == 0) then
    print *, '  --> PROBLEM: Classified as ON EDGE (will return NaN)'
    print *, '  Checking tolerance conditions:'
    print *, '    abs(a) < BARY_TOL:', abs(a) < BARY_TOL, ' (abs(a)=', abs(a), ')'
    print *, '    abs(b) < BARY_TOL:', abs(b) < BARY_TOL, ' (abs(b)=', abs(b), ')'
    print *, '    abs(c) < BARY_TOL:', abs(c) < BARY_TOL, ' (abs(c)=', abs(c), ')'
    print *, '    abs(z_td) > Z_TOL:', abs(z_td) > Z_TOL, ' (abs(z_td)=', abs(z_td), ')'
  else if (trimode == 1) then
    print *, '  --> OK: Configuration I (inside)'
  else if (trimode == -1) then
    print *, '  --> OK: Configuration II (outside)'
  end if
  print *, ''

  ! Test point 9: (-3.0, 3.0, -3.0)
  print *, '========== POINT 9 =========='
  x9 = -3.0_DP
  y9 = 3.0_DP
  z9 = -3.0_DP
  print *, 'Original coords: x=', x9, ', y=', y9, ', z=', z9

  ! Transform to TDCS
  x_td = A(1,1)*x9 + A(1,2)*y9 + A(1,3)*z9
  y_td = A(2,1)*x9 + A(2,2)*y9 + A(2,3)*z9
  z_td = A(3,1)*x9 + A(3,2)*y9 + A(3,3)*z9

  print *, 'TDCS coords: x_td=', x_td, ', y_td=', y_td, ', z_td=', z_td

  ! Calculate barycentric coordinates
  call calc_trimode(y_td, z_td, x_td, p1_td, p2_td, p3_td, a, b, c, trimode)

  print *, 'Barycentric coords: a=', a, ', b=', b, ', c=', c
  print *, 'Sum (should be 1): ', a+b+c
  print *, 'trimode = ', trimode

  if (trimode == 0) then
    print *, '  --> PROBLEM: Classified as ON EDGE (will return NaN)'
    print *, '  Checking tolerance conditions:'
    print *, '    abs(a) < BARY_TOL:', abs(a) < BARY_TOL, ' (abs(a)=', abs(a), ')'
    print *, '    abs(b) < BARY_TOL:', abs(b) < BARY_TOL, ' (abs(b)=', abs(b), ')'
    print *, '    abs(c) < BARY_TOL:', abs(c) < BARY_TOL, ' (abs(c)=', abs(c), ')'
    print *, '    abs(z_td) > Z_TOL:', abs(z_td) > Z_TOL, ' (abs(z_td)=', abs(z_td), ')'
  else if (trimode == 1) then
    print *, '  --> OK: Configuration I (inside)'
  else if (trimode == -1) then
    print *, '  --> OK: Configuration II (outside)'
  end if
  print *, ''

  print *, '============================================================'

contains

  subroutine calc_trimode(x, y, z, p1, p2, p3, a, b, c, trimode)
    real(DP), intent(in) :: x, y, z
    real(DP), dimension(3), intent(in) :: p1, p2, p3
    real(DP), intent(out) :: a, b, c
    integer, intent(out) :: trimode
    real(DP) :: denominator

    ! Calculate barycentric coordinates
    denominator = (p2(2) - p3(2)) * (p1(2) - p3(2)) + (p3(2) - p2(2)) * (p1(3) - p3(3))

    if (abs(denominator) < 1.0e-15_DP) then
      trimode = 1
      a = 0.0_DP
      b = 0.0_DP
      c = 1.0_DP
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

    ! Check for points on triangle sides (0) - WITH TOLERANCE AND BOUNDS CHECK
    ! Also verify point is within triangle (all coords in [0,1])
    if (abs(a) < BARY_TOL .and. b >= -BARY_TOL .and. b <= 1.0_DP + BARY_TOL .and. &
        c >= -BARY_TOL .and. c <= 1.0_DP + BARY_TOL) then
      trimode = 0
    else if (abs(b) < BARY_TOL .and. a >= -BARY_TOL .and. a <= 1.0_DP + BARY_TOL .and. &
             c >= -BARY_TOL .and. c <= 1.0_DP + BARY_TOL) then
      trimode = 0
    else if (abs(c) < BARY_TOL .and. a >= -BARY_TOL .and. a <= 1.0_DP + BARY_TOL .and. &
             b >= -BARY_TOL .and. b <= 1.0_DP + BARY_TOL) then
      trimode = 0
    end if

    ! Special case: if on triangle edge but z != 0, use first configuration
    if (trimode == 0 .and. abs(z) > Z_TOL) then
      trimode = 1
    end if
  end subroutine calc_trimode

  function cross_product(a, b) result(c)
    real(DP), dimension(3), intent(in) :: a, b
    real(DP), dimension(3) :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product

  function norm_vec(v) result(n)
    real(DP), dimension(3), intent(in) :: v
    real(DP) :: n
    n = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
  end function norm_vec

end program debug_trimode_detailed
