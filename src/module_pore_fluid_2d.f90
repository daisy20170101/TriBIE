!===============================================================================
! module_pore_fluid_2d.f90
!
! 2D pore-pressure diffusion on the fault plane for SEAS Benchmark BP8-QD-GS.
!
! Solves (BP8 Eq. 17, with the point source replaced by the Gaussian smear of
! Eq. 19):
!
!   dp/dt = alpha * Laplacian(p) + [qinj(t)/(beta*phi)] *
!           1/(2*pi*Lgauss^2) * exp(-(x2^2+x3^2)/(2*Lgauss^2))
!
! on a uniform structured grid covering the frictional domain
! Omega_f = (-lf,lf) x (-lf,lf), with zero-flux (Neumann) boundaries at the
! edges of that domain (Eq. 18) and qinj(t) = q0 for 0<=t<toff, 0 after
! (Eq. 20).
!
! This is genuinely different physics from src/phy3d_module_bp6.f90's
! compute_pf/compute_G, which solve a 1D line-source problem (pressure a
! function of a single distance coordinate) -- the wrong geometry for BP8's
! point injection at the center of a 2D fault. See NikkhooWalter2015-style
! debug/ scripts for why: BP8's own analytic short-time solution (Eq. 21)
! uses the exponential integral E1 of a *radial* argument x2^2+x3^2, not the
! erfc of a single coordinate.
!
! The grid is deliberately chosen at exactly the fault cell size (10 m, BP8
! Table 1) so profile/station outputs at BP8's required 10 m node spacing
! need no spatial interpolation for pressure; interpolation is only needed
! to couple this structured grid to the (generally off-grid) BEM triangle
! centroids of 3dtri_BP8.f90, via pf2d_interpolate (bilinear).
!
! Time stepping is explicit (forward Euler) with internal sub-stepping at a
! stable dt (dx^2/(4*alpha), which for BP8's parameters is ~500 s -- a small
! fraction of both toff (100 hr) and tf (30 days), so this remains cheap;
! see the verification note in example4/ for the accuracy check against the
! analytic short-time solution, Eq. 21).
!===============================================================================
module pore_fluid_2d
  implicit none
  private
  public :: DP2, pf2d_stable_dt, pf2d_step, pf2d_interpolate, pf2d_source, &
            pf2d_analytic, pf2d_expint_e1

  integer, parameter :: DP2 = kind(1.d0)
  real(DP2), parameter :: PI2 = 3.141592653589793238462643383279502884197_DP2

contains

  !-----------------------------------------------------------------------
  ! Gaussian source term, Eq. (19)'s RHS divided by (beta*phi) -- i.e. the
  ! contribution to dp/dt at (x2,x3), for injection rate qinj(t).
  !-----------------------------------------------------------------------
  real(DP2) function pf2d_source(x2, x3, qinj, beta, phi, Lgauss) result(src)
    implicit none
    real(DP2), intent(in) :: x2, x3, qinj, beta, phi, Lgauss

    src = qinj / (beta * phi) * (1.0_DP2 / (2.0_DP2 * PI2 * Lgauss**2)) * &
          exp(-(x2**2 + x3**2) / (2.0_DP2 * Lgauss**2))
  end function pf2d_source

  !-----------------------------------------------------------------------
  ! Injection rate at time t (Eq. 20): q0 for 0<=t<toff, 0 after.
  !-----------------------------------------------------------------------
  real(DP2) function pf2d_qinj(t, q0, toff) result(q)
    implicit none
    real(DP2), intent(in) :: t, q0, toff
    if (t >= 0.0_DP2 .and. t < toff) then
      q = q0
    else
      q = 0.0_DP2
    end if
  end function pf2d_qinj

  !-----------------------------------------------------------------------
  ! Largest forward-Euler timestep stable for the 5-point Laplacian on a
  ! grid with spacing dx and diffusivity alpha (CFL-type limit dx^2/4alpha
  ! for 2D explicit diffusion); returns a value with a small safety factor.
  !-----------------------------------------------------------------------
  real(DP2) function pf2d_stable_dt(dx, alpha) result(dt_max)
    implicit none
    real(DP2), intent(in) :: dx, alpha
    dt_max = 0.9_DP2 * dx**2 / (4.0_DP2 * alpha)
  end function pf2d_stable_dt

  !-----------------------------------------------------------------------
  ! Advance p(nx,ny) from t_start to t_start+dt_outer in-place, using
  ! internal forward-Euler sub-steps at (or below) the stable dt. Grid
  ! node (i,j), i=1..nx, j=1..ny, sits at x2 = x2min+(i-1)*dx,
  ! x3 = x3min+(j-1)*dx. Zero-flux (Neumann) boundaries via ghost-node
  ! reflection (dp/dn=0 <=> ghost value equals first interior value).
  !-----------------------------------------------------------------------
  subroutine pf2d_step(p, nx, ny, dx, x2min, x3min, t_start, dt_outer, &
                       alpha, beta, phi, q0, toff, Lgauss)
    implicit none
    integer, intent(in) :: nx, ny
    real(DP2), intent(inout) :: p(nx, ny)
    real(DP2), intent(in) :: dx, x2min, x3min, t_start, dt_outer
    real(DP2), intent(in) :: alpha, beta, phi, q0, toff, Lgauss

    real(DP2) :: dt_stable, t, dt, t_remaining, qinj
    real(DP2) :: p_new(nx, ny)
    real(DP2) :: lap, x2, x3, pl, pr, pd, pu
    integer :: i, j

    dt_stable = pf2d_stable_dt(dx, alpha)
    t = t_start
    t_remaining = dt_outer

    do while (t_remaining > 0.0_DP2)
      dt = min(dt_stable, t_remaining)
      qinj = pf2d_qinj(t + 0.5_DP2 * dt, q0, toff)  ! midpoint-in-time injection rate

      do j = 1, ny
        do i = 1, nx
          ! Neumann (zero-flux) boundary via ghost = nearest interior value
          if (i == 1) then
            pl = p(i, j)
          else
            pl = p(i-1, j)
          end if
          if (i == nx) then
            pr = p(i, j)
          else
            pr = p(i+1, j)
          end if
          if (j == 1) then
            pd = p(i, j)
          else
            pd = p(i, j-1)
          end if
          if (j == ny) then
            pu = p(i, j)
          else
            pu = p(i, j+1)
          end if

          lap = (pl + pr + pd + pu - 4.0_DP2 * p(i, j)) / dx**2

          x2 = x2min + real(i-1, DP2) * dx
          x3 = x3min + real(j-1, DP2) * dx

          p_new(i, j) = p(i, j) + dt * (alpha * lap + pf2d_source(x2, x3, qinj, beta, phi, Lgauss))
        end do
      end do

      p = p_new
      t = t + dt
      t_remaining = t_remaining - dt
    end do
  end subroutine pf2d_step

  !-----------------------------------------------------------------------
  ! Bilinear interpolation of the structured grid p(nx,ny) at (xq,yq).
  ! Clamps to the grid boundary rather than extrapolating.
  !-----------------------------------------------------------------------
  real(DP2) function pf2d_interpolate(p, nx, ny, dx, x2min, x3min, xq, yq) result(pq)
    implicit none
    integer, intent(in) :: nx, ny
    real(DP2), intent(in) :: p(nx, ny), dx, x2min, x3min, xq, yq

    real(DP2) :: fx, fy, tx, ty
    integer :: i0, j0, i1, j1

    fx = (xq - x2min) / dx
    fy = (yq - x3min) / dx

    i0 = max(1, min(nx - 1, int(floor(fx)) + 1))
    j0 = max(1, min(ny - 1, int(floor(fy)) + 1))
    i1 = i0 + 1
    j1 = j0 + 1

    tx = max(0.0_DP2, min(1.0_DP2, fx - real(i0 - 1, DP2)))
    ty = max(0.0_DP2, min(1.0_DP2, fy - real(j0 - 1, DP2)))

    pq = (1.0_DP2 - tx) * (1.0_DP2 - ty) * p(i0, j0) + &
         tx * (1.0_DP2 - ty) * p(i1, j0) + &
         (1.0_DP2 - tx) * ty * p(i0, j1) + &
         tx * ty * p(i1, j1)
  end function pf2d_interpolate

  !-----------------------------------------------------------------------
  ! Exponential integral E1(x) = integral_x^inf exp(-t)/t dt, x>0.
  ! Series (Abramowitz & Stegun 5.1.11) for x<=1, continued-fraction
  ! (5.1.56-style asymptotic rational approximation) for x>1. Used only for
  ! the analytic verification solution (Eq. 21), not in the main solver.
  !-----------------------------------------------------------------------
  real(DP2) function pf2d_expint_e1(x) result(e1)
    implicit none
    real(DP2), intent(in) :: x
    real(DP2) :: term, s, euler_gamma
    integer :: k

    euler_gamma = 0.5772156649015329_DP2

    if (x <= 0.0_DP2) then
      e1 = huge(1.0_DP2)
      return
    end if

    if (x <= 1.0_DP2) then
      ! E1(x) = -gamma - ln(x) + sum_{k=1}^inf (-1)^(k+1) x^k / (k*k!)
      s = 0.0_DP2
      term = x
      do k = 1, 40
        s = s + term / k
        term = -term * x / real(k + 1, DP2)
        if (abs(term / (k+1)) < 1.0e-16_DP2) exit
      end do
      e1 = -euler_gamma - log(x) + s
    else
      ! Continued-fraction evaluation (Lentz's algorithm), accurate to
      ! double precision for x>1.
      e1 = pf2d_expint_e1_cf(x)
    end if
  end function pf2d_expint_e1

  ! Continued-fraction evaluation of E1 for x>1 (Lentz's algorithm),
  ! accurate to double precision.
  real(DP2) function pf2d_expint_e1_cf(x) result(e1)
    implicit none
    real(DP2), intent(in) :: x
    real(DP2), parameter :: FPMIN = 1.0e-300_DP2, EPS = 1.0e-16_DP2
    real(DP2) :: a, b, c, d, h, del
    integer :: i

    b = x + 1.0_DP2
    c = 1.0_DP2 / FPMIN
    d = 1.0_DP2 / b
    h = d
    do i = 1, 200
      a = -real(i, DP2) * real(i, DP2)
      b = b + 2.0_DP2
      d = 1.0_DP2 / (a * d + b)
      c = b + a / c
      del = c * d
      h = h * del
      if (abs(del - 1.0_DP2) < EPS) exit
    end do
    e1 = exp(-x) * h
  end function pf2d_expint_e1_cf

  !-----------------------------------------------------------------------
  ! Analytic short-time solution, BP8 Eq. (21), valid for
  ! t << lf^2/(4*alpha) (~220 hours for the Table 1 parameters). For
  ! verification only.
  !-----------------------------------------------------------------------
  real(DP2) function pf2d_analytic(x2, x3, t, q0, alpha, beta, phi, Lgauss) result(p)
    implicit none
    real(DP2), intent(in) :: x2, x3, t, q0, alpha, beta, phi, Lgauss
    real(DP2) :: r2, denom1, denom2

    if (t <= 0.0_DP2) then
      p = 0.0_DP2
      return
    end if

    r2 = x2**2 + x3**2
    if (r2 < 1.0e-12_DP2) then
      p = q0 / (4.0_DP2 * PI2 * alpha * beta * phi) * &
          log((2.0_DP2 * Lgauss**2 + 4.0_DP2 * alpha * t) / (2.0_DP2 * Lgauss**2))
    else
      denom1 = 2.0_DP2 * Lgauss**2 + 4.0_DP2 * alpha * t
      denom2 = 2.0_DP2 * Lgauss**2
      p = q0 / (4.0_DP2 * PI2 * alpha * beta * phi) * &
          (pf2d_expint_e1(r2 / denom1) - pf2d_expint_e1(r2 / denom2))
    end if
  end function pf2d_analytic

end module pore_fluid_2d
