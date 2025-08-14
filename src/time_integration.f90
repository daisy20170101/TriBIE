!===============================================================================
! MODULE: Time Integration
! 
! This module provides time integration methods for the earthquake simulation,
! including Runge-Kutta solvers with error control and adaptive time stepping.
!
! Author: Refactored from original 3dtri_BP5.f90
! Last Modified: Current refactoring
!===============================================================================

module time_integration
  use physical_constants, only: DP
  implicit none
  
  ! Make all functions public
  public
  
  ! Integration method parameters
  integer, parameter :: RK4_ORDER = 4
  integer, parameter :: RK5_ORDER = 5
  integer, parameter :: MAX_STEPS = 1000000
  
  ! Error control parameters
  real(DP), parameter :: DEFAULT_TOLERANCE = 1.0e-6_DP
  real(DP), parameter :: MIN_SCALE_FACTOR = 0.1_DP
  real(DP), parameter :: MAX_SCALE_FACTOR = 10.0_DP
  real(DP), parameter :: SAFETY_FACTOR = 0.9_DP
  
  ! Runge-Kutta coefficients (Cash-Karp method)
  real(DP), parameter :: A2 = 0.2_DP
  real(DP), parameter :: A3 = 0.3_DP
  real(DP), parameter :: A4 = 0.6_DP
  real(DP), parameter :: A5 = 1.0_DP
  real(DP), parameter :: A6 = 0.875_DP
  
  real(DP), parameter :: B21 = 0.2_DP
  real(DP), parameter :: B31 = 3.0_DP/40.0_DP
  real(DP), parameter :: B32 = 9.0_DP/40.0_DP
  real(DP), parameter :: B41 = 0.3_DP
  real(DP), parameter :: B42 = -0.9_DP
  real(DP), parameter :: B43 = 1.2_DP
  real(DP), parameter :: B51 = -11.0_DP/54.0_DP
  real(DP), parameter :: B52 = 2.5_DP
  real(DP), parameter :: B53 = -70.0_DP/27.0_DP
  real(DP), parameter :: B54 = 35.0_DP/27.0_DP
  real(DP), parameter :: B61 = 1631.0_DP/55296.0_DP
  real(DP), parameter :: B62 = 175.0_DP/512.0_DP
  real(DP), parameter :: B63 = 575.0_DP/13824.0_DP
  real(DP), parameter :: B64 = 44275.0_DP/110592.0_DP
  real(DP), parameter :: B65 = 253.0_DP/4096.0_DP
  
  ! Fifth order weights
  real(DP), parameter :: C1 = 37.0_DP/378.0_DP
  real(DP), parameter :: C3 = 250.0_DP/621.0_DP
  real(DP), parameter :: C4 = 125.0_DP/594.0_DP
  real(DP), parameter :: C6 = 512.0_DP/1771.0_DP
  
  ! Fourth order weights
  real(DP), parameter :: D1 = 2825.0_DP/27648.0_DP
  real(DP), parameter :: D3 = 18575.0_DP/48384.0_DP
  real(DP), parameter :: D4 = 13525.0_DP/55296.0_DP
  real(DP), parameter :: D5 = 277.0_DP/14336.0_DP
  real(DP), parameter :: D6 = 0.25_DP
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Runge-Kutta step with error control
  !===============================================================================
  subroutine rkqs_with_error_control(y, dydx, n, t, htry, hdid, hnext)
    implicit none
    real(DP), dimension(:), intent(inout) :: y
    real(DP), dimension(:), intent(in) :: dydx
    integer, intent(in) :: n
    real(DP), intent(in) :: t, htry
    real(DP), intent(out) :: hdid, hnext
    
    real(DP), dimension(n) :: yerr, ytemp
    real(DP) :: h, h_temp, errmax, h_scale
    integer :: i
    
    ! Initialize
    h = htry
    hdid = 0.0_DP
    hnext = 0.0_DP
    
    ! Main integration loop with error control
    do
      ! Take a step
      call rkck_step(dydx, h, n, y, ytemp, yerr)
      
      ! Calculate error
      errmax = 0.0_DP
      do i = 1, n
        errmax = max(errmax, abs(yerr(i) / max(abs(ytemp(i)), 1.0_DP)))
      end do
      
      ! Scale error by tolerance
      errmax = errmax / DEFAULT_TOLERANCE
      
      ! If error is acceptable, accept step
      if (errmax <= 1.0_DP) then
        hdid = h
        y = ytemp
        
        ! Calculate next step size
        if (errmax > 1.89e-4_DP) then
          h_scale = SAFETY_FACTOR * (errmax ** 0.2_DP)
        else
          h_scale = 5.0_DP
        end if
        
        hnext = h * min(MAX_SCALE_FACTOR, max(MIN_SCALE_FACTOR, h_scale))
        exit
      else
        ! Reduce step size and try again
        h_temp = h * SAFETY_FACTOR * (errmax ** 0.25_DP)
        h = sign(max(abs(h_temp), 0.1_DP * abs(h)), h)
        
        ! Check if step size is too small
        if (abs(h) < 1.0e-12_DP) then
          call handle_integration_error('Step size too small in RKQS')
        end if
      end if
    end do
    
  end subroutine rkqs_with_error_control
  
  !===============================================================================
  ! SUBROUTINE: Runge-Kutta Cash-Karp step
  !===============================================================================
  subroutine rkck_step(dydx, h, n, y, yout, yerr)
    implicit none
    real(DP), dimension(:), intent(in) :: dydx
    real(DP), intent(in) :: h
    integer, intent(in) :: n
    real(DP), dimension(:), intent(in) :: y
    real(DP), dimension(:), intent(out) :: yout, yerr
    
    real(DP), dimension(n) :: ak2, ak3, ak4, ak5, ak6
    real(DP), dimension(n) :: ytemp
    integer :: i
    
    ! First step
    do i = 1, n
      ytemp(i) = y(i) + B21 * h * dydx(i)
    end do
    
    ! Second step
    call evaluate_derivatives(ytemp, ak2)
    do i = 1, n
      ytemp(i) = y(i) + h * (B31 * dydx(i) + B32 * ak2(i))
    end do
    
    ! Third step
    call evaluate_derivatives(ytemp, ak3)
    do i = 1, n
      ytemp(i) = y(i) + h * (B41 * dydx(i) + B42 * ak2(i) + B43 * ak3(i))
    end do
    
    ! Fourth step
    call evaluate_derivatives(ytemp, ak4)
    do i = 1, n
      ytemp(i) = y(i) + h * (B51 * dydx(i) + B52 * ak2(i) + B53 * ak3(i) + B54 * ak4(i))
    end do
    
    ! Fifth step
    call evaluate_derivatives(ytemp, ak5)
    do i = 1, n
      ytemp(i) = y(i) + h * (B61 * dydx(i) + B62 * ak2(i) + B63 * ak3(i) + B64 * ak4(i) + B65 * ak5(i))
    end do
    
    ! Sixth step
    call evaluate_derivatives(ytemp, ak6)
    
    ! Calculate fifth order solution
    do i = 1, n
      yout(i) = y(i) + h * (C1 * dydx(i) + C3 * ak3(i) + C4 * ak4(i) + C6 * ak6(i))
    end do
    
    ! Calculate fourth order solution for error estimate
    do i = 1, n
      yerr(i) = h * ((C1 - D1) * dydx(i) + (C3 - D3) * ak3(i) + &
                      (C4 - D4) * ak4(i) + (C6 - D6) * ak6(i))
    end do
    
  end subroutine rkck_step
  
  !===============================================================================
  ! SUBROUTINE: Evaluate derivatives (placeholder - should be implemented by user)
  !===============================================================================
  subroutine evaluate_derivatives(y, dydt)
    implicit none
    real(DP), dimension(:), intent(in) :: y
    real(DP), dimension(:), intent(out) :: dydt
    
    ! This is a placeholder - the actual implementation should be provided
    ! by the user based on their specific physics equations
    dydt = 0.0_DP
    
    ! Example implementation for testing:
    ! dydt(1) = -y(2)
    ! dydt(2) = y(1)
    
  end subroutine evaluate_derivatives
  
  !===============================================================================
  ! SUBROUTINE: Simple RK4 integration (for comparison)
  !===============================================================================
  subroutine rk4_simple(y, dydx, n, h, yout)
    implicit none
    real(DP), dimension(:), intent(in) :: y, dydx
    integer, intent(in) :: n
    real(DP), intent(in) :: h
    real(DP), dimension(:), intent(out) :: yout
    
    real(DP), dimension(n) :: k1, k2, k3, k4, ytemp
    integer :: i
    
    ! First step
    k1 = h * dydx
    
    ! Second step
    do i = 1, n
      ytemp(i) = y(i) + 0.5_DP * k1(i)
    end do
    call evaluate_derivatives(ytemp, k2)
    k2 = h * k2
    
    ! Third step
    do i = 1, n
      ytemp(i) = y(i) + 0.5_DP * k2(i)
    end do
    call evaluate_derivatives(ytemp, k3)
    k3 = h * k3
    
    ! Fourth step
    do i = 1, n
      ytemp(i) = y(i) + k3(i)
    end do
    call evaluate_derivatives(ytemp, k4)
    k4 = h * k4
    
    ! Final result
    do i = 1, n
      yout(i) = y(i) + (k1(i) + 2.0_DP * k2(i) + 2.0_DP * k3(i) + k4(i)) / 6.0_DP
    end do
    
  end subroutine rk4_simple
  
  !===============================================================================
  ! SUBROUTINE: Adaptive time stepping
  !===============================================================================
  subroutine adaptive_time_step(current_time, current_dt, target_dt, min_dt, max_dt, dt_new)
    implicit none
    real(DP), intent(in) :: current_time, current_dt, target_dt, min_dt, max_dt
    real(DP), intent(out) :: dt_new
    
    real(DP) :: dt_candidate
    
    ! Calculate candidate time step
    dt_candidate = target_dt
    
    ! Ensure time step is within bounds
    dt_candidate = max(min_dt, min(max_dt, dt_candidate))
    
    ! Ensure we don't exceed maximum time
    if (current_time + dt_candidate > current_time + max_dt) then
      dt_candidate = max_dt
    end if
    
    ! Ensure we don't go below minimum time step
    if (dt_candidate < min_dt) then
      dt_candidate = min_dt
    end if
    
    dt_new = dt_candidate
    
  end subroutine adaptive_time_step
  
  !===============================================================================
  ! SUBROUTINE: Check integration stability
  !===============================================================================
  subroutine check_integration_stability(y, y_old, tolerance, is_stable)
    implicit none
    real(DP), dimension(:), intent(in) :: y, y_old
    real(DP), intent(in) :: tolerance
    logical, intent(out) :: is_stable
    
    real(DP) :: max_change
    integer :: i, n
    
    n = size(y)
    max_change = 0.0_DP
    
    ! Calculate maximum relative change
    do i = 1, n
      if (abs(y_old(i)) > 1.0e-12_DP) then
        max_change = max(max_change, abs((y(i) - y_old(i)) / y_old(i)))
      end if
    end do
    
    ! Check if change is within tolerance
    is_stable = (max_change <= tolerance)
    
  end subroutine check_integration_stability
  
  !===============================================================================
  ! SUBROUTINE: Handle integration errors
  !===============================================================================
  subroutine handle_integration_error(error_message)
    implicit none
    character(len=*), intent(in) :: error_message
    
    write(*, *) 'Integration error: ', trim(error_message)
    stop
    
  end subroutine handle_integration_error
  
  !===============================================================================
  ! SUBROUTINE: Print integration statistics
  !===============================================================================
  subroutine print_integration_statistics(total_steps, successful_steps, failed_steps, total_time)
    implicit none
    integer, intent(in) :: total_steps, successful_steps, failed_steps
    real(DP), intent(in) :: total_time
    
    write(*, *) '=== Integration Statistics ==='
    write(*, *) 'Total steps attempted: ', total_steps
    write(*, *) 'Successful steps: ', successful_steps
    write(*, *) 'Failed steps: ', failed_steps
    write(*, *) 'Success rate: ', real(successful_steps, DP) / real(total_steps, DP) * 100.0_DP, '%'
    write(*, *) 'Total integration time: ', total_time, ' seconds'
    write(*, *) 'Average time per step: ', total_time / real(successful_steps, DP), ' seconds'
    write(*, *) '================================'
    
  end subroutine print_integration_statistics
  
end module time_integration
