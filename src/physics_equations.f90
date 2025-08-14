!===============================================================================
! MODULE: Physics Equations
! 
! This module implements the physics equations for the earthquake simulation,
! including rate-and-state friction laws, stress calculations, and constitutive relations.
!
! Author: Refactored from original 3dtri_BP5.f90
! Last Modified: Current refactoring
!===============================================================================

module physics_equations
  use physical_constants, only: DP, PI, REFERENCE_VELOCITY, FRICTION_COEFFICIENT, VISCOSITY
  use physical_variables, only: shear_stress_1, shear_stress_2, effective_normal_stress
  use physical_variables, only: friction_parameter_a, friction_parameter_b, characteristic_length
  use physical_variables, only: state_variable_1, state_variable_2, slip_rate
  use physical_variables, only: stiffness_matrix, stiffness_matrix_2
  use simulation_parameters, only: plate_velocity
  use mpi
  implicit none
  
  ! Make all functions public
  public
  
  ! Friction law parameters
  real(DP), parameter :: AGING_LAW = 1
  real(DP), parameter :: SLIP_LAW = 2
  
  ! Current friction law (can be changed at runtime)
  integer :: current_friction_law = AGING_LAW
  
  ! Numerical parameters
  real(DP), parameter :: SMALL_VELOCITY = 1.0e-6_DP
  real(DP), parameter :: LARGE_VELOCITY = 1.0e6_DP
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Calculate derivatives for rate-and-state friction
  !===============================================================================
  subroutine calculate_friction_derivatives(myid, dydt, nv, total_elements, local_elements, &
                                          current_time, yt, z_all, x_coord)
    implicit none
    integer, intent(in) :: myid, nv, total_elements, local_elements
    real(DP), intent(in) :: current_time
    real(DP), dimension(:), intent(in) :: yt, z_all, x_coord
    real(DP), dimension(:), intent(out) :: dydt
    
    real(DP), dimension(local_elements) :: zz, zz_ds, zz_all_global, zz_ds_all_global
    real(DP), dimension(local_elements) :: zzfric, zzfric2
    real(DP), dimension(local_elements) :: sr_local, vi_local
    real(DP) :: psi, help1, help2, help, deriv1, deriv2, deriv3
    integer :: i, j, ierr, master
    real(DP) :: time_start, time_end
    
    ! Initialize
    master = 0
    
    ! Calculate local slip rates and state variables
    do i = 1, local_elements
      zz(i) = yt(2*i-1) * state_variable_1(i) - plate_velocity
      zz_ds(i) = yt(2*i-1) * state_variable_2(i)
    end do
    
    ! Gather slip rates from all processes
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Gather(zz, local_elements, MPI_DOUBLE_PRECISION, zz_all_global, &
                    local_elements, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(zz_all_global, total_elements, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    
    call MPI_Gather(zz_ds, local_elements, MPI_DOUBLE_PRECISION, zz_ds_all_global, &
                    local_elements, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(zz_ds_all_global, total_elements, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    
    ! Calculate stiffness matrix contributions
    call cpu_time(time_start)
    
    ! Calculate friction forces using stiffness matrix
    do i = 1, local_elements
      zzfric(i) = 0.0_DP
      do j = 1, total_elements
        zzfric(i) = zzfric(i) + stiffness_matrix(i, j) * zz_all_global(j)
      end do
    end do
    
    call cpu_time(time_end)
    
    ! Calculate derivatives for each element
    do i = 1, local_elements
      ! Calculate logarithmic term
      psi = log(REFERENCE_VELOCITY * yt(2*i) / characteristic_length(i))
      
      ! Calculate help terms
      help1 = yt(2*i-1) / (2.0_DP * REFERENCE_VELOCITY)
      help2 = (FRICTION_COEFFICIENT + friction_parameter_b(i) * psi) / friction_parameter_a(i)
      help = sqrt(1.0_DP + (help1 * exp(help2))**2)
      
      ! Calculate derivative components
      deriv1 = (effective_normal_stress(i) * friction_parameter_b(i) / yt(2*i)) * &
               help1 * exp(help2) / help
      deriv2 = (effective_normal_stress(i) * friction_parameter_a(i) / (2.0_DP * REFERENCE_VELOCITY)) * &
               exp(help2) / help
      
      ! Calculate state variable derivative based on friction law
      if (current_friction_law == AGING_LAW) then
        deriv3 = 1.0_DP - yt(2*i-1) * yt(2*i) / characteristic_length(i)
      else if (current_friction_law == SLIP_LAW) then
        deriv3 = -yt(2*i-1) * yt(2*i) / characteristic_length(i) * &
                 log(yt(2*i-1) * yt(2*i) / characteristic_length(i))
      else
        deriv3 = 1.0_DP - yt(2*i-1) * yt(2*i) / characteristic_length(i)  ! Default to aging law
      end if
      
      ! Calculate final derivatives
      dydt(2*i-1) = -(zzfric(i) + deriv1 * deriv3) / (VISCOSITY + deriv2)
      dydt(2*i) = deriv3
    end do
    
  end subroutine calculate_friction_derivatives
  
  !===============================================================================
  ! SUBROUTINE: Calculate depth-dependent parameters
  !===============================================================================
  subroutine calculate_depth_dependent_parameters(total_elements, dip_angle, nucleation_height, &
                                                sigma_difference, effective_stress_fixed, &
                                                length_fixed, transition1, transition2, &
                                                slip1, slip2, perturbation_flag, factor1, &
                                                factor2, factor3, factor4, xi_lock1, xi_lock2, &
                                                friction_parameter_a, friction_parameter_b, &
                                                characteristic_length, effective_normal_stress, &
                                                x_all, z_all, initial_velocity)
    implicit none
    integer, intent(in) :: total_elements, perturbation_flag
    real(DP), intent(in) :: dip_angle, nucleation_height, sigma_difference
    real(DP), intent(in) :: effective_stress_fixed, length_fixed
    integer, dimension(:), intent(in) :: transition1, transition2, slip1, slip2
    real(DP), intent(in) :: factor1, factor2, factor3, factor4
    real(DP), intent(in) :: xi_lock1, xi_lock2
    real(DP), dimension(:), intent(out) :: friction_parameter_a, friction_parameter_b
    real(DP), dimension(:), intent(out) :: characteristic_length, effective_normal_stress
    real(DP), dimension(:), intent(in) :: x_all, z_all
    real(DP), dimension(:), intent(out) :: initial_velocity
    
    real(DP) :: depth, effective_stress, characteristic_slip
    integer :: i
    
    ! Calculate parameters for each element based on depth
    do i = 1, total_elements
      ! Calculate depth from coordinates
      depth = abs(z_all(i)) * cos(dip_angle)
      
      ! Calculate effective normal stress (depth-dependent)
      effective_stress = effective_stress_fixed + sigma_difference * depth
      effective_normal_stress(i) = effective_stress
      
      ! Calculate characteristic slip distance (depth-dependent)
      characteristic_slip = length_fixed * (1.0_DP + factor1 * depth)
      characteristic_length(i) = characteristic_slip
      
      ! Calculate friction parameters (depth-dependent)
      friction_parameter_a(i) = factor2 * (1.0_DP + factor3 * depth)
      friction_parameter_b(i) = factor4 * (1.0_DP + factor3 * depth)
      
      ! Set initial velocity
      initial_velocity(i) = 0.0_DP
    end do
    
    ! Apply perturbations if requested
    if (perturbation_flag > 0) then
      call apply_perturbations(total_elements, transition1, transition2, slip1, slip2, &
                              xi_lock1, xi_lock2, friction_parameter_a, friction_parameter_b, &
                              characteristic_length, effective_normal_stress, x_all)
    end if
    
  end subroutine calculate_depth_dependent_parameters
  
  !===============================================================================
  ! SUBROUTINE: Apply perturbations to parameters
  !===============================================================================
  subroutine apply_perturbations(total_elements, transition1, transition2, slip1, slip2, &
                                xi_lock1, xi_lock2, friction_parameter_a, friction_parameter_b, &
                                characteristic_length, effective_normal_stress, x_all)
    implicit none
    integer, intent(in) :: total_elements
    integer, dimension(:), intent(in) :: transition1, transition2, slip1, slip2
    real(DP), intent(in) :: xi_lock1, xi_lock2
    real(DP), dimension(:), intent(in) :: x_all
    real(DP), dimension(:), intent(inout) :: friction_parameter_a, friction_parameter_b
    real(DP), dimension(:), intent(inout) :: characteristic_length, effective_normal_stress
    
    integer :: i, j
    
    ! Apply perturbations to specific elements
    do i = 1, total_elements
      ! Check if element is in transition zone 1
      do j = 1, size(transition1)
        if (i == transition1(j)) then
          friction_parameter_a(i) = friction_parameter_a(i) * 1.1_DP  ! Increase by 10%
          friction_parameter_b(i) = friction_parameter_b(i) * 1.1_DP
        end if
      end do
      
      ! Check if element is in transition zone 2
      do j = 1, size(transition2)
        if (i == transition2(j)) then
          friction_parameter_a(i) = friction_parameter_a(i) * 0.9_DP  ! Decrease by 10%
          friction_parameter_b(i) = friction_parameter_b(i) * 0.9_DP
        end if
      end do
      
      ! Check if element is in slip zone 1
      do j = 1, size(slip1)
        if (i == slip1(j)) then
          characteristic_length(i) = characteristic_length(i) * 1.2_DP  ! Increase by 20%
        end if
      end do
      
      ! Check if element is in slip zone 2
      do j = 1, size(slip2)
        if (i == slip2(j)) then
          characteristic_length(i) = characteristic_length(i) * 0.8_DP  ! Decrease by 20%
        end if
      end do
    end do
    
    ! Apply locked zone perturbations
    do i = 1, total_elements
      if (abs(x_all(i)) >= xi_lock1 .and. abs(x_all(i)) <= xi_lock2) then
        ! Locked zone: reduce friction parameters
        friction_parameter_a(i) = friction_parameter_a(i) * 0.5_DP
        friction_parameter_b(i) = friction_parameter_b(i) * 0.5_DP
        characteristic_length(i) = characteristic_length(i) * 0.1_DP
      end if
    end do
    
  end subroutine apply_perturbations
  
  !===============================================================================
  ! SUBROUTINE: Calculate stress from slip using stiffness matrix
  !===============================================================================
  subroutine calculate_stress_from_slip(local_elements, total_elements, slip_increment, &
                                       stiffness_matrix, stress_increment)
    implicit none
    integer, intent(in) :: local_elements, total_elements
    real(DP), dimension(:), intent(in) :: slip_increment
    real(DP), dimension(:,:), intent(in) :: stiffness_matrix
    real(DP), dimension(:), intent(out) :: stress_increment
    
    integer :: i, j
    
    ! Calculate stress increment using stiffness matrix
    do i = 1, local_elements
      stress_increment(i) = 0.0_DP
      do j = 1, total_elements
        stress_increment(i) = stress_increment(i) + stiffness_matrix(i, j) * slip_increment(j)
      end do
    end do
    
  end subroutine calculate_stress_from_slip
  
  !===============================================================================
  ! SUBROUTINE: Calculate slip from stress using inverse stiffness
  !===============================================================================
  subroutine calculate_slip_from_stress(local_elements, total_elements, stress_increment, &
                                       stiffness_matrix, slip_increment)
    implicit none
    integer, intent(in) :: local_elements, total_elements
    real(DP), dimension(:), intent(in) :: stress_increment
    real(DP), dimension(:,:), intent(in) :: stiffness_matrix
    real(DP), dimension(:), intent(out) :: slip_increment
    
    ! This is a simplified implementation - in practice, you would need to
    ! solve the linear system or use an iterative method
    
    ! For now, just use a simple approximation
    slip_increment = stress_increment / maxval(stiffness_matrix)
    
  end subroutine calculate_slip_from_stress
  
  !===============================================================================
  ! SUBROUTINE: Update state variables
  !===============================================================================
  subroutine update_state_variables(local_elements, current_time, dt, yt, dydt, yt_new)
    implicit none
    integer, intent(in) :: local_elements
    real(DP), intent(in) :: current_time, dt
    real(DP), dimension(:), intent(in) :: yt, dydt
    real(DP), dimension(:), intent(out) :: yt_new
    
    integer :: i
    
    ! Update state variables using Euler method (can be replaced with RK4)
    do i = 1, 2*local_elements
      yt_new(i) = yt(i) + dt * dydt(i)
    end do
    
  end subroutine update_state_variables
  
  !===============================================================================
  ! SUBROUTINE: Check physical constraints
  !===============================================================================
  subroutine check_physical_constraints(local_elements, yt, dydt, is_valid)
    implicit none
    integer, intent(in) :: local_elements
    real(DP), dimension(:), intent(in) :: yt, dydt
    logical, intent(out) :: is_valid
    
    integer :: i
    logical :: constraints_violated
    
    constraints_violated = .false.
    
    ! Check for negative slip rates
    do i = 1, local_elements
      if (yt(2*i-1) < -LARGE_VELOCITY) then
        constraints_violated = .true.
        write(*, *) 'Warning: Negative slip rate detected at element ', i
      end if
    end do
    
    ! Check for negative state variables
    do i = 1, local_elements
      if (yt(2*i) < 0.0_DP) then
        constraints_violated = .false.
        write(*, *) 'Warning: Negative state variable detected at element ', i
      end if
    end do
    
    ! Check for excessive derivatives
    do i = 1, 2*local_elements
      if (abs(dydt(i)) > LARGE_VELOCITY) then
        constraints_violated = .true.
        write(*, *) 'Warning: Excessive derivative detected at component ', i
      end if
    end do
    
    is_valid = .not. constraints_violated
    
  end subroutine check_physical_constraints
  
  !===============================================================================
  ! SUBROUTINE: Set friction law
  !===============================================================================
  subroutine set_friction_law(law_type)
    implicit none
    integer, intent(in) :: law_type
    
    if (law_type == AGING_LAW .or. law_type == SLIP_LAW) then
      current_friction_law = law_type
      write(*, *) 'Friction law set to: ', current_friction_law
    else
      write(*, *) 'Warning: Invalid friction law type, keeping current law'
    end if
    
  end subroutine set_friction_law
  
  !===============================================================================
  ! SUBROUTINE: Get current friction law
  !===============================================================================
  function get_current_friction_law() result(law_type)
    implicit none
    integer :: law_type
    
    law_type = current_friction_law
    
  end function get_current_friction_law
  
  !===============================================================================
  ! SUBROUTINE: Print physics parameters
  !===============================================================================
  subroutine print_physics_parameters(local_elements)
    implicit none
    integer, intent(in) :: local_elements
    
    integer :: i
    
    write(*, *) '=== Physics Parameters ==='
    write(*, *) 'Friction law: ', current_friction_law
    write(*, *) 'Number of elements: ', local_elements
    write(*, *)
    
    write(*, *) 'Sample element parameters:'
    do i = 1, min(5, local_elements)
      write(*, *) 'Element ', i, ':'
      write(*, *) '  a = ', friction_parameter_a(i)
      write(*, *) '  b = ', friction_parameter_b(i)
      write(*, *) '  Dc = ', characteristic_length(i)
      write(*, *) '  sigma_eff = ', effective_normal_stress(i)
      write(*, *)
    end do
    write(*, *) '========================='
    
  end subroutine print_physics_parameters
  
end module physics_equations
