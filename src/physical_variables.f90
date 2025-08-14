!===============================================================================
! MODULE: Physical Variables
! 
! This module defines all physical variables used in the earthquake simulation
! including stress, slip, friction parameters, and state variables.
!
! Author: Refactored from original phy3d_module_non.f90
! Last Modified: Current refactoring
!===============================================================================

module physical_variables
  use physical_constants, only: DP
  implicit none
  
  ! Make all variables public
  public
  
  ! Physical variables (stress, slip, friction parameters)
  real(DP), dimension(:), allocatable, target :: shear_stress_1, shear_stress_2, total_shear_stress
  real(DP), dimension(:), allocatable, target :: friction_parameter_a, friction_parameter_b, characteristic_length
  real(DP), dimension(:), allocatable, target :: effective_normal_stress, state_variable_1, state_variable_2
  real(DP), dimension(:), allocatable, target :: slip_rate, initial_velocity
  
  ! Slip variables
  real(DP), dimension(:), allocatable, target :: total_slip, slip_increment, dip_slip, dip_slip_increment
  
  ! Stiffness matrices
  real(DP), dimension(:,:), allocatable, target :: stiffness_matrix, stiffness_matrix_2
  
  ! Time variables
  real(DP), target :: time_start, time_end, time_day, time_else, time_midnight, time_multiply
  
  ! Physical parameters (from original code, kept for compatibility)
  real(DP), parameter :: pi = 3.14159265358979323_DP
  real(DP), parameter :: amax = 0.025_DP
  real(DP), parameter :: phi = 60.0_DP / 180.0_DP * pi
  real(DP), parameter :: xmu = 0.32038_DP
  real(DP), parameter :: cs = 10.92407e7_DP
  real(DP), parameter :: xnu = 0.25_DP
  real(DP), parameter :: V0 = 3.15e4_DP
  real(DP), parameter :: f0 = 0.6_DP
  real(DP), parameter :: eta = 0.5_DP * xmu / cs
  real(DP), parameter :: gamma = 2.0_DP / pi
  real(DP), parameter :: p18 = 2.0_DP * pi / 360.0_DP
  real(DP), parameter :: yrs = 365.0_DP * 24.0_DP * 3600.0_DP
  real(DP), parameter :: yrd = 365.0_DP
  
  ! Legacy variable names for compatibility (mapped to new names)
  ! Note: stiff and stiff2 are 2D pointers to match stiffness matrices
  real(DP), dimension(:), pointer :: tau1, tau2, tau0, cca, ccb, seff, xLf, phy1, phy2
  real(DP), dimension(:,:), pointer :: stiff, stiff2
  real(DP), dimension(:), pointer :: slip, slipinc, slipds, slipdsinc, sr, vi
  real(DP), dimension(:), pointer :: yt, yt0, dydt, yt_scale
  real(DP), pointer :: tm1, tm2, tmday, tmelse, tmmidn, tmmult, Vpl
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Initialize physical variables
  !===============================================================================
  subroutine initialize_physical_variables(n_elements, nl)
    implicit none
    integer, intent(in) :: n_elements, nl
    integer :: alloc_stat
    
    ! Allocate arrays
    allocate(shear_stress_1(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('shear_stress_1', alloc_stat)
    
    allocate(shear_stress_2(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('shear_stress_2', alloc_stat)
    
    allocate(total_shear_stress(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('total_shear_stress', alloc_stat)
    
    allocate(effective_normal_stress(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('effective_normal_stress', alloc_stat)
    
    allocate(friction_parameter_a(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('friction_parameter_a', alloc_stat)
    
    allocate(friction_parameter_b(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('friction_parameter_b', alloc_stat)
    
    allocate(characteristic_length(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('characteristic_length', alloc_stat)
    
    allocate(total_slip(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('total_slip', alloc_stat)
    
    allocate(slip_increment(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('slip_increment', alloc_stat)
    
    allocate(dip_slip(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('dip_slip', alloc_stat)
    
    allocate(dip_slip_increment(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('dip_slip_increment', alloc_stat)
    
    allocate(state_variable_1(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('state_variable_1', alloc_stat)
    
    allocate(state_variable_2(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('state_variable_2', alloc_stat)
    
    allocate(slip_rate(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('slip_rate', alloc_stat)
    
    allocate(initial_velocity(n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('initial_velocity', alloc_stat)
    
    ! Allocate stiffness matrices (for legacy compatibility)
    allocate(stiffness_matrix(n_elements, n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('stiffness_matrix', alloc_stat)
    
    allocate(stiffness_matrix_2(n_elements, n_elements), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('stiffness_matrix_2', alloc_stat)
    
    ! Initialize arrays to zero
    call initialize_arrays_to_zero()
    
    ! Initialize stiffness matrices to zero
    stiffness_matrix = 0.0_DP
    stiffness_matrix_2 = 0.0_DP
    
    ! Set up legacy variable mappings
    call setup_legacy_mappings()
    
  end subroutine initialize_physical_variables
  
  !===============================================================================
  ! SUBROUTINE: Initialize arrays to zero
  !===============================================================================
  subroutine initialize_arrays_to_zero()
    implicit none
    
    shear_stress_1 = 0.0_DP
    shear_stress_2 = 0.0_DP
    total_shear_stress = 0.0_DP
    effective_normal_stress = 0.0_DP
    friction_parameter_a = 0.0_DP
    friction_parameter_b = 0.0_DP
    characteristic_length = 0.0_DP
    total_slip = 0.0_DP
    slip_increment = 0.0_DP
    dip_slip = 0.0_DP
    dip_slip_increment = 0.0_DP
    state_variable_1 = 0.0_DP
    state_variable_2 = 0.0_DP
    slip_rate = 0.0_DP
    initial_velocity = 0.0_DP
    
  end subroutine initialize_arrays_to_zero
  
  !===============================================================================
  ! SUBROUTINE: Setup legacy variable mappings
  !===============================================================================
  subroutine setup_legacy_mappings()
    implicit none
    
    ! Map new variable names to legacy names for compatibility
    tau1 => shear_stress_1
    tau2 => shear_stress_2
    tau0 => total_shear_stress
    cca => friction_parameter_a
    ccb => friction_parameter_b
    seff => effective_normal_stress
    xLf => characteristic_length
    phy1 => state_variable_1
    phy2 => state_variable_2
    
    stiff => stiffness_matrix
    stiff2 => stiffness_matrix_2
    
    slip => total_slip
    slipinc => slip_increment
    slipds => dip_slip
    slipdsinc => dip_slip_increment
    sr => slip_rate
    vi => initial_velocity
    
    yt => state_variable_1  ! Note: yt is 2*n_elements in original
    yt0 => state_variable_2 ! Note: yt0 is 2*n_elements in original
    dydt => state_variable_1 ! Note: dydt is 2*n_elements in original
    yt_scale => state_variable_2 ! Note: yt_scale is 2*n_elements in original
    
    ! Map time variables
    tm1 => time_start
    tm2 => time_end
    tmday => time_day
    tmelse => time_else
    tmmidn => time_midnight
    tmmult => time_multiply
    ! Note: Vpl should be set from simulation parameters when available
    ! For now, we'll use a local variable
    Vpl => time_start  ! Placeholder - should be updated when plate_velocity is available
    
  end subroutine setup_legacy_mappings
  
  !===============================================================================
  ! SUBROUTINE: Allocate stiffness matrices
  !===============================================================================
  subroutine allocate_stiffness_matrices(n_local, n_total)
    implicit none
    integer, intent(in) :: n_local, n_total
    integer :: alloc_stat
    
    allocate(stiffness_matrix(n_local, n_total), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('stiffness_matrix', alloc_stat)
    
    allocate(stiffness_matrix_2(n_local, n_total), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_allocation_error('stiffness_matrix_2', alloc_stat)
    
    ! Initialize to zero
    stiffness_matrix = 0.0_DP
    stiffness_matrix_2 = 0.0_DP
    
  end subroutine allocate_stiffness_matrices
  
  !===============================================================================
  ! SUBROUTINE: Deallocate physical variables
  !===============================================================================
  subroutine deallocate_physical_variables()
    implicit none
    
    ! Deallocate all arrays
    if (allocated(shear_stress_1)) deallocate(shear_stress_1)
    if (allocated(shear_stress_2)) deallocate(shear_stress_2)
    if (allocated(total_shear_stress)) deallocate(total_shear_stress)
    if (allocated(effective_normal_stress)) deallocate(effective_normal_stress)
    if (allocated(friction_parameter_a)) deallocate(friction_parameter_a)
    if (allocated(friction_parameter_b)) deallocate(friction_parameter_b)
    if (allocated(characteristic_length)) deallocate(characteristic_length)
    if (allocated(total_slip)) deallocate(total_slip)
    if (allocated(slip_increment)) deallocate(slip_increment)
    if (allocated(dip_slip)) deallocate(dip_slip)
    if (allocated(dip_slip_increment)) deallocate(dip_slip_increment)
    if (allocated(state_variable_1)) deallocate(state_variable_1)
    if (allocated(state_variable_2)) deallocate(state_variable_2)
    if (allocated(slip_rate)) deallocate(slip_rate)
    if (allocated(initial_velocity)) deallocate(initial_velocity)
    if (allocated(stiffness_matrix)) deallocate(stiffness_matrix)
    if (allocated(stiffness_matrix_2)) deallocate(stiffness_matrix_2)
    
  end subroutine deallocate_physical_variables
  
  !===============================================================================
  ! SUBROUTINE: Handle allocation errors
  !===============================================================================
  subroutine handle_allocation_error(array_name, stat)
    implicit none
    character(len=*), intent(in) :: array_name
    integer, intent(in) :: stat
    
    write(*, *) 'Error: Failed to allocate array ', trim(array_name), ' with status ', stat
    stop
    
  end subroutine handle_allocation_error
  
  !===============================================================================
  ! SUBROUTINE: Initialize timing variables
  !===============================================================================
  subroutine initialize_timing()
    implicit none
    
    time_start = 0.0_DP
    time_end = 0.0_DP
    time_day = 0.0_DP
    time_else = 0.0_DP
    time_midnight = 0.0_DP
    time_multiply = 0.0_DP
    
  end subroutine initialize_timing
  
  !===============================================================================
  ! SUBROUTINE: Start timing measurement
  !===============================================================================
  subroutine start_timing()
    implicit none
    
    call cpu_time(time_start)
    
  end subroutine start_timing
  
  !===============================================================================
  ! SUBROUTINE: End timing measurement
  !===============================================================================
  subroutine end_timing()
    implicit none
    
    call cpu_time(time_end)
    
  end subroutine end_timing
  
  !===============================================================================
  ! SUBROUTINE: Get elapsed time
  !===============================================================================
  function get_elapsed_time() result(elapsed)
    implicit none
    real(DP) :: elapsed
    
    elapsed = time_end - time_start
    
  end function get_elapsed_time
  
end module physical_variables
