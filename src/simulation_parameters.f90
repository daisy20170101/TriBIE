!===============================================================================
! MODULE: Simulation Parameters
! 
! This module defines all simulation parameters including grid settings,
! time parameters, and numerical settings for the earthquake simulation.
!
! Author: Refactored from original phy3d_module_non.f90
! Last Modified: Current refactoring
!===============================================================================

module simulation_parameters
  use physical_constants, only: DP
  implicit none
  
  ! Make all parameters public
  public
  
  ! Grid and domain parameters
  integer :: total_elements          ! Total number of elements in the domain
  integer :: local_elements          ! Number of elements per process
  integer :: depth_layers            ! Number of depth layers
  integer :: length_ratio            ! Length ratio parameter
  integer :: number_processes        ! Number of MPI processes
  
  ! Time integration parameters
  real(DP) :: maximum_time          ! Maximum simulation time
  real(DP) :: time_step_min         ! Minimum time step
  real(DP) :: time_step_max         ! Maximum time step
  real(DP) :: time_step_current     ! Current time step
  
  ! Output and restart parameters
  integer :: output_interval        ! Output frequency
  integer :: restart_interval       ! Restart file frequency
  integer :: snapshot_interval      ! Snapshot frequency
  integer :: profile_flag           ! Profile output flag
  
  ! Observation and monitoring parameters
  integer :: number_observations    ! Number of observation points
  integer :: profile_strike_points  ! Number of strike profile points
  integer :: profile_dip_points     ! Number of dip profile points
  
  ! Event counting parameters
  integer :: max_velocity_events    ! Maximum velocity events to track
  integer :: aseismic_slip_events   ! Aseismic slip events to track
  integer :: coseismic_events       ! Coseismic events to track
  integer :: null_events            ! Null events to track
  integer :: slow_slip_events       ! Slow slip events to track
  integer :: null_event_interval    ! Null event interval
  
  ! Plate motion parameters
  real(DP) :: plate_velocity        ! Plate convergence velocity
  
  ! Slip timing parameters
  real(DP) :: slip_average_time     ! Time for slip averaging
  real(DP) :: slip_end_time         ! End time for slip calculations
  real(DP) :: slip_average_interval ! Interval for slip averaging
  
  ! Output timing parameters
  real(DP) :: output_time_interval  ! Time interval between outputs
  real(DP) :: output_time_min       ! Minimum time between outputs
  real(DP) :: coseismic_time_interval ! Time interval for coseismic output
  real(DP) :: slow_slip_time_interval ! Time interval for slow slip output
  
  ! Velocity thresholds
  real(DP) :: coseismic_velocity_threshold    ! Threshold for coseismic events
  real(DP) :: slow_slip_velocity_threshold1   ! First threshold for slow slip
  real(DP) :: slow_slip_velocity_threshold2   ! Second threshold for slow slip
  
  ! Output flags (10-element array from original code)
  integer, dimension(10) :: output_flags
  
  ! MPI process identification
  integer :: master_process_id = 0  ! ID of master process
  
  ! Numerical parameters
  real(DP) :: integration_tolerance = 1.0e-6_DP  ! Integration tolerance
  real(DP) :: convergence_tolerance = 1.0e-8_DP  ! Convergence tolerance
  
  ! Buffer zone parameters
  real(DP) :: buffer_zone_width = 20.0_DP  ! Buffer zone width in km
  
  ! Seismogenic zone parameters
  logical :: locked_seismogenic_zone = .true.  ! Whether seismogenic zone is locked
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Initialize default parameters
  !===============================================================================
  subroutine initialize_default_parameters()
    implicit none
    
    ! Set default values
    total_elements = 1000
    local_elements = 100
    depth_layers = 10
    length_ratio = 1
    number_processes = 1
    
    maximum_time = 100.0_DP
    time_step_min = 1.0e-6_DP
    time_step_max = 1.0e-3_DP
    time_step_current = time_step_min
    
    output_interval = 100
    restart_interval = 1000
    snapshot_interval = 500
    profile_flag = 0
    
    number_observations = 10
    profile_strike_points = 5
    profile_dip_points = 5
    
    max_velocity_events = 100
    aseismic_slip_events = 50
    coseismic_events = 50
    null_events = 50
    slow_slip_events = 50
    null_event_interval = 10
    
    plate_velocity = 3.15e4_DP
    
    slip_average_time = 1.0_DP
    slip_end_time = 10.0_DP
    slip_average_interval = 0.1_DP
    
    output_time_interval = 1.0_DP
    output_time_min = 0.1_DP
    coseismic_time_interval = 0.01_DP
    slow_slip_time_interval = 0.1_DP
    
    coseismic_velocity_threshold = 1.0e-3_DP
    slow_slip_velocity_threshold1 = 1.0e-6_DP
    slow_slip_velocity_threshold2 = 1.0e-5_DP
    
    output_flags = 0
    
  end subroutine initialize_default_parameters
  
  !===============================================================================
  ! SUBROUTINE: Validate parameters
  !===============================================================================
  subroutine validate_parameters()
    implicit none
    logical :: valid = .true.
    character(len=256) :: error_message
    
    ! Check grid parameters
    if (total_elements <= 0) then
      error_message = 'total_elements must be positive'
      valid = .false.
    end if
    
    if (local_elements <= 0) then
      error_message = 'local_elements must be positive'
      valid = .false.
    end if
    
    if (mod(total_elements, number_processes) /= 0) then
      error_message = 'total_elements must be divisible by number_processes'
      valid = .false.
    end if
    
    ! Check time parameters
    if (maximum_time <= 0.0_DP) then
      error_message = 'maximum_time must be positive'
      valid = .false.
    end if
    
    if (time_step_min <= 0.0_DP) then
      error_message = 'time_step_min must be positive'
      valid = .false.
    end if
    
    if (time_step_max <= time_step_min) then
      error_message = 'time_step_max must be greater than time_step_min'
      valid = .false.
    end if
    
    ! Check velocity parameters
    if (plate_velocity <= 0.0_DP) then
      error_message = 'plate_velocity must be positive'
      valid = .false.
    end if
    
    ! If validation fails, stop execution
    if (.not. valid) then
      write(*, *) 'Parameter validation failed: ', trim(error_message)
      stop
    end if
    
  end subroutine validate_parameters
  
  !===============================================================================
  ! SUBROUTINE: Print parameters
  !===============================================================================
  subroutine print_parameters()
    implicit none
    
    write(*, *) '=== Simulation Parameters ==='
    write(*, *) 'Grid parameters:'
    write(*, *) '  Total elements: ', total_elements
    write(*, *) '  Local elements: ', local_elements
    write(*, *) '  Depth layers: ', depth_layers
    write(*, *) '  Number of processes: ', number_processes
    write(*, *)
    
    write(*, *) 'Time parameters:'
    write(*, *) '  Maximum time: ', maximum_time
    write(*, *) '  Time step range: [', time_step_min, ', ', time_step_max, ']'
    write(*, *)
    
    write(*, *) 'Output parameters:'
    write(*, *) '  Output interval: ', output_interval
    write(*, *) '  Restart interval: ', restart_interval
    write(*, *) '  Snapshot interval: ', snapshot_interval
    write(*, *)
    
    write(*, *) 'Physical parameters:'
    write(*, *) '  Plate velocity: ', plate_velocity
    write(*, *) '  Buffer zone width: ', buffer_zone_width, ' km'
    write(*, *) '  Locked seismogenic zone: ', locked_seismogenic_zone
    write(*, *) '=============================='
    
  end subroutine print_parameters
  
end module simulation_parameters
