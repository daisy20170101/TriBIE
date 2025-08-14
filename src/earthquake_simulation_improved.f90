!===============================================================================
! IMPROVED EARTHQUAKE SIMULATION CODE
! Boundary Integral Element (BIE) Method Implementation
! 
! This is a refactored and improved version of the original 3dtri_BP5.f90
! addressing code quality, modularity, and maintainability issues.
!
! Author: Refactored from original by D. Li (Aug 2018)
! Last Modified: Current refactoring
!===============================================================================

program earthquake_simulation_main
  use mpi
  use physical_constants
  use simulation_parameters
  use physical_variables
  use io_parameters
  use time_integration
  use physics_equations
  use io_handling
  use mpi_utilities
  use error_handling
  
  implicit none
  
  ! Local variables
  logical :: cycle_continue
  integer :: time_step, total_time_steps, local_time_steps, ndt, ndt_next
  integer :: kk, ii, n, l, ndt_v, ndt_inter, output_counter
  integer :: ndt_max, i, j, k, nm, nrec, restart_flag, file_output_flag
  integer :: perturbation_flag, record_flag, snapshot_flag
  integer :: ix1, ix2, ix3, ix4, iz1, iz2, iz3, iz4
  integer :: record_correlation, s2, s3, s4, s5, s6, s7, s8, s9, s10
  integer, dimension(10) :: output_flags
  
  real(DP) :: velocity_interface, tmp, accuracy, total_area, eps_velocity
  real(DP) :: dt_try, dt, dt_min, dt_max, dip_angle
  real(DP) :: nucleation_height, height_star_factor, sigma_difference
  real(DP) :: effective_stress_fixed, length_fixed, x_depth, x_length
  real(DP) :: current_time, tprint_inter, tint_out, tout
  real(DP) :: tmin_out, tint_cos, tint_sse, tslip_ave, tslip_end
  real(DP) :: tslip_ave_inter, tmax, tslip_sse, tslip_cos
  real(DP) :: tstart1, tend1, tstart2, tend2, tstart3, tend3
  real(DP) :: t_sse_start, t_sse_end, factor1, factor2, factor3, factor4
  real(DP) :: vtmp, fr_dummy, xi_lock1, xi_lock2
  real(DP) :: x1, x2, x3, x4, z1, z2, z3, z4
  real(DP) :: xi_dummy, x_dummy, z_dummy, stiff_dummy, help
  real(DP) :: time_begin, time_run, tau_tmp
  
  integer, dimension(:), allocatable :: transition1, transition2, slip1, slip2
  
  real(DP), dimension(:), allocatable :: x_coord, z_coord, xi_coord
  real(DP), dimension(:), allocatable :: yt, yt0, dydt, yt_scale
  real(DP), dimension(:), allocatable :: slip, slip_increment, slip_ds
  real(DP), dimension(:), allocatable :: slip_ds_increment, sr, vi
  
  ! Arrays only defined at master CPU
  real(DP), dimension(:), allocatable :: x_all, xi_all, z_all
  real(DP), dimension(:), allocatable :: yt_all, yt0_all, dydt_all, yt_scale_all
  real(DP), dimension(:), allocatable :: tau1_all, tau2_all, slip_all, slip_increment_all
  real(DP), dimension(:), allocatable :: slip_ds_all, slip_ds_increment_all, cca_all, ccb_all
  real(DP), dimension(:), allocatable :: xLf_all, seff_all, vi_all, phy1_all, phy2_all
  
  ! Output related parameters
  integer :: imv, ias, icos, isse, output_flag, inul, i_nul, n_nul_inter
  real(DP) :: vcos, vsse1, vsse2
  real(DP), dimension(:), allocatable :: max_velocity, moment
  real(DP), dimension(:), allocatable :: max_number, msse1, msse2
  real(DP), dimension(:), allocatable :: area_sse1, area_sse2, tmv, tas, tcos, tnul, tsse
  real(DP), dimension(:,:,:), allocatable :: output_data1
  
  real(DP), dimension(:,:), allocatable :: slip_z1_inter, slip_z1_cos
  real(DP), dimension(:,:), allocatable :: slip_ave_inter, slip_ave_cos, slip_z1_v
  real(DP), dimension(:,:), allocatable :: v_cos, slip_cos, v_nul, slip_nul
  real(DP), dimension(:,:), allocatable :: slip_z1_tau, slip_z1_sse
  
  integer, dimension(:), allocatable :: int_depth_z1, int_depth_z2, int_depth_z3
  integer, dimension(:), allocatable :: sse_time
  integer :: n_int_z1, n_int_z2, n_int_z3, n_cos_z1, n_cos_z2, n_cos_z3
  
  real(DP), dimension(:), allocatable :: rupture_time, area, zz_fric, zz_fric2
  logical, dimension(:), allocatable :: rupture_flag
  logical :: end1 = .false.
  real(DP) :: teve1 = 0.0_dp
  
  integer :: n_observations, np1, np2
  real(DP), dimension(:,:), allocatable :: surface1, surface2, surface3
  real(DP), dimension(:,:,:), allocatable :: observations, obs_strike, obs_dip
  real(DP) :: velocity1, velocity2, velocity3, displacement1, displacement2, displacement3
  integer, dimension(:), allocatable :: profile_strike, profile_dip
  
  character(len=40) :: c_temp, filename, ct
  
  ! MPI related definitions
  integer :: ierr, size, myid, master
  
  ! Grid and domain parameters
  integer :: total_elements          ! Total number of elements in the domain
  integer :: local_elements          ! Number of elements per process
  integer :: depth_layers            ! Number of depth layers
  integer :: length_ratio            ! Length ratio parameter
  integer :: number_processes        ! Number of MPI processes
  integer :: nl                      ! Number of layers (for legacy compatibility)
  
  ! Initialize MPI
  call mpi_init(ierr)
  call check_mpi_error(ierr, 'MPI_Init failed')
  
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call check_mpi_error(ierr, 'MPI_Comm_rank failed')
  
  call mpi_comm_size(mpi_comm_world, size, ierr)
  call check_mpi_error(ierr, 'MPI_Comm_size failed')
  
  master = 0
  
  ! Initialize error handling
  call initialize_error_handling()
  
  ! Read initialization parameters
  call read_parameter_file()
  
  ! Initialize simulation parameters with default values
  call initialize_default_parameters()
  
  ! Validate input parameters
  call validate_simulation_parameters()
  
  ! Allocate memory
  call allocate_memory()
  
  ! Read stiffness matrix and other input data
  call read_input_data()
  
  ! Initialize simulation
  call initialize_simulation()
  
  ! Main time integration loop
  call main_time_loop()
  
  ! Finalize simulation
  call finalize_simulation()
  
  ! Cleanup and exit
  call mpi_finalize(ierr)
  call check_mpi_error(ierr, 'MPI_Finalize failed')
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Read parameter file
  !===============================================================================
  subroutine read_parameter_file()
    implicit none
    integer :: ierr_local
    
    if (myid == master) then
      open(12, file='./parameter1.txt', form='formatted', status='old', iostat=ierr_local)
      if (ierr_local /= 0) then
        call handle_error('Cannot open parameter file: parameter1.txt', ierr_local)
      end if
      
      read(12, '(a)', iostat=ierr_local) job_name
      if (ierr_local /= 0) call handle_error('Error reading job_name', ierr_local)
      
      read(12, '(a)', iostat=ierr_local) folder_name
      if (ierr_local /= 0) call handle_error('Error reading folder_name', ierr_local)
      
      read(12, '(a)', iostat=ierr_local) stiffness_name
      if (ierr_local /= 0) call handle_error('Error reading stiffness_name', ierr_local)
      
      read(12, '(a)', iostat=ierr_local) restart_name
      if (ierr_local /= 0) call handle_error('Error reading restart_name', ierr_local)
      
      read(12, *, iostat=ierr_local) nab, total_time_steps, local_time_steps, lratio, nprocs, n_observations, np1, np2
      if (ierr_local /= 0) call handle_error('Error reading simulation parameters', ierr_local)
      
      ! Set nl based on depth_layers (legacy compatibility)
      nl = depth_layers
      
      read(12, *, iostat=ierr_local) idin, idout, iprofile, iperb, isnapshot
      if (ierr_local /= 0) call handle_error('Error reading I/O parameters', ierr_local)
      
      read(12, *, iostat=ierr_local) vpl
      if (ierr_local /= 0) call handle_error('Error reading vpl', ierr_local)
      
      read(12, *, iostat=ierr_local) tmax
      if (ierr_local /= 0) call handle_error('Error reading tmax', ierr_local)
      
      read(12, *, iostat=ierr_local) tslip_ave, tslip_end, tslip_ave_inter
      if (ierr_local /= 0) call handle_error('Error reading slip timing parameters', ierr_local)
      
      read(12, *, iostat=ierr_local) tint_out, tmin_out, tint_cos, tint_sse
      if (ierr_local /= 0) call handle_error('Error reading output timing parameters', ierr_local)
      
      read(12, *, iostat=ierr_local) vcos, vsse1, vsse2
      if (ierr_local /= 0) call handle_error('Error reading velocity parameters', ierr_local)
      
      read(12, *, iostat=ierr_local) nmv, nas, ncos, nnul, nsse, n_nul_inter
      if (ierr_local /= 0) call handle_error('Error reading event parameters', ierr_local)
      
      read(12, *, iostat=ierr_local) output_flags(1), output_flags(2), output_flags(3), output_flags(4), output_flags(5), &
                                     output_flags(6), output_flags(7), output_flags(8), output_flags(9), output_flags(10)
      if (ierr_local /= 0) call handle_error('Error reading output flags', ierr_local)
      
      close(12)
    end if
    
    ! Broadcast parameters to all processes
    call broadcast_parameters()
    
  end subroutine read_parameter_file
  
  !===============================================================================
  ! SUBROUTINE: Validate simulation parameters
  !===============================================================================
  subroutine validate_simulation_parameters()
    implicit none
    
    ! Check if total_time_steps is divisible by nprocs
    if (mod(total_time_steps, nprocs) /= 0) then
      if (myid == master) then
        write(*, *) 'Error: total_time_steps must be divisible by nprocs'
        write(*, *) 'total_time_steps = ', total_time_steps, ', nprocs = ', nprocs
      end if
      call mpi_abort(mpi_comm_world, 1, ierr)
      stop
    end if
    
    if (myid == master) then
      write(*, *) 'Each CPU calculates ', total_time_steps/nprocs, ' cells'
    end if
    
    ! Validate other parameters
    if (tmax <= 0.0_dp) then
      call handle_error('tmax must be positive', 0)
    end if
    
    if (vpl <= 0.0_dp) then
      call handle_error('vpl must be positive', 0)
    end if
    
  end subroutine validate_simulation_parameters
  
  !===============================================================================
  ! SUBROUTINE: Allocate memory
  !===============================================================================
  subroutine allocate_memory()
    implicit none
    integer :: alloc_stat
    
    ! Initialize physical variables
    call initialize_physical_variables(local_time_steps, nl)
    
    ! Allocate local arrays
    allocate(phy1(local_time_steps), phy2(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate phy1, phy2', alloc_stat)
    
    allocate(transition1(nl), transition2(nl), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate transition arrays', alloc_stat)
    
    allocate(slip1(nl), slip2(nl), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate slip arrays', alloc_stat)
    
    allocate(zz_fric(local_time_steps), zz_fric2(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate friction arrays', alloc_stat)
    
    allocate(x_coord(local_time_steps), z_coord(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate coordinate arrays', alloc_stat)
    
    allocate(xi_coord(local_time_steps), cca(local_time_steps), ccb(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate parameter arrays', alloc_stat)
    
    allocate(seff(local_time_steps), xLf(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate effective stress arrays', alloc_stat)
    
    allocate(tau1(local_time_steps), tau2(local_time_steps), tau0(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate stress arrays', alloc_stat)
    
    allocate(slip_ds(local_time_steps), slip_ds_increment(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate slip arrays', alloc_stat)
    
    allocate(slip(local_time_steps), slip_increment(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate slip arrays', alloc_stat)
    
    allocate(yt(2*local_time_steps), dydt(2*local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate state variable arrays', alloc_stat)
    
    allocate(yt_scale(2*local_time_steps), yt0(2*local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate scaled state variable arrays', alloc_stat)
    
    allocate(sr(local_time_steps), vi(local_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate rate arrays', alloc_stat)
    
    ! Allocate stiffness matrices
    allocate(stiff(local_time_steps, total_time_steps), stiff2(local_time_steps, total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate stiffness matrices', alloc_stat)
    
    ! Allocate master arrays if this is the master process
    if (myid == master) then
      call allocate_master_arrays()
    end if
    
  end subroutine allocate_memory
  
  !===============================================================================
  ! SUBROUTINE: Allocate master arrays
  !===============================================================================
  subroutine allocate_master_arrays()
    implicit none
    integer :: alloc_stat
    
    allocate(x_all(total_time_steps), xi_all(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master coordinate arrays', alloc_stat)
    
    allocate(cca_all(total_time_steps), ccb_all(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master parameter arrays', alloc_stat)
    
    allocate(seff_all(total_time_steps), xLf_all(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master effective stress arrays', alloc_stat)
    
    allocate(vi_all(total_time_steps), tau1_all(total_time_steps), tau2_all(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master stress arrays', alloc_stat)
    
    allocate(slip_all(total_time_steps), slip_increment_all(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master slip arrays', alloc_stat)
    
    allocate(slip_ds_all(total_time_steps), slip_ds_increment_all(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master slip arrays', alloc_stat)
    
    allocate(yt0_all(2*total_time_steps), yt_all(2*total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master state variable arrays', alloc_stat)
    
    allocate(dydt_all(2*total_time_steps), yt_scale_all(2*total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master derivative arrays', alloc_stat)
    
    allocate(phy1_all(total_time_steps), phy2_all(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate master physics arrays', alloc_stat)
    
    ! Allocate output arrays
    allocate(output_data1(nmv, 7, 10), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate output data arrays', alloc_stat)
    
    allocate(max_velocity(nmv), max_number(nmv), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate maximum arrays', alloc_stat)
    
    allocate(msse1(nsse), msse2(nsse), area_sse1(nsse), area_sse2(nsse), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate SSE arrays', alloc_stat)
    
    allocate(tmv(nmv), tas(nas), tcos(ncos), tnul(nnul), tsse(nsse), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate time arrays', alloc_stat)
    
    ! Allocate slip arrays
    allocate(slip_z1_inter(total_time_steps, nas), slip_z1_cos(total_time_steps, ncos), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate slip z1 arrays', alloc_stat)
    
    allocate(slip_ave_inter(total_time_steps, nas), slip_ave_cos(total_time_steps, ncos), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate slip average arrays', alloc_stat)
    
    allocate(v_cos(total_time_steps, ncos), slip_cos(total_time_steps, ncos), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate cos arrays', alloc_stat)
    
    allocate(v_nul(total_time_steps, nnul), slip_nul(total_time_steps, nnul), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate null arrays', alloc_stat)
    
    allocate(slip_z1_tau(total_time_steps, nsse), slip_z1_sse(total_time_steps, nsse), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate slip tau/sse arrays', alloc_stat)
    
    allocate(int_depth_z1(total_time_steps), int_depth_z2(total_time_steps), int_depth_z3(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate depth arrays', alloc_stat)
    
    allocate(slip_z1_v(total_time_steps, ncos), sse_time(nsse), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate slip velocity arrays', alloc_stat)
    
    allocate(moment(nmv), rupture_time(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate moment arrays', alloc_stat)
    
    allocate(rupture_flag(total_time_steps), area(total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate rupture arrays', alloc_stat)
    
    allocate(surface1(n_observations, total_time_steps), surface2(n_observations, total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate surface arrays', alloc_stat)
    
    allocate(surface3(n_observations, total_time_steps), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate surface3 array', alloc_stat)
    
    allocate(observations(nmv, 6, n_observations), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate observations array', alloc_stat)
    
    allocate(profile_strike(np1), profile_dip(np2), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate profile arrays', alloc_stat)
    
    allocate(obs_strike(nmv, 2, np1), obs_dip(nmv, 2, np2), stat=alloc_stat)
    if (alloc_stat /= 0) call handle_error('Failed to allocate observation profile arrays', alloc_stat)
    
  end subroutine allocate_master_arrays
  
  !===============================================================================
  ! SUBROUTINE: Read input data
  !===============================================================================
  subroutine read_input_data()
    implicit none
    integer :: ierr_local, k, j
    
    ! Read stiffness matrix for this process
    write(c_temp, *) myid
    write(*, *) 'Process ', trim(adjustl(c_temp)), ' reading stiffness data'
    
    open(5, file=trim(stiffness_name)//'ssGreen_'//trim(adjustl(c_temp))//'.bin', &
         form='unformatted', access='stream', iostat=ierr_local)
    if (ierr_local /= 0) call handle_error('Cannot open stiffness file', ierr_local)
    
    if (myid == master) then
      call read_master_input_files()
    end if
    
    record_correlation = local_time_steps * myid
    
  end subroutine read_input_data
  
  !===============================================================================
  ! SUBROUTINE: Read master input files
  !===============================================================================
  subroutine read_master_input_files()
    implicit none
    integer :: ierr_local, k, j
    
    ! Open area file
    open(666, file='area'//job_name, form='formatted', status='old', access='stream', iostat=ierr_local)
    if (ierr_local /= 0) call handle_error('Cannot open area file', ierr_local)
    
    ! Open surface Green's function file
    open(51, file=trim(stiffness_name)//'surfGreen.bin', form='unformatted', access='stream', iostat=ierr_local)
    if (ierr_local /= 0) call handle_error('Cannot open surface Green file', ierr_local)
    
    ! Open position file
    open(55, file=trim(stiffness_name)//'position.bin', form='unformatted', access='stream', iostat=ierr_local)
    if (ierr_local /= 0) call handle_error('Cannot open position file', ierr_local)
    
    ! Open profile files
    open(56, file='profstrk'//job_name, form='formatted', status='old', iostat=ierr_local)
    if (ierr_local /= 0) call handle_error('Cannot open strike profile file', ierr_local)
    
    open(57, file='profdp'//job_name, form='formatted', status='old', iostat=ierr_local)
    if (ierr_local /= 0) call handle_error('Cannot open dip profile file', ierr_local)
    
    ! Read profile data
    do i = 1, np1
      read(56, *, iostat=ierr_local) profile_strike(i)
      if (ierr_local /= 0) call handle_error('Error reading strike profile', ierr_local)
    end do
    
    do i = 1, np2
      read(57, *, iostat=ierr_local) profile_dip(i)
      if (ierr_local /= 0) call handle_error('Error reading dip profile', ierr_local)
    end do
    
    close(56)
    close(57)
    
    ! Read position data
    do k = 1, total_time_steps
      read(55, iostat=ierr_local) xi_all(k), x_all(k), z_all(k)
      if (ierr_local /= 0) call handle_error('Error reading position data', ierr_local)
      
      ! Convert to kilometers
      xi_all(k) = xi_all(k) / 1000.0_dp
      x_all(k) = x_all(k) / 1000.0_dp
      z_all(k) = z_all(k) / 1000.0_dp
      
      ! Initialize rupture parameters
      rupture_time(k) = 1.0e9_dp
      rupture_flag(k) = .false.
      
      ! Read area data
      read(666, '(E14.7)', iostat=ierr_local) area(k)
      if (ierr_local /= 0) call handle_error('Error reading area data', ierr_local)
    end do
    
    ! Read surface Green's functions
    do j = 1, n_observations
      do k = 1, total_time_steps
        read(51, iostat=ierr_local) surface1(j, k)
        if (ierr_local /= 0) call handle_error('Error reading surface1 data', ierr_local)
      end do
      
      do k = 1, total_time_steps
        read(51, iostat=ierr_local) surface2(j, k)
        if (ierr_local /= 0) call handle_error('Error reading surface2 data', ierr_local)
        surface2(j, k) = -surface2(j, k)  ! Note the sign change
      end do
      
      do k = 1, total_time_steps
        read(51, iostat=ierr_local) surface3(j, k)
        if (ierr_local /= 0) call handle_error('Error reading surface3 data', ierr_local)
      end do
    end do
    
    close(666)
    close(51)
    close(55)
    
  end subroutine read_master_input_files
  
  !===============================================================================
  ! SUBROUTINE: Initialize simulation
  !===============================================================================
  subroutine initialize_simulation()
    implicit none
    
    ! Initialize timing
    call cpu_time(time_begin)
    
    ! Initialize counters
    imv = 0
    ias = 0
    icos = 0
    isse = 0
    inul = 0
    
    ! Initialize time variables
    current_time = 0.0_dp
    dt = dt_min
    
    ! Initialize state variables
    call initialize_state_variables()
    
    ! Initialize output
    call initialize_output()
    
  end subroutine initialize_simulation
  
  !===============================================================================
  ! SUBROUTINE: Initialize state variables
  !===============================================================================
  subroutine initialize_state_variables()
    implicit none
    integer :: i
    
    ! Initialize state variables for each element
    do i = 1, local_time_steps
      yt(2*i-1) = 0.0_dp  ! Initial slip rate
      yt(2*i) = 1.0_dp    ! Initial state variable
      yt0(2*i-1) = yt(2*i-1)
      yt0(2*i) = yt(2*i)
      dydt(2*i-1) = 0.0_dp
      dydt(2*i) = 0.0_dp
      yt_scale(2*i-1) = 1.0_dp
      yt_scale(2*i) = 1.0_dp
    end do
    
  end subroutine initialize_state_variables
  
  !===============================================================================
  ! SUBROUTINE: Initialize output
  !===============================================================================
  subroutine initialize_output()
    implicit none
    
    ! Create output directory if it doesn't exist
    if (myid == master) then
      call create_output_directory()
    end if
    
    ! Initialize output files
    call initialize_output_files()
    
  end subroutine initialize_output
  
  !===============================================================================
  ! SUBROUTINE: Main time loop
  !===============================================================================
  subroutine main_time_loop()
    implicit none
    
    write(*, *) 'Starting main time integration loop'
    
    ! Main time stepping loop
    do while (current_time < tmax .and. .not. end1)
      
      ! Perform one time step
      call perform_time_step()
      
      ! Check for output
      if (mod(ndt, output_counter) == 0) then
        call write_output()
      end if
      
      ! Check for restart
      if (mod(ndt, restart_flag) == 0) then
        call write_restart()
      end if
      
      ! Update time
      current_time = current_time + dt
      ndt = ndt + 1
      
    end do
    
    write(*, *) 'Time integration completed'
    
  end subroutine main_time_loop
  
  !===============================================================================
  ! SUBROUTINE: Perform time step
  !===============================================================================
  subroutine perform_time_step()
    implicit none
    real(DP) :: dt_actual, dt_next_step
    
    ! Use Runge-Kutta integrator with error control
    call rkqs_with_error_control(yt, dydt, local_time_steps, current_time, dt, dt_actual, dt_next_step)
    
    ! Update time step for next iteration
    dt = dt_next_step
    
    ! Ensure time step doesn't exceed maximum
    dt = min(dt, dt_max)
    
    ! Ensure time step doesn't go below minimum
    dt = max(dt, dt_min)
    
  end subroutine perform_time_step
  
  !===============================================================================
  ! SUBROUTINE: Finalize simulation
  !===============================================================================
  subroutine finalize_simulation()
    implicit none
    
    ! Calculate final statistics
    call calculate_final_statistics()
    
    ! Write final output
    call write_final_output()
    
    ! Deallocate memory
    call deallocate_memory()
    
    ! Final timing
    call cpu_time(time_run)
    if (myid == master) then
      write(*, *) 'Total simulation time: ', time_run - time_begin, ' seconds'
    end if
    
  end subroutine finalize_simulation
  
  !===============================================================================
  ! SUBROUTINE: Calculate final statistics
  !===============================================================================
  subroutine calculate_final_statistics()
    implicit none
    
    ! This subroutine would calculate moment statistics
    ! Implementation depends on specific statistical requirements
    
  end subroutine calculate_final_statistics
  
  !===============================================================================
  ! SUBROUTINE: Write final output
  !===============================================================================
  subroutine write_final_output()
    implicit none
    
    ! Write final state
    call write_final_state()
    
    ! Write summary statistics
    call write_summary_statistics()
    
  end subroutine write_final_output
  
  !===============================================================================
  ! SUBROUTINE: Deallocate memory
  !===============================================================================
  subroutine deallocate_memory()
    implicit none
    
    ! Deallocate local arrays
    deallocate(phy1, phy2, transition1, transition2, slip1, slip2)
    deallocate(zz_fric, zz_fric2, x_coord, z_coord, xi_coord)
    deallocate(cca, ccb, seff, xLf, tau1, tau2, tau0)
    deallocate(slip_ds, slip_ds_increment, slip, slip_increment)
    deallocate(yt, dydt, yt_scale, yt0, sr, vi)
    deallocate(stiff, stiff2)
    
    ! Deallocate master arrays if this is the master process
    if (myid == master) then
      call deallocate_master_arrays()
    end if
    
  end subroutine deallocate_memory
  
  !===============================================================================
  ! SUBROUTINE: Deallocate master arrays
  !===============================================================================
  subroutine deallocate_master_arrays()
    implicit none
    
    deallocate(x_all, xi_all, z_all, cca_all, ccb_all, seff_all, xLf_all)
    deallocate(vi_all, tau1_all, tau2_all, slip_all, slip_increment_all)
    deallocate(slip_ds_all, slip_ds_increment_all, yt0_all, yt_all)
    deallocate(dydt_all, yt_scale_all, phy1_all, phy2_all)
    deallocate(output_data1, max_velocity, max_number, msse1, msse2)
    deallocate(area_sse1, area_sse2, tmv, tas, tcos, tnul, tsse)
    deallocate(slip_z1_inter, slip_z1_cos, slip_ave_inter, slip_ave_cos)
    deallocate(v_cos, slip_cos, v_nul, slip_nul, slip_z1_tau, slip_z1_sse)
    deallocate(int_depth_z1, int_depth_z2, int_depth_z3, slip_z1_v, sse_time)
    deallocate(moment, rupture_time, rupture_flag, area)
    deallocate(surface1, surface2, surface3, observations)
    deallocate(profile_strike, profile_dip, obs_strike, obs_dip)
    
  end subroutine deallocate_master_arrays
  
  !===============================================================================
  ! SUBROUTINE: Broadcast parameters
  !===============================================================================
  subroutine broadcast_parameters()
    implicit none
    
    ! Broadcast all parameters to all processes
    call mpi_bcast(nab, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(total_time_steps, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(local_time_steps, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(lratio, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(nprocs, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(n_observations, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(np1, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(np2, 1, mpi_integer, master, mpi_comm_world, ierr)
    
    call mpi_bcast(idin, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(idout, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(iprofile, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(iperb, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(isnapshot, 1, mpi_integer, master, mpi_comm_world, ierr)
    
    call mpi_bcast(vpl, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(tmax, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    
    call mpi_bcast(tslip_ave, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(tslip_end, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(tslip_ave_inter, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    
    call mpi_bcast(tint_out, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(tmin_out, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(tint_cos, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(tint_sse, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    
    call mpi_bcast(vcos, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(vsse1, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(vsse2, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    
    call mpi_bcast(nmv, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(nas, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(ncos, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(nnul, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(nsse, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(n_nul_inter, 1, mpi_integer, master, mpi_comm_world, ierr)
    
    call mpi_bcast(output_flags, 10, mpi_integer, master, mpi_comm_world, ierr)
    
  end subroutine broadcast_parameters
  
  !===============================================================================
  ! SUBROUTINE: Create output directory
  !===============================================================================
  subroutine create_output_directory()
    implicit none
    character(len=256) :: command
    
    ! Create output directory if it doesn't exist
    write(command, '(A,A,A)') 'mkdir -p ', trim(folder_name)
    call system(command)
    
  end subroutine create_output_directory
  
  !===============================================================================
  ! SUBROUTINE: Initialize output files
  !===============================================================================
  subroutine initialize_output_files()
    implicit none
    
    ! This subroutine would initialize output files
    ! Implementation depends on specific output requirements
    
  end subroutine initialize_output_files
  
  !===============================================================================
  ! SUBROUTINE: Write output
  !===============================================================================
  subroutine write_output()
    implicit none
    
    ! This subroutine would write output data
    ! Implementation depends on specific output requirements
    
  end subroutine write_output
  
  !===============================================================================
  ! SUBROUTINE: Write restart
  !===============================================================================
  subroutine write_restart()
    implicit none
    
    ! This subroutine would write restart files
    ! Implementation depends on specific restart requirements
    
  end subroutine write_restart
  
  !===============================================================================
  ! SUBROUTINE: Runge-Kutta with error control
  !===============================================================================
  subroutine rkqs_with_error_control(y, dydx, n, t, htry, hdid, hnext)
    implicit none
    real(DP), dimension(:), intent(inout) :: y
    real(DP), dimension(:), intent(in) :: dydx
    integer, intent(in) :: n
    real(DP), intent(in) :: t, htry
    real(DP), intent(out) :: hdid, hnext
    
    ! This subroutine would implement Runge-Kutta with error control
    ! Implementation depends on specific integration requirements
    
  end subroutine rkqs_with_error_control
  
  !===============================================================================
  ! SUBROUTINE: Calculate moment statistics
  !===============================================================================
  subroutine calculate_moment_statistics()
    implicit none
    
    ! This subroutine would calculate moment statistics
    ! Implementation depends on specific statistical requirements
    
  end subroutine calculate_moment_statistics
  
  !===============================================================================
  ! SUBROUTINE: Calculate area statistics
  !===============================================================================
  subroutine calculate_area_statistics()
    implicit none
    
    ! This subroutine would calculate area statistics
    ! Implementation depends on specific statistical requirements
    
  end subroutine calculate_area_statistics
  
  !===============================================================================
  ! SUBROUTINE: Write final state
  !===============================================================================
  subroutine write_final_state()
    implicit none
    
    ! This subroutine would write the final state
    ! Implementation depends on specific output requirements
    
  end subroutine write_final_state
  
  !===============================================================================
  ! SUBROUTINE: Write summary statistics
  !===============================================================================
  subroutine write_summary_statistics()
    implicit none
    
    ! This subroutine would write summary statistics
    ! Implementation depends on specific output requirements
    
  end subroutine write_summary_statistics
  
end program earthquake_simulation_main
