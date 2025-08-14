!===============================================================================
! MODULE: I/O Parameters
! 
! This module defines all I/O parameters including file names, paths,
! and output settings for the earthquake simulation.
!
! Author: Refactored from original phy3d_module_non.f90
! Last Modified: Current refactoring
!===============================================================================

module io_parameters
  implicit none
  
  ! Make all parameters public
  public
  
  ! File names and paths
  character(len=80) :: job_name           ! Job name identifier
  character(len=80) :: folder_name        ! Output folder name
  character(len=80) :: restart_name       ! Restart file name
  character(len=80) :: stiffness_name     ! Stiffness matrix file name
  character(len=80) :: profile_name       ! Profile output name
  
  ! Input/output unit numbers
  integer, parameter :: PARAMETER_FILE_UNIT = 12
  integer, parameter :: STIFFNESS_FILE_UNIT = 5
  integer, parameter :: AREA_FILE_UNIT = 666
  integer, parameter :: SURFACE_GREEN_FILE_UNIT = 51
  integer, parameter :: POSITION_FILE_UNIT = 55
  integer, parameter :: STRIKE_PROFILE_FILE_UNIT = 56
  integer, parameter :: DIP_PROFILE_FILE_UNIT = 57
  
  ! Output file units
  integer, parameter :: SLIP_INTER_FILE_UNIT = 31
  integer, parameter :: TIME_INTER_FILE_UNIT = 34
  integer, parameter :: SLIP_SSE_FILE_UNIT = 25
  integer, parameter :: SLIP_TAU_FILE_UNIT = 26
  integer, parameter :: TIME_SSE_FILE_UNIT = 28
  integer, parameter :: SLIP_COS_FILE_UNIT = 45
  integer, parameter :: SLIP_VELOCITY_FILE_UNIT = 42
  integer, parameter :: TIME_COS_FILE_UNIT = 48
  integer, parameter :: VELOCITY_NULL_FILE_UNIT = 52
  integer, parameter :: NULL_TIME_FILE_UNIT = 53
  
  ! File status options
  character(len=7), parameter :: STATUS_OLD = 'old'
  character(len=7), parameter :: STATUS_NEW = 'new'
  character(len=7), parameter :: STATUS_UNKNOWN = 'unknown'
  character(len=7), parameter :: STATUS_REPLACE = 'replace'
  
  ! File access options
  character(len=10), parameter :: ACCESS_SEQUENTIAL = 'sequential'
  character(len=10), parameter :: ACCESS_DIRECT = 'direct'
  character(len=10), parameter :: ACCESS_STREAM = 'stream'
  
  ! File form options
  character(len=11), parameter :: FORM_FORMATTED = 'formatted'
  character(len=11), parameter :: FORM_UNFORMATTED = 'unformatted'
  
  ! Default file extensions
  character(len=4), parameter :: DEFAULT_EXTENSION = '.dat'
  character(len=4), parameter :: BINARY_EXTENSION = '.bin'
  character(len=4), parameter :: TEXT_EXTENSION = '.txt'
  
  ! Default directory structure
  character(len=256), parameter :: DEFAULT_OUTPUT_DIR = './output/'
  character(len=256), parameter :: DEFAULT_RESTART_DIR = './restart/'
  character(len=256), parameter :: DEFAULT_INPUT_DIR = './input/'
  
  ! File naming conventions
  character(len=20), parameter :: STIFFNESS_PREFIX = 'ssGreen_'
  character(len=20), parameter :: SURFACE_GREEN_PREFIX = 'surfGreen'
  character(len=20), parameter :: POSITION_PREFIX = 'position'
  character(len=20), parameter :: AREA_PREFIX = 'area'
  character(len=20), parameter :: STRIKE_PROFILE_PREFIX = 'profstrk'
  character(len=20), parameter :: DIP_PROFILE_PREFIX = 'profdp'
  
  ! Output file naming conventions
  character(len=20), parameter :: SLIP_INTER_PREFIX = 'slipz1-inter'
  character(len=20), parameter :: TIME_INTER_PREFIX = 't-inter'
  character(len=20), parameter :: SLIP_SSE_PREFIX = 'slipz1_sse'
  character(len=20), parameter :: SLIP_TAU_PREFIX = 'slipz1_tau'
  character(len=20), parameter :: TIME_SSE_PREFIX = 't_sse'
  character(len=20), parameter :: SLIP_COS_PREFIX = 'slipz1-cos'
  character(len=20), parameter :: SLIP_VELOCITY_PREFIX = 'slipz1-v'
  character(len=20), parameter :: TIME_COS_PREFIX = 't-cos'
  character(len=20), parameter :: VELOCITY_NULL_PREFIX = 'vs-nul'
  character(len=20), parameter :: NULL_TIME_PREFIX = 'nul-time'
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Initialize I/O parameters
  !===============================================================================
  subroutine initialize_io_parameters()
    implicit none
    
    ! Set default values
    job_name = 'default_job'
    folder_name = './output/'
    restart_name = 'restart.dat'
    stiffness_name = './stiffness/'
    profile_name = 'profile.dat'
    
  end subroutine initialize_io_parameters
  
  !===============================================================================
  ! SUBROUTINE: Build full file path
  !===============================================================================
  function build_file_path(directory, prefix, suffix, extension) result(full_path)
    implicit none
    character(len=*), intent(in) :: directory, prefix, suffix, extension
    character(len=256) :: full_path
    
    if (len_trim(suffix) > 0) then
      full_path = trim(directory) // trim(prefix) // trim(suffix) // trim(extension)
    else
      full_path = trim(directory) // trim(prefix) // trim(extension)
    end if
    
  end function build_file_path
  
  !===============================================================================
  ! SUBROUTINE: Build stiffness file name
  !===============================================================================
  function build_stiffness_filename(process_id) result(filename)
    implicit none
    integer, intent(in) :: process_id
    character(len=256) :: filename
    character(len=20) :: process_str
    
    write(process_str, '(I0)') process_id
    filename = trim(stiffness_name) // STIFFNESS_PREFIX // trim(adjustl(process_str)) // BINARY_EXTENSION
    
  end function build_stiffness_filename
  
  !===============================================================================
  ! SUBROUTINE: Build surface Green's function filename
  !===============================================================================
  function build_surface_green_filename() result(filename)
    implicit none
    character(len=256) :: filename
    
    filename = trim(stiffness_name) // SURFACE_GREEN_PREFIX // BINARY_EXTENSION
    
  end function build_surface_green_filename
  
  !===============================================================================
  ! SUBROUTINE: Build position filename
  !===============================================================================
  function build_position_filename() result(filename)
    implicit none
    character(len=256) :: filename
    
    filename = trim(stiffness_name) // POSITION_PREFIX // BINARY_EXTENSION
    
  end function build_position_filename
  
  !===============================================================================
  ! SUBROUTINE: Build area filename
  !===============================================================================
  function build_area_filename() result(filename)
    implicit none
    character(len=256) :: filename
    
    filename = AREA_PREFIX // job_name // DEFAULT_EXTENSION
    
  end function build_area_filename
  
  !===============================================================================
  ! SUBROUTINE: Build strike profile filename
  !===============================================================================
  function build_strike_profile_filename() result(filename)
    implicit none
    character(len=256) :: filename
    
    filename = STRIKE_PROFILE_PREFIX // job_name // DEFAULT_EXTENSION
    
  end function build_strike_profile_filename
  
  !===============================================================================
  ! SUBROUTINE: Build dip profile filename
  !===============================================================================
  function build_dip_profile_filename() result(filename)
    implicit none
    character(len=256) :: filename
    
    filename = DIP_PROFILE_PREFIX // job_name // DEFAULT_EXTENSION
    
  end function build_dip_profile_filename
  
  !===============================================================================
  ! SUBROUTINE: Build output filename
  !===============================================================================
  function build_output_filename(prefix) result(filename)
    implicit none
    character(len=*), intent(in) :: prefix
    character(len=256) :: filename
    
    filename = trim(folder_name) // trim(prefix) // job_name
    
  end function build_output_filename
  
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
  ! SUBROUTINE: Validate file existence
  !===============================================================================
  function file_exists(filename) result(exists)
    implicit none
    character(len=*), intent(in) :: filename
    logical :: exists
    integer :: unit_number, status
    
    exists = .false.
    
    ! Try to open file
    open(newunit=unit_number, file=filename, status='old', iostat=status)
    if (status == 0) then
      exists = .true.
      close(unit_number)
    end if
    
  end function file_exists
  
  !===============================================================================
  ! SUBROUTINE: Check and create directories
  !===============================================================================
  subroutine ensure_directories_exist()
    implicit none
    character(len=256) :: command
    
    ! Create output directory
    if (.not. file_exists(folder_name)) then
      write(command, '(A,A,A)') 'mkdir -p ', trim(folder_name)
      call system(command)
    end if
    
    ! Create restart directory
    if (.not. file_exists(restart_name)) then
      write(command, '(A,A,A)') 'mkdir -p ', trim(restart_name)
      call system(command)
    end if
    
  end subroutine ensure_directories_exist
  
  !===============================================================================
  ! SUBROUTINE: Print I/O parameters
  !===============================================================================
  subroutine print_io_parameters()
    implicit none
    
    write(*, *) '=== I/O Parameters ==='
    write(*, *) 'Job name: ', trim(job_name)
    write(*, *) 'Output folder: ', trim(folder_name)
    write(*, *) 'Restart file: ', trim(restart_name)
    write(*, *) 'Stiffness directory: ', trim(stiffness_name)
    write(*, *) 'Profile file: ', trim(profile_name)
    write(*, *) '====================='
    
  end subroutine print_io_parameters
  
end module io_parameters
