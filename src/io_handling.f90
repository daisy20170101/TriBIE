!===============================================================================
! MODULE: I/O Handling
! 
! This module provides I/O handling functions for the earthquake simulation,
! including file operations, output formatting, and restart functionality.
!
! Author: New module for improved I/O handling
! Last Modified: Current refactoring
!===============================================================================

module io_handling
  use physical_constants, only: DP
  use io_parameters, only: build_output_filename, build_stiffness_filename
  use io_parameters, only: build_surface_green_filename, build_position_filename
  use io_parameters, only: build_area_filename, build_strike_profile_filename, build_dip_profile_filename
  use io_parameters, only: SLIP_INTER_FILE_UNIT, TIME_INTER_FILE_UNIT, SLIP_SSE_FILE_UNIT
  use io_parameters, only: SLIP_TAU_FILE_UNIT, TIME_SSE_FILE_UNIT, SLIP_COS_FILE_UNIT
  use io_parameters, only: SLIP_VELOCITY_FILE_UNIT, TIME_COS_FILE_UNIT, VELOCITY_NULL_FILE_UNIT
  use io_parameters, only: NULL_TIME_FILE_UNIT, STATUS_UNKNOWN, ACCESS_APPEND
  use io_parameters, only: FORM_UNFORMATTED
  implicit none
  
  ! Make all functions public
  public
  
  ! Output format parameters
  character(len=20), parameter :: OUTPUT_FORMAT_1 = '(E22.14,7(1X,E15.7))'
  character(len=20), parameter :: OUTPUT_FORMAT_2 = '(E20.13,4X,E20.13,4X,I6)'
  character(len=20), parameter :: OUTPUT_FORMAT_3 = '(E22.14,2(1X,E15.7))'
  character(len=20), parameter :: OUTPUT_FORMAT_4 = '(E20.13)'
  character(len=20), parameter :: OUTPUT_FORMAT_5 = '(E22.14,3(1X,E15.7))'
  character(len=20), parameter :: OUTPUT_FORMAT_6 = '(E20.13,1x,E20.13)'
  character(len=20), parameter :: OUTPUT_FORMAT_7 = '(E15.8,1X,E20.13)'
  character(len=20), parameter :: OUTPUT_FORMAT_8 = '(E15.4,1X,E13.6,1X,E15.8)'
  character(len=20), parameter :: OUTPUT_FORMAT_9 = '(E13.6)'
  character(len=20), parameter :: OUTPUT_FORMAT_10 = '(E15.8)'
  
  ! File unit management
  integer, parameter :: MAX_OPEN_FILES = 100
  logical, dimension(MAX_OPEN_FILES) :: file_units_used = .false.
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Initialize I/O handling
  !===============================================================================
  subroutine initialize_io_handling()
    implicit none
    
    ! Reset file unit usage
    file_units_used = .false.
    
  end subroutine initialize_io_handling
  
  !===============================================================================
  ! SUBROUTINE: Get available file unit
  !===============================================================================
  function get_available_file_unit() result(unit)
    implicit none
    integer :: unit
    integer :: i
    
    unit = -1
    
    ! Find first available unit
    do i = 10, MAX_OPEN_FILES
      if (.not. file_units_used(i)) then
        unit = i
        file_units_used(i) = .true.
        exit
      end if
    end do
    
    if (unit == -1) then
      write(*, *) 'Error: No available file units'
      stop
    end if
    
  end function get_available_file_unit
  
  !===============================================================================
  ! SUBROUTINE: Release file unit
  !===============================================================================
  subroutine release_file_unit(unit)
    implicit none
    integer, intent(in) :: unit
    
    if (unit >= 1 .and. unit <= MAX_OPEN_FILES) then
      file_units_used(unit) = .false.
    end if
    
  end subroutine release_file_unit
  
  !===============================================================================
  ! SUBROUTINE: Write slip data to file
  !===============================================================================
  subroutine write_slip_data(filename, slip_data, num_elements, num_events, append)
    implicit none
    character(len=*), intent(in) :: filename
    real(DP), dimension(:,:), intent(in) :: slip_data
    integer, intent(in) :: num_elements, num_events
    logical, intent(in), optional :: append
    
    integer :: unit_number, i, j, ierr
    logical :: append_file
    
    append_file = .false.
    if (present(append)) append_file = append
    
    unit_number = get_available_file_unit()
    
    if (append_file) then
      open(unit_number, file=filename, form=FORM_UNFORMATTED, access=ACCESS_APPEND, &
           status=STATUS_UNKNOWN, iostat=ierr)
    else
      open(unit_number, file=filename, form=FORM_UNFORMATTED, status='replace', iostat=ierr)
    end if
    
    if (ierr /= 0) then
      write(*, *) 'Error opening file: ', trim(filename)
      call release_file_unit(unit_number)
      return
    end if
    
    ! Write slip data
    do j = 1, num_events
      do i = 1, num_elements
        write(unit_number, iostat=ierr) slip_data(i, j)
        if (ierr /= 0) then
          write(*, *) 'Error writing slip data'
          exit
        end if
      end do
    end do
    
    close(unit_number)
    call release_file_unit(unit_number)
    
  end subroutine write_slip_data
  
  !===============================================================================
  ! SUBROUTINE: Write time data to file
  !===============================================================================
  subroutine write_time_data(filename, time_data, num_events, append)
    implicit none
    character(len=*), intent(in) :: filename
    real(DP), dimension(:), intent(in) :: time_data
    integer, intent(in) :: num_events
    logical, intent(in), optional :: append
    
    integer :: unit_number, i, ierr
    logical :: append_file
    
    append_file = .false.
    if (present(append)) append_file = append
    
    unit_number = get_available_file_unit()
    
    if (append_file) then
      open(unit_number, file=filename, access=ACCESS_APPEND, status=STATUS_UNKNOWN, iostat=ierr)
    else
      open(unit_number, file=filename, status='replace', iostat=ierr)
    end if
    
    if (ierr /= 0) then
      write(*, *) 'Error opening file: ', trim(filename)
      call release_file_unit(unit_number)
      return
    end if
    
    ! Write time data
    do i = 1, num_events
      write(unit_number, *, iostat=ierr) time_data(i)
      if (ierr /= 0) then
        write(*, *) 'Error writing time data'
        exit
      end if
    end do
    
    close(unit_number)
    call release_file_unit(unit_number)
    
  end subroutine write_time_data
  
  !===============================================================================
  ! SUBROUTINE: Write restart file
  !===============================================================================
  subroutine write_restart_file(filename, current_time, time_step, yt, slip, num_elements)
    implicit none
    character(len=*), intent(in) :: filename
    real(DP), intent(in) :: current_time
    integer, intent(in) :: time_step
    real(DP), dimension(:), intent(in) :: yt, slip
    integer, intent(in) :: num_elements
    
    integer :: unit_number, i, ierr
    
    unit_number = get_available_file_unit()
    
    open(unit_number, file=filename, form=FORM_UNFORMATTED, status='replace', iostat=ierr)
    if (ierr /= 0) then
      write(*, *) 'Error opening restart file: ', trim(filename)
      call release_file_unit(unit_number)
      return
    end if
    
    ! Write restart data
    write(unit_number, iostat=ierr) current_time
    write(unit_number, iostat=ierr) time_step
    write(unit_number, iostat=ierr) num_elements
    
    ! Write state variables
    do i = 1, 2*num_elements
      write(unit_number, iostat=ierr) yt(i)
      if (ierr /= 0) exit
    end do
    
    ! Write slip data
    do i = 1, num_elements
      write(unit_number, iostat=ierr) slip(i)
      if (ierr /= 0) exit
    end do
    
    close(unit_number)
    call release_file_unit(unit_number)
    
    if (ierr /= 0) then
      write(*, *) 'Error writing restart data'
    end if
    
  end subroutine write_restart_file
  
  !===============================================================================
  ! SUBROUTINE: Read restart file
  !===============================================================================
  subroutine read_restart_file(filename, current_time, time_step, yt, slip, num_elements, success)
    implicit none
    character(len=*), intent(in) :: filename
    real(DP), intent(out) :: current_time
    integer, intent(out) :: time_step
    real(DP), dimension(:), intent(out) :: yt, slip
    integer, intent(out) :: num_elements
    logical, intent(out) :: success
    
    integer :: unit_number, i, ierr, file_num_elements
    
    success = .false.
    unit_number = get_available_file_unit()
    
    open(unit_number, file=filename, form=FORM_UNFORMATTED, status='old', iostat=ierr)
    if (ierr /= 0) then
      write(*, *) 'Error opening restart file: ', trim(filename)
      call release_file_unit(unit_number)
      return
    end if
    
    ! Read restart data
    read(unit_number, iostat=ierr) current_time
    if (ierr /= 0) goto 100
    
    read(unit_number, iostat=ierr) time_step
    if (ierr /= 0) goto 100
    
    read(unit_number, iostat=ierr) file_num_elements
    if (ierr /= 0) goto 100
    
    ! Check array sizes
    if (file_num_elements /= num_elements) then
      write(*, *) 'Error: Restart file has different number of elements'
      goto 100
    end if
    
    ! Read state variables
    do i = 1, 2*num_elements
      read(unit_number, iostat=ierr) yt(i)
      if (ierr /= 0) goto 100
    end do
    
    ! Read slip data
    do i = 1, num_elements
      read(unit_number, iostat=ierr) slip(i)
      if (ierr /= 0) goto 100
    end do
    
    success = .true.
    
100 continue
    close(unit_number)
    call release_file_unit(unit_number)
    
    if (.not. success) then
      write(*, *) 'Error reading restart data'
    end if
    
  end subroutine read_restart_file
  
  !===============================================================================
  ! SUBROUTINE: Write output summary
  !===============================================================================
  subroutine write_output_summary(filename, summary_data, num_events)
    implicit none
    character(len=*), intent(in) :: filename
    real(DP), dimension(:,:), intent(in) :: summary_data
    integer, intent(in) :: num_events
    
    integer :: unit_number, i, j, ierr
    
    unit_number = get_available_file_unit()
    
    open(unit_number, file=filename, status='replace', iostat=ierr)
    if (ierr /= 0) then
      write(*, *) 'Error opening summary file: ', trim(filename)
      call release_file_unit(unit_number)
      return
    end if
    
    ! Write header
    write(unit_number, '(A)', iostat=ierr) '# Output Summary'
    write(unit_number, '(A)', iostat=ierr) '# Event Time Max_Velocity Max_Slip Area Moment'
    
    ! Write summary data
    do i = 1, num_events
      write(unit_number, OUTPUT_FORMAT_1, iostat=ierr) (summary_data(i, j), j = 1, 7)
      if (ierr /= 0) exit
    end do
    
    close(unit_number)
    call release_file_unit(unit_number)
    
  end subroutine write_output_summary
  
  !===============================================================================
  ! SUBROUTINE: Write formatted output
  !===============================================================================
  subroutine write_formatted_output(filename, data_array, format_string, num_rows, num_cols)
    implicit none
    character(len=*), intent(in) :: filename
    real(DP), dimension(:,:), intent(in) :: data_array
    character(len=*), intent(in) :: format_string
    integer, intent(in) :: num_rows, num_cols
    
    integer :: unit_number, i, j, ierr
    
    unit_number = get_available_file_unit()
    
    open(unit_number, file=filename, status='replace', iostat=ierr)
    if (ierr /= 0) then
      write(*, *) 'Error opening output file: ', trim(filename)
      call release_file_unit(unit_number)
      return
    end if
    
    ! Write formatted data
    do i = 1, num_rows
      write(unit_number, format_string, iostat=ierr) (data_array(i, j), j = 1, num_cols)
      if (ierr /= 0) exit
    end do
    
    close(unit_number)
    call release_file_unit(unit_number)
    
  end subroutine write_formatted_output
  
  !===============================================================================
  ! SUBROUTINE: Create output directory structure
  !===============================================================================
  subroutine create_output_directory_structure(base_directory)
    implicit none
    character(len=*), intent(in) :: base_directory
    
    character(len=256) :: command
    integer :: ierr
    
    ! Create base directory
    write(command, '(A,A,A)') 'mkdir -p ', trim(base_directory)
    call system(command)
    
    ! Create subdirectories
    write(command, '(A,A,A)') 'mkdir -p ', trim(base_directory), '/slip'
    call system(command)
    
    write(command, '(A,A,A)') 'mkdir -p ', trim(base_directory), '/velocity'
    call system(command)
    
    write(command, '(A,A,A)') 'mkdir -p ', trim(base_directory), '/stress'
    call system(command)
    
    write(command, '(A,A,A)') 'mkdir -p ', trim(base_directory), '/restart'
    call system(command)
    
    write(command, '(A,A,A)') 'mkdir -p ', trim(base_directory), '/summary'
    call system(command)
    
  end subroutine create_output_directory_structure
  
  !===============================================================================
  ! SUBROUTINE: Write simulation parameters to file
  !===============================================================================
  subroutine write_simulation_parameters(filename, parameters)
    implicit none
    character(len=*), intent(in) :: filename
    real(DP), dimension(:), intent(in) :: parameters
    
    integer :: unit_number, i, ierr
    
    unit_number = get_available_file_unit()
    
    open(unit_number, file=filename, status='replace', iostat=ierr)
    if (ierr /= 0) then
      write(*, *) 'Error opening parameters file: ', trim(filename)
      call release_file_unit(unit_number)
      return
    end if
    
    ! Write parameters
    do i = 1, size(parameters)
      write(unit_number, '(E15.8)', iostat=ierr) parameters(i)
      if (ierr /= 0) exit
    end do
    
    close(unit_number)
    call release_file_unit(unit_number)
    
  end subroutine write_simulation_parameters
  
  !===============================================================================
  ! SUBROUTINE: Write progress information
  !===============================================================================
  subroutine write_progress_info(current_time, total_time, time_step, num_steps)
    implicit none
    real(DP), intent(in) :: current_time, total_time
    integer, intent(in) :: time_step, num_steps
    
    real(DP) :: progress_percentage
    
    progress_percentage = (current_time / total_time) * 100.0_DP
    
    write(*, '(A,F6.2,A,I0,A,I0)') 'Progress: ', progress_percentage, '% complete, ', &
                                     time_step, ' steps, ', num_steps, ' total steps'
    
  end subroutine write_progress_info
  
  !===============================================================================
  ! SUBROUTINE: Finalize I/O handling
  !===============================================================================
  subroutine finalize_io_handling()
    implicit none
    
    ! Release all file units
    file_units_used = .false.
    
  end subroutine finalize_io_handling
  
end module io_handling
