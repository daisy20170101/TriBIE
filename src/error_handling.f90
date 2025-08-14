!===============================================================================
! MODULE: Error Handling
! 
! This module provides comprehensive error handling and validation functions
! for the earthquake simulation, including MPI error checking and user-defined errors.
!
! Author: New module for improved error handling
! Last Modified: Current refactoring
!===============================================================================

module error_handling
  use mpi
  implicit none
  
  ! Make all functions public
  public
  
  ! Error codes
  integer, parameter :: ERROR_SUCCESS = 0
  integer, parameter :: ERROR_MPI_FAILURE = 1
  integer, parameter :: ERROR_FILE_IO = 2
  integer, parameter :: ERROR_INVALID_PARAMETER = 3
  integer, parameter :: ERROR_MEMORY_ALLOCATION = 4
  integer, parameter :: ERROR_NUMERICAL = 5
  integer, parameter :: ERROR_USER_DEFINED = 100
  
  ! Error message buffer size
  integer, parameter :: MAX_ERROR_MESSAGE_LENGTH = 512
  
  ! Global error state
  logical :: error_occurred = .false.
  integer :: last_error_code = ERROR_SUCCESS
  character(len=MAX_ERROR_MESSAGE_LENGTH) :: last_error_message = ''
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Initialize error handling
  !===============================================================================
  subroutine initialize_error_handling()
    implicit none
    
    error_occurred = .false.
    last_error_code = ERROR_SUCCESS
    last_error_message = ''
    
  end subroutine initialize_error_handling
  
  !===============================================================================
  ! SUBROUTINE: Check MPI error
  !===============================================================================
  subroutine check_mpi_error(ierr, operation_name)
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: operation_name
    
    if (ierr /= MPI_SUCCESS) then
      call handle_mpi_error(ierr, operation_name)
    end if
    
  end subroutine check_mpi_error
  
  !===============================================================================
  ! SUBROUTINE: Handle MPI error
  !===============================================================================
  subroutine handle_mpi_error(ierr, operation_name)
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: operation_name
    character(len=MAX_ERROR_MESSAGE_LENGTH) :: error_msg
    integer :: error_class, error_string_len
    character(len=MPI_MAX_ERROR_STRING) :: error_string
    
    ! Get MPI error class and string
    call MPI_Error_class(ierr, error_class)
    call MPI_Error_string(ierr, error_string, error_string_len)
    
    ! Build error message
    write(error_msg, '(A,A,A,I0,A,A)') 'MPI error in ', trim(operation_name), &
                                       ': code=', ierr, ', class=', trim(error_string)
    
    ! Set global error state
    error_occurred = .true.
    last_error_code = ERROR_MPI_FAILURE
    last_error_message = error_msg
    
    ! Print error message
    write(*, *) 'ERROR: ', trim(error_msg)
    
    ! Abort MPI execution
    call MPI_Abort(MPI_COMM_WORLD, ERROR_MPI_FAILURE, ierr)
    
  end subroutine handle_mpi_error
  
  !===============================================================================
  ! SUBROUTINE: Handle general error
  !===============================================================================
  subroutine handle_error(error_message, error_code)
    implicit none
    character(len=*), intent(in) :: error_message
    integer, intent(in), optional :: error_code
    
    integer :: code
    
    ! Set error code
    if (present(error_code)) then
      code = error_code
    else
      code = ERROR_USER_DEFINED
    end if
    
    ! Set global error state
    error_occurred = .true.
    last_error_code = code
    last_error_message = error_message
    
    ! Print error message
    write(*, *) 'ERROR: ', trim(error_message)
    
    ! Stop execution
    stop
    
  end subroutine handle_error
  
  !===============================================================================
  ! SUBROUTINE: Handle file I/O error
  !===============================================================================
  subroutine handle_file_error(filename, operation, iostat)
    implicit none
    character(len=*), intent(in) :: filename, operation
    integer, intent(in) :: iostat
    
    character(len=MAX_ERROR_MESSAGE_LENGTH) :: error_msg
    
    ! Build error message
    write(error_msg, '(A,A,A,A,A,I0)') 'File I/O error: ', trim(operation), &
                                       ' failed for file "', trim(filename), &
                                       '", iostat=', iostat
    
    ! Handle error
    call handle_error(error_msg, ERROR_FILE_IO)
    
  end subroutine handle_file_error
  
  !===============================================================================
  ! SUBROUTINE: Handle parameter validation error
  !===============================================================================
  subroutine handle_parameter_error(parameter_name, expected_value, actual_value)
    implicit none
    character(len=*), intent(in) :: parameter_name
    character(len=*), intent(in) :: expected_value
    character(len=*), intent(in) :: actual_value
    
    character(len=MAX_ERROR_MESSAGE_LENGTH) :: error_msg
    
    ! Build error message
    write(error_msg, '(A,A,A,A,A,A)') 'Parameter validation error: ', trim(parameter_name), &
                                       ' should be ', trim(expected_value), &
                                       ', but got ', trim(actual_value)
    
    ! Handle error
    call handle_error(error_msg, ERROR_INVALID_PARAMETER)
    
  end subroutine handle_parameter_error
  
  !===============================================================================
  ! SUBROUTINE: Handle memory allocation error
  !===============================================================================
  subroutine handle_memory_error(array_name, size, iostat)
    implicit none
    character(len=*), intent(in) :: array_name
    integer, intent(in) :: size
    integer, intent(in) :: iostat
    
    character(len=MAX_ERROR_MESSAGE_LENGTH) :: error_msg
    
    ! Build error message
    write(error_msg, '(A,A,A,I0,A,I0)') 'Memory allocation error: failed to allocate ', &
                                       trim(array_name), ' with size ', size, &
                                       ', iostat=', iostat
    
    ! Handle error
    call handle_error(error_msg, ERROR_MEMORY_ALLOCATION)
    
  end subroutine handle_memory_error
  
  !===============================================================================
  ! SUBROUTINE: Handle numerical error
  !===============================================================================
  subroutine handle_numerical_error(operation, value, threshold)
    implicit none
    character(len=*), intent(in) :: operation
    real(8), intent(in) :: value, threshold
    
    character(len=MAX_ERROR_MESSAGE_LENGTH) :: error_msg
    
    ! Build error message
    write(error_msg, '(A,A,A,E15.7,A,E15.7)') 'Numerical error: ', trim(operation), &
                                       ' failed with value ', value, &
                                       ', threshold ', threshold
    
    ! Handle error
    call handle_error(error_msg, ERROR_NUMERICAL)
    
  end subroutine handle_numerical_error
  
  !===============================================================================
  ! SUBROUTINE: Check for NaN or infinite values
  !===============================================================================
  subroutine check_numerical_stability(array, array_name)
    implicit none
    real(8), dimension(:), intent(in) :: array
    character(len=*), intent(in) :: array_name
    
    integer :: i, n
    logical :: has_nan, has_inf
    
    n = size(array)
    has_nan = .false.
    has_inf = .false.
    
    ! Check for NaN and infinite values
    do i = 1, n
      if (isnan(array(i))) has_nan = .true.
      if (isinf(array(i))) has_inf = .true.
    end do
    
    ! Report issues
    if (has_nan) then
      write(*, *) 'WARNING: NaN values detected in array ', trim(array_name)
    end if
    
    if (has_inf) then
      write(*, *) 'WARNING: Infinite values detected in array ', trim(array_name)
    end if
    
  end subroutine check_numerical_stability
  
  !===============================================================================
  ! SUBROUTINE: Check array bounds
  !===============================================================================
  subroutine check_array_bounds(array, array_name, expected_size)
    implicit none
    real(8), dimension(:), intent(in) :: array
    character(len=*), intent(in) :: array_name
    integer, intent(in) :: expected_size
    
    integer :: actual_size
    
    actual_size = size(array)
    
    if (actual_size /= expected_size) then
      call handle_parameter_error(trim(array_name) // ' size', &
                                'expected ' // trim(int_to_string(expected_size)), &
                                'actual ' // trim(int_to_string(actual_size)))
    end if
    
  end subroutine check_array_bounds
  
  !===============================================================================
  ! SUBROUTINE: Check positive value
  !===============================================================================
  subroutine check_positive_value(value, value_name)
    implicit none
    real(8), intent(in) :: value
    character(len=*), intent(in) :: value_name
    
    if (value <= 0.0_8) then
      call handle_parameter_error(trim(value_name), 'positive', &
                                trim(real_to_string(value)))
    end if
    
  end subroutine check_positive_value
  
  !===============================================================================
  ! SUBROUTINE: Check value in range
  !===============================================================================
  subroutine check_value_in_range(value, value_name, min_val, max_val)
    implicit none
    real(8), intent(in) :: value, min_val, max_val
    character(len=*), intent(in) :: value_name
    
    if (value < min_val .or. value > max_val) then
      call handle_parameter_error(trim(value_name), &
                                'in range [' // trim(real_to_string(min_val)) // &
                                ', ' // trim(real_to_string(max_val)) // ']', &
                                trim(real_to_string(value)))
    end if
    
  end subroutine check_value_in_range
  
  !===============================================================================
  ! SUBROUTINE: Convert integer to string
  !===============================================================================
  function int_to_string(i) result(str)
    implicit none
    integer, intent(in) :: i
    character(len=20) :: str
    
    write(str, '(I0)') i
    
  end function int_to_string
  
  !===============================================================================
  ! SUBROUTINE: Convert real to string
  !===============================================================================
  function real_to_string(r) result(str)
    implicit none
    real(8), intent(in) :: r
    character(len=20) :: str
    
    write(str, '(E15.7)') r
    
  end function real_to_string
  
  !===============================================================================
  ! SUBROUTINE: Get error status
  !===============================================================================
  function has_error_occurred() result(has_error)
    implicit none
    logical :: has_error
    
    has_error = error_occurred
    
  end function has_error_occurred
  
  !===============================================================================
  ! SUBROUTINE: Get last error code
  !===============================================================================
  function get_last_error_code() result(error_code)
    implicit none
    integer :: error_code
    
    error_code = last_error_code
    
  end function get_last_error_code
  
  !===============================================================================
  ! SUBROUTINE: Get last error message
  !===============================================================================
  function get_last_error_message() result(error_message)
    implicit none
    character(len=MAX_ERROR_MESSAGE_LENGTH) :: error_message
    
    error_message = last_error_message
    
  end function get_last_error_message
  
  !===============================================================================
  ! SUBROUTINE: Clear error state
  !===============================================================================
  subroutine clear_error_state()
    implicit none
    
    error_occurred = .false.
    last_error_code = ERROR_SUCCESS
    last_error_message = ''
    
  end subroutine clear_error_state
  
  !===============================================================================
  ! SUBROUTINE: Print error summary
  !===============================================================================
  subroutine print_error_summary()
    implicit none
    
    if (error_occurred) then
      write(*, *) '=== Error Summary ==='
      write(*, *) 'Error occurred: Yes'
      write(*, *) 'Error code: ', last_error_code
      write(*, *) 'Error message: ', trim(last_error_message)
      write(*, *) '===================='
    else
      write(*, *) 'No errors occurred'
    end if
    
  end subroutine print_error_summary
  
end module error_handling
