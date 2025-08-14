!===============================================================================
! MODULE: MPI Utilities
! 
! This module provides MPI utility functions for the earthquake simulation,
! including communication functions, process management, and data distribution.
!
! Author: New module for improved MPI handling
! Last Modified: Current refactoring
!===============================================================================

module mpi_utilities
  use mpi
  use physical_constants, only: DP
  implicit none
  
  ! Make all functions public
  public
  
  ! MPI process information
  integer :: my_rank = -1
  integer :: num_processes = -1
  integer :: master_rank = 0
  
  ! MPI communicator
  integer :: world_comm = MPI_COMM_WORLD
  
  ! MPI data types
  integer :: mpi_dp_type = MPI_DOUBLE_PRECISION
  integer :: mpi_int_type = MPI_INTEGER
  integer :: mpi_logical_type = MPI_LOGICAL
  
  ! Communication buffers
  integer, parameter :: MAX_BUFFER_SIZE = 1000000
  
contains
  
  !===============================================================================
  ! SUBROUTINE: Initialize MPI utilities
  !===============================================================================
  subroutine initialize_mpi_utilities()
    implicit none
    integer :: ierr
    
    ! Get process information
    call MPI_Comm_rank(world_comm, my_rank, ierr)
    call check_mpi_error(ierr, 'MPI_Comm_rank')
    
    call MPI_Comm_size(world_comm, num_processes, ierr)
    call check_mpi_error(ierr, 'MPI_Comm_size')
    
    ! Set MPI data types based on precision
    if (DP == kind(1.0d0)) then
      mpi_dp_type = MPI_DOUBLE_PRECISION
    else if (DP == kind(1.0)) then
      mpi_dp_type = MPI_REAL
    else
      mpi_dp_type = MPI_DOUBLE_PRECISION  ! Default to double precision
    end if
    
  end subroutine initialize_mpi_utilities
  
  !===============================================================================
  ! SUBROUTINE: Check MPI error
  !===============================================================================
  subroutine check_mpi_error(ierr, operation_name)
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: operation_name
    
    if (ierr /= MPI_SUCCESS) then
      write(*, *) 'MPI error in ', trim(operation_name), ': ', ierr
      call MPI_Abort(world_comm, ierr, ierr)
      stop
    end if
    
  end subroutine check_mpi_error
  
  !===============================================================================
  ! SUBROUTINE: Get process information
  !===============================================================================
  subroutine get_process_info(rank, size, is_master)
    implicit none
    integer, intent(out) :: rank, size
    logical, intent(out) :: is_master
    
    rank = my_rank
    size = num_processes
    is_master = (my_rank == master_rank)
    
  end subroutine get_process_info
  
  !===============================================================================
  ! SUBROUTINE: Check if current process is master
  !===============================================================================
  function is_master_process() result(is_master)
    implicit none
    logical :: is_master
    
    is_master = (my_rank == master_rank)
    
  end function is_master_process
  
  !===============================================================================
  ! SUBROUTINE: Synchronize all processes
  !===============================================================================
  subroutine synchronize_processes()
    implicit none
    integer :: ierr
    
    call MPI_Barrier(world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Barrier')
    
  end subroutine synchronize_processes
  
  !===============================================================================
  ! SUBROUTINE: Broadcast scalar integer
  !===============================================================================
  subroutine broadcast_integer_scalar(value, root)
    implicit none
    integer, intent(inout) :: value
    integer, intent(in) :: root
    integer :: ierr
    
    call MPI_Bcast(value, 1, mpi_int_type, root, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Bcast integer')
    
  end subroutine broadcast_integer_scalar
  
  !===============================================================================
  ! SUBROUTINE: Broadcast scalar double precision
  !===============================================================================
  subroutine broadcast_double_scalar(value, root)
    implicit none
    real(DP), intent(inout) :: value
    integer, intent(in) :: root
    integer :: ierr
    
    call MPI_Bcast(value, 1, mpi_dp_type, root, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Bcast double')
    
  end subroutine broadcast_double_scalar
  
  !===============================================================================
  ! SUBROUTINE: Broadcast integer array
  !===============================================================================
  subroutine broadcast_integer_array(array, count, root)
    implicit none
    integer, dimension(:), intent(inout) :: array
    integer, intent(in) :: count, root
    integer :: ierr
    
    call MPI_Bcast(array, count, mpi_int_type, root, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Bcast integer array')
    
  end subroutine broadcast_integer_array
  
  !===============================================================================
  ! SUBROUTINE: Broadcast double precision array
  !===============================================================================
  subroutine broadcast_double_array(array, count, root)
    implicit none
    real(DP), dimension(:), intent(inout) :: array
    integer, intent(in) :: count, root
    integer :: ierr
    
    call MPI_Bcast(array, count, mpi_dp_type, root, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Bcast double array')
    
  end subroutine broadcast_double_array
  
  !===============================================================================
  ! SUBROUTINE: Gather double precision arrays
  !===============================================================================
  subroutine gather_double_arrays(send_buffer, recv_buffer, send_count, root)
    implicit none
    real(DP), dimension(:), intent(in) :: send_buffer
    real(DP), dimension(:), intent(out) :: recv_buffer
    integer, intent(in) :: send_count, root
    integer :: ierr
    
    call MPI_Gather(send_buffer, send_count, mpi_dp_type, recv_buffer, &
                    send_count, mpi_dp_type, root, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Gather double arrays')
    
  end subroutine gather_double_arrays
  
  !===============================================================================
  ! SUBROUTINE: Allgather double precision arrays
  !===============================================================================
  subroutine allgather_double_arrays(send_buffer, recv_buffer, send_count)
    implicit none
    real(DP), dimension(:), intent(in) :: send_buffer
    real(DP), dimension(:), intent(out) :: recv_buffer
    integer, intent(in) :: send_count
    integer :: ierr
    
    call MPI_Allgather(send_buffer, send_count, mpi_dp_type, recv_buffer, &
                       send_count, mpi_dp_type, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Allgather double arrays')
    
  end subroutine allgather_double_arrays
  
  !===============================================================================
  ! SUBROUTINE: Scatter double precision arrays
  !===============================================================================
  subroutine scatter_double_arrays(send_buffer, recv_buffer, send_count, root)
    implicit none
    real(DP), dimension(:), intent(in) :: send_buffer
    real(DP), dimension(:), intent(out) :: recv_buffer
    integer, intent(in) :: send_count, root
    integer :: ierr
    
    call MPI_Scatter(send_buffer, send_count, mpi_dp_type, recv_buffer, &
                     send_count, mpi_dp_type, root, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Scatter double arrays')
    
  end subroutine scatter_double_arrays
  
  !===============================================================================
  ! SUBROUTINE: Reduce double precision scalar
  !===============================================================================
  subroutine reduce_double_scalar(send_value, recv_value, operation, root)
    implicit none
    real(DP), intent(in) :: send_value
    real(DP), intent(out) :: recv_value
    integer, intent(in) :: operation, root
    integer :: ierr
    
    call MPI_Reduce(send_value, recv_value, 1, mpi_dp_type, operation, root, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Reduce double scalar')
    
  end subroutine reduce_double_scalar
  
  !===============================================================================
  ! SUBROUTINE: Allreduce double precision scalar
  !===============================================================================
  subroutine allreduce_double_scalar(send_value, recv_value, operation)
    implicit none
    real(DP), intent(in) :: send_value
    real(DP), intent(out) :: recv_value
    integer, intent(in) :: operation
    integer :: ierr
    
    call MPI_Allreduce(send_value, recv_value, 1, mpi_dp_type, operation, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Allreduce double scalar')
    
  end subroutine allreduce_double_scalar
  
  !===============================================================================
  ! SUBROUTINE: Reduce double precision array
  !===============================================================================
  subroutine reduce_double_array(send_buffer, recv_buffer, count, operation, root)
    implicit none
    real(DP), dimension(:), intent(in) :: send_buffer
    real(DP), dimension(:), intent(out) :: recv_buffer
    integer, intent(in) :: count, operation, root
    integer :: ierr
    
    call MPI_Reduce(send_buffer, recv_buffer, count, mpi_dp_type, operation, root, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Reduce double array')
    
  end subroutine reduce_double_array
  
  !===============================================================================
  ! SUBROUTINE: Allreduce double precision array
  !===============================================================================
  subroutine allreduce_double_array(send_buffer, recv_buffer, count, operation)
    implicit none
    real(DP), dimension(:), intent(in) :: send_buffer
    real(DP), dimension(:), intent(out) :: recv_buffer
    integer, intent(in) :: count, operation
    integer :: ierr
    
    call MPI_Allreduce(send_buffer, recv_buffer, count, mpi_dp_type, operation, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Allreduce double array')
    
  end subroutine allreduce_double_array
  
  !===============================================================================
  ! SUBROUTINE: Send double precision array
  !===============================================================================
  subroutine send_double_array(buffer, count, dest, tag)
    implicit none
    real(DP), dimension(:), intent(in) :: buffer
    integer, intent(in) :: count, dest, tag
    integer :: ierr
    
    call MPI_Send(buffer, count, mpi_dp_type, dest, tag, world_comm, ierr)
    call check_mpi_error(ierr, 'MPI_Send double array')
    
  end subroutine send_double_array
  
  !===============================================================================
  ! SUBROUTINE: Receive double precision array
  !===============================================================================
  subroutine receive_double_array(buffer, count, source, tag)
    implicit none
    real(DP), dimension(:), intent(out) :: buffer
    integer, intent(in) :: count, source, tag
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    
    call MPI_Recv(buffer, count, mpi_dp_type, source, tag, world_comm, status, ierr)
    call check_mpi_error(ierr, 'MPI_Recv double array')
    
  end subroutine receive_double_array
  
  !===============================================================================
  ! SUBROUTINE: Non-blocking send double precision array
  !===============================================================================
  subroutine isend_double_array(buffer, count, dest, tag, request)
    implicit none
    real(DP), dimension(:), intent(in) :: buffer
    integer, intent(in) :: count, dest, tag
    integer, intent(out) :: request
    integer :: ierr
    
    call MPI_Isend(buffer, count, mpi_dp_type, dest, tag, world_comm, request, ierr)
    call check_mpi_error(ierr, 'MPI_Isend double array')
    
  end subroutine isend_double_array
  
  !===============================================================================
  ! SUBROUTINE: Non-blocking receive double precision array
  !===============================================================================
  subroutine ireceive_double_array(buffer, count, source, tag, request)
    implicit none
    real(DP), dimension(:), intent(out) :: buffer
    integer, intent(in) :: count, source, tag
    integer, intent(out) :: request
    integer :: ierr
    
    call MPI_Irecv(buffer, count, mpi_dp_type, source, tag, world_comm, request, ierr)
    call check_mpi_error(ierr, 'MPI_Irecv double array')
    
  end subroutine ireceive_double_array
  
  !===============================================================================
  ! SUBROUTINE: Wait for non-blocking operation
  !===============================================================================
  subroutine wait_for_request(request)
    implicit none
    integer, intent(inout) :: request
    integer :: ierr, status(MPI_STATUS_SIZE)
    
    call MPI_Wait(request, status, ierr)
    call check_mpi_error(ierr, 'MPI_Wait')
    
  end subroutine wait_for_request
  
  !===============================================================================
  ! SUBROUTINE: Wait for all non-blocking operations
  !===============================================================================
  subroutine wait_for_all_requests(requests, num_requests)
    implicit none
    integer, dimension(:), intent(inout) :: requests
    integer, intent(in) :: num_requests
    integer :: ierr, statuses(MPI_STATUS_SIZE, num_requests)
    
    call MPI_Waitall(num_requests, requests, statuses, ierr)
    call check_mpi_error(ierr, 'MPI_Waitall')
    
  end subroutine wait_for_all_requests
  
  !===============================================================================
  ! SUBROUTINE: Print process information
  !===============================================================================
  subroutine print_process_info()
    implicit none
    
    if (is_master_process()) then
      write(*, *) '=== MPI Process Information ==='
      write(*, *) 'Total processes: ', num_processes
      write(*, *) 'Master rank: ', master_rank
      write(*, *) 'MPI data types:'
      write(*, *) '  Double precision: ', mpi_dp_type
      write(*, *) '  Integer: ', mpi_int_type
      write(*, *) '  Logical: ', mpi_logical_type
      write(*, *) '================================'
    end if
    
  end subroutine print_process_info
  
  !===============================================================================
  ! SUBROUTINE: Finalize MPI utilities
  !===============================================================================
  subroutine finalize_mpi_utilities()
    implicit none
    
    ! Reset process information
    my_rank = -1
    num_processes = -1
    
  end subroutine finalize_mpi_utilities
  
  !===============================================================================
  ! SUBROUTINE: Abort MPI execution
  !===============================================================================
  subroutine abort_mpi_execution(error_code)
    implicit none
    integer, intent(in) :: error_code
    integer :: ierr
    
    call MPI_Abort(world_comm, error_code, ierr)
    
  end subroutine abort_mpi_execution
  
end module mpi_utilities
