! Test program for dynamic load balancing algorithm
program test_load_balancing
  implicit none
  
  integer :: n_cell, size, base_cells, extra_cells
  integer :: i, j, local_cells, start_idx
  
  ! Test cases
  integer, parameter :: test_cases = 3
  integer, dimension(test_cases) :: test_n_cell = [100, 101, 102]
  integer, dimension(test_cases) :: test_size = [4, 4, 4]
  
  write(*,*) "Testing Dynamic Load Balancing Algorithm"
  write(*,*) "========================================"
  
  do i = 1, test_cases
    n_cell = test_n_cell(i)
    size = test_size(i)
    
    write(*,*) ""
    write(*,*) "Test case", i, ": n_cell =", n_cell, ", size =", size
    
    if(mod(n_cell,size)/=0)then
      write(*,*)'n_cell (',n_cell,') is not evenly divisible by MPI processes (',size,')'
      write(*,*)'Implementing dynamic load balancing for optimal distribution...'
      
      ! Calculate optimal distribution using ceiling division
      base_cells = n_cell / size
      extra_cells = mod(n_cell, size)
      
      write(*,*)'Base cells per process:', base_cells
      write(*,*)'Extra cells to distribute:', extra_cells
      write(*,*)'Processes 0 to', extra_cells-1, 'will get', base_cells+1, 'cells'
      write(*,*)'Processes', extra_cells, 'to', size-1, 'will get', base_cells, 'cells'
      
      ! Show distribution for each process
      write(*,*) ""
      write(*,*) "Process distribution:"
      do j = 0, size-1
        if (j < extra_cells) then
          local_cells = base_cells + 1
          start_idx = j * (base_cells + 1)
        else
          local_cells = base_cells
          start_idx = extra_cells * (base_cells + 1) + (j - extra_cells) * base_cells
        end if
        write(*,*) "Process", j, ": cells =", local_cells, "start_idx =", start_idx, &
                   "end_idx =", start_idx + local_cells - 1
      end do
      
    else
      write(*,*)'n_cell (',n_cell,') is evenly divisible by MPI processes (',size,')'
      base_cells = n_cell / size
      extra_cells = 0
      write(*,*)'Each process calculates', base_cells, 'cells'
    end if
    
    write(*,*) "========================================"
  end do
  
  write(*,*) ""
  write(*,*) "Load balancing test completed!"
  
end program test_load_balancing
