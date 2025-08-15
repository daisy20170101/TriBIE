! Alternative version without separate precision module
! This version defines DP directly in the main program

program p_calc_green
    use m_calc_green
    use mpi
    use omp_lib
   implicit none
   
   ! Define precision parameter directly in main program
   integer, parameter :: DP = kind(1.0d0)
   
   character(*),parameter        :: fname="triangular_mesh.gts"
   integer ::                   n_vertex,n_edge,n_cell
   real(DP),DIMENSION(:,:),ALLOCATABLE  ::  arr_vertex
   integer,DIMENSION(:,:),ALLOCATABLE  ::  arr_edge
   integer,DIMENSION(:,:),ALLOCATABLE  ::  arr_cell
  
   integer :: ierr,size,myid
   integer :: Nt,Nt_all,master
   real(DP) :: start_time, end_time
   integer :: local_cells,cells_processed
   
   ! Dynamic load balancing variables
   integer :: base_cells, extra_cells, start_idx
   
   ! OpenMP setup
   integer :: num_threads
   
   ! Error handling
   logical :: error_occurred
   character(len=256) :: error_message
   
   ! Initialize error flag
   error_occurred = .false.
   error_message = ""
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Load input mesh data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call MPI_Init(ierr)
   if (ierr /= MPI_SUCCESS) then
     write(*,*) "Error: Failed to initialize MPI"
     stop
   end if
   
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
   master=0
   
   ! Set OpenMP threads per MPI process
   num_threads = omp_get_max_threads()
   call omp_set_num_threads(num_threads)
   
   if (myid == master) then
     start_time = MPI_Wtime()
     write(*,*) "Starting TriGreen calculation with improved parallelization"
     write(*,*) "MPI processes:", size
     write(*,*) "OpenMP threads per process:", num_threads
     write(*,*) "Total parallel threads:", size * num_threads
   endif
    
   ! Error handling wrapper
   if (.not. error_occurred) then
     !!! load n_vertex,n_vertex,n_cell from fname
     call load_name(fname,n_vertex,n_edge,n_cell, error_occurred, error_message)
   end if
   
   if (.not. error_occurred) then
     allocate(arr_vertex(n_vertex,3),arr_edge(n_edge,2),arr_cell(n_cell,3), stat=ierr)
     if (ierr /= 0) then
       error_occurred = .true.
       error_message = "Failed to allocate arrays"
     end if
   end if

   if (.not. error_occurred) then
     !!!!! load arr_vertex,arr_edge,arr_cell from fname
     call load_gts(fname,n_vertex,n_edge,&
                n_cell,arr_vertex,arr_edge,arr_cell, error_occurred, error_message)
   end if
   
   if (.not. error_occurred) then
     write(*,*) myid

             ! OPTIMIZATION: Implement dynamic load balancing for optimal distribution
       Nt_all = n_cell
       
       ! Special handling for single-CPU execution
       if (size == 1) then
          write(*,*)'Single CPU execution detected - all cells assigned to process 0'
          base_cells = n_cell
          extra_cells = 0
          local_cells = n_cell
          start_idx = 0
       else
          ! Dynamic load balancing for multi-CPU execution
          base_cells = n_cell / size
          extra_cells = mod(n_cell, size)
          
          if (myid < extra_cells) then
             local_cells = base_cells + 1
             start_idx = myid * local_cells
          else
             local_cells = base_cells
             start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells
          end if
       end if
       
       cells_processed = local_cells
       
       write(*,*) "Process", myid, "assigned", local_cells, "cells starting from", start_idx
       
       ! Call the main computation subroutine
       call calc_green_allcell_improved(myid, size, Nt_all, arr_vertex, arr_cell, &
                n_vertex, n_cell, cells_processed, base_cells, extra_cells, error_occurred, error_message)
   end if
   
   ! Final error check
   if (error_occurred) then
     write(*,*) "Process", myid, "encountered an error: ", trim(error_message)
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
   end if
   
   ! Performance monitoring
   if (myid == master) then
     end_time = MPI_Wtime()
     call performance_monitoring(myid, start_time, end_time, n_cell)
   end if
   
   ! Cleanup
   deallocate(arr_vertex, arr_edge, arr_cell)
   
   call MPI_Finalize(ierr)
   
contains

! Include all the subroutines here with DP parameter passed from main program
! This avoids the need for a separate precision module

subroutine calc_ss_ds(v1, v2, v3, v_pl, ss, ds, op)
  use m_calc_green 
  implicit none
  
  ! Pass DP from main program
  integer, parameter :: DP = kind(1.0d0)
  
  real(DP), intent(in) ::  v1(3), v2(3), v3(3)
  real(DP), intent(in) ::  v_pl(3)
  real(DP), intent(out) :: ss, ds, op
  
  real(DP) ::  vpl, vpl0, theta
  real(DP) ::  nv(3)
  real(DP) ::  x, y, z
    
  real(DP) ::  a, b
  real(DP) ::  rl,a11,a12,a13,gamma
  integer ::  i
    
  ! get triangle's normal vector
  call tri_normal_vect(v1, v2, v3, nv)

  if ( nv(3) .lt. 0 ) then
    do i=1, 3
      nv(i) = -nv(i)
    end do
  end if

  ! Handle vertical triangle (nv(3) = 0)
  if (abs(nv(3)) .lt. 1.0d-12) then
    write(*,*) "INFO: Vertical triangle detected in calc_ss_ds - using predefined axes"
    write(*,*) "Triangle vertices: v1=", v1, " v2=", v2, " v3=", v3
    write(*,*) "Normal vector: nv=", nv
    
    ! For vertical triangles, use predefined coordinate system
    ! Set strike direction along x-axis, dip direction along y-axis
    ! This follows the convention for vertical fault planes
    ss = vpl1  ! Strike-slip component
    ds = vpl2  ! Dip-slip component  
    op = 0.0d0 ! Opening component (no opening for vertical faults)
    return
  end if

  ! Check for valid vpl1 and vpl2 parameters
  if (isnan(vpl1) .or. isnan(vpl2)) then
    write(*,*) "ERROR: Invalid vpl1 or vpl2 in calc_ss_ds"
    write(*,*) "vpl1 =", vpl1, "vpl2 =", vpl2
    ss = 0.0d0
    ds = 0.0d0
    op = 0.0d0
    return
  end if

  rl=(vpl1**2+vpl2**2)+(nv(1)*vpl1+nv(2)*vpl2)**2/nv(3)**2
  rl=1.d0/rl
  rl=sqrt(rl)
  gamma=-(nv(1)*vpl1+nv(2)*vpl2)*rl/nv(3)

  a11 = vpl1*rl
  a12 = vpl2*rl
  a13 = gamma

  ! calc ss, ds
  ! Check for horizontal triangles (nv(1) and nv(2) both zero)
  if (abs(nv(1)) < 1.0d-12 .and. abs(nv(2)) < 1.0d-12) then
    write(*,*) "WARNING: Horizontal triangle detected in calc_ss_ds"
    write(*,*) "Triangle vertices: v1=", v1, " v2=", v2, " v3=", v3
    write(*,*) "Normal vector: nv=", nv
    write(*,*) "This triangle will be skipped"
    ss = 0.0d0
    ds = 0.0d0
    op = 0.0d0
    return
  end if
  
  ss =(nv(2)*a11-nv(1)*a12)/sqrt(nv(1)**2+nv(2)**2)
  ds =sqrt(nv(1)**2+nv(2)**2-(nv(2)*a11-nv(1)*a12)**2)
  ds =ds/sqrt(nv(1)**2+nv(2)**2)
  op =0.d0
    
end subroutine calc_ss_ds

! Add other subroutines here...

end program p_calc_green
