!$FREEFORM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Calculate mesh's green coef
!
!   Version:
!       2024-12-19  Improved parallelization with dynamic load balancing
!       2007-09-18  Clean this program
!       2006-12-13  Create this program
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program p_calc_green
    use m_calc_green
    use mpi
    use omp_lib
   implicit none
   
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
          cells_processed = n_cell
          Nt = n_cell
          
          ! Call subroutine directly for single CPU
          call calc_green_allcell_improved(myid,size,Nt,arr_vertex,arr_cell,& 
                  n_vertex,n_cell,cells_processed,base_cells,extra_cells, error_occurred, error_message)
       else
          ! Multi-CPU dynamic load balancing
          if(mod(n_cell,size)/=0)then
             write(*,*)'n_cell (',n_cell,') is not evenly divisible by MPI processes (',size,')'
             write(*,*)'Implementing dynamic load balancing for optimal distribution...'
             
             ! Calculate optimal distribution using ceiling division
             base_cells = Nt_all / size
             extra_cells = mod(Nt_all, size)
             
             write(*,*)'Base cells per process:', base_cells
             write(*,*)'Extra cells to distribute:', extra_cells
             write(*,*)'Processes 0 to', extra_cells-1, 'will get', base_cells+1, 'cells'
             write(*,*)'Processes', extra_cells, 'to', size-1, 'will get', base_cells, 'cells'
          else
             write(*,*)'n_cell (',n_cell,') is evenly divisible by MPI processes (',size,')'
             base_cells = Nt_all / size
             extra_cells = 0
             write(*,*)'Each process calculates', base_cells, 'cells'
          end if

          ! Calculate local cell count for this process using optimal distribution
          ! Ensure master process (myid == 0) always gets at least 1 cell
          if (myid < extra_cells) then
             local_cells = base_cells + 1
             start_idx = myid * (base_cells + 1)
          else
             local_cells = base_cells
             start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells
          end if
          
          ! Special handling for master process to ensure it always gets cells
          if (myid == 0 .and. local_cells == 0) then
             ! If master would get 0 cells, take 1 from the last process
             local_cells = 1
             start_idx = 0
             write(*,*) 'Warning: Master process adjusted to get 1 cell'
          end if
          
          write(*,*)'Process', myid, 'gets', local_cells, 'cells starting from index', start_idx
          cells_processed = local_cells
          Nt = local_cells

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! for all element's center point as op
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call calc_green_allcell_improved(myid,size,Nt,arr_vertex,arr_cell,& 
                 n_vertex,n_cell,cells_processed,base_cells,extra_cells, error_occurred, error_message)
       end if
     end if  ! End of if (.not. error_occurred) block

   
   ! Clean up arrays
   if (allocated(arr_vertex)) deallocate(arr_vertex, stat=ierr)
   if (allocated(arr_edge)) deallocate(arr_edge, stat=ierr)
   if (allocated(arr_cell)) deallocate(arr_cell, stat=ierr)
 
   ! Check for errors and report
   if (error_occurred) then
     write(*,*) "Error occurred: ", trim(error_message)
     write(*,*) "Process", myid, "encountered an error"
   else
     call MPI_barrier(MPI_COMM_WORLD,ierr)
     
     ! Performance monitoring
     if (myid == master) then
       end_time = MPI_Wtime()
       call performance_monitoring(myid, start_time, end_time, cells_processed)
     endif
   end if
   
   ! Always finalize MPI
   call MPI_finalize(ierr)
   
   if (error_occurred) then
     stop 1
   end if
end program


subroutine load_name(fname,n_vertex,n_edge,n_cell, error_occurred, error_message)
    implicit none
   integer, parameter :: DP=kind(1.d0)

    character(*), intent(in) ::  fname
    integer, intent(out) :: n_vertex,n_edge,n_cell
    logical, intent(inout) :: error_occurred
    character(len=*), intent(inout) :: error_message
    
    integer :: i, iostat

    OPEN(unit=33, FILE=fname, status='old', iostat=iostat)
    if (iostat /= 0) then
      error_occurred = .true.
      error_message = "Failed to open file: " // trim(fname)
      return
    end if
    
    read(33, *, iostat=iostat) n_vertex, n_edge, n_cell
    if (iostat /= 0) then
      error_occurred = .true.
      error_message = "Failed to read mesh dimensions from file"
      close(33)
      return
    end if
    
    close(33)
    
    ! Validate dimensions
    if (n_vertex <= 0 .or. n_edge < 0 .or. n_cell <= 0) then
      error_occurred = .true.
      error_message = "Invalid mesh dimensions"
      return
    end if
    
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine load_gts(fname,n_vertex,n_edge,n_cell,arr_vertex,arr_edge,&  
                arr_cell, error_occurred, error_message)
    implicit none
integer, parameter :: DP=kind(1.d0)
    
    character(*), intent(in) ::  fname
    integer  :: n_vertex,n_edge,n_cell
    real(DP), intent(out) :: arr_vertex(n_vertex,3)
    integer, intent(out) ::  arr_edge(n_edge,2)
    integer, intent(out) ::  arr_cell(n_cell,3)
    logical, intent(inout) :: error_occurred
    character(len=*), intent(inout) :: error_message
 
    integer :: i, iostat

    OPEN(unit=10, FILE=fname, status='old', iostat=iostat)
    if (iostat /= 0) then
      error_occurred = .true.
      error_message = "Failed to open GTS file: " // trim(fname)
      return
    end if
    
    read(10, *, iostat=iostat) n_vertex, n_edge, n_cell
    if (iostat /= 0) then
      error_occurred = .true.
      error_message = "Failed to read GTS file header"
      close(10)
      return
    end if
    
    write(*,*) ">>> load_gts"
    write(*,*) "    fname   =", fname
    write(*,*) "    n_vertex=", n_vertex
    write(*,*) "    n_edge  =", n_edge
    write(*,*) "    n_cell  =", n_cell
   
    ! Read vertex coordinates
    do i=1, n_vertex
        read(10, *, iostat=iostat) arr_vertex(i, 1), arr_vertex(i, 2), arr_vertex(i, 3)
        if (iostat /= 0) then
          error_occurred = .true.
          error_message = "Failed to read vertex data"
          close(10)
          return
        end if
    end do
    
    
    ! Read cell definitions
    do i=1, n_cell
        read(10, *, iostat=iostat) arr_cell(i, 1), arr_cell(i, 2), arr_cell(i, 3)
        if (iostat /= 0) then
          error_occurred = .true.
          error_message = "Failed to read cell data"
          close(10)
          return
        end if
    end do
    
    CLOSE(10)
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind(1.d0)) function calc_triangle_area(p1, p2, p3)
    implicit none
 integer, parameter :: DP=kind(1.d0)
    
    real(DP)                     ::  p1(3), p2(3), p3(3)
    
    real(DP)                     ::  a, b, c, s
    real(DP)                     ::  dx, dy, dz
  
    ! a
    dx = p1(1) - p2(1)
    dy = p1(2) - p2(2)
    dz = p1(3) - p2(3)
    a = sqrt(dx**2 + dy**2 + dz**2)
    
    ! b
    dx = p2(1) - p3(1)
    dy = p2(2) - p3(2)
    dz = p2(3) - p3(3)
    b = sqrt(dx**2 + dy**2 + dz**2)
    
    ! c
    dx = p3(1) - p1(1)
    dy = p3(2) - p1(2)
    dz = p3(3) - p1(3)
    c = sqrt(dx**2 + dy**2 + dz**2)
    
    s = (a + b + c) / 2.0
    calc_triangle_area = sqrt(s*(s-a)*(s-b)*(s-c))
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_triangle_centroid(p1, p2, p3, c)
    implicit none
integer, parameter :: DP=kind(1.d0)
    real(DP)                     ::  p1(3), p2(3), p3(3), c(3)

    c(1:3) = (p1(1:3) + p2(1:3) + p3(1:3)) / 3.d0
end subroutine calc_triangle_centroid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate strain from displacement gradient
!
! Parameters:
!   dg              [in] displacement gradient
!   e               [out] strain tensors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_strain(dg, e)
    implicit none
integer, parameter :: DP=kind(1.d0)
    
    real(DP)             ::  dg(9), e(9)
    
    e(1) = dg(1)
    e(2) = ( dg(3*1+1) + dg(2) ) / 2.0
    e(3) = ( dg(3) + dg(3*2 + 1) ) / 2.0
    
    e(4) = ( dg(3*1+1) + dg(2) ) / 2.0
    e(5) = dg(5)
    e(6) = ( dg(3*1+3) + dg(3*2+2) ) / 2.0
    
    e(7) = ( dg(3) + dg(3*2+1) ) / 2.0
    e(DP) = ( dg(3*1+3) + dg(3*2+2) ) / 2.0
    e(9) = dg(9)
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate stress tensor from strain
!
! Parameters:
!   e               [in] strain tensor
!   sig             [out] stress tensor
!   l               [in] 
!   miu             [in] 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_stress(e, sig, l, miu)
    implicit none
  
  integer, parameter :: DP=kind(1.d0)  
    real(DP)             ::  e(9), sig(9), l, miu
    
    sig(1) = l * ( e(1) + e(5) + e(9) ) + 2.0 * miu * e(1)
    sig(2) = 2.0 * miu * e(2)
    sig(3) = 2.0 * miu * e(3)
    
    sig(4) = 2.0 * miu * e(4)
    sig(5) = l * ( e(1) + e(5) + e(9) ) + 2.0 * miu * e(5)
    sig(6) = 2.0 * miu * e(6)
    
    sig(7) = 2.0 * miu * e(7)
    sig(DP) = 2.0 * miu * e(DP)
    sig(9) = l * ( e(1) + e(5) + e(9) ) + 2.0 * miu * e(9)
end subroutine

subroutine calc_n_stress(e, sig, l_miu)
    implicit none
   integer, parameter :: DP=kind(1.d0)
    
    real(DP)             ::  e(9), sig(9), l_miu
 
    sig(1) = l_miu * ( e(1) + e(5) + e(9) ) + 2.d0 * e(1)
    sig(2) = 2.d0 * e(2)
    sig(3) = 2.d0 * e(3)
    
    sig(4) = 2.d0 * e(4)
    sig(5) = ( e(1) + e(5) + e(9) ) + 2.d0 * e(5)
    sig(6) = 2.d0 * e(6)
    
    sig(7) = 2.d0 * e(7)
    sig(DP) = 2.d0 * e(DP)
    sig(9) = (e(1) + e(5) + e(9)) + 2.d0 * e(9)
end subroutine calc_n_stress


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate transformered stress tensor
!
! Parameters:
!   sig             [in] input stress tensor
!   v               [in] new axis & old axis's cos value array
!   s               [out] output new stress tensor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine stress_transformation(sig, v, s)
  implicit none
  integer, parameter :: DP=kind(1.d0)

  real(DP)             ::  sig(9), &       ! input stress tensor
                          v(9), &         ! new axis & old axis ' cos value
                          s(9)            ! output new stress tensor
                            
  integer             ::  i, j, id, jd
    
  do i=1, 9
    s(i) = 0.0D0
  end do
    
  do id=1, 3
    do jd=1, 3          
      do i=1, 3
        do j=1, 3
          s((id-1)*3+jd) = s((id-1)*3 + jd) + &
                  v((id-1)*3+i) * v((jd-1)*3+j) * sig((i-1)*3+j)
        end do
      end do
    end do
  end do
end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate cos value between two coordinate's axes
!
! Parameters:
!   c1              [in] old coordinate
!   c2              [in] new coordinate
!   v               [out] cos value of coodinate axes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_coord_cos(c1, c2, v)
  implicit none
  integer, parameter :: DP=kind(1.d0)

  real(DP)             ::  c1(3, 3), &     ! old coordinate
                          c2(3, 3)        ! new coordinate
  real(DP)             ::  v(9)
  
  integer             ::  i, j
  real(DP)             ::  f1, f2, av1, av2
  
  do j=1, 3
    do i=1, 3
      f1  = c1(i, 1) * c2(j, 1) + c1(i, 2) * c2(j, 2) + c1(i, 3) * c2(j, 3)

      av1 = sqrt(c1(i, 1)**2 + c1(i, 2)**2 + c1(i, 3)**2)
      av2 = sqrt(c2(j, 1)**2 + c2(j, 2)**2 + c2(j, 3)**2)
      f2  = av1 * av2
            
      v((j-1)*3 + i) = f1 / f2
    end do
  end do
end subroutine calc_coord_cos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate vector's product
!
! Parameters:
!   p1              [in] vector 1
!   p2              [in] vector 2
!   p               [in] normal vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vector_product(p1, p2, p)
  implicit none
  integer, parameter :: DP=kind(1.d0) 

  real(DP)                     ::  p1(3), p2(3), p(3)
    
  p(1) = p1(2)*p2(3) - p1(3)*p2(2)
  p(2) = p1(3)*p2(1) - p1(1)*p2(3)
  p(3) = p1(1)*p2(2) - p1(2)*p2(1)
end subroutine vector_product

subroutine unit_vect(v)
  implicit none
 integer, parameter :: DP=kind(1.d0)
  
  real(DP)                   :: v(3)
  real(DP)                   :: l

  l = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
  v(1) = v(1) / l
  v(2) = v(2) / l
  v(3) = v(3) / l
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate given trianle's normal vector
!
! Parameters:
!   p1              [in] triangle vertex 1
!   p2              [in] triangle vertex 2
!   p3              [in] triangle vertex 3
!   p               [out] normal vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tri_normal_vect(p1, p2, p3, p)
  implicit none
 integer, parameter :: DP=kind(1.d0)
    
  real(DP)                     ::  p1(3), p2(3), p3(3), p(3)
  real(DP)                     ::  v1(3), v2(3)
    
  v1(1:3) = p1(1:3) - p2(1:3)
  v2(1:3) = p2(1:3) - p3(1:3)
    
  call vector_product(v1, v2, p)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate special parameters (ss, ds, op)
!
! Parameters:
!   v1              [in] triangle vertex 1
!   v2              [in] triangle vertex 2
!   v3              [in] triangle vertex 3
!   v_pl            [in] input direction vector
!   ss              [out] ss scale parameter (for hor vector)
!   ds              [out] ds scale parameter (for downdip vector)
!   op              [out] op scale parameter (for perp vector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_ss_ds(v1, v2, v3, v_pl, ss, ds, op)
 use m_calc_green 
 implicit none
    
  real(DP)                     ::  v1(3), v2(3), v3(3)
  real(DP)                     ::  v_pl(3)
  real(DP)                     ::  ss, ds, op
  
  real(DP)                     ::  vpl, vpl0, theta
  real(DP)                     ::  nv(3)
  real(DP)                     ::  x, y, z
    
  real(DP)                     ::  a, b
  real(DP)                     ::  rl,a11,a12,a13,gamma
  integer                     ::  i
    
  ! get triangle's normal vector
  call tri_normal_vect(v1, v2, v3, nv)

  if ( nv(3) .lt. 0 ) then
    do i=1, 3
      nv(i) = -nv(i)
    end do
  end if

  ! Handle vertical triangle (nv(3) = 0)
  if (abs(nv(3)) .lt. 1.0d-12) then
    ! For vertical triangles, use predefined coordinate system
    ! Set strike direction along x-axis, dip direction along y-axis
    ! This follows the convention for vertical fault planes
    ss = vpl1  ! Strike-slip component
    ds = vpl2  ! Dip-slip component  
    op = 0.0d0 ! Opening component (no opening for vertical faults)
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
  ss =(nv(2)*a11-nv(1)*a12)/sqrt(nv(1)**2+nv(2)**2)
  ds =sqrt(nv(1)**2+nv(2)**2-(nv(2)*a11-nv(1)*a12)**2)
  ds =ds/sqrt(nv(1)**2+nv(2)**2)
  op =0.d0
    
end subroutine calc_ss_ds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate local coordinate
!
! Parameters:
!   v1              [in] triangle vertex 1
!   v2              [in] triangle vertex 2
!   v3              [in] triangle vertex 3
!   v_pl            [in] input direction vector
!   c               [out] local coordinate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_local_coordinate2(v1, v2, v3, v_pl, c)
 use m_calc_green  
  implicit none
    
    real(DP)                     ::  v1(3), v2(3), v3(3)
    real(DP)                     ::  v_pl(3)
    real(DP)                     ::  c(3, 3)
    
    real(DP)                     ::  nv(3)
    real(DP)                     ::  ss, ds, op
    real(DP)                     ::  a1(3), a2(3), a3(3)
    real(DP)                     ::  rl,gamma
    integer                     ::  i
    
    ! get ss, ds, op
    call calc_ss_ds(v1, v2, v3, v_pl, ss, ds, op)
    
    ! get triangle's normal vector
    call tri_normal_vect(v1, v2, v3, nv)
    
    if( nv(3) .lt. 0 ) then
       do i=1, 3
          nv(i) = -nv(i)
       end do
    end if

! Check for vertical triangle (nv(3) = 0)
    if (abs(nv(3)) .lt. 1.0d-12) then
      ! a1 = strike direction (along x-axis)
      a1(1) = -nv(2)/sqrt(nv(1)**2+nv(2)**2)
      a1(2) = nv(1)/sqrt(nv(1)**2+nv(2)**2)
      a1(3) = 0.0
     
      a3(1:3) = nv(1:3)
      ! calc axis 2
    call vector_product(a1, a3, a2)


    call unit_vect(a1)
    call unit_vect(a2)
    call unit_vect(a3)
      
      ! Set output coordinate system
      c(1, 1:3) = a1(1:3)  ! Strike axis
      c(2, 1:3) = a2(1:3)  ! Dip axis  
      c(3, 1:3) = a3(1:3)  ! Normal axis
      return
    end if

    ! Handle horizontal triangles (nv(1) and nv(2) both zero)
    if (abs(nv(1)) < 1.0d-12 .and. abs(nv(2)) < 1.0d-12) then
      write(*,*) "INFO: Horizontal triangle detected in calc_local_coordinate2 - using predefined axes"
      write(*,*) "Triangle vertices: v1=", v1, " v2=", v2, " v3=", v3
      write(*,*) "Normal vector: nv=", nv
      
      ! For horizontal triangles, construct coordinate system by convention
      ! a1 = strike direction (along x-axis)
      a1(1) = 1.0d0
      a1(2) = 0.0d0
      a1(3) = 0.0d0
      
      ! a3 = normal direction (upward z-axis)
      a3(1) = 0.0d0
      a3(2) = 0.0d0
      a3(3) = 1.0d0
      
      ! a2 = dip direction (along y-axis, perpendicular to a1 and a3)
      a2(1) = 0.0d0
      a2(2) = 1.0d0
      a2(3) = 0.0d0
      
      ! Set output coordinate system
      c(1, 1:3) = a1(1:3)  ! Strike axis
      c(2, 1:3) = a2(1:3)  ! Dip axis  
      c(3, 1:3) = a3(1:3)  ! Normal axis
      return
    end if

    !write(6,*) "nv(3)check",nv(3)
    ! Check for valid ss and ds values
    if (isnan(ss) .or. isnan(ds)) then
      write(*,*) "ERROR: Invalid ss or ds in calc_local_coordinate2"
      write(*,*) "ss =", ss, "ds =", ds
      c = 0.0d0
      return
    end if
    !write(6,*) "nv(3)check",nv(3)
    ! calc axis 1
    a1(1) = ss
    a1(2) = ds
    
    a3(1:3) = nv(1:3)
    a1(3) = -( a3(1)*a1(1) + a3(2)*a1(2) ) / a3(3)
        
    rl=(vpl1**2+vpl2**2)+(nv(1)*vpl1+nv(2)*vpl2)**2/nv(3)**2
    rl=1.d0/rl
    rl=sqrt(rl)
    gamma=-(nv(1)*vpl1+nv(2)*vpl2)*rl/nv(3)

    a1(1) = vpl1*rl
    a1(2) = vpl2*rl

  a1(3) = gamma

    ! calc axis 2
    call vector_product(a1, a3, a2)
    if ( a3(3) .lt. 0 ) then
      a3(1:3) = a3(1:3)
    end if


    call unit_vect(a1)
    call unit_vect(a2)
    call unit_vect(a3)
    
    c(1, 1:3) = a1(1:3)
    c(2, 1:3) = a2(1:3)
    c(3, 1:3) = a3(1:3)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_green_allcell_improved(myid,size,Nt,arr_vertex,arr_cell, & 
                n_vertex,n_cell,cells_processed,base_cells,extra_cells, error_occurred, error_message)
  use m_calc_green,only: DP,parm_nu,parm_l,parm_miu,vpl1,vpl2,PI,ZERO 
  use mod_dtrigreen
  implicit none
  
    integer, intent(in) :: cells_processed,base_cells,extra_cells
   integer, intent(in) :: myid, size, Nt, n_cell, n_vertex
   real(DP), intent(in) :: arr_vertex(n_vertex,3)
   integer, intent(in) :: arr_cell(n_cell,3)
   logical, intent(inout) :: error_occurred
   character(len=*), intent(inout) :: error_message
    
  real(DP) ::       u(3), t(9)
  real(DP) ::                ss, ds, op
  
  integer             :: i, j, k
  integer             ::  vj(3)
  real(DP)             ::  p1(3), p2(3), p3(3), co(3)
  real(DP)             ::  sig33(3,3)
  real(DP)             ::  vpl(3)
  character(20) ::        cTemp
  real(DP)             ::  l_miu
  
  real(DP)             ::  c_local2(3, 3), c_global(3, 3), &
                           c_local(3, 3),  c_local_v(9), c_local_v2(9)
  real(DP),allocatable :: arr_co(:,:),arr_trid(:,:),arr_cl_v2(:,:,:)
  real(DP),allocatable :: arr_out(:,:)
  
  ! Dynamic load balancing variables
  integer :: local_cells, start_idx
  integer :: ierr
  
 ! global coordinate
  c_global(:,:) = 0.d0
  do i=1,3
    c_global(i,i) =1.d0
  end do


  ! Get distribution parameters from main program
  local_cells = cells_processed
  
  ! Calculate start_idx based on myid and distribution
  if (myid < extra_cells) then
     start_idx = myid * (base_cells + 1)
  else
     start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells
  end if
   
   ! Initialize variables
   vpl(1:3) = 1.d0
   l_miu = parm_l/parm_miu
   ss = -1.d0
   ds = 0.d0
   op = 0.d0
   
   ! Allocate arrays only for local cells (memory efficient)
   ! Note: arr_trid must be (9,n_cell) because each process needs triangle data for ALL cells
   ! to compute Green's functions between local cells and all other cells
   ! Ensure minimum size of 1 to avoid allocation issues
   allocate(arr_co(3,max(1,local_cells)), arr_trid(9,n_cell), arr_cl_v2(3,3,max(1,local_cells)))
   allocate(arr_out(max(1,local_cells), n_cell))
  
          ! Pre-compute cell properties for local cells only
       if (local_cells > 0) then
         do j = 1, local_cells
           ! Map local index to global cell index
           k = start_idx + j - 1
           if (k >= n_cell) exit
           
           ! Bounds checking for safety
           if (k < 1 .or. k > n_cell) then
             write(*,*) 'Error: Invalid cell index k =', k, 'for j =', j
             cycle
           end if
           
           vj(1:3) = arr_cell(k,1:3)
           
           ! Bounds checking for vertex indices
           if (vj(1) < 1 .or. vj(1) > n_vertex .or. &
               vj(2) < 1 .or. vj(2) > n_vertex .or. &
               vj(3) < 1 .or. vj(3) > n_vertex) then
             write(*,*) 'Error: Invalid vertex indices for cell', k, ':', vj(1), vj(2), vj(3)
             cycle
           end if
           
           p1(1:3) = arr_vertex(vj(1),1:3)
           p2(1:3) = arr_vertex(vj(2),1:3)
           p3(1:3) = arr_vertex(vj(3),1:3)
           
           call calc_triangle_centroid(p1, p2, p3, co)
           arr_co(1:3,j) = co(1:3)
           arr_trid(1:3,k) = p1(1:3)
           arr_trid(4:6,k) = p2(1:3)
           arr_trid(7:9,k) = p3(1:3)
           
           call calc_local_coordinate2(p1, p2, p3, vpl, c_local2)

           call calc_coord_cos(c_global, c_local2, c_local_v2)
           arr_cl_v2(1:3,1,j) = c_local_v2(1:3)
           arr_cl_v2(1:3,2,j) = c_local_v2(4:6)
           arr_cl_v2(1:3,3,j) = c_local_v2(7:9)
           
           arr_out(j,:) = 0.d0
         end do
       else
         ! Initialize arrays for processes with no cells
         arr_out(1,:) = 0.d0
         arr_co(:,1) = 0.d0
         arr_cl_v2(:,:,1) = 0.d0
       end if
   
   ! CRITICAL: We need triangle data for ALL cells, not just local ones
   ! Each process must compute triangle data for all cells to perform Green's function calculations
   ! This is a necessary overhead for the distributed computation
   write(*,*) "Process", myid, "computing triangle data for all", n_cell, "cells"
   do k = 1, n_cell
     vj(1:3) = arr_cell(k,1:3)
     p1(1:3) = arr_vertex(vj(1),1:3)
     p2(1:3) = arr_vertex(vj(2),1:3)
     p3(1:3) = arr_vertex(vj(3),1:3)
     
     arr_trid(1:3,k) = p1(1:3)
     arr_trid(4:6,k) = p2(1:3)
     arr_trid(7:9,k) = p3(1:3)
   end do
   write(*,*) "Process", myid, "completed triangle data computation"

  write(6,*) "Process", myid, "starting hybrid parallel computation"
  
  ! Hybrid MPI+OpenMP parallel computation
  if (local_cells > 0) then
    !$OMP PARALLEL DO PRIVATE(i, j, u, t, sig33, k) SHARED(arr_co, arr_trid, arr_out, arr_cl_v2)
    do j = 1, local_cells
      ! Map local index to global cell index
      k = start_idx + j - 1
      if (k >= n_cell) cycle
      
      do i = 1, n_cell
        ! Calculate strain gradients using Stuart's method
        call dstuart(parm_nu, arr_co(:,j), arr_trid(:,i), ss, ds, op, u, t)
        
        ! Check for invalid results from dstuart
        if (any(isnan(u)) .or. any(isnan(t))) then
          !$OMP CRITICAL
          error_occurred = .true.
          error_message = "Invalid results from dstuart calculation"
          write(*,*) arr_co(:,j),arr_trid(:,i),ss,ds,op,u,t
          !$OMP END CRITICAL
          cycle
        end if
              
        ! Calculate stress tensor components
        sig33(3,3) = l_miu * (t(1) + t(5) + t(9)) 
        sig33(1,1) = sig33(3,3) + 2.d0 * t(1)
        sig33(2,2) = sig33(3,3) + 2.d0 * t(5)
        sig33(3,3) = sig33(3,3) + 2.d0 * t(9)

        sig33(2,1) = t(4) + t(2)
        sig33(3,1) = t(3) + t(7)
        sig33(3,2) = t(6) + t(DP)

        sig33(1,2) = sig33(2,1)
        sig33(1,3) = sig33(3,1)
        sig33(2,3) = sig33(3,2)

        ! Check for invalid stress tensor
        if (any(isnan(sig33))) then
          !$OMP CRITICAL
          error_occurred = .true.
          error_message = "Invalid stress tensor calculated"
          !$OMP END CRITICAL
          cycle
        end if

        ! Calculate local stress in Bar (0.1MPa)
        arr_out(j,i) = -parm_miu/100 * dot_product(arr_cl_v2(:,3,j), matmul(sig33(:,:), arr_cl_v2(:,1,j)))
        
        ! Check final result
        if (isnan(arr_out(j,i))) then
          !$OMP CRITICAL
          error_occurred = .true.
          error_message = "Invalid output value calculated"
          write(*,*) arr_cl_v2(:,3,j),sig33(1,1),sig33(2,2),sig33(3,3)
          !$OMP END CRITICAL
          cycle
        end if
      end do
    end do
    !$OMP END PARALLEL DO
  else
    ! Process with no cells - no computation needed
    write(*,*) 'Process', myid, 'skipping computation (no cells assigned)'
  end if
  
    ! Check for errors before proceeding
  if (error_occurred) then
    write(*,*) "Process", myid, "encountered an error: ", trim(error_message)
    return
  end if
  
  ! Count processed cells after OpenMP loop
 ! cells_processed = local_cells
   
   write(cTemp,*) myid
   write(*,*) "Process", myid, "completed", local_cells, "cells"
   write(*,*) "Process", myid, "processed cells from index", start_idx, "to", start_idx + local_cells - 1

   ! Performance monitoring for load balancing
   if (myid == 0) then
      write(*,*) "=========================================="
      write(*,*) "Dynamic Load Balancing Summary:"
      write(*,*) "Total cells:", n_cell
      write(*,*) "MPI processes:", size
      write(*,*) "Base cells per process:", base_cells
      write(*,*) "Extra cells distributed:", extra_cells
      write(*,*) "Load imbalance:", extra_cells, "cells (max difference between processes)"
      write(*,*) "=========================================="
   end if

  ! Write results to file
  open(14, file='trigreen_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
  
  ! Only master process writes position data
  if(myid == 0) then 
    open(22, file='position.bin', form='unformatted', access='stream')
    do j = 1, n_cell
      vj(1:3) = arr_cell(j,1:3)
     p1(1:3) = arr_vertex(vj(1),1:3)
     p2(1:3) = arr_vertex(vj(2),1:3)
     p3(1:3) = arr_vertex(vj(3),1:3)

     call calc_triangle_centroid(p1, p2, p3, co)
      write(22) co(1), co(2), co(3)
    end do
    close(22)
  endif

  ! Write local results (only if we have cells to write)
  if (local_cells > 0) then
    do i = 1, local_cells
      write(14) arr_out(i,:)
    end do
  else
    ! Write a dummy entry to ensure file is not empty
    ! This ensures compatibility with 3dtri_BP5.f90 which expects all processes to have files
    write(14) 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
    write(*,*) 'Process', myid, 'wrote dummy entry (no cells assigned)'
  end if
  close(14)

   ! Deallocate arrays allocated in this subroutine
   deallocate (arr_co,arr_trid,arr_cl_v2)
   deallocate (arr_out)
 
return 
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Performance monitoring subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine performance_monitoring(myid, start_time, end_time, cells_processed)
  use mpi
  implicit none
  
 integer, parameter :: DP=kind(1.d0)
  integer, intent(in) :: myid, cells_processed
  real(DP), intent(in) :: start_time, end_time
  
  real(DP) :: local_time, total_time, avg_time
  integer :: size, ierr
  
  ! Get the size of MPI communicator
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
  
  local_time = end_time - start_time
  
  ! Gather timing information from all processes
  call MPI_Reduce(local_time, total_time, 1, MPI_DOUBLE_PRECISION,& 
          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  if (myid == 0) then
    avg_time = total_time / size
    write(*,*) "Performance Summary:"
    write(*,*) "  Total cells processed:", cells_processed
    write(*,*) "  Average time per process:", avg_time
    write(*,*) "  Cells per second:", cells_processed / avg_time
  endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! isnan function for checking invalid values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical function isnan(x)
  use m_calc_green,only: DP
  implicit none
  real(DP), intent(in) :: x
  isnan = (x /= x)
end function isnan








