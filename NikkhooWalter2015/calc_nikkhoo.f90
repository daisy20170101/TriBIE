!===============================================================================
! calc_nikkhoo.f90
!
! MPI/OpenMP parallel program for calculating Green's function coefficients
! (stiffness matrix) using the Nikkhoo & Walter (2015) triangular dislocation
! method for elastic half-space.
!
! This program calculates the shear stress change at each triangular element
! centroid due to unit slip on every other element, producing a stiffness
! matrix for boundary element earthquake simulations.
!
! Reference: Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An
! analytical, artefact-free solution. Geophysical Journal International
!
! Usage:
!   mpirun -np <nprocs> ./calc_nikkhoo
!
! Input:
!   triangular_mesh.gts - GTS format mesh file with triangular elements
!
! Output:
!   trigreen_<rank>.bin - Binary stiffness matrix files (one per MPI process)
!   position.bin - Centroid positions of all elements
!===============================================================================

program calc_nikkhoo
  use m_nikkhoo_green
  use mpi
  use omp_lib
  implicit none

  ! File name for mesh input
  character(*), parameter :: fname = "triangular_mesh.gts"

  ! Mesh data arrays
  integer :: n_vertex, n_edge, n_cell
  real(DP), dimension(:,:), allocatable :: arr_vertex
  integer, dimension(:,:), allocatable :: arr_edge
  integer, dimension(:,:), allocatable :: arr_cell

  ! MPI variables
  integer :: ierr, size, myid, master
  real(DP) :: start_time, end_time

  ! Load balancing variables
  integer :: Nt, Nt_all, local_cells, cells_processed
  integer :: base_cells, extra_cells, start_idx

  ! OpenMP setup
  integer :: num_threads

  ! Error handling
  logical :: error_occurred
  character(len=256) :: error_message

  ! Initialize error flag
  error_occurred = .false.
  error_message = ""

  !============================================================================
  ! Initialize MPI
  !============================================================================
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) then
    write(*,*) "Error: Failed to initialize MPI"
    stop
  end if

  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
  master = 0

  ! Set OpenMP threads per MPI process
  num_threads = omp_get_max_threads()
  call omp_set_num_threads(num_threads)

  if (myid == master) then
    start_time = MPI_Wtime()
    write(*,*) "========================================================"
    write(*,*) "Nikkhoo-Walter Triangular Stiffness Calculation"
    write(*,*) "========================================================"
    write(*,*) "MPI processes:", size
    write(*,*) "OpenMP threads per process:", num_threads
    write(*,*) "Total parallel threads:", size * num_threads
    write(*,*) ""
  end if

  !============================================================================
  ! Load mesh data
  !============================================================================
  if (.not. error_occurred) then
    call load_name(fname, n_vertex, n_edge, n_cell, error_occurred, error_message)
  end if

  if (.not. error_occurred) then
    allocate(arr_vertex(n_vertex, 3), arr_edge(n_edge, 2), arr_cell(n_cell, 3), stat=ierr)
    if (ierr /= 0) then
      error_occurred = .true.
      error_message = "Failed to allocate mesh arrays"
    end if
  end if

  if (.not. error_occurred) then
    call load_gts(fname, n_vertex, n_edge, n_cell, arr_vertex, arr_edge, arr_cell, &
                  error_occurred, error_message)
  end if

  !============================================================================
  ! Dynamic load balancing
  !============================================================================
  if (.not. error_occurred) then
    if (myid == master) then
      write(*,*) "Mesh loaded successfully:"
      write(*,*) "  Vertices:", n_vertex
      write(*,*) "  Cells:", n_cell
      write(*,*) ""
    end if

    Nt_all = n_cell

    if (size == 1) then
      ! Single CPU execution
      if (myid == master) write(*,*) "Single CPU execution - all cells assigned to process 0"
      base_cells = n_cell
      extra_cells = 0
      local_cells = n_cell
      start_idx = 1
      cells_processed = n_cell
      Nt = n_cell

      call calc_nikkhoo_allcell(myid, size, Nt, arr_vertex, arr_cell, &
                                n_vertex, n_cell, cells_processed, base_cells, extra_cells, &
                                error_occurred, error_message)
    else
      ! Multi-CPU dynamic load balancing
      base_cells = Nt_all / size
      extra_cells = mod(Nt_all, size)

      if (myid == master) then
        if (mod(n_cell, size) /= 0) then
          write(*,*) "Dynamic load balancing enabled:"
          write(*,*) "  Base cells per process:", base_cells
          write(*,*) "  Extra cells to distribute:", extra_cells
        else
          write(*,*) "Even distribution: each process gets", base_cells, "cells"
        end if
        write(*,*) ""
      end if

      ! Calculate local cell count for this process
      if (myid < extra_cells) then
        local_cells = base_cells + 1
        start_idx = myid * (base_cells + 1) + 1
      else
        local_cells = base_cells
        start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells + 1
      end if

      cells_processed = local_cells
      Nt = local_cells

      call calc_nikkhoo_allcell(myid, size, Nt, arr_vertex, arr_cell, &
                                n_vertex, n_cell, cells_processed, base_cells, extra_cells, &
                                error_occurred, error_message)
    end if
  end if

  !============================================================================
  ! Cleanup and finalize
  !============================================================================
  if (allocated(arr_vertex)) deallocate(arr_vertex)
  if (allocated(arr_edge)) deallocate(arr_edge)
  if (allocated(arr_cell)) deallocate(arr_cell)

  if (error_occurred) then
    write(*,*) "Process", myid, "error:", trim(error_message)
  else
    call MPI_barrier(MPI_COMM_WORLD, ierr)

    if (myid == master) then
      end_time = MPI_Wtime()
      write(*,*) ""
      write(*,*) "========================================================"
      write(*,*) "Calculation completed successfully!"
      write(*,*) "========================================================"
      write(*,*) "Total time:", end_time - start_time, "seconds"
      write(*,*) "Output files: trigreen_<rank>.bin, position.bin"
    end if
  end if

  call MPI_finalize(ierr)

  if (error_occurred) stop 1

end program calc_nikkhoo


!===============================================================================
! Load mesh dimension from GTS file
!===============================================================================
subroutine load_name(fname, n_vertex, n_edge, n_cell, error_occurred, error_message)
  implicit none
  integer, parameter :: DP = kind(1.d0)

  character(*), intent(in) :: fname
  integer, intent(out) :: n_vertex, n_edge, n_cell
  logical, intent(inout) :: error_occurred
  character(len=*), intent(inout) :: error_message

  integer :: iostat

  open(unit=33, file=fname, status='old', iostat=iostat)
  if (iostat /= 0) then
    error_occurred = .true.
    error_message = "Failed to open file: " // trim(fname)
    return
  end if

  read(33, *, iostat=iostat) n_vertex, n_edge, n_cell
  if (iostat /= 0) then
    error_occurred = .true.
    error_message = "Failed to read mesh dimensions"
    close(33)
    return
  end if

  close(33)

  if (n_vertex <= 0 .or. n_cell <= 0) then
    error_occurred = .true.
    error_message = "Invalid mesh dimensions"
  end if

end subroutine load_name


!===============================================================================
! Load mesh data from GTS file
!===============================================================================
subroutine load_gts(fname, n_vertex, n_edge, n_cell, arr_vertex, arr_edge, arr_cell, &
                    error_occurred, error_message)
  implicit none
  integer, parameter :: DP = kind(1.d0)

  character(*), intent(in) :: fname
  integer, intent(in) :: n_vertex, n_edge, n_cell
  real(DP), intent(out) :: arr_vertex(n_vertex, 3)
  integer, intent(out) :: arr_edge(n_edge, 2)
  integer, intent(out) :: arr_cell(n_cell, 3)
  logical, intent(inout) :: error_occurred
  character(len=*), intent(inout) :: error_message

  integer :: i, iostat
  integer :: dummy_v, dummy_e, dummy_c

  open(unit=10, file=fname, status='old', iostat=iostat)
  if (iostat /= 0) then
    error_occurred = .true.
    error_message = "Failed to open GTS file: " // trim(fname)
    return
  end if

  read(10, *, iostat=iostat) dummy_v, dummy_e, dummy_c
  if (iostat /= 0) then
    error_occurred = .true.
    error_message = "Failed to read GTS header"
    close(10)
    return
  end if

  ! Read vertex coordinates
  do i = 1, n_vertex
    read(10, *, iostat=iostat) arr_vertex(i, 1), arr_vertex(i, 2), arr_vertex(i, 3)
    if (iostat /= 0) then
      error_occurred = .true.
      error_message = "Failed to read vertex data"
      close(10)
      return
    end if
    ! Flip z coordinate sign (mesh has positive z, half-space requires negative z)
    arr_vertex(i, 3) = -arr_vertex(i, 3)
  end do

  ! Read cell definitions (vertex indices)
  do i = 1, n_cell
    read(10, *, iostat=iostat) arr_cell(i, 1), arr_cell(i, 2), arr_cell(i, 3)
    if (iostat /= 0) then
      error_occurred = .true.
      error_message = "Failed to read cell data"
      close(10)
      return
    end if
  end do

  close(10)

end subroutine load_gts


!===============================================================================
! Main calculation subroutine - compute stiffness for all cells
!===============================================================================
subroutine calc_nikkhoo_allcell(myid, size, Nt, arr_vertex, arr_cell, &
                                n_vertex, n_cell, cells_processed, base_cells, extra_cells, &
                                error_occurred, error_message)
  use m_nikkhoo_green
  use nikkhoo_walter, only: tdstress_hs
  use omp_lib
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, intent(in) :: myid, size, Nt, n_vertex, n_cell
  integer, intent(in) :: cells_processed, base_cells, extra_cells
  real(DP), intent(in) :: arr_vertex(n_vertex, 3)
  integer, intent(in) :: arr_cell(n_cell, 3)
  logical, intent(inout) :: error_occurred
  character(len=*), intent(inout) :: error_message

  ! Local variables
  integer :: i, j, k
  integer :: vj(3)
  real(DP) :: p1(3), p2(3), p3(3), co(3)
  real(DP) :: src_p1(3), src_p2(3), src_p3(3)
  real(DP) :: stress(6), strain(6)
  real(DP) :: sig33(3, 3)
  real(DP) :: vpl(3)
  character(20) :: cTemp
  real(DP) :: l_miu

  ! Coordinate transformation
  real(DP) :: c_local(3, 3), c_global(3, 3), c_local_v(9)

  ! Work arrays
  real(DP), allocatable :: arr_co(:,:), arr_trid(:,:), arr_cl_v(:,:,:)
  real(DP), allocatable :: arr_out(:,:)

  ! Load balancing
  integer :: local_cells, start_idx

  ! Slip components for unit slip
  real(DP) :: ss, ds, ts

  ! Global coordinate system
  c_global = 0.d0
  do i = 1, 3
    c_global(i, i) = 1.d0
  end do

  local_cells = cells_processed

  ! Calculate start index for this process
  if (myid < extra_cells) then
    start_idx = myid * (base_cells + 1) + 1
  else
    start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells + 1
  end if

  ! Initialize
  vpl = 1.d0
  l_miu = parm_l / parm_miu

  ! Unit strike-slip (ss=1, ds=0, ts=0)
  ss = 1.d0
  ds = 0.d0
  ts = 0.d0

  ! Allocate work arrays
  allocate(arr_co(3, max(1, local_cells)))
  allocate(arr_trid(9, n_cell))
  allocate(arr_cl_v(3, 3, max(1, local_cells)))
  allocate(arr_out(max(1, local_cells), n_cell))

  ! Pre-compute triangle data for all cells (needed for source triangles)
  do k = 1, n_cell
    vj(1:3) = arr_cell(k, 1:3)
    p1(1:3) = arr_vertex(vj(1), 1:3)
    p2(1:3) = arr_vertex(vj(2), 1:3)
    p3(1:3) = arr_vertex(vj(3), 1:3)

    arr_trid(1:3, k) = p1(1:3)
    arr_trid(4:6, k) = p2(1:3)
    arr_trid(7:9, k) = p3(1:3)
  end do

  ! Pre-compute centroids and local coordinate systems for local cells
  if (local_cells > 0) then
    do j = 1, local_cells
      k = start_idx + j - 1
      if (k > n_cell) cycle

      vj(1:3) = arr_cell(k, 1:3)
      p1(1:3) = arr_vertex(vj(1), 1:3)
      p2(1:3) = arr_vertex(vj(2), 1:3)
      p3(1:3) = arr_vertex(vj(3), 1:3)

      ! Calculate centroid
      co = (p1 + p2 + p3) / 3.d0
      arr_co(1:3, j) = co(1:3)

      ! Calculate local coordinate system (strike, dip, normal)
      call calc_local_coordinate(p1, p2, p3, vpl, c_local)

      ! Calculate coordinate transformation matrix
      call calc_coord_cos(c_global, c_local, c_local_v)
      arr_cl_v(1:3, 1, j) = c_local_v(1:3)
      arr_cl_v(1:3, 2, j) = c_local_v(4:6)
      arr_cl_v(1:3, 3, j) = c_local_v(7:9)

      arr_out(j, :) = 0.d0
    end do
  else
    arr_out(1, :) = 0.d0
    arr_co(:, 1) = 0.d0
    arr_cl_v(:, :, 1) = 0.d0
  end if

  write(*,*) "Process", myid, "starting stiffness calculation for", local_cells, "cells"

  !============================================================================
  ! Main computation loop - MPI + OpenMP hybrid parallelization
  !============================================================================
  if (local_cells > 0) then
    !$OMP PARALLEL DO PRIVATE(i, j, k, src_p1, src_p2, src_p3, stress, strain, sig33) &
    !$OMP& SHARED(arr_co, arr_trid, arr_out, arr_cl_v, n_cell, local_cells, start_idx)
    do j = 1, local_cells
      k = start_idx + j - 1
      if (k > n_cell) cycle

      do i = 1, n_cell
        ! Get source triangle vertices
        src_p1(1:3) = arr_trid(1:3, i)
        src_p2(1:3) = arr_trid(4:6, i)
        src_p3(1:3) = arr_trid(7:9, i)

        ! Calculate stress/strain at observation point due to unit slip on source triangle
        ! Using half-space solution
        call tdstress_hs(arr_co(1, j), arr_co(2, j), arr_co(3, j), &
                         src_p1, src_p2, src_p3, ss, ds, ts, &
                         parm_miu, parm_l, stress, strain)

        ! Check for NaN (singular points)
        if (ieee_is_nan(stress(1)) .or. ieee_is_nan(strain(1))) then
          arr_out(j, i) = 0.d0
          cycle
        end if

        ! Build stress tensor matrix
        sig33(1, 1) = stress(1)  ! Sxx
        sig33(2, 2) = stress(2)  ! Syy
        sig33(3, 3) = stress(3)  ! Szz
        sig33(1, 2) = stress(4)  ! Sxy
        sig33(2, 1) = stress(4)  ! Syx = Sxy
        sig33(1, 3) = stress(5)  ! Sxz
        sig33(3, 1) = stress(5)  ! Szx = Sxz
        sig33(2, 3) = stress(6)  ! Syz
        sig33(3, 2) = stress(6)  ! Szy = Syz

        ! Calculate shear stress in local coordinate system
        ! tau = n' * sigma * s (normal dotted with stress dotted with strike)
        ! Output in Bar (0.1 MPa) - negative sign for convention
        arr_out(j, i) = -1.0d0/100.0d0 * &
                        dot_product(arr_cl_v(:, 3, j), matmul(sig33, arr_cl_v(:, 1, j)))

        ! Check for NaN in result
        if (ieee_is_nan(arr_out(j, i))) then
          arr_out(j, i) = 0.d0
        end if
      end do
    end do
    !$OMP END PARALLEL DO
  end if

  write(*,*) "Process", myid, "completed calculation"

  !============================================================================
  ! Write output files
  !============================================================================
  write(cTemp, *) myid

  if (local_cells > 0) then
    open(14, file='trigreen_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    do i = 1, local_cells
      write(14) arr_out(i, :)
    end do
    close(14)
  else
    ! Write dummy entry for compatibility
    open(14, file='trigreen_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    write(14) (0.d0, i=1, n_cell)
    close(14)
  end if

  ! Master writes position data
  if (myid == 0) then
    open(22, file='position.bin', form='unformatted', access='stream')
    do j = 1, n_cell
      vj(1:3) = arr_cell(j, 1:3)
      p1(1:3) = arr_vertex(vj(1), 1:3)
      p2(1:3) = arr_vertex(vj(2), 1:3)
      p3(1:3) = arr_vertex(vj(3), 1:3)

      co = (p1 + p2 + p3) / 3.d0
      write(22) co(1), co(2), co(3)
    end do
    close(22)
  end if

  ! Load balancing summary
  if (myid == 0) then
    write(*,*) ""
    write(*,*) "Load Balancing Summary:"
    write(*,*) "  Total cells:", n_cell
    write(*,*) "  MPI processes:", size
    write(*,*) "  Base cells per process:", base_cells
    write(*,*) "  Extra cells distributed:", extra_cells
  end if

  ! Cleanup
  deallocate(arr_co, arr_trid, arr_cl_v, arr_out)

end subroutine calc_nikkhoo_allcell


!===============================================================================
! Helper: Cross product
!===============================================================================
subroutine cross_product_local(a, b, c)
  use m_nikkhoo_green, only: DP
  implicit none
  real(DP), dimension(3), intent(in) :: a, b
  real(DP), dimension(3), intent(out) :: c

  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)
end subroutine cross_product_local


!===============================================================================
! Helper: Calculate local coordinate system
!===============================================================================
subroutine calc_local_coordinate(v1, v2, v3, v_pl, c)
  use m_nikkhoo_green, only: DP, PI, vpl1, vpl2, EPS
  implicit none

  real(DP), intent(in) :: v1(3), v2(3), v3(3), v_pl(3)
  real(DP), intent(out) :: c(3, 3)

  real(DP) :: nv(3), side1(3), side2(3)
  real(DP) :: a1(3), a2(3), a3(3)
  real(DP) :: denom, rl, gamma
  integer :: i

  ! Calculate normal vector
  side1 = v2 - v1
  side2 = v3 - v1
  nv(1) = side1(2) * side2(3) - side1(3) * side2(2)
  nv(2) = side1(3) * side2(1) - side1(1) * side2(3)
  nv(3) = side1(1) * side2(2) - side1(2) * side2(1)

  denom = sqrt(sum(nv**2))
  if (denom < EPS) then
    c = 0.0_DP
    do i = 1, 3
      c(i, i) = 1.0_DP
    end do
    return
  end if
  nv = nv / denom

  ! Ensure normal points upward (positive z)
  if (nv(3) < 0.0_DP) nv = -nv

  ! Handle vertical triangles
  if (abs(nv(3)) < EPS) then
    denom = sqrt(nv(1)**2 + nv(2)**2)
    if (denom < EPS) then
      a1 = [1.0_DP, 0.0_DP, 0.0_DP]
    else
      a1(1) = -nv(2) / denom
      a1(2) = nv(1) / denom
      a1(3) = 0.0_DP
    end if

    a3 = nv
    call cross_product_local(a1, a3, a2)

    denom = sqrt(sum(a1**2)); if (denom > EPS) a1 = a1 / denom
    denom = sqrt(sum(a2**2)); if (denom > EPS) a2 = a2 / denom
    denom = sqrt(sum(a3**2)); if (denom > EPS) a3 = a3 / denom

    c(1, :) = a1
    c(2, :) = a2
    c(3, :) = a3
    return
  end if

  ! Handle horizontal triangles
  if (abs(nv(1)) < EPS .and. abs(nv(2)) < EPS) then
    a1 = [1.0_DP, 0.0_DP, 0.0_DP]
    a3 = [0.0_DP, 0.0_DP, 1.0_DP]
    a2 = [0.0_DP, 1.0_DP, 0.0_DP]

    c(1, :) = a1
    c(2, :) = a2
    c(3, :) = a3
    return
  end if

  ! General case
  rl = (vpl1**2 + vpl2**2) + (nv(1) * vpl1 + nv(2) * vpl2)**2 / nv(3)**2
  rl = 1.0_DP / rl
  rl = sqrt(rl)
  gamma = -(nv(1) * vpl1 + nv(2) * vpl2) * rl / nv(3)

  a1(1) = vpl1 * rl
  a1(2) = vpl2 * rl
  a1(3) = gamma

  a3 = nv
  call cross_product_local(a1, a3, a2)

  denom = sqrt(sum(a1**2)); if (denom > EPS) a1 = a1 / denom
  denom = sqrt(sum(a2**2)); if (denom > EPS) a2 = a2 / denom
  denom = sqrt(sum(a3**2)); if (denom > EPS) a3 = a3 / denom

  c(1, :) = a1
  c(2, :) = a2
  c(3, :) = a3

end subroutine calc_local_coordinate


!===============================================================================
! Helper: Calculate coordinate transformation cosines
!===============================================================================
subroutine calc_coord_cos(c1, c2, v)
  use m_nikkhoo_green, only: DP, EPS
  implicit none

  real(DP), intent(in) :: c1(3, 3), c2(3, 3)
  real(DP), intent(out) :: v(9)

  integer :: i, j
  real(DP) :: f1, f2, av1, av2

  do j = 1, 3
    do i = 1, 3
      f1 = c1(i, 1) * c2(j, 1) + c1(i, 2) * c2(j, 2) + c1(i, 3) * c2(j, 3)

      av1 = sqrt(c1(i, 1)**2 + c1(i, 2)**2 + c1(i, 3)**2)
      av2 = sqrt(c2(j, 1)**2 + c2(j, 2)**2 + c2(j, 3)**2)
      f2 = av1 * av2

      if (f2 > EPS) then
        v((j - 1) * 3 + i) = f1 / f2
      else
        v((j - 1) * 3 + i) = 0.0_DP
      end if
    end do
  end do
end subroutine calc_coord_cos
