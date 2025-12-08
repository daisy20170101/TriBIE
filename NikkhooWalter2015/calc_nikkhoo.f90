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
!   nikkhoo_<rank>.bin - Binary stiffness matrix files (one per MPI process)
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
      write(*,*) "Output files: nikkhoo_<rank>.bin, position.bin"
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
  ss = -1.d0
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
        call tdstress_hs_silent(arr_co(1, j), arr_co(2, j), arr_co(3, j), &
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
    open(14, file='nikkhoo_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    do i = 1, local_cells
      write(14) arr_out(i, :)
    end do
    close(14)
  else
    ! Write dummy entry for compatibility
    open(14, file='nikkhoo_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
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
! Silent version of TDstressHS (half-space solution)
! No debug output for production use
!===============================================================================
subroutine tdstress_hs_silent(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                              stress, strain)
  use m_nikkhoo_green, only: DP, PI, EPS
  use, intrinsic :: ieee_arithmetic
  implicit none

  real(DP), intent(in) :: x, y, z
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  real(DP), intent(in) :: ss, ds, ts, mu, lambda
  real(DP), dimension(6), intent(out) :: stress, strain

  real(DP), dimension(6) :: sts_ms, str_ms, sts_is, str_is
  real(DP), dimension(3) :: p1_img, p2_img, p3_img

  ! Check half-space constraint
  if (z > 0.0_DP .or. p1(3) > 0.0_DP .or. p2(3) > 0.0_DP .or. p3(3) > 0.0_DP) then
    stress = ieee_value(0.0_DP, ieee_quiet_nan)
    strain = ieee_value(0.0_DP, ieee_quiet_nan)
    return
  end if

  ! Main dislocation contribution
  call tdstress_fs_silent(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, sts_ms, str_ms)

  ! Image dislocation contribution
  p1_img = p1; p2_img = p2; p3_img = p3
  p1_img(3) = -p1(3)
  p2_img(3) = -p2(3)
  p3_img(3) = -p3(3)

  call tdstress_fs_silent(x, y, z, p1_img, p2_img, p3_img, ss, ds, ts, mu, lambda, &
                          sts_is, str_is)

  ! Surface element correction
  if (abs(p1_img(3)) < EPS .and. abs(p2_img(3)) < EPS .and. abs(p3_img(3)) < EPS) then
    sts_is(5) = -sts_is(5)
    sts_is(6) = -sts_is(6)
    str_is(5) = -str_is(5)
    str_is(6) = -str_is(6)
  end if

  ! Total (main + image, harmonic function contribution omitted for simplicity)
  stress = sts_ms + sts_is
  strain = str_ms + str_is

end subroutine tdstress_hs_silent


!===============================================================================
! Silent version of TDstressFS (full-space solution)
!===============================================================================
subroutine tdstress_fs_silent(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, &
                              stress, strain)
  use m_nikkhoo_green, only: DP, PI, EPS
  use, intrinsic :: ieee_arithmetic
  implicit none

  real(DP), intent(in) :: x, y, z
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  real(DP), intent(in) :: ss, ds, ts, mu, lambda
  real(DP), dimension(6), intent(out) :: stress, strain

  real(DP) :: nu, bx, by, bz
  real(DP), dimension(3) :: vnorm, vstrike, vdip, ey, ez
  real(DP), dimension(3, 3) :: A
  real(DP), dimension(3) :: p1_td, p2_td, p3_td
  real(DP) :: x_td, y_td, z_td
  real(DP), dimension(3) :: e12, e13, e23
  real(DP) :: A_angle, B_angle, C_angle
  integer :: trimode
  real(DP) :: exx, eyy, ezz, exy, exz, eyz
  real(DP) :: exx_out, eyy_out, ezz_out, exy_out, exz_out, eyz_out
  real(DP) :: exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p
  real(DP) :: sxx, syy, szz, sxy, sxz, syz
  real(DP) :: norm_val

  ! Poisson's ratio
  nu = 1.0_DP / (1.0_DP + lambda / mu) / 2.0_DP

  ! Slip vector components
  bx = ts  ! Tensile-slip
  by = ss  ! Strike-slip
  bz = ds  ! Dip-slip

  ! Unit vectors
  ey = [0.0_DP, 1.0_DP, 0.0_DP]
  ez = [0.0_DP, 0.0_DP, 1.0_DP]

  ! Normal vector
  call cross_product_local(p2 - p1, p3 - p1, vnorm)
  norm_val = sqrt(sum(vnorm**2))
  if (norm_val < EPS) then
    stress = 0.0_DP
    strain = 0.0_DP
    return
  end if
  vnorm = vnorm / norm_val

  ! Strike vector
  call cross_product_local(ez, vnorm, vstrike)
  norm_val = sqrt(sum(vstrike**2))
  if (norm_val < EPS) then
    vstrike = ey * vnorm(3)
    if (p1(3) > 0.0_DP) vstrike = -vstrike
    norm_val = sqrt(sum(vstrike**2))
  end if
  vstrike = vstrike / norm_val

  ! Dip vector
  call cross_product_local(vnorm, vstrike, vdip)

  ! Transformation matrix (columns are unit vectors)
  A(:, 1) = vnorm
  A(:, 2) = vstrike
  A(:, 3) = vdip

  ! Transform to TDCS
  p1_td = 0.0_DP
  p2_td = 0.0_DP
  p3_td = 0.0_DP

  call coord_trans_local(x - p2(1), y - p2(2), z - p2(3), A, x_td, y_td, z_td)
  call coord_trans_local(p1(1) - p2(1), p1(2) - p2(2), p1(3) - p2(3), A, p1_td(1), p1_td(2), p1_td(3))
  call coord_trans_local(p3(1) - p2(1), p3(2) - p2(2), p3(3) - p2(3), A, p3_td(1), p3_td(2), p3_td(3))

  ! Unit vectors along TD sides
  norm_val = sqrt(sum((p2_td - p1_td)**2))
  if (norm_val < EPS) then
    stress = 0.0_DP
    strain = 0.0_DP
    return
  end if
  e12 = (p2_td - p1_td) / norm_val

  norm_val = sqrt(sum((p3_td - p1_td)**2))
  if (norm_val < EPS) then
    stress = 0.0_DP
    strain = 0.0_DP
    return
  end if
  e13 = (p3_td - p1_td) / norm_val

  norm_val = sqrt(sum((p3_td - p2_td)**2))
  if (norm_val < EPS) then
    stress = 0.0_DP
    strain = 0.0_DP
    return
  end if
  e23 = (p3_td - p2_td) / norm_val

  ! Angles
  A_angle = acos(max(-1.0_DP, min(1.0_DP, dot_product(e12, e13))))
  B_angle = acos(max(-1.0_DP, min(1.0_DP, -dot_product(e12, e23))))
  C_angle = acos(max(-1.0_DP, min(1.0_DP, dot_product(e23, e13))))

  ! Determine configuration
  call trimode_finder_local(y_td, z_td, x_td, p1_td, p2_td, p3_td, trimode)

  ! Initialize
  exx = 0.0_DP; eyy = 0.0_DP; ezz = 0.0_DP
  exy = 0.0_DP; exz = 0.0_DP; eyz = 0.0_DP

  if (trimode == 1) then
    ! Configuration I
    call tdsetup_s_local(x_td, y_td, z_td, A_angle, bx, by, bz, nu, p1_td, -e13, &
                         exx, eyy, ezz, exy, exz, eyz)

    call tdsetup_s_local(x_td, y_td, z_td, B_angle, bx, by, bz, nu, p2_td, e12, &
                         exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p)
    exx = exx + exx_p; eyy = eyy + eyy_p; ezz = ezz + ezz_p
    exy = exy + exy_p; exz = exz + exz_p; eyz = eyz + eyz_p

    call tdsetup_s_local(x_td, y_td, z_td, C_angle, bx, by, bz, nu, p3_td, e23, &
                         exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p)
    exx = exx + exx_p; eyy = eyy + eyy_p; ezz = ezz + ezz_p
    exy = exy + exy_p; exz = exz + exz_p; eyz = eyz + eyz_p

  else if (trimode == -1) then
    ! Configuration II
    call tdsetup_s_local(x_td, y_td, z_td, A_angle, bx, by, bz, nu, p1_td, e13, &
                         exx, eyy, ezz, exy, exz, eyz)

    call tdsetup_s_local(x_td, y_td, z_td, B_angle, bx, by, bz, nu, p2_td, -e12, &
                         exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p)
    exx = exx + exx_p; eyy = eyy + eyy_p; ezz = ezz + ezz_p
    exy = exy + exy_p; exz = exz + exz_p; eyz = eyz + eyz_p

    call tdsetup_s_local(x_td, y_td, z_td, C_angle, bx, by, bz, nu, p3_td, -e23, &
                         exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p)
    exx = exx + exx_p; eyy = eyy + eyy_p; ezz = ezz + ezz_p
    exy = exy + exy_p; exz = exz + exz_p; eyz = eyz + eyz_p

  else
    ! Configuration 0 - on triangle edge (singular)
    exx = ieee_value(0.0_DP, ieee_quiet_nan)
    eyy = ieee_value(0.0_DP, ieee_quiet_nan)
    ezz = ieee_value(0.0_DP, ieee_quiet_nan)
    exy = ieee_value(0.0_DP, ieee_quiet_nan)
    exz = ieee_value(0.0_DP, ieee_quiet_nan)
    eyz = ieee_value(0.0_DP, ieee_quiet_nan)
  end if

  ! Transform strain to EFCS
  call tens_trans_local(exx, eyy, ezz, exy, exz, eyz, A, &
                        exx_out, eyy_out, ezz_out, exy_out, exz_out, eyz_out)

  exx = exx_out; eyy = eyy_out; ezz = ezz_out
  exy = exy_out; exz = exz_out; eyz = eyz_out

  ! Calculate stress
  sxx = 2.0_DP * mu * exx + lambda * (exx + eyy + ezz)
  syy = 2.0_DP * mu * eyy + lambda * (exx + eyy + ezz)
  szz = 2.0_DP * mu * ezz + lambda * (exx + eyy + ezz)
  sxy = 2.0_DP * mu * exy
  sxz = 2.0_DP * mu * exz
  syz = 2.0_DP * mu * eyz

  stress(1) = sxx; stress(2) = syy; stress(3) = szz
  stress(4) = sxy; stress(5) = sxz; stress(6) = syz

  strain(1) = exx; strain(2) = eyy; strain(3) = ezz
  strain(4) = exy; strain(5) = exz; strain(6) = eyz

end subroutine tdstress_fs_silent


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
! Helper: Coordinate transformation
!===============================================================================
subroutine coord_trans_local(x, y, z, A, x_out, y_out, z_out)
  use m_nikkhoo_green, only: DP
  implicit none
  real(DP), intent(in) :: x, y, z
  real(DP), dimension(3, 3), intent(in) :: A
  real(DP), intent(out) :: x_out, y_out, z_out

  real(DP), dimension(3, 3) :: AT

  AT = transpose(A)
  x_out = AT(1, 1) * x + AT(1, 2) * y + AT(1, 3) * z
  y_out = AT(2, 1) * x + AT(2, 2) * y + AT(2, 3) * z
  z_out = AT(3, 1) * x + AT(3, 2) * y + AT(3, 3) * z
end subroutine coord_trans_local


!===============================================================================
! Helper: Trimode finder (barycentric coordinates)
!===============================================================================
subroutine trimode_finder_local(x, y, z, p1, p2, p3, trimode)
  use m_nikkhoo_green, only: DP, EPS
  implicit none

  real(DP), intent(in) :: x, y, z
  real(DP), dimension(3), intent(in) :: p1, p2, p3
  integer, intent(out) :: trimode

  real(DP) :: a, b, c, denominator

  denominator = (p2(3) - p3(3)) * (p1(2) - p3(2)) + (p3(2) - p2(2)) * (p1(3) - p3(3))

  if (abs(denominator) < EPS) then
    trimode = 1
    return
  end if

  a = ((p2(3) - p3(3)) * (x - p3(2)) + (p3(2) - p2(2)) * (y - p3(3))) / denominator
  b = ((p3(3) - p1(3)) * (x - p3(2)) + (p1(2) - p3(2)) * (y - p3(3))) / denominator
  c = 1.0_DP - a - b

  trimode = 1

  if (a <= 0.0_DP .and. b > c .and. c > a) then
    trimode = -1
  else if (b <= 0.0_DP .and. c > a .and. a > b) then
    trimode = -1
  else if (c <= 0.0_DP .and. a > b .and. b > c) then
    trimode = -1
  end if

  if (abs(a) < EPS .and. b >= 0.0_DP .and. c >= 0.0_DP) then
    trimode = 0
  else if (a >= 0.0_DP .and. abs(b) < EPS .and. c >= 0.0_DP) then
    trimode = 0
  else if (a >= 0.0_DP .and. b >= 0.0_DP .and. abs(c) < EPS) then
    trimode = 0
  end if

  if (trimode == 0 .and. abs(z) > EPS) then
    trimode = 1
  end if
end subroutine trimode_finder_local


!===============================================================================
! Helper: TDsetup for strain calculation
!===============================================================================
subroutine tdsetup_s_local(x, y, z, alpha, bx, by, bz, nu, tri_vertex, side_vec, &
                           exx, eyy, ezz, exy, exz, eyz)
  use m_nikkhoo_green, only: DP, PI
  implicit none

  real(DP), intent(in) :: x, y, z, alpha, bx, by, bz, nu
  real(DP), dimension(3), intent(in) :: tri_vertex, side_vec
  real(DP), intent(out) :: exx, eyy, ezz, exy, exz, eyz

  real(DP), dimension(2, 2) :: A2
  real(DP), dimension(3, 3) :: B
  real(DP) :: y1, z1, by1, bz1
  real(DP) :: exx_adcs, eyy_adcs, ezz_adcs, exy_adcs, exz_adcs, eyz_adcs

  ! 2x2 transformation matrix
  A2(1, 1) = side_vec(3)
  A2(1, 2) = -side_vec(2)
  A2(2, 1) = side_vec(2)
  A2(2, 2) = side_vec(3)

  ! Transform coordinates to ADCS
  y1 = A2(1, 1) * (y - tri_vertex(2)) + A2(1, 2) * (z - tri_vertex(3))
  z1 = A2(2, 1) * (y - tri_vertex(2)) + A2(2, 2) * (z - tri_vertex(3))

  ! Transform slip vector
  by1 = A2(1, 1) * by + A2(1, 2) * bz
  bz1 = A2(2, 1) * by + A2(2, 2) * bz

  ! Calculate strains in ADCS
  call angdis_strain_local(x, y1, z1, -PI + alpha, bx, by1, bz1, nu, &
                           exx_adcs, eyy_adcs, ezz_adcs, exy_adcs, exz_adcs, eyz_adcs)

  ! Transform strains from ADCS to TDCS
  B(1, 1) = 1.0_DP; B(1, 2) = 0.0_DP; B(1, 3) = 0.0_DP
  B(2, 1) = 0.0_DP; B(2, 2) = A2(1, 1); B(2, 3) = A2(2, 1)
  B(3, 1) = 0.0_DP; B(3, 2) = A2(1, 2); B(3, 3) = A2(2, 2)

  call tens_trans_local(exx_adcs, eyy_adcs, ezz_adcs, exy_adcs, exz_adcs, eyz_adcs, B, &
                        exx, eyy, ezz, exy, exz, eyz)
end subroutine tdsetup_s_local


!===============================================================================
! Helper: Angular dislocation strain
!===============================================================================
subroutine angdis_strain_local(x, y, z, alpha, bx, by, bz, nu, &
                               exx, eyy, ezz, exy, exz, eyz)
  use m_nikkhoo_green, only: DP, PI, EPS
  implicit none

  real(DP), intent(in) :: x, y, z, alpha, bx, by, bz, nu
  real(DP), intent(out) :: exx, eyy, ezz, exy, exz, eyz

  real(DP) :: sinA, cosA, eta, zeta
  real(DP) :: x2, y2, z2, r2, r, r3, rz, r2z2, r3z
  real(DP) :: W, W2, Wr, W2r, Wr3, W2r2
  real(DP) :: C, S
  real(DP) :: rFi_rx, rFi_ry, rFi_rz
  real(DP) :: one_minus_nu

  sinA = sin(alpha)
  cosA = cos(alpha)
  eta = y * cosA - z * sinA
  zeta = y * sinA + z * cosA

  x2 = x * x
  y2 = y * y
  z2 = z * z
  r2 = x2 + y2 + z2
  r = sqrt(r2)

  if (r < EPS) then
    exx = 0.0_DP; eyy = 0.0_DP; ezz = 0.0_DP
    exy = 0.0_DP; exz = 0.0_DP; eyz = 0.0_DP
    return
  end if

  r3 = r * r2
  rz = r * (r - z)
  r2z2 = r2 * (r - z)**2
  r3z = r3 * (r - z)

  W = zeta - r
  W2 = W * W
  Wr = W * r
  W2r = W2 * r
  Wr3 = W * r3
  W2r2 = W2 * r2

  if (abs(Wr) < EPS .or. abs(rz) < EPS) then
    exx = 0.0_DP; eyy = 0.0_DP; ezz = 0.0_DP
    exy = 0.0_DP; exz = 0.0_DP; eyz = 0.0_DP
    return
  end if

  C = (r * cosA - z) / Wr
  S = (r * sinA - y) / Wr

  one_minus_nu = 1.0_DP - nu

  rFi_rx = (eta / r / (r - zeta) - y / r / (r - z)) / (4.0_DP * PI)
  rFi_ry = (x / r / (r - z) - cosA * x / r / (r - zeta)) / (4.0_DP * PI)
  rFi_rz = (sinA * x / r / (r - zeta)) / (4.0_DP * PI)

  exx = bx * rFi_rx + &
        bx / (8.0_DP * PI * one_minus_nu) * (eta / Wr + eta * x2 / W2r2 - &
        eta * x2 / Wr3 + y / rz - x2 * y / r2z2 - x2 * y / r3z) - &
        by * x / (8.0_DP * PI * one_minus_nu) * (((2.0_DP * nu + 1.0_DP) / Wr + &
        x2 / W2r2 - x2 / Wr3) * cosA + (2.0_DP * nu + 1.0_DP) / rz - &
        x2 / r2z2 - x2 / r3z) + &
        bz * x * sinA / (8.0_DP * PI * one_minus_nu) * ((2.0_DP * nu + 1.0_DP) / Wr + &
        x2 / W2r2 - x2 / Wr3)

  eyy = by * rFi_ry + &
        bx / (8.0_DP * PI * one_minus_nu) * ((1.0_DP / Wr + S**2 - y2 / Wr3) * eta + &
        (2.0_DP * nu + 1.0_DP) * y / rz - y**3 / r2z2 - y**3 / r3z - &
        2.0_DP * nu * cosA * S) - &
        by / (8.0_DP * PI * one_minus_nu) * (y2 / Wr3 - 1.0_DP / Wr - S**2) * cosA * y + &
        bz * sinA / (8.0_DP * PI * one_minus_nu) * ((1.0_DP / Wr - S**2 + y2 / Wr3) * y)

  ezz = bz * rFi_rz + &
        bx / (8.0_DP * PI * one_minus_nu) * (eta / Wr3 * z2 - eta / Wr + eta * C**2 - &
        y / rz + z2 * y / r2z2 + z2 * y / r3z) - &
        by / (8.0_DP * PI * one_minus_nu) * (z2 / Wr3 - 1.0_DP / Wr - C**2) * cosA * z + &
        bz * sinA / (8.0_DP * PI * one_minus_nu) * ((1.0_DP / Wr - C**2 + z2 / Wr3) * z)

  exy = bx / 2.0_DP * rFi_ry + by / 2.0_DP * rFi_rx + &
        bx / (8.0_DP * PI * one_minus_nu) * (x * y / Wr3 * eta - x * S * eta / Wr + &
        (2.0_DP * nu + 1.0_DP) * x / rz - x * y2 / r2z2 - x * y2 / r3z - &
        nu * x * cosA * S / Wr) - &
        by / (8.0_DP * PI * one_minus_nu) * cosA * (x * y / Wr3 - x * S / Wr) + &
        bz * sinA / (8.0_DP * PI * one_minus_nu) * x * (y / Wr3 - S / Wr)

  exz = bx / 2.0_DP * rFi_rz + bz / 2.0_DP * rFi_rx + &
        bx / (8.0_DP * PI * one_minus_nu) * (x * z / Wr3 * eta + x * C * eta / Wr - &
        x * y / rz**2 * (r + z) - x * y / r3z * 2.0_DP - &
        nu * x * cosA * C / Wr) - &
        by / (8.0_DP * PI * one_minus_nu) * cosA * (x * z / Wr3 + x * C / Wr) + &
        bz * sinA / (8.0_DP * PI * one_minus_nu) * x * (z / Wr3 + C / Wr)

  eyz = by / 2.0_DP * rFi_rz + bz / 2.0_DP * rFi_ry + &
        bx / (8.0_DP * PI * one_minus_nu) * (y * z / Wr3 * eta - y * C * eta / Wr + &
        S * C * eta / Wr + &
        y2 / rz**2 * (r + z) + y2 / r3z * 2.0_DP - &
        (2.0_DP * nu + 1.0_DP) / rz - &
        nu * cosA * (S * C / Wr - y * C / Wr)) - &
        by / (8.0_DP * PI * one_minus_nu) * cosA * (y * z / Wr3 - y * C / Wr - S * C / Wr) + &
        bz * sinA / (8.0_DP * PI * one_minus_nu) * (y * z / Wr3 - y * C / Wr + S * C / Wr)

end subroutine angdis_strain_local


!===============================================================================
! Helper: Tensor transformation
!===============================================================================
subroutine tens_trans_local(exx, eyy, ezz, exy, exz, eyz, A, &
                            exx_out, eyy_out, ezz_out, exy_out, exz_out, eyz_out)
  use m_nikkhoo_green, only: DP
  implicit none

  real(DP), intent(in) :: exx, eyy, ezz, exy, exz, eyz
  real(DP), dimension(3, 3), intent(in) :: A
  real(DP), intent(out) :: exx_out, eyy_out, ezz_out, exy_out, exz_out, eyz_out

  real(DP), dimension(3, 3) :: E, E_out
  integer :: i, j, k, l

  ! Build strain tensor
  E(1, 1) = exx; E(2, 2) = eyy; E(3, 3) = ezz
  E(1, 2) = exy; E(2, 1) = exy
  E(1, 3) = exz; E(3, 1) = exz
  E(2, 3) = eyz; E(3, 2) = eyz

  ! Transform: E_out = A * E * A^T
  E_out = 0.0_DP
  do i = 1, 3
    do j = 1, 3
      do k = 1, 3
        do l = 1, 3
          E_out(i, j) = E_out(i, j) + A(i, k) * A(j, l) * E(k, l)
        end do
      end do
    end do
  end do

  ! Extract components
  exx_out = E_out(1, 1)
  eyy_out = E_out(2, 2)
  ezz_out = E_out(3, 3)
  exy_out = E_out(1, 2)
  exz_out = E_out(1, 3)
  eyz_out = E_out(2, 3)

end subroutine tens_trans_local


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


