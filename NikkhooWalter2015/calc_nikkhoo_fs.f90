!===============================================================================
! calc_nikkhoo_fs.f90
!
! MPI/OpenMP parallel program for calculating Green's function coefficients
! (stiffness matrices) for a planar fault in an elastic WHOLE-SPACE, using
! the Nikkhoo & Walter (2015) full-space triangular dislocation solution
! (tdstress_fs). Adapted from calc_nikkhoo.f90 (which uses the half-space
! solution tdstress_hs) for SEAS Benchmark BP8.
!
! Differences from calc_nikkhoo.f90:
!   - Calls tdstress_fs directly (no free-surface image dislocation, no
!     z<0 restriction) -- see NikkhooWalter2015/README_NIKKHOO.md for the
!     validation of tdstress_fs against an independent Python port.
!   - No z-coordinate sign flip on load: the mesh is expected to already be
!     a flat plane at z=0 (mesh x,y map to BP8's fault-plane coordinates
!     x2 (strike), x3 (dip); the whole-space fault normal is BP8's x1,
!     which the code's local-frame convention places along its own z).
!     A flat mesh at z=0 makes every triangle hit calc_local_coordinate's
!     "horizontal triangle" branch, giving an unambiguous local frame
!     a1=(1,0,0) [x2], a2=(0,1,0) [x3], a3=(0,0,1) [x1/normal] -- so no
!     plate-convergence-direction heuristic is needed to pick a frame.
!   - Elastic constants are BP8's own (mu=32.04 GPa, nu=0.25), not the
!     m_nikkhoo_green module's BP5/6 constants -- kept as local parameters
!     here so the shared module used by calc_nikkhoo.f90 is untouched.
!   - Produces the full 2x2 in-plane shear stiffness system: trigreen_22
!     (x2 traction from unit x2 slip), trigreen_23 (x2 traction from unit
!     x3 slip), trigreen_32 (x3 traction from unit x2 slip), trigreen_33
!     (x3 traction from unit x3 slip). The cross terms (23/32) were checked
!     numerically (debug/validate_calc_nikkhoo_fs.py) and are NOT
!     negligible in general -- an earlier draft of this driver assumed
!     planar-fault symmetry would zero them out; it does not, so all four
!     are computed.
!
!     Pointwise, K23(r) and K32(r) are the SAME function of source-receiver
!     offset r (G_23=G_32 by the symmetry of the elasticity tensor -- this
!     is verified for a single pair in debug/validate_calc_nikkhoo_fs.py).
!     That does NOT mean the assembled matrices are symmetric under
!     swapping which triangle is source vs. receiver (K22(i,j) vs
!     K22(j,i)): this driver, like calc_nikkhoo.f90, evaluates the
!     traction at a single point (the receiver centroid) rather than
!     area-averaging it over the receiver triangle, and centroid-point
!     collocation does not preserve exact Betti reciprocity between
!     unequal source/receiver geometry -- only the true area-integrated
!     (Galerkin) formulation does. The resulting asymmetry is a standard,
!     well-understood BEM discretization artifact (largest for close
!     neighbors, ~1% of the self-term at 10 m spacing here, shrinking as
!     element size shrinks relative to source-receiver distance) already
!     present in calc_nikkhoo.f90's HS matrices; it is not a bug.
!
!     No normal-stress matrix is computed: BP8 Eq. (7) fixes total normal
!     stress in time (only pore pressure perturbs the effective normal
!     stress), so it is not needed for this benchmark.
!
!   - On the flat z=0 mesh, traction on the fault plane is (Sxz, Syz, Szz)
!     (index 5, 6, 3 of the (Sxx,Syy,Szz,Sxy,Sxz,Syz) stress vector) --
!     NOT Sxy, which is an in-plane stress component with no meaning as a
!     fault-plane traction. Physical slip in EFCS x (=BP8 x2) or EFCS y
!     (=BP8 x3) must also be translated into tdstress_fs's (ss,ds,ts)
!     arguments, which are defined relative to that subroutine's own
!     internally-computed Vstrike/Vdip -- for a horizontal (z=const)
!     triangle those come out as Vstrike=+y, Vdip=-x (see TDstressFS.m's
!     degenerate-normal branch, Vstrike=cross(eZ,Vnorm) vanishes when
!     Vnorm=+z so it falls back to Vstrike=eY*Vnorm(3)=+y, then
!     Vdip=cross(Vnorm,Vstrike)=-x). So unit EFCS-x slip is ds=-1 (not
!     ss=1), and unit EFCS-y slip is ss=1. Confirmed against
!     debug/validate_calc_nikkhoo_fs.py.
!
! Usage:
!   mpirun -np <nprocs> ./calc_nikkhoo_fs
!
! Input:
!   triangular_mesh.gts - GTS format mesh file, flat at z=0
!
! Output:
!   trigreen_22_<rank>.bin - x2 traction from unit x2 slip (one per rank)
!   trigreen_23_<rank>.bin - x2 traction from unit x3 slip
!   trigreen_32_<rank>.bin - x3 traction from unit x2 slip
!   trigreen_33_<rank>.bin - x3 traction from unit x3 slip
!   position.bin           - centroid positions of all elements
!===============================================================================

module m_nikkhoo_fs_params
  implicit none
  public

  integer, parameter :: DP = kind(1.d0)
  real(DP), parameter :: EPS = 1.0d-15

  ! BP8 elastic constants (Table 1): mu = 32.04 GPa, nu = 0.25.
  ! Working in MPa to match the magnitude convention of m_nikkhoo_green.
  real(DP), parameter :: fs_mu = 32040.d0                                 ! MPa
  real(DP), parameter :: fs_nu = 0.25d0
  real(DP), parameter :: fs_lambda = 2.d0*fs_nu*fs_mu/(1.d0-2.d0*fs_nu)   ! MPa

end module m_nikkhoo_fs_params


program calc_nikkhoo_fs
  use m_nikkhoo_fs_params
  use mpi
  use omp_lib
  implicit none

  character(*), parameter :: fname = "triangular_mesh.gts"

  integer :: n_vertex, n_edge, n_cell
  real(DP), dimension(:,:), allocatable :: arr_vertex
  integer, dimension(:,:), allocatable :: arr_edge
  integer, dimension(:,:), allocatable :: arr_cell

  integer :: ierr, size, myid, master
  real(DP) :: start_time, end_time

  integer :: Nt, Nt_all, local_cells, cells_processed
  integer :: base_cells, extra_cells

  integer :: num_threads

  logical :: error_occurred
  character(len=256) :: error_message

  error_occurred = .false.
  error_message = ""

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) then
    write(*,*) "Error: Failed to initialize MPI"
    stop
  end if

  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
  master = 0

  num_threads = omp_get_max_threads()
  call omp_set_num_threads(num_threads)

  if (myid == master) then
    start_time = MPI_Wtime()
    write(*,*) "========================================================"
    write(*,*) "Nikkhoo-Walter FULL-SPACE Stiffness Calculation (BP8)"
    write(*,*) "========================================================"
    write(*,*) "MPI processes:", size
    write(*,*) "OpenMP threads per process:", num_threads
    write(*,*) "mu (MPa) =", fs_mu, " nu =", fs_nu, " lambda (MPa) =", fs_lambda
    write(*,*) ""
  end if

  if (.not. error_occurred) then
    call load_name_fs(fname, n_vertex, n_edge, n_cell, error_occurred, error_message)
  end if

  if (.not. error_occurred) then
    allocate(arr_vertex(n_vertex, 3), arr_edge(n_edge, 2), arr_cell(n_cell, 3), stat=ierr)
    if (ierr /= 0) then
      error_occurred = .true.
      error_message = "Failed to allocate mesh arrays"
    end if
  end if

  if (.not. error_occurred) then
    call load_gts_fs(fname, n_vertex, n_edge, n_cell, arr_vertex, arr_edge, arr_cell, &
                     error_occurred, error_message)
  end if

  if (.not. error_occurred) then
    if (myid == master) then
      write(*,*) "Mesh loaded successfully:"
      write(*,*) "  Vertices:", n_vertex
      write(*,*) "  Cells:", n_cell
      write(*,*) ""
    end if

    Nt_all = n_cell

    if (size == 1) then
      base_cells = n_cell
      extra_cells = 0
      local_cells = n_cell
      cells_processed = n_cell
      Nt = n_cell
    else
      base_cells = Nt_all / size
      extra_cells = mod(Nt_all, size)
      if (myid < extra_cells) then
        local_cells = base_cells + 1
      else
        local_cells = base_cells
      end if
      cells_processed = local_cells
      Nt = local_cells
    end if

    call calc_nikkhoo_fs_allcell(myid, size, Nt, arr_vertex, arr_cell, &
                                 n_vertex, n_cell, cells_processed, base_cells, extra_cells, &
                                 error_occurred, error_message)
  end if

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
      write(*,*) "Calculation completed successfully! Total time:", end_time - start_time, "seconds"
      write(*,*) "Output files: trigreen_{22,23,32,33}_<rank>.bin, position.bin"
    end if
  end if

  call MPI_finalize(ierr)
  if (error_occurred) stop 1

end program calc_nikkhoo_fs


subroutine load_name_fs(fname, n_vertex, n_edge, n_cell, error_occurred, error_message)
  implicit none
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
end subroutine load_name_fs


subroutine load_gts_fs(fname, n_vertex, n_edge, n_cell, arr_vertex, arr_edge, arr_cell, &
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

  ! No z-sign flip here: unlike the half-space driver, the full-space
  ! solution has no z<0 requirement, and the BP8 mesh is generated flat at
  ! z=0 by make_bp8_mesh.py.
  do i = 1, n_vertex
    read(10, *, iostat=iostat) arr_vertex(i, 1), arr_vertex(i, 2), arr_vertex(i, 3)
    if (iostat /= 0) then
      error_occurred = .true.
      error_message = "Failed to read vertex data"
      close(10)
      return
    end if
  end do

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
end subroutine load_gts_fs


subroutine calc_nikkhoo_fs_allcell(myid, size, Nt, arr_vertex, arr_cell, &
                                   n_vertex, n_cell, cells_processed, base_cells, extra_cells, &
                                   error_occurred, error_message)
  use m_nikkhoo_fs_params
  use nikkhoo_walter, only: tdstress_fs
  use omp_lib
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, intent(in) :: myid, size, Nt, n_vertex, n_cell
  integer, intent(in) :: cells_processed, base_cells, extra_cells
  real(DP), intent(in) :: arr_vertex(n_vertex, 3)
  integer, intent(in) :: arr_cell(n_cell, 3)
  logical, intent(inout) :: error_occurred
  character(len=*), intent(inout) :: error_message

  integer :: i, j, k
  integer :: vj(3)
  real(DP) :: p1(3), p2(3), p3(3), co(3)
  real(DP) :: src_p1(3), src_p2(3), src_p3(3)
  real(DP) :: stress_x2(6), strain_x2(6), stress_x3(6), strain_x3(6)
  character(20) :: cTemp

  real(DP), allocatable :: arr_co(:,:), arr_trid(:,:)
  real(DP), allocatable :: arr_22(:,:), arr_23(:,:), arr_32(:,:), arr_33(:,:)

  integer :: local_cells, start_idx

  ! Unit slip in EFCS x (BP8 x2) and EFCS y (BP8 x3). tdstress_fs's (ss,ds,ts)
  ! arguments are relative to ITS OWN internal Vstrike/Vdip, which for a
  ! horizontal (z=const) triangle come out as Vstrike=+y, Vdip=-x (see
  ! module header and debug/validate_calc_nikkhoo_fs.py) -- so unit EFCS-x
  ! slip is (ss=0,ds=-1,ts=0) and unit EFCS-y slip is (ss=1,ds=0,ts=0).
  real(DP), parameter :: zero_slip = 0.d0

  local_cells = cells_processed

  if (myid < extra_cells) then
    start_idx = myid * (base_cells + 1) + 1
  else
    start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells + 1
  end if

  allocate(arr_co(3, max(1, local_cells)))
  allocate(arr_trid(9, n_cell))
  allocate(arr_22(max(1, local_cells), n_cell))
  allocate(arr_23(max(1, local_cells), n_cell))
  allocate(arr_32(max(1, local_cells), n_cell))
  allocate(arr_33(max(1, local_cells), n_cell))

  do k = 1, n_cell
    vj(1:3) = arr_cell(k, 1:3)
    arr_trid(1:3, k) = arr_vertex(vj(1), 1:3)
    arr_trid(4:6, k) = arr_vertex(vj(2), 1:3)
    arr_trid(7:9, k) = arr_vertex(vj(3), 1:3)
  end do

  if (local_cells > 0) then
    do j = 1, local_cells
      k = start_idx + j - 1
      if (k > n_cell) cycle
      p1(1:3) = arr_trid(1:3, k)
      p2(1:3) = arr_trid(4:6, k)
      p3(1:3) = arr_trid(7:9, k)
      co = (p1 + p2 + p3) / 3.d0
      arr_co(1:3, j) = co(1:3)
      arr_22(j, :) = 0.d0
      arr_23(j, :) = 0.d0
      arr_32(j, :) = 0.d0
      arr_33(j, :) = 0.d0
    end do
  else
    arr_22(1, :) = 0.d0
    arr_23(1, :) = 0.d0
    arr_32(1, :) = 0.d0
    arr_33(1, :) = 0.d0
    arr_co(:, 1) = 0.d0
  end if

  write(*,*) "Process", myid, "starting full-space stiffness calculation for", local_cells, "cells"

  if (local_cells > 0) then
    !$OMP PARALLEL DO PRIVATE(i, j, k, src_p1, src_p2, src_p3, stress_x2, strain_x2, stress_x3, strain_x3) &
    !$OMP& SHARED(arr_co, arr_trid, arr_22, arr_23, arr_32, arr_33, n_cell, local_cells, start_idx)
    do j = 1, local_cells
      k = start_idx + j - 1
      if (k > n_cell) cycle

      do i = 1, n_cell
        src_p1(1:3) = arr_trid(1:3, i)
        src_p2(1:3) = arr_trid(4:6, i)
        src_p3(1:3) = arr_trid(7:9, i)

        ! Response to unit EFCS-x (BP8 x2) slip: ss=0, ds=-1, ts=0.
        ! On the z=0 plane, traction is (Sxz,Syz,Szz) -> x2-traction is
        ! Sxz=stress(5), x3-traction is Syz=stress(6).
        call tdstress_fs(arr_co(1, j), arr_co(2, j), arr_co(3, j), &
                         src_p1, src_p2, src_p3, zero_slip, -1.d0, zero_slip, &
                         fs_mu, fs_lambda, stress_x2, strain_x2)

        ! Response to unit EFCS-y (BP8 x3) slip: ss=1, ds=0, ts=0.
        call tdstress_fs(arr_co(1, j), arr_co(2, j), arr_co(3, j), &
                         src_p1, src_p2, src_p3, 1.d0, zero_slip, zero_slip, &
                         fs_mu, fs_lambda, stress_x3, strain_x3)

        if (ieee_is_nan(stress_x2(5))) then
          arr_22(j, i) = 0.d0
        else
          arr_22(j, i) = stress_x2(5)
        end if

        if (ieee_is_nan(stress_x2(6))) then
          arr_32(j, i) = 0.d0
        else
          arr_32(j, i) = stress_x2(6)
        end if

        if (ieee_is_nan(stress_x3(5))) then
          arr_23(j, i) = 0.d0
        else
          arr_23(j, i) = stress_x3(5)
        end if

        if (ieee_is_nan(stress_x3(6))) then
          arr_33(j, i) = 0.d0
        else
          arr_33(j, i) = stress_x3(6)
        end if
      end do
    end do
    !$OMP END PARALLEL DO
  end if

  write(*,*) "Process", myid, "completed calculation"

  write(cTemp, *) myid

  if (local_cells > 0) then
    open(14, file='trigreen_22_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    do i = 1, local_cells
      write(14) arr_22(i, :)
    end do
    close(14)

    open(15, file='trigreen_23_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    do i = 1, local_cells
      write(15) arr_23(i, :)
    end do
    close(15)

    open(16, file='trigreen_32_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    do i = 1, local_cells
      write(16) arr_32(i, :)
    end do
    close(16)

    open(17, file='trigreen_33_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    do i = 1, local_cells
      write(17) arr_33(i, :)
    end do
    close(17)
  else
    open(14, file='trigreen_22_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    write(14) (0.d0, i=1, n_cell)
    close(14)

    open(15, file='trigreen_23_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    write(15) (0.d0, i=1, n_cell)
    close(15)

    open(16, file='trigreen_32_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    write(16) (0.d0, i=1, n_cell)
    close(16)

    open(17, file='trigreen_33_'//trim(adjustl(cTemp))//'.bin', form='unformatted', access='stream')
    write(17) (0.d0, i=1, n_cell)
    close(17)
  end if

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

  deallocate(arr_co, arr_trid, arr_22, arr_23, arr_32, arr_33)

end subroutine calc_nikkhoo_fs_allcell
