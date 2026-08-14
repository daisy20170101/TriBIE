!===============================================================================
! module_bp8.f90
!
! Shared state/parameters for src/3dtri_BP8.f90 (SEAS Benchmark BP8-QD-GS).
!
! Unlike phy3d_module_bp6.f90 (BP5/BP6: spatially-varying a-b profiles read
! from var-*.dat files, single scalar shear-traction/slip-rate DOF per
! element via the phy1/phy2 direction-cosine trick), BP8's friction
! parameters are spatially uniform (Table 1) and the friction law is
! genuinely 2-component (V2,V3 solved jointly, not decoupled) -- see the
! derivation in 3dtri_BP8.f90's derivs subroutine docstring for why.
!===============================================================================
module bp8_module
  implicit none
  public

  integer, parameter :: DP = kind(1.d0)
  real(DP), parameter :: PI_BP8 = 3.141592653589793238462643383279502884197_DP

  ! MPI / load balancing (mirrors calc_nikkhoo_fs.f90's decomposition so
  ! trigreen_{22,23,32,33}_<rank>.bin chunks line up with this driver's
  ! per-rank element ranges when run with the same -np).
  integer :: nprocs
  integer :: master_id
  integer :: my_start_idx, my_local_cells
  integer :: n_state

  ! Table 1 parameters (read from parameter1.txt, see example4/parameter1.txt)
  real(DP) :: xmu, xnu, rho, cs, eta          ! elastic constants, radiation damping eta=mu/2cs
  real(DP) :: seff0, tauinit                  ! initial effective normal stress, initial shear stress
  real(DP) :: fric_a, fric_b, fric_Drs, fric_Vstar, fric_fstar
  real(DP) :: Vinit, Vzero
  real(DP) :: Q0inj, Lfwid, toff, Lgauss      ! injection
  real(DP) :: alpha_diff, beta_diff, phi_diff ! pore pressure diffusion
  real(DP) :: tf_end, tint_out
  real(DP) :: lf_half, dz_cell                ! frictional half-length, cell size (mesh = exactly Omega_f)
  integer :: n_side                           ! mesh is n_side x n_side cells (2*n_side^2 triangles)

  character(len=200) :: foldername, stiffname

  ! Per-rank element arrays (local_cells long)
  real(DP), dimension(:,:), allocatable :: K22, K23, K32, K33  ! (local_cells, Nt_all)
  real(DP), dimension(:), allocatable :: cx2, cx3              ! local element centroid (x2,x3)

  ! Full arrays (master only, gathered every accepted step for output/interaction)
  real(DP), dimension(:), allocatable :: cx2_all, cx3_all

  ! Precomputed pore-pressure history at every element centroid, hourly
  ! snapshots (see pf_history_init in 3dtri_BP8.f90): p_hist(0:n_hours, Nt_all)
  real(DP), dimension(:,:), allocatable :: p_hist
  integer :: n_hours_hist
  real(DP) :: dt_hist

end module bp8_module
