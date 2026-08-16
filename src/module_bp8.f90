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

  ! BP8 Table 1: half-length of the rate-and-state fault, lf = 400 m. This
  ! is a FIXED benchmark parameter, not a simulation setting -- it must
  ! stay 400 m regardless of the mesh's own domain size or resolution, so
  ! it is a compile-time constant here rather than read from
  ! parameter1.txt (removing any chance of it silently drifting to match
  ! whatever mesh happens to be loaded, which is exactly the bug this
  ! constant fixes: an earlier version tied the frictional domain and the
  ! pore-pressure diffusion domain to the mesh's own extent
  ! (n_side*dz_cell/2), so a mesh bigger than Omega_f silently gave every
  ! extra element rate-and-state friction and diffusion instead of the
  ! locked (V=0, Eq. 13) boundary condition BP8 actually specifies there.
  real(DP), parameter :: lf_fixed = 400.0_DP

  ! BP8 Section 4.3: profile-line output nodes are required at EXACTLY
  ! 10 m spacing from -400 to 400 m (81 nodes), regardless of the mesh's
  ! own cell size -- so these are fixed too, not derived from n_side/dz_cell.
  integer, parameter :: n_nodes_fixed = 81
  real(DP), parameter :: node_dz_fixed = 10.0_DP

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
  real(DP) :: tau0_2, tau0_3  ! fixed reference traction (Eq. 29), set in set_initial_conditions
  real(DP) :: Q0inj, Lfwid, toff, Lgauss      ! injection
  real(DP) :: alpha_diff, beta_diff, phi_diff ! pore pressure diffusion
  real(DP) :: tf_end, tint_out
  real(DP) :: lf_half, dz_cell                ! MESH's own half-domain and cell size (may exceed lf_fixed)
  integer :: n_side                           ! mesh is n_side x n_side cells (2*n_side^2 triangles)

  character(len=200) :: foldername, stiffname

  ! Per-rank element arrays (local_cells long)
  real(DP), dimension(:,:), allocatable :: K22, K23, K32, K33  ! (local_cells, Nt_all)
  logical, dimension(:), allocatable :: is_active              ! true if |cx2|<lf_fixed and |cx3|<lf_fixed (local elements)

  ! Full arrays (every rank has the complete set, broadcast once at startup)
  real(DP), dimension(:), allocatable :: cx2_all, cx3_all

  ! Precomputed pore-pressure history at every element centroid, hourly
  ! snapshots (see pf_history_init in 3dtri_BP8.f90): p_hist(0:n_hours, Nt_all)
  real(DP), dimension(:,:), allocatable :: p_hist
  integer :: n_hours_hist
  real(DP) :: dt_hist

end module bp8_module
