!===============================================================================
! m_nikkhoo_green.f90
!
! Parameter module for Nikkhoo & Walter (2015) triangular dislocation
! Green's function calculation with MPI/OpenMP parallelization.
!
! This module contains global parameters used by calc_nikkhoo.f90
!===============================================================================

module m_nikkhoo_green
  implicit none
  public

  integer, parameter :: DP = kind(1.d0)
  real(DP), parameter :: PI = 3.141592653589793238462643383279502884197_DP
  real(DP), parameter :: ZERO = 1.0d-26
  real(DP), parameter :: EPS = 1.0d-15

  ! Maximum array sizes for mesh data
  integer, parameter :: n_max_vertex_number = 480000
  integer, parameter :: n_max_edge_number = 240000
  integer, parameter :: n_max_cell_number = 240000

  ! Subduction geometry parameters
  real(DP), parameter :: subd_az = 0.d0  ! Subduction azimuth (degree)
  real(DP), parameter :: rot_deg = 0.d0  ! Rotation angle (same as mesh_gen)
  real(DP), parameter :: vpl1 = cos((90.d0 - subd_az + rot_deg) * PI / 180.d0)  ! East component
  real(DP), parameter :: vpl2 = -sin((90.d0 - subd_az + rot_deg) * PI / 180.d0) ! South component

  ! Material parameters (elastic half-space)
  ! These can be modified based on your specific problem
  real(DP), parameter :: parm_nu = (6.1d0**2 - 2.d0*3.5d0**2) / (2.d0*6.1d0**2 - 2.d0*3.5d0**2)  ! Poisson's ratio
  real(DP), parameter :: parm_miu = 30000.d0  ! Rigidity (MPa)
  real(DP), parameter :: parm_l = (2.d0 * parm_nu * parm_miu) / (1.d0 - 2.d0 * parm_nu)  ! Lambda (Lame parameter, MPa)

end module m_nikkhoo_green
