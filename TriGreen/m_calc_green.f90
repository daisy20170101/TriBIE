module m_calc_green
public

    integer, parameter :: DP=kind(1.d0)
    real(DP), parameter          ::  PI   = 3.1415926536d0
    real(DP), parameter          ::  ZERO = 1d-26

    ! gts data
    integer, parameter          ::  n_max_vertex_number    = 480000, &
                                    n_max_edge_number      = 240000, &
                                    n_max_cell_number      = 240000

!    integer                     ::  n_vertex, n_edge, n_cell
 !   real(DP),DIMENSION(:,:),ALLOCATABLE  ::  arr_vertex
 !   integer,DIMENSION(:,:),ALLOCATABLE  ::  arr_edge
 !   integer,DIMENSION(:,:),ALLOCATABLE  ::  arr_cell

    ! global parameters
    real(DP), parameter          :: subd_az = 0.d0  ! Subduction azimuth (degree)
    real(DP), parameter          :: rot_deg = 0.d0   ! Should be the same as used in mesh_gen
    real(DP), parameter          ::  vpl1 = cos((90.d0 - subd_az+rot_deg)*PI/180.d0)   ! coord-1's component (east)
    real(DP),  parameter          ::  vpl2 = -sin((90.d0 - subd_az+rot_deg)*PI/180.d0)  ! coord-2's component (south)

    real(DP), parameter          ::  parm_nu  = (6.1d0**2-2.d0*3.5d0**2)/(2.d0*6.1d0**2-2.d0*3.5d0**2)
    real(DP), parameter          ::  parm_miu = 30000.d0        ! rigidity (MPa)
    real(DP), parameter          ::  parm_l   = ( 2.d0*parm_nu*parm_miu ) / (1.d0 - 2.d0*parm_nu)   ! lambda in Lame parameters (MPa)

end module m_calc_green
