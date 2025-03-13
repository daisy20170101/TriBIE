module m_calc_green
public

    real(8), parameter          ::  PI   = 3.1415926536d0
    real(8), parameter          ::  ZERO = 1d-6

    ! gts data
    integer, parameter          ::  dp=kind(1.d0), n_max_vertex_number    = 480000, &
                                    n_max_edge_number      = 240000, &
                                    n_max_cell_number      = 240000

!    integer                     ::  n_vertex, n_edge, n_cell
 !   real(8),DIMENSION(:,:),ALLOCATABLE  ::  arr_vertex
 !   integer,DIMENSION(:,:),ALLOCATABLE  ::  arr_edge
 !   integer,DIMENSION(:,:),ALLOCATABLE  ::  arr_cell

    ! global parameters
    real(8), parameter          :: subd_az = 0.d0  ! Subduction azimuth (degree)
    real(8), parameter          :: rot_deg = 0.d0   ! Should be the same as used in mesh_gen
    real(8), parameter          ::  vpl1 = cos((90.d0 - subd_az+rot_deg)*PI/180.d0)   ! coord-1's component (east)
    real(8),  parameter          ::  vpl2 = -sin((90.d0 - subd_az+rot_deg)*PI/180.d0)  ! coord-2's component (south)

    real(8), parameter          ::  parm_nu  = (6.1d0**2-2.d0*3.5d0**2)/(2.d0*6.1d0**2-2.d0*3.5d0**2)
    real(8), parameter          ::  parm_miu = 30000.d0        ! rigidity (MPa)
    real(8), parameter          ::  parm_l   = ( 2.d0*parm_nu*parm_miu ) / (1.d0 - 2.d0*parm_nu)   ! lambda in Lame parameters (MPa)

end module m_calc_green
