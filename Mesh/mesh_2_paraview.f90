!! modified for paraview input by Duo Li

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Generate Nankai Model's mesh from gradient depth line
!
!   Version:
!       2013-07-24 Modified to use memory allocation
!       2012-07-17 Mesh generated from grid data
!       2007-09-11 create Tonankai model
!       2007-06-20 Change to model II
!       2006-11-26 Create this program
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine load_gts(fname,n_vertex,n_edge,n_cell,arr_vertex,arr_edge,arr_cell)
  implicit none
    
  character*(*)               ::  fname
  integer                     ::  i,n_vertex,n_edge,n_cell
  real,dimension(:,:) :: arr_vertex(n_vertex,3)
  integer,dimension(:,:) :: arr_edge(n_edge,2),arr_cell(n_cell,3)
  
  open(unit=10, file=fname)
  read(10, *) n_vertex, n_edge, n_cell
 
  do i=1, n_vertex
    read(10, *) arr_vertex(i, 1), arr_vertex(i, 2), arr_vertex(i, 3)
  end do
    
  do i=1, n_edge
    read(10, *) arr_edge(i, 1), arr_edge(i, 2)
  end do
    
  do i=1, n_cell
    read(10, *) arr_cell(i, 1), arr_cell(i, 2), arr_cell(i, 3)
  end do
  return   
  close(10)
end subroutine load_gts


subroutine save_paraview(fname,n_vertex,n_cell,arr_vertex,arr_cell)
  implicit none

  character*(*)               ::  fname
  integer                     ::  i
  integer :: n_cell,n_vertex
  real,dimension(:,:) :: arr_vertex(n_vertex,3)
  integer :: arr_cell(n_cell,3)

  open(unit=10, file=fname)
  write(10, '(A)')"# vtk DataFile Version 4.0"
  write(10,'(A)')"Cascadia Paraview"
  write(10,'(A)')"ASCII"
  write(10,'(A)')"DATASET POLYDATA"
  write(10,*)"POINTS",n_vertex,"float"


  do i=1, n_vertex
    write(10, *) arr_vertex(i, 1)/1000, arr_vertex(i, 2)/1000, 3*arr_vertex(i, 3)/1000
  end do

  write(10,*)"POLYGONS",n_cell, 4*n_cell
  do i=1, n_cell
    write(10,*) 3, arr_cell(i, 1)-1, arr_cell(i, 2)-1, arr_cell(i, 3)-1
  end do

  write(10,*)"CELL_DATA",n_cell
  write(10,*)"SCALARS cell_scalars float 1"
  write(10,*)"LOOKUP_TABLE default"
  do i=1,n_cell
   write(10,*) -arr_vertex(arr_cell(i,1),3)/1000
  end do

  close(10)
end subroutine save_paraview


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Main Program  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program paraview_mesh
  
  implicit none
  character(*),parameter :: fname="400km_1km_smooth.gts"   
  integer  :: i,j, n_vertex,n_cell,n_edge
  real,dimension(:,:),allocatable :: arr_vertex
  integer,dimension(:,:),allocatable :: arr_edge,arr_cell
  ! save to gnuplot data
  
  open(10,file=fname)
  read(10,*) n_vertex,n_edge,n_cell
  close(10) 
 allocate( arr_vertex(n_vertex,3),arr_edge(n_edge,2),arr_cell(n_cell,3))

  ! save 3-D mesh
  call load_gts(fname,n_vertex,n_edge,n_cell,arr_vertex,arr_edge,arr_cell)
  call save_paraview("slab_2_paraview_50km.vtk",n_vertex,n_cell,arr_vertex,arr_cell)
  
end program paraview_mesh 

