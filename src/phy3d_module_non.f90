! Module to define global variables used in 3d_sub.f90 (or 3d_strike.f90)
Module phy3d_module_non
public

integer, parameter :: DP0=kind(1.d0)
integer :: IDin, IDout,Iprofile,Nd,Nl,Nd_all,Lratio,Nab,nprocs
integer :: nmv,nas,ncos,nnul,nsse
real (DP0), parameter :: pi = 3.14159265358979323
real (DP0), parameter :: xmu= 0.3d0, cs= 9.6d+7,xnu = 0.25d0, &
                         V0=3.2d+4,f0=0.6,eta=0.5*xmu/cs,  &
                         gamma=2.0/pi, &
                         p18 = 2.d0*pi/360.d0 
real (DP0), parameter :: yrs=365.*24.*3600.d0, yrd=365.d0

real (DP0) :: tm1,tm2,tmday,tmelse,tmmidn,tmmult,Vpl

real (DP0), DIMENSION(:), ALLOCATABLE :: cca,ccb,seff,xLf
!real (DP0), DIMENSION(:,:,:), ALLOCATABLE :: fr
real (DP0), DIMENSION(:,:), ALLOCATABLE :: stiff,stiff_ns
character(len=80) :: jobname,foldername,restartname,stiffname,profile
end module phy3d_module_non
