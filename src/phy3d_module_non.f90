! Module to define global variables used in 3d_sub.f90 (or 3d_strike.f90)
Module phy3d_module_non
public

integer, parameter :: DP0=kind(1.d0)
integer :: IDin, IDout,Iprofile,Nd,Nl,Nd_all,Lratio,Nab,nprocs
integer :: nmv,nas,ncos,nnul,nsse
real (8), parameter :: pi = 3.14159265358979323, amax=0.025,phi=60.0/180.0*pi
real (8), parameter :: xmu= 0.32038d0, cs=10.92407d+7, xnu = 0.25d0, &
                         V0=3.15d+4,f0=0.6,eta=0.5*xmu/cs,  &
                         gamma=2.0/pi, &
                         p18 = 2.d0*pi/360.d0 
real (8), parameter :: yrs=365.*24.*3600.d0, yrd=365.d0

real (8) :: tm1,tm2,tmday,tmelse,tmmidn,tmmult,Vpl

real (8), DIMENSION(:), ALLOCATABLE :: tau1,tau2,tau0,cca,ccb,seff,xLf,vi,phy1,phy2
!real (DP0), DIMENSION(:,:,:), ALLOCATABLE :: fr
real (8), DIMENSION(:,:), ALLOCATABLE :: stiff,stiff2
character(len=80) :: jobname,foldername,restartname,stiffname,profile
end module phy3d_module_non
