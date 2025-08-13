! Module to define global variables used in 3d_sub.f90 (or 3d_strike.f90)
Module phy3d_module_bp6
public

integer, parameter :: DP0=kind(1.d0)
integer :: IDin, IDout,Iprofile,Nd,Nl,Nd_all,Lratio,Nab,nprocs
integer :: nmv,nas,ncos,nnul,nsse
real (DP0), parameter :: pi = 3.14159265358979323, amax=0.007
real (DP0), parameter :: xmu= 32.038d9, cs=3464, xnu = 0.25d0, &
                         V0=1d-6,f0=0.6,eta=0.5*xmu/cs,  &
                         gamma=2.0/pi, &
                         p18 = 2.d0*pi/360.d0 
real (DP0), parameter :: yrs=365.*24.*3600.d0, yrd=365.d0
! unit: sec, pa-1, m/s,  alpha, 0.1 m^2/s
real (DP0), parameter :: toff =  100.0*24*3600.0, beta = 1d-8, q0 = 1.25d-6,phi = 0.1, alpha= 0.1 ! unit converse
real (DP0), parameter :: kappa=1d-13, eta_diff= 1d-3
real (DP0), parameter :: tauini =29.20d6,tp=100.0,reb=1d-6

real (DP0) ::tsec, tm1,tm2,tmday,tmelse,tmmidn,tmmult,Vpl
real (DP0) ::dipangle

real (DP0), DIMENSION(:), ALLOCATABLE :: dvel,pp1,tau1,tau2,tau0,cca,ccb,seff,xLf,phy1,phy2
!real (DP0), DIMENSION(:,:,:), ALLOCATABLE :: fr
real (DP0), DIMENSION(:,:), ALLOCATABLE :: stiff,stiff2
character(len=80) :: jobname,foldername,restartname,stiffname,profile
real(DP0),parameter :: vini = 31.5 ! initial vel
real(DP0), parameter :: y2=-50.0,y3=-50.0, T0=1.0d0, tau_p0=17.5d0, rr_nuc=150.0d0

contains

 function heavi(x) !Heaviside function, useful in DSP 
 implicit none
 integer, parameter :: DP0=kind(1.d0)
 real(DP0) :: heavi
 real(DP0) :: x
  heavi = 0.5*(sign(1.d0,x)+1.0)
 end function heavi

end module phy3d_module_bp6

