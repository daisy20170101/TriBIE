! Module to define global variables used in 3d_sub.f90 (or 3d_strike.f90)
Module phy3d_module_bp6
public :: compute_G, compute_dGdt, dirac_delta, heavi

integer, parameter :: DP0=kind(1.d0)
integer :: IDin, IDout,Iprofile,Nd,Nl,Nd_all,Lratio,Nab,nprocs
integer :: nmv,nas,ncos,nnul,nsse
real (DP0), parameter :: pi = 3.14159265358979323, sqrt_pi = 1.772453850905516, amax=0.007
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
integer, DIMENSION(:), ALLOCATABLE :: sendcounts, displs
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

  real(8) function compute_G(z, t, alpha)
        implicit none
        real(8), intent(in) :: z, t, alpha
        real(8) :: abs_z, sqrt_4at, exp_term, erfc_term, erfc_arg
        
        if (t <= 0.0d0) then
            compute_G = 0.0d0
            return
        endif
        
        abs_z = abs(z)
        sqrt_4at = sqrt(4.0d0 * alpha * t)
        erfc_arg = abs_z / sqrt_4at
        
        ! First term: exp(-z²/4αt)/√π
        exp_term = exp(-(z**2) / (4.0d0 * alpha * t)) / sqrt_pi
        
        ! Second term: |z|/√(4αt) * erfc(|z|/√(4αt))
        erfc_term = (abs_z / sqrt_4at) * erfc_function(erfc_arg)
        
        ! G(z,t,α) = √t * [exp_term - erfc_term]
        compute_G = sqrt(t) * (exp_term - erfc_term)
        
    end function compute_G

    ! Time derivative of G: ∂G/∂t
    real(8) function compute_dGdt(z, t, alpha)
        implicit none
        real(8), intent(in) :: z, t, alpha
        real(8) :: abs_z, sqrt_4at, sqrt_t, z_squared
        real(8) :: exp_term, erfc_term, erfc_arg
        real(8) :: dexp_dt, derfc_dt, term1, term2, term3
        
        if (t <= 0.0d0) then
            compute_dGdt = 0.0d0
            return
        endif
        
        abs_z = abs(z)
        sqrt_t = sqrt(t)
        sqrt_4at = sqrt(4.0d0 * alpha * t)
        z_squared = z**2
        erfc_arg = abs_z / sqrt_4at
        
        exp_term = exp(-z_squared / (4.0d0 * alpha * t)) / sqrt_pi
        erfc_term = (abs_z / sqrt_4at) * erfc_function(erfc_arg)
        
        ! ∂G/∂t = (1/2√t) * [exp_term - erfc_term] + √t * [∂exp_term/∂t - ∂erfc_term/∂t]
        
        ! First part: (1/2√t) * [exp_term - erfc_term]
        term1 = 0.5d0 / sqrt_t * (exp_term - erfc_term)
        
        ! ∂exp_term/∂t = exp_term * z²/(4αt²)
        dexp_dt = exp_term * z_squared / (4.0d0 * alpha * t**2)
        
        ! ∂erfc_term/∂t is more complex
        derfc_dt = -abs_z / (2.0d0 * sqrt_4at * t) * erfc_function(erfc_arg) + &
                   abs_z * z_squared / (4.0d0 * sqrt_pi * alpha * t**2 * sqrt_4at) * &
                   exp(-erfc_arg**2)
        
        term2 = sqrt_t * dexp_dt
        term3 = -sqrt_t * derfc_dt
        
        compute_dGdt = term1 + term2 + term3
        
    end function compute_dGdt

    ! Dirac delta function approximation
    ! Uses a narrow Gaussian approximation: δ(t) ≈ (1/ε√π) * exp(-t²/ε²)
    real(8) function dirac_delta(t, epsilon)
        implicit none
        real(8), intent(in) :: t
        real(8), intent(in), optional :: epsilon
        real(8) :: eps
        
        ! Default epsilon for numerical approximation
        if (present(epsilon)) then
            eps = epsilon
        else
            eps = 1.0d-6  ! Small value for sharp approximation
        endif
        
        dirac_delta = (1.0d0 / (eps * sqrt_pi)) * exp(-(t**2) / (eps**2))
        
    end function dirac_delta


    ! Complementary error function approximation
    ! Using rational approximation (Abramowitz and Stegun)
    real(8) function erfc_function(x)
        implicit none
        real(8), intent(in) :: x
        real(8) :: t, tau, result
        real(8), parameter :: a1 = 0.254829592d0
        real(8), parameter :: a2 = -0.284496736d0  
        real(8), parameter :: a3 = 1.421413741d0
        real(8), parameter :: a4 = -1.453152027d0
        real(8), parameter :: a5 = 1.061405429d0
        real(8), parameter :: p = 0.3275911d0
        
        t = 1.0d0 / (1.0d0 + p * abs(x))
        tau = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))))
        result = tau * exp(-x**2)
        
        if (x >= 0.0d0) then
            erfc_function = result
        else
            erfc_function = 2.0d0 - result
        endif
        
    end function erfc_function

end module phy3d_module_bp6

