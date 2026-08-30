! Module to define global variables used in 3d_sub.f90 (or 3d_strike.f90)
Module phy3d_module_bp6
public :: compute_G, compute_dGdt, dirac_delta, heavi, compute_pf, compute_dpf_dt

integer, parameter :: DP0=kind(1.d0)
integer :: IDin, IDout,Iprofile,Nd,Nl,Nd_all,Lratio,Nab,nprocs
integer :: nmv,nas,ncos,nnul,nsse
real (DP0), parameter :: pi = 3.14159265358979323d0, sqrt_pi = 1.772453850905516d0, amax=0.007
real (DP0), parameter :: xmu= 32.038d9, cs=3464, xnu = 0.25d0, &
                         V0=1d-6,f0=0.6,eta=0.5*xmu/cs,  &
                         gamma=2.0/pi, &
                         p18 = 2.d0*pi/360.d0 
real (DP0), parameter :: yrs=365.*24.*3600.d0, yrd=365.d0
! unit: sec, pa-1, m/s,  alpha, 0.1 m^2/s
real (DP0), parameter :: toff =  100.0d0*24.0d0*3600.0d0, beta = 1.d-8, q0 = 1.25d-6,phi = 0.1d0, alpha= 0.1d0 ! unit converse
real (DP0), parameter :: kappa=1.d-13, eta_diff= 1.d-3
real (DP0), parameter :: tauini =29.20d6,tp=100.0,reb=1d-6
! Floor on the evolving effective normal stress (Pa). Fault opening is not
! modelled, so sigma must stay positive or the regularised friction law breaks
! (log/sqrt of a non-positive argument). Applied to a LOCAL copy inside derivs
! and in the tau1 diagnostic; the integrated state is never rewritten.
real (DP0), parameter :: sigma_min = 1.0d5

real (DP0) ::tsec, tm1,tm2,tmday,tmelse,tmmidn,tmmult,Vpl
real (DP0) ::dipangle
integer, DIMENSION(:), ALLOCATABLE :: sendcounts, displs
real (DP0), DIMENSION(:), ALLOCATABLE :: dvel,pp1,tau1,tau2,tau0,cca,ccb,seff,xLf,phy1,phy2
! Per-element plate loading rate (m/s), optional 6th column of the var file.
! Falls back to the scalar Vpl from parameter1.txt for 5-column var files,
! so existing inputs (example1/2/3) are unaffected. See resdep in 3dtri_BP5.f90.
real (DP0), DIMENSION(:), ALLOCATABLE :: vplv
! Per-element initial shear stress (Pa), optional 7th column of the var file.
! Falls back to the scalar tauini above for <7-column files, so example1/2/3
! are unaffected. Set to the steady-state stress at the local loading rate by
! prepare_input.py, which removes the artificial initial loading transient.
real (DP0), DIMENSION(:), ALLOCATABLE :: tauv
!real (DP0), DIMENSION(:,:,:), ALLOCATABLE :: fr
real (DP0), DIMENSION(:,:), ALLOCATABLE :: stiff,stiff2
character(len=80) :: jobname,foldername,restartname,stiffname,profile
real(DP0),parameter :: vini = 1.0d-12 ! initial vel
real(DP0), parameter :: y2=-50.0,y3=-50.0, T0=1.0d0, tau_p0=17.5d0, rr_nuc=150.0d0

contains

real(DP0)  function heavi(x) !Heaviside function, useful in DSP 
 implicit none
 integer, parameter :: DP0=kind(1.d0)
 real(DP0) :: x
  heavi = 0.5d0*(sign(1.d0,x)+1.d0)
 end function heavi

  real(DP0) function compute_G(z, t, alpha)
        implicit none
        real(DP0), intent(in) :: z, t, alpha
        real(DP0) :: abs_z, sqrt_4at, exp_term, erfc_term, erfc_arg
 
        if (t <= 0.0d0) then
            compute_G = 0.0d0
            return
        endif
        
        abs_z = dabs(z)
        sqrt_4at = dsqrt(4.0d0 * alpha * t)
        erfc_arg = abs_z / sqrt_4at
        
        ! First term: exp(-z²/4αt)/√π
        exp_term = dexp(-(z**2) / (4.0d0 * alpha * t)) / sqrt_pi
        
        ! Second term: |z|/√(4αt) * erfc(|z|/√(4αt))
        erfc_term = (abs_z / sqrt_4at) * erfc(erfc_arg)
        
        ! G(z,t,α) = √t * [exp_term - erfc_term]
        compute_G = dsqrt(t) * (exp_term - erfc_term)
        
    end function compute_G

    ! Time derivative of G: ∂G/∂t
 ! ∂G/∂t = (1/(2√t)) * [exp(-z²/4αt)/√π - |z|*erfc(|z|/√(4αt))/√(4αt)]
    ! Note: The z² exponential terms cancel, giving a simple form
    
 real(DP0) function compute_dGdt(z, t, alpha)
        implicit none
        real(DP0), intent(in) :: z, t, alpha
        real(DP0) :: abs_z, sqrt_4at, exp_term, erfc_term, erfc_arg
 
        if (t <= 0.0d0) then
            compute_dGdt = 0.0d0
            return
        endif
        
        abs_z = dabs(z)
        sqrt_4at = dsqrt(4.0d0 * alpha * t)
        erfc_arg = abs_z / sqrt_4at
        
        ! exp(-z²/4αt)/√π
        exp_term = dexp(-(z**2) / (4.0d0 * alpha * t)) / SQRT_PI
        
        ! |z|*erfc(|z|/√(4αt))/√(4αt)
        erfc_term = abs_z * erfc(erfc_arg) / sqrt_4at
        
        ! ∂G/∂t = (1/(2√t)) * [exp_term - erfc_term]
        compute_dGdt = (1.0d0 / (2.0d0 * dsqrt(t))) * (exp_term - erfc_term)
        
    end function compute_dGdt

    ! Dirac delta function approximation
    ! Uses a narrow Gaussian approximation: δ(t) ≈ (1/ε√π) * exp(-t²/ε²)
    real(DP0) function dirac_delta(t, epsilon)
        implicit none
        real(DP0), intent(in) :: t
        real(DP0), intent(in), optional :: epsilon
        real(DP0) :: eps
        
        ! Default epsilon for numerical approximation
        if (present(epsilon)) then
            eps = epsilon
        else
            eps = 1.0d-12  ! Small value for sharp approximation
        endif
        
        dirac_delta = (1.0d0 / (eps * sqrt_pi)) * dexp(-(t/eps)**2)
        
    end function dirac_delta


    ! Complementary error function approximation
    ! Using rational approximation (Abramowitz and Stegun)
    real(DP0) function erfc_function_rat(x)
        implicit none
        real(DP0), intent(in) :: x
        real(DP0) :: t, tau, result
        real(DP0), parameter :: a1 = 0.254829592d0
        real(DP0), parameter :: a2 = -0.284496736d0  
        real(DP0), parameter :: a3 = 1.421413741d0
        real(DP0), parameter :: a4 = -1.453152027d0
        real(DP0), parameter :: a5 = 1.061405429d0
        real(DP0), parameter :: p = 0.3275911d0
        
        t = 1.0d0 / (1.0d0 + p * abs(x))
        tau = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))))
        result = tau * dexp(-x**2)
        
        if (x >= 0.0d0) then
            erfc_function_rat = result
        else
            erfc_function_rat = 2.0d0 - result
        endif
        
    end function erfc_function_rat

real(DP0) function erfc_function(x)
    implicit none
    real(DP0), intent(in) :: x
    real(DP0) :: erfcresult,a, b, c, d, h, del
    real(DP0), parameter :: EPS = 1.0d-15
    real(DP0), parameter :: FPMIN = 1.0d-30
    real(DP0), parameter :: SQRT_PI = 1.772453850905516d0
    integer :: i
    
    
    ! Lentz's method for continued fraction
    b = x*x + 0.5d0
    c = 1.0d0 / FPMIN
    d = 1.0d0 / b
    h = d
    
    do i = 1, 100
        a = -i * (i - 0.5d0)
        b = b + 2.0d0
        d = a * d + b
        if (dabs(d) < FPMIN) d = FPMIN
        c = b + a / c
        if (dabs(c) < FPMIN) c = FPMIN
        d = 1.0d0 / d
        del = d * c
        h = h * del
        if (dabs(del - 1.0d0) < EPS) exit
    enddo
    erfcresult = dexp(-x*x) * h / SQRT_PI   

    if (x>=0.d0) then

      erfc_function = erfcresult
    else
      erfc_function = 2.d0 - erfcresult

    end if
end function erfc_function

    ! Compute pore fluid pressure
    ! p_fluid = (q0/(βφ√α)) * [G(t) - G(t-toff)]
    
    real(DP0) function compute_pf(z, t, alpha, beta, phi, q0, toff)
        implicit none
        real(DP0), intent(in) :: z, t, alpha, beta, phi, q0, toff
        real(DP0) :: G_val, G_val_off
        
        ! Compute Green's functions
        G_val = compute_G(z, t, alpha)
        G_val_off = compute_G(z, t - toff, alpha)
        
        ! Compute pore fluid pressure
        compute_pf = q0 / (beta * phi * dsqrt(alpha)) * (G_val * heavi(t) - G_val_off * heavi(t - toff))
        
    end function compute_pf

    ! Compute time derivative of pore fluid pressure
    ! dp_fluid/dt = (q0/(βφ√α)) * [dG/dt * H(t) + G(t) * δ(t) - dG/dt * H(t-toff) - G(t-toff) * δ(t-toff)]
    real(DP0) function compute_dpf_dt(z, t, alpha, beta, phi, q0, toff)
        implicit none
        real(DP0), intent(in) :: z, t, alpha, beta, phi, q0, toff
        real(DP0) :: delta_t,dGdt_val_off,G_val, G_val_off, dGdt_val,term1, term2,term3,term4

        
        ! Compute Green's functions and their time derivatives
        !if(t>0.d0)then
        !
        !  G_val = compute_G(z, t, alpha)
        !  dGdt_val = compute_dGdt(z, t, alpha)
        !  term1 = dGdt_val * heavi(t)
        !  term2 = G_val * dirac_delta(t)
        !else
        !  term1 = 0.d0
        !  term2 = 0.d0
        !end if
  
        !if(t.gt.toff)then
        !   dGdt_val_off = compute_dGdt(z,t-toff,alpha)
        !   G_val_off = compute_G(z,t-toff,alpha)
        !   term3 = dGdt_val_off * heavi(t-toff)
        !   term4 = G_val_off *dirac_delta(t-toff)
        !else
        !   term3=0.d0
        !   term4=0.d0
        !end if
        ! Compute time derivative of pore fluid pressure
        !compute_dpf_dt = q0 / (beta * phi * dsqrt(alpha)) * (term1 + term2 - term3 - term4) 

        delta_t = 0.01d0
        compute_dpf_dt = 1/delta_t *(compute_pf(z,t+delta_t,alpha,beta,phi,q0,toff) - compute_pf(z,t,alpha,beta,phi,q0,toff))
    end function compute_dpf_dt

end module phy3d_module_bp6

