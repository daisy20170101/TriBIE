! Fortran functions for pore pressure evolution with Heaviside injection
! Based on equations (26) and (27) from the document

module pore_pressure_functions
    implicit none
    private
    public :: compute_G, compute_dGdt, dirac_delta, heaviside_function
    
    ! Mathematical constants
    real(8), parameter :: PI = 3.141592653589793d0
    real(8), parameter :: SQRT_PI = 1.772453850905516d0
    
contains

    ! Function G(z,t,alpha) from equation (27)
    ! G(z,t,α) = √t [exp(-z²/4αt)/√π - |z|/√(4αt) * erfc(|z|/√(4αt))]
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
        exp_term = exp(-(z**2) / (4.0d0 * alpha * t)) / SQRT_PI
        
        ! Second term: |z|/√(4αt) * erfc(|z|/√(4αt))
        erfc_term = (abs_z / sqrt_4at) * erfc_function(erfc_arg)
        
        ! G(z,t,α) = √t * [exp_term - erfc_term]
        compute_G = sqrt(t) * (exp_term - erfc_term)
        
    end function compute_G

    ! Time derivative of G: ∂G/∂t
    ! ∂G/∂t = (1/(2√t)) * [exp(-z²/4αt)/√π - |z|*erfc(|z|/√(4αt))/√(4αt)]
    ! Note: The z² exponential terms cancel, giving a simple form
    real(8) function compute_dGdt(z, t, alpha)
        implicit none
        real(8), intent(in) :: z, t, alpha
        real(8) :: abs_z, sqrt_4at, exp_term, erfc_term, erfc_arg
        
        if (t <= 0.0d0) then
            compute_dGdt = 0.0d0
            return
        endif
        
        abs_z = abs(z)
        sqrt_4at = sqrt(4.0d0 * alpha * t)
        erfc_arg = abs_z / sqrt_4at
        
        ! exp(-z²/4αt)/√π
        exp_term = exp(-(z**2) / (4.0d0 * alpha * t)) / SQRT_PI
        
        ! |z|*erfc(|z|/√(4αt))/√(4αt)
        erfc_term = abs_z * erfc_function(erfc_arg) / sqrt_4at
        
        ! ∂G/∂t = (1/(2√t)) * [exp_term - erfc_term]
        compute_dGdt = (1.0d0 / (2.0d0 * sqrt(t))) * (exp_term - erfc_term)
        
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
        
        dirac_delta = (1.0d0 / (eps * SQRT_PI)) * exp(-(t**2) / (eps**2))
        
    end function dirac_delta

    ! Heaviside function
    real(8) function heaviside_function(t)
        implicit none
        real(8), intent(in) :: t
        
        if (t >= 0.0d0) then
            heaviside_function = 1.0d0
        else
            heaviside_function = 0.0d0
        endif
        
    end function heaviside_function

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

end module pore_pressure_functions

! Example usage program
program test_functions
    use pore_pressure_functions
    implicit none
    
    real(8) :: z, t, alpha, G_val, dGdt_val, delta_val
    real(8) :: q0, beta, phi, sqrt_alpha
    real(8) :: pressure, dpdt
    
    ! Example parameters
    z = 10.0d0          ! Distance from injection point (m)
    t = 3600.0d0        ! Time (s)
    alpha = 1.0d-6      ! Hydraulic diffusivity (m²/s)
    q0 = 0.001d0        ! Injection rate
    beta = 1.0d0        ! Coupling parameter
    phi = 0.2d0         ! Porosity
    sqrt_alpha = sqrt(alpha)
    
    ! Test G function
    G_val = compute_G(z, t, alpha)
    write(*,*) 'G(z,t,α) = ', G_val
    
    ! Test dG/dt
    dGdt_val = compute_dGdt(z, t, alpha)
    write(*,*) 'dG/dt = ', dGdt_val
    
    ! Test Dirac delta at t=0
    delta_val = dirac_delta(0.0d0)
    write(*,*) 'δ(0) ≈ ', delta_val
    
    ! Calculate pressure using equation (26) for single injection period
    pressure = (q0 / (beta * phi * sqrt_alpha)) * G_val * heaviside_function(t)
    write(*,*) 'Pressure = ', pressure
    
    ! Calculate time derivative including Dirac delta contribution
    dpdt = (q0 / (beta * phi * sqrt_alpha)) * &
           (dGdt_val * heaviside_function(t) + G_val * dirac_delta(t))
    write(*,*) 'dp/dt = ', dpdt
    
end program test_functions