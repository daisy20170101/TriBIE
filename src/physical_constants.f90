!===============================================================================
! MODULE: Physical Constants
! 
! This module defines all physical constants used in the earthquake simulation
! including material properties, mathematical constants, and conversion factors.
!
! Author: Refactored from original phy3d_module_non.f90
! Last Modified: Current refactoring
!===============================================================================

module physical_constants
  implicit none
  
  ! Make all constants public
  public
  
  ! Precision specification using modern Fortran
  integer, parameter :: DP = selected_real_kind(15, 307)
  
  ! Mathematical constants
  real(DP), parameter :: PI = 4.0_DP * atan(1.0_DP)
  real(DP), parameter :: DEG_TO_RAD = PI / 180.0_DP
  real(DP), parameter :: RAD_TO_DEG = 180.0_DP / PI
  
  ! Material properties (from original code)
  real(DP), parameter :: SHEAR_MODULUS = 0.32038_DP
  real(DP), parameter :: SHEAR_WAVE_SPEED = 10.92407e7_DP
  real(DP), parameter :: POISSON_RATIO = 0.25_DP
  
  ! Friction law parameters
  real(DP), parameter :: REFERENCE_VELOCITY = 3.15e4_DP
  real(DP), parameter :: FRICTION_COEFFICIENT = 0.6_DP
  real(DP), parameter :: VISCOSITY = 0.5_DP * SHEAR_MODULUS / SHEAR_WAVE_SPEED
  real(DP), parameter :: GAMMA_FACTOR = 2.0_DP / PI
  
  ! Geometric parameters
  real(DP), parameter :: MAX_AMPLITUDE = 0.025_DP
  real(DP), parameter :: FAULT_ANGLE = 60.0_DP * DEG_TO_RAD
  
  ! Time conversion factors
  real(DP), parameter :: SECONDS_PER_YEAR = 365.0_DP * 24.0_DP * 3600.0_DP
  real(DP), parameter :: DAYS_PER_YEAR = 365.0_DP
  
  ! Small numbers for numerical stability
  real(DP), parameter :: EPSILON_SMALL = 1.0e-6_DP
  real(DP), parameter :: EPSILON_LARGE = 1.0e-12_DP
  
  ! Array bounds and limits
  integer, parameter :: MAX_ARRAY_SIZE = 1000000
  integer, parameter :: MAX_STRING_LENGTH = 256
  
end module physical_constants
