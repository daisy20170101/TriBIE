module precision_module
  implicit none
  
  ! Define precision parameter for better portability
  integer, parameter :: DP = kind(1.0d0)
  
end module precision_module
