module param_mod
  ! Maximum number of cores <= number of domains = 10 * 2^(2*DOMAIN_LEVEL)
  ! PATCH_LEVEL = MIN_LEVEL - (DOMAIN_LEVEL + 1) >= 2
  ! Total nodes = 10 * 2^(2*MIN_LEVEL) 
  integer, parameter :: DOMAIN_LEVEL = 2 ! determines maximum cores and PATCH_LEVEL
  integer, parameter :: MIN_LEVEL    = 7 ! coarsest grid
end module param_mod
