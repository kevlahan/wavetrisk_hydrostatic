module param_mod
  ! Maximum number of cores = 10 * 2^(2*DOMAIN_LEVEL)
  integer, parameter :: DOMAIN_LEVEL = 2 ! <= 4 (MIN_LEVEL - PATCH_LEVEL - 1 with PATCH_LEVEL >= 2)
  integer, parameter :: MIN_LEVEL    = 7 ! sets coarsest grid
end module param_mod
