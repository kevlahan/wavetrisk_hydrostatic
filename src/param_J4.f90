module param_mod
  ! Jmin = DOMAIN_LEVEL + PATCH_LEVEL + 1
  integer, parameter :: DOMAIN_LEVEL = 2
  integer, parameter :: PATCH_LEVEL = 3 - DOMAIN_LEVEL
end module param_mod
