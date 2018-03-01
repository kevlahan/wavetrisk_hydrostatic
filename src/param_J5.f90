module param_mod
  ! Jmin = DOMAIN_LEVEL + PATCH_LEVEL + 1
  integer, parameter :: DOMAIN_LEVEL = 1
  integer, parameter :: PATCH_LEVEL = 4 - DOMAIN_LEVEL
end module param_mod
