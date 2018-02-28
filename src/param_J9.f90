module param_mod
  ! Jmin = DOMAIN_LEVEL + PATCH_LEVEL + 1
  integer, parameter :: DOMAIN_LEVEL = 5
  integer, parameter :: PATCH_LEVEL = 8 - DOMAIN_LEVEL
end module param_mod
