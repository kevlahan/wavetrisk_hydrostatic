module kind_mod
  ! Sets variable kinds
  implicit none
  integer, parameter :: sp = selected_real_kind(p=6,  r=37 )             ! single precision
  integer, parameter :: dp = selected_real_kind(p=15, r=307)             ! double precision
end module kind_mod
