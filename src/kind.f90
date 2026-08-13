module kind_mod
  ! Sets variable kinds                                                                                       
  use, intrinsic :: iso_fortran_env
  
  implicit none
  
  integer, parameter :: sp = REAL32 ! single precision
  integer, parameter :: dp = REAL64 ! double precision
end module kind_mod
