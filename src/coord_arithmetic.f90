module coord_arithmetic_mod
  ! Defines arithmetic operations for coord type variables
  ! includes vector addition/subtraction, scalar multiplication/division
  use shared_mod
  implicit none
  private

  public operator (+)
  public operator (-)
  public operator (*)
  public operator (/)

  interface operator (+)
     module procedure coord_plus_cc
     module procedure coord_plus_rc
     module procedure coord_plus_cr
  end interface operator (+)

  interface operator (-)
     module procedure coord_minus_cc
     module procedure coord_minus_rc
     module procedure coord_minus_cr
  end interface operator (-)

  interface operator (*)
     module procedure coord_times_cc
     module procedure coord_times_rc
     module procedure coord_times_cr
  end interface operator (*)

  interface operator (/)
     module procedure coord_divide_cr
  end interface operator (/)
contains
  pure type(coord) function coord_plus_cc (c1, c2)
    type(coord), intent (in) :: c1, c2

    coord_plus_cc = coord (c1%x+c2%x, c1%y+c2%y, c1%z+c2%z)
  end function coord_plus_cc

  pure type(coord) function coord_plus_rc (c1, c2)
    real(8),     intent (in) :: c1
    type(coord), intent (in) :: c2

    coord_plus_rc = coord (c1+c2%x, c1+c2%y, c1+c2%z)
  end function coord_plus_rc

  pure type(coord) function coord_plus_cr (c1, c2)
    real(8),     intent (in) :: c2
    type(coord), intent (in) :: c1

    coord_plus_cr = coord (c2+c1%x, c2+c1%y, c2+c1%z)
  end function coord_plus_cr

  pure type(coord) function coord_minus_cc (c1, c2)
    type(coord), intent (in) :: c1, c2

    coord_minus_cc = coord (c1%x-c2%x, c1%y-c2%y, c1%z-c2%z)
  end function coord_minus_cc

  pure type(coord) function coord_minus_rc (c1, c2)
    real(8),     intent (in) :: c1
    type(coord), intent (in) :: c2

    coord_minus_rc = coord (c1-c2%x, c1-c2%y, c1-c2%z)
  end function coord_minus_rc

  pure type(coord) function coord_minus_cr (c1, c2)
    real(8),     intent (in) :: c2
    type(coord), intent (in) :: c1

    coord_minus_cr = coord (c1%x-c2, c1%y-c2, c1%z-c2)
  end function coord_minus_cr

  pure real(8) function coord_times_cc (c1, c2)
    ! Inner product
    implicit none
    type(coord), intent (in) :: c1, c2

    coord_times_cc = c1%x*c2%x + c1%y*c2%y + c1%z*c2%z
  end function coord_times_cc

  pure type(coord) function coord_times_rc (c1, c2)
    implicit none
    real(8),     intent (in) :: c1
    type(coord), intent (in) :: c2

    coord_times_rc = coord (c1*c2%x, c1*c2%y, c1*c2%z)
  end function coord_times_rc

  pure type(coord) function coord_times_cr (c1, c2)
    implicit none
    real(8),     intent (in) :: c2
    type(coord), intent (in) :: c1

    coord_times_cr = coord (c2*c1%x, c2*c1%y, c2*c1%z)
  end function coord_times_cr

  pure type(coord) function coord_divide_cr (c1, c2)
    implicit none
    real(8),     intent (in) :: c2
    type(coord), intent (in) :: c1

    coord_divide_cr = coord (c1%x/c2, c1%y/c2, c1%z/c2)
  end function coord_divide_cr
end module coord_arithmetic_mod
