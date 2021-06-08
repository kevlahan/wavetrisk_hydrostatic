module spline_mod
  ! Routines for cubic spline interpolation
  implicit none
contains
  subroutine splint (xa, ya, y2a, n, x, y)
    ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
    ! (with the xa(i) in order), and given the array y2a(1:n), which is the output
    ! from the subroutine spline, and given a value of x, this routine returns a
    ! cubic spline interpolated value y.
    !
    ! (adopted from Numerical Recipes in FORTRAN 77)
    !
    integer               :: n
    real(8)               :: x, y
    real(8), dimension(n) :: xa, y2a, ya
    
    integer :: k, khi, klo
    real(8) :: a, b, h

    klo = 1
    khi = n
1   if (khi - klo > 1) then
       k = (khi + klo) / 2
       if (xa(k) > x) then
          khi = k
       else
          klo = k
       end if
       goto 1
    end if

    h = xa(khi) - xa(klo)
    if (h == 0d0) write(6,*)"bad xa input in splint"

    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h
    y = a * ya(klo) + b * ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b) * y2a(khi)) * (h**2)/6d0

    return
  end subroutine splint

  subroutine spline (x, y, n, yp1, ypn, y2)
    ! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
    ! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
    ! the first derivative of the interpolating function at points 1 and n,
    ! respectively, this routine returns an array y2(1:n) of length n which
    ! contains the second derivatives of the interpolating function at the
    ! tabulated points x(i).  If yp1 and/or ypn are equal to 1d30 or larger,
    ! the routine is signaled to set the corresponding boundary condition for a
    ! natural spline with zero second derivative on that boundary.
    !
    ! (adopted from Numerical Recipes in FORTRAN 77)
    !
    integer               :: n
    real(8)               :: yp1, ypn
    real(8), dimension(n) :: x, y, y2
    
    integer                  :: i, k
    integer, parameter       :: nmax = 4000 ! largest anticipated value of n
    real(8)                  :: p, qn, sig, un
    real(8), dimension(nmax) :: u

    if (yp1 > 0.99d30) then
       y2(1) = 0d0
       u(1)  = 0d0
    else
       y2(1) = -0.5d0
       u(1) = (3d0 / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
    end if

    do i = 2, n-1
       sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
       p = sig * y2(i-1) + 2d0
       y2(i) = (sig - 1d0) / p
       u(i) = (6d0 * ((y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1)) / &
            (x(i)-x(i-1))) / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
    end do

    if (ypn > 0.99d30) then
       qn = 0d0
       un = 0d0
    else
       qn = 0.5d0
       un = (3d0 / (x(n) - x(n-1))) * (ypn - (y(n) - y(n-1)) / (x(n) - x(n-1)))
    end if

    y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1d0)

    do k = n-1, 1, -1
       y2(k) = y2(k) * y2(k+1) + u(k)
    end do

    return
  end subroutine spline
end module spline_mod
