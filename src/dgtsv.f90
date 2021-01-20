subroutine dgtsv (n, nrhs, dl, d, du, b, ldb, info)
  ! DGTSV  solves the equation
  !
  !  A*X = B,
  ! 
  ! where A is an n by n tridiagonal matrix, by Gaussian elimination with
  ! partial pivoting.
  !
  ! Note that the equation  A**T*X = B  may be solved by interchanging the
  ! order of the arguments DU and DL.
  !
  integer          :: i, j, info, ldb, n, nrhs
  double precision ::  fact, temp
  double precision :: b( ldb, * ), d( * ), dl( * ), du( * )
  double precision, parameter :: zero = 0.0d+0
  intrinsic abs, max
  external  xerbla

  info = 0
  if ( n < 0 ) then
     info = -1
  else if ( nrhs < 0 ) then
     info = -2
  else if ( ldb < max( 1, n ) ) then
     info = -7
  end if
  if ( info /= 0 ) then
     call xerbla ( 'dgtsv ', -info )
     return
  end if

  if( n == 0 ) return

  if( nrhs == 1 ) then
     do i = 1, n - 2
        if ( abs( d( i ) ) >= abs( dl( i ) ) ) then ! no row interchange required
           if( d( i ) /= zero ) then
              fact = dl( i ) / d( i )
              d( i+1 ) = d( i+1 ) - fact*du( i )
              b( i+1, 1 ) = b( i+1, 1 ) - fact*b( i, 1 )
           else
              info = i
              return
           end if
           dl( i ) = zero
        else ! interchange rows i and i+1
           fact = d( i ) / dl( i )
           d( i ) = dl( i )
           temp = d( i+1 )
           d( i+1 ) = du( i ) - fact*temp
           dl( i ) = du( i+1 )
           du( i+1 ) = -fact*dl( i )
           du( i ) = temp
           temp = b( i, 1 )
           b( i, 1 ) = b( i+1, 1 )
           b( i+1, 1 ) = temp - fact*b( i+1, 1 )
        end if
     end do
     
     if (n > 1) then
        i = n - 1
        
        if (abs( d( i ) ) >= abs( dl( i ) )) then
           if( d( i ) /= zero ) then
              fact = dl( i ) / d( i )
              d( i+1 ) = d( i+1 ) - fact*du( i )
              b( i+1, 1 ) = b( i+1, 1 ) - fact*b( i, 1 )
           else
              info = i
              return
           end if
        else
           fact = d( i ) / dl( i )
           d( i ) = dl( i )
           temp = d( i+1 )
           d( i+1 ) = du( i ) - fact*temp
           du( i ) = temp
           temp = b( i, 1 )
           b( i, 1 ) = b( i+1, 1 )
           b( i+1, 1 ) = temp - fact*b( i+1, 1 )
        end if
     end if

     if ( d( n ) == zero ) then
        info = n
        return
     end if
  else
     do i = 1, n - 2
        if( abs( d( i ) ) >= abs( dl( i ) ) ) then ! no row interchange required
           if ( d( i ) /= zero ) then
              fact = dl( i ) / d( i )
              d (i+1) = d(i+1) - fact*du( i )
              do j = 1, nrhs
                 b( i+1, j ) = b( i+1, j ) - fact*b( i, j )
              end do
           else
              info = i
              return
           end if
           dl( i ) = zero
        else ! interchange rows i and i+1
           fact = d( i ) / dl( i )
           d( i ) = dl( i )
           temp = d( i+1 )
           d( i+1 ) = du( i ) - fact*temp
           dl( i ) = du( i+1 )
           du( i+1 ) = -fact*dl( i )
           du( i ) = temp
           do j = 1, nrhs
              temp = b( i, j )
              b( i, j ) = b( i+1, j )
              b( i+1, j ) = temp - fact*b( i+1, j )
           end do
        end if
     end do

     if ( n > 1 ) then
        i = n - 1
        if ( abs( d( i ) ) >= abs( dl( i ) ) ) then
           if( d( i ) /= zero ) then
              fact = dl( i ) / d( i )
              d( i+1 ) = d( i+1 ) - fact*du( i )
              do j = 1, nrhs
                 b( i+1, j ) = b( i+1, j ) - fact*b( i, j )
              end do
           else
              info = i
              return
           end if
        else
           fact = d( i ) / dl( i )
           d( i ) = dl( i )
           temp = d( i+1 )
           d( i+1 ) = du( i ) - fact*temp
           du( i ) = temp
           do j = 1, nrhs
              temp = b( i, j )
              b( i, j ) = b( i+1, j )
              b( i+1, j ) = temp - fact*b( i+1, j )
           end do
        end if
     end if
     if ( d( n ) == zero ) then
        info = n
        return
     end if
  end if

  !   back solve with the matrix u from the factorization.
  if ( nrhs <= 2 ) then
     j = 1
70   continue
     b( n, j ) = b( n, j ) / d( n )
     if ( n > 1 ) b( n-1, j ) = ( b( n-1, j ) - du( n-1 )*b( n, j ) ) / d( n-1 )
     do i = n - 2, 1, -1
        b( i, j ) = ( b( i, j ) - du( i )*b( i+1, j ) - dl( i ) * b( i+2, j ) ) / d( i )
     end do
     if ( j < nrhs ) then
        j = j + 1
        go to 70
     end if
  else
     do j = 1, nrhs
        b( n, j ) = b( n, j ) / d( n )
        if ( n > 1 ) b( n-1, j ) = ( b( n-1, j ) - du( n-1 )*b( n, j ) ) / d( n-1 )
        do i = n - 2, 1, -1
           b( i, j ) = ( b( i, j ) - du( i )*b( i+1, j )-dl( i )* b( i+2, j ) ) / d( i )
        end do
     end do
  end if
  return
end subroutine dgtsv

