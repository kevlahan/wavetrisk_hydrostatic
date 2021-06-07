module projection_mod
  ! Routines to implement projection onto 2D lon-lat grid
  use comm_mpi_mod
  implicit none
  integer, dimension(2)                :: Nx, Ny
  real(8), dimension(2)                :: lon_lat_range
  real(8), dimension(:),   allocatable :: lat, lon
  real(8), dimension(:,:), allocatable :: xcoord_lat, xcoord_lon
  real(8)                              :: dx_export, dy_export, kx_export, ky_export
  real(4), dimension(:,:), allocatable :: field2d
  real(8), dimension(:),   pointer     :: proj_sclr
contains
  subroutine initialize_projection (m)
    ! Initialize 2d projection variables on lon-lat grid of size (-m/2, m/2) x (-m/4, m/4)
    implicit none
    integer :: m
    
    integer :: i

    Nx = (/-m/2, m/2/)
    Ny = (/-m/4, m/4/)

    lon_lat_range = (/2d0*MATH_PI, MATH_PI/)
    dx_export = lon_lat_range(1) / (Nx(2) - Nx(1) + 1)
    dy_export = lon_lat_range(2) / (Ny(2) - Ny(1) + 1)
    
    kx_export = 1d0 / dx_export
    ky_export = 1d0 / dy_export

    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))
    allocate (lat(Ny(1):Ny(2)), lon(Nx(1):Nx(2)))
    allocate (xcoord_lat(Ny(1):Ny(2),1:2), xcoord_lon(Nx(1):Nx(2),1:2))

    do i = Nx(1), Nx(2)
       lon(i) = -180d0 + dx_export * (i - Nx(1))/MATH_PI * 180d0
    end do

    do i = Ny(1), Ny(2)
       lat(i) = -90d0 + dy_export * (i - Ny(1))/MATH_PI * 180d0
    end do
  end subroutine initialize_projection
  
  subroutine project_field_onto_plane (field, l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    implicit none
    integer           :: l, itype
    real(8)           :: default_val
    Type(Float_field) :: field

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = field%data(d)%elts(id+1)
                valN  = field%data(d)%elts(idN+1)
                valE  = field%data(d)%elts(idE+1)
                valNE = field%data(d)%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_field_onto_plane

  subroutine project_array_onto_plane (array, l, default_val)
    ! Projects array from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    implicit none
    integer      :: l, itype
    real(8)      :: default_val
    character(*) :: array

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       ! Determine Float_Array to project
       select case (array)
       case ("coriolis")
          proj_sclr => grid(d)%coriolis%elts
       case ("surf_press")
          proj_sclr => grid(d)%surf_press%elts
       case ("geopot")
          proj_sclr => grid(d)%geopot%elts
       case ("u_zonal")
          proj_sclr => grid(d)%u_zonal%elts
       case ("v_merid")
          proj_sclr => grid(d)%v_merid%elts
       case ("press_lower")
          proj_sclr => grid(d)%press_lower%elts
       case ("geopot_lower")
          proj_sclr => grid(d)%geopot_lower%elts
       case ("vort")
          proj_sclr => grid(d)%vort%elts
       case ("ke")
          proj_sclr => grid(d)%ke%elts
       case ("bernoulli")
          proj_sclr => grid(d)%bernoulli%elts
       case ("divu")
          proj_sclr => grid(d)%divu%elts
       case ("qe")
          proj_sclr => grid(d)%qe%elts
       case ("topo")
          proj_sclr => grid(d)%topo%elts
       end select
       
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = proj_sclr(id+1)
                valN  = proj_sclr(idN+1)
                valE  = proj_sclr(idE+1)
                valNE = proj_sclr(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
       nullify (proj_sclr)
    end do
    
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_array_onto_plane
  
  subroutine interp_tri_to_2d (a, b, c, val)
    implicit none
    real(8), dimension(2) :: a, b, c
    real(8), dimension(3) :: val

    integer               :: id_x, id_y
    real(8)               :: ival, minx, maxx, miny, maxy
    real(8), dimension(2) :: ll
    real(8), dimension(3) :: bac
    logical               :: inside

    minx = min (min (a(1), b(1)), c(1))
    maxx = max (max (a(1), b(1)), c(1))
    miny = min (min (a(2), b(2)), c(2))
    maxy = max (max (a(2), b(2)), c(2))
    if (maxx-minx > MATH_PI/2d0) then
       write (0,'(A,i4,A)') 'ERROR (rank = ', rank, '): io-333 "export"'
       return
    end if

    do id_x = floor (kx_export*minx), ceiling (kx_export*maxx)
       if (id_x < lbound (field2d,1) .or. id_x > ubound (field2d,1)) cycle
       do id_y = floor (ky_export*miny), ceiling (ky_export*maxy)
          if (id_y < lbound (field2d,2) .or. id_y > ubound (field2d,2)) cycle
          ll = (/ dx_export*id_x, dy_export*id_y /)
          call interp_tria (ll, a, b, c, val, ival, inside)
          if (inside) field2d(id_x,id_y) = ival
       end do
    end do
  end subroutine interp_tri_to_2d

  subroutine interp_tri_to_2d_and_fix_bdry (a0, b0, c0, val)
    implicit none
    real(8), dimension(2) :: a0, b0, c0
    real(8), dimension(3) :: val

    integer               :: i
    integer, dimension(3) :: fixed
    real(8), dimension(2) :: a, b, c

    a = a0
    b = b0
    c = c0
    call fix_boundary (a(1), b(1), c(1), fixed(1))
    call fix_boundary (b(1), c(1), a(1), fixed(2))
    call fix_boundary (c(1), a(1), b(1), fixed(3))
    call interp_tri_to_2d (a, b, c, val)

    if (sum(abs(fixed)) > 1) write (0,'(A)') 'ALARM'

    if (sum(fixed) /= 0) then
       a(1) = a(1) - sum(fixed) * 2d0*MATH_PI
       b(1) = b(1) - sum(fixed) * 2d0*MATH_PI
       c(1) = c(1) - sum(fixed) * 2d0*MATH_PI
       call interp_tri_to_2d (a, b, c, val)
    end if
  end subroutine interp_tri_to_2d_and_fix_bdry

  subroutine cart2sph2 (cin, cout)
    implicit none
    type(Coord)                        :: cin
    real(8), dimension(2), intent(out) :: cout

    call cart2sph (cin, cout(1), cout(2))
  end subroutine cart2sph2

  subroutine interp_tria (ll, coord1, coord2, coord3, values, ival, inside)
    implicit none
    real(8), dimension(2) :: coord1, coord2, coord3
    real(8), dimension(3) :: values
    real(8)               :: ival
    logical               :: inside

    real(8), dimension(2) :: ll
    real(8), dimension(3) :: bc

    bc = bary_coord(ll, coord1, coord2, coord3)
    inside = (0d0 < bc(1) .and. bc(1) < 1d0 .and. 0d0 < bc(2) .and. bc(2) < 1d0 .and. 0d0 < bc(3) .and. bc(3) < 1d0)
    if (inside) ival = sum (values*bc)
  end subroutine interp_tria

  function bary_coord (ll, a, b, c)
    implicit none
    real(8), dimension(3) :: bary_coord
    real(8), dimension(2) :: a, b, c, ll

    real(8)               :: det
    real(8), dimension(3) :: bac
    real(8), dimension(2) :: ca, cb, cll

    cb = b - c
    ca = a - c
    cll = ll - c
    det = cb(2)*ca(1) - cb(1)*ca(2)
    bac(1) = ( cb(2)*cll(1) - cb(1)*cll(2)) / det
    bac(2) = (-ca(2)*cll(1) + ca(1)*cll(2)) / det
    bac(3) = 1d0 - bac(1) - bac(2)
    bary_coord = bac
  end function bary_coord

  subroutine fix_boundary (a, b, c, fixed)
    implicit none
    real(8), intent(inout) :: a
    real(8), intent(in)    :: b, c
    integer, intent(out)   :: fixed

    fixed = 0
    if (a < -MATH_PI/2d0 .and. (b > MATH_PI/2d0 .and. c > MATH_PI/2d0)) then
       a = a + MATH_PI*2d0
       fixed = 1
    elseif (a > MATH_PI/2d0 .and. (b < -MATH_PI/2d0 .and. c < -MATH_PI/2d0)) then
       a = a - MATH_PI*2d0
       fixed = -1
    end if
  end subroutine fix_boundary
end module projection_mod
