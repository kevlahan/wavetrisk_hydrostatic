module io_mod
  use geom_mod
  use domain_mod
  use arch_mod
  use adapt_mod
  use smooth_mod
  use comm_mpi_mod
  implicit none
  integer, parameter                  :: N_VAR_OUT = 5
  integer, dimension(2,4)             :: HR_offs
  real(8)                             :: dx_export, dy_export, kx_export, ky_export, vmin, vmax
  real(8), dimension(N_VAR_OUT)       :: minv, maxv
  real, dimension(:,:,:), allocatable :: field2d
  integer                             :: next_fid
  type(Float_Field)                   :: active_level
  data HR_offs /0,0, 1,0, 1,1, 0,1/
contains
  subroutine init_io_mod
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_domain_mod
    next_fid = 100
    initialized = .True.
  end subroutine init_io_mod

  function get_fid()
    integer :: get_fid
    get_fid  = next_fid
    next_fid = next_fid + 1
  end function get_fid

  subroutine write_dual (dom, p, i, j, zlev, offs, dims, fid)
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev, fid
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                   :: d, id, idE, idN, idNE
    integer, dimension(TRIAG) :: leveldual
    real(8), dimension(TRIAG) :: relvort

    d = dom%id + 1

    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    relvort = get_vort (dom, i, j, offs, dims) 

    if (maxval(dom%mask_n%elts((/id, idE, idNE/)+1)) .ge. ADJZONE) then
       ! avoid segfault if pre_levelout not used
       if (allocated(active_level%data)) leveldual(LORT+1) = maxval(active_level%data(d)%elts((/id, idE, idNE/)+1))

       write (fid,'(9(E14.5E2,1X), E14.5E2, 1X, I3)') dom%node%elts((/id, idE, idNE/)+1), relvort(LORT+1), leveldual(LORT+1)
    end if

    if (maxval(dom%mask_n%elts((/id, idNE, idN/)+1)) .ge. ADJZONE) then
       ! avoid segfault if pre_levelout not used
       if (allocated(active_level%data)) leveldual(UPLT+1) = maxval(active_level%data(d)%elts((/id, idNE, idN/)+1))

       write (fid,'(9(E14.5E2,1X), E14.5E2, 1X, I3)') dom%node%elts((/id, idNE, idN/)+1), relvort(UPLT+1), leveldual(UPLT+1)
    end if
  end subroutine write_dual

  subroutine vort_extrema (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id, idN, idE
    real(8) :: vort

    id  = idx(i,   j,   offs, dims)
    idN = idx(i,   j+1, offs, dims)
    idE = idx(i+1, j,   offs, dims)

    if ( dom%mask_e%elts(id*EDGE+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(id*EDGE+UP+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(idN*EDGE+RT+1) .ge. ADJZONE) then

       vort = dom%vort%elts(id*TRIAG+UPLT+1)
       vmin = min(vmin, vort)
       vmax = max(vmax, vort)

    end if

    if ( dom%mask_e%elts(id*EDGE+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(idE*EDGE+UP+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(id*EDGE+RT+1)  .ge. ADJZONE) then

       vort = dom%vort%elts(id*TRIAG+LORT+1)
       vmin = min(vmin, vort)
       vmax = max(vmax, vort)

    end if
  end subroutine vort_extrema

  subroutine write_step (fid, time, k)
    integer ::  fid, k
    real(8) :: time
    
    integer :: l
    real(8) :: tot_mass

    vmin =  1.0e-16
    vmax = -1.0e-16

    do l = level_start, level_end
       call apply_onescale(vort_extrema, l, z_null, 0, 0)
    end do

    tot_mass = integrate_hex(mu, level_start, k)

    if (rank .eq. 0) write(fid,'(E16.9, I3, 2(1X, I9), 7(1X, E16.8), 1X, F16.7)') &
         time, level_end, n_active, tot_mass, get_timing()
  end subroutine write_step

  function integrate_hex (fun, l, k)
    ! Integrate function defined by fun over hexagons
    real(8)  :: integrate_hex
    external :: fun
    integer  :: l, k
    
    integer                        :: d, ll, p, i, j, c, id
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: s, fun

    s = 0.0_8
    do d = 1, size(grid)
       do ll = 1, grid(d)%lev(l)%length
          p = grid(d)%lev(l)%elts(ll)
          call get_offs_Domain (grid(d), p, offs, dims)
          do j = 1, PATCH_SIZE
             do i = 1, PATCH_SIZE
                id = idx(i-1,j-1,offs,dims)
                s = s + fun(grid(d), i-1, j-1, k, offs, dims)/grid(d)%areas%elts(id+1)%hex_inv
             end do
          end do
       end do

       do c = SOUTHEAST, NORTHWEST, 2
          if (.not. grid(d)%pole_master(c/2-2) .or. .not. grid(d)%penta(c)) cycle
          p = 1
          do while (grid(d)%patch%elts(p+1)%level .lt. l)
             p = grid(d)%patch%elts(p+1)%children(c-4)
             if (p .eq. 0) then
                write(6,*) "ERROR(rank=", rank, "):integrate_hex: level incomplete"
                return
             end if
          end do
          call get_offs_Domain (grid(d), p, offs, dims)
          if (c .eq. NORTHWEST) then
             s = s + fun(grid(d), 0, PATCH_SIZE, k, offs, dims)/ &
                  grid(d)%areas%elts(idx(0,PATCH_SIZE,offs,dims)+1)%hex_inv
          else
             s = s + fun(grid(d), PATCH_SIZE, 0, k, offs, dims)/ &
                  grid(d)%areas%elts(idx(PATCH_SIZE,0,offs,dims)+1)%hex_inv
          end if
       end do
    end do
    integrate_hex = sum_real(s)
  end function integrate_hex

  function integrate_tri (fun, k)
    real(8)  :: integrate_tri
    external :: fun
    integer  :: k
    
    integer                        :: d, ll, p, i, j, t, id
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: s, fun

    s = 0.0_8
    do d = 1, size(grid)
       do ll = 1, grid(d)%lev(level_start)%length
          p = grid(d)%lev(level_start)%elts(ll)
          call get_offs_Domain (grid(d), p, offs, dims)
          do j = 1, PATCH_SIZE
             do i = 1, PATCH_SIZE
                id = idx(i-1, j-1, offs, dims)
                do t = LORT, UPLT
                   s = s + fun(grid(d), i-1, j-1, k, t, offs, dims) &
                        *grid(d)%triarea%elts(id*TRIAG+t+1)
                end do
             end do
          end do
       end do
    end do
    integrate_tri = sum_real(s)
  end function integrate_tri

  function only_area (dom, i, j, offs, dims)
    real(8)                        :: only_area
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    only_area = 1.0_8
  end function only_area

  function mu (dom, i, j, zlev, offs, dims)
    real(8)                        :: mu
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx(i, j, offs, dims)
    mu = sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)
  end function mu

  function pot_energy (dom, i, j, zlev, offs, dims)
    real(8)                        :: pot_energy
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx(i, j, offs, dims)
    pot_energy = sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)**2
  end function pot_energy

  function energy (dom, i, j, zlev, offs, dims)
    real(8)                        :: energy
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id

    d = dom%id+1
    id = idx(i, j, offs, dims)
    energy = bernoulli(id+1) * sol(S_MASS,zlev)%data(d)%elts(id+1)
  end function energy

  function tri_only_area (dom, i, j, t, offs, dims)
    real(8)                        :: tri_only_area
    type(Domain)                   :: dom
    integer                        :: i, j, t
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    tri_only_area = 1.0_8
  end function tri_only_area

  function only_coriolis (dom, i, j, t, offs, dims)
    real(8)                        :: only_coriolis
    type(Domain)                   :: dom
    integer                        :: i, j, t
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx(i, j, offs, dims)
    only_coriolis = (dom%coriolis%elts(TRIAG*id+t+1)/dom%triarea%elts(TRIAG*id+t+1))**2
  end function only_coriolis

  subroutine export_2d (proj, values, n_val, fid, l, Nx, Ny, valrange, default_val)
    external                            :: proj
    type(Float_Field), dimension(n_val) :: values
    integer, dimension(2)               :: Nx, Ny
    integer                             :: l, n_val, fid
    real(8), dimension(2)               :: valrange
    real(8), dimension(n_val)           :: default_val
    
    integer                        :: d, k, i, j, v, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    character(5)                   :: fidv
    character(130)                 :: command

    dx_export = valrange(1)/(Nx(2)-Nx(1)+1)
    dy_export = valrange(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1./dx_export
    ky_export = 1./dy_export

    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2),n_val))

    do v = 1, n_val
       field2d(:,:,v) = default_val(v)
    end do

    do d = 1, size(grid)
       do k = 1, grid(d)%lev(l)%length
          call get_offs_Domain(grid(d), grid(d)%lev(l)%elts(k), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,     j,     offs, dims)
                idN  = idx(i,     j + 1, offs, dims)
                idE  = idx(i + 1, j,     offs, dims)
                idNE = idx(i + 1, j + 1, offs, dims)

                call proj(grid(d)%node%elts(id+1),   cC)
                call proj(grid(d)%node%elts(idN+1),  cN)
                call proj(grid(d)%node%elts(idE+1),  cE)
                call proj(grid(d)%node%elts(idNE+1), cNE)

                do v = 1, n_val
                   val   = values(v)%data(d)%elts(id+1)
                   valN  = values(v)%data(d)%elts(idN+1)
                   valE  = values(v)%data(d)%elts(idE+1)
                   valNE = values(v)%data(d)%elts(idNE+1)

                   if (abs(cN(2) - MATH_PI/2) .lt. sqrt(1.0e-15)) then
                      call interp_tri_to_2d_and_fix_bdry(cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/), v)
                      call interp_tri_to_2d_and_fix_bdry((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/), v)
                   else
                      call interp_tri_to_2d_and_fix_bdry(cNE, cN, cC, (/valNE, valN, val/), v)
                   end if
                   if (abs(cE(2) + MATH_PI/2) .lt. sqrt(1.0e-15)) then
                      call interp_tri_to_2d_and_fix_bdry(cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/), v)
                      call interp_tri_to_2d_and_fix_bdry((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/), v)
                   else
                      call interp_tri_to_2d_and_fix_bdry(cC, cE, cNE, (/val, valE, valNE/), v)
                   end if
                end do
             end do
          end do
       end do
    end do

    do v = 1, n_val
       sync_val = default_val(v)
       call sync_array(field2d(Nx(1),Ny(1),v), size(field2d(:,:,v)))
    end do

    if (rank .eq. 0) then
       do v = 1, n_val
          open(fid+v, recl=32768)
          do i = Ny(1),Ny(2)
             write(fid+v,'(2047(E15.6, 1X))') field2d(:,i,v)
          end do
          close(fid+v)
          write(fidv, '(i5)') fid+v
          command = 'bzip2 fort.' // fidv // ' &'
          !call system(command) !JEMF
       end do
    end if

    deallocate(field2d)
  end subroutine export_2d

  subroutine fix_boundary (a, b, c, fixed)
    real(8), intent(inout) :: a
    real(8), intent(in)    :: b, c
    integer, intent(out)   :: fixed
    
    fixed = 0
    if (a .lt. -MATH_PI/2.0_8 .and. (b .gt. MATH_PI/2.0_8 .and. c .gt. MATH_PI/2.0_8)) then
       a = a + MATH_PI*2.0_8
       fixed = 1
    elseif (a .gt. MATH_PI/2.0_8 .and. (b .lt. -MATH_PI/2.0_8 .and. c .lt. -MATH_PI/2.0_8)) then
       a = a - MATH_PI*2.0_8
       fixed = -1
    end if
  end subroutine fix_boundary

  subroutine interp_tri_to_2d_and_fix_bdry (a0, b0, c0, val, v)
    real(8), dimension(2) :: a0, b0, c0
    real(8), dimension(3) :: val
    integer :: v

    integer               :: i
    integer, dimension(3) :: fixed
    real(8), dimension(2) :: a, b, c
    
    a = a0
    b = b0
    c = c0
    call fix_boundary(a(1), b(1), c(1), fixed(1))
    call fix_boundary(b(1), c(1), a(1), fixed(2))
    call fix_boundary(c(1), a(1), b(1), fixed(3))
    call interp_tri_to_2d(a, b, c, val, v)

    if (sum(abs(fixed)) .gt. 1) write(0,*) 'ALARM'

    if (sum(fixed) .ne. 0) then
       a(1) = a(1) - sum(fixed)*MATH_PI*2.0_8
       b(1) = b(1) - sum(fixed)*MATH_PI*2.0_8
       c(1) = c(1) - sum(fixed)*MATH_PI*2.0_8
       call interp_tri_to_2d(a, b, c, val, v)
    end if
  end subroutine interp_tri_to_2d_and_fix_bdry

  subroutine interp_tri_to_2d(a, b, c, val, v)
    real(8), dimension(2) :: a, b, c
    real(8), dimension(3) :: val
    integer :: v
    real(8) minx
    real(8) maxx
    real(8) miny
    real(8) maxy
    integer id_x, id_y
    real(8) ll(2), bac(3), ival
    logical inside

    minx = min(min(a(1), b(1)), c(1))
    maxx = max(max(a(1), b(1)), c(1))
    miny = min(min(a(2), b(2)), c(2))
    maxy = max(max(a(2), b(2)), c(2))
    if (maxx-minx .gt. MATH_PI/2.0_8) then
       write(0,*) 'ERROR(rank', rank, '):io-333 "export"'
       return
    end if

    do id_x = floor(kx_export*minx), ceiling(kx_export*maxx)
       if (id_x .lt. lbound(field2d,1) .or. id_x .gt. ubound(field2d,1)) cycle
       do id_y = floor(ky_export*miny), ceiling(ky_export*maxy)
          if (id_y .lt. lbound(field2d,2) .or. id_y .gt. ubound(field2d,2)) cycle
          ll = (/dx_export*id_x, dy_export*id_y/)
          call interp_tria(ll, a, b, c, val, ival, inside)
          if (inside) field2d(id_x,id_y,v) = ival
       end do
    end do
  end subroutine interp_tri_to_2d

  subroutine interp_tria (ll, coord1, coord2, coord3, values, ival, inside)
    real(8), dimension(2) :: coord1, coord2, coord3
    real(8), dimension(3) :: values
    real(8)               :: ival
    logical               :: inside
    
    real(8), dimension(2) :: ll
    real(8), dimension(3) :: bc

    bc = bary_coord(ll, coord1, coord2, coord3)
    inside = (0.0_8 .lt. bc(1) .and. bc(1) .lt. 1.0_8 .and. &
         0.0_8 .lt. bc(2) .and. bc(2) .lt. 1.0_8 .and. &
         0.0_8 .lt. bc(3) .and. bc(3) .lt. 1.0_8)
    if (inside) ival = sum(values*bc)
  end subroutine interp_tria

  function bary_coord (ll, a, b, c)
    real(8), dimension(3) :: bary_coord
    real(8), dimension(2) :: a, b, c, ll

    real(8)               :: det
    real(8), dimension(3) :: bac
    real(8), dimension(2) :: ca, cb, cll

    cb = b - c
    ca = a - c
    cll = ll - c
    det = cb(2)*ca(1) - cb(1)*ca(2)
    bac(1) = ( cb(2)*cll(1) - cb(1)*cll(2))/det
    bac(2) = (-ca(2)*cll(1) + ca(1)*cll(2))/det
    bac(3) = 1 - bac(1) - bac(2)
    bary_coord = bac
  end function bary_coord

  subroutine write_primal (dom, p, i, j, zlev, offs, dims, fid)
    ! Write primal grid for vertical level zlev
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                       :: d, fid, id, idW, idSW, idS, outl
    real(4), dimension(N_VAR_OUT) :: outv
    real(8), dimension(2)         :: vel_latlon

    d = dom%id + 1

    id   = idx(i,     j,     offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)

    ! Temperature in layer zlev
    outv(1) = sol(S_TEMP,zlev)%data(d)%elts(id+1)/sol(S_MASS,zlev)%data(d)%elts(id+1)*(dom%press%elts(id+1)/ref_press)**kappa

    ! Zonal and meridional velocities
    outv(2) = dom%u_zonal%elts(id+1)
    outv(3) = dom%v_merid%elts(id+1)

    ! Geopotential height at level zlev
    outv(4) = dom%geopot%elts(id+1)/grav_accel

    ! Mass
    outv(5) = sol(S_MASS,zlev)%data(d)%elts(id+1)

    if (allocated(active_level%data)) then ! avoid segfault pre_levelout not used
       outl = nint(active_level%data(dom%id+1)%elts(id+1))
    else
       outl = 0
    end if

    if (dom%mask_n%elts(id+1) .gt. 0) then
       write (fid,'(18(E14.5E2, 1X), 5(E14.5E2, 1X), I3, 1X, I3)') &
            dom%ccentre%elts(TRIAG*id   +LORT+1), dom%ccentre%elts(TRIAG*id   +UPLT+1), &
            dom%ccentre%elts(TRIAG*idW  +LORT+1), dom%ccentre%elts(TRIAG*idSW +UPLT+1), &
            dom%ccentre%elts(TRIAG*idSW +LORT+1), dom%ccentre%elts(TRIAG*idS  +UPLT+1), &
            outv, dom%mask_n%elts(id+1), outl
       where (minv .gt. outv) minv = outv
       where (maxv .lt. outv) maxv = outv
    end if
  end subroutine write_primal

  subroutine zonal_meridional_vel (dom, i, j, offs, dims, zlev, vel_latlon)
    ! Finds lat-lon velocity (with components in zonal and meridional directions) given index information of node
    ! using lapack least squares routine dgels
    type (Domain)                :: dom
    integer                      :: i, j, zlev
    integer, dimension(N_BDRY+1) :: offs
    integer, dimension(2,N_BDRY+1)      :: dims
    real(8), dimension (3)       :: uvw
    real(8), dimension (2)       :: vel_latlon

    integer                     :: d, e, id, id_i, id_e, idN, idE, idNE
    type (Coord)                :: co_node, co_east, co_north, co_northeast, e_merid, e_zonal
    type (Coord), dimension (3) :: dir 
    real(8)                     :: lon, lat

    ! For least squares solver dgels
    integer                    :: info
    real(8), dimension (3,2)   :: A
    integer, parameter         :: lwork = 2*3*2
    real(8), dimension (lwork) :: work

    d = dom%id+1

    id = idx(i, j, offs, dims)

    id_i = id + 1
    id_e = EDGE*id + 1
    idN  = idx(i, j + 1,     offs, dims) + 1
    idE  = idx(i + 1, j,     offs, dims) + 1
    idNE = idx(i + 1, j + 1, offs, dims) + 1

    uvw(1) = sol(S_VELO,zlev)%data(d)%elts(id_e+RT) ! RT velocity
    uvw(2) = sol(S_VELO,zlev)%data(d)%elts(id_e+DG) ! DG velocity
    uvw(3) = sol(S_VELO,zlev)%data(d)%elts(id_e+UP) ! UP velocity

    ! Calculate velocity directions
    co_node      = dom%node%elts(id_i) 
    co_east      = dom%node%elts(idE)
    co_northeast = dom%node%elts(idNE)
    co_north     = dom%node%elts(idN)

    dir(1) = direction (co_node,      co_east)  ! RT direction
    dir(2) = direction (co_northeast, co_node)  ! DG direction
    dir(3) = direction (co_node,      co_north) ! UP direction

    ! Find longitude and latitude coordinates of node
    call cart2sph (co_node, lon, lat)

    e_zonal = Coord (-sin(lon),          cos(lon),           0.0_8)    ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    ! Least squares overdetermined matrix 
    do e = 1, EDGE
       A(e,1) = inner(dir(e), e_zonal)
       A(e,2) = inner(dir(e), e_merid)
    end do

    ! Solve least squares problem to find zonal and meridional velocities
    call dgels ('N', 3, 2, 1, A, 3, uvw, 3, work, lwork, info)

    vel_latlon = uvw(1:2)
  end subroutine zonal_meridional_vel

  function get_vort (dom, i, j, offs, dims)
    ! Averages vorticity to get smooth field for visualization
    real(8), dimension(TRIAG)      :: get_vort
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idN, idE, idS, idW

    id  = idx(i,   j,   offs, dims)
    idE = idx(i+1, j,   offs, dims)
    idN = idx(i,   j+1, offs, dims)
    idW = idx(i-1, j,   offs, dims)
    idS = idx(i,   j-1, offs, dims)

    get_vort(UPLT+1) = ( &
         interp(dom%vort%elts(TRIAG*idW+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)+ &
         
         interp(dom%vort%elts(TRIAG*id+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(id*EDGE+DG+1)+ &
         
         interp(dom%vort%elts(TRIAG*idN+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*idN+RT+1)*dom%pedlen%elts(EDGE*idN+RT+1)) &
         
         / (dom%len%elts(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)+ &
         dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)+ &
         dom%len%elts(EDGE*idN+RT+1)*dom%pedlen%elts(EDGE*idN+RT+1))

    get_vort(LORT+1) = ( &
         interp(dom%vort%elts(TRIAG*idS+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)+ &
         
         interp(dom%vort%elts(TRIAG*id+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(id*EDGE+DG+1)+ &
         
         interp(dom%vort%elts(TRIAG*idE+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*idE+UP+1)*dom%pedlen%elts(EDGE*idE+UP+1)) &
         
         / (dom%len%elts(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)+ &
         dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)+ &
         dom%len%elts(EDGE*idE+UP+1)*dom%pedlen%elts(EDGE*idE+UP+1))
  end function get_vort

  subroutine write_u_wc(dom, p, i, j, offs, dims, fid)
    !write wavelet coefficients of velocity
    type(Domain) dom
    integer p
    integer i, j, k
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer fid
    integer e
    integer id

    id = idx(i, j, offs, dims)

    do k = 1, zlevels
       do e = 1, EDGE
          write(fid,*) wav_coeff(S_VELO,k)%data(dom%id+1)%elts(EDGE*id+e)
       end do
    end do
  end subroutine write_u_wc

  subroutine write_velo (dom, p, i, j, offs, dims, fid)
    type(Domain)                   :: dom
    integer                        :: fid, i, j, p
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, k

    do k = 1, zlevels
       do e = 1, EDGE
          id = idx(i, j, offs, dims)
          write(fid,*) sol(S_VELO,k)%data(dom%id+1)%elts(EDGE*id+e)
       end do
    end do
  end subroutine write_velo

  subroutine read_u_wc_and_mask (dom, p, i, j, offs, dims, fid)
    !read in wavelet coefficients of velocity (JEMF: not mask though??)
    type(Domain)                   :: dom
    integer                        :: fid, i, j, p
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, k

    id = idx(i, j, offs, dims)
    do k = 1, zlevels
       do e = 1, EDGE
          read(fid,*) wav_coeff(S_VELO,k)%data(dom%id+1)%elts(EDGE*id+e)
       end do
    end do
  end subroutine read_u_wc_and_mask

  subroutine write_scalar (dom, p, i, j, zlev, offs, dims, fid)
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, fid, v

    d = dom%id+1
    id = idx(i, j, offs, dims)

    do v = S_MASS, S_TEMP
       write(fid) sol(v,zlev)%data(d)%elts(id+1) ! for pole
    end do
  end subroutine write_scalar

  subroutine read_scalar (dom, p, i, j, zlev, offs, dims, fid)
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, fid, v

    d  = dom%id+1
    id = idx(i, j, offs, dims)

    do v = S_MASS, S_TEMP
       read(fid) sol(v,zlev)%data(d)%elts(id+1) ! for pole
    end do
  end subroutine read_scalar

  subroutine read_mt_wc_and_mask (dom, p, i, j, offs, dims, fid)
    !read in wavelet coefficients of mass and potential temperature (not mask though??)
    type(Domain)                   :: dom
    integer                        :: fid, p, i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, k, v

    id = idx(i, j, offs, dims)

    do k = 1, zlevels
       write(0,*) 'reading k', k, 'id', id
       do v = S_MASS, S_TEMP
          read(fid,*) wav_coeff(v,k)%data(dom%id+1)%elts(id+1)
       end do
    end do
  end subroutine read_mt_wc_and_mask

  subroutine write_mt_wc (dom, p, i, j, offs, dims, fid)
    ! Write wavelet coefficients of mass and potential temperature
    type(Domain)                   :: dom
    integer                        :: fid, p, i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, k, v

    id = idx(i, j, offs, dims)

    do k = 1, zlevels
       write(0,*) 'writing k', k, 'id', id
       do v = S_MASS, S_TEMP
          write (fid,*) wav_coeff(v,k)%data(dom%id+1)%elts(id+1)
       end do
    end do
  end subroutine write_mt_wc

  function dump_adapt_mpi (node_out_rout, edge_out_rout, id, custom_dump)
    ! One file per domain
    integer  :: dump_adapt_mpi
    external :: node_out_rout, edge_out_rout, custom_dump
    integer  :: id

    character(5+4+1+5) :: filename_no, filename_ed, fname_gr
    integer            :: fid_no, l, fid_grid, fid_ed, d, j, k, v, i
    logical            :: child_required(N_CHDRN)
    integer            :: p_par, p_chd, c, p_lev

    dump_adapt_mpi = 0
    fid_no   = id+1000000
    fid_grid = id+3000000

    call update_array_bdry(wav_coeff(S_MASS:S_TEMP,:), NONE)

    do k = 1, zlevels
       do d = 1, size(grid)
          mass => sol(S_MASS,k)%data(d)%elts
          temp => sol(S_TEMP,k)%data(d)%elts
          wc_m => wav_coeff(S_MASS,k)%data(d)%elts
          wc_t => wav_coeff(S_TEMP,k)%data(d)%elts
          call apply_interscale_d (restrict_scalar, grid(d), min_level-1, k, 0, 1) ! +1 to include poles
          nullify (mass, temp, wc_m, wc_t)
       end do
    end do

    do d = 1, size(grid)
       write(filename_no, '(A,I4.4,A,I5.5)')  "coef.", id, "_", glo_id(rank+1,d)
       write(fname_gr,    '(A,I4.4,A,I5.5)')  "grid.", id, "_", glo_id(rank+1,d)

       open(unit=fid_no,   file=filename_no, form="UNFORMATTED", action='WRITE')
       open(unit=fid_grid, file=fname_gr,    form="UNFORMATTED", action='WRITE')

       write(fid_no) istep
       write(fid_no) time

       call custom_dump (fid_no)

       do k = 1, zlevels
          call apply_to_pole_d (write_scalar, grid(d), min_level-1, k, fid_no, .True.)

          do i = MULT(S_VELO)*grid(d)%patch%elts(1+1)%elts_start+1, &
               MULT(S_VELO)*(grid(d)%patch%elts(1+1)%elts_start+PATCH_SIZE**2)

             if (sol(S_VELO,k)%data(d)%elts(i) .ne. sol(S_VELO,k)%data(d)%elts(i)) then
                write(*,*) d, i, 'writeout NaN scal'
                dump_adapt_mpi = 1
                close(fid_no); close(fid_grid);
                return
             end if
          end do

          do v = S_MASS, S_VELO
             write(fid_no) (sol(v,k)%data(d)%elts(i), i=MULT(v)*grid(d)%patch%elts(1+1)%elts_start+1, &
                  MULT(v)*(grid(d)%patch%elts(1+1)%elts_start+PATCH_SIZE**2))
          end do
       end do

       do l = min_level, level_end
          p_lev = 0
          do j = 1, grid(d)%lev(l)%length
             p_par = grid(d)%lev(l)%elts(j)
             if (grid(d)%patch%elts(p_par+1)%deleted) then
                do c = 1, N_CHDRN
                   p_chd = grid(d)%patch%elts(p_par+1)%children(c)
                   if (p_chd .gt. 0) grid(d)%patch%elts(p_chd+1)%deleted = .True.
                end do
                cycle
             end if

             do k = 1, zlevels
                do i = MULT(S_VELO)*grid(d)%patch%elts(p_par+1)%elts_start+1, &
                     MULT(S_VELO)*(grid(d)%patch%elts(p_par+1)%elts_start+PATCH_SIZE**2)

                   if (wav_coeff(S_VELO,k)%data(d)%elts(i) .ne. wav_coeff(S_VELO,k)%data(d)%elts(i)) then
                      write(0,*) grid(d)%patch%elts(p_par+1)%level, 'writeout NaN wav'
                      dump_adapt_mpi = 1
                      close(fid_no); close(fid_grid);
                      return
                   end if
                end do

                do v = S_MASS, S_VELO
                   write(fid_no) (wav_coeff(v,k)%data(d)%elts(i),  &
                        i=MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start+1, &
                        MULT(v)*(grid(d)%patch%elts(p_par+1)%elts_start+PATCH_SIZE**2))
                end do
             end do

             do c = 1, N_CHDRN
                p_chd = grid(d)%patch%elts(p_par+1)%children(c)

                if (p_chd .gt. 0) then

                   child_required(c) = check_child_required(grid(d), p_par, c-1)
                   grid(d)%patch%elts(p_chd+1)%deleted = .not. child_required(c)

                   if (child_required(c)) then
                      p_lev = p_lev + 1
                      grid(d)%lev(l+1)%elts(p_lev) = p_chd
                   end if
                else
                   child_required(c) = .False.
                end if

             end do
             write(fid_grid) child_required
          end do
          if (l+1 .le. max_level) grid(d)%lev(l+1)%length = p_lev
       end do
       close(fid_no); close(fid_grid)
    end do
  end function dump_adapt_mpi

  subroutine load_adapt_mpi (node_in_rout, edge_in_rout, id, custom_load)
    ! One file per domain
    external                             :: node_in_rout, edge_in_rout, custom_load
    integer, dimension(n_domain(rank+1)) :: fid_no, fid_grid !ASCII, fid_ed 
    integer                              :: id

    character(5+4+1+5)          :: filename_no, filename_ed, fname_gr
    integer                     :: c, d, i, j, k, l, v, old_n_patch, p_par, p_chd
    logical, dimension(N_CHDRN) :: child_required

    do d = 1, size(grid)
       fid_no(d)   = id*1000 + 1000000 + d
       fid_grid(d) = id*1000 + 3000000 + d

       write (filename_no, '(A,I4.4,A,I5.5)')  "coef.", id, "_", glo_id(rank+1,d)
       write (fname_gr,    '(A,I4.4,A,I5.5)')  "grid.", id, "_", glo_id(rank+1,d)

       open (unit=fid_no(d),   file=filename_no, form="UNFORMATTED", action='READ')
       open (unit=fid_grid(d), file=fname_gr,    form="UNFORMATTED", action='READ')

       read (fid_no(d)) istep
       read (fid_no(d)) time

       call custom_load (fid_no(d))

       do k = 1, zlevels
          call apply_to_pole_d (read_scalar, grid(d), min_level-1, k, fid_no(d), .True.)

          do v = S_MASS, S_VELO
             read (fid_no(d)) (sol(v,k)%data(d)%elts(i),i = MULT(v)* grid(d)%patch%elts(1+1)%elts_start+1, &
                  MULT(v)*(grid(d)%patch%elts(1+1)%elts_start+PATCH_SIZE**2) )
          end do

          do i = MULT(S_VELO) * grid(d)%patch%elts(1+1)%elts_start+1, &
               MULT(S_VELO)*(grid(d)%patch%elts(1+1)%elts_start+PATCH_SIZE**2)

             if (sol(S_VELO,k)%data(d)%elts(i) .ne. sol(S_VELO,k)%data(d)%elts(i)) then
                write (0,*) d, i, 'Attempt reading in NaN scal -> corrupted checkpoint', id
                stop
             end if
          end do
       end do
    end do

    l = 1
    do while (level_end .gt. l) ! new level was added -> proceed to it
       l = level_end 
       if (rank .eq. 0) write (6,*) 'loading level', l

       do d = 1, size(grid)
          old_n_patch = grid(d)%patch%length
          do j = 1, grid(d)%lev(l)%length
             p_par = grid(d)%lev(l)%elts(j)

             do k = 1, zlevels
                do v = S_MASS, S_VELO
                   read (fid_no(d)) (wav_coeff(v,k)%data(d)%elts(i), &
                        i=MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start+1, &
                        MULT(v)*(grid(d)%patch%elts(p_par+1)%elts_start+PATCH_SIZE**2))
                end do

                do i = MULT(S_VELO)*grid(d)%patch%elts(p_par+1)%elts_start+1, &
                     MULT(S_VELO)*(grid(d)%patch%elts(p_par+1)%elts_start+PATCH_SIZE**2)
                   if (wav_coeff(S_VELO,k)%data(d)%elts(i) .ne. wav_coeff(S_VELO,k)%data(d)%elts(i)) then
                      write (0,*) d, i, 'Attempt reading in NaN wav -> corrupted checkpoint', id
                      stop
                   end if
                end do
             end do

             read (fid_grid(d)) child_required
             do c = 1, N_CHDRN
                if (child_required(c)) then
                   call refine_patch1 (grid(d), p_par, c-1)
                end if
             end do
          end do

          do p_par = 2, old_n_patch
             do c = 1, N_CHDRN
                p_chd = grid(d)%patch%elts(p_par)%children(c)
                if (p_chd+1 .gt. old_n_patch) then
                   call refine_patch2 (grid(d), p_par - 1, c - 1)
                end if
             end do
          end do
       end do
       call post_refine
    end do

    do d = 1, size(grid)
       close(fid_no(d)); close(fid_grid(d))
    end do

    wav_coeff%bdry_uptodate = .False.
  end subroutine load_adapt_mpi

  subroutine default_dump (fid)
    integer :: fid
  end subroutine default_dump

  subroutine default_load (fid)
    integer :: fid
  end subroutine default_load

  subroutine proj_xz_plane (cin, cout)
    type(Coord)          :: cin
    real(8), intent(out) :: cout(2)

    if (cin%y .gt. 0) then
       cout = (/cin%x-radius, cin%z/)
    else
       cout = (/cin%x+radius, cin%z/)
    end if
  end subroutine proj_xz_plane

  subroutine error (msg)
    character(*) :: msg
    write(0,*) "ERROR: ", msg
  end subroutine error

  subroutine read_lonlat_from_binary (arr, n, fid)
    ! use: 
    !     real(8) arr(n_lon,n_lat)
    !     call read_lonlat_from_binary(arr(1,1),n_lon*n_lat,fid)
    integer               :: n, fid
    real(8), dimension(n) :: arr
    
    integer :: i

    read(fid) (arr(i),i=1,n)

  end subroutine read_lonlat_from_binary

  subroutine read_HR_optim_grid
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: d_HR, p, d_glo, d_sub, fid, loz
    character(19+1)                :: filename
    
    maxerror = 0.0_8
    l2error = 0.0_8

    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call apply_onescale2 (ccentre, level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,   level_end-1, z_null, -1, 1)
    call apply_onescale (check_d,  level_end-1, z_null,  0, 0)

    l2error = sqrt(sum_real(l2error))
    maxerror = sync_max_d(maxerror)

    if (rank .eq. 0) write(6,'(A,2(es12.4,1x))') 'Grid quality before optimization:', maxerror, l2error

    fid = get_fid()
    if (level_start .ne. level_end) then
       write(0,*) level_end, level_start
       write(0,*) "Reading HR grid points for `level_start` unequal `level_end` not implemented"
       return
    end if

    write(filename, '(A,I1)')  "../extern/grid_HR/J", level_start-1
    open(unit=fid,file=filename)

    p = 1
    do d_HR = 1, N_ICOSAH_LOZANGE
       loz = dom_id_from_HR_id(d_HR)
       do d_sub = 1, N_SUB_DOM
          d_glo = loz*N_SUB_DOM + sub_dom_id_from_HR_sub_id(d_sub)
          if (owner(d_glo+1) .eq. rank) &
               call get_offs_Domain(grid(loc_id(d_glo+1)+1), p, offs, dims)
          call coord_from_file(d_glo, PATCH_LEVEL, fid, offs, dims, (/0, 0/))
       end do
    end do
    close(fid)

    call comm_nodes3_mpi (get_coord, set_coord, NONE)

    call apply_onescale2 (ccentre,    level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,      level_end-1, z_null, -1, 1)
    call apply_onescale2 (check_grid, level_end-1, z_null,  0, 0)

    maxerror = 0.0_8
    l2error = 0.0_8

    call comm_nodes3_mpi (get_coord, set_coord, NONE)

    call apply_onescale2 (ccentre, level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,   level_end-1, z_null, -1, 1)
    call apply_onescale (check_d,  level_end-1, z_null,  0, 0)

    l2error = sqrt(sum_real(l2error))
    maxerror = sync_max_d(maxerror)
    if (rank .eq. 0) write(*,'(A,2(es12.4,1x))') 'Grid quality (max. diff. primal dual edge bisection [m]):', maxerror, l2error
  end subroutine read_HR_optim_grid

  function dom_id_from_HR_id (d_HR)
    integer :: dom_id_from_HR_id
    ! d_HR: lozange id as used by Heikes & Randall (starts from 1)
    ! results: domain id (starts from 0)
    integer :: d_HR
    
    dom_id_from_HR_id = modulo(d_HR,2)*5 + modulo(d_HR/2-1,5)
  end function dom_id_from_HR_id

  function sub_dom_id_from_HR_sub_id (sub_id)
    integer ::  sub_dom_id_from_HR_sub_id
    ! sub_id: lozange sub id as used by Heikes & Randall (starts from 1)
    ! results: sub domain id (starts from 0)
    integer :: sub_id
    
    integer :: id, i, j, halv_sub_dom, l, jdiv, idiv

    i = 0
    j = 0
    id = sub_id - 1
    halv_sub_dom = N_SUB_DOM/2
    do l = DOMAIN_LEVEL-1, 0, -1
       jdiv = id/halv_sub_dom
       j = j + jdiv*2**l
       id = modulo(id+4**l,4**(l+1))
       idiv = id/halv_sub_dom
       i = i + idiv*2**l
       halv_sub_dom = halv_sub_dom/4
       id = modulo(id,4**l)
    end do
    sub_dom_id_from_HR_sub_id = j*N_SUB_DOM_PER_DIM + i
  end function sub_dom_id_from_HR_sub_id

  subroutine zrotate (c_in, c_out, angle)
    real(8),      intent(in) :: angle
    type(Coord),  intent(in) :: c_in
    type(Coord), intent(out) :: c_out

    c_out%x =  c_in%x*cos(angle) - c_in%y*sin(angle)
    c_out%y =  c_in%x*sin(angle) + c_in%y*cos(angle)
    c_out%z =  c_in%z
  end subroutine zrotate

  recursive subroutine coord_from_file (d_glo, l, fid, offs, dims, ij0)
    integer,                        intent(in) :: d_glo, l, fid
    integer, dimension(2),          intent(in) :: ij0
    integer, dimension(N_BDRY+1),   intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims

    integer :: d_loc, k, ij(2)
    type(Coord) node, node_r

    d_loc = loc_id(d_glo+1)
    do k = 1, 4
       ij = ij0 + HR_offs(:,k)*2**(l-1)
       if (l .eq. 1) then
          ! if domain on other process still read to get to correct position in file
          if (owner(d_glo+1) .eq. rank) then
             read(fid,*) node
             call zrotate(node, node_r, -0.5_8)  ! icosahedron orientation good for tsunami
             grid(d_loc+1)%node%elts(idx(ij(1), ij(2), offs, dims)+1) = project_on_sphere(node_r)
          else
             read(fid,*)
          end if
       else
          call coord_from_file(d_glo, l-1, fid, offs, dims, ij)
       end if
    end do
  end subroutine coord_from_file

  subroutine pre_levelout
    integer :: d, l, max_output_level, num

    ! FIXME cleaner would be to use init_io routine
    call init_Float_Field (active_level, AT_NODE)

    do d = 1, size(grid)
       num = grid(d)%node%length
       call init(active_level%data(d), num)
       active_level%data(d)%elts(1:num) = grid(d)%level%elts(1:num)
    end do

    do l = level_end-1, level_start, -1
       call apply_interscale(restrict_level, l, z_null, 0, 1)
    end do
  end subroutine pre_levelout

  ! now active_level can be used

  subroutine post_levelout
    integer :: d
    
    do d = 1, size(grid)
       deallocate(active_level%data(d)%elts)
    end do
    deallocate(active_level%data)
  end subroutine post_levelout

  subroutine restrict_level (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx(i_par, j_par, offs_par, dims_par)

    if (dom%mask_n%elts(id_chd+1) .ge. ADJZONE) &
         active_level%data(dom%id+1)%elts(id_par+1) = active_level%data(dom%id+1)%elts(id_chd+1)
  end subroutine restrict_level
end module io_mod
