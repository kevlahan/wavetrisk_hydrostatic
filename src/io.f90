module io_mod
  use domain_mod
  use ops_mod
  use smooth_mod
  use comm_mpi_mod
  use adapt_mod
  implicit none

  integer, parameter                   :: N_VAR_OUT = 7
  integer, dimension(2,4)              :: HR_offs
  real(8)                              :: vmin, vmax
  real(8), dimension(N_VAR_OUT)        :: minv, maxv
  integer                              :: next_fid
  type(Float_Field)                    :: active_level
  data HR_offs /0,0, 1,0, 1,1, 0,1/
contains
  subroutine init_io_mod
    implicit none
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_domain_mod
    next_fid = 100
    initialized = .True.
  end subroutine init_io_mod

  integer function get_fid ()
    implicit none

    get_fid  = next_fid
    next_fid = next_fid + 1
  end function get_fid

  subroutine vort_extrema (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idN, idE
    real(8) :: vort

    id  = idx(i,   j,   offs, dims)
    idN = idx(i,   j+1, offs, dims)
    idE = idx(i+1, j,   offs, dims)

    if ( dom%mask_e%elts(id*EDGE+DG+1)  >= ADJZONE .or. &
         dom%mask_e%elts(id*EDGE+UP+1)  >= ADJZONE .or. &
         dom%mask_e%elts(idN*EDGE+RT+1) >= ADJZONE) then

       vort = dom%vort%elts(id*TRIAG+UPLT+1)
       vmin = min(vmin, vort)
       vmax = max(vmax, vort)
    end if

    if ( dom%mask_e%elts(id*EDGE+DG+1)  >= ADJZONE .or. &
         dom%mask_e%elts(idE*EDGE+UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(id*EDGE+RT+1)  >= ADJZONE) then
       vort = dom%vort%elts(id*TRIAG+LORT+1)
       vmin = min (vmin, vort)
       vmax = max (vmax, vort)
    end if
  end subroutine vort_extrema

  subroutine sum_total_mass (initialgo)
    ! Total mass over all vertical layers
    implicit none
    integer :: errcode, ierr
    logical :: initialgo

    integer :: k
    
    if (initialgo) then
       initotalmass = 0.0_8
       do k = 1, zlevels
          initotalmass = initotalmass + integrate_hex (mu, level_start, k)
       end do
    else
       totalmass = 0.0_8
       do k = 1, zlevels
          totalmass = totalmass + integrate_hex (mu, level_start, k)
       end do
       mass_error = abs (totalmass-initotalmass)/initotalmass
    end if
  end subroutine sum_total_mass

  real(8) function mu (dom, i, j, zlev, offs, dims)
    ! Defines mass for total mass integration
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims)
    mu = sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)
  end function mu

  real(8) function integrate_hex (fun, l, k)
    ! Integrate function defined by fun over hexagons
    implicit none
    integer  :: l, k

    integer                        :: d, ll, p, i, j, c, id
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: s, val

    interface
       real(8) function fun (dom, i, j, zlev, offs, dims)
         import
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function fun
    end interface

    s = 0.0_8
    do d = 1, size(grid)
       do ll = 1, grid(d)%lev(l)%length
          p = grid(d)%lev(l)%elts(ll)
          call get_offs_Domain (grid(d), p, offs, dims)
          do j = 1, PATCH_SIZE
             do i = 1, PATCH_SIZE
                id = idx (i-1, j-1, offs, dims)
                val = fun (grid(d), i-1, j-1, k, offs, dims)
                s = s + val/grid(d)%areas%elts(id+1)%hex_inv
             end do
          end do
       end do

       do c = SOUTHEAST, NORTHWEST, 2
          if (.not. grid(d)%pole_master(c/2-2) .or. .not. grid(d)%penta(c)) cycle
          p = 1
          do while (grid(d)%patch%elts(p+1)%level < l)
             p = grid(d)%patch%elts(p+1)%children(c-4)
             if (p == 0) then
                write (6,'(A, i4, A)') "ERROR(rank = ", rank, "): integrate_hex: level incomplete"
                return
             end if
          end do
          call get_offs_Domain (grid(d), p, offs, dims)
          if (c == NORTHWEST) then
             id = idx (0, PATCH_SIZE, offs, dims)
             val = fun (grid(d), 0, PATCH_SIZE, k, offs, dims)
             s = s + val/grid(d)%areas%elts(id+1)%hex_inv
          else
             id = idx (PATCH_SIZE, 0, offs, dims)
             val = fun (grid(d), PATCH_SIZE, 0, k, offs, dims)
             s = s + val/grid(d)%areas%elts(id+1)%hex_inv
          end if
       end do
    end do
    integrate_hex = sum_real (s)
  end function integrate_hex

  real(8) function integrate_tri (fun, k)
    implicit none
    integer  :: k

    integer                        :: d, l, ll, p, i, j, t, id
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: s

    interface
       real(8) function fun (dom, i, j, zlev, t, offs, dims)
         import 
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, t, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function fun
    end interface

    s = 0.0_8
    do l = level_start, level_end
       do d = 1, size(grid)
          do ll = 1, grid(d)%lev(l)%length
             p = grid(d)%lev(l)%elts(ll)
             call get_offs_Domain (grid(d), p, offs, dims)
             do j = 1, PATCH_SIZE
                do i = 1, PATCH_SIZE
                   id = idx (i-1, j-1, offs, dims)
                   do t = LORT, UPLT
                      s = s + fun (grid(d), i-1, j-1, k, t, offs, dims) * grid(d)%triarea%elts(TRIAG*id+t+1)
                   end do
                end do
             end do
          end do
       end do
    end do
    integrate_tri = sum_real (s)
  end function integrate_tri
  
  real(8) function adaptive_area (routine)
    ! Integrates value define in routine over adaptive grid
    implicit none
    real(8), external :: routine

    integer ::  l
        
    hex_int = 0.0_8
    do l = level_start, level_end-1
       call fine_hex_area (only_area, l)
       call coarse_hex_area (only_area, l)
    end do
    call fine_hex_area (only_area, level_end)
    adaptive_area = sum_real (hex_int)
  end function adaptive_area
  
  subroutine coarse_hex_area (routine, l)
    ! Remove cells that are not at locally finest scale
    implicit none
    integer :: l
    real(8), external :: routine
    
    integer                          :: c, d, i, id, i_par, j, j_par, jj,  p_chd, p_par
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry

    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          p_par = grid(d)%lev(l)%elts(jj)
          call get_offs_Domain (grid(d), p_par, offs_par, dims_par)

          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par+1)%children(c)
             if (p_chd == 0) cycle

             call get_offs_Domain (grid(d), p_chd, offs_chd, dims_chd, inner_bdry)

             bdry = (/0, 0, 0, 0/)

             where (inner_bdry) bdry = 0

             do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
                j_par = j-1 + chd_offs(2,c)
                do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
                   i_par = i-1 + chd_offs(1,c)
                   id = idx (i_par, j_par, offs_par, dims_par)
                   hex_int = hex_int - routine(id)/grid(d)%areas%elts(id+1)%hex_inv
                end do
             end do
          end do
       end do
    end do
  end subroutine coarse_hex_area

  subroutine fine_hex_area (routine, l)
    implicit none
    integer :: l
    real(8), external :: routine
    
    integer                          :: d, i, id, j, jj, p
    integer, dimension(N_BDRY+1)     :: offs
    integer, dimension(2,N_BDRY+1)   :: dims
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry

    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          p = grid(d)%lev(l)%elts(jj)
          call get_offs_Domain (grid(d), p, offs, dims, inner_bdry)

          bdry = (/0, 0, 0, 0/)

          where (inner_bdry) bdry = 0

          do j = bdry(JMINUS) + 1, PATCH_SIZE + bdry(JPLUS)
             do i = bdry(IMINUS) + 1, PATCH_SIZE + bdry(IPLUS)
                id = idx (i, j, offs, dims)
                hex_int = hex_int + routine(id)/grid(d)%areas%elts(id+1)%hex_inv
             end do
          end do
       end do
    end do
  end subroutine fine_hex_area

  real(8) function only_area (id)
    implicit none
    integer :: id
    
    only_area = 1.0_8
  end function only_area

  function pot_energy (dom, i, j, zlev, offs, dims)
    implicit none
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
    implicit none
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

  real(8) function tri_only_area (dom, i, j, zlev, t, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev, t
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    tri_only_area = 1.0_8
  end function tri_only_area

  real(8) function only_coriolis (dom, i, j, t, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, t
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx(i, j, offs, dims)
    only_coriolis = (dom%coriolis%elts(TRIAG*id+t+1)/dom%triarea%elts(TRIAG*id+t+1))**2
  end function only_coriolis

  subroutine interp_tria (ll, coord1, coord2, coord3, values, ival, inside)
    implicit none
    real(8), dimension(2) :: coord1, coord2, coord3
    real(8), dimension(3) :: values
    real(8)               :: ival
    logical               :: inside

    real(8), dimension(2) :: ll
    real(8), dimension(3) :: bc

    bc = bary_coord(ll, coord1, coord2, coord3)
    inside = (0.0_8 < bc(1) .and. bc(1) < 1.0_8 .and. 0.0_8 < bc(2) .and. bc(2) < 1.0_8 .and. 0.0_8 < bc(3) .and. bc(3) < 1.0_8)
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
    bac(1) = ( cb(2)*cll(1) - cb(1)*cll(2))/det
    bac(2) = (-ca(2)*cll(1) + ca(1)*cll(2))/det
    bac(3) = 1 - bac(1) - bac(2)
    bary_coord = bac
  end function bary_coord

  subroutine fix_boundary (a, b, c, fixed)
    implicit none
    real(8), intent(inout) :: a
    real(8), intent(in)    :: b, c
    integer, intent(out)   :: fixed

    fixed = 0
    if (a < -MATH_PI/2.0_8 .and. (b > MATH_PI/2.0_8 .and. c > MATH_PI/2.0_8)) then
       a = a + MATH_PI*2.0_8
       fixed = 1
    elseif (a > MATH_PI/2.0_8 .and. (b < -MATH_PI/2.0_8 .and. c < -MATH_PI/2.0_8)) then
       a = a - MATH_PI*2.0_8
       fixed = -1
    end if
  end subroutine fix_boundary

  subroutine cart2sph2 (cin, cout)
    implicit none
    type(Coord)                        :: cin
    real(8), dimension(2), intent(out) :: cout

    call cart2sph (cin, cout(1), cout(2))
  end subroutine cart2sph2

  subroutine cal_temp (dom, i, j, zlev, offs, dims)
    ! Compute temperature in compressible case
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, d, k
    real(8), dimension(zlevels) :: p

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    ! Integrate the pressure upwards
    p(1) = dom%surf_press%elts(id+1) - 0.5*grav_accel*sol(S_MASS,1)%data(d)%elts(id+1)
    do k = 2, zlevels
       p(k) = p(k-1) - grav_accel*interp(sol(S_MASS,k)%data(d)%elts(id+1), sol(S_MASS,k-1)%data(d)%elts(id+1))
    end do

    ! Temperature at all vertical levels (saved in exner_fun)
    do k = 1, zlevels
       exner_fun(k)%data(d)%elts(id+1) = sol(S_TEMP,k)%data(d)%elts(id+1)/sol(S_MASS,k)%data(d)%elts(id+1) * (p(k)/p_0)**kappa
    end do

    ! temperature at save levels (saved in horiz_flux)
    do k = 1, save_levels
       horiz_flux(k)%data(d)%elts(id+1) = sol_save(S_TEMP,k)%data(d)%elts(id+1)/sol_save(S_MASS,k)%data(d)%elts(id+1) * &
            (pressure_save(k)/p_0)**kappa
    end do
  end subroutine cal_temp

  subroutine cal_geopot (dom, i, j, zlev, offs, dims)
    ! Compute geopotential in compressible case
    ! Assumes that temperature has already been calculated and stored in exner_fun
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, d, k
    real(8) :: pressure_lower, pressure_upper

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    ! Integrate geopotential upwards from surface
    pressure_lower = dom%surf_press%elts(id+1)
    pressure_upper = pressure_lower - grav_accel*sol(S_MASS,1)%data(d)%elts(id+1)
    dom%geopot_lower%elts(id+1) = surf_geopot (dom%node%elts(id+1))/grav_accel

    k = 1
    do while (pressure_upper > pressure_save(1))
       dom%geopot_lower%elts(id+1) = dom%geopot_lower%elts(id+1) + &
            R_d/grav_accel * exner_fun(k)%data(d)%elts(id+1) * (log(pressure_lower)-log(pressure_upper))

       k = k+1
       pressure_lower = pressure_upper
       pressure_upper = pressure_lower - grav_accel*sol(S_MASS,k+1)%data(d)%elts(id+1)
    end do

    ! Add additional contribution up to pressure level pressure_save(1)
    dom%geopot_lower%elts(id+1) = dom%geopot_lower%elts(id+1) &
         + R_d/grav_accel * exner_fun(k)%data(d)%elts(id+1) * (log(pressure_lower)-log(pressure_save(1)))
  end subroutine cal_geopot

  subroutine statistics
    ! Calculates zonal statistics
    use domain_mod
    implicit none
    integer :: d, k, p

    if (rank == 0) write (6,'(a)') 'Incrementing zonal averages'

    call cal_surf_press (sol)

    do k = 1, zlevels
       do d = 1, size (grid)
          mass => sol(S_MASS,k)%data(d)%elts
          temp => sol(S_TEMP,k)%data(d)%elts
          velo => sol(S_VELO,k)%data(d)%elts
          do p = 2, grid(d)%patch%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), p-1, k, 0, 0)
             call apply_onescale_to_patch (cal_pressure,   grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (cal_zonal_avg,  grid(d), p-1, k, 0, 0)
          end do
          nullify (mass, temp, velo)
       end do
    end do
  end subroutine statistics

  subroutine cal_zonal_avg (dom, i, j, zlev, offs, dims)
    ! Zonal average means and covariances over all checkpoints using stable online algorithm
    ! Uses Welford's stable onlne algorithm
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: bin, d, id_i
    real(8) :: lat, lon, temperature, Tprime, Uprime, Vprime, Tprime_new, Uprime_new, Vprime_new

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    call cart2sph (dom%node%elts(id_i), lon, lat)

    temperature = (temp(id_i)/mass(id_i)) * (dom%press%elts(id_i)/p_0)**kappa

    do bin = 1, nbins
       if (lat*1.8d2/MATH_PI <= bounds(bin)) then
          Nstats(zlev,bin) = Nstats(zlev,bin) + 1

          Tprime = temperature            - zonal_avg(zlev,bin,1)
          Uprime = dom%u_zonal%elts(id_i) - zonal_avg(zlev,bin,3)
          Vprime = dom%v_merid%elts(id_i) - zonal_avg(zlev,bin,4)

          ! Mean values
          zonal_avg(zlev,bin,1) = zonal_avg(zlev,bin,1) + Tprime/Nstats(zlev,bin)
          zonal_avg(zlev,bin,3) = zonal_avg(zlev,bin,3) + Uprime/Nstats(zlev,bin)
          zonal_avg(zlev,bin,4) = zonal_avg(zlev,bin,4) + Vprime/Nstats(zlev,bin)
          zonal_avg(zlev,bin,5) = zonal_avg(zlev,bin,5) +  0.5 * (Uprime**2 + Vprime**2)/Nstats(zlev,bin)

          Tprime_new = temperature            - zonal_avg(zlev,bin,1)
          Uprime_new = dom%u_zonal%elts(id_i) - zonal_avg(zlev,bin,3)
          Vprime_new = dom%v_merid%elts(id_i) - zonal_avg(zlev,bin,4)

          ! Update sums of squares (for variance calculation)

          ! Temperature 
          zonal_avg(zlev,bin,2) = zonal_avg(zlev,bin,2) + Tprime * Tprime_new

          ! Eddy momentum flux (covariance)
          zonal_avg(zlev,bin,6) = zonal_avg(zlev,bin,6) + Uprime * Vprime_new

          ! Zonal velocity variance
          zonal_avg(zlev,bin,7) = zonal_avg(zlev,bin,7) + Uprime * Uprime_new 

          ! Meridional velocity variance
          zonal_avg(zlev,bin,8) = zonal_avg(zlev,bin,8) + Vprime * Vprime_new 

          ! Eddy heat flux (covariance)
          zonal_avg(zlev,bin,9) = zonal_avg(zlev,bin,9) + Vprime * Tprime_new

          exit
       end if
    end do
  end subroutine cal_zonal_avg

  subroutine write_out_stats
    ! Writes out zonal average statistics
    implicit none
    integer            :: ibin, k, v
    integer, parameter :: funit = 400
    character(2)       :: var_file
    character(130)     :: command

    write (6,'(a)') 'Saving statistics'

    ! Find sample covariances from sums of squares
    zonal_avg_glo(:,:,2) = zonal_avg_glo(:,:,2) / (Nstats_glo - 1)
    zonal_avg_glo(:,:,6) = zonal_avg_glo(:,:,6) / (Nstats_glo - 1)
    zonal_avg_glo(:,:,7) = zonal_avg_glo(:,:,7) / (Nstats_glo - 1)
    zonal_avg_glo(:,:,8) = zonal_avg_glo(:,:,8) / (Nstats_glo - 1)
    zonal_avg_glo(:,:,9) = zonal_avg_glo(:,:,9) / (Nstats_glo - 1)

    ! Save number of data points
    write (var_file, '(i2.2)') 00
    open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="UNFORMATTED", action='WRITE')
    write (funit) Nstats_glo
    close (funit)

    ! Save zonal statistics (means and covariances)
    do v = 1, nvar_zonal
       write (var_file, '(i2)') v+10
       open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='WRITE')
       do k = zlevels,1,-1
          write (funit,'(2047(E15.6, 1X))') zonal_avg_glo(k,:,v)
       end do
       close (funit)
    end do

    ! Longitude values (dummy)
    write (var_file, '(i2)') 20
    open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='WRITE') 
    write (funit,'(2047(E15.6, 1X))') bounds
    close (funit)

    ! Latitude values
    write (var_file, '(i2)') 21
    open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='WRITE') 
    write (funit,'(2047(E15.6, 1X))') bounds - dbin/2
    close (funit)

    ! Non-dimensional pressure based vertical coordinates p_k/p_s
    write (var_file, '(i2)') 22
    open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='WRITE') 
    write (funit,'(2047(E15.6, 1X))') (0.5*((a_vert(k)+a_vert(k+1))/p_0 + b_vert(k)+b_vert(k+1)), k = zlevels, 1, -1)
    close (funit)

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.3.?? > tmp' 
    call system (command)
    command = 'tar czf '//trim(run_id)//'.3.tgz -T tmp --remove-files &'
    call system (command)
  end subroutine write_out_stats

  subroutine vort_triag_to_hex (dom, i, j, zlev, offs, dims)
    ! Approximate vorticity at hexagon points
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idW, idSW, idS, d

    d = dom%id + 1
    id   = idx(i,   j,   offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)

    dom%press_lower%elts(id+1) = ( &
         dom%areas%elts(id+1)%part(1)*dom%vort%elts(TRIAG*id+LORT+1)   + &
         dom%areas%elts(id+1)%part(2)*dom%vort%elts(TRIAG*id+UPLT+1)   + &
         dom%areas%elts(id+1)%part(3)*dom%vort%elts(TRIAG*idW+LORT+1)  + &
         dom%areas%elts(id+1)%part(4)*dom%vort%elts(TRIAG*idSW+UPLT+1) + &
         dom%areas%elts(id+1)%part(5)*dom%vort%elts(TRIAG*idSW+LORT+1) + &
         dom%areas%elts(id+1)%part(6)*dom%vort%elts(TRIAG*idS+UPLT+1)    &
         ) * dom%areas%elts(id+1)%hex_inv
  end subroutine vort_triag_to_hex

  subroutine write_primal (dom, p, i, j, zlev, offs, dims, funit)
    ! Write primal grid for vertical level zlev
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev, funit
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                       :: d, id, idW, idSW, idS, outl
    real(4), dimension(N_VAR_OUT) :: outv
    real(8), dimension(2)         :: vel_latlon

    d = dom%id + 1

    id   = idx(i,   j,   offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)

    ! Temperature in layer zlev
    outv(1) = sol(S_TEMP,zlev)%data(d)%elts(id+1)/sol(S_MASS,zlev)%data(d)%elts(id+1)*(dom%press%elts(id+1)/p_0)**kappa

    ! Zonal and meridional velocities
    outv(2) = dom%u_zonal%elts(id+1)
    outv(3) = dom%v_merid%elts(id+1)

    ! Geopotential height at level zlev
    outv(4) = dom%geopot%elts(id+1)/grav_accel

    ! Mass
    outv(5) = sol(S_MASS,zlev)%data(d)%elts(id+1)

    ! Surface pressure
    outv(6) = dom%surf_press%elts(id+1)

    ! Vorticity at hexagon points
    outv(7) =  dom%press_lower%elts(id+1)
    ! (dom%areas%elts(id+1)%part(1)*dom%vort%elts(TRIAG*id+LORT+1) + &
    ! dom%areas%elts(id+1)%part(2)*dom%vort%elts(TRIAG*id+UPLT+1) + &
    ! dom%areas%elts(idW+1)%part(3)*dom%vort%elts(TRIAG*idW+LORT+1) + &
    ! dom%areas%elts(idSW+1)%part(4)*dom%vort%elts(TRIAG*idSW+UPLT+1) + &
    ! dom%areas%elts(idSW+1)%part(5)*dom%vort%elts(TRIAG*idSW+LORT+1) + &
    ! dom%areas%elts(idS+1)%part(6)*dom%vort%elts(TRIAG*idS+UPLT+1)) * dom%areas%elts(id+1)%hex_inv

    if (allocated(active_level%data)) then ! avoid segfault pre_levelout not used
       outl = nint(active_level%data(d)%elts(id+1))
    else
       outl = 0
    end if

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       write (funit,'(18(E14.5E2, 1X), 7(E14.5E2, 1X), I3, 1X, I3)') &
            dom%ccentre%elts(TRIAG*id  +LORT+1), dom%ccentre%elts(TRIAG*id  +UPLT+1), &
            dom%ccentre%elts(TRIAG*idW +LORT+1), dom%ccentre%elts(TRIAG*idSW+UPLT+1), &
            dom%ccentre%elts(TRIAG*idSW+LORT+1), dom%ccentre%elts(TRIAG*idS +UPLT+1), &
            outv, dom%mask_n%elts(id+1), outl
       where (minv > outv) minv = outv
       where (maxv < outv) maxv = outv
    end if
  end subroutine write_primal

  subroutine write_dual (dom, p, i, j, zlev, offs, dims, funit)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev, funit
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                   :: d, id, idE, idN, idNE
    integer, dimension(TRIAG) :: leveldual
    real(8), dimension(TRIAG) :: relvort

    d = dom%id + 1

    id   = idx(i,   j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    relvort = get_vort (dom, i, j, offs, dims)

    if (maxval(dom%mask_n%elts((/id, idE, idNE/)+1)) >= ADJZONE) then
       ! avoid segfault if pre_levelout not used
       if (allocated(active_level%data)) leveldual(LORT+1) = maxval(active_level%data(d)%elts((/id, idE, idNE/)+1))

       write (funit,'(9(E14.5E2,1X), E14.5E2, 1X, I3)') dom%node%elts((/id, idE, idNE/)+1), relvort(LORT+1), leveldual(LORT+1)
    end if

    if (maxval(dom%mask_n%elts((/id, idNE, idN/)+1)) >= ADJZONE) then
       ! avoid segfault if pre_levelout not used
       if (allocated(active_level%data)) leveldual(UPLT+1) = maxval(active_level%data(d)%elts((/id, idNE, idN/)+1))

       write (funit,'(9(E14.5E2,1X), E14.5E2, 1X, I3)') dom%node%elts((/id, idNE, idN/)+1), relvort(UPLT+1), leveldual(UPLT+1)
    end if
  end subroutine write_dual

  function get_vort (dom, i, j, offs, dims)
    ! Averages vorticity to get smooth field for visualization
    implicit none
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

  subroutine write_u_wc (dom, p, i, j, offs, dims, fid)
    ! Write wavelet coefficients of velocity
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, p
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, fid, id, k

    id = idx(i, j, offs, dims)

    do k = 1, zlevels
       do e = 1, EDGE
          write (fid,*) wav_coeff(S_VELO,k)%data(dom%id+1)%elts(EDGE*id+e)
       end do
    end do
  end subroutine write_u_wc

  subroutine write_velo (dom, p, i, j, offs, dims, fid)
    implicit none
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

  subroutine write_scalar (dom, p, i, j, zlev, offs, dims, fid)
    implicit none
    type(Domain)                   :: dom
    integer                        :: fid, i, j, p, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, v

    d = dom%id+1
    id = idx(i, j, offs, dims)

    do v = S_MASS, S_TEMP
       write (fid) sol(v,zlev)%data(d)%elts(id+1) ! for pole
       write (fid) trend(v,zlev)%data(d)%elts(id+1) ! for pole`
    end do
  end subroutine write_scalar

  subroutine read_scalar (dom, p, i, j, zlev, offs, dims, fid)
    implicit none
    type(Domain)                   :: dom
    integer                        :: fid, i, j, p, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, v

    d  = dom%id+1
    id = idx(i, j, offs, dims)

    do v = S_MASS, S_TEMP
       read (fid) sol(v,zlev)%data(d)%elts(id+1) ! for pole
       read (fid) trend(v,zlev)%data(d)%elts(id+1) ! for pole
    end do
  end subroutine read_scalar

  subroutine dump_adapt_mpi (id, custom_dump, run_id)
    ! Save data in check point files for restart
    ! One file per domain
    ! NOTE: modifies grid structure
    implicit none
    integer      :: id
    external     :: custom_dump
    character(*) :: run_id

    integer                               :: c, d, i, ibeg, iend, j, k, l, p_chd, p_lev, p_par, v
    integer, dimension(1:size(grid))      :: fid_no, fid_gr
    character(255)                        :: filename_gr, filename_no
    logical, dimension(1:N_CHDRN)         :: required
    type(Domain), dimension(1:size(grid)) :: grid_tmp

    interface
       subroutine custom_dump (fid)
         implicit none
         integer :: fid
       end subroutine custom_dump
    end interface

    ! Make copy of grid since grid is modified when saving checkpoint
    grid_tmp = grid
    
    fid_no = id+1000000
    fid_gr = id+3000000

    sol%bdry_uptodate             = .false.
    trend%bdry_uptodate           = .false.
    wav_coeff%bdry_uptodate       = .false.
    trend_wav_coeff%bdry_uptodate = .false.
    call update_array_bdry (sol,             NONE)
    call update_array_bdry (trend,           NONE)
    call update_array_bdry (wav_coeff,       NONE)
    call update_array_bdry (trend_wav_coeff, NONE)

    do k = 1, zlevels
       do d = 1, size(grid)
          mass => sol(S_MASS,k)%data(d)%elts
          temp => sol(S_TEMP,k)%data(d)%elts
          wc_m => wav_coeff(S_MASS,k)%data(d)%elts
          wc_t => wav_coeff(S_TEMP,k)%data(d)%elts
          call apply_interscale_d (restrict_scalar, grid(d), min_level-1, k, 0, 1) ! +1 to include poles
          nullify (mass, temp, wc_m, wc_t)

          mass => trend(S_MASS,k)%data(d)%elts
          temp => trend(S_TEMP,k)%data(d)%elts
          wc_m => trend_wav_coeff(S_MASS,k)%data(d)%elts
          wc_t => trend_wav_coeff(S_TEMP,k)%data(d)%elts
          call apply_interscale_d (restrict_scalar, grid(d), min_level-1, k, 0, 1) ! +1 to include poles
          nullify (mass, temp, wc_m, wc_t)
       end do
    end do

    do d = 1, size(grid)
       write (filename_no, '(A,A,I4.4,A,I5.5)') trim (run_id), "_coef.", id, "_", glo_id(rank+1,d)
       write (filename_gr, '(A,A,I4.4,A,I5.5)') trim (run_id), "_grid.", id, "_", glo_id(rank+1,d)

       open (unit=fid_no(d), file=trim(filename_no), form="UNFORMATTED", action='WRITE')
       open (unit=fid_gr(d), file=trim(filename_gr), form="UNFORMATTED", action='WRITE')

       write (fid_no(d)) istep
       write (fid_no(d)) time

       call custom_dump (fid_no(d))

       ! Write data at coarsest scale (scaling functions)
       p_par = 1
       do k = 1, zlevels
          call apply_to_pole_d (write_scalar, grid_tmp(d), min_level-1, k, fid_no(d), .true.)
          do v = S_MASS, S_VELO
             ibeg = MULT(v)*grid_tmp(d)%patch%elts(p_par+1)%elts_start + 1
             iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1
             write (fid_no(d))   sol(v,k)%data(d)%elts(ibeg:iend)
             write (fid_no(d)) trend(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do

       ! Write wavelets at finer scales
       do l = min_level, level_end
          p_lev = 0
          do j = 1, grid_tmp(d)%lev(l)%length
             p_par = grid_tmp(d)%lev(l)%elts(j)
             if (grid_tmp(d)%patch%elts(p_par+1)%deleted) then
                do c = 1, N_CHDRN
                   p_chd = grid_tmp(d)%patch%elts(p_par+1)%children(c)
                   if (p_chd > 0) grid_tmp(d)%patch%elts(p_chd+1)%deleted = .true. 
                end do
                cycle ! No data to write
             end if
 
            do k = 1, zlevels
                do v = S_MASS, S_VELO
                   ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
                   iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1
                   write (fid_no(d))       wav_coeff(v,k)%data(d)%elts(ibeg:iend)
                   write (fid_no(d)) trend_wav_coeff(v,k)%data(d)%elts(ibeg:iend)
                end do
             end do

             ! Record whether patch needs to be refined
             do c = 1, N_CHDRN
                p_chd = grid_tmp(d)%patch%elts(p_par+1)%children(c)
                if (p_chd > 0) then
                   required(c) = check_child_required(grid_tmp(d), p_par, c-1)
                   grid_tmp(d)%patch%elts(p_chd+1)%deleted = .not. required(c) 
                   if (required(c)) then
                      p_lev = p_lev + 1
                      grid_tmp(d)%lev(l+1)%elts(p_lev) = p_chd ! ** grid modified **
                   end if
                else
                   required(c) = .false.
                end if
             end do
             write (fid_gr(d)) required
          end do
          if (l+1 <= max_level) grid_tmp(d)%lev(l+1)%length = p_lev ! ** grid modified **
       end do
       close (fid_no(d))
       close (fid_gr(d))
    end do
  end subroutine dump_adapt_mpi

  subroutine load_adapt_mpi (id, custom_load, run_id)
    ! Read data from check point files for restart
    ! One file per domain
    implicit none
    integer      :: id
    external     :: custom_load
    character(*) :: run_id

    character(255)                       :: filename_gr, filename_no
    integer                              :: c, d, i, ibeg, iend, j, k, l, old_n_patch, p_chd, p_par, v
    integer, dimension(1:size(grid))     :: fid_no, fid_gr
    logical, dimension(N_CHDRN)          :: required

     interface
       subroutine custom_load (fid)
         implicit none
         integer :: fid
       end subroutine custom_load
    end interface

    ! Load coarsest scale solution (scaling functions)
    do d = 1, size(grid)
       fid_no(d) = id*1000 + 1000000 + d
       fid_gr(d) = id*1000 + 3000000 + d

       write (filename_no, '(A,A,I4.4,A,I5.5)') trim (run_id), "_coef.", id, "_", glo_id(rank+1,d)
       write (filename_gr, '(A,A,I4.4,A,I5.5)') trim (run_id), "_grid.", id, "_", glo_id(rank+1,d)

       open (unit=fid_no(d), file=trim(filename_no), form="UNFORMATTED", action='READ')
       open (unit=fid_gr(d), file=trim(filename_gr), form="UNFORMATTED", action='READ')

       read (fid_no(d)) istep
       read (fid_no(d)) time

       call custom_load (fid_no(d))

       p_par = 1
       do k = 1, zlevels
          call apply_to_pole_d (read_scalar, grid(d), min_level-1, k, fid_no(d), .true.)
          do v = S_MASS, S_VELO
             ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
             iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1
             read (fid_no(d))   sol(v,k)%data(d)%elts(ibeg:iend)
             read (fid_no(d)) trend(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do
    end do
    
    ! Load finer scales (wavelets) if present
    ! (level_end is initially level_start and is incremented by refine_patch1 if children are present)
    l = 1
    do while (level_end > l) ! New level was added -> proceed to it
       l = level_end 
       if (rank == 0) write (6,'(a,i2)') 'Loading level ', l
       do d = 1, size(grid)
          old_n_patch = grid(d)%patch%length
          do j = 1, grid(d)%lev(l)%length
             p_par = grid(d)%lev(l)%elts(j)
             do k = 1, zlevels
                do v = S_MASS, S_VELO
                   ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
                   iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1
                   read (fid_no(d))       wav_coeff(v,k)%data(d)%elts(ibeg:iend)
                   read (fid_no(d)) trend_wav_coeff(v,k)%data(d)%elts(ibeg:iend)
                end do
             end do

             read (fid_gr(d)) required
             do c = 1, N_CHDRN
                if (required(c)) call refine_patch1 (grid(d), p_par, c-1)
             end do
          end do
          do p_par = 2, old_n_patch
             do c = 1, N_CHDRN
                p_chd = grid(d)%patch%elts(p_par)%children(c)
                if (p_chd+1 > old_n_patch) call refine_patch2 (grid(d), p_par-1, c-1)
             end do
          end do
       end do
       call post_refine
    end do

    do d = 1, size(grid)
       close(fid_no(d)); close(fid_gr(d))
    end do

    sol%bdry_uptodate             = .false.
    trend%bdry_uptodate           = .false.
    wav_coeff%bdry_uptodate       = .false.
    trend_wav_coeff%bdry_uptodate = .false.
    call update_array_bdry (sol,             NONE)
    call update_array_bdry (trend,           NONE)
    call update_array_bdry (wav_coeff,       NONE)
    call update_array_bdry (trend_wav_coeff, NONE)
  end subroutine load_adapt_mpi

  subroutine proj_xz_plane (cin, cout)
    implicit none
    type(Coord)                        :: cin
    real(8), dimension(2), intent(out) :: cout

    if (cin%y > 0) then
       cout = (/cin%x-radius, cin%z/)
    else
       cout = (/cin%x+radius, cin%z/)
    end if
  end subroutine proj_xz_plane

  subroutine error (msg)
    implicit none
    character(*) :: msg
    write(0,*) "ERROR: ", msg
  end subroutine error

  subroutine read_lonlat_from_binary (arr, n, fid)
    !     Use: real(8) arr(n_lon,n_lat)
    !     call read_lonlat_from_binary(arr(1,1),n_lon*n_lat,fid)
    implicit none
    integer               :: n, fid
    real(8), dimension(n) :: arr

    integer :: i

    read(fid) (arr(i),i=1,n)
  end subroutine read_lonlat_from_binary

  subroutine read_HR_optim_grid
    ! Reads in Heikes & Randall (1995) optimized grid from file in directory grid_HR
    ! Need to provide a symbolic link to grid_HR in working directory
    implicit none
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

    l2error = sqrt (sum_real (l2error))
    maxerror = sync_max_real (maxerror)

    if (rank == 0) then
       write (6,'(A)') '-------------------------------------------------------&
            --------------------------------------------------------------------------'
       write (6,'(A,i2,A,/)') 'Heikes-Randall optimizations of level ', level_start-1, ' grid:'
       write (6,'(A,2(es8.2,A))') 'Grid quality before optimization = ', maxerror, ' m (linf) ', l2error, ' m (l2)'
    end if
    
    fid = get_fid()
    if (level_start /= level_end) then
       write (0,'(i2,1x,i2)') level_end, level_start
       write (0,'(A)') "Reading HR grid points for level_start not equal to level_end not implemented"
       return
    end if

    write (filename, '(A,I1)')  "grid_HR/J", level_start-1
    open (unit=fid, file=filename)

    p = 1
    do d_HR = 1, N_ICOSAH_LOZENGE
       loz = dom_id_from_HR_id(d_HR)
       do d_sub = 1, N_SUB_DOM
          d_glo = loz*N_SUB_DOM + sub_dom_id_from_HR_sub_id(d_sub)
          if (owner(d_glo+1) == rank) call get_offs_Domain (grid(loc_id(d_glo+1)+1), p, offs, dims)
          call coord_from_file (d_glo, PATCH_LEVEL, fid, offs, dims, (/0, 0/))
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

    l2error = sqrt (sum_real(l2error))
    maxerror = sync_max_real (maxerror)
    if (rank == 0) then
       write (6,'(A,2(es8.2,A))') 'Grid quality after optimization  = ', maxerror, ' m (linf) ', l2error, ' m (l2)'
       write (6,'(A)') '(distance between midpoints of primal and dual edges)'
       write (6,'(A,/)') '-------------------------------------------------------&
            --------------------------------------------------------------------------'
    end if
  end subroutine read_HR_optim_grid

  integer function dom_id_from_HR_id (d_HR)
    ! d_HR: lozenge id as used by Heikes & Randall (starts from 1)
    ! results: domain id (starts from 0)
    implicit none
    integer :: d_HR

    dom_id_from_HR_id = modulo(d_HR,2)*5 + modulo(d_HR/2-1,5)
  end function dom_id_from_HR_id

  integer function sub_dom_id_from_HR_sub_id (sub_id)
    ! sub_id: lozenge sub id as used by Heikes & Randall (starts from 1)
    ! results: sub domain id (starts from 0)
    implicit none
    integer :: sub_id

    integer :: id, i, j, halv_sub_dom, l, jdiv, idiv

    i = 0
    j = 0
    id = sub_id - 1
    halv_sub_dom = N_SUB_DOM/2
    do l = DOMAIN_LEVEL-1, 0, -1
       jdiv = id/halv_sub_dom
       j = j + jdiv*2**l
       id = modulo (id+4**l,4**(l+1))
       idiv = id/halv_sub_dom
       i = i + idiv*2**l
       halv_sub_dom = halv_sub_dom/4
       id = modulo (id,4**l)
    end do
    sub_dom_id_from_HR_sub_id = j*N_SUB_DOM_PER_DIM + i
  end function sub_dom_id_from_HR_sub_id

  subroutine zrotate (c_in, c_out, angle)
    implicit none
    real(8),      intent(in) :: angle
    type(Coord),  intent(in) :: c_in
    type(Coord), intent(out) :: c_out

    c_out%x =  c_in%x*cos(angle) - c_in%y*sin(angle)
    c_out%y =  c_in%x*sin(angle) + c_in%y*cos(angle)
    c_out%z =  c_in%z
  end subroutine zrotate

  recursive subroutine coord_from_file (d_glo, l, fid, offs, dims, ij0)
    implicit none
    integer,                        intent(in) :: d_glo, l, fid
    integer, dimension(2),          intent(in) :: ij0
    integer, dimension(N_BDRY+1),   intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims

    integer :: d_loc, k, ij(2)
    type(Coord) node, node_r

    d_loc = loc_id(d_glo+1)
    do k = 1, 4
       ij = ij0 + HR_offs(:,k)*2**(l-1)
       if (l == 1) then
          ! if domain on other process still read to get to correct position in file
          if (owner(d_glo+1) == rank) then
             read(fid,*) node
             call zrotate (node, node_r, -0.5_8)  ! icosahedron orientation good for tsunami
             grid(d_loc+1)%node%elts(idx(ij(1), ij(2), offs, dims)+1) = project_on_sphere(node_r)
          else
             read(fid,*)
          end if
       else
          call coord_from_file (d_glo, l-1, fid, offs, dims, ij)
       end if
    end do
  end subroutine coord_from_file

  subroutine pre_levelout
    implicit none
    integer :: d, l, max_output_level, num

    ! FIXME cleaner would be to use init_io routine
    call init_Float_Field (active_level, AT_NODE)

    do d = 1, size(grid)
       num = grid(d)%node%length
       call init (active_level%data(d), num)
       active_level%data(d)%elts(1:num) = grid(d)%level%elts(1:num)
    end do

    do l = level_end-1, level_start, -1
       call apply_interscale (restrict_level, l, z_null, 0, 1)
    end do
  end subroutine pre_levelout

  ! now active_level can be used

  subroutine post_levelout
    implicit none
    integer :: d

    do d = 1, size(grid)
       deallocate (active_level%data(d)%elts)
    end do
    deallocate (active_level%data)
  end subroutine post_levelout

  subroutine restrict_level (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: d, id_par, id_chd

    d = dom%id+1

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx(i_par, j_par, offs_par, dims_par)

    if (dom%mask_n%elts(id_chd+1) >= ADJZONE) active_level%data(d)%elts(id_par+1) = active_level%data(d)%elts(id_chd+1)
  end subroutine restrict_level

  subroutine write_and_export (iwrite)
    implicit none
    integer :: iwrite

    integer      :: d, i, j, k, l, p, u
    character(7) :: var_file

    if (rank == 0) write(6,'(/,A,i4/)') 'Saving fields ', iwrite

    sol%bdry_uptodate = .false.
    call update_array_bdry (sol, NONE)

    call pre_levelout

    ! Compute surface pressure
    call cal_surf_press (sol)

    do l = level_start, level_end
       minv = 1.0d63; maxv = -1.0d63
       u = 1000000+100*iwrite

       ! Calculate pressure, exner and geopotential at vertical level save_zlev and scale l
       do k = 1, save_zlev
          do d = 1, size(grid)
             mass  => sol(S_MASS,k)%data(d)%elts
             temp  => sol(S_TEMP,k)%data(d)%elts
             exner => exner_fun(k)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do
             nullify (mass, temp, exner)
          end do
       end do

       ! Calculate zonal and meridional velocities and vorticity for vertical level save_zlev
       do d = 1, size(grid)
          velo => sol(S_VELO,save_zlev)%data(d)%elts
          vort => grid(d)%vort%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,       grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), l, z_null)
          nullify (velo, vort)
       end do

       ! Calculate vorticity at hexagon points (stored in press_lower)
       call apply_onescale (vort_triag_to_hex, l, z_null, 0, 1)

       call write_level_mpi (write_primal, u+l, l, save_zlev, .true., run_id)

       do i = 1, N_VAR_OUT
          minv(i) = -sync_max_real (-minv(i))
          maxv(i) =  sync_max_real ( maxv(i))
       end do
       if (rank == 0) then
          write (var_file, '(i7)') u
          open(unit=50, file=trim(run_id)//'.'//var_file)
          write(50,'(A, 7(E15.5E2, 1X), I3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
          write(50,'(A, 7(E15.5E2, 1X), I3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
          close(50)
       end if
       u = 2000000+100*iwrite
       call write_level_mpi (write_dual, u+l, l, save_zlev, .False., run_id)
    end do

    call post_levelout
    call barrier
    if (rank == 0) call compress_files (iwrite, run_id)
    call barrier
  end subroutine write_and_export

  subroutine compress_files (iwrite, run_id)
    implicit none
    integer      :: iwrite
    character(*) :: run_id

    character(4)   :: s_time
    character(130) :: command

    write (s_time, '(i4.4)') iwrite

    command = 'ls -1 '//trim(run_id)//'.1'//s_time//'?? > tmp1'
    
    call system (command)

    command = 'tar cfz '//trim(run_id)//'.1'//s_time//'.tgz -T tmp1 --remove-files'
    call system (command)

    command = 'ls -1 '//trim(run_id)//'.2'//s_time //'?? > tmp2' 
    call system (command)

    command = 'tar cfz '//trim(run_id)//'.2'//s_time//'.tgz -T tmp2 --remove-files'
    call system (command)
  end subroutine compress_files
  
  subroutine zonal_meridional_vel (dom, i, j, zlev, offs, dims)
    ! Finds lat-lon velocity (with components in zonal and meridional directions) given index information of node
    ! using lapack least squares routine dgels
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    real(8), dimension (3)         :: uvw
    real(8), dimension (2)         :: vel_latlon

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

    id = idx (i, j, offs, dims)

    id_i = id + 1
    id_e = EDGE*id + 1
    idN  = idx (i,   j+1, offs, dims) + 1
    idE  = idx (i+1, j,   offs, dims) + 1
    idNE = idx (i+1, j+1, offs, dims) + 1

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

    e_zonal = Coord (-sin(lon),           cos(lon),             0.0_8) ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    ! Least squares overdetermined matrix 
    do e = 1, EDGE
       A(e,1) = inner (dir(e), e_zonal)
       A(e,2) = inner (dir(e), e_merid)
    end do

    ! Solve least squares problem to find zonal and meridional velocities
    call dgels ('N', 3, 2, 1, A, 3, uvw, 3, work, lwork, info)

    dom%u_zonal%elts(id_i) = uvw(1)
    dom%v_merid%elts(id_i) = uvw(2)
  end subroutine zonal_meridional_vel
end module io_mod
