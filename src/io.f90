module io_mod
  use geom_mod
  use domain_mod
  use arch_mod
  use adapt_mod
  use smooth_mod
  use comm_mpi_mod
  use multi_level_mod
  use remap_mod
  implicit none

  integer, parameter                   :: N_VAR_OUT = 7
  integer, dimension(2,4)              :: HR_offs
  real(8)                              :: dx_export, dy_export, kx_export, ky_export, vmin, vmax
  real(8), dimension(N_VAR_OUT)        :: minv, maxv
  real(8), pointer                     :: press_save
  real(4), dimension(:,:), allocatable :: field2d
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
       vmin = min(vmin, vort)
       vmax = max(vmax, vort)

    end if
  end subroutine vort_extrema

  subroutine write_step (fid, time, k)
    implicit none
    integer ::  fid, k
    real(8) :: time

    integer :: l
    real(8) :: tot_mass

    vmin =  1d-16
    vmax = -1d-16

    do l = level_start, level_end
       call apply_onescale(vort_extrema, l, z_null, 0, 0)
    end do

    tot_mass = integrate_hex (mu, level_start, k)

    if (rank == 0) write(fid,'(E16.9, I3, 2(1X, I9), 7(1X, E16.8), 1X, F16.7)') &
         time, level_end, n_active, tot_mass, get_timing()
  end subroutine write_step

  real(8) function integrate_hex (fun, l, k)
    ! Integrate function defined by fun over hexagons
    implicit none
    external :: fun
    integer  :: l, k

    integer                        :: d, ll, p, i, j, c, id
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: s, fun

    s = 0.0
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
          do while (grid(d)%patch%elts(p+1)%level < l)
             p = grid(d)%patch%elts(p+1)%children(c-4)
             if (p == 0) then
                write(6,*) "ERROR(rank=", rank, "):integrate_hex: level incomplete"
                return
             end if
          end do
          call get_offs_Domain (grid(d), p, offs, dims)
          if (c == NORTHWEST) then
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

  real(8) function integrate_tri (fun, k)
    implicit none
    external :: fun
    integer  :: k

    integer                        :: d, ll, p, i, j, t, id
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: s, fun

    s = 0.0
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

  real(8) function only_area (dom, i, j, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    only_area = 1.0_8
  end function only_area

  real(8) function mu (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx(i, j, offs, dims)
    mu = sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)
  end function mu

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

  real(8) function tri_only_area (dom, i, j, t, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, t
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

  subroutine export_2d (proj, fid, Nx, Ny, lon_lat_range, test_case)
    ! Interpolate variables defined in valrange onto lon-lat grid of size (Nx(1):Nx(2), Ny(1):Ny(2), zlevels),
    ! save zonal average and horizontal grid at vertical level zlevel
    implicit none
    external               :: proj
    integer, dimension(2)  :: Nx, Ny
    integer                :: fid
    real(8), dimension(2)  :: lon_lat_range
    character(*)           :: test_case

    integer                              :: d, i, id, ix, j, k, level_end_old, p, v
    integer, parameter                   :: funit = 400
    integer, parameter                   :: nvar_save=7, nvar_zonal=8 ! Number of variables to save

    real(8)                              :: N_zonal
    real(8), parameter                   :: TT=200.0_8 ! shift for stable variance calculation
    real, dimension(:,:,:), allocatable  :: field2d_save, zonal_av
    real, dimension(:,:),   allocatable  :: uprime, vprime, Tprime

    character(5)                         :: s_time
    character(7)                         :: var_file
    character(130)                       :: command

    N_zonal = real(Nx(2)-Nx(1)+1)

    ! Fill up grid to level l and do inverse wavelet transform onto the uniform grid at level l
    call fill_up_grid_and_IWT (level_save)

    call cal_surf_press (sol)

    ! Remap to pressure_save vertical levels for saving data
    sol_save = sol(:,1:save_levels)
    call apply_onescale (interp_save, level_save, z_null, -1, 2)

    ! Calculate temperature at all vertical levels (saved in exner_fun) and temperature at interpolated saved vertical levels
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    call update_vector_bdry (exner_fun, NONE)

    dx_export = lon_lat_range(1)/(Nx(2)-Nx(1)+1); dy_export = lon_lat_range(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1.0_8/dx_export; ky_export = 1.0_8/dy_export
    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))
    allocate (uprime(Nx(1):Nx(2),Ny(1):Ny(2)), vprime(Nx(1):Nx(2),Ny(1):Ny(2)), Tprime(Nx(1):Nx(2),Ny(1):Ny(2)))
    allocate (field2d_save(Nx(1):Nx(2),Ny(1):Ny(2),nvar_save*save_levels))
    allocate (zonal_av(1:zlevels,Ny(1):Ny(2),nvar_zonal))

    do k = 1, save_levels
       ! Mass density
       call project_onto_plane (sol_save(S_MASS,k), Nx, Ny, level_save, proj, 0.0_8)
       field2d_save(:,:,1+k-1) = field2d

       ! Temperature
       call project_onto_plane (horiz_flux(k), Nx, Ny, level_save, proj, 0.0_8)
       field2d_save(:,:,2+k-1) = field2d

       ! Calculate zonal and meridional velocities and vorticity
       do d = 1, size(grid)
          velo => sol_save(S_VELO,k)%data(d)%elts
          vort => grid(d)%vort%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(level_save)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,       grid(d), grid(d)%lev(level_save)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), level_save, z_null)
          nullify (velo, vort)
       end do

       ! Zonal velocity
       call project_uzonal_onto_plane (Nx, Ny, level_save, proj, 0.0_8)
       field2d_save(:,:,3+k-1) = field2d

       ! Meridional velocity
       call project_vmerid_onto_plane (Nx, Ny, level_save, proj, 0.0_8)
       field2d_save(:,:,4+k-1) = field2d

       ! Geopotential
       press_save => pressure_save(k)
       call apply_onescale (cal_geopot, level_save, z_null, 0, 1)
       call project_geopot_onto_plane (Nx, Ny, level_save, proj, 1.0_8)
       field2d_save(:,:,5+k-1) = field2d

       ! Vorticity
       call apply_onescale (vort_triag_to_hex, level_save, z_null, 0, 1)
       call project_vorticity_onto_plane (Nx, Ny, level_save, proj, 1.0_8)
       field2d_save(:,:,6+k-1) = field2d

       ! Surface pressure
       call project_surf_press_onto_plane (Nx, Ny, level_save, proj, 1.0_8)
       field2d_save(:,:,7+k-1) = field2d
    end do

    ! Zonal averages
    do k = 1, zlevels
       ! Mass density
       call project_onto_plane (sol(S_MASS,k), Nx, Ny, level_save, proj, 0.0_8)
       zonal_av(k,:,1) = sum(field2d,DIM=1)/N_zonal

       ! Temperature
       call project_onto_plane (exner_fun(k), Nx, Ny, level_save, proj, 1.0_8)
       zonal_av(k,:,2) = sum(field2d,DIM=1)/N_zonal

       ! Peturbation Temperature
       do ix = Nx(1), Nx(2)
          Tprime(ix,:) = field2d(ix,:) - zonal_av(k,:,2)
       end do

       ! Variance of temperature (stable calculation)
       zonal_av(k,:,3) = (sum((field2d-TT)**2,DIM=1) - sum(field2d-TT,DIM=1)**2/N_zonal)/(N_zonal-1)

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo => sol(S_VELO,k)%data(d)%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(level_save)%elts(j), k, 0, 1)
          end do
          nullify (velo)
       end do

       ! Zonal velocity
       call project_uzonal_onto_plane (Nx, Ny, level_save, proj, 0.0_8)
       zonal_av(k,:,4) = sum(field2d,DIM=1)/N_zonal

       ! Peturbation zonal velocity
       do ix = Nx(1), Nx(2)
          uprime(ix,:) = field2d(ix,:) - zonal_av(k,:,4)
       end do

       ! Meridional velocity
       call project_vmerid_onto_plane (Nx, Ny, level_save, proj, 0.0_8)
       zonal_av(k,:,5) = sum(field2d,DIM=1)/N_zonal

       ! Peturbation meridional velocity
       do ix = Nx(1), Nx(2)
          vprime(ix,:) = field2d(ix,:) - zonal_av(k,:,5)
       end do

       ! Eddy momentum flux
       zonal_av(k,:,6) = sum(uprime*vprime,DIM=1)/N_zonal

       ! Eddy kinetic energy
       zonal_av(k,:,7) = sum(0.5_8*(uprime*uprime+vprime*vprime),DIM=1)/N_zonal

       ! Eddy heat flux
       zonal_av(k,:,8) = sum(Tprime*vprime,DIM=1)/N_zonal
    end do

    if (rank == 0) then
       ! Solution at level zlev
       do v = 1, nvar_save*save_levels
          write (var_file, '(i7)') fid+v
          open (unit=funit, file=trim(test_case)//'.'//var_file, recl=32768)
          do i = Ny(1), Ny(2)
             write (funit,'(2047(E15.6, 1X))') field2d_save(:,i,v)
          end do
          close (funit)
       end do

       ! Zonal average of solution over all vertical levels
       do v = 1, nvar_zonal
          write (var_file, '(i7)') fid+v+10
          open (unit=funit, file=trim(test_case)//'.'//var_file, recl=32768)
          do k = zlevels,1,-1
             write (funit,'(2047(E15.6, 1X))') zonal_av(k,:,v)
          end do
          close (funit)
       end do

       ! Coordinates

       ! Longitude values
       write (var_file, '(i7)') fid+20
       open (unit=funit, file=trim(test_case)//'.'//var_file, recl=32768) 
       write (funit,'(2047(E15.6, 1X))') (-180.0_8+dx_export*(i-1)/MATH_PI*180.0_8, i=1,Nx(2)-Nx(1)+1)
       close (funit)

       ! Latitude values
       write (var_file, '(i7)') fid+21
       open (unit=funit, file=trim(test_case)//'.'//var_file, recl=32768) 
       write (funit,'(2047(E15.6, 1X))') (-90.0_8+dy_export*(i-1)/MATH_PI*180.0_8, i=1,Ny(2)-Ny(1)+1)
       close (funit)

       ! Pressure vertical coordinates
       write (var_file, '(i7)') fid+22
       open (unit=funit, file=trim(test_case)//'.'//var_file, recl=32768) 
       write (funit,'(2047(E15.6, 1X))') (a_vert(k)*ref_press/ref_surf_press + b_vert(k), k=zlevels+1,2,-1)
       close (funit)

       write (s_time, '(i5)') fid/100
       command = 'ls -1 '//trim(test_case)//'.'//s_time//'?? > tmp' 
       call system (command)

       command = 'tar czf '//trim(test_case)//'.'//s_time//'.tgz -T tmp --remove-files &'
       call system (command)
    end if
    deallocate (field2d, field2d_save, zonal_av)
  end subroutine export_2d

  subroutine project_onto_plane (field, Nx, Ny, l, proj, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    implicit none
    external              :: proj
    integer               :: l, itype
    integer, dimension(2) :: Nx, Ny
    real(8)               :: default_val
    Type(Float_field)     :: field

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

                call proj (grid(d)%node%elts(id+1),   cC)
                call proj (grid(d)%node%elts(idN+1),  cN)
                call proj (grid(d)%node%elts(idE+1),  cE)
                call proj (grid(d)%node%elts(idNE+1), cNE)

                val   = field%data(d)%elts(id+1)
                valN  = field%data(d)%elts(idN+1)
                valE  = field%data(d)%elts(idE+1)
                valNE = field%data(d)%elts(idNE+1)

                if (abs(cN(2) - MATH_PI/2.0_8) < sqrt(1.0d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs(cE(2) + MATH_PI/2.0_8) < sqrt(1.0d-15)) then
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
  end subroutine project_onto_plane

  subroutine project_geopot_onto_plane (Nx, Ny, l, proj, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    implicit none
    external              :: proj
    integer               :: l, itype
    integer, dimension(2) :: Nx, Ny
    real(8)               :: default_val

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

                call proj (grid(d)%node%elts(id+1),   cC)
                call proj (grid(d)%node%elts(idN+1),  cN)
                call proj (grid(d)%node%elts(idE+1),  cE)
                call proj (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%adj_geopot%elts(id+1)
                valN  = grid(d)%adj_geopot%elts(idN+1)
                valE  = grid(d)%adj_geopot%elts(idE+1)
                valNE = grid(d)%adj_geopot%elts(idNE+1)

                if (abs(cN(2) - MATH_PI/2.0_8) < sqrt(1.0d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs(cE(2) + MATH_PI/2.0_8) < sqrt(1.0d-15)) then
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
  end subroutine project_geopot_onto_plane

  subroutine project_surf_press_onto_plane (Nx, Ny, l, proj, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    implicit none
    external              :: proj
    integer               :: l, itype
    integer, dimension(2) :: Nx, Ny
    real(8)               :: default_val

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

                call proj (grid(d)%node%elts(id+1),   cC)
                call proj (grid(d)%node%elts(idN+1),  cN)
                call proj (grid(d)%node%elts(idE+1),  cE)
                call proj (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%surf_press%elts(id+1)/1.0d2
                valN  = grid(d)%surf_press%elts(idN+1)/1.0d2
                valE  = grid(d)%surf_press%elts(idE+1)/1.0d2
                valNE = grid(d)%surf_press%elts(idNE+1)/1.0d2

                if (abs(cN(2) - MATH_PI/2.0_8) < sqrt(1.0d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs(cE(2) + MATH_PI/2.0_8) < sqrt(1.0d-15)) then
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
  end subroutine project_surf_press_onto_plane

  subroutine project_vorticity_onto_plane (Nx, Ny, l, proj, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    implicit none
    external              :: proj
    integer               :: l, itype
    integer, dimension(2) :: Nx, Ny
    real(8)               :: default_val

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

                call proj (grid(d)%node%elts(id+1),   cC)
                call proj (grid(d)%node%elts(idN+1),  cN)
                call proj (grid(d)%node%elts(idE+1),  cE)
                call proj (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%adj_mass%elts(id+1)
                valN  = grid(d)%adj_mass%elts(idN+1)
                valE  = grid(d)%adj_mass%elts(idE+1)
                valNE = grid(d)%adj_mass%elts(idNE+1)

                if (abs(cN(2) - MATH_PI/2.0_8) < sqrt(1.0d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs(cE(2) + MATH_PI/2.0_8) < sqrt(1.0d-15)) then
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
  end subroutine project_vorticity_onto_plane

  subroutine project_uzonal_onto_plane (Nx, Ny, l, proj, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    implicit none
    external              :: proj
    integer, dimension(2) :: Nx, Ny
    integer               :: l, itype
    real(8)               :: default_val

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

                call proj (grid(d)%node%elts(id+1),   cC)
                call proj (grid(d)%node%elts(idN+1),  cN)
                call proj (grid(d)%node%elts(idE+1),  cE)
                call proj (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%u_zonal%elts(id+1)
                valN  = grid(d)%u_zonal%elts(idN+1)
                valE  = grid(d)%u_zonal%elts(idE+1)
                valNE = grid(d)%u_zonal%elts(idNE+1)

                if (abs(cN(2) - MATH_PI/2.0_8) < sqrt(1.0d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs(cE(2) + MATH_PI/2.0_8) < sqrt(1.0d-15)) then
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
  end subroutine project_uzonal_onto_plane

  subroutine project_vmerid_onto_plane (Nx, Ny, l, proj, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    implicit none
    external              :: proj
    integer, dimension(2) :: Nx, Ny
    integer               :: l, itype
    real(8)               :: default_val

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

                call proj (grid(d)%node%elts(id+1),   cC)
                call proj (grid(d)%node%elts(idN+1),  cN)
                call proj (grid(d)%node%elts(idE+1),  cE)
                call proj (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%v_merid%elts(id+1)
                valN  = grid(d)%v_merid%elts(idN+1)
                valE  = grid(d)%v_merid%elts(idE+1)
                valNE = grid(d)%v_merid%elts(idNE+1)

                if (abs(cN(2) - MATH_PI/2.0_8) < sqrt(1.0d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs(cE(2) + MATH_PI/2.0_8) < sqrt(1.0d-15)) then
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
  end subroutine project_vmerid_onto_plane

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

    if (sum(abs(fixed)) > 1) write(0,*) 'ALARM'

    if (sum(fixed) /= 0) then
       a(1) = a(1) - sum(fixed)*MATH_PI*2.0_8
       b(1) = b(1) - sum(fixed)*MATH_PI*2.0_8
       c(1) = c(1) - sum(fixed)*MATH_PI*2.0_8
       call interp_tri_to_2d (a, b, c, val)
    end if
  end subroutine interp_tri_to_2d_and_fix_bdry

  subroutine interp_tri_to_2d (a, b, c, val)
    real(8), dimension(2) :: a, b, c
    real(8), dimension(3) :: val

    integer               :: id_x, id_y
    real(8)               :: ival, minx, maxx, miny, maxy
    real(8), dimension(2) :: ll
    real(8), dimension(3) :: bac
    logical               :: inside

    minx = min(min(a(1), b(1)), c(1))
    maxx = max(max(a(1), b(1)), c(1))
    miny = min(min(a(2), b(2)), c(2))
    maxy = max(max(a(2), b(2)), c(2))
    if (maxx-minx > MATH_PI/2.0_8) then
       write(0,*) 'ERROR(rank', rank, '):io-333 "export"'
       return
    end if

    do id_x = floor(kx_export*minx), ceiling(kx_export*maxx)
       if (id_x < lbound(field2d,1) .or. id_x > ubound(field2d,1)) cycle
       do id_y = floor(ky_export*miny), ceiling(ky_export*maxy)
          if (id_y < lbound(field2d,2) .or. id_y > ubound(field2d,2)) cycle
          ll = (/dx_export*id_x, dy_export*id_y/)
          call interp_tria (ll, a, b, c, val, ival, inside)
          if (inside) field2d(id_x,id_y) = ival
       end do
    end do
  end subroutine interp_tri_to_2d

  subroutine interp_tria (ll, coord1, coord2, coord3, values, ival, inside)
    implicit none
    real(8), dimension(2) :: coord1, coord2, coord3
    real(8), dimension(3) :: values
    real(8)               :: ival
    logical               :: inside

    real(8), dimension(2) :: ll
    real(8), dimension(3) :: bc

    bc = bary_coord(ll, coord1, coord2, coord3)
    inside = (0.0 < bc(1) .and. bc(1) < 1.0_8 .and. &
         0.0 < bc(2) .and. bc(2) < 1.0_8 .and. &
         0.0 < bc(3) .and. bc(3) < 1.0_8)
    if (inside) ival = sum(values*bc)
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
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, d, k
    real(8), dimension(zlevels) :: pressure

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    ! Integrate the pressure upwards
    pressure(1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*sol(S_MASS,1)%data(d)%elts(id+1)
    do k = 2, zlevels
       pressure(k) = pressure(k-1) - grav_accel*interp(sol(S_MASS,k)%data(d)%elts(id+1), sol(S_MASS,k-1)%data(d)%elts(id+1))
    end do

    ! Temperature at all vertical levels (saved in exner_fun)
    do k = 1, zlevels
       exner_fun(k)%data(d)%elts(id+1) = sol(S_TEMP,k)%data(d)%elts(id+1)/sol(S_MASS,k)%data(d)%elts(id+1) &
            * (pressure(k)/ref_press)**kappa
    end do

    ! temperature at save levels (saved in horiz_flux)
    do k = 1, save_levels
       horiz_flux(k)%data(d)%elts(id+1) = sol_save(S_TEMP,k)%data(d)%elts(id+1)/sol_save(S_MASS,k)%data(d)%elts(id+1) * &
            (pressure_save(k)/ref_press)**kappa
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
    dom%adj_geopot%elts(id+1) = dom%surf_geopot%elts(id+1)/grav_accel

    k = 1
    do while (pressure_upper > press_save)
       dom%adj_geopot%elts(id+1) = dom%adj_geopot%elts(id+1) + &
            R_d/grav_accel * exner_fun(k)%data(d)%elts(id+1) * (log(pressure_lower)-log(pressure_upper))

       k = k+1
       pressure_lower = pressure_upper
       pressure_upper = pressure_lower - grav_accel*sol(S_MASS,k+1)%data(d)%elts(id+1)
    end do

    ! Add additional contribution up to pressure level pressure_save
    dom%adj_geopot%elts(id+1) = dom%adj_geopot%elts(id+1) &
         + R_d/grav_accel * exner_fun(k)%data(d)%elts(id+1) * (log(pressure_lower)-log(press_save))
  end subroutine cal_geopot

  subroutine interp_save (dom, i, j, zlev, offs, dims)
    ! Linear interpolation to save levels
    ! Assumes variables have been remapped to original vertical grid                                                                                                 
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, d, e, k, kk
    real(8) :: dpressure, pressure_lower, pressure_upper, spress

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    spress = 0.0_8
    do k = 1, zlevels
       spress = spress + sol(S_MASS,k)%data(d)%elts(id+1)
    end do
    spress = press_infty + grav_accel*spress

    do kk = 1, save_levels
       ! Find pressure at current levels (not interfaces)                                                                                                            
       pressure_lower = spress
       pressure_upper = 0.5_8*(a_vert(1)+a_vert(2))*ref_press + 0.5_8*(b_vert(1)+b_vert(2))*spress
       k = 1
       do while (pressure_upper > pressure_save(kk))
          k = k+1
          pressure_lower = pressure_upper
          pressure_upper = 0.5_8*(a_vert(k)+a_vert(k+1))*ref_press + 0.5_8*(b_vert(k)+b_vert(k+1))*spress
       end do
       if (k==1) return ! Skip incorrect values for pentagons
       dpressure =  (pressure_save(kk)-pressure_upper)/(pressure_lower-pressure_upper)

       ! Linear interpolation                                                                                                                                        
       sol_save(S_MASS,kk)%data(d)%elts(id+1) = sol(S_MASS,k+1)%data(d)%elts(id+1) + &
            dpressure * (sol(S_MASS,k)%data(d)%elts(id+1) - sol(S_MASS,k+1)%data(d)%elts(id+1))
       sol_save(S_TEMP,kk)%data(d)%elts(id+1) = sol(S_TEMP,k+1)%data(d)%elts(id+1) + &
            dpressure * (sol(S_TEMP,k)%data(d)%elts(id+1) - sol(S_TEMP,k+1)%data(d)%elts(id+1))
       do e = 1, EDGE
          sol_save(S_VELO,kk)%data(d)%elts(EDGE*id+e) = sol(S_VELO,k+1)%data(d)%elts(EDGE*id+e)  + &
               dpressure * (sol(S_VELO,k)%data(d)%elts(EDGE*id+e) - sol(S_VELO,k+1)%data(d)%elts(EDGE*id+e))
       end do
    end do
  end subroutine interp_save

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

    dom%adj_mass%elts(id+1) = ( &
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
    outv(1) = sol(S_TEMP,zlev)%data(d)%elts(id+1)/sol(S_MASS,zlev)%data(d)%elts(id+1)*(dom%press%elts(id+1)/ref_press)**kappa

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
    outv(7) =  dom%adj_mass%elts(id+1)
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

  subroutine zonal_meridional_vel (dom, i, j, offs, dims, zlev, vel_latlon)
    ! Finds lat-lon velocity (with components in zonal and meridional directions) given index information of node
    ! using lapack least squares routine dgels
    implicit none
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

    e_zonal = Coord (-sin(lon),           cos(lon),             0.0_8) ! Zonal direction
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
          write(fid,*) wav_coeff(S_VELO,k)%data(dom%id+1)%elts(EDGE*id+e)
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
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, fid, v

    d = dom%id+1
    id = idx(i, j, offs, dims)

    do v = S_MASS, S_TEMP
       write(fid) sol(v,zlev)%data(d)%elts(id+1) ! for pole
       write(fid) trend(v,zlev)%data(d)%elts(id+1) ! for pole`
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
       read(fid) sol(v,zlev)%data(d)%elts(id+1) ! for pole
       read(fid) trend(v,zlev)%data(d)%elts(id+1) ! for pole
    end do
  end subroutine read_scalar

  subroutine dump_adapt_mpi (id, custom_dump)
    ! Save data in check point files for restart
    ! One file per domain
    implicit none
    external :: custom_dump
    integer  :: id

    character(255)                          :: filename_gr, filename_no
    integer                                 :: c, d, fid_gr, fid_no, i, j, k, l, p_chd, p_lev, p_par, v
    logical, dimension(N_CHDRN)             :: child_required
    type(Domain), dimension(:), allocatable :: grid_tmp

    allocate (grid_tmp, source=grid)
    
    fid_no = id+1000000
    fid_gr = id+3000000

    call update_array_bdry (wav_coeff(S_MASS:S_TEMP,:),       NONE)
    call update_array_bdry (trend_wav_coeff(S_MASS:S_TEMP,:), NONE)

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
       write (filename_no, '(A,I4.4,A,I5.5)') "coef.", id, "_", glo_id(rank+1,d)
       write (filename_gr, '(A,I4.4,A,I5.5)') "grid.", id, "_", glo_id(rank+1,d)

       open (unit=fid_no, file=trim(filename_no), form="UNFORMATTED", action='WRITE')
       open (unit=fid_gr, file=trim(filename_gr), form="UNFORMATTED", action='WRITE')

       write (fid_no) istep
       write (fid_no) time

       call custom_dump (fid_no)

       do k = 1, zlevels
          call apply_to_pole_d (write_scalar, grid(d), min_level-1, k, fid_no, .True.)
          do v = S_MASS, S_VELO
             write (fid_no) (sol(v,k)%data(d)%elts(i), i=MULT(v)*grid(d)%patch%elts(1+1)%elts_start+1, &
                  MULT(v)*(grid(d)%patch%elts(1+1)%elts_start+PATCH_SIZE**2))
             write (fid_no) (trend(v,k)%data(d)%elts(i), i=MULT(v)*grid(d)%patch%elts(1+1)%elts_start+1, &
                  MULT(v)*(grid(d)%patch%elts(1+1)%elts_start+PATCH_SIZE**2))
          end do
       end do

       do l = min_level, level_end
          p_lev = 0
          do j = 1, grid_tmp(d)%lev(l)%length
             p_par = grid_tmp(d)%lev(l)%elts(j)
             if (grid_tmp(d)%patch%elts(p_par+1)%deleted) then
                do c = 1, N_CHDRN
                   p_chd = grid_tmp(d)%patch%elts(p_par+1)%children(c)
                   if (p_chd > 0) grid_tmp(d)%patch%elts(p_chd+1)%deleted = .True.
                end do
                cycle
             end if

             do k = 1, zlevels
                do v = S_MASS, S_VELO
                   write (fid_no) (wav_coeff(v,k)%data(d)%elts(i), i = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start+1, &
                        MULT(v)*(grid(d)%patch%elts(p_par+1)%elts_start+PATCH_SIZE**2))
                   write (fid_no) (trend_wav_coeff(v,k)%data(d)%elts(i), i = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start+1, &
                        MULT(v)*(grid(d)%patch%elts(p_par+1)%elts_start+PATCH_SIZE**2))
                end do
             end do

              do c = 1, N_CHDRN
                p_chd = grid_tmp(d)%patch%elts(p_par+1)%children(c)
                if (p_chd > 0) then
                   child_required(c) = check_child_required(grid_tmp(d), p_par, c-1)
                   grid_tmp(d)%patch%elts(p_chd+1)%deleted = .not. child_required(c)
                   if (child_required(c)) then
                      p_lev = p_lev + 1
                      grid_tmp(d)%lev(l+1)%elts(p_lev) = p_chd
                   end if
                else
                   child_required(c) = .False.
                end if
             end do

             write (fid_gr) child_required
          end do
          if (l+1 <= max_level) grid_tmp(d)%lev(l+1)%length = p_lev
       end do
       close (fid_no); close (fid_gr)
    end do
    deallocate (grid_tmp)
  end subroutine dump_adapt_mpi

  subroutine load_adapt_mpi (id, custom_load)
    ! Read data from check point files for restart
    ! One file per domain
    implicit none
    external                             :: custom_load
    integer                              :: id

    character(255)                       :: filename_gr, filename_no
    integer                              :: c, d, i, j, k, l, old_n_patch, p_chd, p_par, v
    integer, dimension(n_domain(rank+1)) :: fid_no, fid_gr
    logical, dimension(N_CHDRN)          :: child_required

    ! Loop over domains at coarsest scale
    do d = 1, size(grid)
       fid_no(d) = id*1000 + 1000000 + d
       fid_gr(d) = id*1000 + 3000000 + d

       write (filename_no, '(A,I4.4,A,I5.5)') "coef.", id, "_", glo_id(rank+1,d)
       write (filename_gr, '(A,I4.4,A,I5.5)') "grid.", id, "_", glo_id(rank+1,d)

       open (unit=fid_no(d), file=trim(filename_no), form="UNFORMATTED", action='READ')
       open (unit=fid_gr(d), file=trim(filename_gr), form="UNFORMATTED", action='READ')

       read (fid_no(d)) istep
       read (fid_no(d)) time

       call custom_load (fid_no(d))

       do k = 1, zlevels
          call apply_to_pole_d (read_scalar, grid(d), min_level-1, k, fid_no(d), .true.)

          do v = S_MASS, S_VELO
             read (fid_no(d)) (sol(v,k)%data(d)%elts(i),i = MULT(v)* grid(d)%patch%elts(1+1)%elts_start+1, &
                  MULT(v)*(grid(d)%patch%elts(1+1)%elts_start+PATCH_SIZE**2) )
             read (fid_no(d)) (trend(v,k)%data(d)%elts(i),i = MULT(v)* grid(d)%patch%elts(1+1)%elts_start+1, &
                  MULT(v)*(grid(d)%patch%elts(1+1)%elts_start+PATCH_SIZE**2) )
          end do
       end do
    end do

    ! Load finer scales if present, looping over each domain
    l = 1
    do while (level_end > l) ! New level was added -> proceed to it
       l = level_end 
       if (rank == 0) write (6,'(A,i2)') 'Loading level ', l
       do d = 1, size(grid)
          old_n_patch = grid(d)%patch%length
          do j = 1, grid(d)%lev(l)%length
             p_par = grid(d)%lev(l)%elts(j)
             do k = 1, zlevels
                do v = S_MASS, S_VELO
                   read (fid_no(d)) (wav_coeff(v,k)%data(d)%elts(i), i=MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start+1, &
                        MULT(v)*(grid(d)%patch%elts(p_par+1)%elts_start+PATCH_SIZE**2))
                   read (fid_no(d)) (trend_wav_coeff(v,k)%data(d)%elts(i), i=MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start+1, &
                        MULT(v)*(grid(d)%patch%elts(p_par+1)%elts_start+PATCH_SIZE**2))
                end do
             end do

             read (fid_gr(d)) child_required
             do c = 1, N_CHDRN
                if (child_required(c)) call refine_patch1 (grid(d), p_par, c - 1)
             end do
          end do

          do p_par = 2, old_n_patch
             do c = 1, N_CHDRN
                p_chd = grid(d)%patch%elts(p_par)%children(c)
                if (p_chd+1 > old_n_patch) call refine_patch2 (grid(d), p_par - 1, c - 1)
             end do
          end do
       end do
       call post_refine
    end do

    do d = 1, size(grid)
       close(fid_no(d)); close(fid_gr(d))
    end do

    wav_coeff%bdry_uptodate       = .False.
    trend_wav_coeff%bdry_uptodate = .False.
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

    l2error = sqrt(sum_real(l2error))
    maxerror = sync_max_d(maxerror)

    if (rank == 0) then
       write (6,'(A)') '-------------------------------------------------------------------------------------------------'
       write (6,'(A,i2,A,/)') 'Heikes-Randall optimizations of level ', level_start-1, ' grid:'
       write (6,'(A,2(es10.4,A))') 'Grid quality before optimization = ', maxerror, ' [m] (linf) ', l2error, ' [m] (l2)'
    end if
    
    fid = get_fid()
    if (level_start /= level_end) then
       write (0,'(i2,1x,i2)') level_end, level_start
       write (0,'(A)') "Reading HR grid points for level_start not equal to level_end not implemented"
       return
    end if

    write (filename, '(A,I1)')  "../extern/grid_HR/J", level_start-1
    open (unit=fid, file=filename)

    p = 1
    do d_HR = 1, N_ICOSAH_LOZANGE
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
    maxerror = sync_max_d(maxerror)
    if (rank == 0) then
       write (6,'(A,2(es10.4,A))') 'Grid quality after optimization  = ', maxerror, ' [m] (linf) ', l2error, ' [m] (l2)'
       write (6,'(A)') '(distance between midpoints of primal and dual edges)'
       write (6,'(A,/)') '-------------------------------------------------------------------------------------------------'
    end if
  end subroutine read_HR_optim_grid

  integer function dom_id_from_HR_id (d_HR)
    ! d_HR: lozange id as used by Heikes & Randall (starts from 1)
    ! results: domain id (starts from 0)
    implicit none
    integer :: d_HR

    dom_id_from_HR_id = modulo(d_HR,2)*5 + modulo(d_HR/2-1,5)
  end function dom_id_from_HR_id

  integer function sub_dom_id_from_HR_sub_id (sub_id)
    ! sub_id: lozange sub id as used by Heikes & Randall (starts from 1)
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
    implicit none
    integer :: d, l, max_output_level, num

    ! FIXME cleaner would be to use init_io routine
    call init_Float_Field (active_level, AT_NODE)

    do d = 1, size(grid)
       num = grid(d)%node%length
       call init(active_level%data(d), num)
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
end module io_mod
