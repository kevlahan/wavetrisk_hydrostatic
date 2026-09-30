!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: physics_simple.f90
! Author: Gabrielle Ching-Johnson, Nicholas Kevlahan
!
! Date Revised: 2025-05
!
! Contains all subroutines needed to compute a Backwards Euler step using simple physics package
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module physics_simple_mod
  
  use, intrinsic :: iso_c_binding, only : C_BOOL
  
  use kind_mod,   only : dp, sp
  use shared_mod, only : DAY, EDGE, N_BDRY, N_VARIABLE, zlevels, S_MASS, S_TEMP, S_VELO, RT, DG, UP, &
       Nsoil, dt, grav_accel, p_top, time, zlevels, zmin, z_null 

  use diagnostics_mod,    only : cal_surf_press, kinetic_energy, integrate_pressure_up
  use domain_mod,         only : Domain, exner, exner_fun, mass, temp, mean_m, mean_t, grid, sol, sol_mean, topography, id_edge, idx
  use domain_ops_mod,     only : apply_no_bdry2
  use geom_mod,           only : cart2sph
  use init_physics_mod,   only : physics_firstcall_flag
  use single_column_mod,  only : change_latitude_longitude, physics_call_single_col
  use utils_mod,          only : interp

#ifdef PHYSICS_COLUMN_LOCAL
  use iso_fortran_env, only: int64
  use shared_mod, only: compressible,c_p,p_0,kappa
  use arch_mod, only: rank
  use domain_ops_mod, only: apply_onescale_to_patch
  use init_mod, only: surf_geopot
  use parallel_block_primitive_kernels_mod, only: pressure_layer
#endif
  implicit none

  private
  public :: physics_simple_step, Rayleigh_friction
  
#ifdef PHYSICS_COLUMN_LOCAL
  ! Integer visitation metadata only. Numerical hydrostatic state belongs to
  ! the active column; no replacement Domain-shaped pressure image is kept.
  type :: Physics_Visits
     integer,allocatable :: surface(:),seen(:)
  end type
  type(Physics_Visits),allocatable :: column_visits(:)
  integer(int64) :: column_count=0,layer_checks=0
#endif
  real(sp) :: phys_dt
contains
  subroutine physics_simple_step (h)
    ! Uses simple physics modules to take a Backwards Euler step for physics using time step dt set by dynamics
    implicit none
    real(dp) :: h ! time step
#ifdef PHYSICS_COLUMN_LOCAL
    integer :: d,p
#endif
    
    phys_dt = real (h, kind=sp)
    
#ifdef PHYSICS_COLUMN_LOCAL
    if(.not.compressible) error stop "local Simple-physics columns require compressible dynamics"
    if(allocated(column_visits)) error stop "unfinished Simple-physics column epoch"
    allocate(column_visits(size(grid)))
    column_count=0;layer_checks=0
    do d=1,size(grid)
       allocate(column_visits(d)%surface(grid(d)%node%length),column_visits(d)%seen(grid(d)%node%length))
       column_visits(d)%surface=0;column_visits(d)%seen=0
       do p=3,grid(d)%patch%length
          call apply_onescale_to_patch(register_surface_visit,grid(d),p-1,0,0,1)
       end do
    end do
#endif
#if !defined(PHYSICS_COLUMN_LOCAL) || defined(PHYSICS_COLUMN_ORACLE)
    call cal_surf_press (sol(1:N_VARIABLE,1:zlevels))
#endif

    ! Compute Simple Physics split step on all columns
    call apply_no_bdry2 (physics_call, z_null)

    sol%bdry_uptodate      = .false.
    physics_firstcall_flag = .false.
#ifdef PHYSICS_COLUMN_LOCAL
    write(*,'(a,i0,a,2(1x,i0))') 'Local physics columns rank=',rank,' columns/layer checks=',column_count,layer_checks
    deallocate(column_visits)
#endif
  end subroutine physics_simple_step

#ifdef PHYSICS_COLUMN_LOCAL
  subroutine register_surface_visit(dom,i,j,k,offs,dims)
    type(Domain),intent(inout)::dom
    integer,intent(in)::i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer::id
    id=idx(i,j,offs,dims)+1
    column_visits(dom%id+1)%surface(id)=column_visits(dom%id+1)%surface(id)+1
  end subroutine
#ifdef PHYSICS_COLUMN_ORACLE
  subroutine check_physics_value(actual,expected)
    real(dp),intent(in)::actual,expected
    if(transfer(actual,0_int64)/=transfer(expected,0_int64)) then
       write(*,*) 'Physics hydrostatic mismatch rank, local, reference:',rank,actual,expected
       error stop "independent physics hydrostatic mismatch"
    end if
  end subroutine
#endif
#endif

  subroutine Rayleigh_friction (dom, i, j, zlev, offs, dims)
    ! Rayleigh friction in atmospheric boundary layer
    ! (based on Held-Suarez)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    character(5) :: smth = "tropo"          ! where to apply Rayleigh friction (tropo, strat)
    real(dp)     :: k_f     = 1.0_dp  / DAY ! Rayleigh friction
    real(dp)     :: sigma_b1 = 0.7_dp       ! normalized tropopause pressure height
    real(dp)     :: sigma_b2 = 0.1_dp       ! normalized stratosphere friction pressure height

    integer                        :: d, id, id_i, k, l
    integer,  dimension(1:EDGE)    :: id_e
    real(dp)                       :: k_v, rho_dz, sigma
    real(dp), dimension(0:zlevels) :: Pl
    real(dp), dimension(1:zlevels) :: Pk

    d = dom%id + 1
    
    id   = idx (i, j, offs, dims)
    id_i = id + 1
    id_e = id_edge (id)
    
    Pl(zlevels) = p_top
    do l = zlevels-1, 0, -1
       k = l + 1
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(id_i) 
       Pl(l)  = Pl(l+1) + grav_accel * rho_dz
       Pk(k) = interp (Pl(l), Pl(l+1))
    end do

    do k = 1, zlevels
       sigma = (Pk(k) - p_top) / (Pl(0) - p_top)

       select case (smth)
       case ("tropo")
          k_v = k_f * max (0.0_dp,   (sigma - sigma_b1) / (1.0_dp - sigma_b1))
       case ("strat")
          k_v = k_f * max (0.0_dp, - (sigma - sigma_b2) / sigma_b2)
       case default
          k_v = 0.0_dp
       end select

       sol(S_VELO,k)%data(d)%elts(id_e) = (1.0_dp - dt * k_v) * sol(S_VELO,k)%data(d)%elts(id_e)
    end do
  end subroutine Rayleigh_friction

  subroutine physics_call (dom, p_null, i, j, zlev, offs, dims, is_pole)
    !-----------------------------------------------------------------------------------
    !
    !   Backwards Euler physics step on a single element/column
    !
    !-----------------------------------------------------------------------------------

    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, is_pole, p_null, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, id_i, k, mask

    ! Variables needed for simple physics
    real(sp)                       :: nth_day, day_fraction   ! day in simulation, fraction of the day
    real(sp), dimension(1:zlevels) :: phys_Play               ! pressure at layer centres
    real(sp), dimension(0:zlevels) :: phys_Pint               ! pressure at layer interfaces
    real(sp), dimension(1:1)       :: phys_Phisurf            ! surface geopotential

    real(sp), dimension(1:zlevels) :: phys_Theta              ! potential temperature
    real(sp), dimension(1:zlevels) :: phys_Phi                ! geopotential
    real(sp), dimension(1:zlevels) :: phys_U, phys_V, phys_W  ! velocities at edges
    real(sp), dimension(1:zlevels) :: phys_Umag               ! speed at nodes
    real(sp), dimension(1:Nsoil+1) :: Tsoil                   ! surface temp and soil temperatures

    real(dp)                       :: latitude, longitude     ! coordinates of the column
    real(dp), dimension(1:zlevels) :: rho_dz
    
    logical(kind=C_BOOL)           :: lastcall_flag = .false.

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    day_fraction = (real (time, kind=sp) - phys_dt) / real (DAY, kind=sp)
    nth_day      = floor (day_fraction)
    day_fraction = real (day_fraction - nth_day, kind=sp)

    ! Save prognostic variables in physics data structure
    call pack_physics_vars

    ! Get latitude and longitude of the column
    call cart2sph (dom%node%elts(id_i), longitude, latitude)

    ! Update physics latitude and longitude
    call change_latitude_longitude (real (latitude, kind=sp), real (longitude, kind=sp))

    ! Surface geopotential for this column
    phys_Phisurf = real (grav_accel * topography%data(d)%elts(id_i), kind=sp)

    ! Mask for adaptive grid
    mask = grid(d)%mask_n%elts(id_i)

    ! Backwards Euler step on current column
    call physics_call_single_col (1, zlevels, mask, &
         physics_firstcall_flag, lastcall_flag, nth_day, day_fraction, phys_dt, &
         phys_Play, phys_Pint, phys_Phi, phys_Phisurf, phys_Umag, &
         phys_U, phys_V, phys_W, phys_Theta, Tsoil) ! updated values

    ! Assign solution at t+h
    do k = 1, zlevels
       if (is_pole /= 1) sol(S_VELO,k)%data(d)%elts(id_edge(id)) = real ((/ phys_U(k), phys_V(k), phys_W(k) /), kind=dp)
       sol(S_TEMP,k)%data(d)%elts(id_i) = rho_dz(k) * real (phys_Theta(k), kind=dp) - sol_mean(S_TEMP,k)%data(d)%elts(id_i)
    end do

    ! Assign soil column solution at t+phys_dt to WAVETRISK data structure
    do k = zmin, 0
       sol(S_TEMP,k)%data(d)%elts(id_i) = real (Tsoil(abs(k)+1), kind=dp)
    end do
  contains
    subroutine pack_physics_vars
      ! Gathers all required variables from wavetrisk data structure for all layers of column into physics data structure 
      integer  :: k
      real(dp) :: rho_dz_theta
#ifdef PHYSICS_COLUMN_LOCAL
      integer :: visit
      real(dp) :: column_surface,local_lower,local_geopot,local_pressure,local_exner,local_geopot_lower
      if(.not.allocated(column_visits)) error stop "Simple-physics column epoch absent"
      if(column_visits(d)%seen(id_i)/=0) error stop "duplicate Simple-physics column consumer"
      column_visits(d)%seen(id_i)=1
      column_count=column_count+1
      column_surface=0.0_dp
      do k=1,zlevels
         do visit=1,column_visits(d)%surface(id_i)
            column_surface=column_surface+(sol_mean(S_MASS,k)%data(d)%elts(id_i)+sol(S_MASS,k)%data(d)%elts(id_i))
         end do
      end do
      column_surface=grav_accel*column_surface+p_top
      local_lower=column_surface
      local_geopot=surf_geopot(d,id_i)
#ifdef PHYSICS_COLUMN_ORACLE
      call check_physics_value(column_surface,dom%surf_press%elts(id_i))
#endif
      phys_Pint(0)=real(column_surface,kind=sp)
#else

      phys_Pint(0)  = real (dom%surf_press%elts(id_i), kind=4)
#endif

      do k = 1, zlevels
         mass   =>      sol(S_MASS,k)%data(d)%elts
         temp   =>      sol(S_TEMP,k)%data(d)%elts
         mean_m => sol_mean(S_MASS,k)%data(d)%elts
         mean_t => sol_mean(S_TEMP,k)%data(d)%elts
#if !defined(PHYSICS_COLUMN_LOCAL) || defined(PHYSICS_COLUMN_ORACLE)
         exner  =>       exner_fun(k)%data(d)%elts
#endif

         rho_dz(k)    = mass(id_i) + mean_m(id_i)
         rho_dz_theta = temp(id_i) + mean_t(id_i)

         ! Pressure at layers and interface and geopotential at next interface (set in dom%geopot)
#ifdef PHYSICS_COLUMN_LOCAL
         call pressure_layer(mass(id_i),mean_m(id_i),temp(id_i),mean_t(id_i),grav_accel,c_p,p_0,kappa, &
              local_lower,local_geopot,local_pressure,local_exner,local_geopot_lower)
#ifdef PHYSICS_COLUMN_ORACLE
         call integrate_pressure_up (dom, i, j, zlev, offs, dims)
         call check_physics_value(local_lower,dom%press_lower%elts(id_i))
         call check_physics_value(local_pressure,dom%press%elts(id_i))
         call check_physics_value(local_geopot,dom%geopot%elts(id_i))
         call check_physics_value(local_geopot_lower,dom%geopot_lower%elts(id_i))
         call check_physics_value(local_exner,exner(id_i))
         layer_checks=layer_checks+1
#endif
#else
         call integrate_pressure_up (dom, i, j, zlev, offs, dims)
#endif

         ! Input variables for simple physics module
#ifdef PHYSICS_COLUMN_LOCAL
         phys_Pint(k)=real(local_lower,kind=sp)
         phys_Play(k)=real(local_pressure,kind=sp)
         phys_Phi(k)=real(interp(local_geopot_lower,local_geopot),kind=sp)
#else
         phys_Pint(k)  = real (dom%press_lower%elts(id_i),                                  kind=sp)
         phys_Play(k)  = real (dom%press%elts(id_i),                                        kind=sp)
         phys_Phi(k)   = real (interp (dom%geopot_lower%elts(id_i), dom%geopot%elts(id_i)), kind=sp)
#endif
         phys_Umag(k)  = real (sqrt (2 * kinetic_energy (dom, i, j, k, offs, dims)),        kind=sp)
         phys_U(k)     = real (sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1),                    kind=sp)
         phys_V(k)     = real (sol(S_VELO,k)%data(d)%elts(EDGE*id+DG+1),                    kind=sp)
         phys_W(k)     = real (sol(S_VELO,k)%data(d)%elts(EDGE*id+UP+1),                    kind=sp)
         phys_Theta(k) = real (rho_dz_theta / rho_dz(k),                                    kind=sp) 

         nullify (mass, temp, mean_m, mean_t, exner)
      end do

      ! Retrieve surface temperature and soil column temperatures from dynamics
      ! (Tsoil(1) is surface temperature)
      do k = zmin, 0
         Tsoil(abs(k)+1) = real (sol(S_TEMP,k)%data(d)%elts(id_i), kind=sp)
      end do
    end subroutine pack_physics_vars
  end subroutine physics_call
end module physics_simple_mod
