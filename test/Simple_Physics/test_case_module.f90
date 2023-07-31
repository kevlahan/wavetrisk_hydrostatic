module test_case_mod
   ! Description: Test case module for the Simple physics
   ! Author: Gabrielle Ching-Johnson
   ! Date Revised: February 4th
   use comm_mpi_mod
   use utils_mod
   use init_mod
   ! Use case for c booleans for physics
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_BOOL
   implicit none

   ! Standard variables
   integer :: CP_EVERY, resume_init
   real(8) :: dt_cfl, total_cpu_time, dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim

   ! Test case variables
   real(8) :: T_0 ! Reference Temperature
   real(8) :: u_0, e_thick ! geostrophic wind speed, ekman layer thickness
   integer :: Nsoil
   ! Physics variables
   LOGICAL(KIND=C_BOOL) :: physics_firstcall_flag = .true. ! flag for the physics package, true if call physics for 1st time f
   type(Float_Field), dimension(:), pointer :: dzonal, dmerid ! pointers for the tendencies

contains
   subroutine assign_functions
      ! Assigns generic pointer functions to functions defined in test cases
      implicit none

      ! Standard functions
      apply_initial_conditions => apply_initial_conditions_case
      dump                     => dump_case
      load                     => load_case
      initialize_a_b_vert      => initialize_a_b_vert_case
      initialize_dt_viscosity  => initialize_dt_viscosity_case
      initialize_thresholds    => initialize_thresholds_case
      physics_scalar_flux      => physics_scalar_flux_case
      physics_velo_source      => physics_velo_source_case
      set_save_level           => set_save_level_case
      set_thresholds           => set_thresholds_case
      surf_geopot              => surf_geopot_case
      update                   => update_case
      z_coords                 => z_coords_case
   end subroutine assign_functions

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Not sure what these functions do/ how to update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function physics_scalar_flux_case (q, dom, id, idE, idNE, idN, v, zlev, type)
      ! Additional physics for the flux term of the scalar trend
      ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
      !
      ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
      use domain_mod
      implicit none

      real(8), dimension(1:EDGE)                           :: physics_scalar_flux_case
      type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
      type(domain)                                         :: dom
      integer                                              :: d, id, idE, idNE, idN, v, zlev
      logical, optional                                    :: type

      integer                    :: id_i
      real(8), dimension(1:EDGE) :: d_e, grad, l_e
      logical                    :: local_type

      if (present(type)) then
         local_type = type
      else
         local_type = .false.
      end if

      id_i = id + 1
      d = dom%id + 1

      if (Laplace_order == 0) then
         physics_scalar_flux_case = 0.0_8
      else
         if (.not.local_type) then ! usual flux at edges E, NE, N
            l_e =  dom%pedlen%elts(EDGE*id+1:EDGE*id_i)
            d_e =  dom%len%elts(EDGE*id+1:EDGE*id_i)
         else ! flux at SW corner
            l_e(RT+1) = dom%pedlen%elts(EDGE*idE+RT+1)
            l_e(DG+1) = dom%pedlen%elts(EDGE*idNE+DG+1)
            l_e(UP+1) = dom%pedlen%elts(EDGE*idN+UP+1)
            d_e(RT+1) = -dom%len%elts(EDGE*idE+RT+1)
            d_e(DG+1) = -dom%len%elts(EDGE*idNE+DG+1)
            d_e(UP+1) = -dom%len%elts(EDGE*idN+UP+1)
         end if

         ! Calculate gradients
         if (Laplace_order == 1) then
            grad = grad_physics (q(v,zlev)%data(d)%elts)
         elseif (Laplace_order == 2) then
            grad = grad_physics (Laplacian_scalar(v)%data(d)%elts)
         end if

         ! Complete scalar diffusion
         physics_scalar_flux_case = (-1)**Laplace_order * visc_sclr(v) * grad * l_e
      end if
   contains
      function grad_physics (scalar)
         implicit none
         real(8), dimension(1:EDGE) :: grad_physics
         real(8), dimension(:)      :: scalar

         grad_physics(RT+1) = (scalar(idE+1) - scalar(id+1))   / d_e(RT+1)
         grad_physics(DG+1) = (scalar(id+1)  - scalar(idNE+1)) / d_e(DG+1)
         grad_physics(UP+1) = (scalar(idN+1) - scalar(id+1))   / d_e(UP+1)
      end function grad_physics
   end function physics_scalar_flux_case

   function physics_velo_source_case (dom, i, j, zlev, offs, dims)
      ! Additional physics for the source term of the velocity trend
      use domain_mod
      implicit none

      real(8), dimension(1:EDGE)     :: physics_velo_source_case
      type(domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer                    :: id, visc_scale
      real(8), dimension(1:EDGE) :: diffusion

      id = idx (i, j, offs, dims)

      visc_scale = 1.0_8!(max_level/dom%level%elts(id+1))**(2*Laplace_order_init-1)

      if (Laplace_order == 0) then
         diffusion = 0.0_8
      else
         ! Calculate Laplacian of velocity
         diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu()) * visc_scale
      end if

      ! Total physics for source term of velocity trend
      physics_velo_source_case =  diffusion
   contains
      function grad_divu()
         implicit none
         real(8), dimension(3) :: grad_divu

         integer :: idE, idN, idNE

         idE  = idx (i+1, j,   offs, dims)
         idN  = idx (i,   j+1, offs, dims)
         idNE = idx (i+1, j+1, offs, dims)

         grad_divu(RT+1) = (divu(idE+1) - divu(id+1))  /dom%len%elts(EDGE*id+RT+1)
         grad_divu(DG+1) = (divu(id+1)  - divu(idNE+1))/dom%len%elts(EDGE*id+DG+1)
         grad_divu(UP+1) = (divu(idN+1) - divu(id+1))  /dom%len%elts(EDGE*id+UP+1)
      end function grad_divu

      function curl_rotu()
         implicit none
         real(8), dimension(3) :: curl_rotu

         integer :: idS, idW

         idS  = idx (i,   j-1, offs, dims)
         idW  = idx (i-1, j,   offs, dims)

         curl_rotu(RT+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*idS+UPLT+1))/dom%pedlen%elts(EDGE*id+RT+1)
         curl_rotu(DG+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+DG+1)
         curl_rotu(UP+1) = (vort(TRIAG*idW+LORT+1) - vort(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+UP+1)
      end function curl_rotu
   end function physics_velo_source_case
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Dynamics Initialization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_sol (dom, i, j, zlev, offs, dims)
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims

      type(Coord) :: x_i, x_E, x_N, x_NE
      integer     :: id, d, idN, idE, idNE
      real(8)     :: phi, p, p_s, pot_temp, sigma
      real(8), dimension (1:EDGE)     :: uvw

      d = dom%id+1

      id   = idx(i,   j,   offs, dims)
      idN  = idx(i,   j+1, offs, dims)
      idE  = idx(i+1, j,   offs, dims)
      idNE = idx(i+1, j+1, offs, dims)

      x_i  = dom%node%elts(id+1)
      x_E  = dom%node%elts(idE+1)
      x_N  = dom%node%elts(idN+1)
      x_NE = dom%node%elts(idNE+1)

      ! Surface pressure
      !dom%surf_press%elts(id+1) = surf_pressure (x_i)
      p_s = dom%surf_press%elts(id+1)

      ! Pressure at level zlev (layer center)
      p = 0.5 * (a_vert(zlev)+a_vert(zlev+1) + (b_vert(zlev)+b_vert(zlev+1))*p_s)

      ! Mass/Area = rho*dz at level zlev
      sol(S_MASS,zlev)%data(d)%elts(id+1) = a_vert_mass(zlev) + b_vert_mass(zlev)*p_s/grav_accel

      ! Potential temperature
      pot_temp = T_0 * (p/p_0)**(-kappa)
      !     call cal_theta_eq (p, p_s, lat, theta_equil, k_T)

      ! Mass-weighted potential temperature
      sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

      ! Set initial velocity field (of edges, zonal and meridonal should already be calculated for the entire domain)
      ! Calculate u,v,w edge components fromm zonal and meridional
      phi = dom%u_zonal%elts(id + 1)
      phi = dom%v_merid%elts(id+1)
      call interp_latlon_UVW (dom, i, j, zlev, offs, dims, uvw)
      ! Set initial edge velocity
      sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = uvw(RT+1)
      sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) = uvw(DG+1)
      sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = uvw(UP+1)

      ! Means are zero
      sol_mean(S_MASS,zlev)%data(d)%elts(id+1)                      = 0d0
      sol_mean(S_TEMP,zlev)%data(d)%elts(id+1)                      = 0d0
      sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
   end subroutine init_sol

   subroutine init_vels_and_surf(dom, i, j, zlev, offs, dims)
      ! Set the initial_geopotential and zonal & meridional velcity
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims
      integer :: id
      type(Coord) :: x_i

      id   = idx(i,   j,   offs, dims)
      x_i  = dom%node%elts(id+1)

      ! Set Surface pressure
      dom%surf_press%elts(id+1) = surf_pressure (x_i)

      call geopot_initialize(dom, i, j, zlev, offs, dims)

      call vel_initialize(dom, i, j, zlev, offs, dims)

   end subroutine init_vels_and_surf

   real(8) function surf_geopot_case (x_i)
      ! Surface geopotential set to zero
      implicit none
      Type(Coord) :: x_i

      ! Geopotential at sea level across sphere
      surf_geopot_case = 0

   end function surf_geopot_case

   real(8) function surf_pressure (x_i)
      ! Surface pressure
      implicit none
      type(Coord) :: x_i

      surf_pressure = p_0
   end function surf_pressure

   subroutine geopot_initialize(dom, i, j, zlev, offs, dims)
      ! Set the initial_geopotential
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims

      real(8) :: p, p_s, phi
      integer :: id_i

      id_i = idx(i, j, offs, dims) + 1

      p_s = dom%surf_press%elts(id_i)
      p = 0.5 * (a_vert(zlev)+a_vert(zlev+1) + (b_vert(zlev)+b_vert(zlev+1))*p_s)

      ! Geopotential
      phi = -R_d*T_0*LOG((p)/p_s)
      dom%geopot%elts(id_i) = phi

   end subroutine geopot_initialize

   subroutine vel_initialize (dom, i, j, zlev, offs, dims)
      ! Velocities intialized using theory of ekman layer
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims
      integer :: id_i
      real(8) :: phi, u, v, lat, lon
      type(Coord) :: x_i

      id_i = idx(i, j, offs, dims) + 1
      x_i  = dom%node%elts(id_i)
      call cart2sph (x_i, lon, lat)
      phi = dom%geopot%elts(id_i)

      u = u_0 * (1-(exp(-phi/(e_thick*grav_accel))*cos(phi/ (e_thick*grav_accel)))) ! Zonal velocity component
      dom%u_zonal%elts(id_i) = u*cos(lat)
      v = u_0 * (1-(exp(-phi/(e_thick*grav_accel))*sin(phi/ (e_thick*grav_accel)))) ! Meridional velocity component
      dom%v_merid%elts(id_i) = v*cos(lat)
   end subroutine vel_initialize

   subroutine apply_initial_conditions_case
      implicit none
      integer :: k, l

      do l = level_start, level_end
         do k = 1, zlevels
            call apply_onescale(init_vels_and_surf, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
            call apply_onescale (init_sol, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
         end do
      end do
   end subroutine apply_initial_conditions_case

   subroutine initialize_a_b_vert_case
      implicit none
      integer :: k

      ! Allocate vertical grid parameters
      allocate (a_vert(1:zlevels+1),    b_vert(1:zlevels+1))
      allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

      do k = 1, zlevels+1
         a_vert(k) = dble(k-1)/dble(zlevels) * p_top
         b_vert(k) = 1.0_8 - dble(k-1)/dble(zlevels)
      end do

      ! LMDZ grid
      !call cal_AB

      ! Set mass coefficients
      a_vert_mass = (a_vert(1:zlevels) - a_vert(2:zlevels+1))/grav_accel
      b_vert_mass =  b_vert(1:zlevels) - b_vert(2:zlevels+1)
   end subroutine initialize_a_b_vert_case

   subroutine initialize_thresholds_case
      ! Set default thresholds based on dimensional scalings of norms
      implicit none

      integer :: k

      allocate (threshold(1:N_VARIABLE,zmin:zlevels));     threshold     = 0.0_8
      allocate (threshold_def(1:N_VARIABLE,zmin:zlevels)); threshold_def = 0.0_8

      lnorm(S_MASS,:) = dPdim/grav_accel
      do k = 1, zlevels
         lnorm(S_TEMP,k) = (a_vert_mass(k) + b_vert_mass(k)*Pdim/grav_accel)*dTempdim
      end do
      lnorm(S_TEMP,:) = lnorm(S_TEMP,:) + Tempdim*lnorm(S_MASS,:) ! Add component due to tendency in mass
      lnorm(S_VELO,:) = Udim
      threshold_def = tol * lnorm
   end subroutine initialize_thresholds_case

   subroutine initialize_dt_viscosity_case
      ! Initializes viscosity
      implicit none
      real(8) :: area, C_divu, C_sclr, C_rotu, tau_divu, tau_rotu, tau_sclr

      area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
      dx_min = sqrt (4/sqrt(3.0_8) * area)         ! edge length of average triangle

      ! Diffusion constants
      C_sclr = 2d-3       ! <= 1.75e-2 for hyperdiffusion (lower than exact limit 1/6^2 = 2.8e-2 due to non-uniform grid)
      C_divu = 2d-3    ! <= 1.75e-2 for hyperdiffusion (lower than exact limit 1/6^2 = 2.8e-2 due to non-uniform grid)
      C_rotu = C_sclr / 4**Laplace_order_init ! <= 1.09e-3 for hyperdiffusion (lower than exact limit 1/24^2 = 1.7e-3 due to non-uniform grid)

      ! CFL limit for time step
      dt_cfl = cfl_num*dx_min/(wave_speed+Udim) * 0.85 ! corrected for dynamic value
      dt_init = dt_cfl

      tau_sclr = dt_cfl / C_sclr
      tau_divu = dt_cfl / C_divu
      tau_rotu = dt_cfl / C_rotu

      if (Laplace_order_init == 0) then
         visc_sclr = 0.0_8
         visc_divu = 0.0_8
         visc_rotu = 0.0_8
      elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
         visc_sclr = dx_min**(2*Laplace_order_init) / tau_sclr
         visc_rotu = dx_min**(2*Laplace_order_init) / tau_rotu
         visc_divu = dx_min**(2*Laplace_order_init) / tau_divu
      elseif (Laplace_order_init > 2) then
         if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
         stop
      end if

      if (rank == 0) then
         write (6,'(/,3(a,es8.2),a,/)') "dx_min  = ",dx_min/KM, " [km] dt_cfl = ", dt_cfl, " [s] tau_sclr = ", tau_sclr/HOUR, " [h]"
         write (6,'(3(a,es8.2),/)') "C_sclr = ", C_sclr, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
         write (6,'(4(a,es8.2))') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
            " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
            " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
      end if
   end subroutine initialize_dt_viscosity_case

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine set_thresholds_case
      ! Set thresholds dynamically (trend or sol must be known)
      use lnorms_mod
      use wavelet_mod
      implicit none
      integer                                     :: k
      real(8), dimension(1:N_VARIABLE,zmin:zlevels) :: threshold_new
      character(3), parameter                     :: order = "inf"

      if (default_thresholds) then ! Initialize once
         threshold_new = threshold_def
      else
         call cal_lnorm_sol (sol, order)
         !threshold_new = max (tol*lnorm, threshold_def) ! Avoid very small thresholds before instability develops
         threshold_new = tol*lnorm
      end if

      if (istep >= 10) then
         threshold = 0.01*threshold_new + 0.99*threshold
      else
         threshold = threshold_new
      end if
   end subroutine set_thresholds_case

   subroutine set_save_level_case
      ! Determines closest vertical level to desired pressure
      implicit none
      integer :: k
      real(8) :: dpress, p, save_press

      dpress = 1d16; save_zlev = 0
      do k = 1, zlevels
         p = 0.5 * (a_vert(k)+a_vert(k+1) + (b_vert(k)+b_vert(k+1))*p_0)
         if (abs(p-pressure_save(1)) < dpress) then
            dpress = abs (p-pressure_save(1))
            save_zlev = k
            save_press = p
         end if
      end do
      if (rank==0) write (6,'(/,A,i2,A,f5.1,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate pressure = ", save_press/100, " hPa)"
   end subroutine set_save_level_case

   subroutine initialize_seed
      implicit none

      integer                            :: k
      integer, dimension(1:8)            :: values
      integer, dimension(:), allocatable :: seed

      call date_and_time (values=values)
      call random_seed(size=k)
      allocate (seed(1:k))
      seed = values(8)
      call random_seed (put=seed)
   end subroutine initialize_seed

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Reading/Printing!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_test_case_parameters
      implicit none
      integer            :: k, v
      integer, parameter :: fid = 500, funit = 400
      real(8)            :: press_save
      character(255)     :: command, filename, varname
      character(2)       :: var_file
      logical            :: file_exists

      ! Find input parameters file name
      if (command_argument_count() >= 1) then
         CALL getarg (1, filename)
      else
         filename = 'test_case.in'
      end if
      if (rank == 0) write (6,'(a,A)') "Input file = ", trim (filename)

      open(unit=fid, file=filename, action='READ')
      read (fid,*) varname, test_case
      read (fid,*) varname, run_id
      read (fid,*) varname, compressible
      read (fid,*) varname, max_level
      read (fid,*) varname, zlevels
      read (fid,*) varname, uniform
      read (fid,*) varname, remap
      read (fid,*) varname, remapscalar_type
      read (fid,*) varname, remapvelo_type
      read (fid,*) varname, iremap
      read (fid,*) varname, default_thresholds
      read (fid,*) varname, tol
      read (fid,*) varname, optimize_grid
      read (fid,*) varname, adapt_dt
      read (fid,*) varname, cfl_num
      read (fid,*) varname, timeint_type
      read (fid,*) varname, press_save
      read (fid,*) varname, Laplace_order_init
      read (fid,*) varname, dt_write
      read (fid,*) varname, CP_EVERY
      read (fid,*) varname, rebalance
      read (fid,*) varname, time_end
      read (fid,*) varname, resume_init
      close(fid)

      allocate (pressure_save(1))
      pressure_save(1) = 1.0d2*press_save
      dt_write = dt_write * DAY
      time_end = time_end * DAY
      resume   = resume_init
      Laplace_order = Laplace_order_init

      ! Bins for zonal statistics
      !nbins = sqrt (10*4**max_level/2) ! consistent with maximum resolution
      nbins = 300
      allocate (Nstats(zlevels,nbins), Nstats_glo(zlevels,nbins)) ; Nstats = 0 ; Nstats_glo = 0
      allocate (zonal_avg(zlevels,nbins,nvar_zonal), zonal_avg_glo(zlevels,nbins,nvar_zonal))
      zonal_avg = 0.0_8; zonal_avg_glo = 0.0_8
      allocate (bounds(1:nbins-1))
      dbin = 1.8d2/nbins
      bounds = -90+dbin + dbin*(/ (ibin, ibin = 0, nbins-1) /)

      ! Initialize rank 0 with saved statistics data if present
      if (rank == 0) then
         inquire (file = trim(run_id)//'.3.tgz', exist = file_exists)
         if (file_exists) then
            command = 'tar xzf '//trim(run_id)//'.3.tgz'
            call system (command)

            write (var_file, '(i2.2)') 00
            open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="UNFORMATTED", action='READ')
            read (funit) Nstats
            close (funit)

            do v = 1, nvar_zonal
               write (var_file, '(i2)') v+10
               open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='READ')
               do k = zlevels, 1, -1
                  read (funit,*) zonal_avg(k,:,v)
               end do
               close (funit)
            end do

            ! Convert variances to sums of squares for statistics computation
            zonal_avg(:,:,2) = zonal_avg(:,:,2) * (Nstats - 1)
            do v = 6, nvar_zonal
               zonal_avg(:,:,v) = zonal_avg(:,:,v) * (Nstats - 1)
            end do
         end if
      end if
   end subroutine read_test_case_parameters

   subroutine print_test_case_parameters
      implicit none

      if (rank==0) then
         write (6,'(a)') &
            '********************************************************** Parameters &
            ************************************************************'
         write (6,'(a)')        "RUN PARAMETERS"
         write (6,'(a,a)')      "test_case           = ", trim (test_case)
         write (6,'(a,a)')      "run_id              = ", trim (run_id)
         write (6,'(a,l1)')     "compressible        = ", compressible
         write (6,'(a,i3)')     "min_level           = ", min_level
         write (6,'(a,i3)')     "max_level           = ", max_level
         write (6,'(a,i5)')     "number of domains   = ", N_GLO_DOMAIN
         write (6,'(a,i5)')     "number of processors = ", n_process
         write (6,'(a,i5)')     "DOMAIN_LEVEL        = ", DOMAIN_LEVEL
         write (6,'(a,i5)')     "PATCH_LEVEL         = ", PATCH_LEVEL
         write (6,'(a,i3)')     "zlevels             = ", zlevels
         write (6,'(a,l1)')     "uniform             = ", uniform
         write (6,'(a,l1)')     "remap               = ", remap
         write (6,'(a,a)')      "remapscalar_type    = ", trim (remapscalar_type)
         write (6,'(a,a)')      "remapvelo_type      = ", trim (remapvelo_type)
         write (6,'(a,i3)')     "iremap              = ", iremap
         write (6,'(a,l1)')     "default_thresholds  = ", default_thresholds
         write (6,'(a,es10.4)') "tolerance           = ", tol
         write (6,'(a,i1)')     "optimize_grid       = ", optimize_grid
         write (6,'(a,l1)')     "adapt_dt            = ", adapt_dt
         write (6,'(a,es10.4)') "cfl_num             = ", cfl_num
         write (6,'(a,a)')      "timeint_type        = ", trim (timeint_type)
         write (6,'(a,es10.4)') "pressure_save (hPa) = ", pressure_save(1)/100
         write (6,'(a,i1)')     "Laplace_order       = ", Laplace_order_init
         write (6,'(a,i2)')     "n_diffuse           = ", n_diffuse
         write (6,'(a,es10.4)') "dt_write (day)      = ", dt_write/DAY
         write (6,'(a,i6)')     "CP_EVERY            = ", CP_EVERY
         write (6,'(a,l1)')     "rebalance           = ", rebalance
         write (6,'(a,es10.4)') "time_end (day)      = ", time_end/DAY
         write (6,'(a,i6)')     "resume              = ", resume_init

         write (6,'(/,a)')      "STANDARD PARAMETERS"
         write (6,'(a,es10.4)') "radius              = ", radius
         write (6,'(a,es10.4)') "omega               = ", omega
         write (6,'(a,es10.4)') "p_0   (hPa)         = ", p_0/100
         write (6,'(a,es10.4)') "p_top (hPa)         = ", p_top/100
         write (6,'(a,es10.4)') "R_d                 = ", R_d
         write (6,'(a,es10.4)') "c_p                 = ", c_p
         write (6,'(a,es10.4)') "c_v                 = ", c_v
         write (6,'(a,es10.4)') "gamma               = ", gamma
         write (6,'(a,es10.4)') "kappa               = ", kappa

         write (6,'(/,a)')      "TEST CASE PARAMETERS"
         write (6,'(a,es10.4)') "T_0                 = ", T_0
         write (6,'(a)') &
            '*********************************************************************&
            ************************************************************'
      end if
   end subroutine print_test_case_parameters

   subroutine print_log
      ! Prints out and saves logged data to a file
      implicit none

      integer :: min_load, max_load
      real(8) :: avg_load, rel_imbalance, timing

      timing = get_timing(); total_cpu_time = total_cpu_time + timing

      call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)

      if (rank == 0) then
         write (6,'(a,es12.6,4(a,es8.2),a,i2,a,i12,4(a,es8.2,1x))') &
            'time [d] = ', time/DAY, &
            ' dt [s] = ', dt, &
            '  mass tol = ', sum (threshold(S_MASS,1:zlevels))/zlevels, &
            ' temp tol = ', sum (threshold(S_TEMP,1:zlevels))/zlevels, &
            ' velo tol = ', sum (threshold(S_VELO,1:zlevels))/zlevels, &
            ' Jmax = ', level_end, &
            ' dof = ', sum (n_active), &
            ' min rel mass = ', min_mass, &
            ' mass error = ', mass_error, &
            ' balance = ', rel_imbalance, &
            ' cpu = ', timing

         write (12,'(5(es15.9,1x),i2,1x,i12,1x,4(es15.9,1x))')  &
            time/DAY, dt, sum (threshold(S_MASS,1:zlevels))/zlevels, sum (threshold(S_TEMP,1:zlevels))/zlevels, &
            sum (threshold(S_VELO,1:zlevels))/zlevels, level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
      end if
   end subroutine print_log
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine update_case
      ! Update means, bathymetry and penalization mask
      ! not needed in this test case
      use wavelet_mod
      implicit none
      integer :: d, k, l, p

   end subroutine update_case

   subroutine dump_case (fid)
      implicit none
      integer :: fid

      write (fid) itime
      write (fid) iwrite
      write (fid) threshold
   end subroutine dump_case

   subroutine load_case (fid)
      implicit none
      integer :: fid

      read (fid) itime
      read (fid) iwrite
      read (fid) threshold
   end subroutine load_case
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Physics Initialization and Plugins !!!!!!!!!!!!!!!!!!!!!
   !! Should be in its own module
   subroutine init_physics
      !-----------------------------------------------------------------------------------
      !
      !   Description: Sets necessary function pointers for the physics and initializes the
      !                main physics parameters.
      !
      !   Assumptions: Assumes that physics grid parameters has been initialized (ie called
      !                init_soild_grid(_default)).
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------

      !Any use cases from either physics or dynamics (ie zlevels)
      !Physics
      use iniphyparam_mod, ONLY : iniphyparam
      use single_column_mod, ONLY : initialize_extra_levels
      use read_param_mod
      !Dynamics
      use main_mod, ONLY : dt_new
      implicit none
      character(255) :: command, param_file

      !Set physics function pointers (if using) ! it is here where I use soil_mod flag to set soil
      read_paramr_plugin => read_paramr
      read_parami_plugin => read_parami
      read_paramb_plugin => read_paramb

      !write physics read in parameters to file
      write(param_file, '(A,I4.4)') trim(run_id)//'.physics_params.', rank
      call write_physics_params(9*rank,param_file)

      !Call initialization of physics parameters
      open(unit=9*rank, file=trim(param_file), form="FORMATTED", action='READ')
      call iniphyparam(dt_new,  DAY, radius, grav_accel, R_d, c_p)
      close(9*rank)

      !Delete extra files
      write(command, '(A,A)') '\rm ', trim(param_file)
      call system(command)

      ! physics single column module extra levels initialization (as needs soil flag set in iniphyparam)
      call initialize_extra_levels(Nsoil+1)
   end subroutine init_physics

   subroutine init_soil_grid_default
      !-----------------------------------------------------------------------------------
      !
      !   Description: Initialize physics package with dummy longitude & latitude values
      !                and grid parameters. Also set number of soil levels and zmin.
      !
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      ! Physics use cases
      use comgeomfi, ONLY : init_comgeomfi, nsoilmx

      ! local variables
      real(8) :: lat(1), long(1)

      if (.not. soil_mod) then
         PRINT*, 'STOP!! Cannot use default (init_soil_grid_default) to set zmin when soil_mod flag indicates false'
         STOP
      end if

      ! Dummy latatiude and longitude for initialization
      lat(1) = 0
      long(1) = 0

      ! call grid initialization for the physics
      call init_comgeomfi(1, zlevels, long, lat)

      ! set the number of soil levels and zmin (lowest vertical level index)
      Nsoil = nsoilmx
      zmin = -Nsoil
   end subroutine init_soil_grid_default

   subroutine init_soil_grid
      !-----------------------------------------------------------------------------------
      !
      !   Description: Initialize physics package with dummy longitude & latitude values
      !                and grid parameters. Also set number of soil levels and zmin.
      !
      !   Assumptions: Nsoil is set in Simple_Physics.f90 under the test case parameters.
      !                for case where soil model is turned on (ie soil_mod flag = true)
      !
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      ! Physics use cases
      use comgeomfi, ONLY : init_comgeomfi

      ! local variables
      real(8) :: lat(1), long(1)

      ! Dummy latatiude and longitude for initialization
      lat(1) = 0
      long(1) = 0

      ! set the zmin (lowest vertical level index) (! See Assumptions)
      if (soil_mod) then
         zmin = -Nsoil
         ! call grid initialization for the physics
         call init_comgeomfi(1, zlevels, long, lat, Nsoil)
      else
         Nsoil = 0
         zmin = 0
         ! call grid initialization for the physics
         call init_comgeomfi(1, zlevels, long, lat)
      end if

   end subroutine init_soil_grid

   subroutine physics_checkpoint_restart
      !-----------------------------------------------------------------------------------
      !
      !   Description: Initialize physics cases for when checkpointing is used. Only to be
      !                called when cp_idx > 0.
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      ! Physics use cases
      use phyparam_mod, ONLY : alloc, precompute, zday_last, icount

      ! local variables
      real(8) :: day_fraction, nth_day

      ! set flag for first call to physics false (thus the soil levels will get updated)
      physics_firstcall_flag = .false.

      !call allocation for physics call (usually done on the physics first call)
      call alloc(1, zlevels)

      !call precompute for physics call(usually done on the physics first call)
      call precompute()

      !set the previous day in physics
      day_fraction = (time-dt)/DAY
      nth_day = FLOOR(day_fraction)
      day_fraction = day_fraction-nth_day

      zday_last = nth_day + day_fraction - dt/DAY

   end subroutine physics_checkpoint_restart

   subroutine read_paramr(name, defval, val, comment)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Physics plugin to read reals from a file
      !
      !   Assumption: File is already open with unit number: 9*rank
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      character(*), INTENT(IN) :: name, comment
      real(8), INTENT(IN)         :: defval
      real(8), INTENT(OUT)        :: val

      character(15) :: param_name
      character :: equal_sign

      ! read line from input file
      read(9*rank,*) param_name, equal_sign, val

      if (trim(param_name) .ne. trim(name)) val = defval

   end subroutine

   subroutine read_parami(name, defval, val, comment)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Physics plugin to read integers from a file
      !
      !   Assumption: File is already open with unit number: 9*rank
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      character(*), INTENT(IN) :: name, comment
      integer, INTENT(IN)         :: defval
      integer, INTENT(OUT)        :: val

      character(15) :: param_name
      character :: equal_sign

      ! read line from input file
      read(9*rank,*) param_name, equal_sign, val

      if (trim(param_name) .ne. trim(name)) val = defval

   end subroutine

   subroutine read_paramb(name, defval, val, comment)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Physics plugin to read logicals from a file
      !
      !   Assumption: File is already open with unit number: 9*rank
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      character(*), INTENT(IN) :: name, comment
      logical, INTENT(IN)         :: defval
      logical, INTENT(OUT)        :: val

      character(15) :: param_name
      character :: equal_sign

      ! read line from input file
      read(9*rank,*) param_name, equal_sign, val

      if (trim(param_name) .ne. trim(name)) val = defval

   end subroutine

   subroutine write_physics_params(file_unit,file_params)
      integer :: file_unit
      character(*) :: file_params
      logical :: physics_write
      !-----------------------------------------------------------------------------------
      !
      !   Description: Write desired physics parameters to a file.
      !
      !   Notes: Takes into account if mpi is being used, so file_params contains
      !           the file name
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------

      if (rank == 0) then
         physics_write = .true.
      else
         physics_write = .false.
      end if

      open(unit=file_unit, file=trim(file_params), form="FORMATTED", action='WRITE', status='REPLACE')

      !write the desired physics parameters
      write(file_unit,*) "planet_rat = ", radius
      write(file_unit,*) "g = ", grav_accel
      write(file_unit,*) "cpp = ", c_p
      write(file_unit,*) "mugaz = ", 28.9702532_8  !8314.46261815324/R_d
      write(file_unit,*) "unjours = ", DAY
      write(file_unit,*) "year_day = ", 365
      write(file_unit,*) "periheli = ", 150
      write(file_unit,*) "aphelie = ", 150
      write(file_unit,*) "peri_day = ", 0.
      write(file_unit,*) "obliquit = ", 0
      write(file_unit,*) "Cd_mer = ", 0.01_8
      write(file_unit,*) "Cd_ter = ", 0.01_8
      write(file_unit,*) "I_mer = ", 3000.
      write(file_unit,*) "I_ter = ", 3000.
      write(file_unit,*) "alb_ter = ", 0.112
      write(file_unit,*) "alb_mer = ", 0.112
      write(file_unit,*) "emi_mer = ", 1.
      write(file_unit,*) "emi_ter = ", 1.
      write(file_unit,*) "emin_turb = ", 1.e-16
      write(file_unit,*) "lmixmin = ", 100
      write(file_unit,*) "coefvis = ", 0.99_8
      write(file_unit,*) "coefir = ", 0.08_8
      write(file_unit,*) "callrad = ", .true.
      write(file_unit,*) "calldifv = ", .true.
      write(file_unit,*) "calladj = ", .true.
      write(file_unit,*) "callsoil = ", soil_mod
      write(file_unit,*) "season = ", .false.
      write(file_unit,*) "diurnal = ", .true.
      write(file_unit,*) "lverbose = ", physics_write
      write(file_unit,*) "period_sort = ", 1.

      close(file_unit)
   end subroutine write_physics_params

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Routines Needed to Call the Physics !!!!!!!!!!!!!!!!!!!!
   subroutine trend_physics(q, dq)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Trend used to call the physics for each column on each domain. Also
      !                updates velocity edge tendencies once all elements on a domain have
      !                been set in the placeholder trends.
      !
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------

      implicit none
      type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target  :: q !solutions
      type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target  :: dq !tendencies
      integer :: p, d, dom_length, k
      type(Float_Field), dimension(1:zlevels), target :: trend_zonal, trend_merid

      ! Initialize zonal and meridional velocity placeholder trends
      trend_zonal = dq(S_MASS,1:zlevels)
      trend_merid = dq(S_MASS,1:zlevels)
      do k = 1, zlevels
         call zero_float_field (trend_zonal(k), S_MASS)
         call zero_float_field (trend_merid(k), S_MASS)
      end do
      dzonal => trend_zonal
      dmerid => trend_merid

      call update_array_bdry (sol, NONE, 27)

      !Get Surface pressure of all domains on rank
      call cal_surf_press_phys(q)

      ! Loop through each domain on rank to call physics
      do d = 1,size(grid)
         ! Loop through each patch on a domain
         do p = 3, grid(d)%patch%length
            call apply_onescale_to_patch(physics_call, grid(d), p-1, z_null, 0, 1)
         end do
      end do
      ! Update flags to false, once 1st call for all columns finished
      physics_firstcall_flag = .false.

      ! Update the edge tendencies of a domain (Can only occur once all velocities on a domain set)
      do k=1,zlevels
         do d = 1,size(grid)
            grid(d)%u_zonal%elts = dzonal(k)%data(d)%elts
            grid(d)%v_merid%elts = dmerid(k)%data(d)%elts
            !once all columns on domain has been updated, go patch by patch to change tendencies and update them
            do p = 3, grid(d)%patch%length
               call apply_onescale_to_patch(update_velocity_tendencies, grid(d), p-1, k, 0, 0)
            end do
         end do
      end do

      !nullify pointers
      nullify (dzonal, dmerid)

      ! Boundary update ! setting only the zleveld 1:zlevs float fields boundary to be updated
      dq%bdry_uptodate = .false.
   end subroutine trend_physics

   subroutine update_velocity_tendencies(dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Subroutine used to update the trend hybrid structure for a single element/column
      !                with the converted edge tendencies from the zonal and meridional tendencies.
      !
      !   Expectation: Expected that this subroutine is called/used by apply_one_scale_to_patch
      !                (or similar routines) routine of Wavetrisk to update the hybrid structure.
      !
      !   Assumptions: dom%u_zonal and dom%v_merid have been populated at all gird points for
      !                domain dom.
      !
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------

      ! Arguments
      type(Domain)                         :: dom ! domain the column is on
      integer                              :: i, j, zlev
      integer, dimension(N_BDRY+1)         :: offs
      integer, dimension(2,N_BDRY+1)       :: dims
      ! Local variables
      integer :: id, id_i
      real(8), dimension (1:EDGE)     :: uvw

      ! Get id of element
      id = idx(i, j, offs, dims)
      id_i = id + 1

      !Convert zonal and meridional tendencies to edge tendencies
      call interp_latlon_UVW(dom, i, j, zlev, offs, dims, uvw)

      !Update the trend
      trend(S_VELO, zlev)%data(dom%id+1)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = uvw

   end subroutine update_velocity_tendencies

   subroutine physics_call(dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Subroutine used to update calculate the trend for a single element/column.
      !                Converts the dynamics progrnostic vars to needed physics structure, calls
      !                the physics package to get tendencies for the column and sets temp
      !                tendencies back to hybrid structure. With the physics package, the surface
      !                temp and soil temp (if on) will be saved in the pot temp hybrid structure,
      !                for when load balancing occurs.
      !
      !   Expectation: Expected that this subroutine is called/used by apply_one_scale_to_patch
      !                (or similar routines) routine of Wavetrisk to calculate trend for the column
      !
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------

      ! Physics use cases
      USE callkeys,     ONLY : lverbose
      USE single_column_mod, ONLY : change_latitude_longitude, physics_call_single_col
      implicit none

      ! Arguments
      type(Domain)                     :: dom                  ! domain the column is on
      integer                          :: i, j, zlev
      integer, dimension(N_BDRY+1)     :: offs
      integer, dimension(2,N_BDRY+1)   :: dims

      ! Local variables for physics
      integer :: d,&   ! domain id
         id, &  ! id of the column in the domain -1
         id_i   ! id of column
      real(8), dimension(1:zlevels) :: phys_temp,&             ! Column temperature
         phys_zonal,&            ! Column zonal velocity
         phys_meridional, &      ! Column meridional velocity
         phys_geopot,&           ! Column geopotential
         phys_pres_lay, &        ! Column layer pressures
         phys_dtemp,&            ! Column temperature tendencies
         phys_dzonal,&           ! Column zonal velolcity tendencies
         phys_dmeridional, &     ! Column meridonal velcity tendencies
         phys_dsurfpres,&        ! Column surface pressure tendency
         phys_dtheta             ! Column potential temperature tendency
      real(8), dimension(1:Nsoil+1) :: surf_soil_temp           ! Column surface temp and soil temperatures
      real(8), dimension (0:zlevels) :: phys_pres_lev          ! Column interfaces pressure
      real(8) :: latitude, longitude                           ! Coordinates of the column
      real(8) :: nth_day, day_fraction                         ! Day in simulation, fraction of the day
      LOGICAL(KIND=C_BOOL) :: lastcall_flag=.false.

      ! Get domain id
      d = dom%id +1
      ! Get the id of the column in the domain
      id = idx (i, j, offs, dims)
      id_i = id + 1

      !Get the number of days and fraction of day
      day_fraction = (time-dt)/DAY
      nth_day = FLOOR(day_fraction)
      day_fraction = day_fraction-nth_day

      ! Set Surface pressure to what was initially calculated
      phys_pres_lev(0) = dom%press_lower%elts(id_i)

      ! Retrieve prognositc variables and save in physics data structure
      call retrieve_prog_vars

      ! Intialize tendencies to 0
      phys_dtemp = 0
      phys_dzonal = 0
      phys_dmeridional = 0
      phys_dsurfpres = 0

      ! Get latitude and longitude of the column
      call cart2sph(dom%node%elts(id_i), longitude, latitude)

      ! Update physics latitude and longitude
      call change_latitude_longitude(latitude, longitude)

      ! Call physics for current column
      call physics_call_single_col(1, zlevels, physics_firstcall_flag, lastcall_flag, nth_day, day_fraction, dt, &
         phys_pres_lev, phys_pres_lay, phys_geopot, phys_zonal,  phys_meridional, phys_temp, surf_soil_temp, &
         phys_dzonal, phys_dmeridional, phys_dtemp, phys_dsurfpres)

      lverbose = .false.
      ! Convert data received from physics back to hybrid structure
      call save_tendencies

   contains
      subroutine retrieve_prog_vars
         ! Description: Gathers all prognostic variables for all levels of the column into physics 2D data structure
         !              from the dynamics hybrid data structure
         !Arguments
         integer :: k
         real(8) :: theta, full_mass_theta, full_mass

         do k= 1,zlevels
            ! Get pressure at layer centers and interfaces of the column and geopotential at next interface (set in dom%geopot)
            call cal_press_geopot_layer (dom, i, j, k, offs, dims)
            phys_pres_lev(k) = dom%press_lower%elts(id_i)
            phys_pres_lay(k) = dom%press%elts(id_i)

            ! Get geopotential at layer center
            phys_geopot(k) = 0.5*(dom%geopot_lower%elts(id_i)+dom%geopot%elts(id_i))

            ! Convert the dynamics temp (pot temp ie theta) to temperature  (for entire column)
            ! use T = theta(p/p0)^kappa
            full_mass = sol(S_MASS, k)%data(d)%elts(id_i) + sol_mean(S_MASS, k)%data(d)%elts(id_i)
            full_mass_theta = sol(S_TEMP, k)%data(d)%elts(id_i) + sol_mean(S_TEMP, k)%data(d)%elts(id_i)
            theta = full_mass_theta/full_mass
            phys_temp(k) = (((phys_pres_lay(k)/dom%surf_press%elts(id_i))**kappa)*theta)

            ! Convert the edge velocities to zonal and meridional velocities
            velo => sol(S_VELO,k)%data(d)%elts
            velo1 => grid(d)%u_zonal%elts
            velo2 => grid(d)%v_merid%elts
            call interp_UVW_latlon(dom, i, j, k, offs, dims)
            phys_zonal(k) = velo1(id_i)
            phys_meridional(k) = velo2(id_i)
            nullify (velo, velo1, velo2)
         end do

         ! Retrieve column surface and soil temperatures from hybrid structure
         do k = 0,zmin,-1
            surf_soil_temp(abs(k) + 1) = sol(S_TEMP, k)%data(d)%elts(id_i)
         end do
      end subroutine retrieve_prog_vars

      subroutine save_tendencies
         ! Description: Saves tendencies calculated by the physics to dynamics hybrid data structure.
         ! Arguments
         integer :: k
         real(8) :: full_mass

         do k= 1, zlevels
            ! updated mass trend - none due to hydrostatic balance and get full_mass
            trend(S_MASS,k)%data(d)%elts(id_i) = 0.0_8
            full_mass = sol(S_MASS, k)%data(d)%elts(id_i) + sol_mean(S_MASS, k)%data(d)%elts(id_i)

            !Convert temperature tendency
            phys_dtheta(k) = phys_dtemp(k)*(dom%surf_press%elts(id_i)/phys_pres_lay(k))**kappa
            trend(S_TEMP,k)%data(d)%elts(id_i) = phys_dtheta(k)*full_mass

            !Save Zonal and Meridional Velocities, will be converted to edges once entire domain finished
            dzonal(k)%data(d)%elts(id_i) = phys_dzonal(k)
            dmerid(k)%data(d)%elts(id_i) = phys_dmeridional(k)
         end do

         ! Save column surface soil temp in temperature hybrid structure
         do k = 0,zmin,-1
            sol(S_TEMP,k)%data(d)%elts(id_i) = surf_soil_temp(abs(k) + 1)
         end do
      end subroutine save_tendencies
   end subroutine physics_call

   subroutine cal_surf_press_phys (q)
      implicit none
      ! Description : Compute surface pressure of all domains and save in domain press_lower element for upward integration
      !               Set geopotential to surface geopotential for upward integration
      !               Orginal
      ! Date Revised: Feb 1st
      type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q

      integer :: d, k, mass_type, p

      ! Set surface geopotential
      call apply (set_surf_geopot_phys, z_null)

      ! Calculate surface pressure
      do d = 1, size(grid)
         grid(d)%surf_press%elts = 0.0_8
         ! Get total mass of the column
         do k = 1, zlevels
            mass   => q(S_MASS,k)%data(d)%elts
            temp   => q(S_TEMP,k)%data(d)%elts
            mean_m => sol_mean(S_MASS,k)%data(d)%elts
            mean_t => sol_mean(S_TEMP,k)%data(d)%elts
            do p = 3, grid(d)%patch%length
               call apply_onescale_to_patch (column_mass_phys, grid(d), p-1, k, 0, 1)
            end do
            nullify (mass, mean_m, mean_t, temp)
         end do
         ! using hydrostatic approx get surface pressure (P-bottom - p_top = -row*g*dz) (recall mass = reference density*dz)
         grid(d)%surf_press%elts = grav_accel*grid(d)%surf_press%elts + p_top

         grid(d)%press_lower%elts = grid(d)%surf_press%elts
      end do
   end subroutine cal_surf_press_phys

   subroutine column_mass_phys (dom, i, j, zlev, offs, dims)
      ! Description: Sum up mass (ie add mass to previous sum) and save in surface pressure array
      implicit none
      type (Domain)                  :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: id_i
      ! Get id of column
      id_i = idx (i, j, offs, dims) + 1

      dom%surf_press%elts(id_i) = dom%surf_press%elts(id_i) + mass(id_i)
   end subroutine column_mass_phys

   subroutine set_surf_geopot_phys (dom, i, j, zlev, offs, dims)
      ! Set initial geopotential to surface geopotential
      implicit none
      type (Domain)                  :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: id_i

      id_i = idx (i, j, offs, dims) + 1

      dom%geopot%elts(id_i) = surf_geopot_case (dom%node%elts(id_i))
   end subroutine set_surf_geopot_phys

   subroutine cal_press_geopot_layer (dom, i, j, zlev, offs, dims)
      ! Integrate pressure up from surface to top layer
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: id_i, d
      real(8) :: full_mass, pressure_upper, pressure_lower, exner

      id_i = idx (i, j, offs, dims) + 1
      d = dom%id + 1

      ! Calculate mass
      full_mass = sol(S_MASS,zlev)%data(d)%elts(id_i) + sol_mean(S_MASS,zlev)%data(d)%elts(id_i)

      ! Get pressure at next zlevel interface
      pressure_lower = dom%press_lower%elts(id_i)
      pressure_upper = pressure_lower - grav_accel * full_mass
      ! Get pressure at the center of current zlevel
      dom%press%elts(id_i) = 0.5 * (pressure_lower + pressure_upper)

      ! Get geopotential at next zlevel interface
      dom%geopot_lower%elts(id_i) = dom%geopot%elts(id_i)
      exner = c_p * (dom%press%elts(id_i)/p_0)**kappa
      dom%geopot%elts(id_i) = dom%geopot_lower%elts(id_i) + &
         grav_accel*kappa*sol(S_TEMP,zlev)%data(d)%elts(id_i)*exner/dom%press%elts(id_i)

      ! Update pressure lower to next interface for next integration
      dom%press_lower%elts(id_i) = pressure_upper
   end subroutine cal_press_geopot_layer

   function z_coords_case (eta_surf, z_s)
      ! Dummy routine
      ! (see upwelling test case for example)
      implicit none
      real(8)                       :: eta_surf, z_s ! free surface and bathymetry
      real(8), dimension(0:zlevels) :: z_coords_case

      z_coords_case = 0.0_8
   end function z_coords_case
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Place holder subroutine for dynamics time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   subroutine timestep_placeholder(align_time, aligned)
      ! Description: Subroutine to call in place of dynamics time step if only calling the physics by itself.
      !              This will update the time and flags needed to know when to checkpoint.
      use main_mod, ONLY : time_mult, dt_new
      implicit none
      real(8)              :: align_time
      logical, intent(out) :: aligned

      integer(8) :: idt, ialign

      istep       = istep+1
      istep_cumul = istep_cumul+1

      dt = dt_new
      !dt = 1200.0_8 !hard_coded case

      !  Check to see if need to write checkpoint
      idt    = nint (dt*time_mult, 8)
      ialign = nint (align_time*time_mult, 8)
      if (ialign > 0 .and. istep > 20) then
         aligned = (modulo (itime+idt,ialign) < modulo (itime,ialign))
         ! For case when dt is exact divisor of align_time:
         !aligned = modulo(time + dt, align_time) < modulo(time,align_time)      !real(8) case
      else
         aligned = .false.
      end if

      ! Update time
      itime = itime + idt
      time = time + dt

   end subroutine timestep_placeholder

   subroutine timestep_placeholder_real(align_time, aligned,dt_real)
      ! Description: Subroutine to call in place of dynamics time step if only calling the physics by itself.
      !              This will update the time and flags needed to know when to checkpoint.
      use main_mod, ONLY : time_mult, dt_new
      implicit none
      !! REAL version

      real              :: align_time, dt_real
      logical, intent(out) :: aligned

      integer :: idt, ialign

      istep       = istep+1
      istep_cumul = istep_cumul+1

      dt = dt_new
      !dt = 1200.0 !hard coding dt

      !  Check to see if need to write checkpoint
      idt    = dt_real*time_mult
      ialign = align_time*time_mult
      if (ialign > 0 .and. istep > 20) then
         ! For case when dt is exact divisor of align_time
         !aligned = modulo(time + dt_real, align_time) < modulo(time,align_time) !real case
         aligned = (modulo (itime+idt,ialign) < modulo (itime,ialign))
      else
         aligned = .false.
      end if

      ! Update time
      itime = itime + idt
      time = time + dt_real

   end subroutine timestep_placeholder_real
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Post Processing Calculations - Mean Values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mean_values (iwrt)
      ! Saves  Temperature, Velocity and Geopotential averages over the sphere
      ! (assumes non-adaptive grid)
      use io_mod
      implicit none
      integer :: iwrt

      integer                       :: k
      real(8)                       :: area, high_lat_area, mid_lat_area, low_lat_area, low_lat_upper_lim, high_lat_lower_lim
      real(8), dimension(1:zlevels) :: dz, T_avg, Total_Tavg, Pressure_avg, Geopot_avg, Zonal_vel_avg, Meridional_avg, &
         Zonal_KE_avg, Merid_KE_avg, Low_Lat_avg, Mid_lat_avg, High_lat_avg
      real(8), dimension(1:zlevels) ::  z
      real(8), dimension(2) :: mid_lat_range
      real(8), dimension(0:Nsoil) :: tsurf_soil_avg
      character(4)                  :: s_time

      ! Set the ranges
      low_lat_upper_lim = (23.5*MATH_PI/180)
      mid_lat_range = (/(23.5*MATH_PI/180), (66.5*MATH_PI/180)/)
      high_lat_lower_lim = (66.5*MATH_PI/180)

      !Set latatide averages and area to zero
      Low_Lat_avg = 0
      Mid_lat_avg = 0
      High_lat_avg = 0
      high_lat_area = 0
      mid_lat_area = 0
      low_lat_area = 0
      tsurf_soil_avg = 0

      ! Calculate surface pressure, as pressure is calculated in temp_fun to get temperature
      call cal_surf_press_phys(sol(1:N_VARIABLE,1:zmax))

      !Calculate Sums of temp, geopotential, and velocities, Kinetic Energies, Temp of Latitude zones
      area = integrate_hex (area_fun, 1, .true.)
      high_lat_area = sum_real(high_lat_area)
      mid_lat_area = sum_real(mid_lat_area)
      low_lat_area = sum_real(low_lat_area)
      do k = 1, zlevels
         T_avg(k)  = integrate_hex (temp_fun, k, .true.) ! This will take care if mpi is used
         Pressure_avg(k) = integrate_hex(pressure_fun, k, .true.)
         Geopot_avg(k) = integrate_hex (geopot_fun, k, .true.)
         Zonal_vel_avg(k) = integrate_hex(zonal_fun, k, .true.)
         Meridional_avg(k) = integrate_hex(merid_fun, k, .true.)
         Zonal_KE_avg(k) = integrate_hex(zonal_KE, k, .true.)
         Merid_KE_avg(k) = integrate_hex(merid_KE, k, .true.)
         ! Temps of zones on each rank summed in temp_fun call for T_avg
         High_lat_avg(k) = sum_real(High_lat_avg(k)) ! Get temp of high lat zones from all ranks
         Mid_lat_avg(k) = sum_real(Mid_lat_avg(k)) ! Get temp of mid lat zones from all ranks
         Low_Lat_avg(k) = sum_real(Low_Lat_avg(k)) ! Get temp of low lat zones from all ranks
      end do
      ! Calculate sums of the surface temp and soil temps
      do k = 0,zmin,-1
         tsurf_soil_avg(abs(k)) = integrate_hex(surf_soil_temp_fun,k,.true.)
      end do

      if (rank == 0) then
         ! Get average
         Total_Tavg  = T_avg  / area
         Pressure_avg = Pressure_avg / area
         Geopot_avg = Geopot_avg / area
         Zonal_vel_avg = Zonal_vel_avg / area
         Meridional_avg = Meridional_avg / area
         Zonal_KE_avg = 0.5 * Zonal_KE_avg / area
         Merid_KE_avg = 0.5 * Merid_KE_avg / area
         Low_Lat_avg = Low_Lat_avg / low_lat_area
         Mid_lat_avg = Mid_lat_avg / mid_lat_area
         High_lat_avg = High_lat_avg / high_lat_area
         tsurf_soil_avg = tsurf_soil_avg / area

         ! Write values to terminal and file
         write (s_time, '(i4.4)') iwrt
         open (unit=20, file=trim(run_id)//'.6.'//s_time, form="FORMATTED", action='WRITE', status='REPLACE')

         write (6,'(a, f4.1, a)') "Temperature profile at time ", time/HOUR, " h"
         write (6,'(a)') 'Level    z_k     pres_k     geopot_k     T(z_k)       u(z_k)       v(z_k)        1/2u(z_k)^2 '// &
            '  1/2v(z_k)^2  Low_lat_T(z_k) Mid_lat_T(z_k) High_lat_T(z_k)'

         write(20,*) time ! write the time
         do k = 1, zlevels
            z(k) = log(Pressure_avg(k)/p_0)*(-R_d*Total_Tavg(k)/grav_accel)
            write (6,'(i3, 3x, f9.2, 1x, f9.2, 1x, f9.2, 9(1x, es13.4))') &
               k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), Zonal_vel_avg(k), Meridional_avg(k),&
               Zonal_KE_avg(k), Merid_KE_avg(k), Low_Lat_avg(k), Mid_lat_avg(k), High_lat_avg(k)
            ! write (20,'(i3, 1x, 6(es13.6,1x))') k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), &
            !    Zonal_vel_avg(k), Meridional_avg(k)
            write (20,*) k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), &
               Zonal_vel_avg(k), Meridional_avg(k), Zonal_KE_avg(k), Merid_KE_avg(k), Low_Lat_avg(k), Mid_lat_avg(k),&
               High_lat_avg(k)
         end do
         ! Write the surface temperature and soil temp (if soil included)
         do k = 0,zmin,-1
            write(20,*) k, tsurf_soil_avg(abs(k))
         end do

         close (20)
      end if
   contains
      real(8) function temp_fun (dom, i, j, zlev, offs, dims)
         ! Calculates the temperature at center of zlev layer of an element in a domain
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: d, id_i
         real(8) :: full_mass, full_theta, potential_temp, temperature, lat, long

         d = dom%id + 1
         id_i = idx (i, j, offs, dims) + 1
         call cart2sph(dom%node%elts(id_i), long, lat)
         lat = abs(lat)

         ! Calculate the pressure
         call cal_press_geopot_layer(dom, i, j, zlev, offs, dims)

         full_mass  = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + sol(S_MASS,zlev)%data(d)%elts(id_i)
         full_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) + sol(S_TEMP,zlev)%data(d)%elts(id_i)
         potential_temp = full_theta/full_mass

         ! Calculate the temperature at the center of zlayer
         temperature = potential_temp * ((dom%press%elts(id_i)/dom%surf_press%elts(id_i))**kappa)

         ! Gather the zonal lat temps
         if (lat .ge. mid_lat_range(1) .and. lat .le. mid_lat_range(2)) then
            Mid_lat_avg(zlev) = Mid_lat_avg(zlev) + (temperature/dom%areas%elts(id_i)%hex_inv)
         end if
         if (lat .ge. high_lat_lower_lim) High_lat_avg(zlev) = High_lat_avg(zlev) + (temperature/dom%areas%elts(id_i)%hex_inv)
         if (lat .le. low_lat_upper_lim) Low_Lat_avg(zlev) = Low_lat_avg(zlev) + (temperature/dom%areas%elts(id_i)%hex_inv)

         ! Save density in Ke dom arg for KE
         dom%ke%elts(id_i) = dom%press%elts(id_i)/(temperature*R_d)

         !Return temp
         temp_fun = temperature

      end function temp_fun

      real(8) function zonal_fun(dom, i, j, zlev, offs, dims)
         ! Calculates the velocity of an element at the center of a zlev layer in a domain, but returns zonal velocity
         ! meridional velocity is saved in the dom%v_merid%elts for when meridional integration is called
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: d, id_i

         d = dom%id + 1
         id_i = idx (i, j, offs, dims) + 1

         velo => sol(S_VELO,zlev)%data(d)%elts
         velo1 => grid(d)%u_zonal%elts
         velo2 => grid(d)%v_merid%elts
         call interp_UVW_latlon(dom, i, j, zlev, offs, dims)
         zonal_fun = velo1(id_i)
         nullify (velo, velo1, velo2)
      end function zonal_fun

      real(8) function zonal_KE(dom, i, j, zlev, offs, dims)
         ! Calculates the Zonal Kinetic energy of an element at the center of a zlev layer in a domain
         ! Assumes that the velocity was already calculated and stored in dom%u_zonal, from zonal velocity call
         ! Assumes that the density was already calculated and stored in dom%ke, from the temp call
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: d, id_i

         d = dom%id + 1
         id_i = idx (i, j, offs, dims) + 1

         zonal_KE = dom%ke%elts(id_i)*((dom%u_zonal%elts(id_i))**2)

      end function zonal_KE

      real(8) function merid_fun(dom, i, j, zlev, offs, dims)
         ! Calculates the meridional velocity of an element at the center of a zlev layer in a domain
         ! Assumes that the velocity was already calculated and stored in dom%v_merid, from zonal velocity call
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i

         id_i = idx (i, j, offs, dims) + 1

         merid_fun = dom%v_merid%elts(id_i)

      end function merid_fun

      real(8) function merid_KE(dom, i, j, zlev, offs, dims)
         ! Calculates the Meridional kinetic energy of an element at the center of a zlev layer in a domain
         ! Assumes that the velocity was already calculated and stored in dom%v_merid, from zonal velocity call
         ! Assumes that the density was already calculated and stored in dom%ke, from the temp call
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i

         id_i = idx (i, j, offs, dims) + 1

         merid_KE = dom%ke%elts(id_i)*((dom%v_merid%elts(id_i))**2)

      end function merid_KE

      real(8) function pressure_fun (dom, i, j, zlev, offs, dims)
         ! Calculates the pressure of an element at the center of a zlev layer in a domain
         ! Assumes pressure at the zlevel interfaces are already calculated and stored in dom%press
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i

         id_i = idx (i, j, offs, dims) + 1

         pressure_fun = dom%press%elts(id_i)

      end function pressure_fun
      real(8) function geopot_fun (dom, i, j, zlev, offs, dims)
         ! Calculates the geopotential of an element at the center of a zlev layer in a domain
         ! Assumes geopot at the zlevel interfaces above and below already calculated and stored in dom%geopot and dom%geopot_lower
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i

         id_i = idx (i, j, offs, dims) + 1

         geopot_fun = 0.5*(dom%geopot_lower%elts(id_i)+dom%geopot%elts(id_i))

      end function geopot_fun

      real(8) function area_fun (dom, i, j, zlev, offs, dims)
         ! Defines mass for total mass integration
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i
         real(8) :: lat, long

         id_i = idx (i, j, offs, dims) + 1
         call cart2sph(dom%node%elts(id_i), long, lat)
         lat = abs(lat)

         if (lat .ge. mid_lat_range(1) .and. lat .le. mid_lat_range(2)) then
            mid_lat_area = mid_lat_area + (1/dom%areas%elts(id_i)%hex_inv)
         end if
         if (lat .ge. high_lat_lower_lim) high_lat_area = high_lat_area + (1/dom%areas%elts(id_i)%hex_inv)
         if (lat .le. low_lat_upper_lim) low_lat_area = low_lat_area + (1/dom%areas%elts(id_i)%hex_inv)

         area_fun = 1.0_8

      end function area_fun

      real(8) function surf_soil_temp_fun(dom, i, j, zlev, offs, dims)
         ! Get the temperature of the surface or the soil at zlev
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i,d

         id_i = idx (i, j, offs, dims) + 1
         d = dom%id + 1

         surf_soil_temp_fun = sol(S_TEMP,zlev)%data(d)%elts(id_i)

      end function surf_soil_temp_fun

   end subroutine mean_values

   subroutine get_coordinates()
      ! Description: Subroutine to retrieve all latitude and longitude coordinates and save in a file named: run_id_coordinates

      !Arguments
      integer   :: d, p
      real(8)  :: lat, long
      type(Coord) :: x_i

      ! Open a file to write the latitude and longitude
      open(90, file=trim(run_id)//'_coordinates', form="FORMATTED", action='WRITE', status='REPLACE')

      !Retrieve longitude and latitude elements of each domain and combine into single array
      do d = 1, size(grid)
         do p = 3, grid(d)%patch%length
            call apply_onescale_to_patch(get_lat_long, grid(d), p-1, z_null, 0, 0)
         end do
      end do

      close(90)
   end subroutine

   subroutine get_lat_long(dom, i, j, zlev, offs, dims)

      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims
      integer :: id_i
      real(8) :: lat, long

      id_i = idx (i, j, offs, dims) + 1

      call cart2sph(dom%node%elts(id_i), long, lat)

      write(90, '(2(es13.6,1x))') lat, long

   end subroutine get_lat_long


end module test_case_mod
