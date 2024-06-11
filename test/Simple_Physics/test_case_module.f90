module test_case_mod
   ! Description: Test case module for the Simple physics
   ! Author: Gabrielle Ching-Johnson
   ! Date Revised: November 1 2023
   use comm_mpi_mod
   use utils_mod
   use init_mod
   implicit none

   ! Standard variables
   integer :: CP_EVERY, resume_init
   real(8) :: dt_cfl, total_cpu_time, dPdim, R_ddim, specvoldim, dTempdim

   ! Test case variables
   real(8) :: T_0 ! Reference Temperature
   real(8) :: u_0, e_thick ! geostrophic wind speed, ekman layer thickness

contains
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Dynamics test case routines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

   subroutine init_sol (dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Initialization routine to be called for each zlev and grid point.
      !                Initializes the initial conditions:
      !                    - potenital temp -> mass weighted potential temp
      !                    - pressure at the layer center
      !                    - Intepolates velocites from nodes to edges
      !
      !-----------------------------------------------------------------------------------
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
      !-----------------------------------------------------------------------------------
      !
      !   Description: Initialize surface pressure, geoptentials and velocities
      !                (zonal and meridional) at the node.
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
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

   real(8) function surf_geopot_case (d, id)
      ! Surface geopotential set to zero
      implicit none
      integer :: d, id

      ! Geopotential at sea level across sphere
      surf_geopot_case = 0d0
   end function surf_geopot_case

   real(8) function surf_pressure (x_i)
      ! Surface pressure
      implicit none
      type(Coord) :: x_i

      surf_pressure = p_0
   end function surf_pressure

   subroutine geopot_initialize(dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Initialize geopotential at element node
      !
      !   Assumption: Surface pressure is set and grid has been set
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims

      real(8) :: p, p_s, phi
      integer :: id_i

      id_i = idx(i, j, offs, dims) + 1

      p_s = dom%surf_press%elts(id_i)
      ! Calculate pressure at the center of vertical layer
      p = 0.5 * (a_vert(zlev)+a_vert(zlev+1) + (b_vert(zlev)+b_vert(zlev+1))*p_s)

      ! Geopotential
      phi = -R_d*T_0*LOG((p)/p_s)
      dom%geopot%elts(id_i) = phi

   end subroutine geopot_initialize

   subroutine vel_initialize (dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Initialize velocities (zonal and meridional) at the element node
      !                 using the theory of the ekman layer.
      !
      !   Assumption: Geopotential at the element has been set in dom%geopot%elts
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
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
            ! Set velocities at the nodes for all domains
            call apply_onescale(init_vels_and_surf, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
            ! Initialize prognostic variables
            call apply_onescale (init_sol, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
         end do
      end do
   end subroutine apply_initial_conditions_case

   subroutine initialize_a_b_vert_case
      !-----------------------------------------------------------------------------------
      !
      !   Description: Initialize the vertical grid using equal pressure layers.
      !
      !-----------------------------------------------------------------------------------
      implicit none
      integer :: k

      ! Allocate vertical grid parameters
      allocate (a_vert(1:zlevels+1),    b_vert(1:zlevels+1))
      allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

      do k = 1, zlevels+1
         a_vert(k) = dble(k-1)/dble(zlevels) * p_top
         b_vert(k) = 1.0_8 - dble(k-1)/dble(zlevels)
      end do

      ! Set mass coefficients
      a_vert_mass = (a_vert(1:zlevels) - a_vert(2:zlevels+1))/grav_accel
      b_vert_mass =  b_vert(1:zlevels) - b_vert(2:zlevels+1)
   end subroutine initialize_a_b_vert_case

   subroutine initialize_thresholds_case
      ! Set default thresholds based on dimensional scalings of norms
      implicit none
      integer :: k

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

   function z_coords_case (eta_surf, z_s)
      ! Dummy routine
      ! (see upwelling test case for example)
      implicit none
      real(8)                       :: eta_surf, z_s ! free surface and bathymetry
      real(8), dimension(0:zlevels) :: z_coords_case

      z_coords_case = 0.0_8
   end function z_coords_case

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
         call cal_lnorm (order)
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
            call system (trim(command))

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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Place holder subroutine for dynamics time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine timestep_placeholder(align_time, aligned)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Subroutine to call in place of dynamics time step if only calling
      !                  the physics by itself. This will update the time and flags
      !                  needed to know when to checkpoint.
      !
      !   Key Note: This is the REAL(8) case, used for physics only cases.
      !             It includes hard coded dt case, commented out if desired
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
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
      !-----------------------------------------------------------------------------------
      !
      !   Description: Subroutine to call in place of dynamics time step if only calling
      !                  the physics by itself. This will update the time and flags
      !                  needed to know when to checkpoint.
      !
      !   Key Note: This is the REAL case, used for testing and physics only cases.
      !             It includes hard coded dt case, commented out if desired
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module test_case_mod
