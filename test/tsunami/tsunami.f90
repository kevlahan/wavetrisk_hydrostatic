program Tsunami
  ! Simple tsunami test case
  !
  ! The permeability penalization parameter eta = dt_cfl and the porosity parameter alpha is set in the input file.
  use main_mod
  use test_case_mod
  use io_mod
  use lin_solve_mod
  implicit none
  integer                                :: N
  integer, dimension(2)                  :: Nx, Ny
  real(4), dimension(:,:),   allocatable :: field2d
  real(8)                                :: dx_export, dy_export, kx_export, ky_export
  real(8), dimension(2)                  :: lon_lat_range
  real(8), dimension(:,:,:), allocatable :: field2d_save
  logical                                :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6371.229   * KM              ! mean radius of the Earth
  grav_accel     = 9.80616    * METRE/SECOND**2 ! gravitational acceleration 
  omega          = 7.29211d-5 * RAD/SECOND      ! Earth’s angular velocity 
  p_top          = 0.0_8      * hPa             ! pressure at free surface
  ref_density    = 1027       * KG/METRE**3     ! reference density (seawater)

  ! Local test case parameters
  min_depth      = -50  * METRE                 ! minimum allowed depth (must be negative)
  max_depth      = -4   * KM                    ! maximum allowed depth (must be negative)

  dH             =  7   * METRE                 ! initial perturbation to the free surface
  pert_radius    =  1e3 * KM                    ! radius of Gaussian free surface perturbation
  lon_c          = -50  * DEG                   ! longitude location of perturbation
  lat_c          =  25  * DEG                   ! latitude  location of perturbation

  ! Dimensional scaling
  wave_speed     = sqrt (grav_accel*abs(max_depth)) ! inertia-gravity wave speed based on maximum allowed depth
  Udim           = wave_speed                   ! velocity scale
  Ldim           = 2 * pert_radius              ! length scale (free surface perturbation width)
  Tdim           = Ldim/Udim                    ! time scale
  Hdim           = abs (max_depth)              ! vertical length scale

  ! Parameters for 2D projection
  N              = 1024                         ! size of lat-lon grid in 2D projection
  lon_lat_range  = (/2*MATH_PI, MATH_PI/)       ! region to save in 2D projection

  mode_split     = .false.                      ! use explicit time step for accuracy
  compressible   = .false.                      ! incompressible
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize 2D projection grid
  call init_2d_grid

  ! Read test case parameters
  call read_test_case_parameters

  ! Initialize variables
  call initialize (run_id)

  ! Initialize tsunami diagnostic variables
  call init_diagnostics

  ! Save initial conditions
  call print_test_case_parameters
  call write_and_export (iwrite)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  total_cpu_time = 0.0_8
  do while (time < time_end)
     call start_timing;  call time_step (dt_write, aligned); call stop_timing

     call update_diagnostics
     
     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) then
           call deallocate_diagnostics
           call write_checkpoint (run_id, rebalance)
           call init_diagnostics !! resets diagnostics !!
        end if

        ! Save fields
        call write_and_export (iwrite)
        call save_diagnostics
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
contains
  subroutine save_diagnostics 
    ! Interpolate diagnostics onto lon-lat and save
    use test_case_mod
    implicit none
    
    integer            :: fid, i, v
    integer, parameter :: funit = 400
    character(130)     :: command
    character(4)       :: s_time
    character(7)       :: var_file

    ! Project diagnostics onto latitude-longitude plane
    call update_bdry (max_wave_height, min_level, 80)
    call update_bdry (arrival_time,    min_level, 81)

    call project_onto_plane (max_wave_height, min_level, 0.0_8)
    field2d_save(:,:,1) = field2d

    call project_onto_plane (arrival_time,    min_level, 0.0_8)
    field2d_save(:,:,2) = field2d

    ! Write out diagnostics
    if (rank == 0) then
       do v = 1, 2
          fid = 6000000+100*iwrite + 10 + v
          write (var_file, '(I7)') fid
          open (unit=funit, file=trim(run_id)//'.'//var_file)
          do i = Ny(1), Ny(2)
             write (funit,'(2047(E15.6, 1X))') field2d_save(:,i,v)
          end do
          close (funit)
       end do

       ! Coordinates

       ! Longitude values
       fid = 6000000+100*iwrite + 20
       write (var_file, '(I7)') fid
       open (unit=funit, file=trim(run_id)//'.'//var_file)
       write (funit,'(2047(E15.6, 1X))') (-180+dx_export*(i-1)/MATH_PI*180, i=1,Nx(2)-Nx(1)+1)
       close (funit)

       ! Latitude values
       fid = 6000000+100*iwrite + 21
       write (var_file, '(I7)') fid
       open (unit=funit, file=trim(run_id)//'.'//var_file)
       write (funit,'(2047(E15.6, 1X))') (-90+dy_export*(i-1)/MATH_PI*180, i=1,Ny(2)-Ny(1)+1)
       close (funit)

       ! Compress files
       write (s_time, '(i4.4)') iwrite
       command = 'ls -1 '//trim(run_id)//'.6'//s_time//'?? > tmp3'
       call system (command)
       command = 'tar cfz '//trim(run_id)//'.6'//s_time//'.tgz -T tmp3 --remove-files'
       call system (command)
    end if
  end subroutine save_diagnostics

  subroutine init_2d_grid
    implicit none

    Nx = (/-N/2, N/2/)
    Ny = (/-N/4, N/4/)

    dx_export = lon_lat_range(1)/(Nx(2)-Nx(1)+1); dy_export = lon_lat_range(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1.0_8/dx_export; ky_export = 1.0_8/dy_export
    
    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))
    allocate (field2d_save(Nx(1):Nx(2),Ny(1):Ny(2),1:2))
  end subroutine init_2d_grid

  subroutine project_onto_plane (field, l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    implicit none
    integer               :: l, itype
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
  end subroutine project_onto_plane

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
    if (maxx-minx > MATH_PI/2) then
       write (0,'(A,i4,A)') 'ERROR (rank = ', rank, '): io-333 "export"'
       return
    end if

    do id_x = floor (kx_export*minx), ceiling (kx_export*maxx)
       if (id_x < lbound (field2d,1) .or. id_x > ubound (field2d,1)) cycle
       do id_y = floor (ky_export*miny), ceiling (ky_export*maxy)
          if (id_y < lbound (field2d,2) .or. id_y > ubound (field2d,2)) cycle
          ll = (/dx_export*id_x, dy_export*id_y/)
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
       a(1) = a(1) - sum(fixed) * 2*MATH_PI
       b(1) = b(1) - sum(fixed) * 2*MATH_PI
       c(1) = c(1) - sum(fixed) * 2*MATH_PI
       call interp_tri_to_2d (a, b, c, val)
    end if
  end subroutine interp_tri_to_2d_and_fix_bdry
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics routines for this test case (including diffusion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function physics_scalar_flux (q, dom, id, idE, idNE, idN, v, zlev, type)
  ! Additional physics for the flux term of the scalar trend
  ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
  !
  ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
  use domain_mod
  implicit none

  real(8), dimension(1:EDGE)                           :: physics_scalar_flux
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
     physics_scalar_flux = 0.0_8
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
     physics_scalar_flux = (-1)**Laplace_order * visc_sclr(v) * grad * l_e
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
end function physics_scalar_flux

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  use domain_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                    :: d, id, id_i, visc_scale
  real(8), dimension(1:EDGE) :: diffusion

  d = dom%id + 1
  id = idx (i, j, offs, dims)
  id_i = id + 1

  visc_scale = 1.0_8!(max_level/dom%level%elts(id+1))**(2*Laplace_order_init-1)

  if (Laplace_order == 0) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu()) * visc_scale
  end if

  ! Total physics for source term of velocity trend including volume penalization
  physics_velo_source = diffusion 
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
end function physics_velo_source





