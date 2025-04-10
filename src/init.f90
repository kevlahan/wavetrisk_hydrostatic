module init_mod
  use geom_mod
  use domain_mod
  use arch_mod
  implicit none
  real(8), parameter :: YANGLE = 0d0
  
  abstract interface
     real(8) function fun1 (eta, ri, z)
       implicit none
       real(8) :: eta, ri, z
     end function fun1
     function fun2 (eta, ri, z)
       use shared_mod
       implicit none
       real(8), dimension(1:EDGE) :: fun2, eta, z
       real(8)                    :: ri
     end function fun2
     real(8) function fun3 (dom, i, j, z_lev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, z_lev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end function fun3
     function fun4 (dom, i, j, z_lev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, z_lev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
       real(8), dimension(1:EDGE)     :: fun4
     end function fun4
     subroutine sub4 (dom, i, j, z_lev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, z_lev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end subroutine sub4
     function fun5 (q, dom, id, idE, idNE, idN, v, zlev, type)
       use domain_mod
       implicit none
       real(8), dimension(1:EDGE)                           :: fun5
       type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
       type(domain)                                         :: dom
       integer                                              :: d, id, idE, idNE, idN, v, zlev
       logical, optional                                    :: type
     end function fun5
     real(8) function coord_fun (p)
       use geom_mod
       implicit none
       type(Coord) :: p
     end function coord_fun
     real(8) function int2_fun (d, id)
       implicit none
       integer :: d, id
     end function int2_fun
     subroutine io_fun (fid)
       implicit none
       integer :: fid
     end subroutine io_fun
     subroutine noarg_fun
       implicit none
       integer :: fid
     end subroutine noarg_fun
      function zcoords_fun (eta_surf, z_s)
        use shared_mod
        implicit none
        real(8)                       :: eta_surf, z_s 
        real(8), dimension(0:zlevels) :: zcoords_fun
      end function zcoords_fun
      subroutine solver (u, f, Lu, Lu_diag)
        use domain_mod
        implicit none
        type(Float_Field), intent(in)    :: f
        type(Float_Field), intent(inout) :: u
        interface
           function Lu (u, l)
             use domain_mod
             implicit none
             integer                   :: l
             type(Float_Field), target :: Lu, u
           end function Lu
           function Lu_diag (u, l)
             use domain_mod
             implicit none
             integer                   :: l
             type(Float_Field), target :: Lu_diag, u
           end function Lu_diag
        end interface
      end subroutine solver
      subroutine physics_fun (q, dq)
        use domain_mod
        implicit none
        type(Float_Field), dimension(1:N_VARIABLE,zmin:zmax), target :: q
        type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: dq
      end subroutine physics_fun
   end interface

  ! General procedure pointers
  procedure (sub4),        pointer :: u_source                 => null ()

  ! Pointers to variables and procedures that may be defined in test cases
  real(8),                 pointer :: bottom_friction
  procedure (noarg_fun),   pointer :: apply_initial_conditions => null ()
  procedure (fun3),        pointer :: bottom_buoy_flux         => null ()
  procedure (io_fun),      pointer :: dump                     => null ()
  procedure (io_fun),      pointer :: load                     => null ()
  procedure (noarg_fun),   pointer :: initialize_a_b_vert      => null ()
  procedure (noarg_fun),   pointer :: initialize_dt_viscosity  => null ()
  procedure (noarg_fun),   pointer :: initialize_thresholds    => null ()
  procedure (noarg_fun),   pointer :: set_thresholds           => null ()
  procedure (int2_fun),    pointer :: surf_geopot              => null ()
  procedure (coord_fun),   pointer :: tau_mag                  => null ()
  procedure (fun3),        pointer :: top_buoy_flux            => null ()
  procedure (noarg_fun),   pointer :: update                   => null ()
  procedure (fun4),        pointer :: wind_flux                => null ()
  procedure (fun4),        pointer :: physics_velo_source      => null ()
  procedure (fun5),        pointer :: physics_scalar_flux      => null ()
  procedure (zcoords_fun), pointer :: z_coords                 => null ()

  ! Elliptic solver 
  procedure (solver),      pointer :: elliptic_solver          => null ()

  ! Physics trend
  procedure (physics_fun), pointer :: trend_physics            => null ()
contains
  subroutine init_init_mod
    implicit none
    logical :: initialized = .false.
    
    if (initialized) return ! initialize only once
    call init_domain_mod
    call init_arch_mod
    initialized = .true.
  end subroutine init_init_mod

  subroutine init_grid
    implicit none
    integer :: b, d, i_, k, loz, p, s, v

    if (mode_split) then ! separate free surface level at zlev = zlevels+1
       zmax = zlevels+1
    else
       zmax = zlevels
    end if

    if (ref_density == 0d0) then ! ref_density not set in test case: choose correct default value
       if (compressible) then
          ref_density = ref_density_air
       else
          ref_density = ref_density_water
       end if
    end if

    allocate (grid(n_domain(rank+1)))
    allocate (sol(1:N_VARIABLE,zmin:zmax), sol_mean(1:N_VARIABLE,zmin:zmax), trend(1:N_VARIABLE,1:zmax))
    allocate (wav_coeff(1:N_VARIABLE,zmin:zmax))
    allocate (exner_fun(zmin:zmax+1))
    allocate (penal_node(zmin:zmax), penal_edge(zmin:zmax))
    allocate (horiz_flux(scalars(1):scalars(2)), Laplacian_scalar(scalars(1):scalars(2)))
    allocate (Laplacian_vector(S_DIVU:S_ROTU))
    allocate (lnorm(1:N_VARIABLE,zmin:zmax))
    if (vert_diffuse) allocate (Kt(0:zlevels), Kv(0:zlevels), tke(1:zlevels), wav_tke(1:zlevels))

    call init_Float_Field (topography, AT_NODE)
    if (sso) then
       do k = 1, 4
          call init_Float_Field (sso_param(k), AT_NODE)
       end do
    end if
    
    do k = zmin, zmax
       call init_Float_Field (penal_node(k), AT_NODE)
       call init_Float_Field (penal_edge(k), AT_EDGE)
       call init_Float_Field (exner_fun(k),  AT_NODE)
       do v = 1, N_VARIABLE
          call init_Float_Field (sol(v,k),      POSIT(v))
          call init_Float_Field (sol_mean(v,k), POSIT(v))
          if (k > 0) call init_Float_Field (trend(v,k), POSIT(v))
       end do
    end do
    call init_Float_Field (exner_fun(zmax+1), AT_NODE)
    
    if (vert_diffuse) then
       call init_Float_Field (Kv(0), AT_NODE)
       call init_Float_Field (Kt(0), AT_NODE)
       do k = 1, zlevels
          call init_Float_Field (Kv(k),  AT_NODE)
          call init_Float_Field (Kt(k),  AT_NODE)
          call init_Float_Field (tke(k), AT_NODE)
       end do
    end if

    call init_Float_Field (Laplacian_vector(S_DIVU), AT_NODE)
    call init_Float_Field (Laplacian_vector(S_ROTU), AT_EDGE)
    do v = scalars(1), scalars(2)
       call init_Float_Field (horiz_flux(v),       AT_EDGE)
       call init_Float_Field (Laplacian_scalar(v), AT_NODE)
    end do

    do d = 1, n_domain(rank+1)
       call init_Domain (grid(d))

       do k = zmin, zmax
          do v = scalars(1), scalars(2)
             call init (sol(v,k)%data(d),      1)
             call init (sol_mean(v,k)%data(d), 1)
          end do
          call init (sol(S_VELO,k)%data(d),      EDGE)
          call init (sol_mean(S_VELO,k)%data(d), EDGE)
       end do
    end do

    !  Initializes grid for icosahedron
    do d = 1, n_domain(rank+1)

       grid(d)%id = d-1

       do s = 1, N_BDRY
          b = add_bdry_patch_Domain (grid(d), s)
       end do

       grid(d)%penta = .false.

       p = add_patch_Domain (grid(d), min_level-1)

       grid(d)%patch%elts(p+1)%children = 0
       grid(d)%patch%elts(p+1)%neigh    = (/ ( i_ , i_ = -1, -N_BDRY, -1 ) /)
    end do

    do loz = 1, N_ICOSAH_LOZENGE
       call set_penta (N_SUB_DOM*(loz - 1),                                 SOUTHWEST)
       call set_penta (N_SUB_DOM*(loz - 1) + N_SUB_DOM_PER_DIM - 1,         SOUTHEAST)
       call set_penta (N_SUB_DOM*(loz - 1) + N_SUB_DOM - N_SUB_DOM_PER_DIM, NORTHWEST)
       call set_penta (N_SUB_DOM*(loz - 1) + N_SUB_DOM - 1,                 NORTHEAST)
    end do

    call init_connections
    call init_coordinates
  end subroutine init_grid

  subroutine init_coordinates
    implicit none
    integer                     :: d, d_glo, i, ii, j, jj, k, loz
    real(8), dimension(4)       :: lat
    real(8), dimension(10)      :: lon
    type(Coord)                 :: ne, se, sw, nw
    type(Coord), dimension(2,2) :: cnr

    lat = (/ - MATH_PI/2d0, - atan (1d0/2d0), atan (1d0/2d0), MATH_PI/2d0 /)
    lon = (/ ((MATH_PI * dble(k))/5d0, k = 0, 10-1) /)

    do ii = 1, 2
       do jj = 1, 5
          loz = 5*ii - 5 + (jj-1)

          ne = sph2cart (lon(modulo(ii + 2*jj - 2, 10) + 1), lat(ii+1)) / radius
          se = sph2cart (lon(ii+2*jj-2),                     lat(ii)  ) / radius
          sw = sph2cart (lon(modulo(ii + 2*jj - 4, 10) + 1), lat(ii+1)) / radius
          nw = sph2cart (lon(ii+2*jj-2),                     lat(ii+2)) / radius

          call yrotate (nw, cnr(1,2), YANGLE); call yrotate (ne, cnr(2,2), YANGLE)
          call yrotate (sw, cnr(1,1), YANGLE); call yrotate (se, cnr(2,1), YANGLE)

          do j = 1, N_SUB_DOM_PER_DIM
             do i = 1, N_SUB_DOM_PER_DIM
                d_glo = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j-1) + (i-1)

                if (.not. owner(d_glo+1) == rank) cycle

                d = loc_id(d_glo+1)
                call assign_coord (grid(d+1), 1, &
                     get_J0_coord(i,   j,   DOMAIN_LEVEL), &
                     get_J0_coord(i,   j-1, DOMAIN_LEVEL), &
                     get_J0_coord(i-1, j-1, DOMAIN_LEVEL), &
                     get_J0_coord(i-1, j,   DOMAIN_LEVEL))
             end do
          end do
       end do
    end do
  contains
    subroutine  yrotate (c_in, c_out, angle)
      implicit none
      real(8),     intent(in)  :: angle
      type(Coord), intent(in)  :: c_in
      type(Coord), intent(out) :: c_out

      c_out%x =  c_in%x*cos(angle) - c_in%z*sin(angle)
      c_out%y =  c_in%y
      c_out%z =  c_in%x*sin(angle) + c_in%z*cos(angle)
    end subroutine yrotate

    type(Coord) recursive function get_J0_coord (i, j, l) result(c)
      implicit none
      integer :: i, j, l

      if (l > 0) then
         c = mid_pt (get_J0_coord(i/2, j/2, l-1), get_J0_coord((i+1)/2, (j+1)/2, l-1))
         return
      else
         c = cnr (i+1,j+1)
         return
      end if
    end function get_J0_coord

  end subroutine init_coordinates

  subroutine ccentre_penta (dom, p)
    implicit none
    type(Domain) :: dom
    integer      :: p
    
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    type(Coord)                    :: pt3
    integer                        :: i, id, idN, idNE, idW, idE, idS, id_mm, id_0m1,id_m10, id_0m2,  id_m20, id_m12, id_m21, j

    call get_offs_Domain (dom, p, offs, dims)

    if (is_penta(dom, p, IJMINUS - 1)) then
       id_mm  = idx (-1, -1, offs, dims)
       id_0m1 = idx ( 0, -1, offs, dims)
       id_m10 = idx (-1,  0, offs, dims)

       dom%ccentre%elts(TRIAG*id_mm+LORT+1) = &
            circumcentre(dom%node%elts(idx(0,0,offs,dims)+1), dom%node%elts(id_0m1+1), dom%node%elts(id_m10+1))

       dom%ccentre%elts(TRIAG*id_mm+UPLT+1) = dom%ccentre%elts(TRIAG*id_mm+LORT+1)

       id_0m2 = idx( 0, -2, offs, dims)
       id_m20 = idx(-2,  0, offs, dims)

       pt3 = mid_pt (dom%node%elts(id_0m2+1), dom%node%elts(id_m20+1))

       id_m12 = idx(-1, -2, offs, dims)
       id_m21 = idx(-2, -1, offs, dims)

       dom%ccentre%elts(TRIAG*id_m12+LORT+1) = circumcentre (dom%node%elts(id_0m1+1), dom%node%elts(id_0m2+1), pt3)
       dom%ccentre%elts(TRIAG*id_m21+UPLT+1) = circumcentre (pt3, dom%node%elts(id_m20+1), dom%node%elts(id_m10+1))
    end if

    if (is_penta(dom, p, IPLUSJMINUS - 1)) then
       i = PATCH_SIZE
       j = -1

       id   = idx(i,   j,   offs, dims)
       idN  = idx(i,   j+1, offs, dims)
       idNE = idx(i+1, j+1, offs, dims)
       idW  = idx(i-1, j,   offs, dims)

       dom%ccentre%elts(TRIAG*idW+LORT+1) = circumcentre(dom%node%elts(idW+1), dom%node%elts(idN+1), dom%node%elts(idNE+1))
       dom%ccentre%elts(TRIAG*id +UPLT+1) = dom%ccentre%elts(TRIAG*idW+LORT+1)
    end if

    if (is_penta(dom, p, IMINUSJPLUS - 1)) then
       i = -1
       j = PATCH_SIZE

       id   = idx (i,   j,   offs, dims)
       idE  = idx (i+1, j,   offs, dims)
       idNE = idx (i+1, j+1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       dom%ccentre%elts(TRIAG*idS+UPLT+1) = circumcentre (dom%node%elts(idS+1), dom%node%elts(idNE+1), dom%node%elts(idE+1))
       dom%ccentre%elts(TRIAG*id+LORT+1) = dom%ccentre%elts(TRIAG*idS+UPLT+1)
    end if

    if (is_penta(dom, p, IJPLUS - 1)) then
       i = PATCH_SIZE
       j = PATCH_SIZE

       id  = idx (i,   j,   offs, dims)
       idN = idx (i,   j+1, offs, dims)
       idE = idx (i+1, j,   offs, dims)

       dom%ccentre%elts(TRIAG*id+LORT+1) = circumcentre (dom%node%elts(id+1), dom%node%elts(idN+1), dom%node%elts(idE+1))
       dom%ccentre%elts(TRIAG*id+UPLT+1) = dom%ccentre%elts(TRIAG*id+LORT+1)
    end if
  end subroutine ccentre_penta

  subroutine lengths (dom, p, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        ::  p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id, idS, idW, idN, idE, idNE

    id   = idx (i,   j,   offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    if (j >= PATCH_SIZE + 1) then
       dom%len%elts(EDGE*id+RT+1) = dist (dom%node%elts(id+1), dom%node%elts(idE+1))
       return
    end if

    if (i >= PATCH_SIZE + 1) then
       dom%len%elts(EDGE*id+UP+1) = dist (dom%node%elts(id+1), dom%node%elts(idN+1))
       return
    end if

    dom%len%elts(EDGE*id+RT+1)    = dist (dom%node%elts(id+1),   dom%node%elts(idE+1))
    dom%len%elts(EDGE*id+DG+1)    = dist (dom%node%elts(idNE+1), dom%node%elts(id+1))
    dom%len%elts(EDGE*id+UP+1)    = dist (dom%node%elts(id+1),   dom%node%elts(idN+1))
    dom%pedlen%elts(EDGE*id+RT+1) = dist (dom%ccentre%elts(TRIAG*idS+UPLT+1), dom%ccentre%elts(LORT+TRIAG*id+1))
    dom%pedlen%elts(EDGE*id+DG+1) = dist (dom%ccentre%elts(TRIAG*id+UPLT+1),  dom%ccentre%elts(LORT+TRIAG*id+1))
    dom%pedlen%elts(EDGE*id+UP+1) = dist (dom%ccentre%elts(TRIAG*id+UPLT+1),  dom%ccentre%elts(LORT+TRIAG*idW+1))

    if (i == -1 .and. j == -1 .and. is_penta (dom, p, IJMINUS - 1)) then
       dom%pedlen%elts(EDGE*id+DG+1) = 0
       dom%len%elts(EDGE*id+1) = dist (dom%node%elts(idE+1), dom%node%elts(idN+1))
    end if

    if (i == PATCH_SIZE .and. j == PATCH_SIZE .and. is_penta(dom, p, IJPLUS - 1)) &
       dom%len%elts(EDGE*id+DG+1) = dist (dom%node%elts(idE+1), dom%node%elts(idN+1))
  end subroutine lengths

  subroutine init_geometry
    implicit none
    integer :: d, i, v

    do d = 1, size(grid)
       call init (grid(d)%ccentre, grid(d)%node%length * TRIAG)

       do i = 1, grid(d)%node%length*TRIAG
          call init_Coord (grid(d)%ccentre%elts(i), 0d0, 0d0, 0d0)
       end do

       call init (grid(d)%midpt, grid(d)%node%length * EDGE)

       do i = 1, grid(d)%node%length*EDGE
          call init_Coord (grid(d)%midpt%elts(i), 0d0, 0d0, 0d0)
       end do

       call init (grid(d)%areas,    grid(d)%node%length)
       call init (grid(d)%pedlen,   grid(d)%node%length * EDGE)
       call init (grid(d)%len,      grid(d)%node%length * EDGE)
       call init (grid(d)%triarea,  grid(d)%node%length * TRIAG)
       call init (grid(d)%coriolis, grid(d)%node%length * TRIAG)
    end do
  end subroutine init_geometry

  subroutine precompute_geometry
    implicit none
    integer :: d, k, v

    call apply_onescale2 (ccentre, min_level-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)

    do d = 1, size(grid)
       call ccentre_penta (grid(d), 1)
    end do

    call apply_onescale2 (midpt,      min_level-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)
    call apply_onescale2 (cpt_areas,  min_level-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)
    call apply_onescale2 (lengths,    min_level-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)

    call apply_onescale (cpt_triarea, min_level-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)
    call apply_onescale (coriolis,    min_level-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)
    
    do d = 1, size(grid)
       call init (grid(d)%surf_press,     grid(d)%node%length)
       call init (grid(d)%press,          grid(d)%node%length)
       call init (grid(d)%geopot,         grid(d)%node%length)
       call init (grid(d)%u_zonal,        grid(d)%node%length)
       call init (grid(d)%v_merid,        grid(d)%node%length)
       call init (grid(d)%press_lower,    grid(d)%node%length)
       call init (grid(d)%geopot_lower,   grid(d)%node%length)
       call init (grid(d)%bernoulli,      grid(d)%node%length)
       call init (grid(d)%ke,             grid(d)%node%length)
       call init (grid(d)%divu,           grid(d)%node%length)
       call init (grid(d)%vort,   TRIAG * grid(d)%node%length)
       call init (grid(d)%qe,      EDGE * grid(d)%node%length)
       
       call init (Laplacian_vector(S_DIVU)%data(d),        grid(d)%node%length)
       call init (Laplacian_vector(S_ROTU)%data(d), EDGE * grid(d)%node%length)
       
       do v = scalars(1), scalars(2)
          call init (horiz_flux(v)%data(d), EDGE * grid(d)%node%length)
          call init (Laplacian_scalar(v)%data(d),  grid(d)%node%length)
       end do
       
       call init (topography%data(d), grid(d)%node%length)
       if (sso) then
          do k = 1, 4
             call init (sso_param(k)%data(d), grid(d)%node%length)
          end do
       end if
       
       do k = zmin, zmax
          call init (penal_node(k)%data(d),        grid(d)%node%length)
          call init (penal_edge(k)%data(d), EDGE * grid(d)%node%length)
          call init (exner_fun(k)%data(d),         grid(d)%node%length)
          if (k > 0) then
             do v = scalars(1), scalars(2)
                call init (trend(v,k)%data(d), grid(d)%node%length)
             end do
             call init (trend(S_VELO,k)%data(d), EDGE * grid(d)%node%length)
          end if
       end do
       call init (exner_fun(zmax+1)%data(d),  grid(d)%node%length)
    end do

    if (vert_diffuse) then
       do d = 1, size(grid)
          call init (Kt(0)%data(d), grid(d)%node%length)
          call init (Kv(0)%data(d), grid(d)%node%length)
          do k = 1, zlevels
             call init (Kt(k)%data(d),  grid(d)%node%length)
             call init (Kv(k)%data(d),  grid(d)%node%length)
             call init (tke(k)%data(d), grid(d)%node%length)
          end do
       end do
    end if
  end subroutine precompute_geometry

  subroutine init_connections
    implicit none
    logical, dimension(2) :: pole_assigned
    integer, dimension(2) :: neigh_over_pole
    integer               :: d, d_glo, i, ii, j, jj, k, loz, ngb_loz, rot, s, split

    pole_assigned = .false.

    do ii = 1, 2
       do jj = 1, 5
          loz = 5*ii - 5 + (jj-1)
          do i = 1, N_SUB_DOM_PER_DIM
             do j = 1, N_SUB_DOM_PER_DIM
                split = N_SUB_DOM_PER_DIM * (j-1) + (i-1)
                d_glo = N_SUB_DOM*loz + split

                if (.not. owner(d_glo+1) == rank) then
                   ! check if pole master on other rank
                   if (ii-1 == 1 .and. i-1 == 0 .and. j-1 == N_SUB_DOM_PER_DIM - 1) pole_assigned(1) = .true.
                   if (ii-1 == 0 .and. j-1 == 0 .and. i-1 == N_SUB_DOM_PER_DIM - 1) pole_assigned(2) = .true.
                   cycle
                end if

                d = loc_id(d_glo+1)
                grid(d+1)%neigh = NONE

                neigh_over_pole = (/5*ii - 5 + modulo(jj+1, 5), 5*ii - 5 + modulo(jj - 3, 5)/)

                if (ii-1 == 1 .and. i-1 == 0 .and. j-1 == N_SUB_DOM_PER_DIM - 1) then
                   grid(d+1)%neigh(NORTHWEST) = POLE
                   grid(d+1)%neigh_over_pole = N_SUB_DOM*neigh_over_pole + split

                   if (.not. pole_assigned(1)) then
                      grid(d+1)%pole_master(2) = .True.
                      pole_assigned(1) = .True.
                   end if
                else
                   if (ii-1 == 0 .and. j-1 == 0 .and. i-1 == N_SUB_DOM_PER_DIM - 1) then
                      grid(d+1)%neigh(SOUTHEAST) = POLE
                      grid(d+1)%neigh_over_pole = N_SUB_DOM*neigh_over_pole + split
                      if (.not. pole_assigned(2)) then
                         grid(d+1)%pole_master(1) = .True.
                         pole_assigned(2) = .True.
                      end if
                   end if
                end if

                call init (grid(d+1)%neigh_pa_over_pole, min_level*2)

                grid(d+1)%neigh_pa_over_pole%elts(min_level*2-1:min_level*2) = 1
                grid(d+1)%neigh_rot = 0

                grid(d+1)%neigh(NORTH) = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*j + (i-1)
                grid(d+1)%neigh(EAST)  = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j-1) + i-1 + 1
                grid(d+1)%neigh(SOUTH) = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j - 2) + (i-1)
                grid(d+1)%neigh(WEST) = (N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j-1) + i-1) - 1

                if (i < N_SUB_DOM_PER_DIM) then
                   if (j < N_SUB_DOM_PER_DIM) grid(d+1)%neigh(IJPLUS)      = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*j     + i-1 + 1
                   if (j-1 > 0)               grid(d+1)%neigh(IPLUSJMINUS) = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j-2) + i-1 + 1
                end if

                if (i-1 > 0) then
                   if (j < N_SUB_DOM_PER_DIM) grid(d+1)%neigh(IMINUSJPLUS) = (N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*j+i-1)     - 1
                   if (j-1 > 0)               grid(d+1)%neigh(IJMINUS)     = (N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j-2)+i-1) - 1
                end if

                if (j == N_SUB_DOM_PER_DIM) then
                   rot = ii-1
                   ngb_loz = (5 + modulo(ii + jj - 2, 5))*N_SUB_DOM
                   call set_dom_neigh (d, NORTH, ngb_loz, i-1, 0, NORTH, rot)

                   if (i < N_SUB_DOM_PER_DIM) then
                      call set_dom_neigh (d, IJPLUS, ngb_loz, i, 0, NORTH, rot)
                   end if

                   if (i-1 > 0) then
                      call set_dom_neigh (d, IMINUSJPLUS, ngb_loz, i - 2, 0, NORTH, rot)
                   end if
                end if

                if (i == N_SUB_DOM_PER_DIM) then
                   rot = 1 - (ii-1)
                   ngb_loz = (0 + modulo(jj, 5))*N_SUB_DOM
                   call set_dom_neigh (d, EAST, ngb_loz, 0, j-1, EAST, rot)

                   if (j < N_SUB_DOM_PER_DIM) call set_dom_neigh (d, IJPLUS,      ngb_loz, 0, j,   EAST, rot)
                   if (j-1 > 0)               call set_dom_neigh (d, IPLUSJMINUS, ngb_loz, 0, j-2, EAST, rot)
                end if

                if (j-1 == 0) then
                   rot = 1 - (ii-1)
                   ngb_loz = (0 + modulo(ii + jj - 3, 5))*N_SUB_DOM

                   call set_dom_neigh (d, SOUTH, ngb_loz, i-1, N_SUB_DOM_PER_DIM - 1, SOUTH, rot)

                   if (i < N_SUB_DOM_PER_DIM) call set_dom_neigh (d, IPLUSJMINUS, ngb_loz, i,   N_SUB_DOM_PER_DIM-1, SOUTH, rot)
                   if (i-1 > 0)               call set_dom_neigh (d, IJMINUS,     ngb_loz, i-2, N_SUB_DOM_PER_DIM-1, SOUTH, rot)
                end if

                if (i-1 == 0) then
                   rot = ii-1
                   ngb_loz = (5 + modulo(jj - 2, 5))*N_SUB_DOM

                   call set_dom_neigh(d, WEST, ngb_loz, N_SUB_DOM_PER_DIM - 1, j-1, WEST, rot)

                   if (j < N_SUB_DOM_PER_DIM) call set_dom_neigh (d, IMINUSJPLUS, ngb_loz, N_SUB_DOM_PER_DIM-1, j,   WEST, rot)
                   if (j-1 > 0)             call set_dom_neigh (d, IJMINUS,     ngb_loz, N_SUB_DOM_PER_DIM-1, j-2, WEST, rot)
                end if
                do s = 1, N_BDRY
                   grid(d+1)%bdry_patch%elts(s+1)%neigh = 1
                end do
             end do
          end do
       end do
    end do
  end subroutine init_connections

  subroutine set_penta (d_glo, side)
    implicit none
    integer :: d_glo, side

    if (owner(d_glo+1) == rank) grid(loc_id(d_glo+1)+1)%penta(side) = .true.
  end subroutine set_penta

  subroutine assign_coord (dom, p, ne, se, sw, nw)
    implicit none
    type(Domain) :: dom
    integer      :: p
    type(Coord)  :: ne, se, sw, nw

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: i, j, l, pc_incr, pf_offs

    call get_offs_Domain (dom, p, offs, dims)

    dom%node%elts(offs(NORTHEAST+1) + 1) = project_on_sphere (ne)
    dom%node%elts(offs(EAST+1) + 1)      = project_on_sphere (se)
    dom%node%elts(offs(1) + 1)           = project_on_sphere (sw)
    dom%node%elts(offs(NORTH+1) + 1)     = project_on_sphere (nw)

    do l = 1, PATCH_LEVEL
       pc_incr = 2**(PATCH_LEVEL - l + 1)
       pf_offs = pc_incr/2

       do j = 1, PATCH_SIZE, pc_incr
          do i = 1, PATCH_SIZE, pc_incr
             dom%node%elts(idx(i+pf_offs-1, j-1, offs, dims) + 1) = mid_pt(dom%node%elts(idx(i-1, j-1, offs, dims) + 1), &
                  dom%node%elts(idx(i + pc_incr - 1, j-1, offs, dims) + 1))

             dom%node%elts(idx(i-1, j+pf_offs-1, offs, dims) + 1) = mid_pt(dom%node%elts(idx(i-1, j-1, offs, dims) + 1), &
                  dom%node%elts(idx(i-1, j+pc_incr-1, offs, dims) + 1))

             dom%node%elts(idx(i+pf_offs-1, j+pf_offs-1, offs, dims) + 1) = &
                  mid_pt(dom%node%elts(idx(i+pc_incr-1, j+pc_incr-1, offs, dims) + 1), &
                  dom%node%elts(idx(i-1, j-1, offs, dims) + 1))
          end do

          i = PATCH_SIZE
          dom%node%elts(idx(i, j+pf_offs-1, offs, dims) + 1) = mid_pt(dom%node%elts(idx(i, j-1, offs, dims) + 1), &
               dom%node%elts(idx(i, j+pc_incr-1, offs, dims) + 1))
       end do

       j = PATCH_SIZE
       do i = 1, PATCH_SIZE, pc_incr
          dom%node%elts(idx(i+pf_offs - 1, j, offs, dims) + 1) = mid_pt(dom%node%elts(idx(i-1, j, offs, dims) + 1), &
               dom%node%elts(idx(i+pc_incr-1, j, offs, dims) + 1))
       end do
    end do
  end subroutine assign_coord

  subroutine cpt_areas (dom, p, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer     :: idNW, idN, idNE, idW, id, idE, idSW, idS, idSE
    type(Areas) :: area

    idNW = idx(i-1, j+1, offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    id   = idx(i,   j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idSE = idx(i+1, j-1, offs, dims)

    call init_Areas (area, dom%node%elts(id+1), &
         
         (/ dom%ccentre%elts(TRIAG*id  +LORT+1), &
            dom%ccentre%elts(TRIAG*id  +UPLT+1), &
            dom%ccentre%elts(TRIAG*idW +LORT+1), &
            dom%ccentre%elts(TRIAG*idSW+UPLT+1), &
            dom%ccentre%elts(TRIAG*idSW+LORT+1), &
            dom%ccentre%elts(TRIAG*idS +UPLT+1) /), &
         
         (/ dom%midpt%elts(EDGE*id  +RT+1), &
            dom%midpt%elts(EDGE*id  +DG+1), &
            dom%midpt%elts(EDGE*id  +UP+1), &
            dom%midpt%elts(EDGE*idW +RT+1), &
            dom%midpt%elts(EDGE*idSW+DG+1), &
            dom%midpt%elts(EDGE*idS +UP+1)/) )

    if (j >= PATCH_SIZE + 1) then
       dom%areas%elts(id+1)%part(4:6) = area%part(4:6)
       return
    end if

    if (i >= PATCH_SIZE + 1) then
       dom%areas%elts(id+1)%part(3:5) = area%part(3:5)
       return
    end if

    dom%areas%elts(id+1) = area
  end subroutine cpt_areas

  integer function sub_dom_id (i, j, s, rot)
    implicit none
    integer :: i, j, s, rot
    
    integer, dimension(2) :: ij

    ij = (/i, j/)
    if (rot == 1) then
       ij(modulo(s, 2) + 1) = N_SUB_DOM_PER_DIM - 1 - ij(modulo(s, 2) + 1)
       ij = (/ij(2), ij(1)/)
    end if

    sub_dom_id = ij(2) * N_SUB_DOM_PER_DIM + ij(1)
  end function sub_dom_id

  subroutine ccentre (dom, p, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id, idN, idE, idNE, idS, idW

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)

    dom%ccentre%elts(TRIAG*id+LORT+1) = circumcentre (dom%node%elts(id+1), dom%node%elts(idNE+1), dom%node%elts(idE+1))
    dom%ccentre%elts(TRIAG*id+UPLT+1) = circumcentre (dom%node%elts(id+1), dom%node%elts(idN+1),  dom%node%elts(idNE+1))
  end subroutine ccentre

  subroutine cpt_triarea (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id, idN, idE, idNE

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    dom%triarea%elts(TRIAG*id+LORT+1) = dom%areas%elts(id+1)%part(1) + dom%areas%elts(idE+1)%part(3) &
         + dom%areas%elts(idNE+1)%part(5)
    dom%triarea%elts(TRIAG*id+UPLT+1) = dom%areas%elts(id+1)%part(2) + dom%areas%elts(idNE+1)%part(4) &
         + dom%areas%elts(idN+1)%part(6)
  end subroutine cpt_triarea

  subroutine coriolis (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idN, idE, idNE

    id   = idx (i,   j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    dom%coriolis%elts(TRIAG*id+LORT+1) = dom%ccentre%elts(TRIAG*id+LORT+1)%z/radius * 2d0*omega * &
         (dom%areas%elts(id+1)%part(1) + dom%areas%elts(idE+1)%part(3) + dom%areas%elts(idNE+1)%part(5))

    dom%coriolis%elts(TRIAG*id+UPLT+1) = dom%ccentre%elts(TRIAG*id+UPLT+1)%z/radius * 2d0*omega * &
         (dom%areas%elts(id+1)%part(2) + dom%areas%elts(idNE+1)%part(4) + dom%areas%elts(idN+1)%part(6))
  end subroutine coriolis

  subroutine set_dom_neigh (d, s, ngb_loz, i, j, s1, rot)
    implicit none
    integer :: d, s, ngb_loz, i, j, s1, rot

    call set_neigh_Domain (grid(d+1), s, ngb_loz + sub_dom_id(i, j, s1 - 1, rot), rot)
  end subroutine set_dom_neigh

  subroutine set_level (dom, p, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    dom%level%elts(idx(i,j,offs,dims)+1) = dom%patch%elts(p+1)%level
  end subroutine set_level

  subroutine midpt (dom, p, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id, idN, idE, idNE

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    if (j >= PATCH_SIZE + 1) then
       dom%midpt%elts(EDGE*id+RT+1) = mid_pt(dom%node%elts(id+1), dom%node%elts(idE+1))
       return
    end if

    if (i >= PATCH_SIZE + 1) then
       dom%midpt%elts(EDGE*id+UP+1) = mid_pt(dom%node%elts(id+1), dom%node%elts(idN+1))
       return
    end if

    if (j == -1) then
       if (i == -1 .and. is_penta(dom, p, IJMINUS - 1)) then
          dom%midpt%elts(EDGE*id+DG+1) = dom%ccentre%elts(TRIAG*id+LORT+1)
          return
       else
          if (i == PATCH_SIZE .and. is_penta(dom, p, IPLUSJMINUS - 1)) then
             dom%midpt%elts(EDGE*id+UP+1) = dom%ccentre%elts(TRIAG*id+UPLT+1)
             return
          end if
       end if
    else
       if (j == PATCH_SIZE .and. is_penta(dom, p, IMINUSJPLUS - 1)) then
          if (i == -1) then
             dom%midpt%elts(EDGE*id+RT+1) = dom%ccentre%elts(TRIAG*id+LORT+1)
             return
          end if
       end if
    end if

    dom%midpt%elts(EDGE*id+RT+1) = mid_pt(dom%node%elts(id+1),   dom%node%elts(idE+1))
    dom%midpt%elts(EDGE*id+DG+1) = mid_pt(dom%node%elts(idNE+1), dom%node%elts(id+1))
    dom%midpt%elts(EDGE*id+UP+1) = mid_pt(dom%node%elts(id+1),   dom%node%elts(idN+1))

    if (j == PATCH_SIZE) then
       if (i == PATCH_SIZE .and. is_penta(dom, p, IJPLUS - 1)) dom%midpt%elts(EDGE*id+DG+1) = dom%ccentre%elts(TRIAG*id+LORT+1)
    end if
  end subroutine midpt
end module init_mod
