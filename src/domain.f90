 module domain_mod
  use param_mod
  use shared_mod
  use geom_mod
  use patch_mod
  use dyn_arrays
  use arch_mod
  implicit none

  integer, dimension(2,N_BDRY+1) :: sides_dims
  integer, dimension(2,4)        :: chd_offs

  ! Objects same on all zlevels
  type Domain
     ! Geometry
     integer, dimension(2)      :: neigh_over_pole
     integer, dimension(N_BDRY) :: neigh, neigh_rot
     logical, dimension(2)      :: pole_master
     logical, dimension(N_BDRY) :: penta
     type(Int_Array)            :: neigh_pa_over_pole
     type(Coord_Array)          :: ccentre     ! circumcentres
     type(Coord_Array)          :: midpt       ! midpoints of edges
     type(Coord_Array)          :: node        ! coordinates of hexagon nodes
     type(Areas_Array)          :: areas       ! hexagon areas
     type(Float_Array)          :: len         ! primal edge lengths
     type(Float_Array)          :: pedlen      ! dual edge lengths
     type(Float_Array)          :: triarea     ! triangle areas
     type(Overl_Area_Array)     :: overl_areas ! overlapping areas

     ! Multiscale and data structure
     integer                                    :: id
     type(Int_Array)                            :: level
     type(Int_Array)                            :: mask_n
     type(Int_Array)                            :: mask_e
     type(Int_Array), dimension(:), allocatable :: lev
     type(Iu_Wgt_Array)                         :: I_u_wgt
     type(RF_Wgt_Array)                         :: R_F_wgt
     type(Patch_Array)                          :: patch
     type(Bdry_Patch_Array)                     :: bdry_patch

     ! Communication
     type(Int_Array)                                          :: send_pa_all
     type(Int_Array), dimension(N_GLO_DOMAIN)                 :: recv_pa, send_conn
     type(Int_Array), dimension(AT_NODE:AT_EDGE,N_GLO_DOMAIN) :: pack, unpk
     type(Int_Array), dimension(:,:), allocatable             :: src_patch
     
     ! Physical quantities (should evaluate on -(BDRY_THICKNESS-1), BDRY_THICKNESS since cannot use bdry_update)
     type(Float_Array) :: coriolis    ! Coriolis force
     type(Float_Array) :: surf_press  ! surface pressure (compressible) or surface Lagrange multiplier (incompressible)
     type(Float_Array) :: press       ! pressure (compressible case) or Lagrange multiplier (incompressible case)
     type(Float_Array) :: geopot      ! geopotential
     type(Float_Array) :: u_zonal     ! zonal velocity
     type(Float_Array) :: v_merid     ! meridional velocity
     type(Float_Array) :: press_lower ! mass in adjacent vertical cell
     type(Float_Array) :: geopot_lower! geopotential in adjacent vertical cell
     type(Float_Array) :: vort        ! vorticity on triangles
     type(Float_Array) :: ke          ! kinetic energy
     type(Float_Array) :: bernoulli   ! Bernoulli function
     type(Float_Array) :: divu        ! divergence of velocity
     type(Float_Array) :: qe          !
  end type Domain

  type Float_Field
     integer                                      :: pos
     logical                                      :: bdry_uptodate
     type(Float_Array), dimension(:), allocatable :: data
  end type Float_Field

  type(Domain), dimension(:), allocatable, target        :: grid

  type(Float_Field),                              target :: topography
  type(Topo_Array),  dimension(:,:), allocatable         :: topography_data
  type(Float_Field), dimension(4),                target :: sso_param
  type(Float_Field), dimension(:),   allocatable, target :: exner_fun, horiz_flux, penal_node, penal_edge
  type(Float_Field), dimension(:),   allocatable, target :: Kv, Kt, tke, wav_tke
  type(Float_Field), dimension(:),   allocatable, target :: Laplacian_scalar, Laplacian_vector
  type(Float_Field), dimension(:,:), allocatable, target :: sol, sol_mean, trend
  type(Float_Field), dimension(:,:), allocatable, target :: wav_coeff

  real(8), dimension(:), pointer :: diag, mass, mass1, h_flux, h_mflux
  real(8), dimension(:), pointer :: dmass, dtemp, dscalar, scalar, scalar_2d, temp, temp1
  real(8), dimension(:), pointer :: velo, velo1, velo2, velo_2d, dvelo, dvelo_2d
  real(8), dimension(:), pointer :: mean_m, mean_t
  real(8), dimension(:), pointer :: Laplacian
  real(8), dimension(:), pointer :: bernoulli, divu, exner, ke, qe, vort
  real(8), dimension(:), pointer :: wc_s, wc_u
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Initialization subroutines !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_Float_Field (self, pos)
    implicit none
    type(Float_Field) :: self
    integer           :: pos

    self%bdry_uptodate = .false.
    self%pos = pos

    allocate (self%data(n_domain(rank+1)))
  end subroutine init_Float_Field

  subroutine init_domain_mod
    implicit none
    integer :: i, v
    logical :: initialized = .false.

    if (initialized) return ! initialize only once

    call init_shared_mod
    call init_patch_mod
    call init_arch_mod

    sides_dims = reshape ((/PATCH_SIZE, PATCH_SIZE, PATCH_SIZE, &
         BDRY_THICKNESS, BDRY_THICKNESS, PATCH_SIZE, PATCH_SIZE, &
         BDRY_THICKNESS, BDRY_THICKNESS, PATCH_SIZE, BDRY_THICKNESS, &
         BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS, &
         BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS/), (/2, 9/))

    chd_offs = reshape ((/PATCH_SIZE/2, PATCH_SIZE/2, PATCH_SIZE/2, 0, 0, 0, 0, PATCH_SIZE/2/), (/2, 4/))

    initialized = .true.
  end subroutine init_domain_mod

  subroutine init_Domain (self)
    implicit none
    type(Domain) :: self
    
    integer :: i, k, l, r

    call init (self%patch, 1)
    call init (self%bdry_patch, 1)
    call init (self%node, 1)

    allocate (self%src_patch(n_process,min_level:max_level))

    do l = min_level, max_level
       do r = 1, n_process
          call init (self%src_patch(r,l), 0)
       end do
    end do

    allocate (self%lev(min_level-1:max_level))

    do i = lbound(self%lev,1), ubound(self%lev,1)
       call init (self%lev(i), 0)
    end do

    do i = 1, N_GLO_DOMAIN
       call init (self%send_conn(i), 0)
    end do

    call init (self%send_pa_all, 0)

    do i = 1, N_GLO_DOMAIN
       call init (self%recv_pa(i), 0)
    end do

    do k = AT_NODE, AT_EDGE
       do i = 1, N_GLO_DOMAIN
          call init (self%pack(k,i), 0)
          call init (self%unpk(k,i), 0)
       end do
    end do

    self%pole_master = .false.
  end subroutine init_Domain

  integer function add_patch_Domain (self, level)
    ! Add new patch to the domain
    implicit none
    type(Domain) :: self
    
    integer :: level, p

    p = self%patch%length

    call append (self%lev(level), p)
    call append (self%patch, Patch (self%node%length, level, 0, 0, 0, .false.))
    
    call extend_Domain (self, PATCH_SIZE**2)
    
    add_patch_Domain = p
  end function add_patch_Domain

  integer function add_bdry_patch_Domain (self, side)
    ! Add boundary patch to the domain
    implicit none
    type(Domain) :: self
    integer      :: side

    integer :: p

    p = self%bdry_patch%length

    call append (self%bdry_patch, Bdry_Patch (self%node%length, side, 0))

    call extend_Domain (self, BDRY_THICKNESS * PATCH_SIZE)

    add_bdry_patch_Domain = p
  end function add_bdry_patch_Domain

  subroutine extend_Domain (self, num)
    implicit none
    type(Domain) :: self
    integer      :: num
    
    integer :: d, k, v
    real(8) :: def_val, dz, z

    d = self%id + 1

    call extend (self%node, num, ORIGIN)

    ! Atmosphere/ocean layers
    do k = 1, zmax
       ! Set reasonable default values for new boundary patches to avoid NaN if variable undefined in boundary
       if (compressible) then
          def_val = a_vert_mass(k) + b_vert_mass(k) * p_0 / grav_accel
       else
          dz     = b_vert_mass(k) * max_depth
          z      = 0.5d0 * (b_vert(k) + b_vert(k-1)) * max_depth
          def_val = ref_density * dz
       end if

       do v = scalars(1), scalars(2)
          if (split_mean_perturbation) then 
             call extend (sol(v,k)%data(d),      num, 0d0)     
             call extend (sol_mean(v,k)%data(d), num, def_val) 
          else
             call extend (sol(v,k)%data(d),      num, def_val) 
             call extend (sol_mean(v,k)%data(d), num, 0d0)     
          end if
       end do
       call extend (sol(S_VELO,k)%data(d),      EDGE * num, 0d0)
       call extend (sol_mean(S_VELO,k)%data(d), EDGE * num, 0d0)
    end do

    ! Soil layers
    do k = zmin, 0
       do v = scalars(1), scalars(2)
          if (split_mean_perturbation) then 
             call extend (sol(v,k)%data(d),      num, 0d0)     
             call extend (sol_mean(v,k)%data(d), num, 0d0) 
          else
             call extend (sol(v,k)%data(d),      num, 0d0) 
             call extend (sol_mean(v,k)%data(d), num, 0d0)     
          end if
       end do
       call extend (sol(S_VELO,k)%data(d),      EDGE * num, 0d0)
       call extend (sol_mean(S_VELO,k)%data(d), EDGE * num, 0d0)
    end do
  end subroutine extend_Domain

  subroutine set_neigh_Domain (self, s, id, rot)
    implicit none
    type(Domain) :: self
    integer      :: id, rot, s

    self%neigh(s) = id
    self%neigh_rot(s) = rot
  end subroutine set_neigh_Domain

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Subroutines for applying a routine to domains, patches, ... !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine apply (routine, zlev)
    ! Applies routine over all levels and over entire boundary
    implicit none
    integer  :: zlev
    external :: routine
    
    integer :: l

    do l = level_start, level_end
       call apply_onescale (routine, l, zlev, -(BDRY_THICKNESS-1), BDRY_THICKNESS)
    end do
  end subroutine apply

  subroutine apply_bdry (routine, zlev, st, en)
    ! Applies routine to nodes/edges at all levels and including boundary cells specified by (st,en)
    implicit none
    integer  :: en, st, zlev
    external :: routine
    
    integer :: l

    do l = level_start, level_end
       call apply_onescale (routine, l, zlev, st, en)
    end do
  end subroutine apply_bdry

  subroutine apply_onescale (routine, l, zlev, st, en)
     ! Applies routine to nodes/edges at all level l including boundary cells specified by (st,en)
    implicit none
    external :: routine
    integer  :: en, l, st, zlev

    integer :: d, j

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (routine, grid(d), grid(d)%lev(l)%elts(j), zlev, st, en)
       end do
    end do
  end subroutine apply_onescale

  subroutine apply_onescale__int (routine, l, zlev, st, en, ival)
    ! Applies routine to nodes/edges at all level l including boundary cells specified by (st,en)
    ! and passes integer ival to routine
    implicit none
    external :: routine
    integer  :: en, ival, l, st, zlev

    integer               :: d, j
    logical, dimension(2) :: pole_done

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch__int (routine, grid(d), grid(d)%lev(l)%elts(j), zlev, st, en, ival)
       end do
    end do
  end subroutine apply_onescale__int

  subroutine apply_no_bdry (routine, zlev)
    ! Applies routine to nodes/edges at all levels excluding all boundary cells
    !
    ! The integer flag is_pole allows the routine to not update values for edges if is_pole = 1
    !
    ! Routine needs to accept an integer which indicates pole (if equals 1) or non-pole (if equals 0)
    !
    implicit none
    integer  :: zlev
    external :: routine

    integer :: is_pole, l

    do l = level_start, level_end
       is_pole = 0; call apply_onescale__int (routine, l, zlev, 0, 0, is_pole)
       is_pole = 1; call apply_to_pole       (routine, l, zlev,       is_pole, .true.)  
    end do
  end subroutine apply_no_bdry

  subroutine apply_d (routine, dom, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, st, zlev

    integer :: j, l
    
    do l = level_start, level_end
       do j = 1, dom%lev(l)%length
          call apply_onescale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
       end do
    end do
  end subroutine apply_d
  
  subroutine apply_onescale_to_patch (routine, dom, p, zlev, st, en)
    ! Applies routine to nodes/edges at all level associated with p.
    ! Includes boundary cells specified by (st,en).
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, p, st, zlev

    integer                          :: i, j
    integer, dimension(N_BDRY+1)     :: offs
    integer, dimension(2,N_BDRY+1)   :: dims
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry

    call get_offs_Domain (dom, p, offs, dims, inner_bdry)

    bdry = (/en, en, st, st/)

    where (inner_bdry) bdry = 0

    do j = bdry(JMINUS) + 1, PATCH_SIZE + bdry(JPLUS)
       do i = bdry(IMINUS) + 1, PATCH_SIZE + bdry(IPLUS)
          call routine (dom, i - 1, j - 1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch

  subroutine apply_onescale_d (routine, dom, l, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, l, st, zlev

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_onescale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_onescale_d

  subroutine apply_onescale_to_patch__int (routine, dom, p, zlev, st, en, ival)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, ival, p, st, zlev

    integer :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call get_offs_Domain(dom, p, offs, dims)

    do j = st + 1, PATCH_SIZE + en
       do i = st + 1, PATCH_SIZE + en
          call routine (dom, p, i - 1, j - 1, zlev, offs, dims, ival)
       end do
    end do
  end subroutine apply_onescale_to_patch__int

  subroutine apply_onescale2 (routine, l, zlev, st, en)
    implicit none
    external :: routine
    integer  :: en, l, st, zlev

    integer :: d, j

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch2 (routine, grid(d), grid(d)%lev(l)%elts(j), zlev, st, en)
       end do
    end do
  end subroutine apply_onescale2

  subroutine apply_onescale_to_patch5 (routine, dom, p, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, p, st, zlev

    integer                          :: i, j
    integer, dimension(N_BDRY+1)     :: offs
    integer, dimension(2,N_BDRY+1)   :: dims
    logical, dimension(JPlUS:IMINUS) :: inner_patch
    integer, dimension(JPlUS:IMINUS) :: bdry

    call get_offs_Domain5 (dom, p, offs, dims, inner_patch)

    bdry = (/en, en, st, st/)

    where (inner_patch) bdry = 0

    do j = bdry(JMINUS) + 1, PATCH_SIZE + bdry(JPLUS)
       do i = bdry(IMINUS) + 1, PATCH_SIZE + bdry(IPLUS)
          call routine (dom, p, i - 1, j - 1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch5

  subroutine apply_onescale_to_patch2 (routine, dom, p, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, p, st, zlev

    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call get_offs_Domain(dom, p, offs, dims)

    do j = st + 1, PATCH_SIZE + en
       do i = st + 1, PATCH_SIZE + en
          call routine (dom, p, i - 1, j - 1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch2

  subroutine apply_interscale (routine, l, zlev, st, en)
    ! Applies interscale routine to coarse scale l and fine scale l+1
    implicit none
    external :: routine
    integer  :: en, l, st, zlev

    integer :: d

    do d = 1, size(grid)
       call apply_interscale_d (routine, grid(d), l, zlev, st, en)
    end do
  end subroutine apply_interscale

  subroutine apply_interscale_d (routine, dom, l, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, l, st, zlev

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_interscale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_interscale_d
  
  subroutine apply_interscale_d2 (routine, dom, l, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, l, st, zlev

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_interscale_to_patch22 (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_interscale_d2

  subroutine apply_interscale_to_patch (routine, dom, p_par, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, p_par, st, zlev

    integer                          :: c, i, i_chd, i_par, j, j_chd, j_par, p_chd
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do c = 1, N_CHDRN
       p_chd = dom%patch%elts(p_par+1)%children(c)
       if (p_chd == 0) cycle

       call get_offs_Domain (dom, p_chd, offs_chd, dims_chd, inner_bdry)

       bdry = (/en, en, st, st/)

       where (inner_bdry) bdry = 0

       do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
          j_chd = (j - 1)*2
          j_par = j - 1 + chd_offs(2,c)
          do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
             i_chd = (i - 1)*2
             i_par = i - 1 + chd_offs(1,c)
             call routine (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
          end do
       end do
    end do
  end subroutine apply_interscale_to_patch

  subroutine apply_interscale_to_patch2 (routine, dom, p_par, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, st, p_par, zlev

    integer                          :: c, i, i_chd, i_par, j, j_chd, j_par, p_chd
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry
    !TODO{uncomment & test}  if (dom%patch%elts(p_par+1)%active == NONE) return

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do c = 1, N_CHDRN
       p_chd = dom%patch%elts(p_par+1)%children(c)

       if (p_chd == 0) cycle

       call get_offs_Domain (dom, p_chd, offs_chd, dims_chd, inner_bdry)

       bdry = (/en, en, st, st/)

       where (inner_bdry) bdry = 0

       do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
          j_chd = (j - 1)*2
          j_par = j - 1 + chd_offs(2,c)
          do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
             i_chd = (i - 1)*2
             i_par = i - 1 + chd_offs(1,c)
             call routine (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
          end do
       end do
    end do
  end subroutine apply_interscale_to_patch2

  subroutine apply_interscale_to_patch22 (routine, dom, p_par, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: en, p_par, st, zlev
    
    integer                          :: c, i, i_chd, i_par, j, j_chd, j_par, p_chd
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do c = 1, N_CHDRN
       p_chd = dom%patch%elts(p_par+1)%children(c)
       if (p_chd == 0) cycle

       call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

       bdry = (/en, en, st, st/)

       do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
          j_chd = (j - 1)*2
          j_par = j - 1 + chd_offs(2,c)
          do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
             i_chd = (i - 1)*2
             i_par = i - 1 + chd_offs(1,c)
             call routine (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
          end do
       end do
    end do
  end subroutine apply_interscale_to_patch22
  
  subroutine apply_interscale_to_patch3 (routine, dom, p_par, c, zlev, st, en)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: c, en, p_par, st, zlev
    
    integer                          :: i, i_chd, i_par, j, j_chd, j_par, p_chd
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry
    !TODO{uncomment & test}  if (dom%patch%elts(p_par+1)%active == NONE) return

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    p_chd = dom%patch%elts(p_par+1)%children(c)

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd, inner_bdry)

    bdry = (/en, en, st, st/)

    do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
       j_chd = (j - 1)*2
       j_par = j - 1 + chd_offs(2,c)
       do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
          i_chd = (i - 1)*2
          i_par = i - 1 + chd_offs(1,c)
          call routine (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
       end do
    end do
  end subroutine apply_interscale_to_patch3

  subroutine apply_to_pole (routine, l, zlev, ival, to_all)
    implicit none
    integer  :: ival, l, zlev
    logical  :: to_all
    external :: routine
    
    integer :: d

    do d = 1, size(grid)
       call apply_to_pole_d (routine, grid(d), l, zlev, ival, to_all)
    end do
  end subroutine apply_to_pole

  subroutine apply_to_pole_d (routine, dom, l, zlev, ival, to_all)
    implicit none
    type(Domain) :: dom
    integer      :: ival, l, zlev
    logical      :: to_all
    external     :: routine

    integer                        :: c, l_cur, p, p_par
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do c = SOUTHEAST, NORTHWEST, 2
       if (.not. dom%pole_master(c/2-2) .and. .not. to_all) cycle
       if (.not. dom%penta(c)) cycle
       if (.not. dom%neigh(c) == POLE) cycle

       p = 1
       do while (p > 0)
          p_par = p
          p     = dom%patch%elts(p_par+1)%children(c-4)

          if (.not. l == NONE) then
             l_cur = dom%patch%elts(p_par+1)%level
             if (l_cur < l) then
                cycle
             else
                if (l_cur > l) exit
             end if
          end if

          call get_offs_Domain (dom, p_par, offs, dims)

          if (c == NORTHWEST) then     ! north pole
             call routine (dom, p_par, 0, PATCH_SIZE, zlev, offs, dims, ival) 
          elseif (c == SOUTHEAST) then ! south pole
             call routine (dom, p_par, PATCH_SIZE, 0, zlev, offs, dims, ival) 
          end if
          
       end do
    end do
  end subroutine apply_to_pole_d

  subroutine apply_to_penta (routine, l, zlev)
    implicit none
    external :: routine
    integer  :: l, zlev

    integer                        :: c, d, l_cur, p, p_par
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do d = 1, size(grid)
       call apply_to_penta_d(routine, grid(d), l, zlev)
    end do
  end subroutine apply_to_penta

  subroutine apply_to_penta_d (routine, dom, l, zlev)
    implicit none
    external     :: routine
    type(Domain) :: dom
    integer      :: l, zlev

    integer                        :: c, d, l_cur, p, p_par
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do c = NORTHEAST, NORTHWEST
       if (.not. dom%penta(c)) cycle

       p = 1
       do while (p > 0)
          p_par = p
          p = dom%patch%elts(p_par+1)%children(c-4)
          if (.not. l == NONE) then
             l_cur = dom%patch%elts(p_par+1)%level
             if (l_cur < l) then
                cycle
             else
                if (l_cur > l) exit
             end if
          end if

          call get_offs_Domain (dom, p_par, offs, dims)
          call routine (dom, p_par, c, offs, dims, zlev)
       end do
    end do
  end subroutine apply_to_penta_d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Index and neighbour routines !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function idx (i0, j0, offs, dims)
    ! Given regular array coordinates (i0,j0), offset array offs and domain array dims returns associated grid element as elts(idx+1)
    implicit none
    integer                        :: i0, j0
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: i, j

    i = i0
    j = j0
    if (i < 0) then
       if (j < 0) then
          idx = offs(IJMINUS+1) + (j0 + dims(2,IJMINUS+1))*dims(1,IJMINUS+1) + i + dims(1,IJMINUS+1)
          return
       else
          if (j >= PATCH_SIZE) then
             idx = offs(IMINUSJPLUS+1) + (-PATCH_SIZE + j)*dims(1,IMINUSJPLUS+1) + i + dims(1,IMINUSJPLUS+1)
             return
          else
             idx = offs(WEST+1) + j*dims(1,WEST+1) + i + dims(1,IMINUS+1)
             return
          end if
       end if
    else
       if (i >= PATCH_SIZE) then
          i = i - PATCH_SIZE
          if (j >= PATCH_SIZE) then
             j = j - PATCH_SIZE
             idx = offs(NORTHEAST+1) + j*dims(1,NORTHEAST+1) + i
             return
          else
             if (j < 0) then
                idx = offs(SOUTHEAST+1) + (j0 + dims(2,SOUTHEAST+1))*dims(1,SOUTHEAST+1) + (i0 - PATCH_SIZE) 
                return
             else
                idx = offs(EAST+1) + j*dims(1,EAST+1) + i
                return
             end if
          end if
       else
          if (j < 0) then
             idx = offs(JMINUS+1) + (j + dims(2,JMINUS+1))*dims(1,JMINUS+1) + i
             return
          else
             if (j >= PATCH_SIZE) then
                idx = offs(NORTH+1) + (-PATCH_SIZE + j)*dims(1,NORTH+1) + i
                return
             else
                idx = offs(1) + PATCH_SIZE*j + i
                return
             end if
          end if
       end if
    end if
  end function idx

  integer function idx__fast (i, j, offs)
    implicit none
    integer :: i, j, offs

    idx__fast = PATCH_SIZE*j + i + offs
  end function idx__fast

  function id_edge (id) 
    ! Returns vector with the indices of the three edges associated to node id
    implicit none
    integer                    :: id
    integer, dimension(1:EDGE) :: id_edge
    
    id_edge = EDGE*id + (/ RT, DG, UP /) + 1
  end function id_edge

  function idx_hex (dom, i, j, offs, dims)
    ! Returns vector with the indices of the hexagon and its six neighbours hexagons at node (i,j,offs,dims)
    implicit none
    type(Domain)                     :: dom
    integer                          :: i, j
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims 
    integer, dimension(1:2*EDGE+1)   :: idx_hex

    integer :: id, idE, idNE, idN, idW, idSW, idS

    id   = idx (i,   j,   offs, dims)
    
    idE  = idx (i+1, j,   offs, dims) 
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    idS  = idx (i-1, j,   offs, dims) 
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    idx_hex = (/ id, idE, idNE, idN, idS, idSW, idS /)
  end function idx_hex

  function idx_hex_LORT (dom, i, j, offs, dims)
    ! Returns vector with the indices of the six child hexagons overlapping with the LORT parent triangle
    implicit none
    type(Domain)                     :: dom
    integer                          :: i, j
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims 
    integer, dimension(1:2*EDGE)     :: idx_hex_LORT

    integer :: id, idE, idNE, id2E, id2NE, id2EN

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims) 
    idNE = idx (i+1, j+1, offs, dims)
    
    id2E  = idx (i+2, j,   offs, dims) 
    id2NE = idx (i+2, j+2, offs, dims)
    id2EN = idx (i+2, j+1, offs, dims)

    idx_hex_LORT = (/ id, idE, idNE, id2E, id2NE, id2EN /)
  end function idx_hex_LORT

  function idx_hex_LORT2 (dom, i, j, offs, dims)
    ! Returns vector with the indices of the three hexagons overlapping with the LORT triangle
    implicit none
    type(Domain)                     :: dom
    integer                          :: i, j
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims 
    integer, dimension(1:EDGE)       :: idx_hex_LORT2

    integer :: id, idE, idNE

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims) 
    idNE = idx (i+1, j+1, offs, dims)
    
    idx_hex_LORT2 = (/ id, idE, idNE /)
  end function idx_hex_LORT2

  function idx_hex_UPLT (dom, i, j, offs, dims)
    ! Returns vector with the indices of the six child hexagons overlapping with the LORT parent triangle
    implicit none
    type(Domain)                     :: dom
    integer                          :: i, j
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims 
    integer, dimension(2*EDGE)       :: idx_hex_UPLT

    integer :: id, idNE, idN, id2NE, idE2N, id2N

    id    = idx (i,   j,   offs, dims)
    idNE  = idx (i+1, j+1, offs, dims)
    idN   = idx (i,   j+1, offs, dims)
    
    idE2N = idx (i+1, j+2, offs, dims)
    id2NE = idx (i+2, j+2, offs, dims)
    id2N  = idx (i,   j+2, offs, dims)

    idx_hex_UPLT = (/ id, idNE, idN, id2N, id2NE, idE2N/)
  end function idx_hex_UPLT

  function idx_hex_UPLT2 (dom, i, j, offs, dims)
    ! Returns vector with the indices of the three hexagons overlapping with the LORT triangle
    implicit none
    type(Domain)                     :: dom
    integer                          :: i, j
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims 
    integer, dimension(1:EDGE)       :: idx_hex_UPLT2

    integer :: id, idNE, idN

    id   = idx (i,   j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i+1, j,   offs, dims)
    
    idx_hex_UPLT2 = (/ id, idNE, idN /)
  end function idx_hex_UPLT2

  integer function tri_idx (i, j, tri, offs, dims)
    implicit none
    integer                        :: i, j
    integer, dimension(3)          :: tri
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    tri_idx = TRIAG * idx(i + tri(1), j + tri(2), offs, dims) + tri(3)
  end function tri_idx

  integer function nidx (i, j, s, offs, dims)
    implicit none
    integer                        :: i, j, s
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    nidx = offs(s+1) + j*dims(1,s+1) + i
  end function nidx

  integer function idx2 (i, j, noffs, offs, dims)
    implicit none
    integer                        :: i, j
    integer, dimension(2)          :: noffs
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    idx2 = idx(i + noffs(1), j + noffs(2), offs, dims)
  end function idx2

  integer function ed_idx (i, j, ed, offs, dims)
    ! Return edge index
    implicit none
    integer                        :: i, j
    integer, dimension(3)          :: ed
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    ed_idx = EDGE*idx(i + ed(1), j + ed(2), offs, dims) + ed(3)
  end function ed_idx

  logical function is_penta (dom, p, s)
    implicit none
    type(Domain) :: dom
    integer      :: p, s

    integer :: n, side
    logical :: penta

    penta = .false.

    n = dom%patch%elts(p+1)%neigh(s+1)

    if (n < 0) then
       n = -n
       side = dom%bdry_patch%elts(n+1)%side

       if (side > 0) then
          is_penta = dom%penta(side)
          return
       end if
    end if

    is_penta = penta
  end function is_penta

  integer function find_neigh_bdry_patch_Domain (self, p_par, c, s)
    implicit none
    type(Domain) :: self
    integer      :: p_par, c, s
    
    integer :: p_chd, p_par1, c1

    find_neigh_bdry_patch_Domain = 0

    if (p_par <= 0) return

    p_chd = self%patch%elts(p_par+1)%children(c+1)
    if (p_chd > 0) find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh(s+1)

    if (s >= 4) return

    if (find_neigh_bdry_patch_Domain == 0) then
       call find_neigh_patch2_Domain (self, p_par, c, modulo(s+1,4), p_par1, c1)
       if (p_par1 > 0) then
          p_chd = self%patch%elts(p_par1+1)%children(c1+1)
          if (p_chd > 0) find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh((modulo(s-1,4)+4)+1)
       end if
    end if

    if (find_neigh_bdry_patch_Domain == 0) then
       call find_neigh_patch2_Domain (self, p_par, c, modulo(s-1,4), p_par1, c1)
       if (p_par1 > 0) then
          p_chd = self%patch%elts(p_par1+1)%children(c1+1)
          if (p_chd > 0) find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh((s+4)+1)
       end if
    end if
  end function find_neigh_bdry_patch_Domain

  subroutine find_neigh_patch2_Domain (self, p_par0, c0, s0, p_par, c)
    ! For patch given as `c0`-th child of `p_par0` find neighbour with respect to side `s0`
    ! result as `c`-th child of patch `p_par`
    implicit none
    type(Domain)         :: self
    integer, intent(in)  :: p_par0, c0, s0
    integer, intent(out) :: p_par, c
    
    integer :: s_par

    s_par = par_side(c0, s0)
    if (s0 == c0 .or. s0 == modulo(c0 + 1, 4) .or. (s0 == c0 + 4)) then
       p_par = self%patch%elts(p_par0+1)%neigh(s0+1)
    else if (s0 == modulo(c0 + 1, 4) + 4 .or. s0 == modulo(c0 - 1, 4) + 4) then
       p_par = self%patch%elts(p_par0+1)%neigh(s_par+1)
    else ! neighbour patch on same parent
       p_par = p_par0
    end if

    c = ngb_chd_idx(c0, s0)
  contains
    integer function ngb_chd_idx(c, s)
      implicit none
      integer :: c, s

      if (s < 4) then
         ngb_chd_idx = modulo(-c + 2*s + 1, 4)
         return
      else
         ngb_chd_idx = modulo(c + 2, 4)
         return
      end if
    end function ngb_chd_idx
  end subroutine find_neigh_patch2_Domain

  integer function find_neigh_patch_Domain (self, p_par, c, s_chd)
    ! Finds the neighbour at side s of c-th child of parent patch p_par through p_par
    implicit none
    type(Domain) :: self
    integer      :: c, p_par, s_chd

    integer :: c1, n_par, n_side, p_chd, p_par1, s_help, s_par, typ

    call find_neigh_patch2_Domain (self, p_par, c, s_chd, n_par, c1)

    s_par = par_side(c, s_chd)

    if (n_par > 0) then
       find_neigh_patch_Domain = self%patch%elts(n_par+1)%children(c1+1)
    else if (n_par < 0) then ! bdry patch
       n_par = -n_par
       typ = self%bdry_patch%elts(n_par+1)%side
       if (typ < 1) then
          find_neigh_patch_Domain = 0
          return
       end if

       if (s_chd + 1 == typ) then ! side
          find_neigh_patch_Domain = find_neigh_bdry_patch_Domain(self, p_par, c, s_chd)
       else ! patch corner, but domain side
          if (s_chd + 1 == typ + 4) then ! s_chd after typ (clockwise sense)
             s_help = modulo(typ, 4)
          else if (s_chd + 1 == modulo(typ-2, 4) + 5) then ! s_chd before typ
             s_help = modulo(typ-2, 4)
          end if
          call find_neigh_patch2_Domain (self, p_par, c, s_help, p_par1, c1)
          find_neigh_patch_Domain = find_neigh_bdry_patch_Domain(self, p_par1, c1, typ-1)
       end if
    end if
  end function find_neigh_patch_Domain

  integer function par_side (c, s)
    implicit none
    integer :: c, s

    if (s == modulo(c+1, 4) + 4) then
       par_side = modulo(c+1, 4)
       return
    else
       if (s == modulo(c-1, 4) + 4) then
          par_side = c
          return
       else
          par_side = s
          return
       end if
    end if
  end function par_side

  subroutine get_offs_Domain5 (self, p, offs, dims, inner_patch)
    implicit none
    type(Domain)                    :: self
    integer                         :: p
    integer, dimension(N_BDRY+1)    :: offs
    integer, dimension(2,N_BDRY+1)  :: dims
    logical, dimension(4), optional :: inner_patch

    integer :: i, n

    offs = -1
    dims =  0
    offs(1) = self%patch%elts(p+1)%elts_start

    if (present(inner_patch)) inner_patch = .false.

    do i = 1, N_BDRY
       n = self%patch%elts(p+1)%neigh(i)
       if (n > 0) then
          offs(i+1) = self%patch%elts(n+1)%elts_start
          dims(:,i+1) = PATCH_SIZE
          if (present(inner_patch) .and. i <= 4) inner_patch(i) = .true.
       else
          if (n < 0) then
             offs(i+1) = self%bdry_patch%elts(-n+1)%elts_start
             dims(:,i+1) = sides_dims(:,abs(self%bdry_patch%elts(-n+1)%side) + 1)
          end if
       end if
    end do
  end subroutine get_offs_Domain5

  subroutine get_offs_Domain (self, p, offs, dims, inner_bdry)
    implicit none
    type(Domain)                    :: self
    integer                         :: p
    integer, dimension(N_BDRY+1)    :: offs
    integer, dimension(2,N_BDRY+1)  :: dims
    logical, dimension(4), optional :: inner_bdry

    integer :: i, n

    offs = -1
    dims = 0
    offs(1) = self%patch%elts(p+1)%elts_start

    do i = 1, N_BDRY
       n = self%patch%elts(p+1)%neigh(i)
       if (n > 0) then
          offs(i+1) = self%patch%elts(n+1)%elts_start
          dims(:,i+1) = PATCH_SIZE
          if (present(inner_bdry) .and. i <= 4) inner_bdry(i) = .true. 
       else
          if (n < 0) then
             offs(i+1) = self%bdry_patch%elts(-n+1)%elts_start
             dims(:,i+1) = sides_dims(:,abs(self%bdry_patch%elts(-n+1)%side) + 1)
             if (present(inner_bdry) .and. i <= 4) inner_bdry(i) = self%bdry_patch%elts(-n+1)%side < 0
          end if
       end if
    end do
  end subroutine get_offs_Domain
end module domain_mod
