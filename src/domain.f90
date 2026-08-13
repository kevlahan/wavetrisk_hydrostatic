module domain_mod

  use kind_mod,   only : dp

  use shared_mod, only : AT_NODE, AT_EDGE, BDRY_THICKNESS, b_vert, b_vert_mass, compressible, &
       EAST, EDGE, grav_accel, IMINUS, IJMINUS, IMINUSJPLUS, JMINUS, &
       N_GLO_DOMAIN, N_BDRY, N_DOMAIN, NORTH, NORTHEAST, ORIGIN, RT, DG, SOUTHEAST, S_VELO, TRIAG, UP, WEST, &
       a_vert, a_vert_mass, p_0, ref_density, max_level, max_depth, min_level, mode_split, &
       split_mean_perturbation, scalars, zlevels, zmin

  use patch_mod, only : BDRY_PATCH, Patch, PATCH_SIZE
  use arch_mod,  only : n_glo_block, n_process, rank

  use dyn_arrays, only : Areas_Array, Bdry_Patch_Array, Coord_Array, Float_Array, Int_Array, &
       Iu_Wgt_Array, Logical_Array, Overl_Area_Array, Patch_Array, RF_Wgt_Array, Topo_Array, &
       append, extend, init

  implicit none

  private
  public :: Domain, Float_Field, Int_Field, Logical_Field
  public :: init_field, sides_dims, chd_offs
  public :: grid, topography, topography_data, sso_param, exner_fun, horiz_flux, penal_node, penal_edge, Kv, Kt, tke
  public :: wav_tke, Laplacian_scalar, Laplacian_vector, sol, sol_mean, trend, wav_coeff
  public :: diag, mass, mass1, h_flux, h_mflux, dmass, dtemp, dscalar, scalar, scalar_2d, temp, temp1
  public :: velo, velo1, velo2, velo_2d, dvelo, dvelo_2d
  public :: mean_m, mean_t, Laplacian,  bernoulli, divu, exner, ke, qe, vort,  wc_s, wc_u
  public :: init_Domain, add_bdry_patch_Domain, add_patch_Domain, extend_Domain, set_neigh_Domain
  public :: idx, idx2, idx__fast, idx_hex, idx_hex_LORT, idx_hex_LORT2, idx_hex_UPLT, idx_hex_UPLT2, ed_idx, id_edge, tri_idx
  public :: nidx, is_penta, find_neigh_bdry_patch_domain, find_neigh_patch2_domain, find_neigh_patch_domain
  public :: par_side, get_offs_Domain, get_offs_Domain5

  
  integer, parameter :: sides_dims(2,N_BDRY+1) = reshape ( [PATCH_SIZE, PATCH_SIZE, PATCH_SIZE, &
       BDRY_THICKNESS, BDRY_THICKNESS, PATCH_SIZE, PATCH_SIZE, &
       BDRY_THICKNESS, BDRY_THICKNESS, PATCH_SIZE, BDRY_THICKNESS, &
       BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS, &
       BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS], [2, 9])

  integer, parameter :: chd_offs(2,4) = reshape ([PATCH_SIZE/2, PATCH_SIZE/2, PATCH_SIZE/2, 0, 0, 0, 0, PATCH_SIZE/2], [2, 4])


  interface init_field
     procedure :: init_Float_Field, init_Int_Field, init_Logical_Field
  end interface init_field

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
     type(Int_Array)              :: send_pa_all
     type(Int_Array), allocatable :: recv_pa(:), send_conn(:)
     type(Int_Array), allocatable :: pack(:,:), unpk(:,:)
     type(Int_Array), allocatable :: src_patch(:,:)

     ! Physical quantities
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
     integer                        :: pos
     integer                        :: bdry_tag = -1
     logical                        :: bdry_uptodate
     type(Float_Array), allocatable :: data(:)
  end type Float_Field

  type Int_Field
     integer                      :: pos
     integer                      :: bdry_tag = -1
     logical                      :: bdry_uptodate
     type(Int_Array), allocatable :: data(:)
  end type Int_Field

  type Logical_Field
     integer                          :: pos
     integer                          :: bdry_tag = -1
     logical                          :: bdry_uptodate
     type(Logical_Array), allocatable :: data(:)
  end type Logical_Field

  type(Domain), allocatable, target :: grid(:)

  type(Float_Field),                              target :: topography
  type(Topo_Array),  dimension(:,:), allocatable         :: topography_data
  type(Float_Field), dimension(4),                target :: sso_param
  type(Float_Field), dimension(:),   allocatable, target :: exner_fun, horiz_flux, penal_node, penal_edge
  type(Float_Field), dimension(:),   allocatable, target :: Kv, Kt, tke, wav_tke
  type(Float_Field), dimension(:),   allocatable, target :: Laplacian_scalar, Laplacian_vector
  type(Float_Field), dimension(:,:), allocatable, target :: sol, sol_mean, trend
  type(Float_Field), dimension(:,:), allocatable, target :: wav_coeff

  real(dp), dimension(:), pointer :: diag, mass, mass1, h_flux, h_mflux
  real(dp), dimension(:), pointer :: dmass, dtemp, dscalar, scalar, scalar_2d, temp, temp1
  real(dp), dimension(:), pointer :: velo, velo1, velo2, velo_2d, dvelo, dvelo_2d
  real(dp), dimension(:), pointer :: mean_m, mean_t
  real(dp), dimension(:), pointer :: Laplacian
  real(dp), dimension(:), pointer :: bernoulli, divu, exner, ke, qe, vort
  real(dp), dimension(:), pointer :: wc_s, wc_u

  
contains

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialization subroutines !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_Float_Field (self, pos)
    implicit none
    type(Float_Field), intent(out) :: self
    integer,           intent(in) :: pos

    self%bdry_uptodate = .false.
    self%pos = pos

    allocate (self%data(n_domain(rank+1)))
  end subroutine init_Float_Field


  subroutine init_Int_Field (self, pos)
    implicit none
    type(Int_Field), intent(out) :: self
    integer,         intent(in)  :: pos

    self%bdry_uptodate = .false.
    self%pos = pos

    allocate (self%data(n_domain(rank+1)))
  end subroutine init_Int_Field

  
  subroutine init_Logical_Field (self, pos)
    implicit none
    type(Logical_Field), intent(out) :: self
    integer,             intent(in)  :: pos

    self%bdry_uptodate = .false.
    self%pos = pos

    allocate (self%data(n_domain(rank+1)))
  end subroutine init_Logical_Field
  
  
  subroutine init_Domain (self)
    implicit none
    type(Domain), intent(out) :: self

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

    allocate (self%recv_pa(n_glo_block))
    allocate (self%send_conn(n_glo_block))

    allocate (self%pack(AT_NODE:AT_EDGE,n_glo_block))
    allocate (self%unpk(AT_NODE:AT_EDGE,n_glo_block))

    do i = 1, n_glo_block
       call init (self%send_conn(i), 0)
    end do

    call init (self%send_pa_all, 0)

    do i = 1, n_glo_block
       call init (self%recv_pa(i), 0)
    end do

    do k = AT_NODE, AT_EDGE
       do i = 1, n_glo_block
          call init (self%pack(k,i), 0)
          call init (self%unpk(k,i), 0)
       end do
    end do

    self%pole_master = .false.
  end subroutine init_Domain


  integer function add_patch_Domain (self, level)
    ! Add new patch to the domain
    implicit none
    type(Domain), intent(inout) :: self

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
    type(Domain), intent(inout) :: self
    integer,      intent(in)    :: side

    integer :: p

    p = self%bdry_patch%length

    call append (self%bdry_patch, Bdry_Patch (self%node%length, side, 0))

    call extend_Domain (self, BDRY_THICKNESS * PATCH_SIZE)

    add_bdry_patch_Domain = p
  end function add_bdry_patch_Domain

  
  subroutine extend_Domain (self, num)
    implicit none
    type(Domain), intent(inout) :: self
    integer,      intent(in)    :: num

    integer  :: d, k, v
    real(dp) :: def_val, dz, z

    d = self%id + 1

    call extend (self%node, num, ORIGIN)

    ! Atmosphere/ocean layers
    do k = 1, zlevels
       ! Set reasonable default values for new boundary patches to avoid NaN if variable undefined in boundary
       if (compressible) then
          def_val = a_vert_mass(k) + b_vert_mass(k) * p_0 / grav_accel
       else
          dz     = b_vert_mass(k) * max_depth
          z      = 0.5 * (b_vert(k) + b_vert(k-1)) * max_depth
          def_val = ref_density * dz
       end if

       do v = scalars(1), scalars(2)
          if (split_mean_perturbation) then 
             call extend (sol(v,k)%data(d),      num, 0.0_dp)     
             call extend (sol_mean(v,k)%data(d), num, def_val) 
          else
             call extend (sol(v,k)%data(d),      num, def_val) 
             call extend (sol_mean(v,k)%data(d), num, 0.0_dp)     
          end if
       end do
       call extend (sol(S_VELO,k)%data(d),      EDGE * num, 0.0_dp)
       call extend (sol_mean(S_VELO,k)%data(d), EDGE * num, 0.0_dp)
    end do

    ! Free surface
    if (mode_split) then
       do v = scalars(1), scalars(2)
          if (split_mean_perturbation) then 
             call extend (sol(v,k)%data(d),      num, 0.0_dp)     
             call extend (sol_mean(v,k)%data(d), num, 0.0_dp) 
          else
             call extend (sol(v,k)%data(d),      num, 0.0_dp) 
             call extend (sol_mean(v,k)%data(d), num, 0.0_dp)     
          end if
       end do
       call extend (sol(S_VELO,k)%data(d),      EDGE * num, 0.0_dp)
       call extend (sol_mean(S_VELO,k)%data(d), EDGE * num, 0.0_dp)
    end if

    ! Soil layers
    do k = zmin, 0
       do v = scalars(1), scalars(2)
          if (split_mean_perturbation) then 
             call extend (sol(v,k)%data(d),      num, 0.0_dp)     
             call extend (sol_mean(v,k)%data(d), num, 0.0_dp) 
          else
             call extend (sol(v,k)%data(d),      num, 0.0_dp) 
             call extend (sol_mean(v,k)%data(d), num, 0.0_dp)     
          end if
       end do
       call extend (sol(S_VELO,k)%data(d),      EDGE * num, 0.0_dp)
       call extend (sol_mean(S_VELO,k)%data(d), EDGE * num, 0.0_dp)
    end do
  end subroutine extend_Domain

  
  subroutine set_neigh_Domain (self, s, id, rot)
    implicit none
    type(Domain), intent(inout) :: self
    integer,      intent(in)    :: id, rot, s

    self%neigh(s)     = id
    self%neigh_rot(s) = rot
  end subroutine set_neigh_Domain

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Index and neighbour routines !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function idx (i0, j0, offs, dims) result(id)
    ! Given regular-array coordinates (i0,j0), offset array offs,
    ! and domain dimensions dims, return the corresponding grid index.
    implicit none

    integer, intent(in) :: i0, j0
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: id
    integer :: i, j

    i = i0
    j = j0

    if (i < 0) then
       if (j < 0) then
          id = offs(IJMINUS+1)                                      &
               + (j + dims(2,IJMINUS+1))*dims(1,IJMINUS+1)          &
               + i + dims(1,IJMINUS+1)

       else if (j >= PATCH_SIZE) then
          id = offs(IMINUSJPLUS+1)                                 &
               + (j - PATCH_SIZE)*dims(1,IMINUSJPLUS+1)            &
               + i + dims(1,IMINUSJPLUS+1)

       else
          id = offs(WEST+1)                                        &
               + j*dims(1,WEST+1)                                  &
               + i + dims(1,WEST+1)
       end if

    else if (i >= PATCH_SIZE) then
       i = i - PATCH_SIZE

       if (j >= PATCH_SIZE) then
          j = j - PATCH_SIZE

          id = offs(NORTHEAST+1)                                   &
               + j*dims(1,NORTHEAST+1)                             &
               + i

       else if (j < 0) then
          id = offs(SOUTHEAST+1)                                   &
               + (j + dims(2,SOUTHEAST+1))*dims(1,SOUTHEAST+1)     &
               + i

       else
          id = offs(EAST+1)                                        &
               + j*dims(1,EAST+1)                                  &
               + i
       end if

    else
       if (j < 0) then
          id = offs(JMINUS+1)                                      &
               + (j + dims(2,JMINUS+1))*dims(1,JMINUS+1)           &
               + i

       else if (j >= PATCH_SIZE) then
          id = offs(NORTH+1)                                       &
               + (j - PATCH_SIZE)*dims(1,NORTH+1)                  &
               + i

       else
          id = offs(1) + PATCH_SIZE*j + i
       end if
    end if
  end function idx

  
  pure function idx__fast(i, j, offs) result(id)
    implicit none

    integer, intent(in) :: i, j, offs
    integer             :: id

    id = PATCH_SIZE*j + i + offs
  end function idx__fast

  
  pure function id_edge (id) result(edges)
    ! Returns indices of the three edges associated with node id
    implicit none

    integer, intent(in) :: id
    integer             :: edges(EDGE)        ! EDGE must be 3 here

    edges = EDGE*id + [ RT, DG, UP ] + 1
  end function id_edge

  
  pure function idx_hex (i, j, offs, dims) result(hex_idx)
    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)
    integer             :: hex_idx(2*EDGE+1)

    hex_idx = [ &
         idx(i,   j,   offs, dims), &
         idx(i+1, j,   offs, dims), &
         idx(i+1, j+1, offs, dims), &
         idx(i,   j+1, offs, dims), &
         idx(i-1, j,   offs, dims), &
         idx(i-1, j-1, offs, dims), &
         idx(i,   j-1, offs, dims)  &
         ]
  end function idx_hex

  
  pure function idx_hex_LORT (i, j, offs, dims) result(hex_idx)
    ! Return the indices of the six child hexagons overlapping
    ! the LORT parent triangle.
    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: hex_idx(2*EDGE)

    hex_idx = [ &
         idx(i,   j,   offs, dims), &
         idx(i+1, j,   offs, dims), &
         idx(i+1, j+1, offs, dims), &
         
         idx(i+2, j,   offs, dims), &
         idx(i+2, j+2, offs, dims), &
         idx(i+2, j+1, offs, dims)  &
         ]
  end function idx_hex_LORT

  
  pure function idx_hex_LORT2 (i, j, offs, dims) result(hex_idx)
    ! Return the indices of the three hexagons overlapping the LORT triangle.
    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: hex_idx(EDGE)

    hex_idx = [ &
         idx(i,   j,   offs, dims), &
         idx(i+1, j,   offs, dims), &
         idx(i+1, j+1, offs, dims)  &
         ]
  end function idx_hex_LORT2

  
  pure function idx_hex_UPLT (i, j, offs, dims) result(hex_idx)
    ! Return the indices of the six child hexagons overlapping
    ! the LORT parent triangle.
    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: hex_idx(2*EDGE)

    hex_idx = [ &
         idx(i,   j,   offs, dims), &
         idx(i+1, j+1, offs, dims), &
         idx(i,   j+1, offs, dims), &
         idx(i,   j+2, offs, dims), &
         idx(i+2, j+2, offs, dims), &
         idx(i+1, j+2, offs, dims)  &
         ]
  end function idx_hex_UPLT

  
  pure function idx_hex_UPLT2 (i, j, offs, dims) result(hex_idx)
    ! Return the indices of the three hexagons overlapping the UPLT triangle.
    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: hex_idx(EDGE)

    hex_idx = [ &
         idx(i,   j,   offs, dims), &
         idx(i+1, j+1, offs, dims), &
         idx(i,   j+1, offs, dims)  &
         ]
  end function idx_hex_UPLT2

  
  pure function tri_idx (i, j, tri, offs, dims) result(id)
    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: tri(3)
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: id

    id = TRIAG*idx(i + tri(1), j + tri(2), offs, dims) + tri(3)
  end function tri_idx

  
  pure function nidx (i, j, s, offs, dims) result(id)
    implicit none

    integer, intent(in) :: i, j, s
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: id

    id = offs(s+1) + j*dims(1,s+1) + i
  end function nidx

  
  pure function idx2 (i, j, noffs, offs, dims) result(id)
    ! Return the index of node (i+noffs(1), j+noffs(2)).
    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: noffs(2)
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: id

    id = idx(i + noffs(1), j + noffs(2), offs, dims)
  end function idx2

  
  pure function ed_idx (i, j, ed, offs, dims) result(id)
    ! Return the edge index.
    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: ed(3)
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: id

    id = EDGE*idx(i + ed(1), j + ed(2), offs, dims) + ed(3)
  end function ed_idx

  
  pure logical function is_penta (dom, p, s) result(penta)
    implicit none

    type(Domain), intent(in) :: dom
    integer,      intent(in) :: p, s

    integer :: n, side

    penta = .false.

    n = dom%patch%elts(p+1)%neigh(s+1)

    if (n < 0) then
       n = -n
       side = dom%bdry_patch%elts(n+1)%side

       if (side > 0) penta = dom%penta(side)
    end if
  end function is_penta

  
  pure function find_neigh_patch_Domain (self, p_par, c, s_chd) result(neigh)
    ! Find the neighbour of the c-th child of parent patch p_par
    ! across child side s_chd.
    implicit none

    type(Domain), intent(in) :: self
    integer,      intent(in) :: p_par, c, s_chd

    integer :: neigh
    integer :: c1, n_par, p_par1
    integer :: s_help, typ

    neigh = 0

    call find_neigh_patch2_Domain( &
         self, p_par, c, s_chd, n_par, c1)

    if (n_par > 0) then
       neigh = self%patch%elts(n_par+1)%children(c1+1)

    else if (n_par < 0) then
       ! Boundary patch.
       n_par = -n_par
       typ = self%bdry_patch%elts(n_par+1)%side

       if (typ < 1) return

       if (s_chd+1 == typ) then
          ! Domain side.
          neigh = find_neigh_bdry_patch_Domain( &
               self, p_par, c, s_chd)

       else
          ! Patch corner lying on a domain side.
          if (s_chd+1 == typ+4) then
             ! s_chd follows typ in the clockwise direction.
             s_help = modulo(typ, 4)

          else if (s_chd+1 == modulo(typ-2, 4)+5) then
             ! s_chd precedes typ.
             s_help = modulo(typ-2, 4)

          else
             ! No recognized corner relationship.
             return
          end if

          call find_neigh_patch2_Domain( &
               self, p_par, c, s_help, p_par1, c1)

          neigh = find_neigh_bdry_patch_Domain( &
               self, p_par1, c1, typ-1)
       end if
    end if
  end function find_neigh_patch_Domain

  
  pure function find_neigh_bdry_patch_Domain (self, p_par, c, s) result(neigh)
    implicit none

    type(Domain), intent(in) :: self
    integer,      intent(in) :: p_par, c, s

    integer :: neigh
    integer :: p_chd, p_par1, c1

    neigh = 0

    if (p_par <= 0) return

    p_chd = self%patch%elts(p_par+1)%children(c+1)

    if (p_chd > 0) then
       neigh = self%patch%elts(p_chd+1)%neigh(s+1)
    end if

    if (s >= 4) return

    if (neigh == 0) then
       call find_neigh_patch2_Domain ( &
            self, p_par, c, modulo(s+1, 4), p_par1, c1)

       if (p_par1 > 0) then
          p_chd = self%patch%elts(p_par1+1)%children(c1+1)

          if (p_chd > 0) then
             neigh = self%patch%elts(p_chd+1)%neigh( &
                  modulo(s-1, 4) + 5)
          end if
       end if
    end if

    if (neigh == 0) then
       call find_neigh_patch2_Domain ( &
            self, p_par, c, modulo(s-1, 4), p_par1, c1)

       if (p_par1 > 0) then
          p_chd = self%patch%elts(p_par1+1)%children(c1+1)

          if (p_chd > 0) then
             neigh = self%patch%elts(p_chd+1)%neigh(s+5)
          end if
       end if
    end if
  end function find_neigh_bdry_patch_Domain

  
  pure subroutine find_neigh_patch2_Domain (self, p_par0, c0, s0, p_par, c)
    ! For the patch given by the c0-th child of p_par0, find the
    ! neighbour across side s0. Return it as the c-th child of p_par.
    implicit none

    type(Domain), intent(in)  :: self
    integer,      intent(in)  :: p_par0, c0, s0
    integer,      intent(out) :: p_par, c

    integer :: s_par

    s_par = par_side(c0, s0)

    if (s0 == c0 .or.                            &
         s0 == modulo(c0+1, 4) .or.               &
         s0 == c0+4) then

       p_par = self%patch%elts(p_par0+1)%neigh(s0+1)

    else if (s0 == modulo(c0+1, 4)+4 .or.        &
         s0 == modulo(c0-1, 4)+4) then

       p_par = self%patch%elts(p_par0+1)%neigh(s_par+1)

    else
       ! Neighbour patch has the same parent.
       p_par = p_par0
    end if

    c = ngb_chd_idx(c0, s0)

  contains

    pure integer function ngb_chd_idx(c0, s0) result(c_neigh)
      implicit none

      integer, intent(in) :: c0, s0

      if (s0 < 4) then
         c_neigh = modulo(-c0 + 2*s0 + 1, 4)
      else
         c_neigh = modulo(c0 + 2, 4)
      end if
    end function ngb_chd_idx

  end subroutine find_neigh_patch2_Domain

  
  pure function par_side(c, s) result(side)
    implicit none

    integer, intent(in) :: c, s
    integer             :: side
    
    if (s == modulo(c+1, 4) + 4) then
       side = modulo(c+1, 4)
    else if (s == modulo(c-1, 4) + 4) then
       side = c
    else
       side = s
    end if
  end function par_side

  
  pure subroutine get_offs_Domain5 (self, p, offs, dims, inner_patch)
    implicit none

    type(Domain), intent(in)            :: self
    integer,      intent(in)            :: p
    integer,      intent(out)           :: offs(N_BDRY+1)
    integer,      intent(out)           :: dims(2,N_BDRY+1)
    logical,      intent(out), optional :: inner_patch(4)

    integer :: i, n
    integer :: bdry, side

    offs = -1
    dims = 0

    offs(1) = self%patch%elts(p+1)%elts_start

    do i = 1, N_BDRY
       n = self%patch%elts(p+1)%neigh(i)

       if (n > 0) then
          offs(i+1)   = self%patch%elts(n+1)%elts_start
          dims(:,i+1) = PATCH_SIZE

       else if (n < 0) then
          bdry = -n
          side = self%bdry_patch%elts(bdry+1)%side

          offs(i+1)   = self%bdry_patch%elts(bdry+1)%elts_start
          dims(:,i+1) = sides_dims(:,abs(side)+1)
       end if
    end do

    if (present(inner_patch)) then
       inner_patch = .false.

       do i = 1, min(N_BDRY, size(inner_patch))
          n = self%patch%elts(p+1)%neigh(i)
          inner_patch(i) = n /= 0
       end do
    end if
  end subroutine get_offs_Domain5

  
  pure subroutine get_offs_Domain(self, p, offs, dims, inner_bdry)
    implicit none

    type(Domain), intent(in)            :: self
    integer,      intent(in)            :: p
    integer,      intent(out)           :: offs(N_BDRY+1)
    integer,      intent(out)           :: dims(2,N_BDRY+1)
    logical,      intent(out), optional :: inner_bdry(4)

    integer :: i, n
    integer :: bdry, side

    offs = -1
    dims = 0
    offs(1) = self%patch%elts(p+1)%elts_start

    if (present(inner_bdry)) inner_bdry = .false.

    do i = 1, N_BDRY
       n = self%patch%elts(p+1)%neigh(i)

       if (n > 0) then
          offs(i+1)   = self%patch%elts(n+1)%elts_start
          dims(:,i+1) = PATCH_SIZE

       else if (n < 0) then
          bdry = -n
          side = self%bdry_patch%elts(bdry+1)%side

          offs(i+1)   = self%bdry_patch%elts(bdry+1)%elts_start
          dims(:,i+1) = sides_dims(:,abs(side)+1)
       end if
    end do

    if (present(inner_bdry)) then
       do i = 1, min(N_BDRY, size(inner_bdry))
          n = self%patch%elts(p+1)%neigh(i)

          if (n > 0) then
             inner_bdry(i) = .true.
          else if (n < 0) then
             side = self%bdry_patch%elts(-n+1)%side
             inner_bdry(i) = side < 0
          end if
       end do
    end if
  end subroutine get_offs_Domain

  
end module domain_mod
