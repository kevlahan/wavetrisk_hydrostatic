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
     integer                    :: id
     integer, dimension(N_BDRY) :: neigh, neigh_rot
     integer, dimension(2)      :: neigh_over_pole

     logical, dimension(2)      :: pole_master
     logical, dimension(N_BDRY) :: penta
     
     type(Patch_Array) :: patch
     
     type(Bdry_Patch_Array) :: bdry_patch
     
     type(Coord_Array) :: ccentre, midpt, node
     
     type(Int_Array)                                          :: level, mask_e, mask_n, neigh_pa_over_pole, send_pa_all
     type(Int_Array), dimension(:), allocatable               :: lev
     type(Int_Array), dimension(N_GLO_DOMAIN)                 :: recv_pa, send_conn
     type(Int_Array), dimension(AT_NODE:AT_EDGE,N_GLO_DOMAIN) :: pack, unpk
     type(Int_Array), dimension(:,:), allocatable             :: src_patch

     type(Float_Array) :: coriolis, len, pedlen, triarea

     type(Areas_Array) :: areas

     type(Float_Array) :: surf_press  ! surface pressure (compressible) or surface Lagrange multiplier (incompressible)
     type(Float_Array) :: press       ! pressure (compressible case) or Lagrange multiplier (incompressible case)
     type(Float_Array) :: surf_geopot ! surface geopotential
     type(Float_Array) :: geopot      ! geopotential
     type(Float_Array) :: u_zonal     ! zonal velocity
     type(Float_Array) :: v_merid     ! meridional velocity
     type(Float_Array) :: adj_mass    ! mass in adjacent vertical cell
     type(Float_Array) :: adj_temp    ! temp in adjacent vertical cell
     type(Float_Array) :: adj_geopot  ! geopotential in adjacent vertical cell
     type(Float_Array) :: vort        ! vorticity
     type(Float_Array) :: bernoulli   ! Bernoulli function
     type(Float_Array) :: bern_fast   ! Bernoulli function
     type(Float_Array) :: bern_slow   ! Bernoulli function (slow part)
     type(Float_Array) :: divu        ! divergence of velocity
     type(Float_Array) :: qe          ! 

     type(Overl_Area_Array) :: overl_areas
     
     type(Iu_Wgt_Array)     :: I_u_wgt
     type(RF_Wgt_Array)     :: R_F_wgt
  end type Domain

  type Float_Field
     integer                                      :: pos
     logical                                      :: bdry_uptodate
     type(Float_Array), dimension(:), allocatable :: data
  end type Float_Field

  type(Domain), dimension(:), allocatable, target        :: grid
  type(Float_Field), dimension(:,:), allocatable, target :: sol, trend, wav_coeff, trend_wav_coeff
  type(Float_Field), dimension(:), allocatable, target   :: bernoulli_fast, exner_fun, horiz_flux

  !note that the theta in the DYNAMICO paper is in fact theta^b (buoyancy)
  !we have theta^b=(theta_r-theta_k)/theta_r where theta_r is the reference potential temperature
  !and theta_k is the potential temperature
  !therefore 1-theta in the DYNAMICO paper is really 1-theta^b=theta_k/theta_r=theta'_k which
  !is what they are solving for in the code (and changes some of the equations in the paper, more specifically
  !the equation below (18) and the equation below (24))
  !we will solve for theta'_k and Theta'_k
  real(8), dimension(:), pointer :: mass, dmass, h_mflux
  real(8), dimension(:), pointer :: temp, dtemp, h_tflux
  real(8), dimension(:), pointer :: velo, dvelo
  real(8), dimension(:), pointer :: bernoulli, divu, exner, qe, vort
  real(8), dimension(:), pointer :: wc_u, wc_m, wc_t
contains
  subroutine init_Float_Field (self, pos)
    type(Float_Field) :: self
    integer           :: pos

    self%bdry_uptodate = .False.
    self%pos = pos

    allocate(self%data(n_domain(rank+1)))
  end subroutine init_Float_Field

  subroutine init_domain_mod
    integer :: i, v
    logical :: initialized = .False.

    if (initialized) return ! initialize only once

    call init_shared_mod
    call init_sphere_mod
    call init_patch_mod
    call init_arch_mod

    sides_dims = reshape((/PATCH_SIZE, PATCH_SIZE, PATCH_SIZE, &
         BDRY_THICKNESS, BDRY_THICKNESS, PATCH_SIZE, PATCH_SIZE, &
         BDRY_THICKNESS, BDRY_THICKNESS, PATCH_SIZE, BDRY_THICKNESS, &
         BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS, &
         BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS/), (/2, 9/))

    chd_offs = reshape((/PATCH_SIZE/2, PATCH_SIZE/2, PATCH_SIZE/2, 0, 0, 0, &
         0, PATCH_SIZE/2/), (/2, 4/))

    initialized = .True.
  end subroutine init_domain_mod

  subroutine apply_onescale_d (routine, dom, l, zlev, st, en)
    external     :: routine
    type(Domain) :: dom
    integer      :: en, l, st, zlev

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_onescale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_onescale_d
  
  subroutine apply_onescale_to_patch (routine, dom, p, zlev, st, en)
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

  subroutine apply_onescale_to_patch__int (routine, dom, p, zlev, st, en, ival)
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
    external :: routine
    integer  :: en, l, st, zlev

    integer :: d, k

    do d = 1, size(grid)
       do k = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch2(routine, grid(d), grid(d)%lev(l)%elts(k), zlev, st, en)
       end do
    end do
  end subroutine apply_onescale2

  subroutine apply_onescale_to_patch5 (routine, dom, p, zlev, st, en)
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
          call routine(dom, p, i - 1, j - 1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch5

  subroutine apply_onescale_to_patch2 (routine, dom, p, zlev, st, en)
    external     :: routine
    type(Domain) :: dom
    integer      :: en, p, st, zlev

    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call get_offs_Domain(dom, p, offs, dims)

    do j = st + 1, PATCH_SIZE + en
       do i = st + 1, PATCH_SIZE + en
          call routine(dom, p, i - 1, j - 1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch2

  subroutine apply_interscale_d (routine, dom, l, zlev, st, en)
    external     :: routine
    type(Domain) :: dom
    integer      :: en, l, st, zlev

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_interscale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_interscale_d
  
  subroutine apply_interscale (routine, l, zlev, st, en)
    external :: routine
    integer  :: en, l, st, zlev

    integer :: d

    do d = 1, size(grid)
       call apply_interscale_d (routine, grid(d), l, zlev, st, en)
    end do
  end subroutine apply_interscale

   subroutine apply_interscale_d2 (routine, dom, l, zlev, st, en)
    external     :: routine
    type(Domain) :: dom
    integer      :: en, l, st, zlev

    integer :: k

    do k = 1, dom%lev(l)%length
       call apply_interscale_to_patch22 (routine, dom, dom%lev(l)%elts(k), zlev, st, en)
    end do
  end subroutine apply_interscale_d2

  subroutine apply_interscale_to_patch (routine, dom, p_par, zlev, st, en)
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
       if (p_chd .eq. 0) cycle

       call get_offs_Domain (dom, p_chd, offs_chd, dims_chd, inner_bdry)

       bdry = (/en, en, st, st/)

       where (inner_bdry) bdry = 0

       do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
          j_chd = (j - 1)*2
          j_par = j - 1 + chd_offs(2,c)
          do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
             i_chd = (i - 1)*2
             i_par = i - 1 + chd_offs(1,c)
             call routine(dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
          end do
       end do
    end do
  end subroutine apply_interscale_to_patch

   subroutine apply_interscale_to_patch2 (routine, dom, p_par, zlev, st, en)
    external     :: routine
    type(Domain) :: dom
    integer      :: en, st, p_par, zlev

    integer                          :: c, i, i_chd, i_par, j, j_chd, j_par, p_chd
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry
    !TODO{uncomment & test}  if (dom%patch%elts(p_par+1)%active .eq. NONE) return

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do c = 1, N_CHDRN
       p_chd = dom%patch%elts(p_par+1)%children(c)

       if (p_chd .eq. 0) cycle

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
       if (p_chd .eq. 0) cycle

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
    external     :: routine
    type(Domain) :: dom
    integer      :: c, en, p_par, st, zlev
    
    integer                          :: i, i_chd, i_par, j, j_chd, j_par, p_chd
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry
    !TODO{uncomment & test}  if (dom%patch%elts(p_par+1)%active .eq. NONE) return

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

  function ed_idx (i, j, ed, offs, dims)
    !return edge index
    integer                        :: ed_idx
    integer                        :: i, j
    integer, dimension(3)          :: ed
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    ed_idx = EDGE*idx(i + ed(1), j + ed(2), offs, dims) + ed(3)
  end function ed_idx

  subroutine apply_to_pole_d (routine, dom, l, zlev, ival, to_all)
    external     :: routine
    type(Domain) :: dom
    integer      :: ival, l, zlev
    logical      :: to_all

    integer                        :: c, l_cur, p, p_par
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do c = SOUTHEAST, NORTHWEST, 2

       if (.not. dom%pole_master(c/2-2) .and. .not. to_all) cycle
       if (.not. dom%penta(c)) cycle
       if (.not. dom%neigh(c) .eq. POLE) cycle

       p=1
       do while (p .gt. 0)
          p_par = p
          p = dom%patch%elts(p_par+1)%children(c-4)

          if (.not. l .eq. NONE) then
             l_cur = dom%patch%elts(p_par+1)%level
             if (l_cur .lt. l) then
                cycle
             else
                if (l_cur .gt. l) exit
             end if
          end if

          call get_offs_Domain (dom, p_par, offs, dims)

          if (c .eq. NORTHWEST) then
             call routine (dom, p_par, 0, PATCH_SIZE, zlev, offs, dims, ival) ! NORTHPOLE
          else
             call routine (dom, p_par, PATCH_SIZE, 0, zlev, offs, dims, ival) ! SOUTHPOLE
          end if

       end do
    end do
  end subroutine apply_to_pole_d

  subroutine apply_to_pole (routine, l, zlev, ival, to_all)
    external :: routine
    integer  :: ival, l, zlev
    logical  :: to_all

    integer                        :: c, d, l_cur, p, p_par
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do d = 1, size(grid)
       call apply_to_pole_d (routine, grid(d), l, zlev, ival, to_all)
    end do
  end subroutine apply_to_pole

  subroutine apply_to_penta_d (routine, dom, l, zlev)
    external     :: routine
    type(Domain) :: dom
    integer      :: l, zlev

    integer                        :: c, d, l_cur, p, p_par
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do c = NORTHEAST, NORTHWEST
       if (.not. dom%penta(c)) cycle

       p = 1
       do while (p .gt. 0)
          p_par = p
          p = dom%patch%elts(p_par+1)%children(c-4)
          if (.not. l .eq. NONE) then
             l_cur = dom%patch%elts(p_par+1)%level
             if (l_cur .lt. l) then
                cycle
             else
                if (l_cur .gt. l) exit
             end if
          end if

          call get_offs_Domain (dom, p_par, offs, dims)
          call routine (dom, p_par, c, offs, dims, zlev)
       end do
    end do
  end subroutine apply_to_penta_d

  subroutine apply_to_penta (routine, l, zlev)
    external :: routine
    integer  :: l, zlev
    
    integer                        :: c, d, l_cur, p, p_par
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do d = 1, size(grid)
       call apply_to_penta_d(routine, grid(d), l, zlev)
    end do
  end subroutine apply_to_penta

  function tri_idx (i, j, tri, offs, dims)
    integer                        :: tri_idx
    integer                        :: i, j
    integer, dimension(3)          :: tri
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    tri_idx = TRIAG * idx(i + tri(1), j + tri(2), offs, dims) + tri(3)
  end function tri_idx

  function nidx (i, j, s, offs, dims)
    integer                        :: nidx
    integer                        :: i, j, s
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    nidx = offs(s+1) + j*dims(1,s+1) + i
  end function nidx

  function idx2 (i, j, noffs, offs, dims)
    integer                        :: idx2
    integer                        :: i, j
    integer, dimension(2)          :: noffs
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    idx2 = idx(i + noffs(1), j + noffs(2), offs, dims)
  end function idx2

  function is_penta (dom, p, s)
    logical      :: is_penta
    type(Domain) :: dom
    integer      :: p, s

    integer :: n, side
    logical :: penta

    penta = .False.

    n = dom%patch%elts(p+1)%neigh(s+1)

    if (n .lt. 0) then
       n = -n
       side = dom%bdry_patch%elts(n+1)%side

       if (side .gt. 0) then
          is_penta = dom%penta(side)
          return
       end if
    end if

    is_penta = penta
  end function is_penta

  subroutine apply (routine, zlev)
    external :: routine
    integer  :: zlev

    integer :: l

    do l = level_start, level_end
       call apply_onescale(routine, l, zlev, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
  end subroutine apply

  subroutine apply_onescale (routine, l, zlev, st, en)
    external :: routine
    integer  :: en, l, st, zlev

    integer :: d, k

    do d = 1, size(grid)
       do k = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch(routine, grid(d), grid(d)%lev(l)%elts(k), zlev, st, en)
       end do
    end do
  end subroutine apply_onescale

  subroutine apply_onescale__int (routine, l, zlev, st, en, ival)
    external :: routine
    integer  :: en, ival, l, st, zlev

    integer               :: d, k
    logical, dimension(2) :: pole_done

    do d = 1, size(grid)
       do k = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch__int(routine, grid(d), grid(d)%lev(l)%elts(k), zlev, st, en, ival)
       end do
    end do
  end subroutine apply_onescale__int

  function idx__fast (i, j, offs)
    integer ::  idx__fast
    integer :: i, j, offs

    idx__fast = PATCH_SIZE*j + i + offs
  end function idx__fast

  function idx (i0, j0, offs, dims)
    ! given regular array coordinates (i0,j0), offset array offs and domain array dims returns associated grid element as elts(idx+1)
    integer                        :: idx
    integer                        :: i0, j0
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: i, j

    i = i0
    j = j0
    if (i .lt. 0) then
       if (j .lt. 0) then
          idx = offs(IJMINUS+1) + (j0 + dims(2,IJMINUS+1))*dims(1,IJMINUS+1) + i + dims(1,IJMINUS+1)
          return
       else
          if (j .ge. PATCH_SIZE) then
             idx = offs(IMINUSJPLUS+1) + (-PATCH_SIZE + j)*dims(1,IMINUSJPLUS+1) + i + dims(1,IMINUSJPLUS+1)
             return
          else
             idx = offs(WEST+1) + j*dims(1,WEST+1) + i + dims(1,IMINUS+1)
             return
          end if
       end if
    else
       if (i .ge. PATCH_SIZE) then
          i = i - PATCH_SIZE
          if (j .ge. PATCH_SIZE) then
             j = j - PATCH_SIZE
             idx = offs(NORTHEAST+1) + j*dims(1,NORTHEAST+1) + i
             return
          else
             if (j .lt. 0) then
                idx = offs(SOUTHEAST+1) + (j0 + dims(2,SOUTHEAST+1))*dims(1,SOUTHEAST+1) + (i0 - PATCH_SIZE) 
                return
             else
                idx = offs(EAST+1) + j*dims(1,EAST+1) + i
                return
             end if
          end if
       else
          if (j .lt. 0) then
             idx = offs(JMINUS+1) + (j + dims(2,JMINUS+1))*dims(1,JMINUS+1) + i
             return
          else
             if (j .ge. PATCH_SIZE) then
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

  function add_bdry_patch_Domain (self, side)
    integer      :: add_bdry_patch_Domain
    type(Domain) :: self
    integer      :: side
    
    integer :: p

    p = self%bdry_patch%length

    call append(self%bdry_patch, Bdry_Patch(self%node%length, side, 0))

    call extend_Domain(self, BDRY_THICKNESS*PATCH_SIZE)

    add_bdry_patch_Domain = p
  end function add_bdry_patch_Domain

  subroutine get_offs_Domain5 (self, p, offs, dims, inner_patch)
    type(Domain)                    :: self
    integer                         :: p
    integer, dimension(N_BDRY+1)    :: offs
    integer, dimension(2,N_BDRY+1)  :: dims
    logical, dimension(4), optional :: inner_patch(4)

    integer :: i, n

    offs = -1
    dims = 0
    offs(1) = self%patch%elts(p+1)%elts_start

    if (present(inner_patch)) inner_patch = .False.

    do i = 1, N_BDRY
       n = self%patch%elts(p+1)%neigh(i)
       if (n .gt. 0) then
          offs(i+1) = self%patch%elts(n+1)%elts_start
          dims(:,i+1) = PATCH_SIZE
          if (present(inner_patch) .and. i .le. 4) inner_patch(i) = .True.
       else
          if (n .lt. 0) then
             offs(i+1) = self%bdry_patch%elts(-n+1)%elts_start
             dims(:,i+1) = sides_dims(:,abs(self%bdry_patch%elts(-n+1)%side) + 1)
          end if
       end if
    end do
  end subroutine get_offs_Domain5

  subroutine get_offs_Domain (self, p, offs, dims, inner_bdry)
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
       if (n .gt. 0) then
          offs(i+1) = self%patch%elts(n+1)%elts_start
          dims(:,i+1) = PATCH_SIZE
          if (present(inner_bdry) .and. i .le. 4) inner_bdry(i) = .True.
       else
          if (n .lt. 0) then
             offs(i+1) = self%bdry_patch%elts(-n+1)%elts_start
             dims(:,i+1) = sides_dims(:,abs(self%bdry_patch%elts(-n+1)%side) + 1)
             if (present(inner_bdry) .and. i .le. 4) &
                  inner_bdry(i) = self%bdry_patch%elts(-n+1)%side .lt. 0
          end if
       end if
    end do
  end subroutine get_offs_Domain

  subroutine set_neigh_Domain (self, s, id, rot)
    type(Domain) :: self
    integer      :: id, rot, s

    self%neigh(s) = id
    self%neigh_rot(s) = rot
  end subroutine set_neigh_Domain

  function find_neigh_bdry_patch_Domain (self, p_par, c, s)
    integer      :: find_neigh_bdry_patch_Domain
    type(Domain) :: self
    integer      :: p_par, c, s
    
    integer :: p_chd, p_par1, c1

    find_neigh_bdry_patch_Domain = 0

    if (p_par .le. 0) return

    p_chd = self%patch%elts(p_par+1)%children(c+1)
    if (p_chd .gt. 0) find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh(s+1)

    if (s .ge. 4) return

    if (find_neigh_bdry_patch_Domain .eq. 0) then
       call find_neigh_patch2_Domain (self, p_par, c, modulo(s+1,4), p_par1, c1)
       if (p_par1 .gt. 0) then
          p_chd = self%patch%elts(p_par1+1)%children(c1+1)
          if (p_chd .gt. 0) find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh((modulo(s-1,4)+4)+1)
       end if
    end if

    if (find_neigh_bdry_patch_Domain .eq. 0) then
       call find_neigh_patch2_Domain (self, p_par, c, modulo(s-1,4), p_par1, c1)
       if (p_par1 .gt. 0) then
          p_chd = self%patch%elts(p_par1+1)%children(c1+1)
          if (p_chd .gt. 0) find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh((s+4)+1)
       end if
    end if
  end function find_neigh_bdry_patch_Domain

  subroutine find_neigh_patch2_Domain (self, p_par0, c0, s0, p_par, c)
    ! for patch given as `c0`-th child of `p_par0`
    ! find neighbour with respect to side `s0`
    ! result as `c`-th child of patch `p_par`
    type(Domain)         :: self
    integer, intent(in)  :: p_par0, c0, s0
    integer, intent(out) :: p_par, c
    
    integer :: s_par

    s_par = par_side(c0, s0)
    if (s0 .eq. c0 .or. s0 .eq. modulo(c0 + 1, 4) .or. (s0 .eq. c0 + 4)) then
       p_par = self%patch%elts(p_par0+1)%neigh(s0+1)
    else if (s0 .eq. modulo(c0 + 1, 4) + 4 .or. s0 .eq. modulo(c0 - 1, 4) + 4) then
       p_par = self%patch%elts(p_par0+1)%neigh(s_par+1)
    else ! neighbour patch on same parent
       p_par = p_par0
    end if

    c = ngb_chd_idx(c0, s0)
  contains
     function ngb_chd_idx(c, s)
      integer :: ngb_chd_idx
      integer :: c, s

      if (s .lt. 4) then
         ngb_chd_idx = modulo(-c + 2*s + 1, 4)
         return
      else
         ngb_chd_idx = modulo(c + 2, 4)
         return
      end if
    end function ngb_chd_idx
  end subroutine find_neigh_patch2_Domain

  function find_neigh_patch_Domain (self, p_par, c, s_chd)
    !  finds the neighbour at side s of c-th child
    !             of parent patch p_par through p_par
    integer      :: find_neigh_patch_Domain
    type(Domain) :: self
    integer      :: c, p_par, s_chd

    integer :: c1, n_par, n_side, p_chd, p_par1, s_help, s_par, typ

    call find_neigh_patch2_Domain (self, p_par, c, s_chd, n_par, c1)

    s_par = par_side(c, s_chd)

    if (n_par .gt. 0) then
       find_neigh_patch_Domain = self%patch%elts(n_par+1)%children(c1+1)
    else if (n_par .lt. 0) then ! bdry patch
       n_par = -n_par
       typ = self%bdry_patch%elts(n_par+1)%side
       if (typ .lt. 1) then
          find_neigh_patch_Domain = 0
          return
       end if

       if (s_chd + 1 .eq. typ) then ! side
          find_neigh_patch_Domain = find_neigh_bdry_patch_Domain(self, p_par, c, s_chd)
       else ! patch corner, but domain side
          if (s_chd + 1 .eq. typ + 4) then ! s_chd after typ (clockwise sense)
             s_help = modulo(typ, 4)
          else if (s_chd + 1 .eq. modulo(typ-2, 4) + 5) then ! s_chd before typ
             s_help = modulo(typ-2, 4)
          end if
          call find_neigh_patch2_Domain(self, p_par, c, s_help, p_par1, c1)
          find_neigh_patch_Domain = find_neigh_bdry_patch_Domain(self, p_par1, c1, typ-1)
       end if
    end if
  end function find_neigh_patch_Domain

  function par_side (c, s)
    integer :: par_side
    integer :: c, s

    if (s .eq. modulo(c + 1, 4) + 4) then
       par_side = modulo(c + 1, 4)
       return
    else
       if (s .eq. modulo(c - 1, 4) + 4) then
          par_side = c
          return
       else
          par_side = s
          return
       end if
    end if
  end function par_side

  subroutine extend_Domain (self, num)
    type(Domain) :: self
    integer      :: num
    
    integer :: d, k, v

    d = self%id + 1

    call extend(self%node, num, ORIGIN)

    do k = 1, zlevels
       call extend(sol(S_MASS,k)%data(d), num, 1.0_8) ! set 1.0 so PV computation does not raise float pt exception if undefined
       call extend(sol(S_TEMP,k)%data(d), num, 0.0_8)
       call extend(sol(S_VELO,k)%data(d), EDGE*num, 0.0_8)
    end do
  end subroutine extend_Domain

  subroutine init_Domain (self)
    type(Domain) :: self
    
    integer :: i, k, l, r

    call init(self%patch, 1)
    call init(self%bdry_patch, 1)
    call init(self%node, 1)

    allocate(self%src_patch(n_process,min_level:max_level))

    do l = min_level, max_level
       do r = 1, n_process
          call init(self%src_patch(r,l), 0)
       end do
    end do

    allocate(self%lev(min_level-1:max_level))

    do i = lbound(self%lev,1), ubound(self%lev,1)
       call init(self%lev(i), 0)
    end do

    do i = 1, N_GLO_DOMAIN
       call init(self%send_conn(i), 0)
    end do

    call init(self%send_pa_all, 0)

    do i = 1, N_GLO_DOMAIN
       call init(self%recv_pa(i), 0)
    end do

    do k = AT_NODE, AT_EDGE
       do i = 1, N_GLO_DOMAIN
          call init(self%pack(k,i), 0)
          call init(self%unpk(k,i), 0)
       end do
    end do

    self%pole_master = .False.
  end subroutine init_Domain

  function add_patch_Domain (self, level)
    integer      :: add_patch_Domain
    type(Domain) :: self
    
    integer :: level, p

    p = self%patch%length

    call append(self%lev(level), p)
    call append(self%patch, Patch(self%node%length, level, 0, 0, 0, .False.))
    call extend_Domain(self, PATCH_SIZE**2)

    add_patch_Domain = p
  end function add_patch_Domain

end module domain_mod
