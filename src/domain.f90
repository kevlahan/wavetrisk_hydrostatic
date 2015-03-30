module domain_mod
  use param_mod
  use shared_mod
  use geom_mod
  use patch_mod
  use dyn_arrays
  use arch_mod
  implicit none
  integer, dimension(2,9) :: sides_dims
  integer, dimension(2,4) :: chd_offs
  type Domain
      type(Patch_Array) patch
      type(Bdry_Patch_Array) bdry_patch
      type(Coord_Array) node
      type(Int_Array), allocatable :: lev(:)
      type(Int_Array), dimension(N_GLO_DOMAIN) :: send_conn
      type(Int_Array) send_pa_all
      type(Int_Array), dimension(N_GLO_DOMAIN) :: recv_pa
      type(Int_Array), dimension(AT_NODE:AT_EDGE,N_GLO_DOMAIN) :: pack
      type(Int_Array), dimension(AT_NODE:AT_EDGE,N_GLO_DOMAIN) :: unpk
      type(Int_Array), allocatable :: src_patch(:,:)
      logical, dimension(N_BDRY) :: penta
      integer id
      integer, dimension(N_BDRY) :: neigh
      integer, dimension(2) :: neigh_over_pole
      type(Int_Array) neigh_pa_over_pole
      integer, dimension(N_BDRY) :: neigh_rot
      type(Coord_Array) ccentre
      type(Coord_Array) midpt
      type(Float_Array) pedlen
      type(Float_Array) len
      type(Float_Array) corolis
      type(Float_Array) topo
      type(Areas_Array) areas
      type(Float_Array) triarea
      type(Float_Array) bernoulli
      type(Float_Array) kin_energy
      type(Float_Array) windstress
      type(Float_Array) qe
      type(Float_Array) divu
      type(Float_Array) vort
      type(Overl_Area_Array) overl_areas
      type(Iu_Wgt_Array) I_u_wgt
      type(RF_Wgt_Array) R_F_wgt
      type(Int_Array) mask_p
      type(Int_Array) mask_u
      type(Int_Array) level
      logical pole_master(2)
  end type
  type Float_Field
      type(Float_Array), allocatable :: data(:)
      logical bdry_uptodate
      integer pos
  end type
  type(Domain), allocatable :: grid(:)
  type(Float_Field), TARGET :: sol(S_HEIGHT:S_VELO), trend(S_HEIGHT:S_VELO), thickflux
  real(8), pointer :: height(:), velo(:), dheight(:), dvelo(:), wc_u(:), wc_h(:), tflux(:)
  ! for penalization boundary condition
  type(Float_Field), target :: penal
  real(8), pointer :: chi(:), chiflux(:)

contains
  subroutine init_Float_Field(self, pos)
      type(Float_Field) self
      integer pos
      self%bdry_uptodate = .False.
      self%pos = pos
      allocate(self%data(n_domain(rank+1)))
  end subroutine

  subroutine init_domain_mod()
      integer i, v
      logical :: initialized = .False.
      if (initialized) return ! initialize only once
      call init_shared_mod()
      call init_sphere_mod()
      call init_patch_mod()
      call init_arch_mod()
      sides_dims = reshape((/PATCH_SIZE, PATCH_SIZE, PATCH_SIZE, &
              BDRY_THICKNESS, BDRY_THICKNESS, PATCH_SIZE, PATCH_SIZE, &
              BDRY_THICKNESS, BDRY_THICKNESS, PATCH_SIZE, BDRY_THICKNESS, &
              BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS, &
              BDRY_THICKNESS, BDRY_THICKNESS, BDRY_THICKNESS/), (/2, 9/))
      chd_offs = reshape((/PATCH_SIZE/2, PATCH_SIZE/2, PATCH_SIZE/2, 0, 0, 0, &
              0, PATCH_SIZE/2/), (/2, 4/))
      initialized = .True.
  end subroutine init_domain_mod

  subroutine apply_onescale_to_patch__int(routine, dom, p, st, en, ival)
      external routine
      type(Domain) dom
      integer p
      integer st, en
      integer ival
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer j
      integer i
      call get_offs_Domain(dom, p, offs, dims)
      do j = st + 1, PATCH_SIZE + en
          do i = st + 1, PATCH_SIZE + en
              call routine(dom, p, i - 1, j - 1, offs, dims, ival)
          end do
      end do
  end subroutine

  subroutine apply_onescale_to_patch(routine, dom, p, st, en)
      external routine
      type(Domain) dom
      integer p
      integer st
      integer en
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer bdry(JPlUS:IMINUS)
      logical inner_bdry(JPlUS:IMINUS)
      integer j
      integer i
      call get_offs_Domain(dom, p, offs, dims, inner_bdry)
      bdry = (/en, en, st, st/)
      where (inner_bdry) bdry = 0
      do j = bdry(JMINUS) + 1, PATCH_SIZE + bdry(JPLUS)
          do i = bdry(IMINUS) + 1, PATCH_SIZE + bdry(IPLUS)
              call routine(dom, i - 1, j - 1, offs, dims)
          end do
      end do
  end subroutine

  integer function idx__fast(i, j, offs)
      integer i
      integer j
      integer offs
      idx__fast = PATCH_SIZE*j + i + offs
      return
  end function

  integer function idx(i0, j0, offs, dims)
      integer i0
      integer j0
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer i
      integer j
      i = i0
      j = j0
      if (i .lt. 0) then
          if (j .lt. 0) then
              idx = offs(IJMINUS+1) + (j0 + &
                      dims(2,IJMINUS+1))*dims(1,IJMINUS+1) + i + &
                      dims(1,IJMINUS+1)
              return
          else
              if (j .ge. PATCH_SIZE) then
                  idx = offs(IMINUSJPLUS+1) + (-PATCH_SIZE + &
                          j)*dims(1,IMINUSJPLUS+1) + i + dims(1,IMINUSJPLUS+1)
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
                      idx = offs(SOUTHEAST+1) + (j0 + dims(2,SOUTHEAST+1))*dims(1,SOUTHEAST+1) &
                                              + (i0 - PATCH_SIZE) 
                      return
                  else
                      idx = offs(EAST+1) + j*dims(1,EAST+1) + i
                      return
                  end if
              end if
          else
              if (j .lt. 0) then
                  idx = offs(JMINUS+1) + (j + &
                          dims(2,JMINUS+1))*dims(1,JMINUS+1) + i
                  return
              else
                  if (j .ge. PATCH_SIZE) then
                      idx = offs(NORTH+1) + (-PATCH_SIZE + j)*dims(1,NORTH+1) + &
                              i
                      return
                  else
                      idx = offs(1) + PATCH_SIZE*j + i
                      return
                  end if
              end if
          end if
      end if
  end function

  subroutine apply_interscale_to_patch2(routine, dom, p_par, st, en)
      external routine
      type(Domain) dom
      integer p_par
      integer st
      integer en
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,N_BDRY + 1) :: dims_par
      integer c
      integer p_chd
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,N_BDRY + 1) :: dims_chd
      integer j
      integer j_chd
      integer j_par
      integer i
      integer i_chd
      integer i_par
      integer bdry(JPlUS:IMINUS)
      logical inner_bdry(JPlUS:IMINUS)
      !TODO{uncomment & test}  if (dom%patch%elts(p_par+1)%active .eq. NONE) return
      call get_offs_Domain(dom, p_par, offs_par, dims_par)
      do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .eq. 0) then
              cycle
          end if
          call get_offs_Domain(dom, p_chd, offs_chd, dims_chd, inner_bdry)
          bdry = (/en, en, st, st/)
          where (inner_bdry) bdry = 0
          do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
              j_chd = (j - 1)*2
              j_par = j - 1 + chd_offs(2,c)
              do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
                  i_chd = (i - 1)*2
                  i_par = i - 1 + chd_offs(1,c)
                  call routine(dom, p_chd, i_par, j_par, i_chd, j_chd, &
                          offs_par, dims_par, offs_chd, dims_chd)
              end do
          end do
      end do
  end subroutine

  integer function ed_idx(i, j, ed, offs, dims)
      integer i
      integer j
      integer, dimension(3) :: ed
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      ed_idx = idx(i + ed(1), j + ed(2), offs, dims)*EDGE + ed(3)
  end function

  subroutine apply_onescale_d(routine, dom, l, st, en)
      external routine
      type(Domain) dom
      integer l
      integer st
      integer en
      integer k
      do k = 1, dom%lev(l)%length
          call apply_onescale_to_patch(routine, dom, &
                  dom%lev(l)%elts(k), st, en)
      end do
  end subroutine

  subroutine apply_interscale_d(routine, dom, l, st, en)
      external routine
      type(Domain) dom
      integer l
      integer st
      integer en
      integer k
      do k = 1, dom%lev(l)%length
          call apply_interscale_to_patch(routine, dom, &
                  dom%lev(l)%elts(k), st, en)
      end do
  end subroutine

  subroutine apply_interscale_d2(routine, dom, l, st, en)
      external routine
      type(Domain) dom
      integer l
      integer st
      integer en
      integer k
      do k = 1, dom%lev(l)%length
          call apply_interscale_to_patch22(routine, dom, &
                  dom%lev(l)%elts(k), st, en)
      end do
  end subroutine

  subroutine apply_to_pole_d(routine, dom, l, ival, to_all)
      external routine
      type(Domain) dom
      integer l
      integer ival
      logical to_all
      integer c
      integer p
      integer p_par
      integer l_cur
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims

      do c = SOUTHEAST, NORTHWEST, 2
          if (.not. dom%pole_master(c/2-2) .and. .not. to_all) cycle
          if ( .not. dom%penta(c)) then
              cycle
          end if
          p = 1
          if (.not. dom%neigh(c) .eq. POLE) cycle
          do while (p .gt. 0)
              p_par = p
              p = dom%patch%elts(p_par+1)%children(c-4)
              if ( .not. l .eq. NONE) then
                  l_cur = dom%patch%elts(p_par+1)%level
                  if (l_cur .lt. l) then
                      cycle
                  else
                      if (l_cur .gt. l) then
                          exit
                      end if
                  end if
              end if
              call get_offs_Domain(dom, p_par, offs, dims)
              if (c .eq. NORTHWEST) then
                  call routine(dom, p_par, 0, PATCH_SIZE, offs, dims, ival) ! NORTHPOLE
              else
                  call routine(dom, p_par, PATCH_SIZE, 0, offs, dims, ival) ! SOUTHPOLE
              end if
          end do
      end do
  end subroutine

  subroutine apply_to_pole(routine, l, ival, to_all)
      external routine
      integer l
      integer ival
      logical to_all
      integer d
      integer c
      integer p
      integer p_par
      integer l_cur
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      do d = 1, size(grid)
          call apply_to_pole_d(routine, grid(d), l, ival, to_all)
      end do
  end subroutine

  subroutine apply_to_penta_d(routine, dom, l)
      external routine
      type(Domain) dom
      integer l
      integer d
      integer c
      integer p
      integer p_par
      integer l_cur
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      do c = NORTHEAST, NORTHWEST
          p = 1
          if ( .not. dom%penta(c)) then
              cycle
          end if
          do while (p .gt. 0)
              p_par = p
              p = dom%patch%elts(p_par+1)%children(c-4)
              if (.not. l .eq. NONE) then
                  l_cur = dom%patch%elts(p_par+1)%level
                  if (l_cur .lt. l) then
                      cycle
                  else
                      if (l_cur .gt. l) then
                          exit
                      end if
                  end if
              end if
              call get_offs_Domain(dom, p_par, offs, dims)
              call routine(dom, p_par, c, offs, dims)
          end do
      end do
  end subroutine

  subroutine apply_to_penta(routine, l)
      external routine
      integer l
      integer d
      integer c
      integer p
      integer p_par
      integer l_cur
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      do d = 1, size(grid)
          call apply_to_penta_d(routine, grid(d), l)
      end do
  end subroutine

  integer function tri_idx(i, j, tri, offs, dims)
      integer i
      integer j
      integer, dimension(3) :: tri
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      tri_idx = idx(i + tri(1), j + tri(2), offs, dims)*TRIAG + tri(3)
      return
  end function

  integer function nidx(i, j, s, offs, dims)
      integer i
      integer j
      integer s
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      nidx = offs(s+1) + j*dims(1,s+1) + i
      return
  end function

  subroutine apply_interscale_to_patch3(routine, dom, p_par, c, st, en)
      external routine
      type(Domain) dom
      integer p_par
      integer st
      integer en
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,N_BDRY + 1) :: dims_par
      integer c
      integer p_chd
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,N_BDRY + 1) :: dims_chd
      integer j
      integer j_chd
      integer j_par
      integer i
      integer i_chd
      integer i_par
      integer bdry(JPlUS:IMINUS)
      logical inner_bdry(JPlUS:IMINUS)
      !TODO{uncomment & test}  if (dom%patch%elts(p_par+1)%active .eq. NONE) return
      call get_offs_Domain(dom, p_par, offs_par, dims_par)
      p_chd = dom%patch%elts(p_par+1)%children(c)
      call get_offs_Domain(dom, p_chd, offs_chd, dims_chd, inner_bdry)
      bdry = (/en, en, st, st/)
      do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
          j_chd = (j - 1)*2
          j_par = j - 1 + chd_offs(2,c)
          do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
              i_chd = (i - 1)*2
              i_par = i - 1 + chd_offs(1,c)
              call routine(dom, p_chd, i_par, j_par, i_chd, j_chd, offs_par, &
                      dims_par, offs_chd, dims_chd)
          end do
      end do
  end subroutine

  subroutine apply_interscale_to_patch22(routine, dom, p_par, st, en)
      external routine
      type(Domain) dom
      integer p_par
      integer st
      integer en
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,N_BDRY + 1) :: dims_par
      integer c
      integer p_chd
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,N_BDRY + 1) :: dims_chd
      integer j
      integer j_chd
      integer j_par
      integer i
      integer i_chd
      integer i_par
      integer bdry(JPlUS:IMINUS)
      logical inner_bdry(JPlUS:IMINUS)
      call get_offs_Domain(dom, p_par, offs_par, dims_par)
      do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .eq. 0) then
              cycle
          end if
          call get_offs_Domain(dom, p_chd, offs_chd, dims_chd)
          bdry = (/en, en, st, st/)
          do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
              j_chd = (j - 1)*2
              j_par = j - 1 + chd_offs(2,c)
              do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
                  i_chd = (i - 1)*2
                  i_par = i - 1 + chd_offs(1,c)
                  call routine(dom, i_par, j_par, i_chd, j_chd, offs_par, &
                          dims_par, offs_chd, dims_chd)
              end do
          end do
      end do
  end subroutine

  subroutine apply_interscale_to_patch(routine, dom, p_par, st, en)
      external routine
      type(Domain) dom
      integer p_par
      integer st
      integer en
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,N_BDRY + 1) :: dims_par
      integer c
      integer p_chd
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,N_BDRY + 1) :: dims_chd
      integer j
      integer j_chd
      integer j_par
      integer i
      integer i_chd
      integer i_par
      integer bdry(JPlUS:IMINUS)
      logical inner_bdry(JPlUS:IMINUS)
      call get_offs_Domain(dom, p_par, offs_par, dims_par)
      do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .eq. 0) then
              cycle
          end if
          call get_offs_Domain(dom, p_chd, offs_chd, dims_chd, inner_bdry)
          bdry = (/en, en, st, st/)
          where (inner_bdry) bdry = 0
          do j = bdry(JMINUS) + 1, PATCH_SIZE/2 + bdry(JPLUS)
              j_chd = (j - 1)*2
              j_par = j - 1 + chd_offs(2,c)
              do i = bdry(IMINUS) + 1, PATCH_SIZE/2 + bdry(IPLUS)
                  i_chd = (i - 1)*2
                  i_par = i - 1 + chd_offs(1,c)
                  call routine(dom, i_par, j_par, i_chd, j_chd, offs_par, &
                          dims_par, offs_chd, dims_chd)
              end do
          end do
      end do
  end subroutine

  subroutine apply_onescale2(routine, l, st, en)
      external routine
      integer l
      integer st
      integer en
      integer d
      integer k
      do d = 1, size(grid)
          do k = 1, grid(d)%lev(l)%length
              call apply_onescale_to_patch2(routine, grid(d), &
                      grid(d)%lev(l)%elts(k), st, en)
          end do
      end do
  end subroutine

  subroutine apply_onescale_to_patch5(routine, dom, p, st, en)
      external routine
      type(Domain) dom
      integer p
      integer st
      integer en
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      logical inner_patch(JPlUS:IMINUS)
      integer bdry(JPlUS:IMINUS)
      integer j
      integer i
      call get_offs_Domain5(dom, p, offs, dims, inner_patch)
      bdry = (/en, en, st, st/)
      where (inner_patch) bdry = 0
      do j = bdry(JMINUS) + 1, PATCH_SIZE + bdry(JPLUS)
          do i = bdry(IMINUS) + 1, PATCH_SIZE + bdry(IPLUS)
              call routine(dom, p, i - 1, j - 1, offs, dims)
          end do
      end do
  end subroutine

  subroutine apply_onescale_to_patch2(routine, dom, p, st, en)
      external routine
      type(Domain) dom
      integer p
      integer st
      integer en
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer j
      integer i
      call get_offs_Domain(dom, p, offs, dims)
      do j = st + 1, PATCH_SIZE + en
          do i = st + 1, PATCH_SIZE + en
              call routine(dom, p, i - 1, j - 1, offs, dims)
          end do
      end do
  end subroutine

  integer function idx2(i, j, noffs, offs, dims)
      integer i
      integer j
      integer, dimension(2) :: noffs
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      idx2 = idx(i + noffs(1), j + noffs(2), offs, dims)
  end function

  logical function is_penta(dom, p, s)
      type(Domain) dom
      integer p
      integer s
      logical penta
      integer n
      integer side
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
      return
  end function

  subroutine apply_interscale(routine, l, st, en)
      external routine
      integer l
      integer st
      integer en
      integer d
      do d = 1, size(grid)
          call apply_interscale_d(routine, grid(d), l, st, en)
      end do
  end subroutine

  subroutine apply(routine)
      external routine
      integer l
      do l = level_start, level_end
          call apply_onescale(routine, l, -BDRY_THICKNESS, BDRY_THICKNESS)
      end do
  end subroutine

  subroutine apply_onescale(routine, l, st, en)
      external routine
      integer l
      integer st
      integer en
      integer d
      integer k
      do d = 1, size(grid)
          do k = 1, grid(d)%lev(l)%length
              call apply_onescale_to_patch(routine, grid(d), &
                      grid(d)%lev(l)%elts(k), st, en)
          end do
      end do
  end subroutine

  subroutine apply_onescale__int(routine, l, st, en, ival)
      external routine
      integer l
      integer st
      integer en
      integer d
      integer k
      integer ival
      logical pole_done(2)
      do d = 1, size(grid)
          do k = 1, grid(d)%lev(l)%length
              call apply_onescale_to_patch__int(routine, grid(d), &
                      grid(d)%lev(l)%elts(k), st, en, ival)
          end do
      end do
  end subroutine

  integer function add_bdry_patch_Domain(self, side)
      type(Domain) self
      integer side
      integer p
      p = self%bdry_patch%length
      call append(self%bdry_patch, Bdry_Patch(self%node%length, side, 0))
      call extend_Domain(self, BDRY_THICKNESS*PATCH_SIZE)
      add_bdry_patch_Domain = p
  end function

  subroutine get_offs_Domain5(self, p, offs, dims, inner_patch)
      type(Domain) self
      integer p
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      logical, optional :: inner_patch(4)
      integer i
      integer n
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
  end subroutine

  subroutine get_offs_Domain(self, p, offs, dims, inner_bdry)
      type(Domain) self
      integer p
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      logical, optional :: inner_bdry(4)
      integer i
      integer n
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
  end subroutine

  subroutine set_neigh_Domain(self, s, id, rot)
      type(Domain) self
      integer s
      integer id
      integer rot
      self%neigh(s) = id
      self%neigh_rot(s) = rot
  end subroutine

  integer function find_neigh_bdry_patch_Domain(self, p_par, c, s)
      type(Domain) self
      integer p_par, c, s
      integer p_chd, p_par1, c1
      find_neigh_bdry_patch_Domain = 0
      if (p_par .le. 0) then
          return
      end if
      p_chd = self%patch%elts(p_par+1)%children(c+1)
      if (p_chd .gt. 0) &
          find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh(s+1)
      if (s .ge. 4) return
      if (find_neigh_bdry_patch_Domain .eq. 0) then
          call find_neigh_patch2_Domain(self, p_par, c, modulo(s+1,4), p_par1, c1)
          if (p_par1 .gt. 0) then
              p_chd = self%patch%elts(p_par1+1)%children(c1+1)
              if (p_chd .gt. 0) &
                  find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh((modulo(s-1,4)+4)+1)
          end if
      end if
      if (find_neigh_bdry_patch_Domain .eq. 0) then
          call find_neigh_patch2_Domain(self, p_par, c, modulo(s-1,4), p_par1, c1)
          if (p_par1 .gt. 0) then
              p_chd = self%patch%elts(p_par1+1)%children(c1+1)
              if (p_chd .gt. 0) &
                  find_neigh_bdry_patch_Domain = self%patch%elts(p_chd+1)%neigh((s+4)+1)
          end if
      end if
  end function

  subroutine find_neigh_patch2_Domain(self, p_par0, c0, s0, p_par, c)
      ! for patch given as `c0`-th child of `p_par0`
      ! find neighbour with respect to side `s0`
      ! result as `c`-th child of patch `p_par`
      type(Domain) self
      integer, intent(in) :: p_par0, c0, s0
      integer, intent(out) :: p_par, c
      integer s_par
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
      integer function ngb_chd_idx(c, s)
          integer c
          integer s
          if (s .lt. 4) then
              ngb_chd_idx = modulo(-c + 2*s + 1, 4)
              return
          else
              ngb_chd_idx = modulo(c + 2, 4)
              return
          end if
      end function
  end subroutine

  integer function find_neigh_patch_Domain(self, p_par, c, s_chd)
      type(Domain) self
      integer p_par, p_par1
      integer c, c1
      integer s_chd
      integer s_par
      integer n_par
      integer typ
      integer p_chd
      integer s_help
      integer n_side
      !  finds the neighbour at side s of c-th child
      !             of parent patch p_par through p_par

      call find_neigh_patch2_Domain(self, p_par, c, s_chd, n_par, c1)
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
  end function

  integer function par_side(c, s)
      integer c
      integer s
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
  end function

  subroutine extend_Domain(self, num)
      type(Domain) self
      integer num, d
      d = self%id + 1
      call extend(self%node, num, ORIGIN)
      call extend(sol(S_HEIGHT)%data(d), num, 1.0_8) ! set 1.0 so PV cpt does not raise float pt exception if undef
      call extend(sol(S_VELO)%data(d), EDGE*num, 0.0_8)
  end subroutine

  subroutine init_Domain(self)
      type(Domain) self
      integer i, k, l, r
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
  end subroutine

  integer function add_patch_Domain(self, level)
      type(Domain) self
      integer level
      integer p
      p = self%patch%length
      call append(self%lev(level), p)
      call append(self%patch, Patch(self%node%length, level, 0, 0, 0, .False.))
      call extend_Domain(self, PATCH_SIZE**2)
      add_patch_Domain = p
      return
  end function
end module domain_mod
