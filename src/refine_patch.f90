module refine_patch_mod
  use domain_mod
  use init_mod
  use wavelet_mod
  use mask_mod
  implicit none
  integer, parameter :: DOF_PER_PATCH = PATCH_SIZE * PATCH_SIZE * (EDGE + 1)
  integer, parameter :: FILLED_AND_FROZEN = DOF_PER_PATCH + 1
  
  logical :: max_level_exceeded
contains
  logical function refine ()
    ! Determines where new patches are needed when grid is refineds
    implicit none
    integer :: c, d, did_refine, old_n_patch, p_chd, p_par
    logical :: required
    
    ! Use threshold mask to call refine patch where necessary
    did_refine = FALSE
    do d = 1, size(grid)
       old_n_patch = grid(d)%patch%length
       do p_par = 3, grid(d)%patch%length
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par)%children(c)
             
             required = check_child_required (grid(d), p_par-1, c-1)   ! is child (finer grid) required?
             
             if (required .and. p_chd <= 0) then                       ! new patch required that does not yet exist
                if (grid(d)%patch%elts(p_par)%level == max_level) then ! cannot refine further
                   max_level_exceeded = .true.
                else
                   call refine_patch1 (grid(d), p_par-1, c-1)
                   did_refine = TRUE
                end if
             end if
          end do
       end do
       
       do p_par = 2, old_n_patch
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par)%children(c)
             if (p_chd+1 > old_n_patch) call refine_patch2 (grid(d), p_par-1, c-1)
          end do
       end do
       
    end do
    refine = sync_max_int (did_refine) == TRUE
    return
  end function refine

  subroutine post_refine
    implicit none
    integer :: d, p

    call comm_masks_mpi (NONE)

    level_end = sync_max_int (level_end)

    do d = 1, n_domain(rank+1)
       do p = 3, grid(d)%patch%length
          call connect_children (grid(d), p-1)
       end do
    end do

    call comm_patch_conn_mpi

    do d = 1, size(grid)
       call update_comm (grid(d))
    end do

    call comm_communication_mpi
    call comm_nodes9_mpi (get_areas, set_areas, NONE)
    call apply_to_penta (area_post_comm, NONE, z_null)
  end subroutine post_refine

  subroutine refine_patch1 (dom, p_par, c0)
    implicit none
    type(Domain) :: dom
    integer      :: p_par, c0
    
    integer                   :: c, d, i, i_chd, i_par, id_par, j, j_chd, j_par, k, lev, num, p_chd, s
    type(Coord), dimension(6) :: tmp
    
    !  Main difficulty: in order to precompute geometry, weights etc
    !         nodes outside the patch are needed that might not be part of the grid yet
    !         Better compute them temporarily than adding additional patches
    !         that are empty in terms of evaluating operators and slow down the inner loops
    !         Stategy A: use temporary bdry_patch

    lev = dom%patch%elts(p_par+1)%level
    
    if (lev       == max_level) return
    if (level_end == lev) level_end = level_end + 1

    num = dom%node%length 
    p_chd = add_patch_Domain (dom, lev+1)
    
    c = c0 + 1

    dom%patch%elts(p_par+1)%children(c) = p_chd
    num = dom%node%length - num

    call extend (dom%level, num, dom%patch%elts(p_chd+1)%level)
  end subroutine refine_patch1

  subroutine refine_patch2 (dom, p_par, c0)
    implicit none
    type(Domain) :: dom
    integer      :: p_par, c0
    
    type(Coord), dimension(6)      :: tmp
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd
    integer                        :: c, d, i, id_par, i_chd, i_par, j, j_chd, j_par, k, lev, num, p_chd, s, v

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    lev = dom%patch%elts(p_par+1)%level
    c = c0 + 1
    p_chd = dom%patch%elts(p_par+1)%children(c)
    do k = 1, 2
       call attach_bdry (dom, p_par, c-1, modulo(c + k - 2, N_CHDRN), side (dom, p_par, modulo(c + k - 2, N_CHDRN)))
       call attach_bdry (dom, p_par, c-1, modulo(c + k, N_CHDRN), -(modulo(c + k, N_CHDRN) + 1))
       call connect_cousin (dom, p_par, p_chd, modulo(c + k - 2, N_CHDRN), modulo(c + k - 2, N_CHDRN), modulo(c - 2*k + 2, N_CHDRN))
    end do

    call attach_bdry (dom, p_par, c-1, c + 3, side (dom, p_par, c + 3))
    call connect_cousin (dom, p_par, p_chd, c + 3, c + 3, modulo(c + 1, N_CHDRN))

    call attach_bdry (dom, p_par, c-1, modulo(c + 1, N_CHDRN) + 4, -(modulo(c + 1, N_CHDRN) + 4 + 1))
    call attach_bdry (dom, p_par, c-1, modulo(c,     N_CHDRN) + 4, side (dom, p_par, modulo(c, N_CHDRN)))
    call attach_bdry (dom, p_par, c-1, modulo(c + 2, N_CHDRN) + 4, side (dom, p_par, modulo(c-1, N_CHDRN)))

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

    do j = 0, PATCH_SIZE/2 + 1
       j_chd = (j - 1)*2
       j_par = j - 1 + chd_offs(2,c)
       do i = 0, PATCH_SIZE/2 + 1
          i_chd = (i - 1)*2
          i_par = i - 1 + chd_offs(1,c)
          id_par = idx(i_par, j_par, offs_par, dims_par)

          dom%node%elts(idx(i_chd, j_chd, offs_chd, dims_chd) + 1) = dom%node%elts(id_par+1)

          dom%node%elts(idx(i_chd + 1, j_chd,     offs_chd, dims_chd) + 1) = dom%midpt%elts(EDGE*id_par+RT+1)
          dom%node%elts(idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd) + 1) = dom%midpt%elts(EDGE*id_par+DG+1)
          dom%node%elts(idx(i_chd,     j_chd + 1, offs_chd, dims_chd) + 1) = dom%midpt%elts(EDGE*id_par+UP+1)
       end do
    end do

    if (is_penta (dom, p_chd, IPLUSJMINUS-1)) &
         dom%node%elts(idx(PATCH_SIZE, -1, offs_chd, dims_chd) + 1) = &
         mid_pt ( &
         dom%node%elts(idx(PATCH_SIZE + 1, 0, offs_par, dims_par) + 1), &
         dom%node%elts(idx(PATCH_SIZE,     0, offs_par, dims_par) + 1) )

    if (is_penta (dom, p_chd, IMINUSJPLUS-1)) &
         dom%node%elts(idx(-1, PATCH_SIZE, offs_chd, dims_chd) + 1) = &
         mid_pt ( &
         dom%node%elts(idx(0, PATCH_SIZE+1, offs_par, dims_par)+1), &
         dom%node%elts(idx(0, PATCH_SIZE,   offs_par, dims_par)+1) )

    num = dom%node%length - dom%areas%length
    d = dom%id + 1

    call extend (dom%ccentre, TRIAG * num, ORIGIN)
    call apply_onescale_to_patch2 (ccentre, dom, p_chd, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    
    call ccentre_penta (dom, p_chd)
    call extend (dom%midpt, EDGE * num, ORIGIN)
    
    call apply_onescale_to_patch2 (midpt, dom, p_chd, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)
    
    call extend (dom%pedlen, EDGE * num, 0d0)
    call extend (dom%len,    EDGE * num, 0d0)
    
    call apply_onescale_to_patch2 (lengths, dom, p_chd, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)

    tmp = ORIGIN
    
    call extend (dom%areas, num, Areas (0d0, 0d0))
    call apply_onescale_to_patch2 (cpt_areas, dom, p_chd, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)
    
    call extend (dom%triarea, EDGE * num, 1d0)
    call apply_onescale_to_patch (cpt_triarea, dom, p_chd, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)
    
    call extend (dom%coriolis, TRIAG * num, 0d0)
    call apply_onescale_to_patch (coriolis, dom, p_chd, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS)

    ! Initialize domain variables to zero
    call extend (dom%surf_press,    num, 0d0)
    call extend (dom%press,         num, 0d0)
    call extend (dom%geopot,        num, 0d0)
    call extend (dom%u_zonal,       num, 0d0)
    call extend (dom%v_merid,       num, 0d0)
    call extend (dom%press_lower,   num, 0d0)
    call extend (dom%geopot_lower,  num, 0d0)
    call extend (dom%bernoulli,     num, 0d0)
    call extend (dom%ke,            num, 0d0)
    call extend (dom%divu,          num, 0d0)
    
    call extend (dom%qe,      EDGE * num, 0d0)
    call extend (dom%vort,   TRIAG * num, 0d0)

    call extend (topography%data(d), num, 0d0)

    ! Initialize float fields and float arrays to zero
    if (sso) then
       do k = 1, 4
          call extend (sso_param(k)%data(d), num, 0d0)
       end do
    end if
    
    do k = zmin, zmax
       call extend (penal_node(k)%data(d),        num, 0d0)
       call extend (penal_edge(k)%data(d), EDGE * num, 0d0)
       call extend (exner_fun(k)%data(d),         num, 0d0)
       
       do v = scalars(1), scalars(2)
          if (k > 0) call extend (trend(v,k)%data(d), num, 0d0)
          call extend (wav_coeff(v,k)%data(d), num, 0d0)
       end do
       if (k > 0) call extend (trend(S_VELO,k)%data(d),     EDGE*num, 0d0)
       call extend (wav_coeff(S_VELO,k)%data(d), EDGE * num, 0d0)
    end do
    call extend (exner_fun(zmax+1)%data(d), num, 0d0)

    ! Initialize vertical diffusion variables to zero
    if (vert_diffuse) then
       call extend (Kt(0)%data(d), num, 0d0)
       call extend (Kv(0)%data(d), num, 0d0)
       do k = 1, zlevels
          call extend (Kt(k)%data(d),      num, 0d0)
          call extend (Kv(k)%data(d),      num, 0d0)
          call extend (tke(k)%data(d),     num, 0d0)
          call extend (wav_tke(k)%data(d), num, 0d0)
       end do
    end if

    ! Initialize Laplacian diffusion variables to zero
    call extend (Laplacian_vector(S_DIVU)%data(d),     num,  0d0)
    call extend (Laplacian_vector(S_ROTU)%data(d), EDGE * num, 0d0)
    do v = scalars(1), scalars(2)
       call extend (horiz_flux(v)%data(d),       EDGE * num, 0d0)
       call extend (Laplacian_scalar(v)%data(d),      num, 0d0)
    end do
    
    ! Initialize mask and wavelet variables to zero
    call extend (dom%overl_areas, EDGE * num, Overl_Area(0d0, 0d0))
    call extend (dom%I_u_wgt,     EDGE * num, Iu_Wgt (0d0))
    call extend (dom%R_F_wgt,            num, RF_Wgt (0d0))
    call extend (dom%mask_n,             num, ZERO)
    call extend (dom%mask_e,      EDGE * num, ZERO)
    
    call apply_interscale_to_patch3 (set_WT_wgts, dom, p_par, c, z_null, 0, 0)
    call apply_interscale_to_patch3 (set_RF_wgts, dom, p_par, c, z_null, 0, 0)

    ! Extend patch level variable
    num = dom%node%length - dom%level%length
    call extend (dom%level, num, dom%patch%elts(p_chd+1)%level)
  end subroutine refine_patch2

  subroutine add_second_level
    implicit none
    integer :: d, c

    do d = 1, size(grid)
       do c = 1, N_CHDRN
          call refine_patch1 (grid(d), 1, c-1)
       end do
       do c = 1, N_CHDRN
          call refine_patch2 (grid(d), 1, c-1)
       end do
       call connect_children (grid(d), 1)
    end do

    call comm_patch_conn_mpi

    do d = 1, size(grid)
       call update_comm (grid(d))
    end do

    call comm_communication_mpi
    call comm_nodes9_mpi (get_areas, set_areas, NONE)
    call apply_to_penta (area_post_comm, NONE, z_null)
  end subroutine add_second_level

  logical function check_child_required (dom, p_par, c)
    ! Determines where child (finer grid) is required based on mask value for parent
    ! child is required if parent node is in the adjacent zone or parent edge can be obtained by restriction
    implicit none
    type(Domain) :: dom
    integer      :: p_par, c

    integer                        :: e, j0, j_par, i0, i_par, id_e, id_par
    integer                        :: st, en
    integer, dimension(N_BDRY+1)   :: offs_par
    integer, dimension(2,N_BDRY+1) :: dims_par
    logical                        :: required

    st = -BDRY_THICKNESS
    en =  BDRY_THICKNESS

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do j0 = st + 1, PATCH_SIZE/2 + en
       j_par = j0 - 1 + chd_offs(2,c+1)
       do i0 = st + 1, PATCH_SIZE/2 + en
          i_par = i0 - 1 + chd_offs(1,c+1)
          
          id_par = idx (i_par, j_par, offs_par, dims_par) 

          required = dom%mask_n%elts(id_par+1) >= ADJZONE               ! child is required if parent is in adjacent zone
          do e = 1, EDGE
             id_e = EDGE*id_par + e
             required = required .or. dom%mask_e%elts(id_e) >= RESTRCT  ! child is required if parent edge can be obtained by restriction
          end do
          
          if (required) then
             check_child_required = .true.
             return
          end if
       end do
    end do
    check_child_required = .false.
  end function check_child_required

  function get_child_and_neigh_patches (dom, p_par, c)
    implicit none
    integer, dimension(4) :: get_child_and_neigh_patches
    type(Domain)          :: dom
    integer               :: p_par, c
    
    integer :: n

    get_child_and_neigh_patches = 0
    get_child_and_neigh_patches(1) = dom%patch%elts(p_par+1)%children(c)
    n = dom%patch%elts(p_par+1)%neigh(c) ! side
    if (n > 0) then
       get_child_and_neigh_patches(2) = dom%patch%elts(n+1)%children(modulo((c+1)-1,4)+1) 
       get_child_and_neigh_patches(3) = dom%patch%elts(n+1)%children(modulo((c+2)-1,4)+1) 
    endif

    n = dom%patch%elts(p_par+1)%neigh(c+4) ! corner
    if (n > 0) get_child_and_neigh_patches(4) = dom%patch%elts(n+1)%children(modulo((c+2)-1,4)+1) 
  end function get_child_and_neigh_patches  

  integer function side (dom, p, s)
    implicit none
    type(Domain) :: dom
    integer      :: p, s

    integer :: n

    n = dom%patch%elts(p+1)%neigh(s+1)
    if (n >= 0) then
       side = -s - 1
       return
    else
       side = dom%bdry_patch%elts(-n+1)%side
       return
    end if
  end function side

  subroutine connect_cousin (dom, p_par, p_chd, s_par, s_chd, c)
    implicit none
    type(Domain) :: dom
    integer      :: c, p_par, p_chd, s_par, s_chd

    integer :: n, typ
    !  c: which child on neighbour

    n = dom%patch%elts(p_par+1)%neigh(s_par+1)
    if (n < 0) then
       n = -n
       typ = dom%bdry_patch%elts(n+1)%side
       if (.not. s_chd == typ - 1) return
       if (dom%neigh(typ) == POLE) then
          call connect_pole (dom, n, p_chd, s_par)
          return
       end if
       call append (dom%send_pa_all, n)
       call append (dom%send_pa_all, c)
       call append (dom%send_pa_all, p_chd)
       call append (dom%send_pa_all, s_par)
    end if
  end subroutine connect_cousin

  subroutine connect_pole (dom, n, p_chd, s_par)
    implicit none
    type(Domain) :: dom
    integer      :: n, p_chd, s_par

    integer :: i
    do i = 1, 2
       call append (dom%send_pa_all, n)
       call append (dom%send_pa_all, i - 1)
       call append (dom%send_pa_all, p_chd)
       call append (dom%send_pa_all, s_par)
    end do
  end subroutine connect_pole

  subroutine connect_children (dom, p_par)
    ! children of patch `p_par` are connected to neighbours on same level if they exist
    !        and temporary boundaries are removed
    !        considers the case that not all four children are present
    ! update: still used to connect old patches to new patches (new patches already connected now)
    implicit none
    type(Domain) :: dom
    integer      :: p_par

    integer, dimension(N_CHDRN) :: children
    integer                     :: c, n_chd, n_tmp, p_chd, s

    children = dom%patch%elts(p_par+1)%children
    do c = 1, N_CHDRN
       p_chd = children(c)
       if (p_chd == 0) cycle
       do s = 1, N_BDRY
          n_tmp = dom%patch%elts(p_chd+1)%neigh(s)
          if (n_tmp >= 1) cycle ! already connected
          n_chd = find_neigh_patch_Domain(dom, p_par, c-1, s-1)
          if (n_chd == 0) then
             cycle
          else
             if (.not. n_chd == n_tmp) then
                dom%bdry_patch%elts(-n_tmp+1)%side = 0
                dom%patch%elts(p_chd+1)%neigh(s) = n_chd
             end if
          end if
       end do
    end do
  end subroutine connect_children

  subroutine set_areas (dom, id, val)
    implicit none
    type(Domain)          :: dom
    integer               :: id
    real(8), dimension(7) :: val

    real(8), dimension(4) :: area

    area = val(1:4)
    if (id < 0) area = (/area(2), area(1), area(4), area(3)/)

    dom%overl_areas%elts(abs(id) + 1)%a     = area
    dom%overl_areas%elts(abs(id) + 1)%split = val(5:6)
    dom%areas%elts(abs(id) + 1)%hex_inv     = val(7)
  end subroutine set_areas

  subroutine get_areas (dom, id, val)
    implicit none
    real(8), dimension(7), intent(out) :: val
    type(Domain)                       :: dom
    integer                            :: id

    real(8), dimension(7) :: area

    area = 0d0
    
    area(1:4) = dom%overl_areas%elts(id+1)%a
    area(5:6) = dom%overl_areas%elts(id+1)%split
    area(7)   = dom%areas%elts(id+1)%hex_inv
    val       = area
    return
  end subroutine get_areas

  subroutine area_post_comm (dom, p, c, offs, dims, zlev)
    implicit none
    type(Domain)                   :: dom
    integer                        :: c, p, zlev
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    integer, dimension(N_BDRY+1)   :: offs

    if (c == IPLUSJMINUS) then
       id = idx (PATCH_SIZE, -1, offs, dims)
       dom%overl_areas%elts(id+1)%a = 0d0
    end if
    if (c == IMINUSJPLUS) then
       id = idx (-1, PATCH_SIZE, offs, dims)
       dom%overl_areas%elts(id+1)%a = 0d0
    end if
  end subroutine area_post_comm

  subroutine attach_bdry (dom, p_par, c, s, side)
    implicit none
    type(Domain) :: dom
    integer      :: p_par, c, s, side

    integer :: n_chd, p_chd

    n_chd = find_neigh_patch_Domain(dom, p_par, c, s)

    if (n_chd == 0) n_chd = - add_bdry_patch_Domain (dom, side)

    p_chd = dom%patch%elts(p_par+1)%children(c+1)
    dom%patch%elts(p_chd+1)%neigh(s+1) = n_chd
  end subroutine attach_bdry

  subroutine init_refine_patch_mod
    implicit none
    logical :: initialized = .false.
    
    if (initialized) return ! initialize only once

    call init_shared_mod
    call init_domain_mod
    call init_init_mod
    call init_wavelet_mod
    call init_mask_mod
    
    initialized = .true.
  end subroutine init_refine_patch_mod

  subroutine fill_up_level
    ! Fills up level level_start+1 and increases level_start
    implicit none
    integer :: d, j, p_par, c, p_chd

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(level_start)%length
          p_par = grid(d)%lev(level_start)%elts(j)
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par+1)%children(c)
             if (p_chd == 0) call refine_patch (grid(d), p_par, c-1)
          end do
       end do
    end do
    call post_refine
    level_start = level_start+1
  end subroutine fill_up_level

  subroutine refine_patch (dom, p, c0)
    implicit none
    type(Domain) :: dom
    integer      :: p, c0

    call refine_patch1 (dom, p, c0)
    call refine_patch2 (dom, p, c0)
  end subroutine refine_patch
end module refine_patch_mod
