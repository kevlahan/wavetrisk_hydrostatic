module comm_mod
  use arch_mod
  use domain_mod
  implicit none
  integer, dimension(4,4)            :: shift_arr
  integer, dimension(:), allocatable :: n_active_edges, n_active_nodes
  real(8)                            :: beta_sclr_loc, beta_divu_loc, beta_rotu_loc, dt_loc, min_mass_loc, sync_val
contains
  subroutine init_comm_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_arch_mod
    call init_domain_mod
    shift_arr = reshape ((/0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0/), (/4, 4/))
    initialized = .true.
  end subroutine init_comm_mod

  subroutine init_comm
    implicit none
    integer :: k, s, d

    allocate(n_active_edges(min_level-1:max_level), n_active_nodes(min_level-1:max_level))

    n_active_edges = 0
    n_active_nodes = 0
    do d = 1, size(grid)
       do s = 1, N_BDRY
          if (.not. is_penta(grid(d), 1, s - 1)) then
             call create_comm (grid(d), 1, s)
          else
             do k = 1, 2
                call create_comm_pole (grid(d), 1, s, k - 1)
             end do
          end if
       end do
    end do
  end subroutine init_comm

  subroutine comm_nodes9 (get, set)
    implicit none
    external :: get, set
    
    integer               :: i, dest_glo, dest_id, dest_loc, src_glo, src_id, src_loc
    real(8), dimension(7) :: val

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
             src_id = grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)
             call get (grid(src_loc), src_id, val)
             call set (grid(dest_loc), dest_id, val)
          end do
       end do
    end do
  end subroutine comm_nodes9

  subroutine create_pack_st2 (dom, src, i, j, pa, e, id, e_pack, orient)
    implicit none
    type(Domain) :: dom
    integer      :: src, i, j, u, pa, e, id, e_pack, orient

    ! if PATCH_SIZE == 4 (i.e. PATCH_LEVEL == 2) and BDRY_THICKNESS == 3
    ! then some halo cells have their partner not on a neighbouring patch
    ! luckily these (very few) halo cells are not used and can be skipped:

    if (i < 0 .or. j < 0 .or. i >= PATCH_SIZE .or. j >= PATCH_SIZE) return

    if (e == NODE) then
       call create_pack_st (dom, AT_NODE, src, i, j, pa, e, id*orient)
    else
       call create_pack_st (dom, AT_EDGE, src, i, j, pa, e_pack, orient*(EDGE*id + e))
    end if
  end subroutine create_pack_st2

  subroutine comm_nodes3 (get, set)
    implicit none
    external    :: get, set
    type(Coord) :: get
    integer     :: i, dest_id, dest_loc, dest_glo, src_glo, src_id, src_loc

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
             src_id = grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)
             call set (grid(dest_loc), dest_id, get(grid(src_loc), src_id))
          end do
       end do
    end do
  end subroutine comm_nodes3

  subroutine create_comm_e (dom, p, s, e)
    implicit none
    type(Domain)                   :: dom
    integer                        :: b, e, e_pack, i, i_pack, i_recv, j, j_pack, j_recv, ngb_pa
    integer                        :: orient, p, rot, s, s_adj, shift, src, t_last, t_next, typ
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call get_offs_Domain (dom, p, offs, dims)

    b = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side
    if (typ < 1)  return

    src = dom%neigh(typ)

    if (src == NONE) return

    rot = dom%neigh_rot(typ)
    ngb_pa = dom%bdry_patch%elts(b+1)%neigh
    orient = (-1)**rot
    if (s > 4) then
       t_last = typ - 4
       t_next = modulo(typ - 4, 4) + 1
       if (dom%neigh_rot(t_last) == 1) then
          s_adj = s - 4
          shift = shift_arr(e+1,s_adj)
          e_pack = modulo(e + rot*(modulo(s_adj, 2) + 1), 3)
       else
          if (dom%neigh_rot(t_next) == 1) then
             s_adj = modulo(s - 4, 4) + 1
             shift = shift_arr(e+1,s_adj)
             e_pack = modulo(e + rot*(modulo(s_adj, 2) + 1), 3)
          else
             shift = 0
             e_pack = e
          end if
       end if
    else
       shift = shift_arr(e+1,s)
       e_pack = modulo(e + rot*(modulo(s, 2) + 1), 3)
    end if

    if (s == NORTH) then
       do j = 1, BDRY_THICKNESS
          do i = 1, PATCH_SIZE
             call pack_idx (i - 1, j - 1, rot, s, shift, i_recv, j_recv, i_pack, j_pack)
             call create_pack_st2 (dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do
    end if

    if (s == EAST) then
       do i = 1, BDRY_THICKNESS
          do j = 1, PATCH_SIZE
             call pack_idx (i - 1, j - 1, rot, s, shift, i_recv, j_recv, i_pack, j_pack)
             call create_pack_st2 (dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do
    end if

    if (s == SOUTH) then
       do j = 1, BDRY_THICKNESS
          do i = 1, PATCH_SIZE
             call pack_idx(i - 1, j - 1, rot, s, shift, i_recv, j_recv, i_pack, j_pack)
             call create_pack_st2 (dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do
    end if

    if (s == WEST) then
       do i = 1, BDRY_THICKNESS
          do j = 1, PATCH_SIZE
             call pack_idx (i - 1, j - 1, rot, s, shift, i_recv, j_recv, &
                  i_pack, j_pack)
             call create_pack_st2 (dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do
    end if

    if (s == NORTHEAST) then
       if (dom%neigh_rot(t_last) == 1) then
          do j = 1, BDRY_THICKNESS
             do i = 1, BDRY_THICKNESS - rot*(j + shift - 1)
                call pack_idx(i - 1, j - 1, rot, s_adj, shift, i_recv, j_recv, i_pack, j_pack)
                i_recv = i_recv + PATCH_SIZE
                call create_pack_st2 (dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), e_pack, orient)
             end do
          end do
       else
          if (dom%neigh_rot(t_next) == 1) then
             do i = 1, BDRY_THICKNESS
                do j = 1, BDRY_THICKNESS - rot*(i + shift - 1)
                   call pack_idx (i - 1, j - 1, rot, s_adj, shift, &
                        i_recv, j_recv, i_pack, j_pack)
                   j_recv = j_recv + PATCH_SIZE
                   call create_pack_st2 (dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          else
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS
                   i_recv = i - 1 + PATCH_SIZE
                   j_recv = j - 1 + PATCH_SIZE
                   i_pack = i - 1
                   j_pack = j - 1
                   call create_pack_st2 (dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          end if
       end if
    end if

    if (s == SOUTHEAST) then
       if (dom%neigh_rot(t_last) == 1) then
          do i = 1, BDRY_THICKNESS
             do j = 1, BDRY_THICKNESS + rot*(i + shift - 1)
                i_recv = i - 1 + PATCH_SIZE
                j_recv = -j + rot*(i + shift - 1)
                i_pack = j - 1
                j_pack = i - 1
                call create_pack_st2 (dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), e_pack, orient)
             end do
          end do
       else
          if (dom%neigh_rot(t_next) == 1) then
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS + rot*(j + shift - 1)
                   i_recv = (PATCH_SIZE + i - 1) - rot*(j + shift - 1)
                   j_recv = -1 - (j - 1)
                   i_pack = LAST - (j - 1)
                   j_pack = LAST - (i - 1)
                   call create_pack_st2(dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), e_pack, orient)
                end do
             end do
          else
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS
                   i_recv = i - 1 + PATCH_SIZE
                   j_recv = -1 - (j - 1)
                   i_pack = i - 1
                   j_pack = LAST - (j - 1)
                   call create_pack_st2 (dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), e_pack, orient)
                end do
             end do
          end if
       end if
    end if

    if (s == SOUTHWEST) then
       if (dom%neigh_rot(t_last) == 1) then
          do j = 1, BDRY_THICKNESS
             do i = 1, BDRY_THICKNESS - rot*(j + shift - 1)
                i_recv = -i - rot*(j + shift - 1)
                j_recv = -1 - (j - 1)
                i_pack = LAST - (j - 1)
                j_pack = i - 1
                call create_pack_st2 (dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), e_pack, orient)
             end do
          end do
       else
          if (dom%neigh_rot(t_next) == 1) then
             do i = 1, BDRY_THICKNESS
                do j = 1, BDRY_THICKNESS - rot*(i + shift - 1)
                   i_recv = -1 - (i - 1)
                   j_recv = -j - rot*(i + shift - 1)
                   i_pack = j - 1
                   j_pack = LAST - (i - 1)
                   call create_pack_st2 (dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), e_pack, orient)
                end do
             end do
          else
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS - rot*(j + shift - 1)
                   i_recv = -1 - (i - 1)
                   j_recv = -1 - (j - 1)
                   i_pack = LAST - (i - 1)
                   j_pack = LAST - (j - 1)
                   call create_pack_st2 (dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), e_pack, orient)
                end do
             end do
          end if
       end if
    end if

    if (s == NORTHWEST) then
       if (dom%neigh_rot(t_last) == 1) then
          do i = 1, BDRY_THICKNESS
             do j = 1, BDRY_THICKNESS + rot*(i + shift - 1)
                i_recv = -1 - (i - 1)
                j_recv = (PATCH_SIZE + j - 1) - rot*(i + shift - 1)
                i_pack = LAST - (j - 1)
                j_pack = LAST - (i - 1)
                call create_pack_st2 (dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), e_pack, orient)
             end do
          end do
       else
          if (dom%neigh_rot(t_next) == 1) then
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS + rot*(j + shift - 1)
                   i_recv = -i + rot*(j + shift - 1)
                   j_recv = j - 1 + PATCH_SIZE
                   i_pack = j - 1
                   j_pack = i - 1
                   call create_pack_st2 (dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), e_pack, orient)
                end do
             end do
          else
             do i = 1, BDRY_THICKNESS
                do j = 1, BDRY_THICKNESS
                   i_recv = -1 - (i - 1)
                   j_recv = j - 1 + PATCH_SIZE
                   i_pack = LAST - (i - 1)
                   j_pack = j - 1
                   call create_pack_st2 (dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), e_pack, orient)
                end do
             end do
          end if
       end if
    end if
  end subroutine create_comm_e

  integer function rot_direction (dom, typ)
    implicit none
    type(Domain) :: dom
    integer      :: typ

    integer :: t_last, t_next

    if (typ <= 4) then
       rot_direction = modulo(typ, 2)
       return
    else
       t_last = typ - 4
       if (dom%neigh_rot(t_last) == 1) then
          rot_direction = modulo(t_last, 2)
          return
       else
          t_next = (modulo(typ, 4) + 1) - 4
          rot_direction = modulo(t_next, 2)
          return
       end if
    end if
  end function rot_direction

  subroutine comm_communication
    implicit none
    integer               :: i, dest_glo, dest_loc, src_glo, src_loc
    integer, dimension(4) :: st

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%send_conn(dest_glo+1)%length, 4
             st = grid(src_loc)%send_conn(dest_glo+1)%elts(i:i + 3)
             call unpack_comm_struct (grid(dest_loc), src_glo, st(1), st(2), st(3), st(4))
          end do
          grid(src_loc)%send_conn(dest_glo+1)%length = 0
       end do
    end do
  end subroutine comm_communication

  type(Coord) function get_coord (dom, id)
    type(Domain) :: dom
    integer      :: id

    get_coord = dom%node%elts(id+1)
  end function get_coord

  subroutine create_comm (dom, p, s)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, s
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: b, e, pa, src, rot, typ
    
    call get_offs_Domain (dom, p, offs, dims)

    b = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side

    if (typ < 1) return

    src = dom%neigh(typ)
    rot = dom%neigh_rot(typ)
    pa = dom%bdry_patch%elts(b+1)%neigh

    do e = 1, EDGE + 1
       call create_comm_e (dom, p, s, e - 1)
    end do

    if (rot == 1 .and. typ == WEST) then
       if (is_penta(dom, p, SOUTHWEST - 1)) then
          call create_pack_st (dom, AT_EDGE, src, LAST, LAST, pa, RT, &
               nidx(LAST_BDRY, LAST_BDRY, SOUTHWEST, offs, dims)*EDGE)
       end if
    end if

    if (rot == 1 .and. typ == SOUTH) then
       if (is_penta(dom, p, SOUTHWEST - 1)) then
          call create_pack_st(dom, AT_EDGE, src, LAST, LAST, pa, UP, &
               -nidx(LAST_BDRY, LAST_BDRY, SOUTHWEST, offs, dims)*EDGE)
       end if
    end if

    if (rot == 0 .and. typ == EAST) then
       if (is_penta(dom, p, SOUTHEAST - 1)) then
          call create_pack_st(dom, AT_NODE, src, 1, 0, pa, NODE, nidx(0, LAST_BDRY, SOUTHEAST, offs, dims))
       end if
    end if

    if (rot == 0 .and. typ == NORTH) then
       if (is_penta(dom, p, NORTHWEST - 1)) then
          call create_pack_st(dom, AT_NODE, src, 0, 1, pa, NODE, nidx(LAST_BDRY, 0, NORTHWEST, offs, dims))
       end if
    end if

    if (rot == 1 .and. typ == EAST) then
       if (is_penta(dom, p, NORTHEAST - 1)) then
          call create_pack_st (dom, AT_NODE, src, 0, 1, pa, NODE,  nidx(0, 1, NORTHEAST, offs, dims))
          call create_pack_st (dom, AT_NODE, src, 1, 1, pa, NODE, -nidx(1, 0, NORTHEAST, offs, dims))
          call create_pack_st (dom, AT_EDGE, src, 0, 1, pa, RT,    nidx(0, 1, NORTHEAST, offs, dims)*EDGE + RT)
          call create_pack_st (dom, AT_EDGE, src, 0, 0, pa, DG,  -(nidx(0, 0, NORTHEAST, offs, dims)*EDGE + RT))
          call create_pack_st (dom, AT_EDGE, src, 0, 0, pa, UP,    nidx(0, 0, NORTHEAST, offs, dims)*EDGE + UP)
       end if
    end if

    if (rot == 1 .and. typ == NORTH) then
       if (is_penta(dom, p, NORTHEAST - 1)) then
          call create_pack_st (dom, AT_NODE, src, 1, 0, pa, NODE,  nidx(1, 0, NORTHEAST, offs, dims))
          call create_pack_st (dom, AT_NODE, src, 1, 1, pa, NODE, -nidx(0, 1, NORTHEAST, offs, dims))
          call create_pack_st (dom, AT_EDGE, src, 1, 0, pa, UP,  -(nidx(0, 1, NORTHEAST, offs, dims)*EDGE + RT))
          call create_pack_st (dom, AT_EDGE, src, 0, 0, pa, DG,  -(nidx(0, 0, NORTHEAST, offs, dims)*EDGE + UP))
          call create_pack_st (dom, AT_EDGE, src, 0, 0, pa, RT,    nidx(0, 0, NORTHEAST, offs, dims)*EDGE + RT)
       end if
    end if
  end subroutine create_comm

  subroutine create_comm_pole (dom, p, s, i)
    implicit none
    type(Domain) :: dom
    integer      :: i, p, s

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: b, lev, pa, s_side, src, typ
    integer, dimension(2)          :: ij_node, ij_send

    !  side s of patch p is a pole and i-th connection over this pole is available
    b = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side

    if (typ < 1) return

    if (.not. dom%neigh(typ) == POLE) return

    lev = dom%patch%elts(p+1)%level
    src = dom%neigh_over_pole(i+1)
    pa = dom%neigh_pa_over_pole%elts(i+2*lev+1)

    call get_offs_Domain (dom, p, offs, dims)

    if (s == NORTHWEST) then
       s_side = NORTH
    else
       if (s == SOUTHEAST) s_side = EAST
    end if

    if (i == 0) then
       ij_node = (/0, 1/)
       if (s == SOUTHEAST) ij_node = (/ij_node(2), ij_node(1)/)
       call create_pack_st(dom, AT_EDGE, src, 0, LAST, pa, DG, &
            nidx(ij_node(1), ij_node(2), s_side, offs, dims)*EDGE + 2*s_side - 2)
       ij_node = (/LAST_BDRY, 1/)
       ij_send = (/1, LAST/)
       if (s == SOUTHEAST) then
          ij_node = (/ij_node(2), ij_node(1)/)
          ij_send = (/ij_send(2), ij_send(1)/)
       end if
       call create_pack_st (dom, AT_NODE, src, ij_send(1), ij_send(2), pa, NODE, -nidx(ij_node(1), ij_node(2), s, offs, dims))
    end if

    if (i == 1) then
       ij_node = (/LAST_BDRY, 0/)
       ij_send = (/0, LAST/)
       if (s == SOUTHEAST) then
          ij_node = (/ij_node(2), ij_node(1)/)
          ij_send = (/ij_send(2), ij_send(1)/)
       end if
       call create_pack_st (dom, AT_NODE, src, ij_send(1), ij_send(2), pa, NODE, nidx(ij_node(1), ij_node(2), s, offs, dims))
    end if

    ij_node = (/-i + 1, 1/)
    ij_send = (/0, LAST/)

    if (s == SOUTHEAST) then
       ij_node = (/ij_node(2), ij_node(1)/)
       ij_send = (/ij_send(2), ij_send(1)/)
    end if

    call create_pack_st (dom, AT_NODE, src, ij_send(1), ij_send(2), pa, &
         NODE, (-1)**i*nidx(ij_node(1), ij_node(2), s_side, offs, dims))

    call create_pack_st (dom, AT_EDGE, src, ij_send(1), ij_send(2), pa, UP &
         - 2*s_side + 2, (-1)**i*(nidx(0, 0, s_side, offs, dims)*EDGE + DG + i*(-2*s_side + 3)))

    if (i == 0) return

    if (s == NORTHWEST) then
       s_side = WEST
    else
       if (s == SOUTHEAST) s_side = SOUTH
    end if

    ij_node = (/LAST_BDRY - 1, LAST/)
    ij_send = (/1, LAST/)

    if (s == SOUTHEAST) then
       ij_node = (/ij_node(2), ij_node(1)/)
       ij_send = (/ij_send(2), ij_send(1)/)
    end if

    call create_pack_st (dom, AT_NODE, src, ij_send(1), ij_send(2), pa, NODE, nidx(ij_node(1), ij_node(2), s_side, offs, dims))

    ij_node = (/LAST, LAST_BDRY/)

    if (s == NORTHWEST) then
       ij_node = (/ij_node(2), ij_node(1)/)
    end if

    call create_pack_st (dom, AT_EDGE, src, LAST*(-s_side + 4), &
         LAST*(s_side - 3), pa, DG, nidx(ij_node(1), ij_node(2), s_side, offs, dims)*EDGE + 2*s_side - 6)
  end subroutine create_comm_pole

  subroutine get_areas (dom, id, val)
    implicit none
    real(8), dimension(7), intent(out) :: val
    type(Domain)                       :: dom
    integer                            :: id
    
    real(8), dimension(7) :: area

    area = 0.0_8
    area(1:4) = dom%overl_areas%elts(id+1)%a
    area(5:6) = dom%overl_areas%elts(id+1)%split
    area(7) = dom%areas%elts(id+1)%hex_inv
    val = area
    return
  end subroutine get_areas

  subroutine comm_masks
    implicit none
    integer dest_glo, dest_id, dest_loc, i, k, src_glo, src_id, src_loc

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
             src_id  = grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)
             grid(dest_loc)%mask_n%elts(abs(dest_id)+1) = grid(src_loc)%mask_n%elts(abs(src_id)+1) 
          end do
          do i = 1, grid(src_loc)%pack(AT_EDGE,dest_glo+1)%length
             src_id  = grid(src_loc)%pack(AT_EDGE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_EDGE,src_glo+1)%elts(i)
             grid(dest_loc)%mask_e%elts(abs(dest_id)+1) = grid(src_loc)%mask_e%elts(abs(src_id)+1) 
          end do
       end do
    end do
  end subroutine comm_masks

  subroutine comm_edges (get, set)
    implicit none
    external :: get, set
    
    real(8)  :: get
    integer  :: dest_glo, dest_id, dest_loc, i, src_id, src_glo, src_loc

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_EDGE,dest_glo+1)%length
             src_id  = grid(src_loc)%pack(AT_EDGE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_EDGE,src_glo+1)%elts(i)
             call set (grid(dest_loc), dest_id, get(grid(src_loc), src_id))
          end do
       end do
    end do
  end subroutine comm_edges

  subroutine unpack_comm_struct (dom, src, i, j, p, e)
    implicit none
    type(Domain) :: dom
    integer      :: e, i, j, p, src

    integer :: offs

    if (e == NODE) then
       offs = dom%patch%elts(p+1)%elts_start
       call append (dom%pack(AT_NODE,src+1), idx__fast(i, j, offs))
    end if

    if (.not. e == NODE) then
       offs = dom%patch%elts(p+1)%elts_start
       call append (dom%pack(AT_EDGE,src+1), idx__fast(i, j, offs)*EDGE + e)
    end if
  end subroutine unpack_comm_struct

  subroutine pack_idx (i, j, rot, s, shift, i_recv, j_recv, i_pack, j_pack)
    implicit none
    integer :: i, i_pack, i_recv, j, j_pack, j_recv, rot, shift, s

    if (s == NORTH) then
       i_recv = i + rot*(j + shift)
       j_recv = j + PATCH_SIZE
       if (rot == 1) then
          i_pack = j
          j_pack = LAST - i
       else
          i_pack = i
          j_pack = j
       end if
    end if

    if (s == EAST) then
       i_recv = i + PATCH_SIZE
       j_recv = j + rot*(i + shift)
       if (rot == 1) then
          i_pack = LAST - j
          j_pack = i
       else
          i_pack = i
          j_pack = j
       end if
    end if

    if (s == SOUTH) then
       i_recv = i - rot*(j + shift)
       j_recv = -1 - j
       if (rot == 1) then
          i_pack = LAST - j
          j_pack = LAST - i
       else
          i_pack = i
          j_pack = LAST - j
       end if
    end if

    if (s == WEST) then
       i_recv = -1 - i
       j_recv = j - rot*(i + shift)
       if (rot == 1) then
          i_pack = LAST - j
          j_pack = LAST - i
       else
          i_pack = LAST - i
          j_pack = j
       end if
    end if
  end subroutine pack_idx

  subroutine update_comm (dom)
    implicit none
    type(Domain) :: dom
    
    integer c, ii, lev, n_chd, n_par, ngb_pa, p_chd, p_par, s, src_glo, typ,  unused_elements

    do src_glo = 1, N_GLO_DOMAIN
       unused_elements = 0
       do ii = 1, dom%recv_pa(src_glo)%length, 4
          ngb_pa = dom%recv_pa(src_glo)%elts(ii)
          c = dom%recv_pa(src_glo)%elts(ii+1)
          p_par = dom%recv_pa(src_glo)%elts(ii+2)
          s = dom%recv_pa(src_glo)%elts(ii+3)
          n_par = -dom%patch%elts(p_par+1)%neigh(s+1)
          typ = dom%bdry_patch%elts(n_par+1)%side
          if (typ < 1) return
          if (dom%neigh(typ) == POLE) then
             p_chd = dom%patch%elts(p_par+1)%children(s-3)
          else
             p_chd = dom%patch%elts(p_par+1)%children(c+1)
          end if
          if (p_chd == 0) then
             dom%recv_pa(src_glo)%elts(unused_elements + 1:unused_elements + 4) = (/ngb_pa, c, p_par, s/)
             unused_elements = unused_elements + 4
             cycle
          end if
          if (dom%neigh(typ) == POLE) then
             lev = dom%patch%elts(p_chd+1)%level
             if (2*lev + 2 > dom%neigh_pa_over_pole%length) then
                call extend (dom%neigh_pa_over_pole, 2, 0)
             end if
             dom%neigh_pa_over_pole%elts((1-c)+2*lev+1) = ngb_pa
             call create_comm_pole (dom, p_chd, s + 1, 1-c)
             cycle
          end if
          n_chd = -dom%patch%elts(p_chd+1)%neigh(s+1)
          dom%bdry_patch%elts(n_chd+1)%neigh = ngb_pa
          if (.not. is_penta(dom, p_chd, s)) call create_comm (dom, p_chd, s + 1)
       end do
       dom%recv_pa(src_glo)%length = unused_elements
    end do
  end subroutine update_comm

  subroutine area_post_comm (dom, p, c, offs, dims, zlev)
    implicit none
    type(Domain)                   :: dom
    integer                        :: c, p, zlev
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    integer, dimension(N_BDRY+1)   :: offs
    
    if (c == IPLUSJMINUS) then
       id = idx(PATCH_SIZE, -1, offs, dims)
       dom%overl_areas%elts(id+1)%a = 0.0
    end if
    if (c == IMINUSJPLUS) then
       id = idx(-1, PATCH_SIZE, offs, dims)
       dom%overl_areas%elts(id+1)%a = 0.0
    end if
  end subroutine area_post_comm

  subroutine comm_patch_conn
    implicit none
    integer :: b, c, dest, i, l_par, ngh_pa, p_chd, rot, rot_shift, s, src, src_glo, typ, unused_elements
    logical :: is_pole

    do src = 1, size(grid)
       unused_elements = 0
       src_glo = glo_id(rank+1,src)
       do i = 1, grid(src)%send_pa_all%length, 4
          b = grid(src)%send_pa_all%elts(i)
          c = grid(src)%send_pa_all%elts(i+1)
          p_chd = grid(src)%send_pa_all%elts(i+2)
          s = grid(src)%send_pa_all%elts(i+3)
          typ = grid(src)%bdry_patch%elts(b+1)%side

          if (typ < 1) return

          dest = grid(src)%neigh(typ)
          is_pole = dest == POLE

          if (is_pole) then
             dest = grid(src)%neigh_over_pole(c+1)
             l_par = grid(src)%patch%elts(p_chd+1)%level - 1
             if (grid(src)%neigh_pa_over_pole%length < l_par*2 + c + 1) then
                ngh_pa = 0
             else
                ngh_pa = grid(src)%neigh_pa_over_pole%elts(l_par*2 + c + 1)
             end if
          else
             ngh_pa = grid(src)%bdry_patch%elts(b+1)%neigh
          end if

          if (ngh_pa == 0) then
             grid(src)%send_pa_all%elts(unused_elements + &
                  1:unused_elements + 4) = (/b, c, p_chd, s/)
             unused_elements = unused_elements + 4
             cycle
          else 
             if (dest == NONE) cycle
          end if

          ! cycle if dest is not on rank
          if (.not. owner(dest+1) == rank) cycle

          dest = loc_id(dest+1)
          rot = grid(src)%neigh_rot(typ)
          rot_shift = (rot_direction(grid(src), typ)*2 - 1)*rot

          call append (grid(dest+1)%recv_pa(src_glo+1), p_chd)

          if (is_pole) then
             call append (grid(dest+1)%recv_pa(src_glo+1), c)
          else
             call append (grid(dest+1)%recv_pa(src_glo+1), modulo(c + rot_shift, 4))
          end if

          call append (grid(dest+1)%recv_pa(src_glo+1), ngh_pa)

          if (is_pole) then
             call append (grid(dest+1)%recv_pa(src_glo+1), s)
          else
             call append (grid(dest+1)%recv_pa(src_glo+1), modulo(rot_shift + s + 2, 4) + 4*(s/4))
          end if
       end do
       grid(src)%send_pa_all%length = unused_elements
    end do
  end subroutine comm_patch_conn

  subroutine create_pack_st (dom, unpk_pos, src, i, j, pa, e, id)
    implicit none
    type(Domain) :: dom
    integer :: e, i, id, j, pa, src, unpk_pos

    call append (dom%send_conn(src+1),      i)
    call append  (dom%send_conn(src+1),     j)
    call append (dom%send_conn(src+1),     pa)
    call append (dom%send_conn(src+1),      e)
    call append (dom%unpk(unpk_pos,src+1), id)
  end subroutine create_pack_st

  subroutine set_coord (dom, id, val)
    implicit none
    type(Domain) :: dom
    integer      :: id
    type(Coord)  :: val

    dom%node%elts(abs(id) + 1) = val
  end subroutine set_coord

  subroutine cp_bdry_inside (field)
    implicit none
    type(Float_Field) :: field
    
    integer ::  dest_id, dest_loc, dest_glo, i, pos, src_id, src_loc, src_glo

    pos = field%pos
    
    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
             src_id  = grid(src_loc)%pack(pos,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)
             field%data(dest_loc)%elts(abs(dest_id)+1) = field%data(src_loc)%elts(src_id+1)
             if (dest_id < 0 .and. pos == AT_EDGE) &
                  field%data(dest_loc)%elts(abs(dest_id)+1) = - field%data(dest_loc)%elts(abs(dest_id)+1)
          end do
       end do
    end do
  end subroutine cp_bdry_inside

  subroutine cp_bdry_inside_vector (field)
    implicit none
    type(Float_Field), dimension(:) :: field
    
    integer :: i, i1, dest_id, dest_loc, dest_glo, pos, src_id, src_loc, src_glo

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i1 = 1, size(field)
             pos = field(i1)%pos
             do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
                src_id  = grid(src_loc)%pack(pos,dest_glo+1)%elts(i)
                dest_id = grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)
                field(i1)%data(dest_loc)%elts(abs(dest_id)+1) = field(i1)%data(src_loc)%elts(src_id+1)
                if (dest_id < 0 .and. pos == AT_EDGE) &
                     field(i1)%data(dest_loc)%elts(abs(dest_id)+1) = - field(i1)%data(dest_loc)%elts(abs(dest_id)+1)
             end do
          end do
       end do
    end do
  end subroutine cp_bdry_inside_vector

  subroutine cp_bdry_inside_array (field)
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: i, i1, i2, dest_id, dest_loc, dest_glo, pos, src_id, src_loc, src_glo

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i2 = 1, size(field,2)
             do i1 = 1, size(field,1)
                pos = field(i1,i2)%pos
                do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
                   src_id  = grid(src_loc)%pack(pos,dest_glo+1)%elts(i)
                   dest_id = grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)
                   field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1) = field(i1,i2)%data(src_loc)%elts(src_id+1)
                   if (dest_id < 0 .and. pos == AT_EDGE) &
                        field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1) = - field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1)
                end do
             end do
          end do
       end do
    end do
  end subroutine cp_bdry_inside_array

  subroutine comm_nodes (get, set)
    implicit none
    external :: get, set

    real(8) get
    integer :: i, dest_id, dest_glo, dest_loc, src_id, src_glo, src_loc

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
             src_id = grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)
             call set (grid(dest_loc), dest_id, get(grid(src_loc), src_id))
          end do
       end do
    end do
  end subroutine comm_nodes

  subroutine set_areas (dom, id, val)
    implicit none
    type(Domain)          :: dom
    integer               :: id
    real(8), dimension(7) :: val
    
    real(8), dimension(4) :: area

    area = val(1:4)
    if (id < 0) area = (/area(2), area(1), area(4), area(3)/)

    dom%overl_areas%elts(abs(id) + 1)%a = area
    dom%overl_areas%elts(abs(id) + 1)%split = val(5:6)
    dom%areas%elts(abs(id) + 1)%hex_inv = val(7)
  end subroutine set_areas

  subroutine min_dt (dom, i, j, zlev, offs, dims)
    ! Calculates time step and number of active nodes and edges
    ! time step is smallest of barotropic time step, advective time step and internal wave time step
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer            :: d, e, id, id_e, id_i, k, l
    real(8)            :: dx, v_mag
    real(8), parameter :: cfl = 0.8

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d  = dom%id + 1
    l  = dom%level%elts(id_i)
        
    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       n_active_nodes(l) = n_active_nodes(l) + 1 
       if (adapt_dt) then 
          dx = minval (dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1))
          do k = 1, zlevels
             v_mag = velo_mag (dom, i, j, k, offs, dims)
             if (mode_split) then
                dt_loc = min (dt_loc, dt_init, cfl_num*dx/wave_speed, cfl*dx/v_mag, dx/c1)
             else
                dt_loc = min (dt_loc, cfl_num*dx/(v_mag + wave_speed))
             end if
          end do
       end if
    end if

    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) n_active_edges(l) = n_active_edges(l) + 1
    end do
  end subroutine min_dt

  real(8) function velo_mag (dom, i, j, zlev, offs, dims)
    ! Calculate magnitude of velocity as sqrt(2 kinetic_energy)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, idS, idSW, idW
    real(8) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(8) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    d  = dom%id + 1

    u_prim_RT = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
    u_prim_DG = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
    u_prim_UP = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

    u_dual_RT = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_DG = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

    u_prim_UP_S  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
    u_prim_DG_SW = sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
    u_prim_RT_W  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
    
    u_dual_RT_W  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
    u_dual_DG_SW = sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
    u_dual_UP_S  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

    velo_mag = sqrt ((u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT +  &
                      u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
                      * dom%areas%elts(id+1)%hex_inv/2)
  end function velo_mag

  subroutine cal_min_mass (dom, i, j, zlev, offs, dims)
    ! Calculates minimum mass and diffusion stability limits
    use, intrinsic :: ieee_arithmetic
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, id_e, id_i, k, l
    real(8) :: col_mass, d_e, fac, full_mass, init_mass

    id   = idx (i, j, offs, dims)
    id_i = id + 1
    d    = dom%id + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       col_mass = 0.0_8
       do k = 1, zlevels
          full_mass = sol(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)
          if (full_mass < 0.0_8 .or. ieee_is_nan (full_mass)) then
             write (6,'(A,i8,A,i2,2(A,i2),A)') "Mass negative at id = ", id_i, " with scale j = ", dom%level%elts(id_i),  &
                  " vertical level k = ", k, " and mask = ", dom%mask_n%elts(id_i), " ... aborting"
             call abort
          end if
          col_mass = col_mass + full_mass
       end do
       
       ! Measure relative change in mass
       do k = 1, zlevels
          if (compressible) then
             init_mass = a_vert_mass(k) + b_vert_mass(k)*col_mass
             min_mass_loc = min (min_mass_loc, sol(S_MASS,k)%data(d)%elts(id_i)/init_mass)
          else
             init_mass = sol_mean(S_MASS,k)%data(d)%elts(id_i)
             full_mass = sol(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)
             min_mass_loc = min (min_mass_loc, full_mass/init_mass)
          end if
       end do
       
       ! Check diffusion stability
       do e = 1, EDGE
          id_e = EDGE*id + e
          if (dom%mask_e%elts(id_e) >= ADJZONE) then
             d_e = dom%len%elts(id_e) ! triangle edge length
             fac = dt/d_e**(2*Laplace_order)
             beta_sclr_loc = max (beta_sclr_loc, maxval(visc_sclr) * fac)
             beta_divu_loc = max (beta_divu_loc, visc_divu * fac)
             beta_rotu_loc = max (beta_rotu_loc, visc_rotu * fac)
          end if
       end do
    end if
  end subroutine cal_min_mass

  integer function domain_load (dom)
    implicit none
    type(Domain) :: dom

    domain_load = count(abs(dom%mask_n%elts(1+1:dom%node%length)) > ADJZONE) &
         + count(abs(dom%mask_e%elts(EDGE+1:dom%midpt%length)) > ADJZONE)
  end function domain_load

  subroutine write_load_conn1 (fid)
    implicit none
    integer :: fid

    integer :: d, n_active_d

    do d = 1, size(grid)
       ! The following includes load for boundaries, but that seem just fair
       n_active_d = domain_load(grid(d))
       write(fid,'(I10, 99999(1X,I8))') n_active_d, ( &
            grid(d)%pack(AT_NODE,:)%length + grid(d)%pack(AT_EDGE,:)%length + &
            grid(d)%unpk(AT_NODE,:)%length + grid(d)%unpk(AT_EDGE,:)%length)/2
    end do
  end subroutine write_load_conn1
end module comm_mod
