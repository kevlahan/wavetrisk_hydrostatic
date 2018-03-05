module comm_mod
  use arch_mod
  use domain_mod
  implicit none
  integer, dimension(4,4)            :: shift_arr
  integer, dimension(:), allocatable :: n_active_edges, n_active_nodes
  real(8)                            :: dt_loc, sync_val
contains
  subroutine init_comm_mod
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_arch_mod()
    call init_domain_mod()
    shift_arr = reshape((/0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0/), &
         (/4, 4/))
    initialized = .True.
  end subroutine init_comm_mod

  subroutine init_comm
    integer :: k, s, d

    allocate(n_active_edges(min_level-1:max_level), n_active_nodes(min_level-1:max_level))

    n_active_edges = 0
    n_active_nodes = 0
    do d = 1, size(grid)
       do s = 1, N_BDRY
          if (.not. is_penta(grid(d), 1, s - 1)) then
             call create_comm(grid(d), 1, s)
          else
             do k = 1, 2
                call create_comm_pole(grid(d), 1, s, k - 1)
             end do
          end if
       end do
    end do
  end subroutine init_comm

  subroutine comm_nodes9 (get, set)
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
             call get(grid(src_loc), src_id, val)
             call set(grid(dest_loc), dest_id, val)
          end do
       end do
    end do
  end subroutine comm_nodes9

  subroutine create_pack_st2(dom, src, i, j, pa, e, id, e_pack, orient)
    type(Domain) dom
    integer src
    integer i
    integer j
    integer pa
    integer e
    integer id
    integer e_pack
    integer orient

    ! if PATCH_SIZE == 4 (i.e. PATCH_LEVEL == 2) and BDRY_THICKNESS == 3
    ! then some halo cells have their partner not on a neighbouring patch
    ! luckily these (very few) halo cells are not used and can be skipped:

    if (i .lt. 0 .or. j .lt. 0 .or. i .ge. PATCH_SIZE .or. j .ge. PATCH_SIZE) return

    if (e .eq. NODE) then
       call create_pack_st(dom, AT_NODE, src, i, j, pa, e, id*orient)
    else
       call create_pack_st(dom, AT_EDGE, src, i, j, pa, e_pack, orient*(EDGE*id + e))
    end if
  end subroutine create_pack_st2

  subroutine comm_nodes3(get, set)
    external get
    type(Coord) get
    external set
    integer src_loc
    integer src_glo
    integer dest_loc
    integer dest_glo
    integer i
    integer src_id
    integer dest_id

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
             src_id = grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)
             call set(grid(dest_loc), dest_id, get(grid(src_loc), src_id))
          end do
       end do
    end do
  end subroutine comm_nodes3

  subroutine create_comm_e(dom, p, s, e)
    type(Domain) dom
    integer p
    integer s
    integer e
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer b
    integer typ
    integer src
    integer rot
    integer ngb_pa
    integer orient
    integer t_last
    integer t_next
    integer s_adj
    integer shift
    integer e_pack
    integer j
    integer i
    integer i_recv
    integer j_recv
    integer i_pack
    integer j_pack

    call get_offs_Domain(dom, p, offs, dims)

    b = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side
    if (typ .lt. 1) then
       return
    end if

    src = dom%neigh(typ)

    if (src .eq. NONE) then
       return
    end if

    rot = dom%neigh_rot(typ)
    ngb_pa = dom%bdry_patch%elts(b+1)%neigh
    orient = (-1)**rot
    if (s .gt. 4) then
       t_last = typ - 4
       t_next = modulo(typ - 4, 4) + 1
       if (dom%neigh_rot(t_last) .eq. 1) then
          s_adj = s - 4
          shift = shift_arr(e+1,s_adj)
          e_pack = modulo(e + rot*(modulo(s_adj, 2) + 1), 3)
       else
          if (dom%neigh_rot(t_next) .eq. 1) then
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

    if (s .eq. NORTH) then
       do j = 1, BDRY_THICKNESS
          do i = 1, PATCH_SIZE
             call pack_idx(i - 1, j - 1, rot, s, shift, i_recv, j_recv, &
                  i_pack, j_pack)
             call create_pack_st2(dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do
    end if

    if (s .eq. EAST) then
       do i = 1, BDRY_THICKNESS
          do j = 1, PATCH_SIZE
             call pack_idx(i - 1, j - 1, rot, s, shift, i_recv, j_recv, &
                  i_pack, j_pack)
             call create_pack_st2(dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do
    end if

    if (s .eq. SOUTH) then
       do j = 1, BDRY_THICKNESS
          do i = 1, PATCH_SIZE
             call pack_idx(i - 1, j - 1, rot, s, shift, i_recv, j_recv, &
                  i_pack, j_pack)
             call create_pack_st2(dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do
    end if

    if (s .eq. WEST) then
       do i = 1, BDRY_THICKNESS
          do j = 1, PATCH_SIZE
             call pack_idx(i - 1, j - 1, rot, s, shift, i_recv, j_recv, &
                  i_pack, j_pack)
             call create_pack_st2(dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do
    end if

    if (s .eq. NORTHEAST) then
       if (dom%neigh_rot(t_last) .eq. 1) then
          do j = 1, BDRY_THICKNESS
             do i = 1, BDRY_THICKNESS - rot*(j + shift - 1)
                call pack_idx(i - 1, j - 1, rot, s_adj, shift, i_recv, &
                     j_recv, i_pack, j_pack)
                i_recv = i_recv + PATCH_SIZE
                call create_pack_st2(dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), e_pack, orient)
             end do
          end do
       else
          if (dom%neigh_rot(t_next) .eq. 1) then
             do i = 1, BDRY_THICKNESS
                do j = 1, BDRY_THICKNESS - rot*(i + shift - 1)
                   call pack_idx(i - 1, j - 1, rot, s_adj, shift, &
                        i_recv, j_recv, i_pack, j_pack)
                   j_recv = j_recv + PATCH_SIZE
                   call create_pack_st2(dom, src, i_pack, j_pack, &
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
                   call create_pack_st2(dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          end if
       end if
    end if

    if (s .eq. SOUTHEAST) then
       if (dom%neigh_rot(t_last) .eq. 1) then
          do i = 1, BDRY_THICKNESS
             do j = 1, BDRY_THICKNESS + rot*(i + shift - 1)
                i_recv = i - 1 + PATCH_SIZE
                j_recv = -j + rot*(i + shift - 1)
                i_pack = j - 1
                j_pack = i - 1
                call create_pack_st2(dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), e_pack, orient)
             end do
          end do
       else
          if (dom%neigh_rot(t_next) .eq. 1) then
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS + rot*(j + shift - 1)
                   i_recv = (PATCH_SIZE + i - 1) - rot*(j + shift - 1)
                   j_recv = -1 - (j - 1)
                   i_pack = LAST - (j - 1)
                   j_pack = LAST - (i - 1)
                   call create_pack_st2(dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          else
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS
                   i_recv = i - 1 + PATCH_SIZE
                   j_recv = -1 - (j - 1)
                   i_pack = i - 1
                   j_pack = LAST - (j - 1)
                   call create_pack_st2(dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          end if
       end if
    end if

    if (s .eq. SOUTHWEST) then
       if (dom%neigh_rot(t_last) .eq. 1) then
          do j = 1, BDRY_THICKNESS
             do i = 1, BDRY_THICKNESS - rot*(j + shift - 1)
                i_recv = -i - rot*(j + shift - 1)
                j_recv = -1 - (j - 1)
                i_pack = LAST - (j - 1)
                j_pack = i - 1
                call create_pack_st2(dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), e_pack, orient)
             end do
          end do
       else
          if (dom%neigh_rot(t_next) .eq. 1) then
             do i = 1, BDRY_THICKNESS
                do j = 1, BDRY_THICKNESS - rot*(i + shift - 1)
                   i_recv = -1 - (i - 1)
                   j_recv = -j - rot*(i + shift - 1)
                   i_pack = j - 1
                   j_pack = LAST - (i - 1)
                   call create_pack_st2(dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          else
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS - rot*(j + shift - 1)
                   i_recv = -1 - (i - 1)
                   j_recv = -1 - (j - 1)
                   i_pack = LAST - (i - 1)
                   j_pack = LAST - (j - 1)
                   call create_pack_st2(dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          end if
       end if
    end if

    if (s .eq. NORTHWEST) then
       if (dom%neigh_rot(t_last) .eq. 1) then
          do i = 1, BDRY_THICKNESS
             do j = 1, BDRY_THICKNESS + rot*(i + shift - 1)
                i_recv = -1 - (i - 1)
                j_recv = (PATCH_SIZE + j - 1) - rot*(i + shift - 1)
                i_pack = LAST - (j - 1)
                j_pack = LAST - (i - 1)
                call create_pack_st2(dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), e_pack, orient)
             end do
          end do
       else
          if (dom%neigh_rot(t_next) .eq. 1) then
             do j = 1, BDRY_THICKNESS
                do i = 1, BDRY_THICKNESS + rot*(j + shift - 1)
                   i_recv = -i + rot*(j + shift - 1)
                   j_recv = j - 1 + PATCH_SIZE
                   i_pack = j - 1
                   j_pack = i - 1
                   call create_pack_st2(dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          else
             do i = 1, BDRY_THICKNESS
                do j = 1, BDRY_THICKNESS
                   i_recv = -1 - (i - 1)
                   j_recv = j - 1 + PATCH_SIZE
                   i_pack = LAST - (i - 1)
                   j_pack = j - 1
                   call create_pack_st2(dom, src, i_pack, j_pack, &
                        ngb_pa, e, idx(i_recv, j_recv, offs, dims), &
                        e_pack, orient)
                end do
             end do
          end if
       end if
    end if
  end subroutine create_comm_e

  function rot_direction (dom, typ)
    integer      :: rot_direction
    type(Domain) :: dom
    integer      :: typ

    integer :: t_last, t_next

    if (typ .le. 4) then
       rot_direction = modulo(typ, 2)
       return
    else
       t_last = typ - 4
       if (dom%neigh_rot(t_last) .eq. 1) then
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
    integer               :: i, dest_glo, dest_loc, src_glo, src_loc
    integer, dimension(4) :: st

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%send_conn(dest_glo+1)%length, 4
             st = grid(src_loc)%send_conn(dest_glo+1)%elts(i:i + 3)
             call unpack_comm_struct(grid(dest_loc), src_glo, st(1), &
                  st(2), st(3), st(4))
          end do
          grid(src_loc)%send_conn(dest_glo+1)%length = 0
       end do
    end do
  end subroutine comm_communication

  type(Coord) function get_coord(dom, id)
    type(Domain) dom
    integer id

    get_coord = dom%node%elts(id+1)
  end function get_coord

  subroutine create_comm(dom, p, s)
    type(Domain) dom
    integer p
    integer s
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer b
    integer typ
    integer src
    integer rot
    integer pa
    integer e

    call get_offs_Domain(dom, p, offs, dims)

    b = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side

    if (typ .lt. 1) then
       return
    end if

    src = dom%neigh(typ)
    rot = dom%neigh_rot(typ)
    pa = dom%bdry_patch%elts(b+1)%neigh

    do e = 1, EDGE + 1
       call create_comm_e(dom, p, s, e - 1)
    end do

    if (rot .eq. 1 .and. typ .eq. WEST) then
       if (is_penta(dom, p, SOUTHWEST - 1)) then
          call create_pack_st(dom, AT_EDGE, src, LAST, LAST, pa, RT, &
               nidx(LAST_BDRY, LAST_BDRY, SOUTHWEST, offs, dims)*EDGE)
       end if
    end if

    if (rot .eq. 1 .and. typ .eq. SOUTH) then
       if (is_penta(dom, p, SOUTHWEST - 1)) then
          call create_pack_st(dom, AT_EDGE, src, LAST, LAST, pa, UP, &
               -nidx(LAST_BDRY, LAST_BDRY, SOUTHWEST, offs, dims)*EDGE)
       end if
    end if

    if (rot .eq. 0 .and. typ .eq. EAST) then
       if (is_penta(dom, p, SOUTHEAST - 1)) then
          call create_pack_st(dom, AT_NODE, src, 1, 0, pa, NODE, nidx(0, &
               LAST_BDRY, SOUTHEAST, offs, dims))
       end if
    end if

    if (rot .eq. 0 .and. typ .eq. NORTH) then
       if (is_penta(dom, p, NORTHWEST - 1)) then
          call create_pack_st(dom, AT_NODE, src, 0, 1, pa, NODE, &
               nidx(LAST_BDRY, 0, NORTHWEST, offs, dims))
       end if
    end if

    if (rot .eq. 1 .and. typ .eq. EAST) then
       if (is_penta(dom, p, NORTHEAST - 1)) then
          call create_pack_st(dom, AT_NODE, src, 0, 1, pa, NODE, nidx(0, &
               1, NORTHEAST, offs, dims))
          call create_pack_st(dom, AT_NODE, src, 1, 1, pa, NODE, &
               -nidx(1, 0, NORTHEAST, offs, dims))
          call create_pack_st(dom, AT_EDGE, src, 0, 1, pa, RT, nidx(0, &
               1, NORTHEAST, offs, dims)*EDGE + RT)
          call create_pack_st(dom, AT_EDGE, src, 0, 0, pa, DG, -(nidx(0, &
               0, NORTHEAST, offs, dims)*EDGE + RT))
          call create_pack_st(dom, AT_EDGE, src, 0, 0, pa, UP, nidx(0, &
               0, NORTHEAST, offs, dims)*EDGE + UP)
       end if
    end if

    if (rot .eq. 1 .and. typ .eq. NORTH) then
       if (is_penta(dom, p, NORTHEAST - 1)) then
          call create_pack_st(dom, AT_NODE, src, 1, 0, pa, NODE, nidx(1, &
               0, NORTHEAST, offs, dims))
          call create_pack_st(dom, AT_NODE, src, 1, 1, pa, NODE, &
               -nidx(0, 1, NORTHEAST, offs, dims))
          call create_pack_st(dom, AT_EDGE, src, 1, 0, pa, UP, -(nidx(0, &
               1, NORTHEAST, offs, dims)*EDGE + RT))
          call create_pack_st(dom, AT_EDGE, src, 0, 0, pa, DG, -(nidx(0, &
               0, NORTHEAST, offs, dims)*EDGE + UP))
          call create_pack_st(dom, AT_EDGE, src, 0, 0, pa, RT, nidx(0, &
               0, NORTHEAST, offs, dims)*EDGE + RT)
       end if
    end if
  end subroutine create_comm

  subroutine create_comm_pole(dom, p, s, i)
    type(Domain) dom
    integer p
    integer s
    integer i
    integer b
    integer typ
    integer lev
    integer src
    integer pa
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer s_side
    integer, dimension(2) :: ij_node
    integer, dimension(2) :: ij_send

    !  side s of patch p is a pole and i-th connection over this pole is available
    b = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side

    if (typ .lt. 1) then
       return
    end if

    if (.not. dom%neigh(typ) .eq. POLE) then
       return
    end if

    lev = dom%patch%elts(p+1)%level
    src = dom%neigh_over_pole(i+1)
    pa = dom%neigh_pa_over_pole%elts(i+2*lev+1)

    call get_offs_Domain(dom, p, offs, dims)

    if (s .eq. NORTHWEST) then
       s_side = NORTH
    else
       if (s .eq. SOUTHEAST) then
          s_side = EAST
       end if
    end if

    if (i .eq. 0) then
       ij_node = (/0, 1/)
       if (s .eq. SOUTHEAST) then
          ij_node = (/ij_node(2), ij_node(1)/)
       end if
       call create_pack_st(dom, AT_EDGE, src, 0, LAST, pa, DG, &
            nidx(ij_node(1), ij_node(2), s_side, offs, dims)*EDGE + &
            2*s_side - 2)
       ij_node = (/LAST_BDRY, 1/)
       ij_send = (/1, LAST/)
       if (s .eq. SOUTHEAST) then
          ij_node = (/ij_node(2), ij_node(1)/)
          ij_send = (/ij_send(2), ij_send(1)/)
       end if
       call create_pack_st(dom, AT_NODE, src, ij_send(1), ij_send(2), pa, &
            NODE, -nidx(ij_node(1), ij_node(2), s, offs, dims))
    end if

    if (i .eq. 1) then
       ij_node = (/LAST_BDRY, 0/)
       ij_send = (/0, LAST/)
       if (s .eq. SOUTHEAST) then
          ij_node = (/ij_node(2), ij_node(1)/)
          ij_send = (/ij_send(2), ij_send(1)/)
       end if
       call create_pack_st(dom, AT_NODE, src, ij_send(1), ij_send(2), pa, &
            NODE, nidx(ij_node(1), ij_node(2), s, offs, dims))
    end if

    ij_node = (/-i + 1, 1/)
    ij_send = (/0, LAST/)

    if (s .eq. SOUTHEAST) then
       ij_node = (/ij_node(2), ij_node(1)/)
       ij_send = (/ij_send(2), ij_send(1)/)
    end if

    call create_pack_st(dom, AT_NODE, src, ij_send(1), ij_send(2), pa, &
         NODE, (-1)**i*nidx(ij_node(1), ij_node(2), s_side, offs, dims))

    call create_pack_st(dom, AT_EDGE, src, ij_send(1), ij_send(2), pa, UP &
         - 2*s_side + 2, (-1)**i*(nidx(0, 0, s_side, offs, dims)*EDGE + DG &
         + i*(-2*s_side + 3)))

    if (i .eq. 0) then
       return
    end if

    if (s .eq. NORTHWEST) then
       s_side = WEST
    else
       if (s .eq. SOUTHEAST) then
          s_side = SOUTH
       end if
    end if

    ij_node = (/LAST_BDRY - 1, LAST/)
    ij_send = (/1, LAST/)

    if (s .eq. SOUTHEAST) then
       ij_node = (/ij_node(2), ij_node(1)/)
       ij_send = (/ij_send(2), ij_send(1)/)
    end if

    call create_pack_st(dom, AT_NODE, src, ij_send(1), ij_send(2), pa, &
         NODE, nidx(ij_node(1), ij_node(2), s_side, offs, dims))

    ij_node = (/LAST, LAST_BDRY/)

    if (s .eq. NORTHWEST) then
       ij_node = (/ij_node(2), ij_node(1)/)
    end if

    call create_pack_st(dom, AT_EDGE, src, LAST*(-s_side + 4), &
         LAST*(s_side - 3), pa, DG, nidx(ij_node(1), ij_node(2), s_side, &
         offs, dims)*EDGE + 2*s_side - 6)
  end subroutine create_comm_pole

  subroutine get_areas(dom, id, val)
    real(8), dimension(7), intent(out) :: val
    type(Domain) dom
    integer id
    real(8), dimension(7) :: area

    area = 0.0
    area(1:4) = dom%overl_areas%elts(id+1)%a
    area(5:6) = dom%overl_areas%elts(id+1)%split
    area(7) = dom%areas%elts(id+1)%hex_inv
    val = area
    return
  end subroutine get_areas

  subroutine comm_masks()
    integer src_loc
    integer src_glo
    integer dest_loc
    integer dest_glo
    integer i, k
    integer src_id
    integer dest_id

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
             src_id = grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)
             grid(dest_loc)%mask_n%elts(abs(dest_id)+1) = grid(src_loc)%mask_n%elts(abs(src_id)+1) 
          end do
          do i = 1, grid(src_loc)%pack(AT_EDGE,dest_glo+1)%length
             src_id = grid(src_loc)%pack(AT_EDGE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_EDGE,src_glo+1)%elts(i)
             grid(dest_loc)%mask_e%elts(abs(dest_id)+1) = grid(src_loc)%mask_e%elts(abs(src_id)+1) 
          end do
       end do
    end do
  end subroutine comm_masks

  subroutine comm_edges(get, set)
    external get
    real(8) get
    external set
    integer src_loc
    integer src_glo
    integer dest_loc
    integer dest_glo
    integer i
    integer src_id
    integer dest_id

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_EDGE,dest_glo+1)%length
             src_id = grid(src_loc)%pack(AT_EDGE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_EDGE,src_glo+1)%elts(i)
             call set(grid(dest_loc), dest_id, get(grid(src_loc), src_id))
          end do
       end do
    end do
  end subroutine comm_edges

  subroutine unpack_comm_struct(dom, src, i, j, p, e)
    type(Domain) dom
    integer src
    integer i
    integer j
    integer p
    integer e
    integer offs

    if (e .eq. NODE) then
       offs = dom%patch%elts(p+1)%elts_start
       call append(dom%pack(AT_NODE,src+1), idx__fast(i, j, offs))
    end if

    if (.not. e .eq. NODE) then
       offs = dom%patch%elts(p+1)%elts_start
       call append(dom%pack(AT_EDGE,src+1), idx__fast(i, j, offs)*EDGE + e)
    end if
  end subroutine unpack_comm_struct

  subroutine pack_idx(i, j, rot, s, shift, i_recv, j_recv, i_pack, j_pack)
    integer i
    integer j
    integer rot
    integer s
    integer shift
    integer i_recv
    integer j_recv
    integer i_pack
    integer j_pack

    if (s .eq. NORTH) then
       i_recv = i + rot*(j + shift)
       j_recv = j + PATCH_SIZE
       if (rot .eq. 1) then
          i_pack = j
          j_pack = LAST - i
       else
          i_pack = i
          j_pack = j
       end if
    end if

    if (s .eq. EAST) then
       i_recv = i + PATCH_SIZE
       j_recv = j + rot*(i + shift)
       if (rot .eq. 1) then
          i_pack = LAST - j
          j_pack = i
       else
          i_pack = i
          j_pack = j
       end if
    end if

    if (s .eq. SOUTH) then
       i_recv = i - rot*(j + shift)
       j_recv = -1 - j
       if (rot .eq. 1) then
          i_pack = LAST - j
          j_pack = LAST - i
       else
          i_pack = i
          j_pack = LAST - j
       end if
    end if

    if (s .eq. WEST) then
       i_recv = -1 - i
       j_recv = j - rot*(i + shift)
       if (rot .eq. 1) then
          i_pack = LAST - j
          j_pack = LAST - i
       else
          i_pack = LAST - i
          j_pack = j
       end if
    end if
  end subroutine pack_idx

  subroutine update_comm(dom)
    type(Domain) dom
    integer unused_elements
    integer src_glo
    integer ii
    integer ngb_pa
    integer c
    integer p_par
    integer s
    integer n_par
    integer typ
    integer p_chd
    integer lev
    integer n_chd

    do src_glo = 1, N_GLO_DOMAIN
       unused_elements = 0
       do ii = 1, dom%recv_pa(src_glo)%length, 4
          ngb_pa = dom%recv_pa(src_glo)%elts(ii)
          c = dom%recv_pa(src_glo)%elts(ii+1)
          p_par = dom%recv_pa(src_glo)%elts(ii+2)
          s = dom%recv_pa(src_glo)%elts(ii+3)
          n_par = -dom%patch%elts(p_par+1)%neigh(s+1)
          typ = dom%bdry_patch%elts(n_par+1)%side
          if (typ .lt. 1) then
             return
          end if
          if (dom%neigh(typ) .eq. POLE) then
             p_chd = dom%patch%elts(p_par+1)%children(s-3)
          else
             p_chd = dom%patch%elts(p_par+1)%children(c+1)
          end if
          if (p_chd .eq. 0) then
             dom%recv_pa(src_glo)%elts(unused_elements + &
                  1:unused_elements + 4) = (/ngb_pa, c, p_par, s/)
             unused_elements = unused_elements + 4
             cycle
          end if
          if (dom%neigh(typ) .eq. POLE) then
             lev = dom%patch%elts(p_chd+1)%level
             if (2*lev + 2 .gt. dom%neigh_pa_over_pole%length) then
                call extend(dom%neigh_pa_over_pole, 2, 0)
             end if
             dom%neigh_pa_over_pole%elts((1-c)+2*lev+1) = ngb_pa
             call create_comm_pole(dom, p_chd, s + 1, 1-c)
             cycle
          end if
          n_chd = -dom%patch%elts(p_chd+1)%neigh(s+1)
          dom%bdry_patch%elts(n_chd+1)%neigh = ngb_pa
          if (.not. is_penta(dom, p_chd, s)) then
             call create_comm(dom, p_chd, s + 1)
          end if
       end do
       dom%recv_pa(src_glo)%length = unused_elements
    end do
  end subroutine update_comm

  subroutine area_post_comm(dom, p, c, offs, dims, zlev)
    type(Domain) dom
    integer p
    integer c
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer id

    if (c .eq. IPLUSJMINUS) then
       id = idx(PATCH_SIZE, -1, offs, dims)
       dom%overl_areas%elts(id+1)%a = 0.0
    end if
    if (c .eq. IMINUSJPLUS) then
       id = idx(-1, PATCH_SIZE, offs, dims)
       dom%overl_areas%elts(id+1)%a = 0.0
    end if
  end subroutine area_post_comm

  subroutine comm_patch_conn()
    integer src
    integer unused_elements
    integer src_glo
    integer i
    integer b
    integer c
    integer p_chd, l_par
    integer s
    integer typ
    integer dest
    logical is_pole
    integer ngh_pa
    integer rot
    integer rot_shift

    do src = 1, size(grid)
       unused_elements = 0
       src_glo = glo_id(rank+1,src)
       do i = 1, grid(src)%send_pa_all%length, 4
          b = grid(src)%send_pa_all%elts(i)
          c = grid(src)%send_pa_all%elts(i+1)
          p_chd = grid(src)%send_pa_all%elts(i+2)
          s = grid(src)%send_pa_all%elts(i+3)
          typ = grid(src)%bdry_patch%elts(b+1)%side

          if (typ .lt. 1) then
             return
          end if

          dest = grid(src)%neigh(typ)
          is_pole = dest .eq. POLE

          if (is_pole) then
             dest = grid(src)%neigh_over_pole(c+1)
             l_par = grid(src)%patch%elts(p_chd+1)%level - 1
             if (grid(src)%neigh_pa_over_pole%length .lt. l_par*2 + c + 1) then
                ngh_pa = 0
             else
                ngh_pa = grid(src)%neigh_pa_over_pole%elts(l_par*2 + c + 1)
             end if
          else
             ngh_pa = grid(src)%bdry_patch%elts(b+1)%neigh
          end if

          if (ngh_pa .eq. 0) then
             grid(src)%send_pa_all%elts(unused_elements + &
                  1:unused_elements + 4) = (/b, c, p_chd, s/)
             unused_elements = unused_elements + 4
             cycle
          else 
             if (dest .eq. NONE) cycle
          end if

          ! cycle if dest is not on rank
          if (.not. owner(dest+1) .eq. rank) cycle

          dest = loc_id(dest+1)
          rot = grid(src)%neigh_rot(typ)
          rot_shift = (rot_direction(grid(src), typ)*2 - 1)*rot

          call append(grid(dest+1)%recv_pa(src_glo+1), p_chd)

          if (is_pole) then
             call append(grid(dest+1)%recv_pa(src_glo+1), c)
          else
             call append(grid(dest+1)%recv_pa(src_glo+1), modulo(c + &
                  rot_shift, 4))
          end if

          call append(grid(dest+1)%recv_pa(src_glo+1), ngh_pa)

          if (is_pole) then
             call append(grid(dest+1)%recv_pa(src_glo+1), s)
          else
             call append(grid(dest+1)%recv_pa(src_glo+1), modulo(rot_shift &
                  + s + 2, 4) + 4*(s/4))
          end if
       end do
       grid(src)%send_pa_all%length = unused_elements
    end do
  end subroutine comm_patch_conn

  subroutine create_pack_st(dom, unpk_pos, src, i, j, pa, e, id)
    type(Domain) dom
    integer unpk_pos
    integer src
    integer i
    integer j
    integer pa
    integer e
    integer id

    call append(dom%send_conn(src+1), i)
    call append(dom%send_conn(src+1), j)
    call append(dom%send_conn(src+1), pa)
    call append(dom%send_conn(src+1), e)
    call append(dom%unpk(unpk_pos,src+1), id)
  end subroutine create_pack_st

  subroutine set_coord(dom, id, val)
    type(Domain) dom
    integer id
    type(Coord) val

    dom%node%elts(abs(id) + 1) = val
  end subroutine set_coord

  subroutine cp_bdry_inside(field)
    type(Float_Field) field
    integer src_loc, src_glo
    integer dest_loc, dest_glo
    integer src_id, dest_id
    integer i, pos

    pos = field%pos
    
    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
             src_id = grid(src_loc)%pack(pos,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)
             
             field%data(dest_loc)%elts(abs(dest_id)+1) = field%data(src_loc)%elts(src_id+1)
             
             if (dest_id .lt. 0 .and. pos .eq. AT_EDGE) &
                  field%data(dest_loc)%elts(abs(dest_id)+1) = &
                  -field%data(dest_loc)%elts(abs(dest_id)+1)
          end do
       end do
    end do
  end subroutine cp_bdry_inside

  subroutine cp_bdry_inside_vector(field)
    type(Float_Field), dimension(:) :: field
    integer src_loc, src_glo
    integer dest_loc, dest_glo
    integer src_id, dest_id
    integer i, pos, i1, sz

    ! Find size of field
    sz = size(field)

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i1 = 1, sz
             pos = field(i1)%pos
             do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
                src_id = grid(src_loc)%pack(pos,dest_glo+1)%elts(i)
                dest_id = grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)

                field(i1)%data(dest_loc)%elts(abs(dest_id)+1) = field(i1)%data(src_loc)%elts(src_id+1)

                if (dest_id .lt. 0 .and. pos .eq. AT_EDGE) &
                     field(i1)%data(dest_loc)%elts(abs(dest_id)+1) = &
                     -field(i1)%data(dest_loc)%elts(abs(dest_id)+1)
             end do
          end do
       end do
    end do
  end subroutine cp_bdry_inside_vector

   subroutine cp_bdry_inside_array(field)
    type(Float_Field), dimension(:,:) :: field
    integer :: src_loc, src_glo
    integer :: dest_loc, dest_glo
    integer :: src_id, dest_id
    integer :: i, pos, i1, i2
    integer, dimension(2) :: sz

    ! Find size of field
    sz = shape(field)

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i2 = 1, sz(2)
             do i1 = 1, sz(1)
                pos = field(i1,i2)%pos
                do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
                   src_id = grid(src_loc)%pack(pos,dest_glo+1)%elts(i)
                   dest_id = grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)

                   field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1) = field(i1,i2)%data(src_loc)%elts(src_id+1)

                   if (dest_id .lt. 0 .and. pos .eq. AT_EDGE) &
                        field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1) = &
                        -field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1)
                end do
             end do
          end do
       end do
    end do
  end subroutine cp_bdry_inside_array

  subroutine comm_nodes(get, set)
    external get
    real(8) get
    external set
    integer src_loc
    integer src_glo
    integer dest_loc
    integer dest_glo
    integer i
    integer src_id
    integer dest_id

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1,src_loc)
       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1,dest_loc)
          do i = 1, grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
             src_id = grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)
             dest_id = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)
             call set(grid(dest_loc), dest_id, get(grid(src_loc), src_id))
          end do
       end do
    end do
  end subroutine comm_nodes

  subroutine set_areas (dom, id, val)
    type(Domain)          :: dom
    integer               :: id
    real(8), dimension(7) :: val
    
    real(8), dimension(4) :: area

    area = val(1:4)
    if (id .lt. 0) then
       area = (/area(2), area(1), area(4), area(3)/)
    end if

    dom%overl_areas%elts(abs(id) + 1)%a = area
    dom%overl_areas%elts(abs(id) + 1)%split = val(5:6)
    dom%areas%elts(abs(id) + 1)%hex_inv = val(7)
  end subroutine set_areas

  subroutine min_dt (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: d, e, id, k, l
    real(8) :: C_visc, csq, d_e, v_e, viscosity

    C_visc = 0.25_8
    
    id = idx(i, j, offs, dims)
    d  = dom%id + 1
    l  = dom%level%elts(id+1)

    if (dom%mask_n%elts(id+1) .ge. ADJZONE) then
       n_active_nodes(l) = n_active_nodes(l) + 1

       do e = 1, EDGE
          if (dom%mask_e%elts(EDGE*id+e) .ge. ADJZONE) then
             n_active_edges(l) = n_active_edges(l) + 1
             if (adapt_dt) then
                d_e = dom%len%elts(EDGE*id+e) ! Triangle edge length
                do k = 1, zlevels
                   v_e = abs(sol(S_VELO,k)%data(d)%elts(EDGE*id+e))
                   if (v_e.ne.0.0_8) dt_loc =  min(dt_loc, cfl_num*d_e/(v_e+wave_speed))
                end do
                dt_loc = min (dt_loc, C_visc*d_e**2/viscosity)
             end if
          end if
       end do
    end if
  end subroutine min_dt

  function domain_load (dom)
    integer      :: domain_load
    type(Domain) :: dom

    domain_load = count(abs(dom%mask_n%elts(1+1:dom%node%length)) .gt. ADJZONE) + &
         count(abs(dom%mask_e%elts(EDGE+1:dom%midpt%length)) .gt. ADJZONE)
  end function domain_load

  subroutine write_load_conn1 (fid)
    integer :: fid

    integer :: d, n_active_d

    do d = 1, size(grid)
       ! the following adds load for boundaries, but that seem just fair
       n_active_d = domain_load(grid(d))
       write(fid,'(I10, 99999(1X,I8))') n_active_d, ( &
            grid(d)%pack(AT_NODE,:)%length + grid(d)%pack(AT_EDGE,:)%length + &
            grid(d)%unpk(AT_NODE,:)%length + grid(d)%unpk(AT_EDGE,:)%length)/2
    end do
  end subroutine write_load_conn1

end module comm_mod
