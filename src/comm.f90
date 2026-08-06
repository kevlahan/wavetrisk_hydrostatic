module comm_mod
  use kind_mod,   only : dp

  use shared_mod, only : ADJZONE, AT_EDGE, AT_NODE, Coord, EDGE, N_GLO_DOMAIN, NODE, NONE, &
       EAST, NORTH, NORTHEAST, NORTHWEST, SOUTH, SOUTHEAST, SOUTHWEST, WEST, POLE, RT, DG, UP, N_BDRY, BDRY_THICKNESS

  use arch_mod,   only : glo_id, init_arch_mod, loc_id, owner, rank
  use domain_mod, only : get_offs_domain, grid, Domain, Float_Field, idx, idx__fast, init_domain_mod, nidx, is_penta
  use dyn_arrays, only : append, extend
  use patch_mod,  only : LAST, LAST_BDRY, PATCH_SIZE

  implicit none

  private
  public :: comm_communication, comm_edges, comm_nodes, comm_nodes3, comm_nodes9, comm_patch_conn
  public :: init_comm, init_comm_mod, update_comm
  public :: coord_get, coord_set, get_coord, set_coord, get9, set9, unpack, unpack_comm_struct
  public :: domain_load, rot_direction, sync_val
  public :: comm_masks, cp_bdry_inside

  integer, dimension(4,4) :: shift_arr
  real(dp)                :: sync_val

  abstract interface

     type(Coord) function coord_get (dom, id) 
       use shared_mod, only : Coord
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(in) :: dom
       integer,      intent(in) :: id
     end function coord_get

     subroutine coord_set (dom, id, val)
       use shared_mod, only : Coord
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: id
       type(Coord),  intent(in)    :: val
     end subroutine coord_set

     real(dp) function real_get (dom, id)
       use kind_mod,   only : dp
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(in) :: dom
       integer,      intent(in) :: id
     end function real_get

     subroutine real_set (dom, id, val)
       use kind_mod,   only : dp
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: id
       real(dp),     intent(in)    :: val
     end subroutine real_set

     subroutine get9 (dom, id, val)
       use kind_mod,   only : dp
       use domain_mod, only : Domain
       implicit none

       type(Domain),           intent(in)  :: dom
       integer,                intent(in)  :: id
       real(dp), dimension(7), intent(out) :: val
     end subroutine get9

     subroutine set9 (dom, id, val)
       use kind_mod,   only : dp
       use domain_mod, only : Domain
       implicit none

       type(Domain),           intent(inout) :: dom
       integer,                intent(in)    :: id
       real(dp), dimension(7), intent(in)    :: val
     end subroutine set9

     subroutine unpack (dom, src, i, j, p, e)
       use domain_mod, only : Domain
       implicit none
       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: e, i, j, p, src
     end subroutine unpack

  end interface

  interface cp_bdry_inside
     procedure :: cp_bdry_inside_0, cp_bdry_inside_1, cp_bdry_inside_2
  end interface cp_bdry_inside

contains

  
  subroutine init_comm_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_arch_mod
    call init_domain_mod

    shift_arr = reshape ([ 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0 ], [ 4, 4 ])

    initialized = .true.
  end subroutine init_comm_mod


  subroutine init_comm
    implicit none

    integer :: d, k, s

    do d = 1, size(grid)
       do s = 1, N_BDRY
          if (is_penta (grid(d), 1, s-1)) then
             do k = 0, 1
                call create_comm_pole (grid(d), 1, s, k)
             end do
          else
             call create_comm (grid(d), 1, s)
          end if
       end do
    end do
  end subroutine init_comm

  
  subroutine comm_nodes (get, set)
    implicit none

    procedure(real_get) :: get
    procedure(real_set) :: set

    integer :: i
    integer :: dest_id, dest_glo, dest_loc
    integer :: src_id, src_glo, src_loc
    integer :: n_pack, n_unpk

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          n_pack = grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
          n_unpk = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%length

          if (n_pack /= n_unpk) then
             error stop "comm_nodes: pack/unpack length mismatch"
          end if

          do i = 1, n_pack
             src_id = &
                  grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)

             dest_id = &
                  grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)

             call set ( &
                  grid(dest_loc), dest_id, &
                  get(grid(src_loc), src_id))
          end do
       end do
    end do
  end subroutine comm_nodes

  
  subroutine comm_nodes3 (get, set)
    implicit none

    procedure(coord_get) :: get
    procedure(coord_set) :: set

    integer :: i
    integer :: dest_id, dest_glo, dest_loc
    integer :: src_id, src_glo, src_loc
    integer :: n_pack, n_unpk
    type(Coord) :: value

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          n_pack = grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
          n_unpk = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%length

          if (n_pack /= n_unpk) then
             error stop "comm_nodes3: pack/unpack length mismatch"
          end if

          do i = 1, n_pack
             src_id = &
                  grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)

             dest_id = &
                  grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)

             value = get(grid(src_loc), src_id)
             call set(grid(dest_loc), dest_id, value)
          end do
       end do
    end do
  end subroutine comm_nodes3

  
  subroutine comm_nodes9 (get, set)
    implicit none

    procedure(get9) :: get
    procedure(set9) :: set

    integer, parameter :: NVAL = 7

    integer :: i
    integer :: dest_id, dest_glo, dest_loc
    integer :: src_id, src_glo, src_loc
    integer :: n_pack, n_unpk

    real(dp) :: val(NVAL)

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          n_pack = grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
          n_unpk = grid(dest_loc)%unpk(AT_NODE,src_glo+1)%length

          if (n_pack /= n_unpk) then
             error stop "comm_nodes9: pack/unpack length mismatch"
          end if

          do i = 1, n_pack
             src_id = &
                  grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)

             dest_id = &
                  grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)

             call get(grid(src_loc), src_id, val)
             call set(grid(dest_loc), dest_id, val)
          end do
       end do
    end do
  end subroutine comm_nodes9

  
  function get_coord (dom, id) result(val)
    type(Domain), intent(in) :: dom
    integer,      intent(in) :: id
    type(Coord)              :: val

    val = dom%node%elts(id+1)
  end function get_coord

  
  subroutine set_coord (dom, id, val)
    implicit none
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: id
    type(Coord),  intent(in)    :: val

    dom%node%elts(abs(id)+1) = val
  end subroutine set_coord
  

  subroutine create_pack_st2 (dom, src, i, j, pa, e, id, e_pack, orient)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: src, i, j, pa
    integer,      intent(in)    :: e, id, e_pack, orient

    ! For PATCH_SIZE == 4 and BDRY_THICKNESS == 3, a few halo
    ! cells have no partner on an adjacent patch. They are unused.
    if (i < 0 .or. i >= PATCH_SIZE) return
    if (j < 0 .or. j >= PATCH_SIZE) return

    if (e == NODE) then
       call create_pack_st (dom, AT_NODE, src, i, j, pa, e, id*orient)
    else
       call create_pack_st( &
            dom, AT_EDGE, src, i, j, pa, e_pack, &
            orient*(EDGE*id + e))
    end if
  end subroutine create_pack_st2

  
  subroutine create_comm_e (dom, p, s, e)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, s, e

    integer :: b, e_pack
    integer :: i, i_pack, i_recv
    integer :: j, j_pack, j_recv
    integer :: ngb_pa
    integer :: orient, rot
    integer :: s_adj, shift
    integer :: src
    integer :: t_last, t_next, typ

    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    logical :: last_rotated, next_rotated

    call get_offs_Domain(dom, p, offs, dims)

    ! Create_comm_e is meaningful only when side s refers to a
    ! boundary patch, represented by a negative neighbour index.
    if (dom%patch%elts(p+1)%neigh(s) >= 0) then
       error stop "create_comm_e: side is not a boundary patch"
    end if

    b   = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side

    if (typ < 1) return

    src = dom%neigh(typ)
    if (src == NONE) return

    rot    = dom%neigh_rot(typ)
    ngb_pa = dom%bdry_patch%elts(b+1)%neigh

    ! Equivalent to (-1)**rot, but explicit about rotation parity.
    if (modulo(rot,2) == 0) then
       orient = 1
    else
       orient = -1
    end if

    s_adj  = 0
    shift  = 0
    e_pack = e

    t_last = 0
    t_next = 0

    last_rotated = .false.
    next_rotated = .false.

    if (s > WEST) then
       t_last = typ - 4
       t_next = modulo(typ-4,4) + 1

       last_rotated = dom%neigh_rot(t_last) == 1
       next_rotated = dom%neigh_rot(t_next) == 1

       if (last_rotated) then
          s_adj  = s - 4
          shift  = shift_arr(e+1,s_adj)
          e_pack = modulo( &
               e + rot*(modulo(s_adj,2)+1), EDGE)

       else if (next_rotated) then
          s_adj  = modulo(s-4,4) + 1
          shift  = shift_arr(e+1,s_adj)
          e_pack = modulo( &
               e + rot*(modulo(s_adj,2)+1), EDGE)
       end if

    else
       shift  = shift_arr(e+1,s)
       e_pack = modulo(e + rot*(modulo(s,2)+1), EDGE)
    end if

    select case (s)

    case (NORTH)
       do j = 1, BDRY_THICKNESS
          do i = 1, PATCH_SIZE
             call pack_idx( &
                  i-1, j-1, rot, s, shift, &
                  i_recv, j_recv, i_pack, j_pack)

             call create_pack_st2( &
                  dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do

    case (EAST)
       do i = 1, BDRY_THICKNESS
          do j = 1, PATCH_SIZE
             call pack_idx( &
                  i-1, j-1, rot, s, shift, &
                  i_recv, j_recv, i_pack, j_pack)

             call create_pack_st2( &
                  dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do

    case (SOUTH)
       do j = 1, BDRY_THICKNESS
          do i = 1, PATCH_SIZE
             call pack_idx( &
                  i-1, j-1, rot, s, shift, &
                  i_recv, j_recv, i_pack, j_pack)

             call create_pack_st2( &
                  dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do

    case (WEST)
       do i = 1, BDRY_THICKNESS
          do j = 1, PATCH_SIZE
             call pack_idx( &
                  i-1, j-1, rot, s, shift, &
                  i_recv, j_recv, i_pack, j_pack)

             call create_pack_st2( &
                  dom, src, i_pack, j_pack, ngb_pa, e, &
                  idx(i_recv, j_recv, offs, dims), e_pack, orient)
          end do
       end do

    case (NORTHEAST)
       if (last_rotated) then
          do j = 1, BDRY_THICKNESS
             do i = 1, &
                  BDRY_THICKNESS - rot*(j+shift-1)

                call pack_idx( &
                     i-1, j-1, rot, s_adj, shift, &
                     i_recv, j_recv, i_pack, j_pack)

                i_recv = i_recv + PATCH_SIZE

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do

       else if (next_rotated) then
          do i = 1, BDRY_THICKNESS
             do j = 1, &
                  BDRY_THICKNESS - rot*(i+shift-1)

                call pack_idx( &
                     i-1, j-1, rot, s_adj, shift, &
                     i_recv, j_recv, i_pack, j_pack)

                j_recv = j_recv + PATCH_SIZE

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do

       else
          do j = 1, BDRY_THICKNESS
             do i = 1, BDRY_THICKNESS
                i_recv = i-1 + PATCH_SIZE
                j_recv = j-1 + PATCH_SIZE
                i_pack = i-1
                j_pack = j-1

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do
       end if

    case (SOUTHEAST)
       if (last_rotated) then
          do i = 1, BDRY_THICKNESS
             do j = 1, &
                  BDRY_THICKNESS + rot*(i+shift-1)

                i_recv = i-1 + PATCH_SIZE
                j_recv = -j + rot*(i+shift-1)
                i_pack = j-1
                j_pack = i-1

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do

       else if (next_rotated) then
          do j = 1, BDRY_THICKNESS
             do i = 1, &
                  BDRY_THICKNESS + rot*(j+shift-1)

                i_recv = PATCH_SIZE+i-1 - rot*(j+shift-1)
                j_recv = -j
                i_pack = LAST-(j-1)
                j_pack = LAST-(i-1)

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do

       else
          do j = 1, BDRY_THICKNESS
             do i = 1, BDRY_THICKNESS
                i_recv = i-1 + PATCH_SIZE
                j_recv = -j
                i_pack = i-1
                j_pack = LAST-(j-1)

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do
       end if

    case (SOUTHWEST)
       if (last_rotated) then
          do j = 1, BDRY_THICKNESS
             do i = 1, &
                  BDRY_THICKNESS - rot*(j+shift-1)

                i_recv = -i - rot*(j+shift-1)
                j_recv = -j
                i_pack = LAST-(j-1)
                j_pack = i-1

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do

       else if (next_rotated) then
          do i = 1, BDRY_THICKNESS
             do j = 1, &
                  BDRY_THICKNESS - rot*(i+shift-1)

                i_recv = -i
                j_recv = -j - rot*(i+shift-1)
                i_pack = j-1
                j_pack = LAST-(i-1)

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do

       else
          do j = 1, BDRY_THICKNESS
             do i = 1, &
                  BDRY_THICKNESS - rot*(j+shift-1)

                i_recv = -i
                j_recv = -j
                i_pack = LAST-(i-1)
                j_pack = LAST-(j-1)

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do
       end if

    case (NORTHWEST)
       if (last_rotated) then
          do i = 1, BDRY_THICKNESS
             do j = 1, &
                  BDRY_THICKNESS + rot*(i+shift-1)

                i_recv = -i
                j_recv = PATCH_SIZE+j-1 - rot*(i+shift-1)
                i_pack = LAST-(j-1)
                j_pack = LAST-(i-1)

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do

       else if (next_rotated) then
          do j = 1, BDRY_THICKNESS
             do i = 1, &
                  BDRY_THICKNESS + rot*(j+shift-1)

                i_recv = -i + rot*(j+shift-1)
                j_recv = j-1 + PATCH_SIZE
                i_pack = j-1
                j_pack = i-1

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do

       else
          do i = 1, BDRY_THICKNESS
             do j = 1, BDRY_THICKNESS
                i_recv = -i
                j_recv = j-1 + PATCH_SIZE
                i_pack = LAST-(i-1)
                j_pack = j-1

                call create_pack_st2( &
                     dom, src, i_pack, j_pack, ngb_pa, e, &
                     idx(i_recv, j_recv, offs, dims), &
                     e_pack, orient)
             end do
          end do
       end if

    case default
       error stop "create_comm_e: invalid side index"
    end select
  end subroutine create_comm_e

  
  function rot_direction (dom, typ) result(val)
    implicit none
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: typ
    integer                     :: val

    integer :: t_last, t_next

    t_last = 0
    t_next = 0
    if (typ <= 4) then
       val = modulo(typ, 2)
       return
    else
       t_last = typ - 4
       if (dom%neigh_rot(t_last) == 1) then
          val = modulo (t_last, 2)
          return
       else
          t_next = (modulo (typ, 4) + 1) - 4
          val = modulo (t_next, 2)
          return
       end if
    end if
  end function rot_direction

  
  subroutine comm_communication
    implicit none

    integer, parameter :: RECORD_SIZE = 4

    integer :: i
    integer :: dest_glo, dest_loc
    integer :: src_glo, src_loc
    integer :: n_conn
    integer :: st(RECORD_SIZE)

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          n_conn = grid(src_loc)%send_conn(dest_glo+1)%length

          if (mod(n_conn, RECORD_SIZE) /= 0) then
             error stop &
                  "comm_communication: send_conn length not divisible by 4"
          end if

          do i = 1, n_conn, RECORD_SIZE
             st = grid(src_loc)%send_conn(dest_glo+1)%elts(i:i+3)

             call unpack_comm_struct( &
                  grid(dest_loc), src_glo, &
                  st(1), st(2), st(3), st(4))
          end do

          grid(src_loc)%send_conn(dest_glo+1)%length = 0
       end do
    end do
  end subroutine comm_communication

  
  subroutine create_comm (dom, p, s)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, s

    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    integer :: b, e, pa, rot, src, typ

    call get_offs_Domain (dom, p, offs, dims)

    b   = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side

    if (typ < 1) return

    src = dom%neigh(typ)
    rot = dom%neigh_rot(typ)
    pa  = dom%bdry_patch%elts(b+1)%neigh

    do e = 0, EDGE
       call create_comm_e(dom, p, s, e)
    end do

    if (rot == 1 .and. typ == WEST) then
       if (is_penta(dom, p, SOUTHWEST-1)) then
          call create_pack_st ( &
               dom, AT_EDGE, src, LAST, LAST, pa, RT, &
               nidx(LAST_BDRY, LAST_BDRY, SOUTHWEST, offs, dims)*EDGE)
       end if
    end if

    if (rot == 1 .and. typ == SOUTH) then
       if (is_penta(dom, p, SOUTHWEST-1)) then
          call create_pack_st ( &
               dom, AT_EDGE, src, LAST, LAST, pa, UP, &
               -nidx(LAST_BDRY, LAST_BDRY, SOUTHWEST, offs, dims)*EDGE)
       end if
    end if

    if (rot == 0 .and. typ == EAST) then
       if (is_penta(dom, p, SOUTHEAST-1)) then
          call create_pack_st ( &
               dom, AT_NODE, src, 1, 0, pa, NODE, &
               nidx(0, LAST_BDRY, SOUTHEAST, offs, dims))
       end if
    end if

    if (rot == 0 .and. typ == NORTH) then
       if (is_penta(dom, p, NORTHWEST-1)) then
          call create_pack_st ( &
               dom, AT_NODE, src, 0, 1, pa, NODE, &
               nidx(LAST_BDRY, 0, NORTHWEST, offs, dims))
       end if
    end if

    if (rot == 1 .and. typ == EAST) then
       if (is_penta(dom, p, NORTHEAST-1)) then
          call create_pack_st ( &
               dom, AT_NODE, src, 0, 1, pa, NODE, &
               nidx(0, 1, NORTHEAST, offs, dims))

          call create_pack_st ( &
               dom, AT_NODE, src, 1, 1, pa, NODE, &
               -nidx(1, 0, NORTHEAST, offs, dims))

          call create_pack_st ( &
               dom, AT_EDGE, src, 0, 1, pa, RT, &
               nidx(0, 1, NORTHEAST, offs, dims)*EDGE + RT)

          call create_pack_st ( &
               dom, AT_EDGE, src, 0, 0, pa, DG, &
               -(nidx(0, 0, NORTHEAST, offs, dims)*EDGE + RT))

          call create_pack_st ( &
               dom, AT_EDGE, src, 0, 0, pa, UP, &
               nidx(0, 0, NORTHEAST, offs, dims)*EDGE + UP)
       end if
    end if

    if (rot == 1 .and. typ == NORTH) then
       if (is_penta(dom, p, NORTHEAST-1)) then
          call create_pack_st ( &
               dom, AT_NODE, src, 1, 0, pa, NODE, &
               nidx(1, 0, NORTHEAST, offs, dims))

          call create_pack_st ( &
               dom, AT_NODE, src, 1, 1, pa, NODE, &
               -nidx(0, 1, NORTHEAST, offs, dims))

          call create_pack_st ( &
               dom, AT_EDGE, src, 1, 0, pa, UP, &
               -(nidx(0, 1, NORTHEAST, offs, dims)*EDGE + RT))

          call create_pack_st ( &
               dom, AT_EDGE, src, 0, 0, pa, DG, &
               -(nidx(0, 0, NORTHEAST, offs, dims)*EDGE + UP))

          call create_pack_st ( &
               dom, AT_EDGE, src, 0, 0, pa, RT, &
               nidx(0, 0, NORTHEAST, offs, dims)*EDGE + RT)
       end if
    end if
  end subroutine create_comm

  
  subroutine create_comm_pole (dom, p, s, i)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, s, i

    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)
    integer :: ij_node(2), ij_send(2)

    integer :: b, lev, pa, src, typ
    integer :: s_side, s_side_second

    ! This routine is called for every boundary side. Only these two
    ! corner sides correspond to pole communication.
    select case (s)
    case (NORTHWEST)
       s_side        = NORTH
       s_side_second = WEST

    case (SOUTHEAST)
       s_side        = EAST
       s_side_second = SOUTH

    case default
       return
    end select

    ! Side s of patch p is a pole, and the i-th connection over
    ! this pole is available.
    b   = -dom%patch%elts(p+1)%neigh(s)
    typ = dom%bdry_patch%elts(b+1)%side

    if (typ < 1) return
    if (dom%neigh(typ) /= POLE) return

    lev = dom%patch%elts(p+1)%level
    src = dom%neigh_over_pole(i+1)
    pa  = dom%neigh_pa_over_pole%elts(i + 2*lev + 1)

    call get_offs_Domain(dom, p, offs, dims)

    if (i == 0) then
       ij_node = [0, 1]

       if (s == SOUTHEAST) then
          ij_node = [ij_node(2), ij_node(1)]
       end if

       call create_pack_st ( &
            dom, AT_EDGE, src, 0, LAST, pa, DG, &
            nidx(ij_node(1), ij_node(2), s_side, offs, dims)*EDGE &
            + 2*s_side - 2)

       ij_node = [LAST_BDRY, 1]
       ij_send = [1, LAST]

       if (s == SOUTHEAST) then
          ij_node = [ij_node(2), ij_node(1)]
          ij_send = [ij_send(2), ij_send(1)]
       end if

       call create_pack_st ( &
            dom, AT_NODE, src, ij_send(1), ij_send(2), pa, NODE, &
            -nidx(ij_node(1), ij_node(2), s, offs, dims))
    end if

    if (i == 1) then
       ij_node = [LAST_BDRY, 0]
       ij_send = [0, LAST]

       if (s == SOUTHEAST) then
          ij_node = [ij_node(2), ij_node(1)]
          ij_send = [ij_send(2), ij_send(1)]
       end if

       call create_pack_st ( &
            dom, AT_NODE, src, ij_send(1), ij_send(2), pa, NODE, &
            nidx(ij_node(1), ij_node(2), s, offs, dims))
    end if

    ij_node = [1-i, 1]
    ij_send = [0, LAST]

    if (s == SOUTHEAST) then
       ij_node = [ij_node(2), ij_node(1)]
       ij_send = [ij_send(2), ij_send(1)]
    end if

    call create_pack_st ( &
         dom, AT_NODE, src, ij_send(1), ij_send(2), pa, NODE, &
         (-1)**i*nidx( &
         ij_node(1), ij_node(2), s_side, offs, dims))
    
    call create_pack_st ( &
         dom, AT_EDGE, src, ij_send(1), ij_send(2), pa, &
         UP - 2*s_side + 2, &
         (-1)**i*( &
         nidx(0, 0, s_side, offs, dims)*EDGE + DG &
         + i*(-2*s_side + 3)))

    if (i == 0) return

    s_side = s_side_second

    ij_node = [LAST_BDRY-1, LAST]
    ij_send = [1, LAST]

    if (s == SOUTHEAST) then
       ij_node = [ij_node(2), ij_node(1)]
       ij_send = [ij_send(2), ij_send(1)]
    end if

    call create_pack_st ( &
         dom, AT_NODE, src, ij_send(1), ij_send(2), pa, NODE, &
         nidx(ij_node(1), ij_node(2), s_side, offs, dims))

    ij_node = [LAST, LAST_BDRY]

    if (s == NORTHWEST) then
       ij_node = [ij_node(2), ij_node(1)]
    end if

    call create_pack_st ( &
         dom, AT_EDGE, src, &
         LAST*(-s_side + 4), LAST*(s_side - 3), pa, DG, &
         nidx(ij_node(1), ij_node(2), s_side, offs, dims)*EDGE &
         + 2*s_side - 6)
  end subroutine create_comm_pole

  
  subroutine comm_masks
    implicit none

    integer :: dest_glo, dest_id, dest_loc
    integer :: src_glo, src_id, src_loc
    integer :: i

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          do i = 1, grid(src_loc)%pack(AT_NODE,dest_glo+1)%length
             src_id = &
                  grid(src_loc)%pack(AT_NODE,dest_glo+1)%elts(i)

             dest_id = &
                  grid(dest_loc)%unpk(AT_NODE,src_glo+1)%elts(i)

             grid(dest_loc)%mask_n%elts(abs(dest_id)+1) = &
                  grid(src_loc)%mask_n%elts(abs(src_id)+1)
          end do

          do i = 1, grid(src_loc)%pack(AT_EDGE,dest_glo+1)%length
             src_id = &
                  grid(src_loc)%pack(AT_EDGE,dest_glo+1)%elts(i)

             dest_id = &
                  grid(dest_loc)%unpk(AT_EDGE,src_glo+1)%elts(i)

             grid(dest_loc)%mask_e%elts(abs(dest_id)+1) = &
                  grid(src_loc)%mask_e%elts(abs(src_id)+1)
          end do
       end do
    end do
  end subroutine comm_masks

  
  subroutine comm_edges (get, set)
    implicit none

    procedure(real_get) :: get
    procedure(real_set) :: set

    integer :: dest_glo, dest_id, dest_loc
    integer :: src_glo, src_id, src_loc
    integer :: i

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          do i = 1, grid(src_loc)%pack(AT_EDGE,dest_glo+1)%length
             src_id = &
                  grid(src_loc)%pack(AT_EDGE,dest_glo+1)%elts(i)

             dest_id = &
                  grid(dest_loc)%unpk(AT_EDGE,src_glo+1)%elts(i)

             call set ( &
                  grid(dest_loc), dest_id, &
                  get(grid(src_loc), src_id))
          end do
       end do
    end do
  end subroutine comm_edges

  
  subroutine unpack_comm_struct (dom, src, i, j, p, e)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: src, i, j, p, e

    integer :: offs

    offs = dom%patch%elts(p+1)%elts_start

    if (e == NODE) then
       call append ( &
            dom%pack(AT_NODE,src+1), &
            idx__fast(i, j, offs))
    else
       call append ( &
            dom%pack(AT_EDGE,src+1), &
            idx__fast(i, j, offs)*EDGE + e)
    end if
  end subroutine unpack_comm_struct

  
  subroutine pack_idx (i, j, rot, s, shift, i_recv, j_recv, i_pack, j_pack)
    implicit none

    integer, intent(in)  :: i, j, rot, s, shift
    integer, intent(out) :: i_recv, j_recv, i_pack, j_pack

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

    type(Domain), intent(inout) :: dom

    integer :: c, ii, lev
    integer :: n_chd, n_par, ngb_pa
    integer :: p_chd, p_par, s
    integer :: src_glo, typ
    integer :: unused_elements

    do src_glo = 1, N_GLO_DOMAIN
       unused_elements = 0

       do ii = 1, dom%recv_pa(src_glo)%length, 4
          ngb_pa = dom%recv_pa(src_glo)%elts(ii)
          c      = dom%recv_pa(src_glo)%elts(ii+1)
          p_par  = dom%recv_pa(src_glo)%elts(ii+2)
          s      = dom%recv_pa(src_glo)%elts(ii+3)

          n_par = -dom%patch%elts(p_par+1)%neigh(s+1)
          typ   = dom%bdry_patch%elts(n_par+1)%side

          if (typ < 1) return

          if (dom%neigh(typ) == POLE) then
             p_chd = dom%patch%elts(p_par+1)%children(s-3)
          else
             p_chd = dom%patch%elts(p_par+1)%children(c+1)
          end if

          if (p_chd == 0) then
             dom%recv_pa(src_glo)%elts( &
                  unused_elements+1:unused_elements+4) = &
                  [ngb_pa, c, p_par, s]

             unused_elements = unused_elements + 4
             cycle
          end if

          if (dom%neigh(typ) == POLE) then
             lev = dom%patch%elts(p_chd+1)%level

             if (2*lev + 2 > dom%neigh_pa_over_pole%length) then
                call extend (dom%neigh_pa_over_pole, 2, 0)
             end if

             dom%neigh_pa_over_pole%elts((1-c)+2*lev+1) = ngb_pa

             call create_comm_pole (dom, p_chd, s+1, 1-c)
             cycle
          end if

          n_chd = -dom%patch%elts(p_chd+1)%neigh(s+1)
          dom%bdry_patch%elts(n_chd+1)%neigh = ngb_pa

          if (.not. is_penta(dom, p_chd, s)) then
             call create_comm (dom, p_chd, s+1)
          end if
       end do

       dom%recv_pa(src_glo)%length = unused_elements
    end do
  end subroutine update_comm

  
  subroutine comm_patch_conn
    implicit none

    integer :: b, c, dest, i, l_par
    integer :: ngh_pa, p_chd
    integer :: rot, rot_shift, s
    integer :: src, src_glo, typ
    integer :: unused_elements
    logical :: is_pole

    do src = 1, size(grid)
       unused_elements = 0
       src_glo = glo_id(rank+1, src)

       do i = 1, grid(src)%send_pa_all%length, 4
          b     = grid(src)%send_pa_all%elts(i)
          c     = grid(src)%send_pa_all%elts(i+1)
          p_chd = grid(src)%send_pa_all%elts(i+2)
          s     = grid(src)%send_pa_all%elts(i+3)

          typ = grid(src)%bdry_patch%elts(b+1)%side

          if (typ < 1) return

          dest = grid(src)%neigh(typ)
          is_pole = dest == POLE

          if (is_pole) then
             dest = grid(src)%neigh_over_pole(c+1)
             l_par = grid(src)%patch%elts(p_chd+1)%level - 1

             if (grid(src)%neigh_pa_over_pole%length < &
                  2*l_par + c + 1) then
                ngh_pa = 0
             else
                ngh_pa = grid(src)%neigh_pa_over_pole%elts( &
                     2*l_par + c + 1)
             end if
          else
             ngh_pa = grid(src)%bdry_patch%elts(b+1)%neigh
          end if

          if (ngh_pa == 0) then
             grid(src)%send_pa_all%elts( &
                  unused_elements+1:unused_elements+4) = &
                  [b, c, p_chd, s]

             unused_elements = unused_elements + 4
             cycle
          else
             if (dest == NONE) cycle
          end if

          ! Skip if destination is not on this rank.
          if (owner(dest+1) /= rank) cycle

          dest = loc_id(dest+1)
          rot = grid(src)%neigh_rot(typ)
          rot_shift = &
               (2*rot_direction(grid(src),typ)-1)*rot

          call append ( &
               grid(dest+1)%recv_pa(src_glo+1), p_chd)

          if (is_pole) then
             call append ( &
                  grid(dest+1)%recv_pa(src_glo+1), c)
          else
             call append ( &
                  grid(dest+1)%recv_pa(src_glo+1), &
                  modulo(c + rot_shift, 4))
          end if

          call append ( &
               grid(dest+1)%recv_pa(src_glo+1), ngh_pa)

          if (is_pole) then
             call append ( &
                  grid(dest+1)%recv_pa(src_glo+1), s)
          else
             call append ( &
                  grid(dest+1)%recv_pa(src_glo+1), &
                  modulo(rot_shift + s + 2, 4) + 4*(s/4))
          end if
       end do

       grid(src)%send_pa_all%length = unused_elements
    end do
  end subroutine comm_patch_conn

  
  subroutine create_pack_st (dom, unpk_pos, src, i, j, pa, e, id)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: unpk_pos
    integer,      intent(in)    :: src, i, j, pa, e, id

    call append(dom%send_conn(src+1), i)
    call append(dom%send_conn(src+1), j)
    call append(dom%send_conn(src+1), pa)
    call append(dom%send_conn(src+1), e)

    call append(dom%unpk(unpk_pos,src+1), id)
  end subroutine create_pack_st

  
  subroutine cp_bdry_inside_0 (field)
    implicit none

    type(Float_Field), intent(inout) :: field

    integer :: dest_id, dest_loc, dest_glo
    integer :: src_id, src_loc, src_glo
    integer :: i, pos

    pos = field%pos

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
             src_id = &
                  grid(src_loc)%pack(pos,dest_glo+1)%elts(i)

             dest_id = &
                  grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)

             field%data(dest_loc)%elts(abs(dest_id)+1) = &
                  field%data(src_loc)%elts(src_id+1)

             if (dest_id < 0 .and. pos == AT_EDGE) then
                field%data(dest_loc)%elts(abs(dest_id)+1) = &
                     -field%data(dest_loc)%elts(abs(dest_id)+1)
             end if
          end do
       end do
    end do
  end subroutine cp_bdry_inside_0

  
  subroutine cp_bdry_inside_1 (field)
    implicit none

    type(Float_Field), intent(inout) :: field(:)

    integer :: dest_id, dest_loc, dest_glo
    integer :: src_id, src_loc, src_glo
    integer :: i, i1, pos

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          do i1 = 1, size(field)
             pos = field(i1)%pos

             do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
                src_id = &
                     grid(src_loc)%pack(pos,dest_glo+1)%elts(i)

                dest_id = &
                     grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)

                field(i1)%data(dest_loc)%elts(abs(dest_id)+1) = &
                     field(i1)%data(src_loc)%elts(src_id+1)

                if (dest_id < 0 .and. pos == AT_EDGE) then
                   field(i1)%data(dest_loc)%elts(abs(dest_id)+1) = &
                        -field(i1)%data(dest_loc)%elts(abs(dest_id)+1)
                end if
             end do
          end do
       end do
    end do
  end subroutine cp_bdry_inside_1

  
  subroutine cp_bdry_inside_2 (field)
    implicit none

    type(Float_Field), intent(inout) :: field(:,:)

    integer :: dest_id, dest_loc, dest_glo
    integer :: src_id, src_loc, src_glo
    integer :: i, i1, i2, pos

    do src_loc = 1, size(grid)
       src_glo = glo_id(rank+1, src_loc)

       do dest_loc = 1, size(grid)
          dest_glo = glo_id(rank+1, dest_loc)

          do i2 = 1, size(field,2)
             do i1 = 1, size(field,1)
                pos = field(i1,i2)%pos

                do i = 1, grid(src_loc)%pack(pos,dest_glo+1)%length
                   src_id = &
                        grid(src_loc)%pack(pos,dest_glo+1)%elts(i)

                   dest_id = &
                        grid(dest_loc)%unpk(pos,src_glo+1)%elts(i)

                   field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1) = &
                        field(i1,i2)%data(src_loc)%elts(src_id+1)

                   if (dest_id < 0 .and. pos == AT_EDGE) then
                      field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1) = &
                           -field(i1,i2)%data(dest_loc)%elts(abs(dest_id)+1)
                   end if
                end do
             end do
          end do
       end do
    end do
  end subroutine cp_bdry_inside_2

  
  pure function domain_load (dom) result(load)
    implicit none

    type(Domain), intent(in) :: dom
    integer                  :: load

    load = count( &
         abs(dom%mask_n%elts(2:dom%node%length)) > ADJZONE) &
         + count( &
         abs(dom%mask_e%elts(EDGE+1:dom%midpt%length)) > ADJZONE)
  end function domain_load
end module comm_mod
