module parallel_block_primitive_walk_mod
  ! Topology-only enumeration of the exact regular/SW operation order.
  use shared_mod, only: N_BDRY,NORTH,EAST,SOUTH,WEST,NORTHEAST,NORTHWEST,SOUTHEAST,SOUTHWEST
  use domain_mod, only: Domain,sides_dims
  use patch_mod, only: PATCH_SIZE,LAST
  implicit none
  private
  public :: walk_primitive_patch,comp_offs3
  abstract interface
     subroutine Primitive_Visitor(ids,southwest)
       integer, intent(in) :: ids(0:7)
       logical, intent(in) :: southwest
     end subroutine
  end interface
contains
  subroutine walk_primitive_patch(dom,p,visit)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: p
    procedure(Primitive_Visitor) :: visit
    integer :: j,id,n,e,s,w,ne,sw,offs(0:N_BDRY),dims(2,N_BDRY)
    logical :: S_bdry,W_bdry
    call comp_offs3 (dom, p, offs, dims)

    S_bdry = (dom%patch%elts(p+1)%neigh(SOUTH) < 0)
    if (S_bdry) S_bdry = (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(SOUTH)+1)%side > 0)

    W_bdry = (dom%patch%elts(p+1)%neigh(WEST)  < 0)
    if (W_bdry) W_bdry = (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(WEST)+1)%side > 0)

    id = offs(0)

    n = PATCH_SIZE       ; ne = n+1
    w = offs(WEST)       ; e =  +1
    sw = offs(SOUTHWEST) ; s = offs(SOUTH)

    call emit(.false.)
    if (W_bdry .or. S_bdry) call emit(.true.)

    w = -1; sw = s+w
    do id = offs(0)+1, offs(0)+LAST-1
       call emit(.false.)
       if (S_bdry) call emit(.true.)
    end do

    e = offs(EAST); ne = dims(1,EAST) + e
    id = offs(0)+LAST
    call emit(.false.)
    if (S_bdry) call emit(.true.)

    s = -PATCH_SIZE
    do j = 2, PATCH_SIZE-1
       id = offs(0) + PATCH_SIZE*(j-1)
       e = +1; ne = n+e
       w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*(j-1) ! Correct for dimension smaller than patch if boundary
       sw = w-dims(1,WEST)
       call emit(.false.)
       if (W_bdry) call emit(.true.)

       w = -1; sw = s+w; ne = n+e
       do id = offs(0)+PATCH_SIZE*(j-1)+1, offs(0)+PATCH_SIZE*(j-1)+LAST-1
          call emit(.false.)
       end do

       id = offs(0)+PATCH_SIZE*j-1
       e = offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*(j-1); ne = e+dims(1,EAST)
       call emit(.false.)
    end do

    id = offs(0)+PATCH_SIZE*LAST
    n = offs(NORTH)
    w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*LAST
    sw = w-dims(1,WEST)
    e = +1; ne = n+e
    call emit(.false.)
    if (W_bdry) call emit(.true.)

    w = -1; sw = s+w
    do id = offs(0)+PATCH_SIZE*LAST+1, offs(0)+PATCH_SIZE*PATCH_SIZE-2
       call emit(.false.)
    end do

    id = offs(0)+PATCH_SIZE*PATCH_SIZE-1
    ne = offs(NORTHEAST)
    e = offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*LAST
    call emit(.false.)

    if (dom%patch%elts(p+1)%neigh(NORTH) < 0) then ! Neighbour is boundary
       if (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(NORTH)+1)%side > 0) then ! Domain boundary
          id = offs(0)+PATCH_SIZE*LAST + offs(NORTH) ! id + n
          s =                                   - offs(NORTH) ! relative to current id
          w = offs(NORTHWEST)                   - offs(NORTH) 
          sw = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*LAST - offs(NORTH)
          e = 1
          n = dims(1,NORTH)
          ne = n+e
          call emit(.false.)
          if (W_bdry) call emit(.true.)

          w = -1
          sw = -offs(NORTH) + w
          do id = offs(0)+PATCH_SIZE*LAST+1+offs(NORTH), offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTH)-2
             call emit(.false.)
          end do

          e = offs(NORTHEAST)              - offs(NORTH) 
          ne = e+dims(1,NORTHEAST)
          id = offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTH)-1
          call emit(.false.)

          ! NORTHEAST point
          id = offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTHEAST)-1
          s = -offs(NORTHEAST) + offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*LAST
          sw= -offs(NORTHEAST)
          w = -offs(NORTHEAST) + offs(NORTH)
          n = dims(1,NORTHEAST)
          e = 1
          ne = n+e
          call emit(.false.)
       end if
    end if

    if (dom%patch%elts(p+1)%neigh(EAST) < 0) then ! neighbour is boundary
       if (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(EAST)+1)%side <= 0) return ! not domain bdry
       id = offs(0) + LAST + offs(EAST) ! id + n
       w  =                - offs(EAST) ! relative to current id
       s  = offs(SOUTHEAST)- offs(EAST) 
       sw = offs(SOUTH)    - offs(EAST)
       n  = dims(1,EAST)
       e  = 1
       ne = n+e
       call emit(.false.)
       if (S_bdry) call emit(.true.)

       s = -dims(1,EAST)
       do j = 1, LAST-1
          id = offs(0) + LAST + offs(EAST) + j*dims(1,EAST)
          w  =                - offs(EAST) + j*(PATCH_SIZE-dims(1,EAST))
          sw = w-PATCH_SIZE
          call emit(.false.)
       end do

       id = offs(0) + LAST + offs(EAST) + LAST*dims(1,EAST)
       n  =   offs(NORTHEAST) + LAST*(PATCH_SIZE-dims(1,EAST)) - offs(EAST)
       w  = - offs(EAST)      + LAST*(PATCH_SIZE-dims(1,EAST))
       sw = w-PATCH_SIZE
       ne = n+e
       call emit(.false.)
    end if

  contains
    subroutine emit(southwest)
      logical, intent(in) :: southwest
      call visit([id,id+n,id+e,id+s,id+w,id+ne,id,id+sw],southwest)
    end subroutine
  end subroutine
  subroutine comp_offs3 (dom, p, offs, dims)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p
    integer,      intent(out)   :: offs(0:N_BDRY)
    integer,      intent(out)   :: dims(2,N_BDRY)

    integer :: i, n

    offs(0) = dom%patch%elts(p+1)%elts_start
    do i = 1, N_BDRY
       n = dom%patch%elts(p+1)%neigh(i)
       if (n > 0) then ! regular patch
          offs(i)  = dom%patch%elts(n+1)%elts_start
          dims(:,i) = PATCH_SIZE
       elseif (n < 0) then
          offs(i)  = dom%bdry_patch%elts(-n+1)%elts_start
          dims(:,i) = sides_dims(:,abs(dom%bdry_patch%elts(-n+1)%side)+1)
       end if
    end do

    offs(1:N_BDRY)  = offs(1:N_BDRY)  - offs(0) ! make relative
    offs(SOUTH)     = offs(SOUTH)     + dims(1,SOUTH)*(dims(2,SOUTH)-1)
    offs(SOUTHEAST) = offs(SOUTHEAST) + dims(1,SOUTHEAST)*(dims(2,SOUTH)-1)
    offs(WEST)      = offs(WEST)      + (dims(1,WEST)-1)
    offs(NORTHWEST) = offs(NORTHWEST) + (dims(1,NORTHWEST)-1)
    offs(SOUTHWEST) = offs(SOUTHWEST) + dims(1,SOUTHWEST)*dims(2,SOUTHWEST) - 1
    offs(NORTH)     = offs(NORTH)     - PATCH_SIZE*(PATCH_SIZE-1) 
    offs(NORTHWEST) = offs(NORTHWEST) - PATCH_SIZE*(PATCH_SIZE-1) 
    offs(EAST)      = offs(EAST)      - (PATCH_SIZE-1)
    offs(SOUTHEAST) = offs(SOUTHEAST) - (PATCH_SIZE-1)
    offs(NORTHEAST) = offs(NORTHEAST) - (PATCH_SIZE*PATCH_SIZE-1)
  end subroutine comp_offs3
end module
