
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !! Subroutines for applying a routine to domains, patches, ...  !!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module domain_ops_mod

  use kind_mod,    only : dp

  use shared_mod, only : N_BDRY, BDRY_THICKNESS, EDGE, level_start, level_end, IMINUS, IPLUS, JMINUS, JPlUS,  N_CHDRN, &
       NONE, NORTHEAST, NORTHWEST, POLE, SOUTHEAST

  use domain_mod, only : chd_offs, Domain, get_offs_Domain, get_offs_Domain5, grid
  use patch_mod,  only : PATCH_SIZE

  implicit none

  private
  public :: fun3, fun4, sub4, sub5, sub6, sub7, sub8, sub9, sub10
  public :: apply, apply_bdry, apply_onescale, apply_onescale__int, apply_onescale_to_patch__int, apply_no_bdry
  public :: apply_no_bdry_to_patch, apply_no_bdry2, apply_d, apply_onescale_to_patch, apply_onescale_d, apply_onescale2
  public :: apply_onescale_to_patch5, apply_onescale_to_patch2, apply_interscale, apply_interscale_d, apply_interscale_d2
  public :: apply_interscale_to_patch, apply_interscale_to_patch2, apply_interscale_to_patch22, apply_interscale_to_patch3
  public :: apply_to_pole, apply_to_pole_d, apply_to_pole_patch, apply_to_pole2, apply_to_penta, apply_to_penta_d

  
  interface

     function fun3 (dom, i, j, zlev, offs, dims) result(val)
       use kind_mod,   only : dp
       use shared_mod, only : N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: i, j, zlev
       integer,      intent(in)    :: offs(N_BDRY+1)
       integer,      intent(in)    :: dims(2,N_BDRY+1)
       
       real(dp) :: val
     end function fun3

     function fun4 (dom, i, j, zlev, offs, dims) result(val)
       use kind_mod,   only : dp
       use shared_mod, only : EDGE, N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: i, j, zlev
       integer,      intent(in)    :: offs(N_BDRY+1)
       integer,      intent(in)    :: dims(2,N_BDRY+1)

       real(dp) :: val(EDGE)
     end function fun4

     subroutine sub4 (dom, i, j, zlev, offs, dims)
       use shared_mod, only : N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: i, j, zlev
       integer,      intent(in)    :: offs(N_BDRY+1)
       integer,      intent(in)    :: dims(2,N_BDRY+1)
     end subroutine sub4

     subroutine sub5 (dom, p, i, j, zlev, offs, dims)
       use shared_mod, only : N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: p, i, j, zlev
       integer,      intent(in)    :: offs(N_BDRY+1)
       integer,      intent(in)    :: dims(2,N_BDRY+1)
     end subroutine sub5

     subroutine sub6 (  &
          dom, i_par, j_par, i_chd, j_chd, zlev, &
          offs_par, dims_par, offs_chd, dims_chd)

       use shared_mod, only : N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: i_par, j_par
       integer,      intent(in)    :: i_chd, j_chd, zlev
       integer,      intent(in)    :: offs_par(N_BDRY+1)
       integer,      intent(in)    :: dims_par(2,N_BDRY+1)
       integer,      intent(in)    :: offs_chd(N_BDRY+1)
       integer,      intent(in)    :: dims_chd(2,N_BDRY+1)
     end subroutine sub6

     subroutine sub7 ( &
          dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, &
          offs_par, dims_par, offs_chd, dims_chd)

       use shared_mod, only : N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: p_chd
       integer,      intent(in)    :: i_par, j_par
       integer,      intent(in)    :: i_chd, j_chd, zlev
       integer,      intent(in)    :: offs_par(N_BDRY+1)
       integer,      intent(in)    :: dims_par(2,N_BDRY+1)
       integer,      intent(in)    :: offs_chd(N_BDRY+1)
       integer,      intent(in)    :: dims_chd(2,N_BDRY+1)
     end subroutine sub7

     subroutine sub8 (dom, p, i, j, zlev, offs, dims, ival)
       use shared_mod, only : N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: p, i, j, zlev, ival
       integer,      intent(in)    :: offs(N_BDRY+1)
       integer,      intent(in)    :: dims(2,N_BDRY+1)
     end subroutine sub8

     subroutine sub9 (dom, p, c, offs, dims, zlev)
       use shared_mod, only : N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: p, c, zlev
       integer,      intent(in)    :: offs(N_BDRY+1)
       integer,      intent(in)    :: dims(2,N_BDRY+1)
     end subroutine sub9


     subroutine sub10 (dom, p, i, j, zlev, offs, dims)
       use shared_mod, only : N_BDRY
       use domain_mod, only : Domain
       implicit none

       type(Domain), intent(inout) :: dom
       integer,      intent(in)    :: p, i, j, zlev
       integer,      intent(in)    :: offs(N_BDRY+1)
       integer,      intent(in)    :: dims(2,N_BDRY+1)
     end subroutine sub10
     
  end interface

  
contains

  
  subroutine apply (routine, zlev)
    ! Apply routine over all levels and the full boundary range.
    implicit none

    procedure(sub4)     :: routine
    integer, intent(in) :: zlev

    integer :: l

    do l = level_start, level_end
       call apply_onescale (routine, l, zlev, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
  end subroutine apply

  
  subroutine apply_bdry (routine, zlev, st, en)
    ! Apply routine at all levels, including boundary cells
    ! specified by the range (st,en).
    implicit none

    procedure(sub4)     :: routine
    integer, intent(in) :: zlev, st, en

    integer :: l

    do l = level_start, level_end
       call apply_onescale (routine, l, zlev, st, en)
    end do
  end subroutine apply_bdry

  
  subroutine apply_onescale(routine, l, zlev, st, en)
    ! Apply routine at level l, including boundary cells
    ! specified by the range (st,en).
    implicit none

    procedure(sub4)     :: routine
    integer, intent(in) :: l, zlev, st, en

    integer :: d, j

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (routine, grid(d), grid(d)%lev(l)%elts(j), zlev, st, en)
       end do
    end do
  end subroutine apply_onescale

  
  subroutine apply_onescale__int(routine, l, zlev, st, en, ival)
    ! Apply routine at level l, including boundary cells specified
    ! by (st,en), and pass integer ival to routine.
    implicit none

    procedure(sub8)     :: routine
    integer, intent(in) :: l, zlev, st, en, ival

    integer :: d, j

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch__int( &
               routine, grid(d), grid(d)%lev(l)%elts(j), &
               zlev, st, en, ival)
       end do
    end do
  end subroutine apply_onescale__int

  
  subroutine apply_onescale_to_patch__int (routine, dom, p, zlev, st, en, ival)

    implicit none

    procedure(sub8)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, zlev, st, en, ival

    integer :: i, j
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    call get_offs_Domain(dom, p, offs, dims)

    do j = st+1, PATCH_SIZE+en
       do i = st+1, PATCH_SIZE+en
          call routine( &
               dom, p, i-1, j-1, zlev, offs, dims, ival)
       end do
    end do
  end subroutine apply_onescale_to_patch__int


  subroutine apply_no_bdry(routine, zlev)
    ! Apply routine at all levels, including poles but excluding
    ! boundary cells.
    implicit none

    procedure(sub4)     :: routine
    integer, intent(in) :: zlev

    integer :: l

    do l = level_start, level_end
       call apply_onescale(routine, l, zlev, 0, 0)
       call apply_to_pole2(routine, l, zlev)
    end do
  end subroutine apply_no_bdry

  
  subroutine apply_no_bdry_to_patch(routine, dom, p, zlev)
    ! Apply routine to one patch, excluding boundary cells.
    ! Apply the pole point only when patch p contains a pole.
    implicit none

    procedure(sub4)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, zlev

    call apply_onescale_to_patch (routine, dom, p, zlev, 0, 0)
    call apply_to_pole_one_patch (routine, dom, p, zlev)
  end subroutine apply_no_bdry_to_patch

  
  subroutine apply_to_pole_one_patch (routine, dom, p, zlev)
    ! Apply routine to the pole point of patch p, provided p lies
    ! on one of the domain's pole branches.
    implicit none

    procedure(sub4)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, zlev

    integer :: c, p_cur
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    do c = SOUTHEAST, NORTHWEST, 2
       if (.not. dom%penta(c)) cycle
       if (dom%neigh(c) /= POLE) cycle

       ! Follow the refinement branch containing this pole.
       p_cur = 1

       do while (p_cur > 0)
          if (p_cur == p) then
             call get_offs_Domain (dom, p, offs, dims)

             select case (c)
             case (NORTHWEST)
                call routine (dom, 0, PATCH_SIZE, zlev, offs, dims)

             case (SOUTHEAST)
                call routine (dom, PATCH_SIZE, 0, zlev, offs, dims)
             end select

             return
          end if

          p_cur = dom%patch%elts(p_cur+1)%children(c-4)
       end do
    end do
  end subroutine apply_to_pole_one_patch

  
  subroutine apply_no_bdry2 (routine, zlev)
    ! Apply routine at all levels, excluding boundary cells.
    !
    ! is_pole is zero for ordinary points and one for pole points.
    implicit none

    procedure(sub8)     :: routine
    integer, intent(in) :: zlev

    integer :: l

    do l = level_start, level_end
       call apply_onescale__int (routine, l, zlev, 0, 0, 0)
       call apply_to_pole (routine, l, zlev, 1, .true.)
    end do
  end subroutine apply_no_bdry2

  
  subroutine apply_d (routine, dom, zlev, st, en)
    implicit none

    procedure(sub4)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: zlev, st, en

    integer :: j, l

    do l = level_start, level_end
       do j = 1, dom%lev(l)%length
          call apply_onescale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
       end do
    end do
  end subroutine apply_d

  
  subroutine apply_onescale_to_patch (routine, dom, p, zlev, st, en)
    ! Apply routine to the nodes/edges associated with patch p,
    ! including boundary cells specified by (st,en).
    implicit none

    procedure(sub4)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, zlev, st, en

    integer :: i, j
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)
    integer :: bdry(JPLUS:IMINUS)
    logical :: inner_bdry(JPLUS:IMINUS)

    call get_offs_Domain (dom, p, offs, dims, inner_bdry)

    bdry = [en, en, st, st]

    where (inner_bdry)
       bdry = 0
    end where

    do j = bdry(JMINUS)+1, PATCH_SIZE+bdry(JPLUS)
       do i = bdry(IMINUS)+1, PATCH_SIZE+bdry(IPLUS)
          call routine (dom, i-1, j-1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch

  
  subroutine apply_onescale_d (routine, dom, l, zlev, st, en)
    implicit none

    procedure(sub4)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: l, zlev, st, en

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_onescale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_onescale_d

  
  subroutine apply_onescale2 (routine, l, zlev, st, en)
    implicit none

    procedure(sub5)     :: routine
    integer, intent(in) :: l, zlev, st, en

    integer :: d, j

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch2 ( &
               routine, grid(d), grid(d)%lev(l)%elts(j), &
               zlev, st, en)
       end do
    end do
  end subroutine apply_onescale2

  
  subroutine apply_onescale_to_patch5 (routine, dom, p, zlev, st, en)
    implicit none

    procedure(sub5)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, zlev, st, en

    integer :: i, j
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)
    integer :: bdry(JPLUS:IMINUS)
    logical :: inner_patch(JPLUS:IMINUS)

    call get_offs_Domain5 (dom, p, offs, dims, inner_patch)

    bdry = [en, en, st, st]

    where (inner_patch)
       bdry = 0
    end where

    do j = bdry(JMINUS)+1, PATCH_SIZE+bdry(JPLUS)
       do i = bdry(IMINUS)+1, PATCH_SIZE+bdry(IPLUS)
          call routine (dom, p, i-1, j-1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch5

  
  subroutine apply_onescale_to_patch2 (routine, dom, p, zlev, st, en)
    implicit none

    procedure(sub5)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, zlev, st, en

    integer :: i, j
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    call get_offs_Domain( dom, p, offs, dims)

    do j = st+1, PATCH_SIZE+en
       do i = st+1, PATCH_SIZE+en
          call routine (dom, p, i-1, j-1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch2

  
  subroutine apply_interscale (routine, l, zlev, st, en)
    ! Apply an interscale routine to coarse level l and fine level l+1.
    implicit none

    procedure(sub6)     :: routine
    integer, intent(in) :: l, zlev, st, en

    integer :: d

    do d = 1, size(grid)
       call apply_interscale_d (routine, grid(d), l, zlev, st, en)
    end do
  end subroutine apply_interscale

  
  subroutine apply_interscale_d (routine, dom, l, zlev, st, en)
    implicit none

    procedure(sub6)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: l, zlev, st, en

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_interscale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_interscale_d

  
  subroutine apply_interscale_d2 (routine, dom, l, zlev, st, en)
    implicit none

    procedure(sub6)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: l, zlev, st, en

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_interscale_to_patch22 (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_interscale_d2

  
  subroutine apply_interscale_to_patch (routine, dom, p_par, zlev, st, en)
    implicit none

    procedure(sub6)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p_par, zlev, st, en

    integer :: c, i, j
    integer :: i_chd, i_par, j_chd, j_par, p_chd
    integer :: offs_chd(N_BDRY+1), offs_par(N_BDRY+1)
    integer :: dims_chd(2,N_BDRY+1), dims_par(2,N_BDRY+1)
    integer :: bdry(JPLUS:IMINUS)
    logical :: inner_bdry(JPLUS:IMINUS)

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do c = 1, N_CHDRN
       p_chd = dom%patch%elts(p_par+1)%children(c)
       if (p_chd == 0) cycle

       call get_offs_Domain (dom, p_chd, offs_chd, dims_chd, inner_bdry)

       bdry = [en, en, st, st]

       where (inner_bdry)
          bdry = 0
       end where

       do j = bdry(JMINUS)+1, PATCH_SIZE/2+bdry(JPLUS)
          j_chd = 2*(j-1)
          j_par = j-1 + chd_offs(2,c)

          do i = bdry(IMINUS)+1, PATCH_SIZE/2+bdry(IPLUS)
             i_chd = 2*(i-1)
             i_par = i-1 + chd_offs(1,c)

             call routine( &
                  dom, i_par, j_par, i_chd, j_chd, zlev, &
                  offs_par, dims_par, offs_chd, dims_chd)
          end do
       end do
    end do
  end subroutine apply_interscale_to_patch

  
  subroutine apply_interscale_to_patch2 (routine, dom, p_par, zlev, st, en)
    implicit none

    procedure(sub7)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p_par, zlev, st, en

    integer :: c, i, j
    integer :: i_chd, i_par, j_chd, j_par, p_chd
    integer :: offs_chd(N_BDRY+1), offs_par(N_BDRY+1)
    integer :: dims_chd(2,N_BDRY+1), dims_par(2,N_BDRY+1)
    integer :: bdry(JPLUS:IMINUS)
    logical :: inner_bdry(JPLUS:IMINUS)

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do c = 1, N_CHDRN
       p_chd = dom%patch%elts(p_par+1)%children(c)
       if (p_chd == 0) cycle

       call get_offs_Domain( &
            dom, p_chd, offs_chd, dims_chd, inner_bdry)

       bdry = [en, en, st, st]

       where (inner_bdry)
          bdry = 0
       end where

       do j = bdry(JMINUS)+1, PATCH_SIZE/2+bdry(JPLUS)
          j_chd = 2*(j-1)
          j_par = j-1 + chd_offs(2,c)

          do i = bdry(IMINUS)+1, PATCH_SIZE/2+bdry(IPLUS)
             i_chd = 2*(i-1)
             i_par = i-1 + chd_offs(1,c)

             call routine( &
                  dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, &
                  offs_par, dims_par, offs_chd, dims_chd)
          end do
       end do
    end do
  end subroutine apply_interscale_to_patch2

  
  subroutine apply_interscale_to_patch22 (routine, dom, p_par, zlev, st, en)
    implicit none

    procedure(sub6)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p_par, zlev, st, en

    integer :: c, i, j
    integer :: i_chd, i_par, j_chd, j_par, p_chd
    integer :: offs_chd(N_BDRY+1), offs_par(N_BDRY+1)
    integer :: dims_chd(2,N_BDRY+1), dims_par(2,N_BDRY+1)
    integer :: bdry(JPLUS:IMINUS)

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do c = 1, N_CHDRN
       p_chd = dom%patch%elts(p_par+1)%children(c)
       if (p_chd == 0) cycle

       call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

       bdry = [en, en, st, st]

       do j = bdry(JMINUS)+1, PATCH_SIZE/2+bdry(JPLUS)
          j_chd = 2*(j-1)
          j_par = j-1 + chd_offs(2,c)

          do i = bdry(IMINUS)+1, PATCH_SIZE/2+bdry(IPLUS)
             i_chd = 2*(i-1)
             i_par = i-1 + chd_offs(1,c)

             call routine ( &
                  dom, i_par, j_par, i_chd, j_chd, zlev, &
                  offs_par, dims_par, offs_chd, dims_chd)
          end do
       end do
    end do
  end subroutine apply_interscale_to_patch22

  
  subroutine apply_interscale_to_patch3 (routine, dom, p_par, c, zlev, st, en)
    implicit none

    procedure(sub7)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p_par, c, zlev, st, en

    integer :: i, j
    integer :: i_chd, i_par, j_chd, j_par, p_chd
    integer :: offs_chd(N_BDRY+1), offs_par(N_BDRY+1)
    integer :: dims_chd(2,N_BDRY+1), dims_par(2,N_BDRY+1)
    integer :: bdry(JPLUS:IMINUS)
    logical :: inner_bdry(JPLUS:IMINUS)

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    p_chd = dom%patch%elts(p_par+1)%children(c)
    if (p_chd == 0) return

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd, inner_bdry)

    bdry = [en, en, st, st]

    where (inner_bdry)
       bdry = 0
    end where

    do j = bdry(JMINUS)+1, PATCH_SIZE/2+bdry(JPLUS)
       j_chd = 2*(j-1)
       j_par = j-1 + chd_offs(2,c)

       do i = bdry(IMINUS)+1, PATCH_SIZE/2+bdry(IPLUS)
          i_chd = 2*(i-1)
          i_par = i-1 + chd_offs(1,c)

          call routine ( &
               dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, &
               offs_par, dims_par, offs_chd, dims_chd)
       end do
    end do
  end subroutine apply_interscale_to_patch3

  
  subroutine apply_to_pole (routine, l, zlev, ival, to_all)
    implicit none

    procedure(sub8)     :: routine
    integer, intent(in) :: l, zlev, ival
    logical, intent(in) :: to_all

    integer :: d

    do d = 1, size(grid)
       call apply_to_pole_d (routine, grid(d), l, zlev, ival, to_all)
    end do
  end subroutine apply_to_pole

  
  subroutine apply_to_pole_d (routine, dom, l, zlev, ival, to_all)
    implicit none

    procedure(sub8)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: l, zlev, ival
    logical,      intent(in)    :: to_all

    integer :: c, l_cur, p, p_par
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    do c = SOUTHEAST, NORTHWEST, 2
       if (.not. to_all .and. .not. dom%pole_master(c/2-2)) cycle
       if (.not. dom%penta(c)) cycle
       if (dom%neigh(c) /= POLE) cycle

       p = 1

       do while (p > 0)
          p_par = p
          p = dom%patch%elts(p_par+1)%children(c-4)

          if (l /= NONE) then
             l_cur = dom%patch%elts(p_par+1)%level

             if (l_cur < l) then
                cycle
             else if (l_cur > l) then
                exit
             end if
          end if

          call get_offs_Domain (dom, p_par, offs, dims)

          select case (c)
          case (NORTHWEST)
             call routine (dom, p_par, 0, PATCH_SIZE, zlev, offs, dims, ival)

          case (SOUTHEAST)
             call routine (dom, p_par, PATCH_SIZE, 0, zlev, offs, dims, ival)
          end select
       end do
    end do
  end subroutine apply_to_pole_d

  
  subroutine apply_to_pole_patch (routine, dom, p_start, zlev)
    implicit none

    procedure(sub4)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p_start, zlev

    integer :: c, p, p_par
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    do c = SOUTHEAST, NORTHWEST, 2
       if (.not. dom%penta(c)) cycle
       if (dom%neigh(c) /= POLE) cycle

       p = p_start

       do while (p > 0)
          p_par = p
          p = dom%patch%elts(p_par+1)%children(c-4)

          call get_offs_Domain(dom, p_par, offs, dims)

          select case (c)
          case (NORTHWEST)
             call routine (dom, 0, PATCH_SIZE, zlev, offs, dims)

          case (SOUTHEAST)
             call routine (dom, PATCH_SIZE, 0, zlev, offs, dims)
          end select
       end do
    end do
  end subroutine apply_to_pole_patch

  
  subroutine apply_to_pole2 (routine, l, zlev)
    implicit none

    procedure(sub4)     :: routine
    integer, intent(in) :: l, zlev

    integer :: c, d, l_cur, p, p_par
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    do d = 1, size(grid)
       do c = SOUTHEAST, NORTHWEST, 2
          if (.not. grid(d)%penta(c)) cycle
          if (grid(d)%neigh(c) /= POLE) cycle

          p = 1

          do while (p > 0)
             p_par = p
             p = grid(d)%patch%elts(p_par+1)%children(c-4)

             if (l /= NONE) then
                l_cur = grid(d)%patch%elts(p_par+1)%level

                if (l_cur < l) then
                   cycle
                else if (l_cur > l) then
                   exit
                end if
             end if

             call get_offs_Domain (grid(d), p_par, offs, dims)

             select case (c)
             case (NORTHWEST)
                call routine (grid(d), 0, PATCH_SIZE, zlev, offs, dims)

             case (SOUTHEAST)
                call routine (grid(d), PATCH_SIZE, 0, zlev, offs, dims)
             end select
          end do
       end do
    end do
  end subroutine apply_to_pole2

  
  subroutine apply_to_penta (routine, l, zlev)
    implicit none

    procedure(sub9)     :: routine
    integer, intent(in) :: l, zlev

    integer :: d

    do d = 1, size(grid)
       call apply_to_penta_d (routine, grid(d), l, zlev)
    end do
  end subroutine apply_to_penta

  
  subroutine apply_to_penta_d (routine, dom, l, zlev)
    implicit none

    procedure(sub9)             :: routine
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: l, zlev

    integer :: c, l_cur, p, p_par
    integer :: offs(N_BDRY+1)
    integer :: dims(2,N_BDRY+1)

    do c = NORTHEAST, NORTHWEST
       if (.not. dom%penta(c)) cycle

       p = 1

       do while (p > 0)
          p_par = p
          p = dom%patch%elts(p_par+1)%children(c-4)

          if (l /= NONE) then
             l_cur = dom%patch%elts(p_par+1)%level

             if (l_cur < l) then
                cycle
             else if (l_cur > l) then
                exit
             end if
          end if

          call get_offs_Domain (dom, p_par, offs, dims)
          call routine (dom, p_par, c, offs, dims, zlev)
       end do
    end do
  end subroutine apply_to_penta_d

  
end module domain_ops_mod
