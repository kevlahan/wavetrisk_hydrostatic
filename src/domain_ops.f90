
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !! Subroutines for applying a routine to domains, patches, ...  !!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module domain_ops_mod
  use domain_mod
  implicit none
  interface
     real(8) function fun3 (dom, i, j, zlev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end function fun3
     function fun4 (dom, i, j, zlev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
       real(8), dimension(1:EDGE)     :: fun4
     end function fun4
     subroutine sub4 (dom, i, j, zlev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end subroutine sub4
     subroutine sub5 (dom, p, i, j, zlev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: p
       integer                        :: i, j, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end subroutine sub5
     subroutine sub6 (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i_par, j_par, i_chd, j_chd, zlev
       integer, dimension(N_BDRY+1)   :: offs_chd, offs_par
       integer, dimension(2,N_BDRY+1) :: dims_chd, dims_par
     end subroutine sub6
     subroutine sub7 (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev   
       integer, dimension(N_BDRY+1)   :: offs_chd, offs_par
       integer, dimension(2,N_BDRY+1) :: dims_chd, dims_par
     end subroutine sub7
     subroutine sub8 (dom, p, i, j, zlev, offs, dims, ival)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: ival, i, j, p, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end subroutine sub8
     subroutine sub9 (dom, p, c, offs, dims, zlev)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: p, c, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end subroutine sub9
      subroutine sub10 (dom, p, i, j, zlev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, p, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end subroutine sub10
  end interface
contains
  subroutine apply (routine, zlev)
    ! Applies routine over all levels and over entire boundary
    implicit none
    integer          :: zlev
    procedure (sub4) :: routine
    
    integer :: l

    do l = level_start, level_end
       call apply_onescale (routine, l, zlev, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
  end subroutine apply

  subroutine apply_bdry (routine, zlev, st, en)
    ! Applies routine to nodes/edges at all levels and including boundary cells specified by (st,en)
    implicit none
    integer          :: en, st, zlev
    procedure (sub4) :: routine
    
    integer :: l

    do l = level_start, level_end
       call apply_onescale (routine, l, zlev, st, en)
    end do
  end subroutine apply_bdry

  subroutine apply_onescale (routine, l, zlev, st, en)
     ! Applies routine to nodes/edges at all level l including boundary cells specified by (st,en)
    implicit none
    integer          :: en, l, st, zlev
    procedure (sub4) :: routine

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
    integer          :: en, ival, l, st, zlev
    procedure (sub8) :: routine

    integer :: d, j

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch__int (routine, grid(d), grid(d)%lev(l)%elts(j), zlev, st, en, ival)
       end do
    end do
  end subroutine apply_onescale__int

  subroutine apply_onescale_to_patch__int (routine, dom, p, zlev, st, en, ival)
    implicit none
    type(Domain)     :: dom
    integer          :: en, ival, p, st, zlev
    procedure (sub8) :: routine

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
  
  subroutine apply_no_bdry (routine, zlev)
    ! Applies routine to nodes and/or edges at all levels including poles but excluding all boundary cells
    implicit none
    integer          :: zlev 
    procedure (sub4) :: routine

    integer :: l

    do l = level_start, level_end
       call apply_onescale (routine, l, zlev, 0, 0)
       call apply_to_pole2 (routine, l, zlev)  
    end do
  end subroutine apply_no_bdry

  subroutine apply_no_bdry_to_patch (routine, dom, p, zlev)
    ! Applies routine to nodes and/or edges at all levels including poles but excluding all boundary cells
    implicit none
    type(Domain)     :: dom
    integer          :: p, zlev
    procedure (sub4) :: routine

    integer :: l

    l = dom%patch%elts(p+1)%level
    
    call apply_onescale_to_patch (routine, dom, p, zlev, 0, 0)
    call apply_to_pole2          (routine, l, zlev)
  end subroutine apply_no_bdry_to_patch

  subroutine apply_no_bdry2 (routine, zlev)
    ! Applies routine to nodes/edges at all levels excluding all boundary cells
    !
    ! The integer flag is_pole allows the routine to not update values for edges if is_pole = 1
    !
    ! Routine needs to accept an integer which indicates pole (if equals 1) or non-pole (if equals 0)
    !
    implicit none
    integer          :: zlev
    procedure (sub8) :: routine

    integer :: is_pole, l

    do l = level_start, level_end
       is_pole = 0; call apply_onescale__int (routine, l, zlev, 0, 0, is_pole)
       is_pole = 1; call apply_to_pole       (routine, l, zlev,       is_pole, .true.)  
    end do
  end subroutine apply_no_bdry2

  subroutine apply_d (routine, dom, zlev, st, en)
    implicit none
    type(Domain)     :: dom
    integer          :: en, st, zlev
    procedure (sub4) :: routine

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
    type(Domain)     :: dom
    integer          :: en, p, st, zlev
    procedure (sub4) :: routine

    integer                          :: i, j
    integer, dimension(N_BDRY+1)     :: offs
    integer, dimension(2,N_BDRY+1)   :: dims
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry

    call get_offs_Domain (dom, p, offs, dims, inner_bdry)

    bdry = [ en, en, st, st ]

    where (inner_bdry) bdry = 0

    do j = bdry(JMINUS) + 1, PATCH_SIZE + bdry(JPLUS)
       do i = bdry(IMINUS) + 1, PATCH_SIZE + bdry(IPLUS)
          call routine (dom, i - 1, j - 1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch

  subroutine apply_onescale_d (routine, dom, l, zlev, st, en)
    implicit none
    type(Domain)     :: dom
    integer          :: en, l, st, zlev
    procedure (sub4) :: routine

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_onescale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_onescale_d

  subroutine apply_onescale2 (routine, l, zlev, st, en)
    implicit none
    integer          :: en, l, st, zlev
    procedure (sub5) :: routine
    
    integer :: d, j

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch2 (routine, grid(d), grid(d)%lev(l)%elts(j), zlev, st, en)
       end do
    end do
  end subroutine apply_onescale2

  subroutine apply_onescale_to_patch5 (routine, dom, p, zlev, st, en)
    implicit none
    type(Domain)     :: dom
    integer          :: en, p, st, zlev
    procedure (sub5) :: routine

    integer                          :: i, j
    integer, dimension(N_BDRY+1)     :: offs
    integer, dimension(2,N_BDRY+1)   :: dims
    logical, dimension(JPlUS:IMINUS) :: inner_patch
    integer, dimension(JPlUS:IMINUS) :: bdry

    call get_offs_Domain5 (dom, p, offs, dims, inner_patch)

    bdry = [ en, en, st, st ]

    where (inner_patch) bdry = 0

    do j = bdry(JMINUS) + 1, PATCH_SIZE + bdry(JPLUS)
       do i = bdry(IMINUS) + 1, PATCH_SIZE + bdry(IPLUS)
          call routine (dom, p, i - 1, j - 1, zlev, offs, dims)
       end do
    end do
  end subroutine apply_onescale_to_patch5

  subroutine apply_onescale_to_patch2 (routine, dom, p, zlev, st, en)
    implicit none
    type(Domain)     :: dom
    integer          :: en, p, st, zlev
    procedure (sub5) :: routine

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
    integer          :: en, l, st, zlev
    procedure (sub6) :: routine
    
    integer :: d

    do d = 1, size(grid)
       call apply_interscale_d (routine, grid(d), l, zlev, st, en)
    end do
  end subroutine apply_interscale

  subroutine apply_interscale_d (routine, dom, l, zlev, st, en)
    implicit none
    type(Domain)     :: dom
    integer          :: en, l, st, zlev
    procedure (sub6) :: routine

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_interscale_to_patch (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_interscale_d
  
  subroutine apply_interscale_d2 (routine, dom, l, zlev, st, en)
    implicit none
    type(Domain)     :: dom
    integer          :: en, l, st, zlev
    procedure (sub6) :: routine

    integer :: j

    do j = 1, dom%lev(l)%length
       call apply_interscale_to_patch22 (routine, dom, dom%lev(l)%elts(j), zlev, st, en)
    end do
  end subroutine apply_interscale_d2

  subroutine apply_interscale_to_patch (routine, dom, p_par, zlev, st, en)
    implicit none
    type(Domain)     :: dom
    integer          :: p_par, en, st, zlev
    procedure (sub6) :: routine   

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

       bdry = [en, en, st, st]

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
    type(Domain)     :: dom
    integer          :: en, st, p_par, zlev
    procedure (sub7) :: routine

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

       bdry = [en, en, st, st]

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
    type(Domain)     :: dom
    integer          :: en, p_par, st, zlev
    procedure (sub6) :: routine
    
    integer                          :: c, i, i_chd, i_par, j, j_chd, j_par, p_chd
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    do c = 1, N_CHDRN
       p_chd = dom%patch%elts(p_par+1)%children(c)
       if (p_chd == 0) cycle

       call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

       bdry = [en, en, st, st]

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
    type(Domain)     :: dom
    integer          :: c, en, p_par, st, zlev
    procedure (sub7) :: routine
    
    integer                          :: i, i_chd, i_par, j, j_chd, j_par, p_chd
    integer, dimension(N_BDRY+1)     :: offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry
    !TODO{uncomment & test}  if (dom%patch%elts(p_par+1)%active == NONE) return

    call get_offs_Domain (dom, p_par, offs_par, dims_par)

    p_chd = dom%patch%elts(p_par+1)%children(c)

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd, inner_bdry)

    bdry = [en, en, st, st]

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
    integer          :: ival, l, zlev
    logical          :: to_all
    procedure (sub8) :: routine
    
    integer :: d

    do d = 1, size(grid)
       call apply_to_pole_d (routine, grid(d), l, zlev, ival, to_all)
    end do
  end subroutine apply_to_pole

  subroutine apply_to_pole_d (routine, dom, l, zlev, ival, to_all)
    implicit none
    type(Domain)     :: dom
    integer          :: ival, l, zlev
    logical          :: to_all
    procedure (sub8) :: routine

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

  subroutine apply_to_pole_patch (routine, dom, p_par, zlev)
    implicit none
    type(Domain)     :: dom
    integer          :: p_par, zlev
    procedure (sub4) :: routine

    integer                        :: c, p
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do c = SOUTHEAST, NORTHWEST, 2
       if (.not. dom%penta(c))         cycle
       if (.not. dom%neigh(c) == POLE) cycle
       p = 1
       do while (p > 0)
          p_par = p
          p     = dom%patch%elts(p_par+1)%children(c-4)

          call get_offs_Domain (dom, p_par, offs, dims)

          if (c == NORTHWEST) then     ! north pole
             call routine (dom, 0, PATCH_SIZE, zlev, offs, dims) 
          elseif (c == SOUTHEAST) then ! south pole
             call routine (dom, PATCH_SIZE, 0, zlev, offs, dims) 
          end if
       end do
    end do
  end subroutine apply_to_pole_patch

  subroutine apply_to_pole2 (routine, l, zlev)
    implicit none
    integer         :: l, zlev
    procedure(sub4) :: routine

    integer                        :: c, d, l_cur, p, p_par
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do d = 1, size(grid)
       do c = SOUTHEAST, NORTHWEST, 2
          if (.not. grid(d)%penta(c))         cycle
          if (.not. grid(d)%neigh(c) == POLE) cycle

          p = 1
          do while (p > 0)
             p_par = p
             p     = grid(d)%patch%elts(p_par+1)%children(c-4)

             if (.not. l == NONE) then
                l_cur = grid(d)%patch%elts(p_par+1)%level
                if (l_cur < l) then
                   cycle
                else
                   if (l_cur > l) exit
                end if
             end if

             call get_offs_Domain (grid(d), p_par, offs, dims)

             if (c == NORTHWEST) then     ! north pole
                call routine (grid(d), 0, PATCH_SIZE, zlev, offs, dims)
             elseif (c == SOUTHEAST) then ! south pole
                call routine (grid(d), PATCH_SIZE, 0, zlev, offs, dims) 
             end if

          end do
       end do
    end do
  end subroutine apply_to_pole2

  subroutine apply_to_penta (routine, l, zlev)
    implicit none
    integer          :: l, zlev
    procedure (sub9) :: routine

    integer :: d

    do d = 1, size(grid)
       call apply_to_penta_d(routine, grid(d), l, zlev)
    end do
  end subroutine apply_to_penta

  subroutine apply_to_penta_d (routine, dom, l, zlev)
    implicit none
    type(Domain)     :: dom
    integer          :: l, zlev
    procedure (sub9) :: routine

    integer                        :: c, l_cur, p, p_par
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
 end module domain_ops_mod
