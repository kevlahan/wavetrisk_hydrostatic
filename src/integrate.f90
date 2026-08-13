module integrate_mod
  
  ! Routines for integrating over hexagonal and triangular grid and for determing equivalent adaptive triangular grid
  
  use kind_mod,   only : dp
  use shared_mod, only : ADJZONE, AT_NODE, EDGE, N_BDRY, LORT, UPLT, RT, DG, UP, TRIAG, NORTHWEST, SOUTHEAST, &
       level_start, level_end, z_null

  use arch_mod,       only : rank
  use comm_mpi_mod,   only : sum_real
  use domain_mod,     only : Domain, Init_Field, Int_Field, Logical_Field, get_offs_Domain, grid, idx_hex_lort, idx_hex_uplt, idx
  use domain_ops_mod, only : apply_interscale, apply_onescale, fun3
  use dyn_arrays,     only : init
  use patch_mod,      only : PATCH_SIZE

  implicit none

  private
  public :: active_level, pre_levelout, post_levelout, save_tri
  public :: integrate_hex, integrate_tri

  
  real(dp)                                  :: integral
  type(Int_Field)                           :: active_level
  type(Logical_Field), dimension(LORT:UPLT) :: save_tri   


  procedure (fun3), pointer :: hex_fun => null ()

  
contains

  
  function integrate_hex (fun, zlev, level) result(val)
    ! Integrate over adaptive hexagons, where the hex_fun is defined by the routine fun.
    ! If optional variable coarse_only = .true. the integration is carried out over level_start only.
    implicit none
    integer, intent(in)           :: zlev
    integer, intent(in), optional :: level
    real(dp)                      :: val
    procedure(fun3)               :: fun

    hex_fun => fun

    integral = 0.0_dp

    if (present (level)) then ! integrate over a single scales
       call integrate_hex_scale (level, zlev)
    else                      ! integrate over coarsest scale
       call integrate_hex_scale (level_start, zlev)
    end if

    val = sum_real (integral)

    nullify (hex_fun)
  end function integrate_hex


  subroutine integrate_hex_scale (l, zlev)
    ! Integrate function pointer hex_fun over hexagons at a single scale l
    implicit none
    integer, intent(in) :: l, zlev

    integer :: c, d, i, id, j, jj, p
    integer :: offs(N_BDRY+1) 
    integer :: dims(2,N_BDRY+1)

    do d = 1, size(grid)
       ! Regular hexagons/pentagons
       do jj = 1, grid(d)%lev(l)%length
          p = grid(d)%lev(l)%elts(jj)
          call get_offs_Domain (grid(d), p, offs, dims)
          do j = 1, PATCH_SIZE
             do i = 1, PATCH_SIZE
                id = idx (i-1, j-1, offs, dims)
                integral = integral + hex_fun (grid(d), i-1, j-1, zlev, offs, dims) / grid(d)%areas%elts(id+1)%hex_inv
             end do
          end do
       end do

       ! Check domain d to see if its SOUTHEAST or NORTHWEST corners have associated poles
       do c = SOUTHEAST, NORTHWEST, 2
          if (.not. grid(d)%pole_master(c/2-2) .or. .not. grid(d)%penta(c)) cycle
          p = 1
          do while (grid(d)%patch%elts(p+1)%level < l)
             p = grid(d)%patch%elts(p+1)%children(c-4)
             if (p == 0) then
                write (6,'(A, i4, A)') "ERROR(rank = ", rank, "): integrate_hex: level incomplete"
                return
             end if
          end do
          call get_offs_Domain (grid(d), p, offs, dims)

          if (c == NORTHWEST) then     ! north pole
             id = idx (0, PATCH_SIZE, offs, dims)
             integral = integral + hex_fun (grid(d), 0, PATCH_SIZE, zlev, offs, dims) / grid(d)%areas%elts(id+1)%hex_inv
          elseif (c == SOUTHEAST) then ! south pole
             id = idx (PATCH_SIZE, 0, offs, dims)
             integral = integral + hex_fun (grid(d), PATCH_SIZE, 0, zlev, offs, dims) / grid(d)%areas%elts(id+1)%hex_inv
          end if

       end do
    end do
  end subroutine integrate_hex_scale


  function integrate_tri (fun, zlev) result(val)
    ! Integrate a function defined by fun over complete (non-overlapping) adaptive triangular grid.
    implicit none
    integer, intent(in)    :: zlev
    real(dp)               :: val
    procedure(fun3) :: fun
    
    integer :: l

    interface
     
    end interface
    
    hex_fun => fun

    call pre_levelout
    
    integral = 0.0_dp
    do l = level_start, level_end
       call apply_onescale (fdA_tri, l, zlev, 0, 1)
    end do

    val = sum_real (integral)
    nullify (hex_fun)

    call post_levelout
  end function integrate_tri

  
  subroutine fdA_tri (dom, i, j, zlev, offs, dims)
    ! Integrate over adaptive triangle grid.

    implicit none

    type (Domain),                  intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer                        :: d, id, idE, idN, idNE
    real(dp), dimension(LORT:UPLT) :: FdTri, tri_area
    real(dp), dimension(2*EDGE)    :: hex_area
    real(dp), dimension(0:EDGE)    :: val

    d = dom%id + 1

    id   = idx(i,   j,   offs, dims)

    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idN  = idx(i,   j+1, offs, dims)

    tri_area(LORT) = dom%triarea%elts(TRIAG*id+LORT+1)
    tri_area(UPLT) = dom%triarea%elts(TRIAG*id+UPLT+1)

    hex_area(1) = dom%areas%elts(id+1  )%part(1)
    hex_area(2) = dom%areas%elts(id+1  )%part(2)
    hex_area(3) = dom%areas%elts(idE+1 )%part(3)
    hex_area(4) = dom%areas%elts(idNE+1)%part(4)
    hex_area(5) = dom%areas%elts(idNE+1)%part(5)
    hex_area(6) = dom%areas%elts(idN+1 )%part(6)

    val(0)    = hex_fun (dom, i,   j,   zlev, offs, dims)
    val(RT+1) = hex_fun (dom, i+1, j,   zlev, offs, dims)
    val(DG+1) = hex_fun (dom, i+1, j+1, zlev, offs, dims)
    val(UP+1) = hex_fun (dom, i,   j+1, zlev, offs, dims)

    FdTri = hex2tri3 (val, hex_area, tri_area)  
    
    if (save_tri(LORT)%data(d)%elts(id+1)) integral = integral + FdTri(LORT)
    if (save_tri(UPLT)%data(d)%elts(id+1)) integral = integral + FdTri(UPLT)
  end subroutine fdA_tri


  function hex2tri3 (sclr, hex_area, tri_area) result(val)
    ! Integrates sclr given at hexagons over triangles

    implicit none

    real(dp), intent(in) :: sclr(0:EDGE)
    real(dp), intent(in) :: hex_area(2*EDGE)
    real(dp), intent(in) :: tri_area(LORT:UPLT)

    real(dp) :: val(LORT:UPLT) 

    val(LORT) = sclr(0) * hex_area(1) + sclr(1) * hex_area(3) + sclr(2) * hex_area(5)
    val(UPLT) = sclr(0) * hex_area(2) + sclr(2) * hex_area(4) + sclr(3) * hex_area(6) 
  end function hex2tri3


  subroutine pre_levelout
    implicit none
    integer :: d, l, num, t

    ! FIXME cleaner would be to use init_io routine
    call init_Field (active_level,   AT_NODE)
    call init_Field (save_tri(LORT), AT_NODE)
    call init_Field (save_tri(UPLT), AT_NODE)

    do d = 1, size(grid)
       num = grid(d)%node%length
       call init (active_level%data(d), num)

       do t = LORT, UPLT
          call init (save_tri(t)%data(d), num)
          save_tri(t)%data(d)%elts(1:num) = .false.
       end do

       active_level%data(d)%elts = grid(d)%level%elts
    end do

    do l = level_end-1, level_start-1, -1
       call apply_interscale (mark_save_chd, l, z_null, 0, 0)
    end do
    
    do l = level_start, level_end-1
       call apply_interscale (mark_save_par, l, z_null, 0, 0)
    end do
  end subroutine pre_levelout

  subroutine post_levelout
    implicit none
    integer :: d, t

    do d = 1, size(grid)
       deallocate (active_level%data(d)%elts)
       do t = LORT, UPLT
          deallocate (save_tri(t)%data(d)%elts)
       end do
    end do
    deallocate (active_level%data)

    do t = LORT, UPLT
       deallocate (save_tri(t)%data)
    end do
  end subroutine post_levelout

  
  subroutine mark_save_chd (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Mark active child triangles
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer                    :: d, id_chd, idE_chd, idNE_chd, idN_chd
    integer                    :: mask_LORT, mask_UPLT
    integer, dimension(2*EDGE) :: id_LORT, id_UPLT
    
    d = dom%id + 1

    id_chd   = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
   
    id_LORT = idx_hex_LORT (i_chd, j_chd, offs_chd, dims_chd)
    id_UPLT = idx_hex_UPLT (i_chd, j_chd, offs_chd, dims_chd)

    mask_LORT = minval (dom%mask_n%elts(id_LORT+1))
    if (mask_LORT >= ADJZONE) then
       save_tri(LORT)%data(d)%elts([id_chd, idE_chd, idNE_chd]+1) = .true.
       save_tri(UPLT)%data(d)%elts(idE_chd+1)                     = .true.
    end if

    mask_UPLT = minval (dom%mask_n%elts(id_UPLT+1))
    if (mask_UPLT >= ADJZONE) then
       save_tri(LORT)%data(d)%elts(idN_chd+1)                     = .true.
       save_tri(UPLT)%data(d)%elts([id_chd, idNE_chd, idN_chd]+1) = .true.
    end if
  end subroutine mark_save_chd

  
  subroutine mark_save_par (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer :: d, id_par, id_chd

    d = dom%id + 1

    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    if (save_tri(LORT)%data(d)%elts(id_chd+1)) save_tri(LORT)%data(d)%elts(id_par+1) = .false.
    if (save_tri(UPLT)%data(d)%elts(id_chd+1)) save_tri(UPLT)%data(d)%elts(id_par+1) = .false.
  end subroutine mark_save_par

  
end module integrate_mod
