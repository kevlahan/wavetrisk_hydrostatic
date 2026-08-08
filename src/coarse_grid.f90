module coarse_grid_mod
  ! Routines for generating coarsest grid using Xu (2006) Laplacian algorithm or by reading Heikes and Randall (1995) grids
  !
  ! Xu (Int J Comput J , 16(1), 75-93 2006) produces a spherical triangulation that is optimal for the Laplace-Beltrami operator in
  ! the sense that the approximation of the discrete mean curvature is exact and the error of the Laplace-Beltrami operator is minimal.
  !
  ! Heikes and Randall (Monthly Weather Rev 123, 1881-1887) provides spherical triangulations that are consistent for the Laplace and
  ! flux-divergence operators on the sphere. Laplace, Jacobian and flux-divergence operators have approximately second-order accuracy
  !
  ! Both schemes produce grids with similar l2 errors for the mis-match between the midpoints of primal and dual grid edges, but
  ! the linf error is smaller for the Heikes and Randall scheme.

  use kind_mod,   only : dp
  use shared_mod, only : BDRY_THICKNESS, Coord, DOMAIN_LEVEL, EDGE, N_BDRY, N_ICOSAH_LOZENGE, &
       N_SUB_DOM, N_SUB_DOM_PER_DIM, PATCH_LEVEL, nghb_pt, NONE, ORIGIN, &
       RT, DG, UP, LORT, UPLT, TRIAG, SOUTHWEST, level_end, level_start, min_level, grid_type, z_null, theta_grid, dx_avg
  
  use arch_mod,       only : barrier, loc_id, n_process, owner, rank
  use comm_mod,       only : get_coord, set_coord
  use comm_mpi_mod,   only : comm_nodes3_mpi, sum_real, sync_max_real
  use domain_mod,     only : Domain, get_offs_domain, grid, idx, idx2
  use domain_ops_mod, only : apply_onescale, apply_onescale2
  use geom_mod,       only : dist, init_Coord, inner, mid_pt, norm, number_hex, project_on_sphere, vector
  use init_mod,       only : ccentre, ccentre_penta, check_grid, midpt

  use coord_arithmetic_mod
  
  implicit none

  private
  public :: read_optim_grid, smooth_Xu, update_geom_check_grid, zrotate
  
  integer                                  :: next_fid = 100
  integer,     dimension(2,4), parameter   :: HR_offs = reshape ( [0,0, 1,0, 1,1, 0,1], [2,4] ) 
  real(dp)                                 :: dx_coarse, linf_err, l2_err, linf_err_loc, l2_err_loc
  type(coord), dimension(:,:), allocatable :: new_node

  
contains

  
  subroutine read_optim_grid  
    ! Reads in optimized grid from directory grids
    ! Need to provide a symbolic link to grids directory in working directory
    
    implicit none
    
    integer        :: d_glo, d_HR, d_sub, fid, loz, p, r
    integer        :: offs(N_BDRY+1)
    integer        :: dims(2,N_BDRY+1)
    character(999) :: filename

    dx_coarse = dx_avg (min_level-1) ! average edge lengths

    ! Initial error
    call grid_error
    if (rank == 0) then
       write (6,'(a)') '-------------------------------------------------------&
            &---------------------------------------------------------------------------'
       write (6,'(a,2(es8.2,a))') 'Grid quality of non-optimized grid  = ', linf_err, ' (linf) ', l2_err, ' (l2)'

    end if

    ! Read optimized grid
    fid = get_fid()
    if (level_start /= level_end) then
       write (0,'(i2,1x,i2)') level_end, level_start
       write (0,'(a)') "Reading optimized vertices for level_start not equal to level_end not implemented"
       return
    end if

    do r = 1, n_process

       if (r /= rank+1) then ! read only if our turn, otherwise wait at barrier
          call barrier
          cycle 
       end if

       write (filename, '(a,a,a,i3.3)')  "grids/", trim (grid_type), "_WT_", level_start-1
       open (unit = fid, file = trim(filename), status = "old", form = "formatted", access = "stream", action = "read")

       p = 1
       do d_HR = 1, N_ICOSAH_LOZENGE
          loz = dom_id_from_HR_id (d_HR)
          do d_sub = 1, N_SUB_DOM
             d_glo = loz * N_SUB_DOM + sub_dom_id_from_HR_sub_id (d_sub)
             if (owner(d_glo+1) == rank) call get_offs_Domain (grid(loc_id(d_glo+1)+1), p, offs, dims)
             
             call coord_from_file (d_glo, PATCH_LEVEL, fid, offs, dims, [ 0, 0])
          end do
       end do
       close (fid)
    end do

    call update_geom_check_grid

    ! Final error
    call grid_error
    if (rank == 0) then
       write (6,'(a,a,a,2(es8.2,a))') 'Grid quality of ', trim (filename), ' = ', linf_err, ' (linf) ', l2_err, ' (l2)'
       write (6,'(a)') '(relative distance between midpoints of primal and dual grid edges compared to average edge length)'
       write (6,'(a)') '-------------------------------------------------------&
            &---------------------------------------------------------------------------'
    end if
  end subroutine read_optim_grid

  
  subroutine smooth_Xu
    
    implicit none
    
    real(dp) :: tol

    dx_coarse = dx_avg (min_level-1)
    allocate (new_node(maxval(grid(:)%node%length), size(grid)))
    
    ! Initial error
    call grid_error
    if (rank == 0) then
       write (6,'(a)') '-------------------------------------------------------&
            &---------------------------------------------------------------------------'
       write (6,'(a,i2/)') 'Xu (2006) diffusion optimization of level ', level_end-1
       write (6,'(a,2(es8.2,a))') 'Grid quality before optimization = ', linf_err, ' (linf) ', l2_err, ' (l2)'
    end if

    tol = 1e-6_dp
    do while (linf_err > tol)
       linf_err_loc = 0.0_dp
       call apply_onescale (Xu_smooth_cpt,    min_level-1, z_null, 0, 0)
       call apply_onescale (Xu_smooth_assign, min_level-1, z_null, 0, 0)
       linf_err = sync_max_real (linf_err_loc)
       
       call update_geom_check_grid
    end do
  
    ! Final error
    call grid_error
    if (rank == 0) then
       write (6,'(a,2(es8.2,a))') 'Grid quality after optimization  = ', linf_err, ' (linf) ', l2_err, ' (l2)'
       write (6,'(a)') '(relative distance between midpoints of primal and dual grid edges compared to average edge length)'
       write (6,'(a,/)') '-------------------------------------------------&
            &----------------------------------------------------------------------'
    end if
    deallocate (new_node)

    ! Rotate grid
    call apply_onescale (rotate_grid, min_level-1, z_null,  0, 0)
  end subroutine smooth_Xu

  
  recursive subroutine coord_from_file (d_glo, l, fid, offs, dims, ij0)
    
    implicit none

    integer, intent(in) :: d_glo, l, fid
    integer, intent(in) :: ij0(2)
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer     :: d_loc, id, k
    integer     :: ij(2)
    type(Coord) :: node, node_rot

    d_loc = loc_id(d_glo+1)

    do k = 1, 4
       ij = ij0 + HR_offs(:,k)*2**(l-1)
       id = idx(ij(1), ij(2), offs, dims)

       if (l == 1) then
          read(fid,*) node

          call zrotate (node, node_rot, theta_grid)

          if (owner(d_glo+1) == rank) then
             grid(d_loc+1)%node%elts(id+1) = project_on_sphere (node_rot)
          end if
       else
          call coord_from_file (d_glo, l-1, fid, offs, dims, ij)
       end if
    end do
  end subroutine coord_from_file

  
  subroutine update_geom_check_grid
    
    implicit none
    
    integer :: d

    ! Communicate nodes
    call comm_nodes3_mpi (get_coord, set_coord)

    ! Update geometry
    call apply_onescale2 (ccentre, min_level-1, z_null, -2, 1)
    do d = 1, size(grid)
       call ccentre_penta (grid(d), 1)
    end do
    call apply_onescale2 (midpt,   min_level-1, z_null, -1, 2)

    ! Fix grid
    call apply_onescale2 (check_grid, min_level-1, z_null,  0, 0)

    ! Communicate new nodes
    call comm_nodes3_mpi (get_coord, set_coord)   
  end subroutine update_geom_check_grid

  
  subroutine Xu_smooth_cpt (dom, i, j, zlev, offs, dims)
    ! Algorithm 1 of Xu (2006)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in) :: i, j, zlev
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)
    
    integer     :: d, id, n
    real(dp)    :: alpha, beta, cosalpha, cosbeta
    type(Coord) :: s, p_i, p_ip, p_im, p_j, v1, v2
    
    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1
    
    p_i = dom%node%elts(id)

    ! Do not move pentagon nodes 
    if (i == 0 .and. j == 0 .and. dom%penta(SOUTHWEST)) then 
       new_node(id,d) = p_i 
       return
    end if

    call init_Coord (s, 0.0_dp, 0.0_dp, 0.0_dp)
    
    do n = 1, 6 ! sum over hexagon vertices 
       p_j  = dom%node%elts(idx2(i, j, nghb_pt(:,n),                  offs, dims) + 1)
       p_ip = dom%node%elts(idx2(i, j, nghb_pt(:,modulo(n,   6) + 1), offs, dims) + 1)
       p_im = dom%node%elts(idx2(i, j, nghb_pt(:,modulo(n-2, 6) + 1), offs, dims) + 1)
       
       v1 = vector (p_im, p_j)
       v2 = vector (p_im, p_i)
       cosalpha = inner (v1, v2) / (norm (v1) * norm (v2))
       alpha = acos (cosalpha)
       
       v1 = vector (p_ip, p_j)
       v2 = vector (p_ip, p_i)
       cosbeta = inner (v1, v2) / (norm (v1) * norm (v2))
       beta  = acos (cosbeta)

       s = s + (1/tan(alpha) + 1/tan(beta)) * (p_i - p_j)
    end do

    new_node(id,d) = project_on_sphere (s)
  end subroutine Xu_smooth_cpt

  
  subroutine Xu_smooth_assign (dom, i, j, zlev, offs, dims)
    ! Update node position
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in) :: i, j, zlev
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)
    
    integer :: d, id

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    linf_err_loc = max (linf_err_loc, dist (dom%node%elts(id), new_node(id,d))/dx_coarse)

    dom%node%elts(id) = new_node(id,d)
  end subroutine Xu_smooth_assign

  
  subroutine rotate_grid (dom, i, j, zlev, offs, dims)
    ! Rotate entire grid
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in) :: i, j, zlev
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)
    
    integer :: id

    id = idx (i, j, offs, dims) + 1

    call zrotate (dom%node%elts(id), dom%node%elts(id), theta_grid) 
  end subroutine rotate_grid

  
  subroutine check_d (dom, i, j, zlev, offs, dims)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    
    integer                :: id, idS, idW
    real(dp), dimension(3) :: error

    id  = idx(i,   j,   offs, dims)
    idS = idx(i,   j-1, offs, dims)
    idW = idx(i-1, j,   offs, dims)

    error = [ &
         dist (dom%midpt%elts(EDGE*id+RT+1), mid_pt (dom%ccentre%elts(TRIAG*id +LORT+1), dom%ccentre%elts(TRIAG*idS+UPLT+1))), &
         dist (dom%midpt%elts(EDGE*id+DG+1), mid_pt (dom%ccentre%elts(TRIAG*id +LORT+1), dom%ccentre%elts(TRIAG*id +UPLT+1))), &
         dist (dom%midpt%elts(EDGE*id+UP+1), mid_pt (dom%ccentre%elts(TRIAG*idW+LORT+1), dom%ccentre%elts(TRIAG*id +UPLT+1)) ) ]

    linf_err_loc = max (linf_err_loc, maxval (error))
    l2_err_loc = l2_err_loc + sum (error**2)
  end subroutine check_d
  
  function dom_id_from_HR_id (d_HR) result(val)
    ! d_HR: lozenge id as used by Heikes & Randall (starts from 1)
    ! results: domain id (starts from 0)
    
    implicit none
    
    integer, intent(in) :: d_HR
    integer             :: val

    val = modulo (d_HR, 2) * 5 + modulo (d_HR/2 - 1, 5)
  end function dom_id_from_HR_id

  function sub_dom_id_from_HR_sub_id (sub_id) result(val)
    ! sub_id: lozenge sub id as used by Heikes & Randall (starts from 1)
    ! results: sub domain id (starts from 0)

    implicit none

    integer, intent(in) :: sub_id
    integer             :: val

    integer :: id, i, j, halv_sub_dom, l, jdiv, idiv

    i = 0
    j = 0
    id = sub_id - 1
    halv_sub_dom = N_SUB_DOM/2
    
    do l = DOMAIN_LEVEL-1, 0, -1
       jdiv = id/halv_sub_dom
       j = j + jdiv*2**l
       id = modulo (id+4**l,4**(l+1))
       idiv = id/halv_sub_dom
       i = i + idiv*2**l
       halv_sub_dom = halv_sub_dom/4
       id = modulo (id,4**l)
    end do

    val = j * N_SUB_DOM_PER_DIM + i
  end function sub_dom_id_from_HR_sub_id

  function get_fid () result(val)
    
    implicit none
    
    integer :: val

    val  = next_fid
    next_fid = next_fid + 1
  end function get_fid

  subroutine grid_error
    ! Computes error
    
    implicit none

    call comm_nodes3_mpi (get_coord, set_coord)
    call apply_onescale2 (ccentre, level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)
    call apply_onescale2 (midpt,   level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)

    l2_err_loc   = 0.0_dp; linf_err_loc = 0.0_dp
    call  apply_onescale (check_d, level_end-1, z_null, 0, 0)
    l2_err   = sqrt (sum_real (l2_err_loc)) / (dx_coarse * 3 * number_hex (level_end-1))
    linf_err = sync_max_real (linf_err_loc) / dx_coarse
  end subroutine grid_error

  subroutine zrotate (c_in, c_out, angle)
    
    implicit none

    type(Coord), intent(in)  :: c_in
    type(Coord), intent(out) :: c_out
    real(dp),    intent(in)  :: angle

    c_out%x = c_in%x * cos(angle) - c_in%y * sin(angle)
    c_out%y = c_in%x * sin(angle) + c_in%y * cos(angle)
    c_out%z = c_in%z
  end subroutine zrotate

  
end module coarse_grid_mod
