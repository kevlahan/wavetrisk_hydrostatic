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
  use shared_mod
  use geom_mod
  use domain_mod
  use init_mod
  use comm_mpi_mod
  implicit none
  integer                                  :: ncell
  integer                                  :: next_fid = 100
  integer,     dimension(2,4), parameter   :: HR_offs = reshape ( [0,0, 1,0, 1,1, 0,1], [2,4] ) 
  real(dp)                                 :: dx_coarse, linf_err, l2_err
  type(coord), dimension(:,:), allocatable :: new_node
contains
  subroutine read_optim_grid  
    ! Reads in optimized grid from directory grids
    ! Need to provide a symbolic link to grids directory in working directory
    implicit none
    integer                        :: d_glo, d_HR, d_sub, fid, loz, p, r
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(999)                 :: filename

    dx_coarse = dx_avg (min_level) ! average edge lengths

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
#ifdef MPI
       if (r /= rank+1) then ! read only if our turn, otherwise wait at barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if
#endif
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

    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call apply_onescale2 (ccentre,    level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)
    call apply_onescale2 (midpt,      level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)
    call apply_onescale2 (check_grid, level_end-1, z_null,  0, 0)

    ! Final error
    call grid_error
    if (rank == 0) then
       write (6,'(a,a,a,2(es8.2,a))') 'Grid quality of ', trim (filename), ' = ', linf_err, ' (linf) ', l2_err, ' (l2)'
       write (6,'(a)') '(relative distance between midpoints of primal and dual grid edges compared to average edge length)'
       write (6,'(a)') '-------------------------------------------------------&
            &---------------------------------------------------------------------------'
    end if
  end subroutine read_optim_grid

  recursive subroutine coord_from_file (d_glo, l, fid, offs, dims, ij0)
    implicit none
    integer,                        intent(in) :: d_glo, l, fid
    integer, dimension(2),          intent(in) :: ij0
    integer, dimension(N_BDRY+1),   intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims

    integer               :: d_loc, id, k
    integer, dimension(2) :: ij
    real(dp), parameter   :: theta = -0.5_dp ! rotate grid around pole (for backwards compatibility)
    type(Coord)           :: node

    d_loc = loc_id(d_glo+1)
    do k = 1, 4
       ij = ij0 + HR_offs(:,k) * 2**(l-1)
       id = idx (ij(1), ij(2), offs, dims) 
       if (l == 1) then
          read (fid,*) node
          call zrotate (node, node, theta) 
          if (owner(d_glo+1) == rank) grid(d_loc+1)%node%elts(id+1) = project_on_sphere (node)
       else
          call coord_from_file (d_glo, l-1, fid, offs, dims, ij)
       end if
    end do

    do k = 1, 12
       call zrotate (penta_node(k), penta_node(k), theta)
    end do
  end subroutine coord_from_file

  subroutine smooth_Xu
    implicit none
    real(dp) :: tol

    tol = 1e9_dp * eps (radius) ! tolerance in [m], about 1.4 m on Earth or relative error of O(1e-7)

    dx_coarse = dx_avg (min_level)
    
    ! Initial error
    call grid_error
    
    if (rank == 0) then
       write (6,'(a)') '-------------------------------------------------------&
            &---------------------------------------------------------------------------'
       write (6,'(a,i2/)') 'Xu (2006) diffusion optimization of level ', level_end-1
       write (6,'(a,2(es8.2,a))') 'Grid quality before optimization = ', linf_err, ' (linf) ', l2_err, ' (l2)'
    end if
        
    allocate (new_node(maxval(grid(:)%node%length), size(grid)))

    linf_err = 2*tol
    do while (linf_err > tol)
       call comm_nodes3_mpi (get_coord, set_coord, NONE)
       call apply_onescale2 (ccentre,   level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)
       call apply_onescale2 (midpt,     level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)
       call apply_onescale2 (cpt_areas, level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)

       call apply_onescale (Xu_smooth_cpt, level_end-1, z_null, 0, 0)

       linf_err = 0.0_dp
       call apply_onescale (Xu_smooth_assign, level_end-1, z_null, 0, 0)
       linf_err = sync_max_real (linf_err)
    end do
    
    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call apply_onescale2 (ccentre,    level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-2)
    call apply_onescale2 (midpt,      level_end-1, z_null, -(BDRY_THICKNESS-2), BDRY_THICKNESS-1)
    call apply_onescale2 (check_grid, level_end-1, z_null, 0, 0)

    ! Final error
    call grid_error
    if (rank == 0) then
       write (6,'(a,2(es8.2,a))') 'Grid quality after optimization  = ', linf_err, ' (linf) ', l2_err, ' (l2)'
       write (6,'(a)') '(relative distance between midpoints of primal and dual grid edges compared to average edge length)'
       write (6,'(a,/)') '-------------------------------------------------&
            &----------------------------------------------------------------------'
    end if
    deallocate (new_node)
  end subroutine smooth_Xu

  subroutine Xu_smooth_cpt (dom, i, j, zlev, offs, dims)
    ! Algorithm 1 of Xu (2006)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
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

       s = s + (1/tan(alpha) + 1/tan(beta))/4 * (p_i - p_j)
    end do

    new_node(id,d) = project_on_sphere (dom%areas%elts(id)%hex_inv * s) 
  end subroutine Xu_smooth_cpt

  subroutine Xu_smooth_assign (dom, i, j, zlev, offs, dims)
    ! Update node position
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer ::  d, id

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    linf_err = max  (linf_err, dist (dom%node%elts(id), new_node(id,d)))

    dom%node%elts(id) = new_node(id,d)
  end subroutine Xu_smooth_assign

  subroutine init_smooth_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_shared_mod
    call init_domain_mod
    initialized = .true.
  end subroutine init_smooth_mod

  subroutine check_grid (dom, p, i, j, zlev, offs, dims)
    implicit none
    integer                        :: i, j, p, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    type(Domain) :: dom
    integer      :: id, idE, idN, idNE, idS, idW

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    call check_triag (dom, TRIAG*id + LORT, [TRIAG*idE + UPLT, TRIAG*id  + UPLT, TRIAG*idS + UPLT], &
         [id, idE, idNE], [EDGE*idE+UP, EDGE*id+DG, EDGE*id+RT])
    call check_triag (dom, TRIAG*id + UPLT, [TRIAG*idN + LORT, TRIAG*idW + LORT, TRIAG*id  + LORT], &
         [id, idNE, idN], [EDGE*idN+RT, EDGE*id+UP, EDGE*id+DG])
  end subroutine check_grid

  subroutine check_triag (dom, id, id_neigh, id_cnr, id_side)
    implicit none
    type(Domain)          :: dom
    integer               :: id
    integer, dimension(3) :: id_neigh, id_cnr, id_side
    
    integer                   :: i
    type(Coord)               :: cc_fine
    type(Coord), dimension(3) :: inters_pt
    logical, dimension(3)     :: does_inters, troubles

    cc_fine = circumcentre (dom%midpt%elts(id_side(1)+1), dom%midpt%elts(id_side(3)+1), dom%midpt%elts(id_side(2)+1))

    do i = 1, 3
       call arc_inters (dom%ccentre%elts(id+1), dom%ccentre%elts(id_neigh(i)+1), &
            cc_fine, circumcentre (dom%node%elts(id_cnr(i)+1), dom%midpt%elts(id_side(O2(1,i))+1), &
            dom%midpt%elts(id_side(O2(2,i))+1)), &
            inters_pt(i), does_inters(i), troubles(i))
    end do

    if (any (does_inters) .or. any (troubles)) then
       dom%node%elts(id_cnr(1)+1)%x = dom%node%elts(id_cnr(1)+1)%x + 1.0d7*eps(radius)
       dom%node%elts(id_cnr(1)+1) = project_on_sphere(dom%node%elts(id_cnr(1)+1))
    end if
  end subroutine check_triag

  subroutine check_d (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer                :: id, idS, idW
    real(dp), dimension(3) :: error

    id  = idx(i,   j,   offs, dims)
    idS = idx(i,   j-1, offs, dims)
    idW = idx(i-1, j,   offs, dims)

    error = [ &
         dist (dom%midpt%elts(EDGE*id+RT+1), mid_pt (dom%ccentre%elts(TRIAG*id +LORT+1), dom%ccentre%elts(TRIAG*idS+UPLT+1))), &
         dist (dom%midpt%elts(EDGE*id+DG+1), mid_pt (dom%ccentre%elts(TRIAG*id +LORT+1), dom%ccentre%elts(TRIAG*id +UPLT+1))), &
         dist (dom%midpt%elts(EDGE*id+UP+1), mid_pt (dom%ccentre%elts(TRIAG*idW+LORT+1), dom%ccentre%elts(TRIAG*id +UPLT+1)) ) ]

    linf_err = max (linf_err, maxval (error))
    l2_err = l2_err + sum (error**2)
  end subroutine check_d
  
  integer function dom_id_from_HR_id (d_HR)
    ! d_HR: lozenge id as used by Heikes & Randall (starts from 1)
    ! results: domain id (starts from 0)
    implicit none
    integer :: d_HR

    dom_id_from_HR_id = modulo (d_HR, 2) * 5 + modulo (d_HR/2 - 1, 5)
  end function dom_id_from_HR_id

  integer function sub_dom_id_from_HR_sub_id (sub_id)
    ! sub_id: lozenge sub id as used by Heikes & Randall (starts from 1)
    ! results: sub domain id (starts from 0)
    implicit none
    integer :: sub_id

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
    sub_dom_id_from_HR_sub_id = j*N_SUB_DOM_PER_DIM + i
  end function sub_dom_id_from_HR_sub_id

   integer function get_fid ()
    implicit none

    get_fid  = next_fid
    next_fid = next_fid + 1
  end function get_fid

  subroutine grid_error
    ! Computes error
    implicit none

    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call apply_onescale2 (ccentre, level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)
    call apply_onescale2 (midpt,   level_end-1, z_null, -(BDRY_THICKNESS-1), BDRY_THICKNESS-1)

    l2_err   = 0.0_dp; linf_err = 0.0_dp
    call  apply_onescale (check_d, level_end-1, z_null, 0, 0)
    l2_err   = sqrt (sum_real (l2_err)) / (dx_coarse * 3 * number_hex (level_end-1))
    linf_err = sync_max_real (linf_err) / dx_coarse
  end subroutine grid_error

  subroutine zrotate (c_in, c_out, angle)
    ! Rotates a point by longitude angle around pole
    ! (used to rotate entire grid)
    implicit none
    real(dp),    intent(in)  :: angle
    type(Coord), intent(in)  :: c_in
    type(Coord), intent(inout) :: c_out

    c_out%x =  c_in%x * cos (angle) - c_in%y * sin (angle)
    c_out%y =  c_in%x * sin (angle) + c_in%y * cos (angle)
    c_out%z =  c_in%z
  end subroutine zrotate
end module coarse_grid_mod
