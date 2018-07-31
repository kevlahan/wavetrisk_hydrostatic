module smooth_mod
  use shared_mod
  use geom_mod
  use domain_mod
  use init_mod
  use comm_mpi_mod
  implicit none
  real(8)                                  :: maxerror, l2error
  type(Coord), dimension(:,:), allocatable :: sums
contains
  subroutine init_smooth_mod
    implicit none
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_shared_mod
    call init_sphere_mod
    call init_domain_mod
    initialized = .True.
  end subroutine init_smooth_mod

  subroutine Xu_smooth_cpt (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer     :: n
    real(8)     :: alpha, beta, cosalpha, cosbeta, t
    type(Coord) :: s, p_i, p_ip, p_im, p_j, v, v1, v2
    
    call init_Coord (s, 0.0_8, 0.0_8, 0.0_8)

    p_i = dom%node%elts(idx(i, j, offs, dims) + 1)

    if (i == 0 .and. j == 0 .and. dom%penta(SOUTHWEST)) then ! pentagon
       sums(idx(i,j,offs,dims)+1,dom%id+1) = p_i !project_on_sphere(s)
       return
    end if

    do n = 1, 6
       p_j = dom%node%elts(idx2(i, j, nghb_pt(:,n), offs, dims) + 1)
       p_ip = dom%node%elts(idx2(i, j, nghb_pt(:,modulo(n, 6) + 1), offs, &
            dims) + 1)
       p_im = dom%node%elts(idx2(i, j, nghb_pt(:,modulo(n - 2, 6) + 1), &
            offs, dims) + 1)
       v1 = vector(p_im, p_j)
       v2 = vector(p_im, p_i)
       cosalpha = inner(v1, v2)/(norm(v1)*norm(v2))
       v1 = vector(p_ip, p_j)
       v2 = vector(p_ip, p_i)
       cosbeta = inner(v1, v2)/(norm(v1)*norm(v2))
       alpha = acos(cosalpha)
       beta = acos(cosbeta)
       v = vector(p_j, p_i)
       t = 1.0_8/tan(alpha) + (1.0_8/tan(beta))
       s%x = s%x + v%x*t
       s%y = s%y + v%y*t
       s%z = s%z + v%z*t
    end do
    sums(idx(i,j,offs,dims)+1,dom%id+1) = project_on_sphere(s)
  end subroutine Xu_smooth_cpt

  subroutine smooth_Xu (tol)
    implicit none
    real(8) :: tol
    
    integer :: k

    call comm_nodes3_mpi (get_coord, set_coord, NONE)

    call apply_onescale2 (ccentre, level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,   level_end-1, z_null, -1, 1)

    maxerror = 0.0_8
    l2error  = 0.0_8
    call  apply_onescale (check_d, level_end-1, z_null, 0, 0)

    l2error = sqrt(sum_real(l2error))
    maxerror = sync_max_d(maxerror)
    if (rank == 0) then
       write (6,'(A)') '-------------------------------------------------------&
            --------------------------------------------------------------------------'
       write (6,'(A,i2,A,es10.4/)') 'Xu (2006) diffusion optimization of level ', level_end-1, ' grid with tolerance ', tol
       write (6,'(A,2(es10.4,A))') 'Grid quality before optimization = ', maxerror, ' m (linf) ', l2error, ' m (l2)'
    end if
        
    allocate (sums(maxval(grid(:)%node%length), size(grid)))

    k = 0
    maxerror = 2.0_8*tol
    do while(maxerror > tol)
       maxerror = 0.0_8
       call comm_nodes3_mpi (get_coord, set_coord, NONE)

       call apply_onescale (Xu_smooth_cpt,    level_end-1, z_null, 0, 0)
       call apply_onescale (Xu_smooth_assign, level_end-1, z_null, 0, 0)
       maxerror = sync_max_d(maxerror)
       k = k + 1
    end do


    call comm_nodes3_mpi (get_coord, set_coord, NONE)

    call apply_onescale2 (ccentre,    level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,      level_end-1, z_null, -1, 1)
    call apply_onescale2 (check_grid, level_end-1, z_null,  0, 0)

    call comm_nodes3_mpi (get_coord, set_coord, NONE)

    call apply_onescale2 (ccentre, level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,   level_end-1, z_null, -1, 1)

    maxerror = 0.0_8
    l2error = 0.0_8
    call  apply_onescale (check_d, level_end-1, z_null, 0, 0)
    l2error = sqrt (sum_real(l2error))
    maxerror = sync_max_d(maxerror)

    if (rank == 0) then
       write (6,'(A,2(es10.4,A))') 'Grid quality after optimization  = ', maxerror, ' m (linf) ', l2error, ' m (l2)'
       write (6,'(A)') '(distance between midpoints of primal and dual edges)'
       write (6,'(A,/)') '-------------------------------------------------&
            --------------------------------------------------------------------'
    end if
    deallocate(sums)
  end subroutine smooth_Xu

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

    call check_triag (dom, id*TRIAG+LORT, (/TRIAG*idE+UPLT, TRIAG*id+UPLT, TRIAG*idS+UPLT/), &
         (/id, idE, idNE/), (/EDGE*idE+UP, EDGE*id+DG, EDGE*id+RT/))
    call check_triag (dom, id*TRIAG+UPLT, (/TRIAG*idN+LORT, TRIAG*idW+LORT, TRIAG*id+LORT/), &
         (/id, idNE, idN/), (/EDGE*idN+RT, EDGE*id+UP, EDGE*id+DG/))
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

    cc_fine = circumcentre(dom%midpt%elts(id_side(1)+1), dom%midpt%elts(id_side(3)+1), dom%midpt%elts(id_side(2)+1))

    do i = 1, 3
       call arc_inters(dom%ccentre%elts(id+1), dom%ccentre%elts(id_neigh(i)+1), &
            cc_fine, circumcentre(dom%node%elts(id_cnr(i)+1), dom%midpt%elts(id_side(O2(1,i))+1), &
            dom%midpt%elts(id_side(O2(2,i))+1)), &
            inters_pt(i), does_inters(i), troubles(i))
    end do

    if (any(does_inters) .or. any(troubles)) then
       dom%node%elts(id_cnr(1)+1)%x = dom%node%elts(id_cnr(1)+1)%x + 1.0d7*eps()
       dom%node%elts(id_cnr(1)+1) = project_on_sphere(dom%node%elts(id_cnr(1)+1))
    end if
  end subroutine check_triag

  subroutine Xu_smooth_assign (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer id

    id = idx(i, j, offs, dims)

    maxerror = max(maxerror, dist(dom%node%elts(id+1), sums(id+1,dom%id+1)))

    dom%node%elts(id+1) = sums(id+1,dom%id+1)
  end subroutine Xu_smooth_assign

  subroutine check_d (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer                :: id, idS, idW
    real(8), dimension (3) :: error

    id  = idx(i,   j,   offs, dims)
    idS = idx(i,   j-1, offs, dims)
    idW = idx(i-1, j,   offs, dims)

    error = (/dist(dom%midpt%elts(EDGE*id+RT+1), &
         mid_pt(dom%ccentre%elts(TRIAG*id+LORT+1),dom%ccentre%elts(TRIAG*idS+UPLT+1))), &
         dist(dom%midpt%elts(EDGE*id+DG+1), mid_pt(dom%ccentre%elts(TRIAG*id+LORT+1),dom%ccentre%elts(TRIAG*id+UPLT+1))), &
         dist(dom%midpt%elts(EDGE*id+UP+1), mid_pt(dom%ccentre%elts(TRIAG*idW+LORT+1),dom%ccentre%elts(TRIAG*id+UPLT+1)) )/)

    maxerror = max(maxerror, maxval(error))
    l2error = l2error + sum(error**2)
  end subroutine check_d
end module smooth_mod
