module io_vtk_mod
  use domain_mod
  use ops_mod
  use utils_mod
  use multi_level_mod
  implicit none
  integer                                :: ncell, ncoord,  nvertex
  integer                                :: ncell_loc, ncoord_unique_loc, nvertex_unique_loc
  integer,     dimension(:), allocatable :: cell_vert_index,  cell_vert_index_loc
  real(4),     dimension(:), allocatable :: cell_data_loc, vert_coord_unique_loc
  real(4),     dimension(:), allocatable :: cell_data, vert_coord_unique
  type(coord), dimension(:), allocatable :: points_loc
    
  integer, parameter                     :: nvar = 12 ! number of variables to save
contains
  subroutine write_and_export (isave)
    implicit none
    integer :: isave

    integer         :: d, j, k, l
    character(3)    :: layer
    character(1000) :: file_hex, isv, file_tri, file_vtk, scale

    if (rank == 0) then
       write (6,'(a)') '**************************************************************&
            ********************************************************************'
       write (6,'(a,i4)') 'Saving field ', isave
    end if
    
    write (isv, '(i4.4)') isave

    sol%bdry_uptodate = .false.
    call update_bdry (sol, NONE, 883)

    call pre_levelout

    ! Compute surface pressure
    call cal_surf_press (sol(1:N_VARIABLE,1:zlevels))

    ! Compute vertical velocity (stored in trend(S_TEMP,:))
    if (compressible) then
       call omega_velocity    ! vertical velocity in pressure coordinates D_t P = OMEGA [Pa/s]
    else
       call vertical_velocity ! vertical velocity w [m/s]
    end if

    ! if (compressible .or. zlevels /= 2 .or. .not. mode_split) then
    do k = 1, 12
       ! Compute pressure, exner and geopotential
       do d = 1, size(grid)
          mass    =>      sol(S_MASS,k)%data(d)%elts
          temp    =>      sol(S_TEMP,k)%data(d)%elts
          mean_m  => sol_mean(S_MASS,k)%data(d)%elts
          mean_t  => sol_mean(S_TEMP,k)%data(d)%elts
          exner   =>       exner_fun(k)%data(d)%elts
          call apply_d (integrate_pressure_up, grid(d), k, 0, 1)
          nullify (mass, mean_m, temp, mean_t, exner)
       end do

       ! Zonal and meridional velocities 
       do d = 1, size(grid)
          velo  => sol(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          call apply_d (interp_UVW_latlon, grid(d), z_null, 0, 1)
          nullify (velo, velo1, velo2)
       end do

       ! Vorticity
       do d = 1, size(grid)
          velo  => sol(S_VELO,k)%data(d)%elts
          vort  => grid(d)%vort%elts
          call apply_d (cal_vort, grid(d), z_null, -1, 1)
          call apply_to_penta_d (post_vort, grid(d), l, z_null)
          nullify (velo, vort)
       end do

       ! Calculate vorticity at hexagon points (stored in ke)
       do d = 1, size(grid)
          vort => grid(d)%ke%elts
          call apply_d (vort_triag_to_hex, grid(d), k, 0, 1)
          nullify (vort)
       end do
    end do
    
    ! Save data
    call find_vertices
    file_vtk = trim(run_id) // "_tri_" // trim(isv) // ".vtk"
    if (rank == 0) then
       call write_vtk     (file_vtk)
       call compress_file (file_vtk)
    end if
    deallocate (cell_data, cell_vert_index, cell_vert_index_loc, cell_data_loc, vert_coord_unique)

     if (rank == 0) then
       if (write_tri_data) write (6,'(2(a,i8),a,f6.1)') &
            "Number of active cells = ", ncell, " number of unique vertices = ", nvertex, &
            " compression ratio = ", dble (2 * (2 + 10 * 4**max_level)) / dble (ncell)

       write (6,'(a)') '*************************************************************&
            *********************************************************************'
    end if

    call post_levelout
    call barrier
  end subroutine write_and_export
  
  subroutine find_vertices
    ! Find all unique triangle cell vertices on adaptive grid
    use mpi
    implicit none
    integer                       :: i, ibeg, iend, r
    integer, dimension(n_process) :: displs, ncell_glo, ncoord_unique_glo, nvertex_unique_glo, ncell_vert_index_glo

    allocate (cell_data_loc(0), cell_vert_index_loc(0), points_loc(0), vert_coord_unique_loc(0))

    ncoord_unique_loc = 0; nvertex_unique_loc = 0
    
    ncell = 0; ncell_loc = 0

    call apply_bdry (unique_tri_points, z_null, 0, 0)

    ncell   = sum_int (ncell_loc)          ! number of active cells
    ncoord  = sum_int (ncoord_unique_loc)  ! number of active unique vertex coordinates
    nvertex = sum_int (nvertex_unique_loc) ! number of unique vertices

    call gather_int (nvertex_unique_loc, nvertex_unique_glo,       displs)
    call gather_vec (ncoord_unique_loc,  ncoord_unique_glo,        vert_coord_unique_loc, vert_coord_unique)
    call gather_vec (3*ncell_loc,        ncell_vert_index_glo,     cell_vert_index_loc,   cell_vert_index)
    call gather_vec (nvar*ncell_loc,     ncell_glo, cell_data_loc, cell_data)

     ! Shift cell vertex indices
    if (rank == 0) then
       do r = 2, n_process
          do i = sum (ncell_vert_index_glo(1:r-1)) + 1, sum (ncell_vert_index_glo(1:r))
             cell_vert_index(i) = cell_vert_index(i) + sum (nvertex_unique_glo(1:r-1))
          end do
       end do
    end if

    deallocate (points_loc, vert_coord_unique_loc)
  end subroutine find_vertices

  subroutine unique_tri_points (dom, i, j, zlev, offs, dims)
    ! Finds all unique triangle vertices
    use utils_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                                  :: d, id, imin, ivert, t
    integer                                  :: idE, idN, idNE
    integer, dimension(1:EDGE)               :: new_vert_index
    integer, dimension(0:EDGE)               :: neigh_id
    integer, dimension(:), allocatable       :: tmp
    real(8)                                  :: dmin

    type(coord)                              :: p
    type(coord), dimension(LORT:UPLT,1:EDGE) :: vertex
    
    real(4), dimension(nvar)                 :: outv
    
    d    = dom%id + 1
    
    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx( i+1, j+1, offs, dims)

    vertex(LORT,:) = dom%node%elts((/id, idE, idNE/)+1)
    vertex(UPLT,:) = dom%node%elts((/id, idNE, idN/)+1)

    do t = LORT, UPLT
       if (save_tri(t)%data(d)%elts(id+1) == 1d0) then ! cell is active
          ncell_loc = ncell_loc + 1
          
          do ivert = 1, 3
             p = vertex(t,ivert)
             call min_dist (p, points_loc, dmin, imin)

             if (dmin > dx_min/4) then ! add vertex if is not already present
                nvertex_unique_loc = nvertex_unique_loc + 1
                ncoord_unique_loc  = ncoord_unique_loc  + 3
                
                points_loc     = [points_loc,     p]
                vert_coord_unique_loc = [vert_coord_unique_loc, real(p%x), real(p%y), real(p%z)]
                
                new_vert_index(ivert) = nvertex_unique_loc - 1 ! index of new vertex
             else                       
                new_vert_index(ivert) = imin - 1               ! index of existing vertex

             end if
          end do

          ! Add to cell vertices array
          cell_vert_index_loc = [cell_vert_index_loc, new_vert_index]

          ! Compute cell data
          call compute_data (12)

          ! Add to cell data array
          cell_data_loc = [cell_data_loc, outv]
       end if
    end do
  contains
    subroutine compute_data (k)
      use utils_mod
      implicit none
      integer :: k
      
      integer, dimension(0:EDGE)    :: neigh_id
      real(4)                       :: z
      real(4), dimension(0:EDGE)    :: rho_dz, rho_dz_theta
      real(4), dimension(0:EDGE)    :: temperature
      real(4)                       :: Ps, tri_area
      real(4), dimension(2*EDGE)    :: hex_area
      real(4), dimension(LORT:UPLT) :: relvort

      z = z_i (dom, i, j, zlev, offs, dims, sol)

      neigh_id = (/ id, idE, idNE, idN /) + 1

      tri_area = dom%triarea%elts(TRIAG*id+t+1)

      hex_area(1) = dom%areas%elts(id+1  )%part(1)
      hex_area(2) = dom%areas%elts(id+1  )%part(2)
      hex_area(3) = dom%areas%elts(idE+1 )%part(3)
      hex_area(4) = dom%areas%elts(idNE+1)%part(4)
      hex_area(5) = dom%areas%elts(idNE+1)%part(5)
      hex_area(6) = dom%areas%elts(idN+1 )%part(6)

      relvort = get_vort (dom, i, j, offs, dims)

      rho_dz       = sol(S_MASS,k)%data(d)%elts(neigh_id) + sol_mean(S_MASS,k)%data(d)%elts(neigh_id)
      rho_dz_theta = sol(S_TEMP,k)%data(d)%elts(neigh_id) + sol_mean(S_TEMP,k)%data(d)%elts(neigh_id)

      if (compressible) then
         temperature = rho_dz_theta/rho_dz * (dom%press%elts(neigh_id)/p_0)**kappa
      else
         temperature = ref_density * (1d0 - rho_dz_theta/rho_dz)
      end if

      Ps = hex2tri2 (real(dom%surf_press%elts(neigh_id)), hex_area, tri_area, t)                ! surface pressure
      
      outv(1) = hex2tri2 (rho_dz,                                       hex_area, tri_area, t)  ! rho_dz
      outv(2) = hex2tri2 (temperature,                                  hex_area, tri_area, t)  ! temperature (compressible) or density (incompressible)
      outv(3) = hex2tri2 (real(dom%u_zonal%elts(neigh_id)),             hex_area, tri_area, t)  ! zonal velocity
      outv(4) = hex2tri2 (real(dom%v_merid%elts(neigh_id)),             hex_area, tri_area, t)  ! meridional velocity
      outv(5) = hex2tri2 (real(trend(S_TEMP,k)%data(d)%elts(neigh_id)), hex_area, tri_area, t)  ! vertical velocity OMEGA (stored in trend(S_TEMP,:))
      outv(6) = relvort(t)                                                                      ! vorticity

      if (k == 1) then
         outv(7) = hex2tri2 (real(topography%data(d)%elts(neigh_id)),    hex_area, tri_area, t) ! topography
         outv(8) = hex2tri2 (real(penal_node(k)%data(d)%elts(neigh_id)), hex_area, tri_area, t) ! penalization mask
      else
         outv(7:8) = 0.0
      end if

      if (compressible) then
         outv(9) = Ps       
      else                                                                                 
         if (mode_split) then ! free surface perturbation
            outv(9) = hex2tri2 (real(sol(S_MASS,zlevels+1)%data(d)%elts(neigh_id)), hex_area, tri_area, t) 
         else
            outv(9) = hex2tri2 (real(sol(S_MASS,1)%data(d)%elts(neigh_id)),         hex_area, tri_area, t) 
         end if
      end if

      outv(10) = hex2tri2 (real(dom%geopot%elts(neigh_id) / grav_accel), hex_area, tri_area, t)      ! geopotential height
      outv(11) = hex2tri2 (real(dom%press%elts(neigh_id)),               hex_area, tri_area, t) / Ps ! P/Ps

      outv(12) =  nint (active_level%data(d)%elts(id+1))
    end subroutine compute_data
  end subroutine unique_tri_points
  
  subroutine write_vtk (filename)
    ! Each rank writes to vtk file
    implicit none
    character(*) :: filename

    integer                   :: icell, ivar
    integer, dimension(ncell) :: vtk_type
    character(12)             :: str1, str2, str3
    type(coord)               :: p
    
    integer, parameter :: funit        = 300
    integer, parameter :: nvert        = 3 ! triangular cells
    integer, parameter :: VTK_TRIANGLE = 5
    integer, parameter :: VTK_WEDGE    = 13
    real(8), parameter :: dr           = 100 * KM

    vert_coord_unique = shift_vertices (1000*KM)

    vtk_type = VTK_TRIANGLE

    write (str1(1:12),'(i12)') nvertex
    write (str2(1:12),'(i12)') ncell
    write (str3(1:12),'(i12)') ncell * (1 + nvert)
    
    open (unit=funit, file=trim(filename), form="unformatted", access='stream', status='replace', convert='BIG_ENDIAN')
    
    ! Write vtk header
    write (funit) '# vtk DataFile Version 2.0'//lf
    write (funit) 'WAVETRISK adaptive data output'//lf              
    write (funit) 'BINARY'//lf                   
    write (funit) 'DATASET UNSTRUCTURED_GRID'//lf

    ! Write coordinates of unique vertices
    write (funit) 'POINTS ' // trim(str1) // ' float' // lf
    write (funit) vert_coord_unique
    
    ! Write cell vertex indices (refers to POINTS)
    write (funit) 'CELLS ' // trim(str2) // trim(str3) // lf
    do icell = 1, ncell
       write (funit) nvert, cell_vert_index(3*(icell-1) + 1 : 3*(icell-1) + 3)
    end do

    ! Write out type of each cell
    write (funit) 'CELL_TYPES '//trim(str2)//lf
    write (funit) vtk_type

    ! Write out cell data
    write (funit) 'CELL_DATA ' // trim(str2) // lf

    write (funit) 'SCALARS rho_dz float' // lf
    write (funit) 'LOOKUP_TABLE default' // lf
    ivar = 1; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    if (compressible) then
       write (funit) 'SCALARS Temperature float'//lf
    else
       write (funit) 'SCALARS Density float'//lf
    end if
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 2; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS velocity_zonal float'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 3; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS velocity_meridional float'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 4; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS OMEGA float'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 5; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS Vorticity float'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 6; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS Topography float'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 7; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS penalization float'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 8; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    if (compressible) then
       write (funit) 'SCALARS Ps float'//lf
    else
       write (funit) 'SCALARS eta float'//lf
    end if
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 9; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS geopot_height float'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 10; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS P_Ps float'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 11; write (funit) cell_data(ivar : nvar*(ncell-1) + ivar : nvar)

    write (funit) 'SCALARS Level int'//lf
    write (funit) 'LOOKUP_TABLE default'//lf
    ivar = 12; write (funit) int (cell_data(ivar : nvar*(ncell-1) + ivar : nvar))

    close (funit)
  end subroutine write_vtk

  function shift_vertices (dr)
    ! Shifts all vertices radially by dr
    implicit none
    real(8)                    :: dr
    real(4), dimension(ncoord) :: shift_vertices
    
    integer     :: i
    real(8)     :: nrm, r
    type(coord) :: p

    do i = 1, ncoord, 3
       p = coord (vert_coord_unique(i), vert_coord_unique(i+1), vert_coord_unique(i+2))
       
       nrm = sqrt (p%x**2 + p%y**2 + p%z**2)

       r = (radius + dr) / nrm

       p%x = p%x * r
       p%y = p%y * r
       p%z = p%z * r

       shift_vertices(i:i+2) = [real(p%x), real(p%y), real(p%z)]
    end do
  end function shift_vertices

  function get_vort (dom, i, j, offs, dims)
    ! Averages vorticity to get smooth field for visualization
    implicit none
    real(8), dimension(TRIAG)      :: get_vort
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idN, idE, idS, idW

    id  = idx(i,   j,   offs, dims)
    idE = idx(i+1, j,   offs, dims)
    idN = idx(i,   j+1, offs, dims)
    idW = idx(i-1, j,   offs, dims)
    idS = idx(i,   j-1, offs, dims)

    get_vort(UPLT+1) = ( &
         interp(dom%vort%elts(TRIAG*idW+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)+ &
         
         interp(dom%vort%elts(TRIAG*id+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(id*EDGE+DG+1)+ &
         
         interp(dom%vort%elts(TRIAG*idN+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*idN+RT+1)*dom%pedlen%elts(EDGE*idN+RT+1)) &
         
         / (dom%len%elts(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)+ &
         dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)+ &
         dom%len%elts(EDGE*idN+RT+1)*dom%pedlen%elts(EDGE*idN+RT+1))

    get_vort(LORT+1) = ( &
         interp(dom%vort%elts(TRIAG*idS+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)+ &
         
         interp(dom%vort%elts(TRIAG*id+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(id*EDGE+DG+1)+ &
         
         interp(dom%vort%elts(TRIAG*idE+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*idE+UP+1)*dom%pedlen%elts(EDGE*idE+UP+1)) &
         
         / (dom%len%elts(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)+ &
         dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)+ &
         dom%len%elts(EDGE*idE+UP+1)*dom%pedlen%elts(EDGE*idE+UP+1))
  end function get_vort

  subroutine barotropic_velocity (dom, i, j, zlev, offs, dims)
    ! Calculate barotropic velocity in two-layer model
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                       :: d, e, id, id_e, id_i, idE, idNE, idN, k
    real(8)                       :: dz0
    real(8), dimension (1:EDGE,2) :: dz

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d = dom%id + 1

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    do k = 1, 2
       dz0 = sol_mean(S_MASS,k)%data(d)%elts(id+1) + sol(S_MASS,k)%data(d)%elts(id_i)
       dz(RT+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idE+1)  + sol(S_MASS,k)%data(d)%elts(idE+1))
       dz(DG+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idNE+1) + sol(S_MASS,k)%data(d)%elts(idNE+1))
       dz(UP+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idN+1)  + sol(S_MASS,k)%data(d)%elts(idN+1))
    end do

    do e = 1, EDGE
       id_e = EDGE*id + e
       velo(id_e) = (dz(e,1)*sol(S_VELO,1)%data(d)%elts(id_e) + dz(e,2)*sol(S_VELO,2)%data(d)%elts(id_e)) / (dz(e,1) + dz(e,2))
    end do
  end subroutine barotropic_velocity

  subroutine baroclinic_velocity (dom, i, j, zlev, offs, dims)
    ! Calculate baroclinic velocity in top layer
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, id_e

    id = idx (i, j, offs, dims)
    d = dom%id + 1

    do e = 1, EDGE
       id_e = EDGE*id + e
       velo(id_e) = velo2(id_e) - velo1(id_e)
    end do
  end subroutine baroclinic_velocity

  subroutine OMEGA_velocity
    ! Computes vertical velocity in pressure coordinates D_t P = OMEGA [Pa/s]
    ! stored in trend(S_TEMP,1:zlevels)
    ! note that OMEGA > 0 corresponds to negative vertical velocity (w < 0)
    implicit none
    integer :: d, j, k, l, p

    call update_bdry (sol, NONE)

    ! Compute surface pressure
    call cal_surf_press (sol(1:N_VARIABLE,1:zlevels))

    do k = 1, zlevels
       do l = level_end, level_start, -1
          ! Divergence of mass flux at each vertical level
          ! stored in trend(S_MASS,1:zlevels)
          do d = 1, size(grid)
             mass   =>      sol(S_MASS,k)%data(d)%elts
             velo   =>      sol(S_VELO,k)%data(d)%elts
             mean_m => sol_mean(S_MASS,k)%data(d)%elts
             h_flux => horiz_flux(S_MASS)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=7)
             end do
             nullify (mass, mean_m, velo)
             if (l < level_end) then
                dscalar => trend(S_MASS,k)%data(d)%elts
                call cpt_or_restr_flux (grid(d), l)
                nullify (dscalar)
             end if
             nullify (h_flux)
          end do
          horiz_flux(S_MASS)%bdry_uptodate = .false.
          call update_bdry (horiz_flux(S_MASS), l)
          do d = 1, size(grid)
             dscalar => trend(S_MASS,k)%data(d)%elts
             h_flux  => horiz_flux(S_MASS)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (dscalar, h_flux)
          end do

          ! u grad(P_S) at hexagon centres
          ! stored in trend(S_TEMP,1:zlevels)
          do d = 1, size(grid)
             scalar =>         grid(d)%surf_press%elts
             velo   =>      sol(S_VELO,k)%data(d)%elts
             h_flux => horiz_flux(S_TEMP)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=6)
             end do
             if (l < level_end) then
                dscalar => trend(S_TEMP,k)%data(d)%elts
                call cpt_or_restr_flux (grid(d), l)
                nullify (dscalar)
             end if
             nullify (h_flux, scalar, velo)
          end do
          horiz_flux(S_TEMP)%bdry_uptodate = .false.
          call update_bdry (horiz_flux(S_TEMP), l)
          do d = 1, size(grid)
             dscalar =>    trend(S_TEMP,k)%data(d)%elts
             h_flux  => horiz_flux(S_TEMP)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (dscalar, h_flux)
          end do

          trend(S_MASS:S_TEMP,k)%bdry_uptodate = .false.
          call update_bdry (trend(S_MASS:S_TEMP,k), l)
       end do
    end do

    ! Compute OMEGA
    call apply_bdry (cal_omega, z_null, 0, 1)
    
    trend(S_TEMP,1:zlevels)%bdry_uptodate = .false.
    call update_bdry (trend(S_TEMP,1:zlevels), NONE)
  end subroutine omega_velocity

  subroutine cal_omega (dom, i, j, zlev, offs, dims)
    ! Velocity flux across interfaces
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id_i, k
    real(8), dimension(1:zlevels)   :: u_gradP
    real(8), dimension(1:zlevels+1) :: div_mass

    d    = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    ! Vertically integrate div(mass flux) from top to bottom
    ! (results at interfaces)
    div_mass(zlevels+1) = 0d0 ! zero flux at top boundary
    do k = zlevels, 1, -1
       div_mass(k) = div_mass(k+1) + trend(S_MASS,k)%data(d)%elts(id_i)
    end do

    ! u.gradP at layers
    do k = 1, zlevels
       u_gradP(k) = interp (b_vert(k), b_vert(k+1)) * trend(S_TEMP,k)%data(d)%elts(id_i)
    end do

    ! Complete computation of OMEGA
    do k = 1, zlevels
       trend(S_TEMP,k)%data(d)%elts(id_i) = - grav_accel * interp (div_mass(k), div_mass(k+1)) + u_gradP(k) 
    end do
  end subroutine cal_omega

  subroutine vertical_velocity
    ! Computes vertical velocity w [m/s]
    ! stored in trend(S_TEMP,1:zlevels)
    implicit none
    integer :: d, j, k, l, p

    call omega_velocity

    call apply_bdry (cal_w, z_null, 0, 1)

    trend(S_TEMP,1:zlevels)%bdry_uptodate = .false.
    call update_bdry (trend(S_TEMP,1:zlevels), NONE)
  end subroutine vertical_velocity

  subroutine cal_w (dom, i, j, zlev, offs, dims)
    ! Vertical velocity w = - OMEGA / (rho_0 g) + (vertical projection of horizontal velocity)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, k

    d   = dom%id + 1
    id  = idx (i, j, offs, dims)

    do k = 1, zlevels
       trend(S_TEMP,k)%data(d)%elts(id+1) = - trend(S_TEMP,k)%data(d)%elts(id+1) / (ref_density * grav_accel) + proj_vel_vertical () 
    end do
  contains
    real(8) function proj_vel_vertical ()
      ! Computes grad_zonal(z) * u_zonal + grad_merid(z) * u_merid at hexagon centres for vertical velocity computation.
      ! Uses Perot formula as also used for kinetic energy:
      ! u = sum ( u.edge_normal * hexagon_edge_length * (edge_midpoint-hexagon_center) ) / cell_area
      implicit none
      integer :: idE, idN, idNE, idS, idSW, idW

      idE  = idx (i+1, j,   offs, dims)
      idNE = idx (i+1, j+1, offs, dims)
      idN  = idx (i,   j+1, offs, dims)
      
      idW  = idx (i-1, j,   offs, dims)
      idSW = idx (i-1, j-1, offs, dims)
      idS  = idx (i,   j-1, offs, dims)

      velo => sol(S_VELO,k)%data(d)%elts

      proj_vel_vertical =  &
           (vert_vel (i,j,i+1,j,EDGE*id +RT+1) + vert_vel (i+1,j+1,i,j,EDGE*id  +DG+1) + vert_vel (i,j,i,j+1,EDGE*id +UP+1) + &
           (vert_vel (i-1,j,i,j,EDGE*idW+RT+1) + vert_vel (i,j,i-1,j-1,EDGE*idSW+DG+1) + vert_vel (i,j-1,i,j,EDGE*idS+UP+1))) / 6d0

      nullify (velo)
    end function proj_vel_vertical

    real(8) function vert_vel (i1, j1, i2, j2, id_e)
      implicit none
      integer :: i1, j1, i2, j2, id_e

      real(8) :: dl, dz

      dz =  z_i (dom, i2, j2, k, offs, dims, sol) - z_i (dom, i1, j1, k, offs, dims, sol)
      dl = dom%len%elts(id_e)

      vert_vel = dz / sqrt (dl**2 + dz**2) * velo(id_e)
    end function vert_vel
  end subroutine cal_w

  subroutine compress_file (filename)
    implicit none
    character(300) :: filename

    integer         :: info
    character(1300) :: bash_cmd, command
    
    command = 'bash -c "gtar caf ' // trim(filename) // '.tgz '// trim(filename) // '"'
    call system (trim(command), info)
    
    if (info /= 0) then
       write (6,'(a)') 'gtar error info = 0 ... aborting'
       call abort
    end if
  end subroutine compress_file
end module io_vtk_mod
