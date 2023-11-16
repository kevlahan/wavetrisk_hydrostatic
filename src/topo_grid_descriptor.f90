module topo_grid_descriptor_mod
  use comm_mpi_mod
  use utils_mod
  use init_mod
  implicit none
#  include "netcdf.inc"
  !
  !  DATE CODED:  November 2023
  !  DESCRIPTION:
  !
  !  Generates smoothed topography data for WAVETRISK from NCAR topography NetCDF files using cube_to_target
  !  program that remaps topography data from cubed-sphere grid to target grid (non-adaptive WAVETRISK grid at max_level)
  !  grid using rigorous remapping (Lauritzen, Nair and Ullrich, 2010, J. Comput. Phys.).
  !
  !  Workflow
  !
  !  First step: pre-processing of coordinate data
  !
  !                               call write_grid_coords
  !
  !  to produce a grid descriptor file (e.g. J08_topo_grid.nc) for ESMF/SCRIP software in NetCDF file format
  !  for the hexagons on a given non-adaptive WAVETRISK grid (e.g. the grid corresponding to the desired max_level).
  !
  !  Second step:
  !
  !  Use the external program
  !
  !                                    cube_to_target
  !
  !  to generate the NetCDF file that provides the surface geopotential phi_S = z/g corresponding to the hexagons 
  !  saved in Step 1.  It is practical to use a script to specify appropriate parameters e.g. J08.sh.
  !                         
  !  Third step: each time you run the test case:
  !
  !                               call assign_height (fname)
  !
  !  to read in the .nc file generated in Step 2 to assign the topogragraphy data to the float field topography,
  !  which must have the same max_level and domain configuration as the WAVETRISK grid that generated the data in Step 1.
  !
  !  Author: Peter Hjort Lauritzen (pel@ucar.edu)
  !  Modified for inclusion WAVETRISK by Nicholas Kevlahan (kevlahan@mcmaster.ca) 2023-10
  !
contains
  subroutine write_grid_coords
    ! Find grid coordinates on each domain at finest level over entire grid
    ! saves results to netcdf coordinate file
    use mpi
    implicit none
    integer                               :: i, icol, loc_size
    integer                               :: grid_size          ! number of nodes
    integer, parameter                    :: grid_corners = 6   ! hexagons
    integer, dimension(:),    allocatable :: grid_dom, loc_dom
    integer, dimension(:),    allocatable :: grid_imask          
    integer, dimension(:),    allocatable :: grid_id, loc_id            
    real (8), dimension(:),   allocatable :: grid_area, loc_area          
    real (8), dimension(:),   allocatable :: grid_center_lat, loc_center_lat
    real (8), dimension(:),   allocatable :: grid_center_lon, loc_center_lon
    real (8), dimension(:,:), allocatable :: grid_corner_lat, loc_corner_lat
    real (8), dimension(:,:), allocatable :: grid_corner_lon, loc_corner_lon

    integer :: i_node

    character(255) :: grid_name

    loc_size = 0
    call apply_onescale (count_nodes, max_level, z_null, 0, 1)
    allocate (loc_dom(1:loc_size), loc_id(1:loc_size))
    allocate (loc_area(1:loc_size), loc_center_lat(1:loc_size), loc_center_lon(1:loc_size))
    allocate (loc_corner_lat(1:grid_corners,1:loc_size), loc_corner_lon(1:grid_corners,1:loc_size))

    icol = 0
    call apply_onescale (grid_coords, max_level, z_null, 0, 1)
#ifdef MPI
    grid_size = sum_int (loc_size)
#else
    grid_size = loc_size
#endif
    allocate (grid_imask(1:grid_size)); grid_imask = 1
    allocate (grid_dom(1:grid_size), grid_id(1:grid_size))
    allocate (grid_area(1:grid_size), grid_center_lat(1:grid_size), grid_center_lon(1:grid_size))
    allocate (grid_corner_lat(1:grid_corners,1:grid_size), grid_corner_lon(1:grid_corners,1:grid_size))

#ifdef MPI
    call MPI_Gather (loc_dom, loc_size, MPI_INTEGER, grid_dom, loc_size, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_Gather (loc_id,  loc_size, MPI_INTEGER, grid_id,  loc_size, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    call MPI_Gather (loc_area, loc_size, MPI_DOUBLE_PRECISION, grid_area, loc_size, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierror)

    call MPI_Gather (loc_center_lat, loc_size, MPI_DOUBLE_PRECISION, grid_center_lat, loc_size, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierror)

    call MPI_Gather (loc_center_lon, loc_size, MPI_DOUBLE_PRECISION, grid_center_lon, loc_size, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierror)

    do i = 1, grid_corners
       call MPI_Gather (loc_corner_lat(i,:), loc_size, MPI_DOUBLE_PRECISION, grid_corner_lat(i,:), loc_size, &
            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
       call MPI_Gather (loc_corner_lon(i,:), loc_size, MPI_DOUBLE_PRECISION, grid_corner_lon(i,:), loc_size, &
            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    end do
#else
    grid_dom = loc_dom; grid_id = loc_id; grid_area = loc_area; grid_center_lat = loc_center_lat; 
    grid_center_lon = loc_center_lon; grid_corner_lat = loc_corner_lat; grid_corner_lon = loc_corner_lon
#endif
    
    if (rank==0) then
       write (grid_name, '(a,i2.2,a)') "J", max_level, "_topo_grid"
       call wrt_esmf_rll
    end if
  contains
    subroutine count_nodes (dom, i, j, zlev, offs, dims)
      ! Count number of nodes
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      loc_size = loc_size + 1
    end subroutine count_nodes

    subroutine grid_coords (dom, i, j, zlev, offs, dims)
      ! Set grid coordinates
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id, id_i, idS, idSW, idW
      real(8) :: lat, lon

      d = dom%id + 1
      id = idx (i, j, offs, dims)
      id_i = id + 1

      idW  = idx (i-1, j,   offs, dims)
      idSW = idx (i-1, j-1, offs, dims)
      idS  = idx (i,   j-1, offs, dims)

      icol = icol + 1

      loc_id(icol)  = id                   ! id of node
      loc_dom(icol) = glo_id(rank+1,d) + 1 ! global domain label associated to node id

      loc_area(icol) = 1d0 / dom%areas%elts(id_i)%hex_inv / radius**2 ! hexagon area (unit sphere)

      call cart2sph (dom%node%elts(id_i), lon, lat) ! longitude and latitude coordinates of node in radians
      loc_center_lat(icol) = lat
      loc_center_lon(icol) = lon

      call cart2sph (dom%ccentre%elts(TRIAG*id+LORT+1), lon, lat)    
      loc_corner_lat(1,icol) = lat
      loc_corner_lon(1,icol) = lon

      call cart2sph (dom%ccentre%elts(TRIAG*id+UPLT+1), lon, lat)    
      loc_corner_lat(2,icol) = lat
      loc_corner_lon(2,icol) = lon

      call cart2sph (dom%ccentre%elts(TRIAG*idW+LORT+1), lon, lat)    
      loc_corner_lat(3,icol) = lat
      loc_corner_lon(3,icol) = lon

      call cart2sph (dom%ccentre%elts(TRIAG*idSW+UPLT+1), lon, lat)    
      loc_corner_lat(4,icol) = lat
      loc_corner_lon(4,icol) = lon
    
      call cart2sph (dom%ccentre%elts(TRIAG*idSW+LORT+1), lon, lat)    
      loc_corner_lat(5,icol) = lat
      loc_corner_lon(5,icol) = lon
    
      call cart2sph (dom%ccentre%elts(TRIAG*idS+UPLT+1), lon, lat)    
      loc_corner_lat(6,icol) = lat
      loc_corner_lon(6,icol) = lon
    end subroutine grid_coords
    
    subroutine wrt_esmf_rll
      !
      ! write netCDF grid descriptor file for ESMF remapping
      ! 
      implicit none
      !-----------------------------------------------------------------------
      !
      !     grid coordinates and masks
      !
      !-----------------------------------------------------------------------
      integer  :: ncstat             ! general netCDF status variable
      integer  :: nc_grid_id         ! netCDF grid dataset id
      integer  :: nc_lon_id          ! netCDF grid dataset id
      integer  :: nc_lat_id          ! netCDF grid dataset id
      integer  :: nc_gridsize_id     ! netCDF grid size dim id
      integer  :: nc_gridcorn_id     ! netCDF grid corner dim id
      integer  :: nc_gridrank_id     ! netCDF grid rank dim id
      integer  :: nc_griddims_id     ! netCDF grid dimension size id
      integer  :: nc_grddom_id       ! netCDF grid domain id
      integer  :: nc_grdid_id        ! netCDF grid id id
      integer  :: nc_grdarea_id      ! netCDF grid area id
      integer  :: nc_grdcntrlat_id   ! netCDF grid center lat id
      integer  :: nc_grdcntrlon_id   ! netCDF grid center lon id
      integer  :: nc_grdcrnrlat_id   ! netCDF grid corner lat id
      integer  :: nc_grdcrnrlon_id   ! netCDF grid corner lon id
      integer  :: nc_grdimask_id     ! netCDF grid mask id

      integer, dimension(2) :: nc_dims2_id ! netCDF dim id array for 2-d arrays
      integer, dimension(2) :: grid_dims
      integer               :: status      ! return value for error control of netcdf routine

      character (len=255)   :: file_out

      !-----------------------------------------------------------------------
      !
      !     set up attributes for netCDF file
      !
      !-----------------------------------------------------------------------
      !***
      !*** create netCDF dataset for this grid
      !***
      ncstat = nf_create (trim(grid_name)//".nc", NF_CLOBBER, nc_grid_id)
      call handle_err (ncstat)

      ncstat = nf_put_att_text (nc_grid_id, NF_GLOBAL, 'title', len_trim(grid_name), grid_name)
      call handle_err (ncstat)

      !***
      !*** define grid size dimension
      !***
      ncstat = nf_def_dim (nc_grid_id, 'grid_size', grid_size, nc_gridsize_id)
      call handle_err (ncstat)

      !***
      !*** define grid corner dimension
      !***
      ncstat = nf_def_dim (nc_grid_id, 'grid_corners', grid_corners, nc_gridcorn_id)
      call handle_err (ncstat)

      !***
      !*** define grid rank dimension
      !***
      ncstat = nf_def_dim (nc_grid_id, 'grid_rank', 1, nc_gridrank_id)
      call handle_err (ncstat)

      !***
      !*** define grid dimension size array
      !***
      ncstat = nf_def_var (nc_grid_id, 'grid_dims', NF_INT, 1, nc_gridrank_id, nc_griddims_id)
      call handle_err (ncstat)

      !***
      !*** define domain
      !***
      ncstat = nf_def_var (nc_grid_id, 'grid_dom', NF_INT, 1, nc_gridsize_id, nc_grddom_id)
      call handle_err (ncstat)

      !***
      !*** define id
      !***
      ncstat = nf_def_var (nc_grid_id, 'grid_id', NF_INT, 1, nc_gridsize_id, nc_grdid_id)
      call handle_err (ncstat)

      !***
      !*** define grid area array
      !***
      ncstat = nf_def_var (nc_grid_id, 'grid_area', NF_DOUBLE, 1, nc_gridsize_id, nc_grdarea_id)
      call handle_err (ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdarea_id, 'units', 9, 'radians^2')
      call handle_err (ncstat)

      !***
      !*** define grid center latitude array
      !***
      ncstat = nf_def_var (nc_grid_id, 'grid_center_lat', NF_DOUBLE, 1, nc_gridsize_id, nc_grdcntrlat_id)
      call handle_err (ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcntrlat_id, 'units', 7, 'radians')
      call handle_err (ncstat)

      !***
      !*** define grid center longitude array
      !***
      ncstat = nf_def_var (nc_grid_id, 'grid_center_lon', NF_DOUBLE, 1, nc_gridsize_id, nc_grdcntrlon_id)
      call handle_err (ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcntrlon_id, 'units', 7, 'radians')
      call handle_err (ncstat)

      !***
      !*** define grid mask
      !***
      ncstat = nf_def_var (nc_grid_id, 'grid_imask', NF_INT, 1, nc_gridsize_id, nc_grdimask_id)
      call handle_err (ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdimask_id, 'units', 8, 'unitless')
      call handle_err (ncstat)

      !***
      !*** define grid corner latitude array
      !***
      nc_dims2_id(1) = nc_gridcorn_id
      nc_dims2_id(2) = nc_gridsize_id

      ncstat = nf_def_var (nc_grid_id, 'grid_corner_lat', NF_DOUBLE, 2, nc_dims2_id(1), nc_grdcrnrlat_id)
      call handle_err(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcrnrlat_id, 'units', 7, 'radians')
      call handle_err (ncstat)

      !***
      !*** define grid corner longitude array
      !***
      ncstat = nf_def_var (nc_grid_id, 'grid_corner_lon', NF_DOUBLE, 2, nc_dims2_id(1), nc_grdcrnrlon_id)
      call handle_err (ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcrnrlon_id, 'units', 7, 'radians')
      call handle_err (ncstat)

      !***
      !*** end definition stage
      !***

      ncstat = nf_enddef (nc_grid_id)
      call handle_err (ncstat)

      !-----------------------------------------------------------------------
      !
      !     write grid data
      !
      !-----------------------------------------------------------------------
      write (6,'(/,a,a,a)',advance="no") "Writing NCAR topography grid descriptor file ", trim(grid_name)//".nc ..."

      ncstat = nf_put_var_int (nc_grid_id, nc_griddims_id, grid_dims)
      call handle_err (ncstat)

      ncstat = nf_put_var_int (nc_grid_id, nc_grdimask_id, grid_imask)
      call handle_err (ncstat)

      ncstat = nf_put_var_int (nc_grid_id, nc_grddom_id, grid_dom)
      call handle_err (ncstat)

      ncstat = nf_put_var_int (nc_grid_id, nc_grdid_id, grid_id)
      call handle_err (ncstat)

      ncstat = nf_put_var_double (nc_grid_id, nc_grdarea_id, grid_area)
      call handle_err (ncstat)

      ncstat = nf_put_var_double (nc_grid_id, nc_grdcntrlat_id, grid_center_lat)
      call handle_err (ncstat)

      ncstat = nf_put_var_double (nc_grid_id, nc_grdcntrlon_id, grid_center_lon)
      call handle_err (ncstat)

      ncstat = nf_put_var_double (nc_grid_id, nc_grdcrnrlat_id, grid_corner_lat)
      call handle_err (ncstat)

      ncstat = nf_put_var_double (nc_grid_id, nc_grdcrnrlon_id, grid_corner_lon)
      call handle_err (ncstat)

      ncstat = nf_close(nc_grid_id)
      call handle_err (ncstat)

      write (6, '(a)') " finished writing file descriptor data file."
      write (6, '(a,/)') "Use cube_to_target to generate corresponding surface geopotential .nc file."
    end subroutine wrt_esmf_rll
  end subroutine write_grid_coords

  subroutine assign_height (fname)
    ! Reads netcdf geopotential data and saves it to the variable topography
    implicit none
    character(*), intent(in) :: fname

    integer                            :: d, ncol
    integer, dimension(:), allocatable :: grid_dom, grid_id
    real(8), dimension(:), allocatable :: phi_s

    if (rank == 0) write (6,'(/,a,a,a)',advance="no") "Reading NCAR topography data ", trim(fname)//".nc ..."
    
    call read_geopotential

    do d = 1, size(grid)
       call apply_onescale_d (cal_assign_height, grid(d), max_level, z_null, 0, 1)
    end do

    topography%bdry_uptodate = .false.
    call update_bdry (topography, max_level)

    if (rank == 0) write (6,'(a,/)') "finished reading topography data."
  contains
    subroutine read_geopotential
      ! Reads netcdf geopotential data onto a single core and saves the data as its wavelet coefficients
#ifdef MPI
      use mpi
#endif
      implicit none
#     include         <netcdf.inc>

      integer :: dimid, domid, idxid, ierror, ncid, status,  phisid
      !********************************************
      !
      ! Get dimension on rank 0
      !
      !********************************************
      if (rank==0) then
         status = nf_open(TRIM(fname)//".nc", 0, ncid)
         if (STATUS /= NF_NOERR) call handle_err (status) 
         status = NF_INQ_DIMID (ncid, 'ncol', dimid)
         if (status /= NF_NOERR) call handle_err (status)
         status = NF_INQ_DIMLEN (ncid, dimid, ncol)
         if (status /=  NF_NOERR) call handle_err (status)
      end if
#ifdef MPI
      call MPI_Bcast (ncol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif
      allocate (grid_dom(1:ncol), grid_id(1:ncol), phi_s(1:ncol))

      !********************************************
      !
      !  Read variables on rank 0
      !   - all are dimension(ncol)
      !
      !********************************************
      if (rank == 0) then
         status = NF_INQ_VARID (ncid, 'dom', domid)
         if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
         status = NF_GET_VAR_INT (ncid, domid, grid_dom)
         if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)

         status = NF_INQ_VARID (ncid, 'idx', idxid)
         if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
         status = NF_GET_VAR_INT (ncid, idxid, grid_id)
         if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)

         status = NF_INQ_VARID (ncid, 'PHIS', phisid)
         if (status /= NF_NOERR) call handle_err (status)
         status = NF_GET_VAR_DOUBLE (ncid, phisid, phi_s)
         if (status /= NF_NOERR) call handle_err (status)
      end if
#ifdef MPI
      call MPI_Bcast (grid_dom, ncol, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast (grid_id,  ncol, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast (phi_s,    ncol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
#endif
    end subroutine read_geopotential

    subroutine cal_assign_height (dom, i, j, zlev, offs, dims)
      ! Assign topography to multiscale float field variable topography
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: id, d_glo, icol

      d_glo = glo_id(rank+1,d) + 1 ! global domain
      do icol = 1, ncol
         if (grid_dom(icol) == d_glo) then ! data is for this domain
            id = grid_id(icol) + 1
            topography%data(d)%elts(id) = phi_s(icol) / grav_accel
         end if
      end do
    end subroutine cal_assign_height
  end subroutine assign_height

  subroutine handle_err (status)
    implicit         none
    integer          status

    if (status /= nf_noerr) then
       print *, nf_strerror(status)
       stop 'Stopped'
    endif

  end subroutine handle_err
end module topo_grid_descriptor_mod
