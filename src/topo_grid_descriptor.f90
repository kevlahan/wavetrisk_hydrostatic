module topo_grid_descriptor_mod
  use comm_mpi_mod
  use utils_mod
  use init_mod
  implicit none
#  include "netcdf.inc"
  !
  !  DATE CODED:   January, 2012
  !  DESCRIPTION:  This program creates a grid descriptior file for ESMF/SCRIP software in NetCDF file format.
  !
  !  Author: Peter Hjort Lauritzen (pel@ucar.edu)
  !  Modified for inclusion WAVETRISK by Nicholas Kevlahan (kevlahan@mcmaster.ca) 2023-09-27
  !
  integer                             :: i_node, grid_size
  integer, parameter                  :: grid_corners = 6   ! hexagons
  integer, dimension(:),  allocatable :: grid_imask         ! not used

  real (8), dimension(:),   allocatable :: grid_area        ! hexagon area on unit sphere
  real (8), dimension(:),   allocatable :: grid_center_lat  ! lat/lon coordinates
  real (8), dimension(:),   allocatable :: grid_center_lon  ! each grid centre in radians
  real (8), dimension(:,:), allocatable :: grid_corner_lat  ! lat/lon coordinates 
  real (8), dimension(:,:), allocatable :: grid_corner_lon  ! each grid corner in radians
contains

  subroutine write_grid_coords
    ! Find grid coordinates on each domain at coarsest level over entire grid
    ! saves results to netcdf coordinate file
    implicit none
    integer        :: d, i, l
    character(255) :: grid_name

    l = min_level ! coarsest level

    grid_size = 0
    call apply_onescale (count_nodes, l, z_null, 0, 1)

    allocate (grid_imask(1:grid_size)); grid_imask = 1
    allocate (grid_area(1:grid_size), grid_center_lat(1:grid_size), grid_center_lon(1:grid_size))
    allocate (grid_corner_lat(1:grid_corners,1:grid_size), grid_corner_lon(1:grid_corners,1:grid_size))

    i_node = 0
    call apply_onescale (grid_coords, l, z_null, 0, 1) ! need to add poles

    write (grid_name, '(a,a,i3.3)') trim (run_id), "_topo_grid_", rank+1
    call wrt_esmf_rll (grid_name)
  end subroutine write_grid_coords

  subroutine count_nodes (dom, i, j, zlev, offs, dims)
    ! Count number of nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    grid_size = grid_size + 1
  end subroutine count_nodes

  subroutine grid_coords (dom, i, j, zlev, offs, dims)
    ! Set grid coordinates
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_i, idS, idSW, idW
    real(8) :: lat, lon

    id = idx (i, j, offs, dims)
    id_i = id + 1

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    
    i_node = i_node + 1

    grid_area(i_node) = 1d0/dom%areas%elts(id_i)%hex_inv/radius**2 ! hexagon area (unit sphere)
    call cart2sph (dom%node%elts(id_i), lon, lat) ! longitude and latitude coordinates of node in radians
    grid_center_lat(i_node) = lat
    grid_center_lon(i_node) = lon

    call cart2sph (dom%ccentre%elts(TRIAG*id+LORT+1), lon, lat)    
    grid_corner_lat(1,i_node) = lat
    grid_corner_lon(1,i_node) = lon

    call cart2sph (dom%ccentre%elts(TRIAG*id+UPLT+1), lon, lat)    
    grid_corner_lat(2,i_node) = lat
    grid_corner_lon(2,i_node) = lon

    call cart2sph (dom%ccentre%elts(TRIAG*idW+LORT+1), lon, lat)    
    grid_corner_lat(3,i_node) = lat
    grid_corner_lon(3,i_node) = lon

    call cart2sph (dom%ccentre%elts(TRIAG*idSW+UPLT+1), lon, lat)    
    grid_corner_lat(4,i_node) = lat
    grid_corner_lon(4,i_node) = lon

    call cart2sph (dom%ccentre%elts(TRIAG*idSW+LORT+1), lon, lat)    
    grid_corner_lat(5,i_node) = lat
    grid_corner_lon(5,i_node) = lon

    call cart2sph (dom%ccentre%elts(TRIAG*idS+UPLT+1), lon, lat)    
    grid_corner_lat(6,i_node) = lat
    grid_corner_lon(6,i_node) = lon


  end subroutine grid_coords
  
  !************************************************************************
  !handle_err
  !************************************************************************
  !
  !ROUTINE:      handle_err
  !DESCRIPTION:  error handler
  !--------------------------------------------------------------------------

  subroutine handle_err (status)
    implicit         none
    integer          status

    if (status .ne. nf_noerr) then
       print *, nf_strerror(status)
       stop 'Stopped'
    endif

  end subroutine handle_err


  !
  ! write netCDF grid descriptor file for ESMF remapping
  ! 
  subroutine wrt_esmf_rll (grid_name)
    implicit none
    character(35), intent(in) :: grid_name
    
    !-----------------------------------------------------------------------
    !
    !     grid coordinates and masks
    !
    !-----------------------------------------------------------------------
    integer  :: ncstat             ! general netCDF status variable
    integer  :: nc_grid_id         ! netCDF grid dataset id
    integer  :: nc_lon_id          ! netCDF grid dataset id
    integer  :: nc_lat_id          ! netCDF grid dataset id
    integer  :: nc_pole_id         ! netCDF grid dataset id
    integer  :: nc_gridsize_id     ! netCDF grid size dim id
    integer  :: nc_gridcorn_id     ! netCDF grid corner dim id
    integer  :: nc_gridrank_id     ! netCDF grid rank dim id
    integer  :: nc_griddims_id     ! netCDF grid dimension size id
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

    ncstat = nf_put_var_int (nc_grid_id, nc_griddims_id, grid_dims)
    call handle_err (ncstat)

    ncstat = nf_put_var_int (nc_grid_id, nc_grdimask_id, grid_imask)
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

    !
    ! Close output file
    !        
    ncstat = nf_close(nc_grid_id)
    call handle_err (ncstat)

    deallocate (grid_imask, grid_area, grid_center_lat, grid_center_lon, grid_corner_lat, grid_corner_lon)
  end subroutine wrt_esmf_rll
end module topo_grid_descriptor_mod
