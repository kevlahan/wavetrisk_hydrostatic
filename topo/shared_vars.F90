module shared_vars
  use shr_kind_mod, only: r8 => shr_kind_r8

  !---ARH
  !+++ARH
  real(r8), allocatable, dimension(:) :: landm_coslat, landfrac, terr, var30, refine_l
  !real(r8), allocatable, dimension(:) :: landm_coslat, terr, var30, refine_l
  !---ARH
  integer,  allocatable, dimension(:) :: refine_li

  !+++ARH
  real(r8), allocatable, dimension(:) :: landfrac_target, terr_target, sgh30_target, sgh_target
  !real(r8), allocatable, dimension(:) :: terr_target, sgh30_target, sgh_target
  !---ARH
  real(r8), allocatable, dimension(:) :: landm_coslat_target
  !+++ARH
  real(r8), allocatable, dimension(:) :: sumwgts_target
  !---ARH
  real(r8), allocatable, dimension(:) :: terr_uf_target, sgh_uf_target

  real(r8) , allocatable, dimension(:,:,:) :: terr_sm, terr_dev

  real(r8) :: pi, piq, pih, deg2rad, rad2deg, rotate_cube

contains 
  subroutine set_constants
    implicit none
    pi          = 4.D0*DATAN(1.D0)
    piq         = 0.25*pi
    pih         = 0.50*pi
    deg2rad     = pi/180.0
    rad2deg     = 180.0/pi
    rotate_cube = 0.0D0 !default is 0
  end subroutine set_constants

  subroutine progress_bar(txt, n, x)
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    character(*) :: txt
    integer :: n
    real(r8) :: x
    integer, parameter :: s = 5
    character :: c(0:s-1) = (/ achar(127), "/", "|", "/","-" /)
    character, parameter :: CR = achar(13)
    write( *, "((a1,a, t4,i10, f10.2,' percent  done ', a1, '  '))", advance = "NO") CR, txt, n, x, c(mod(n, s))
  end subroutine progress_bar

  subroutine smooth_terrain (lexternal_smooth_terr, ltarget_latlon, terr_target, ntarget, externally_smoothed_topo_file, &
       nlon, nlat)
    implicit none
#     include         <netcdf.inc>
    logical, intent(in) :: lexternal_smooth_terr, ltarget_latlon
    integer, intent(in) :: ntarget,nlon,nlat
    character(len=1024), intent(in) :: externally_smoothed_topo_file
    real(r8), intent(out):: terr_target(ntarget)

    integer :: ncid,status, alloc_error, ntarget_id, ntarget_smooth, phisid
    integer :: i,j,ii
    integer :: nlon_smooth,nlat_smooth
    real(r8), allocatable :: terr_smooth(:,:)

    write(*,*) "smoothing PHIS"
    if (lexternal_smooth_terr) then
       write(*,*) "Using externally generated smoothed topography"

       status = nf_open(externally_smoothed_topo_file, 0, ncid)
       if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)           
       !
       if (.NOT.ltarget_latlon) then
          !
          !*********************************************************
          !
          ! read in smoothed topography
          !
          !*********************************************************
          !
          status = NF_INQ_DIMID (ncid, 'ncol', ntarget_id    )
          status = NF_INQ_DIMLEN(ncid, ntarget_id , ntarget_smooth)
          if (ntarget/=ntarget_smooth) then
             write(*,*) "mismatch in smoothed data-set and target grid specification"
             write(*,*) ntarget, ntarget_smooth
             STOP
          end if
          status = NF_INQ_VARID(ncid, 'PHIS', phisid)
          !
          ! overwrite terr_target with smoothed version
          !
          status = NF_GET_VAR_DOUBLE(ncid, phisid,terr_target)
          terr_target = terr_target/9.80616
       else
          !
          ! read in smoothed lat-lon topography
          !
          status = NF_INQ_DIMID(ncid, 'lon', ntarget_id)
          status = NF_INQ_DIMLEN(ncid, ntarget_id, nlon_smooth)
          status = NF_INQ_DIMID(ncid, 'lat', ntarget_id)
          status = NF_INQ_DIMLEN(ncid, ntarget_id, nlat_smooth)
          if (nlon/=nlon_smooth.OR.nlat/=nlat_smooth) then
             write(*,*) "smoothed topography dimensions do not match target grid dimensions"
             write(*,*) "target grid  : nlon       ,nlat        =",nlon,nlat
             write(*,*) "smoothed topo: nlon_smooth,nlat_smooth =",nlon_smooth,nlat_smooth
             stop
          end if
          ALLOCATE (terr_smooth(nlon_smooth,nlat_smooth),stat=alloc_error)
          status = NF_INQ_VARID(ncid, 'PHIS', phisid)
          status = NF_GET_VAR_DOUBLE(ncid, phisid,terr_smooth)
          !
          ! overwrite terr_target with smoothed version
          !
          ii=1
          DO j=1,nlat
             DO i=1,nlon
                terr_target(ii) = terr_smooth(i,j)/9.80616                  
                ii=ii+1
             end DO
          end DO
          deallocate (terr_smooth)
       end if
    end if

  end subroutine smooth_terrain

  subroutine allocate_target_vars (ntarget)
    implicit none
#     include         <netcdf.inc>
    integer, intent(in) :: ntarget

    integer :: alloc_error
    allocate (terr_target(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for terr_target'
       stop
    end if
    !+++ARH
    allocate (landfrac_target(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for landfrac_target'
       stop
    end if
    !---ARH
    allocate (landm_coslat_target(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for landfrac_target'
       stop
    end if
    allocate (sgh30_target(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for sgh30_target'
       stop
    end if
    allocate (sgh_target(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for sgh_target'
       stop
    end if

    allocate (terr_uf_target(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for terr_uf_target'
       stop
    end if
    allocate (sgh_uf_target(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for sgh_uf_target'
       stop
    end if
    !+++ARH
    allocate (sumwgts_target(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for sgh_target'
       stop
    end if
    !---ARH
  end subroutine allocate_target_vars


  subroutine read_intermediate_cubed_sphere_grid (intermediate_cubed_sphere_fname, ncube, llandfrac)
    implicit none
#     include         <netcdf.inc>
    integer              :: ncid,status, dimid, alloc_error, landid,n
    integer, intent(out) :: ncube
    logical, intent(out) :: llandfrac
    character(len=1024), intent(in) :: intermediate_cubed_sphere_fname

    !
    !****************************************************
    !
    ! get dimension of cubed-sphere grid
    !
    !****************************************************
    !
    write(*,*) "Opening intermediate cubed-sphere file : ",TRIM(intermediate_cubed_sphere_fname)
    status = nf_open(TRIM(intermediate_cubed_sphere_fname), 0, ncid)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)

    status = NF_INQ_DIMID(ncid, 'grid_size', dimid)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, n)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)

    ncube = INT(SQRT(DBLE(n/6)))
    write(*,*) "cubed-sphere dimension: ncube = ",ncube
    write(*,*) "average grid-spacing at the Equator (degrees):" ,90.0/ncube

    if (status .ne. NF_NOERR) call handle_err(status)          
    !
    !****************************************************
    !
    ! read cubed-sphere 3km data
    !
    !****************************************************
    !
    status = NF_INQ_DIMID(ncid, 'grid_size', dimid)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, n)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)

    !
    !  Calculate no. of cells per side of cubed-sphere faces
    !   - 6 faces of ncube x ncube cells
    !  
    ncube = INT(SQRT(DBLE(n/6)))
    write(*,*) "cubed-sphere dimension, ncube: ",ncube


    !********************************************
    !
    !  Begin reading variables from file
    !   - All are dimension(n) where n is the
    !     number of cells in the cubed sphere 
    !
    !********************************************
    !
    !  read LANDM_COSLAT. A smoothed land fraction variable.
    !  (left in for backwards compatibility with CAM4)
    !
    allocate ( landm_coslat(n),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for landm_coslat'
       stop
    end if

    status = NF_INQ_VARID (ncid, 'LANDM_COSLAT', landid)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)

    status = NF_GET_VAR_DOUBLE(ncid, landid,landm_coslat)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)
    write(*,*) "min/max of landm_coslat",MINVAL(landm_coslat),MAXVAL(landm_coslat)

    !+++ARH
    !!
    !! read LANDFRAC
    !!
    allocate ( landfrac(n),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for landfrac'
       stop
    end if
    !
    status = NF_INQ_VARID(ncid, 'LANDFRAC', landid)
    if (status /= NF_NOERR) then
       write(*,*) "LANDFRAC not on file"
       llandfrac = .false.
       landfrac  = 1.0
    else  
       status = NF_GET_VAR_DOUBLE(ncid, landid,landfrac)
       if (status /= NF_NOERR) CALL HANDLE_ERR(status)
       write(*,*) "min/max of landfrac",MINVAL(landfrac),MAXVAL(landfrac)
       llandfrac = .true.
    end if

    !---ARH
    !
    ! read terr - this is the elevation data (meters) on cubed sphere grid
    !
    allocate ( terr(n),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for terr'
       stop
    end if

    status = NF_INQ_VARID(ncid, 'terr', landid)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)

    status = NF_GET_VAR_DOUBLE(ncid, landid,terr)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)
    write(*,*) "min/max of terr",MINVAL(terr),MAXVAL(terr)

    !
    ! read var30 (variance) of 1km topography in each
    ! cubed sphere cell (in meters) 
    !
    allocate ( var30(n),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for var30'
       stop
    end if

    status = NF_INQ_VARID(ncid, 'var30', landid)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)

    status = NF_GET_VAR_DOUBLE(ncid, landid,var30)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)
    write(*,*) "min/max of sgh30",MINVAL(SQRT(var30)),MAXVAL(SQRT(var30))
    print *,"close file"
    status = nf_close (ncid)
    if (status .ne. NF_NOERR) call handle_err(status)

    write(*,*) 'done reading in data from netCDF file'

  end subroutine read_intermediate_cubed_sphere_grid

  subroutine read_target_grid (grid_descriptor_fname, lregional_refinement, ltarget_latlon, lpole, &
       nlat, nlon, ntarget, ncorner, nrank, dom, idx, &
       target_corner_lon, target_corner_lat, target_center_lon, target_center_lat, target_area, target_rrfac)
    implicit none
#     include         <netcdf.inc>
    integer,                               intent(out) :: nlat,nlon,ntarget,ncorner,nrank
    integer,  allocatable, dimension(:),   intent(out) :: dom, idx
    real(r8), allocatable, dimension(:)  , intent(out) :: target_center_lon, target_center_lat, target_area
    real(r8), allocatable, dimension(:)  , intent(out) :: target_rrfac
    real(r8), allocatable, dimension(:,:), intent(out) :: target_corner_lon, target_corner_lat
    logical,                               intent(out) :: lregional_refinement
    logical,                               intent(out) :: lpole, ltarget_latlon
    character(len=1024),                   intent(in)  :: grid_descriptor_fname

    integer :: ncid,status
    integer :: ntarget_id, ncorner_id, nrank_id, nodeCount_id,nodeCoords_id,elementConn_id,numElementConn_id,centerCoords_id
    integer :: alloc_error
    integer :: domid, idxid, lonid, latid, nodeCount
    integer,  allocatable, dimension(:)   :: numElementConn
    integer,  allocatable, dimension(:,:) :: elementConn
    real(r8), allocatable, dimension(:,:) :: centerCoords,nodeCoords
    !+++ARH
    integer :: rrfacid
    !---ARH
    integer                         :: esmf_file = 1  ! assume SCRIP naming convention
    integer                         :: icorner, icell,num
    integer,           dimension(2) :: grid_dims
    character(len=23), dimension(2) :: str_size, str_corners, str_rank, str_corner_lon, str_corner_lat, str_area,str_rrfac
    character(len=23), dimension(2) :: str_center_lon, str_center_lat
    character(len=23), dimension(2) :: str_dom, str_id

    ltarget_latlon = .FALSE.
    !
    !*********************************************************
    !
    ! read in target grid
    !
    !*********************************************************
    !
    !
    ! SCRIP and ESMF naming convention
    ! 
    str_size         = (/"grid_size             ","elementCount          "/)
    str_corners      = (/"grid_corners          ","maxNodePElement       "/)
    str_rank         = (/"grid_rank             ","1                     "/)
    str_corner_lon   = (/"grid_corner_lon       ","nodeCoords            "/) !(nodeCount, coordDim)
    str_corner_lat   = (/"grid_corner_lat       ","nodeCoords            "/) !(nodeCount, coordDim)
    str_center_lon   = (/"grid_center_lon       ","centerCoordsnodeCoords"/) !(elementCount, coordDim)
    str_center_lat   = (/"grid_center_lat       ","centerCoordsnodeCoords"/) !(elementCount, coordDim)
    str_area         = (/"grid_area             ","elementArea           "/)
    str_rrfac        = (/"rrfac                 ","elementRefinementRatio"/)
    str_dom          = (/"grid_dom              ","nodeDomain            "/)
    str_id           = (/"grid_id               ","nodeID                "/)

    write(*,*) "Opening grid descriptor file :  ",TRIM(grid_descriptor_fname)
    status = nf_open (TRIM(grid_descriptor_fname), 0, ncid)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)

    status = NF_INQ_DIMID (ncid, TRIM(str_size(1)), ntarget_id)
    status = NF_INQ_DIMLEN (ncid, ntarget_id, ntarget)
    write(*,*) "dimension of target grid: ntarget = ", ntarget

    status = NF_INQ_DIMID (ncid, TRIM(str_corners(esmf_file)), ncorner_id)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)

    status = NF_INQ_DIMLEN (ncid, ncorner_id, ncorner)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
    write(*,*) "maximum number of corners: ncorner = ",ncorner

    status = NF_INQ_DIMID (ncid, TRIM(str_rank(esmf_file)), nrank_id)
    status = NF_INQ_DIMLEN (ncid, nrank_id, nrank)

    allocate (target_center_lon(ntarget), stat=alloc_error)
    allocate (target_center_lat(ntarget), stat=alloc_error)
    allocate (target_area      (ntarget), stat=alloc_error)

    allocate (target_corner_lon(ncorner,ntarget), stat=alloc_error)
    allocate (target_corner_lat(ncorner,ntarget), stat=alloc_error)

    allocate (dom(ntarget), stat=alloc_error)
    allocate (idx(ntarget), stat=alloc_error)

    write(*,*) "grid rank: nrank = ", nrank

    if (nrank==2) then
       write(*,*) "target grid is a lat-lon grid"
       ltarget_latlon = .TRUE.
       status = NF_INQ_VARID(ncid,'grid_dims', ntarget_id)
       if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
       status = NF_GET_VAR_INT(ncid, ntarget_id,grid_dims)
       nlon = grid_dims(1)
       nlat = grid_dims(2)

       write(*,*) "nlon=",nlon,"nlat=",nlat
    else if (nrank==1) then
       ltarget_latlon = .FALSE.
    else
       write(*,*) "nrank out of range", nrank
       stop
    endif

    status = NF_INQ_VARID (ncid, TRIM(str_dom(esmf_file)), domid)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
    status = NF_GET_VAR_INT (ncid, domid, dom)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)

    status = NF_INQ_VARID (ncid, TRIM(str_id(esmf_file)), idxid)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
    status = NF_GET_VAR_INT (ncid, idxid, idx)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)

    status = NF_INQ_VARID (ncid, TRIM(str_corner_lon(esmf_file)), lonid)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
    status = NF_GET_VAR_DOUBLE (ncid, lonid, target_corner_lon)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
    
    if (maxval(target_corner_lon)>10.0) target_corner_lon = deg2rad*target_corner_lon

    status = NF_INQ_VARID(ncid, TRIM(str_corner_lat(esmf_file)), latid)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
    status = NF_GET_VAR_DOUBLE(ncid, latid,target_corner_lat)

    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
    if (maxval (target_corner_lat) > 10.0) target_corner_lat = deg2rad*target_corner_lat    
    !
    ! for writing remapped data on file at the end of the program
    !    
    status = NF_INQ_VARID(ncid, 'grid_center_lon', lonid)
    status = NF_GET_VAR_DOUBLE(ncid, lonid, target_center_lon)
    if (maxval(target_center_lon)>10.0) target_center_lon = deg2rad*target_center_lon    

    status = NF_INQ_VARID(ncid, 'grid_center_lat', latid)
    status = NF_GET_VAR_DOUBLE(ncid, latid,target_center_lat)
    if (maxval(target_center_lat)>10.0) target_center_lat = deg2rad*target_center_lat    
    if (ltarget_latlon) then
       if (maxval(target_center_lat)>pih-1E-5) then
          lpole=.true.
       end if
    end if
    if (lpole) then      
       write(*,*) "center of most Northern grid cell is lat=90; similarly for South pole"
    else
       write(*,*) "center of most Northern grid cell is NOT lat=90; similarly for South pole"
    end if

    status = NF_INQ_VARID (ncid, TRIM(str_rrfac(esmf_file)), rrfacid)
    if (STATUS /= NF_NOERR) then
       lregional_refinement = .false.
       write(*,*) "rrfac not on file; setting lregional_refinement = .false."
       write(*,*) "  ... allocating target_rrfac ANYWAY"
       ! allocate target_rrfac anyway since it may be invoked
       ! if rrfac write is requested
       allocate ( target_rrfac(ntarget),stat=alloc_error)
    else
       lregional_refinement = .true.
       allocate ( target_rrfac(ntarget),stat=alloc_error)
       status = NF_GET_VAR_DOUBLE(ncid, rrfacid,target_rrfac)
       if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)
       write(*,*) "rrfac on file; setting lregional_refinement = .true."
    end if

    status = NF_INQ_VARID (ncid, TRIM(str_area(esmf_file)), latid)
    status = NF_GET_VAR_DOUBLE (ncid, latid, target_area)

    status = nf_close (ncid)
    if (status .ne. NF_NOERR) call handle_err(status)          
  end subroutine read_target_grid

  subroutine read_refinement_factor(rrfactor_fname,ncube_rr)
    implicit none
#     include         <netcdf.inc>
    character(len=1024), intent(in) :: rrfactor_fname
    integer, intent(out) :: ncube_rr
    integer :: ncid,status, dimid, alloc_error, landid,n

    write(*,*) "Opening file w/ refinement factors : ",TRIM(rrfactor_fname)
    status = nf_open(TRIM(rrfactor_fname) , 0, ncid)
    if (STATUS /= NF_NOERR) CALL HANDLE_ERR(STATUS)

    status = NF_INQ_DIMID(ncid, 'grid_size', dimid)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, n)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)

    ncube_rr = INT(SQRT(DBLE(n/6)))
    write(*,*) "RR_factor dimension: ncube = ",ncube_rr

    !
    ! read  floating point refinement level
    !
    allocate ( refine_l(n),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for refine_l'
       stop
    end if

    status = NF_INQ_VARID(ncid, 'refine_level', landid)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)

    status = NF_GET_VAR_DOUBLE(ncid, landid,refine_l)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)
    write(*,*) "min/max of FLOAT Refine level",MINVAL(refine_l),MAXVAL(refine_l)

    !
    ! read  INTEGER refinement level
    !
    allocate ( refine_li(n),stat=alloc_error )
    if( alloc_error /= 0 ) then
       print*,'Program could not allocate space for refine_l'
       stop
    end if

    status = NF_INQ_VARID(ncid, 'refine_level', landid)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)

    status = NF_GET_VAR_INT(ncid, landid,refine_li)
    if (status /= NF_NOERR) CALL HANDLE_ERR(status)
    write(*,*) "min/max of INT Refine level",MINVAL(refine_li),MAXVAL(refine_li)


    print *,"close file"
    status = nf_close (ncid)
    if (status .ne. NF_NOERR) call handle_err(status)

    write(*,*) 'done reading in data from netCDF file'
  end subroutine read_refinement_factor
end module shared_vars
