!
!  dESCRIPTION:  Remap topo data from cubed-sphere grid to target grid using rigorous remapping
!                (Lauritzen, Nair and Ullrich, 2010, J. Comput. Phys.)
!
!  author: Peter Hjort Lauritzen (pel@ucar.edu), aMP/CGd/NCaR 
!          Julio Bacmeister, aMP/CGd/NCaR 
!          adam Herrington, aMP/CGd/NCaR
!
! ex: ./cube_to_target --help to get list of long and short option names.
!
program convterr
  use shr_kind_mod, only: r8 => shr_kind_r8
  use smooth_topo_cube_sph
  use ridge_ana
  use shared_vars
  use reconstruct
  use f90getopt

  implicit none
#     include         <netcdf.inc>

  integer :: ncube                   ! dimension of intermediate cubed-sphere grid

  integer :: alloc_error 
  !
  ! turn extra debugging on/off
  ! 
  logical :: ldbg=.false.
  real(r8):: wt
  integer :: ii,ip,jx,jy,jp,np,counti !counters,dimensions
  integer :: jmax_segments,jall,jall_anticipated !overlap segments
  integer, parameter :: ngauss = 3               !quadrature for line integrals

  integer                               :: ntarget, ncorner, nrank, nlon, nlat               ! target grid dimensions
  logical                               :: ltarget_latlon,lpole                              ! if target grid lat-lon
  integer , allocatable, dimension(:)   :: dom, idx
  real(r8), allocatable, dimension(:)   :: rrfac_target,target_rrfac
  real(r8), allocatable, dimension(:,:) :: target_corner_lon, target_corner_lat              ! target grid coordinates
  real(r8), allocatable, dimension(:)   :: target_center_lon, target_center_lat, target_area, area_target !target grid coordinates

  real(r8), allocatable, dimension(:,:) :: weights_all                        !overlap weights
  integer , allocatable, dimension(:)   :: weights_lgr_index_all              !overlap index
  integer , allocatable, dimension(:,:) :: weights_eul_index_all              !overlap index
  integer :: ix,iy  , i
  !
  ! volume of topography
  !
  real(r8) :: vol_target, vol_target_un, area_target_total,vol_source,area_source,mea_source

  logical :: lphis_gll=.false.
  logical :: llandfrac=.false. !if landfrac is on the intermediate cubed-sphere file it will be mapped to target grid
  logical :: lzero_negative_peaks  = .TRUE.
  !
  ! for internal filtering
  !
  real(r8), allocatable, dimension(:,:) :: da, terr_2(:,:,:),rrfac(:,:,:)
  integer  :: nreconstruction
  real(r8) :: da_min_ncube, da_min_target ! used to compute jmax_segments
  real(r8) :: volterr, volterr_sm  
  !
  ! namelist variables
  !
  logical :: ldevelopment_diags    = .false.
  logical :: lread_smooth_topofile = .false.
  logical :: luse_prefilter        = .false.
  logical :: lstop_after_smoothing = .false.
  logical :: lrrfac_manipulation   = .false.
  !
  ! Cubed sphere terr is band-pass filtered using circular kernels
  !                             *Radii* of smoothing circles
  integer :: ncube_sph_smooth_coarse = -1
  integer :: ncube_sph_smooth_fine   =  0
  !
  ! namelist variables for detection of sub-grid scale orientation
  ! i.e., "ridge finding"
  !
  logical :: lfind_ridges = .TRUE.
  !                             Ridge analysis takes place on
  !                             squares of 2*NW+1
  integer :: nwindow_halfwidth =  0
  !                             
  !                             for backwards compat with CESM2.0
  !                             Not used, 0 here for naming
  integer :: nridge_subsample = 0 !
  !
  logical :: lridgetiles = .false.

  logical :: lregional_refinement = .false. !set in read_target_grid if rrfac is on file
  integer :: rrfac_max = 1
  logical :: lread_pre_smoothtopo = .false.      !use pre-smoothed (on intermediate cubed-sphere grid) topo file
  logical :: lwrite_rrfac_to_topo_file = .false. !for debugging write rrfac on target grid to topo file
  logical :: linterp_phis = .false.              !interpolate PHIS to grid center (instead of area average remapping; used for GLL grids)
  logical :: lsmoothing_over_ocean = .false.      !default is that no smoothing is applied where landfrac=0; turn off
  logical :: ldistance_weighted_smoother = .false.!use distance weighted smoother instead of Laplacian smoother

  real (r8):: nu_lap = -1
  integer  :: smooth_phis_numcycle=-1
  real (r8):: smoothing_scale=0
  !
  integer :: UNIT, ioptarg

  integer :: NSCL_f, NSCL_c, nhalo,nsw

  !++JTB
  integer :: iopt_ridge_seed = 2
  ! cube quantities for remapping
  real(r8), allocatable, dimension(:) :: uniqiC , uniqwC,  cwghtC, wedgoC
  real(r8), allocatable, dimension(:) :: anglxC,  anisoC,  hwdthC, clngtC, mxdisC
  real(r8), allocatable, dimension(:) :: riseqC,  fallqC,  mxvrxC, mxvryC, nodesC
  integer,  allocatable, dimension(:) :: itrgtC
  !--JTB

  !
  ! namelist filenames
  !
  character(len=1024) :: grid_descriptor_fname,intermediate_cubed_sphere_fname,output_fname=''
  character(len=1024) :: grid_descriptor_fname_gll
  character(len=1024) :: output_grid='', ofile,smooth_topo_fname = '',str_dir=''
  character(len=1024) :: rrfactor_fname, command_line_arguments, str, str_creator, str_source=''

  character(len=8)  :: date
  character(len=10) :: time


  type(option_s):: opts(23)
  !               
  !                     long name                   has     | short | specified    | required
  !                                                 argument| name  | command line | argument
  ! 
  opts(1 ) = option_s( "smoothing_scale"           ,.true.    , 'c'   ,.false.       ,.true.)
  opts(2 ) = option_s( "fine_radius"               ,.true.    , 'f'   ,.false.       ,.false.)!xxx remove
  opts(3 ) = option_s( "grid_descriptor_file"      ,.true.    , 'g'   ,.false.       ,.true.)
  opts(4 ) = option_s( "help"                      ,.false.   , 'h'   ,.false.       ,.false.)
  opts(5 ) = option_s( "intermediate_cs_name"      ,.true.    , 'i'   ,.false.       ,.true.)
  opts(6 ) = option_s( "output_grid"               ,.true.    , 'o'   ,.false.       ,.true.)
  opts(7 ) = option_s( "use_prefilter"             ,.false.   , 'p'   ,.false.       ,.false.)
  opts(8 ) = option_s( "no_ridges"                 ,.false.   , 'r'   ,.false.       ,.false.)
  opts(9)  = option_s( "stop_after_smooth"         ,.false.   , 'x'   ,.false.       ,.false.)
  opts(10) = option_s( "rrfac_max"                 ,.true.    , 'y'   ,.false.       ,.false.)
  opts(11) = option_s( "rrfac_manipulation"        ,.false.   , 'v'   ,.false.       ,.false.)
  opts(12) = option_s( "development_diags"         ,.false.   , 'z'   ,.false.       ,.false.)
  opts(13) = option_s( "ridge2tiles"               ,.false.   , '1'   ,.false.       ,.false.)
  opts(14) = option_s( "smooth_topo_file"          ,.true.    , 't'   ,.false.       ,.false.)
  opts(15) = option_s( "write_rrfac_to_topo_file"  ,.false.   , 'd'   ,.false.       ,.false.)
  opts(16) = option_s( "name_email_of_creator"     ,.true.    , 'u'   ,.false.       ,.true.)
  opts(17) = option_s( "source_data_identifier"    ,.true.    , 'n'   ,.false.       ,.false.)
  opts(18) = option_s( "output_data_directory"     ,.true.    , 'q'   ,.false.       ,.false.)
  opts(19) = option_s( "grid_descriptor_file_gll"  ,.true.    , 'a'   ,.false.       ,.false.)
  opts(20) = option_s( "interpolate_phis"          ,.false.   , 's'   ,.false.       ,.false.)
  opts(21) = option_s( "distance_weighted_smoother",.false.   , 'b'   ,.false.       ,.false.)
  opts(22) = option_s( "smooth_phis_numcycle"      ,.true.    , 'l'   ,.false.       ,.false.)
  opts(23) = option_s( "smoothing_over_ocean"      ,.false.   , 'm'   ,.false.       ,.false.)

  ! end longopts
  ! If no options were committed
  if (command_argument_count() .eq. 0 ) call print_help
  !
  ! collect command line arguments in this string for netCdF meta data
  !
  command_line_arguments    = './cube_to_target'
  grid_descriptor_fname     = ''
  grid_descriptor_fname_gll = ''

  ! Process options one by one
  do
     select case( getopt( "c:f:g:hi:o:prxy:vz1:t:du:n:q:a:sbl:m", opts ) ) ! opts is optional (for longopts only)
     case( char(0) )
        exit
     case( 'c' )
        read (optarg, *) smoothing_scale
        write(str,*) smoothing_scale
        command_line_arguments = TRIM(command_line_arguments)//' --smoothing_scale '//TRIM(adJUSTL(str))
        opts(1)%specified = .true.
     case( 'f' )
        read (optarg, '(i3)') ioptarg
        ncube_sph_smooth_fine = ioptarg
        write(str,*) ioptarg
        command_line_arguments = TRIM(command_line_arguments)//' --fine_radius '//TRIM(adJUSTL(str))
        opts(2)%specified = .true.
     case( 'g' )
        grid_descriptor_fname = optarg
        write(str,*) TRIM(optarg)
        command_line_arguments = TRIM(command_line_arguments)//' --grid_descriptor_file '//TRIM(adJUSTL(str))
        opts(3)%specified = .true.
     case( 'h' )
        call print_help
        opts(4)%specified = .true.
     case( 'i' )
        intermediate_cubed_sphere_fname = optarg
        write(str,*) TRIM(optarg)
        command_line_arguments = TRIM(command_line_arguments)//' --intermediate_cs_name '//TRIM(adJUSTL(str))
        opts(5)%specified = .true.
     case( 'o' )
        output_grid = optarg
        write(str,*) TRIM(optarg)
        command_line_arguments = TRIM(command_line_arguments)//' --output_grid '//TRIM(adJUSTL(str))
        opts(6)%specified = .true.
     case( 'p' )
        luse_prefilter=.TRUE.
        command_line_arguments = TRIM(command_line_arguments)//' --use_prefilter '
        opts(7)%specified = .true.
     case( 'r' )
        lfind_ridges = .false.
        command_line_arguments = TRIM(command_line_arguments)//' --no_ridges '
        opts(8)%specified = .true.
     case( 'x' )
        lstop_after_smoothing = .TRUE.
        command_line_arguments = TRIM(command_line_arguments)//' --stop_after_smooth '//TRIM(adJUSTL(str))
        opts(9)%specified = .true.
     case( 'y' )
        read (optarg, '(i3)') ioptarg
        rrfac_max = ioptarg
        lregional_refinement =.true.
        write(str,*) ioptarg
        command_line_arguments = TRIM(command_line_arguments)//' --rrfac_max '//TRIM(adJUSTL(str))
        opts(10)%specified = .true.
     case( 'v' )
        lrrfac_manipulation= .TRUE.
        command_line_arguments = TRIM(command_line_arguments)//' --rrfac_manipulation '
        opts(11)%specified = .true.
     case( 'z' )
        ldevelopment_diags = .TRUE.
        command_line_arguments = TRIM(command_line_arguments)//' --development_diags '
        opts(12)%specified = .true.
     case( '1' )
        lridgetiles = .TRUE.
        command_line_arguments = TRIM(command_line_arguments)//' --ridge2tiles '
        opts(13)%specified = .true.
     case( 't' )
        smooth_topo_fname = optarg
        write(str,*) TRIM(optarg)
        write (6,*) str
        command_line_arguments = TRIM(command_line_arguments)//' --smooth_topo_file '//TRIM(adJUSTL(str))
        opts(14)%specified = .true.
     case( 'd' )
        lwrite_rrfac_to_topo_file = .TRUE.
        command_line_arguments = TRIM(command_line_arguments)//' --write_rrfac_to_topo_file '
        opts(15)%specified = .true.
     case( 'u' )
        str_creator = optarg
        write(str,*) TRIM(optarg)
        command_line_arguments = TRIM(command_line_arguments)//' --name_email_of_creator '//TRIM(adJUSTL(str))
        opts(16)%specified = .true.
     case( 'n' )
        str_source = optarg
        write(str,*) TRIM(optarg)
        command_line_arguments = TRIM(command_line_arguments)//' --source_data_identifier '//TRIM(adJUSTL(str))
        opts(17)%specified = .true.
     case( 'q' )
        str_dir = optarg
        write(str,*) TRIM(optarg)
        command_line_arguments = TRIM(command_line_arguments)//' --output_data_directory '//TRIM(adJUSTL(str))
        opts(18)%specified = .true.
     case( 'a' )
        lphis_gll=.TRUE.
        grid_descriptor_fname_gll = optarg
        write(str,*) TRIM(optarg)
        command_line_arguments = TRIM(command_line_arguments)//' --grid_descriptor_file_gll '//TRIM(adJUSTL(str))
        opts(19)%specified = .true.
     case( 's' )
        linterp_phis = .TRUE.
        command_line_arguments = TRIM(command_line_arguments)//' --interpolate_phis '
        opts(20)%specified = .true.
     case( 'b' )
        ldistance_weighted_smoother = .true.
        command_line_arguments = TRIM(command_line_arguments)//' --distance_weighted_smoother '//TRIM(adJUSTL(str))
        opts(21)%specified = .true.
     case( 'l' )
        read (optarg, '(i5)') smooth_phis_numcycle
        write(str,*) smooth_phis_numcycle
        write (6,*) trim(str)
        command_line_arguments = TRIM(command_line_arguments)//' --smooth_phis_numcycle '//TRIM(adJUSTL(str))
        opts(22)%specified = .true.
     case( 'm' )
        lsmoothing_over_ocean = .TRUE.
        command_line_arguments = TRIM(command_line_arguments)//' --smoothing_over_ocean '
        opts(23)%specified = .true.
     case default
        write (6,*) "Option unknown: ",char(0)        
        stop
     end select
  end do

  if (TRIM(smooth_topo_fname)/='') then
     lread_smooth_topofile = .TRUE.
     write (6,*) " Use pre-computed smooth topo " 
     write (6,*) " File = ", trim(smooth_topo_fname)
  else 
     write (6,*) " No smoothed topo file specified"
  end if
  !
  ! check that all required arguments are specified/initialized
  ! (stopping after smoothing is only for developers so it
  ! is not checked if required arguments are present)
  !
  if (.not.lstop_after_smoothing) then
     do i=1,SIZE(opts)
        if (.not.opts(i)%specified.and.opts(i)%required) then
           write (6,*) "Required argument not specified: ",opts(i)%name
           stop
        end if
     end do
  end if

  if (LEN(TRIM(str_source))==0) then
     !
     ! default setting for source topography
     !
     str_source = 'gmted2010_bedmachine-ncube0540-220518'
  end if
  if (LEN(TRIM(str_dir))==0) then
     !
     ! default output directory
     !
     str_dir = '.'
  end if

  write (6,*) " "
  write (6,*) "Namelist settings"
  write (6,*) "================="
  write (6,*)
  write (6,*) "smoothing_scale                 = ",smoothing_scale
  write (6,*) "nwindow_halfwidth               = ",nwindow_halfwidth
  write (6,*) "ncube_sph_smooth_fine           = ",ncube_sph_smooth_fine
  write (6,*) "grid_descriptor_fname           = ",trim(grid_descriptor_fname)
  write (6,*) "intermediate_cubed_sphere_fname = ",trim(intermediate_cubed_sphere_fname)
  write (6,*) "output_grid                     = ",trim(output_grid)
  write (6,*) "luse_prefilter                  = ",luse_prefilter
  write (6,*) "lfind_ridges                    = ",lfind_ridges
  write (6,*) "rrfac_max                       = ",rrfac_max
  write (6,*) "ldevelopment_diags              = ",ldevelopment_diags
  write (6,*) "lridgetiles                     = ",lridgetiles
  write (6,*) "smooth_topo_fname               = ",trim(smooth_topo_fname)
  write (6,*) "lwrite_rrfac_to_topo_file       = ",lwrite_rrfac_to_topo_file
  write (6,*) "str_source                      = ",trim(str_source)
  write (6,*) "interpolate_phis                = ",linterp_phis
  write (6,*) "ldistance_weighted_smoother     = ",ldistance_weighted_smoother
  write (6,*) "smooth_phis_numcycle            = ",smooth_phis_numcycle
  write (6,*) "smoothing_over_ocean            = ",lsmoothing_over_ocean

  !*********************************************************

  call  set_constants

  ! Read in target grid
  !------------------------------------------------------------------------------------------------
  if (.not. lstop_after_smoothing) then
     
     call read_target_grid (grid_descriptor_fname, lregional_refinement, ltarget_latlon, lpole, &
          nlat, nlon, ntarget, ncorner, nrank, dom, idx, &
          target_corner_lon, target_corner_lat, target_center_lon, target_center_lat, target_area, target_rrfac)

     if (.not. lregional_refinement .and. rrfac_max /= 1) then
        write (6,*) "User has set rrfac_max =",rrfac_max
        write (6,*) "which turns on regional refinement, however, the refinementfactor is not on grid descriptor file"
        write (6,*) "SCRIP format: rrfac; ESMF format: elementRefinementRatio"
        stop
     end if

     allocate (area_target(ntarget), stat=alloc_error )
     area_target = 0d0
  end if

  if (lregional_refinement) then
     if (lstop_after_smoothing) then
        write (6,*) "stop after smoothing is not supported for variable resolution grids!"
        write (6,*) " "
        write (6,*) "stop after smoothing is intended for efficiency when produing several"
        write (6,*) "topo files which need the same amount of smoothing but different"
        write (6,*) "target grids"
        stop
     end if
     write (6,*) "rrfac_max = ", rrfac_max
     if (rrfac_max.le.1) then
        if (rrfac_max<1) then
           write (6,*) "refinement factor must be >1"
           stop
        end if
        if (rrfac_max==1) then
           write (6,*) "max refinement factor must be specified in namelist for regional refinement"
           write (6,*) " "
           write (6,*) "  --rrfac_max xx"
           write (6,*) " "
           write (6,*) "where xx = coarse resolution / finest resolution"
           stop
        end if
     end if
     if (lrrfac_manipulation) then
        write (6,*) " "
        write (6,*) " lrrfac_manipulation=.TRUE. -> the following manipulation of rrfac_max is taking place:"
        write (6,*) " "
        write (6,*) "   rrfac = real(NinT(rrfac))"
        write (6,*) "   where (rrfac.gt.rrfac_max) rrfac = rrfac_max"
        write (6,*) "   where (rrfac.lt.1.0) rrfac = 1.0"
        write (6,*) "   Laplacian smoother is applied to rrfac"
        write (6,*) "   (same level of smoothing as PHIS)"
        write (6,*) " "
     end if
  else
     if (lrrfac_manipulation) then
        write (6,*) "setting lrrfac_manipulation=.TRUE. (namelist option -v or --rrfac_manipulation)"
        write (6,*) "has no effect when not running regional refinement. Regional refinement is "
        write (6,*) "activated with -y=R --rrfac_max=R where R is the maximum refinement factor"
        write (6,*) " "
        write (6,*) "To keep the user safe - aBORT"
        stop
     end if
  end if

  ! Read in topo data on cubed sphere grid
  !------------------------------------------------------------------------------------------------

  call read_intermediate_cubed_sphere_grid (intermediate_cubed_sphere_fname, ncube, llandfrac)

  !
  ! set derived variables - scaling for smoothing
  !  
  ncube_sph_smooth_coarse = NinT(60.0*(smoothing_scale/100.0)/(3000.0/real(ncube)))
  nu_lap                  = 20.0E7*(smoothing_scale/100.0)**2
  write (6,*) "ncube_sph_smooth_coarse=",ncube_sph_smooth_coarse
  write (6,*) "nu_lap                  =",nu_lap

  if (.not.ldistance_weighted_smoother) then
     if (smooth_phis_numcycle<0) then
        write (6,*) "Recommended setting for stability"
        smooth_phis_numcycle= (nu_lap/20.0E7)*60*(real(ncube)/540.0)**2
        write (6,*) "smooth_phis_numcycle = ",smooth_phis_numcycle
     end if
  end if
  !
  ! calculate some defaults
  !
  if (lfind_ridges) then
     if (nwindow_halfwidth<=0) then
        nwindow_halfwidth = floor(real(ncube_sph_smooth_coarse)/sqrt(2.))
        !
        ! nwindow_halfwidth does NOT actually have to be even (JTB Mar 2022)
        !
        if (nwindow_halfwidth<5) then
           write (6,*) "nwindow_halfwidth can not be < 4"
           write (6,*) "setting nwindow_halfwidth=4"
           nwindow_halfwidth = 4
        end if
     end if
     if (ncube_sph_smooth_coarse<5) then
        write (6,*) "can not find ridges when ncube_sph_smooth_coarse<5"
        stop
     end if
  end if

  if (ncube_sph_smooth_fine > 0) then 
     luse_prefilter=.TRUE.
  else
     luse_prefilter=.false.
  end if

  !
  ! sanity check
  !
  if (.not.lsmoothing_over_ocean.and..not.llandfrac) then
     write (6,*) "landfrac is needed for not smoothing over ocean"
     write (6,*) "LaNdFRaC not found in file: ",intermediate_cubed_sphere_fname
     write (6,*) "aBORT"
     stop
  end if

  allocate ( da(ncube,ncube),stat=alloc_error )
  call Equiangularallareas (ncube, da)

  !*********************************************************
  !
  ! set standard output file name
  !
  !*********************************************************
  if (ldistance_weighted_smoother) then
     if ( ncube_sph_smooth_fine==0) then
        if (lfind_ridges) then
           nsw = nwindow_halfwidth
           write (ofile , "(i0.4, '_Co',i0.3)" ) ncube_sph_smooth_coarse
        else
           write (ofile, "(i0.3)" ) ncube_sph_smooth_coarse
        endif
     else
        if (lfind_ridges) then
           nsw = nwindow_halfwidth
           write (ofile ,"(i0.4,'_Co',i0.3,'_Fi',i0.3 )" ) ncube_sph_smooth_coarse, ncube_sph_smooth_fine
        else
           write( ofile , "(i0.3,'_Fi',i0.3)" ) ncube_sph_smooth_coarse, ncube_sph_smooth_fine
        endif
     end if
  else
     !
     ! Laplacian smoother standard file name
     !
     if (lfind_ridges) then
        nsw = nwindow_halfwidth        
        if (lsmoothing_over_ocean) then
           write (ofile, "('_',i0.4,'km')") nint (smoothing_scale)
        else
           write (ofile ,"('_',i0.4,'km','_noleak')" )  nint (smoothing_scale)
        endif
     else
        if (lsmoothing_over_ocean) then
           write (ofile, "('_',i0.4,'km')" ) NinT (smoothing_scale)
        else
           write (ofile , "(i0.4,'km','_noleak')" ) nint (smoothing_scale)
        end if
     endif
  end if

  output_fname = TRIM(str_dir)//'/'//trim(output_grid)//'_'//trim(str_source)//trim(ofile)//'.nc'
  write (6,*) "Writing topo file to ", output_fname
  
  !+++aRH
  ! Compute overlap weights
  !------------------------------------------------------------------------------------------------

  ! On entry to overlap_weights 'jall' is a generous guess at the number of cells in
  ! in the 'exchange grid'
  allocate ( rrfac(ncube,ncube,6)  )
  rrfac = 0d0
  if (.not. lstop_after_smoothing) then    
     if (nrank == 1) then
        da_min_ncube  = 4.0*pi/(6.0*dBLE(ncube*ncube))
        da_min_target = maxval(target_area)
        if (da_min_target==0) then !bug with MPaS files
           write (6,*) "ERROR: da_min_target =",da_min_target
           stop
        else
           write (6,*) "using dynamic estimate for jmax_segments " 
           jmax_segments = 4 * ncorner * nint (da_min_target/da_min_ncube)
        end if
        write (6,*) "ncorner, da_min_target, da_min_ncube =", ncorner, da_min_target, da_min_ncube
        write (6,*) "jmax_segments",jmax_segments,da_min_target,da_min_ncube
     else
        jmax_segments = 100000   !can be tweaked
     end if
     if (real(ntarget)*real(jmax_segments)>huge(real(jall_anticipated))) then
        jall_anticipated = 1080000000 !huge(jmax_segments) !anticipated number of weights (can be tweaked)
        write (6,*) "truncating jall_anticipated to ",jall_anticipated
     else
        jall_anticipated = ntarget*jmax_segments !anticipated number of weights (can be tweaked)
     end if
     if (jall_anticipated<0) then
        write (6,*) "anticipated number of overlaps likely not representable: jall_anticipated=", jall_anticipated
        jall_anticipated = 1080000000
        write (6,*) "setting to large value = ",jall_anticipated
     else
        write (6,*) "anticipated number of overlaps jall_anticipated=", jall_anticipated
     end if

     jmax_segments = Min( jmax_segments, 10000 )

     nreconstruction = 1
     allocate (weights_all(jall_anticipated,nreconstruction),stat=alloc_error )
     allocate (weights_eul_index_all(jall_anticipated,3),stat=alloc_error )
     allocate (weights_lgr_index_all(jall_anticipated),stat=alloc_error )
     jall=jall_anticipated

     if (.not.lstop_after_smoothing) then
        write (6,*) "Compute overlap weights: "
        call overlap_weights(weights_lgr_index_all,weights_eul_index_all,weights_all,&
             jall,ncube,ngauss,ntarget,ncorner,jmax_segments,target_corner_lon,target_corner_lat,nreconstruction,ldbg)
     end if
     deallocate(target_corner_lon,target_corner_lat)
  end if
  ! On exit from overlap_weights 'jall' is the correct number of cells in the exchange grid. 
  !------------------------------------------------------------------------------------------------

  ! Set-up regional refinement control.
  !------------------------------------------
  ! array rrfac is a refinement factor >= 1.0 
  ! Passed to smooth topo and ridge finder to
  ! control lengthscales used in algorithms. 
  ! RRfac is always used. If output_grid has no 
  ! regional refinement then rrfac(:,:,:)=1.

  if (lregional_refinement) then
     !--- remap rrfac to cube
     !-----------------------------------------------------------------
     do counti=1,jall
        i    = weights_lgr_index_all(counti)!!
        !
        ix  = weights_eul_index_all(counti,1)
        iy  = weights_eul_index_all(counti,2)
        ip  = weights_eul_index_all(counti,3)
        !
        ! convert to 1d indexing of cubed-sphere
        !
        ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
        !
        wt = weights_all(counti,1)
        !
        rrfac(ix,iy,ip) = rrfac(ix,iy,ip) + wt*(target_rrfac(i))/da(ix,iy)
     end do
  else
     write (6,*) " NO refinement: RRFaC = 1. everywhere "
     rrfac = 1d0
  endif
  write (6,*) "MinMaX RRFaC RaW MaPPEd FIELd",minval(rrfac),maxval(rrfac)

  !++jtb
  NSCL_c = 2*ncube_sph_smooth_coarse
  NSCL_f = 2*ncube_sph_smooth_fine
  nhalo  = NSCL_c  !*ncube_sph_smooth_iter ! 120      

  allocate( terr_sm(ncube,ncube,6)  )
  allocate( terr_dev(ncube,ncube,6) )
  allocate( terr_2(ncube,ncube,6)  )
  terr_2 = reshape( terr,    (/ncube,ncube,6/) )

  write (6,*) " SMOOTHING ON CUBED SPHERE "

  if (NSCL_c > 0 .or. .not.ldistance_weighted_smoother) then
     !+++aRH
     !!NSCL_c = 4*2*ncube_sph_smooth_coarse
     nhalo  = NSCL_c


     !---rrfac limiting
     if (rrfac_max>1.and.lrrfac_manipulation) then
        rrfac = real(NinT(rrfac))
        where (rrfac.gt.rrfac_max) rrfac = rrfac_max
        where (rrfac.lt.1.0) rrfac = 1.0
        write (6,*) "RRFaC Massaged .... "
        write (6,*) "MinMaX RRFaC FinaL",minval(rrfac),maxval(rrfac)
     end if


     write (6,*) "Entering smooth_intermediate_topo_wrap ..."

     call  smooth_intermediate_topo_wrap (terr, rrfac,da,  & 
          ncube,nhalo, NSCL_f,NSCL_c, &
          terr_sm, terr_dev ,         &
          smooth_topo_fname,          &
          lread_smooth_topofile,      &
          ldistance_weighted_smoother,&
          luse_prefilter, &
          lstop_after_smoothing, &
          lregional_refinement, rrfac_max, &
          ldevelopment_diags, &
          command_line_arguments,str_dir,str_source,&
          output_grid,&
          nu_lap, smooth_phis_numcycle,landfrac,&
          lsmoothing_over_ocean,lrrfac_manipulation,&
          smooth_topo_fname=smooth_topo_fname&
          )
  else
     terr_dev = terr_2
  endif

  volterr=0.
  volterr_sm=0.
  do np=1,6 
     volterr    =  volterr    + sum( terr_2(:,:,np) * da )
     volterr_sm =  volterr_sm + sum( terr_sm(:,:,np) * da )
  end do

  write (6,*) " Topo volume BEFORE smoother = ",volterr/(6*sum(da))
  write (6,*) " Topo volume  aFTER smoother = ",volterr_sm/(6*sum(da))
  write (6,*) "            difference       = ",(volterr - volterr_sm)/(6*sum(da))

  if (ldistance_weighted_smoother) then
     terr_sm = (volterr/volterr_sm)*terr_sm! should we do this?
  end if
  volterr_sm = 0d0
  do np=1,6 
     volterr_sm =  volterr_sm + sum( terr_sm(:,:,np) * da )
  end do

  write (6,*) " Topo volume  aFTER smoother aNd fixer = ",volterr_sm/(6*sum(da))


  if (lfind_ridges) then
     nsw = nwindow_halfwidth
     nhalo=2*nsw

     call find_local_maxes ( terr_dev, ncube, nhalo, nsw, iopt_ridge_seed )

     call find_ridges ( terr_dev, terr, ncube, nhalo, nsw,&
          ncube_sph_smooth_coarse   , ncube_sph_smooth_fine,   &
          ldevelopment_diags, lregional_refinement=lregional_refinement,&
          rr_factor = rrfac  )

  endif

  !*********************************************************
  !
  ! Begin actual remapping calculations
  !
  !*********************************************************

  call allocate_target_vars (ntarget)
  if (lwrite_rrfac_to_topo_file) allocate (rrfac_target(ntarget))

  !*********************************************************************
  !      In the following loops "counti" is the index of a piece of 
  !      the "exchange grid"
  !********************************************************************

  !
  ! Sum exchange grid cells within each target
  ! grid cell
  !
  do counti=1,jall
     i    = weights_lgr_index_all(counti)
     wt = weights_all(counti,1)
     area_target        (i) = area_target(i) + wt
  end do
  write (6,*) "Min/MaX area_target",minval(area_target),MaXVal(area_target)
  write (6,*) "Min/MaX target_area",minval(target_area),MaXVal(target_area)

  !+++aRH
  if (llandfrac) then
     write (6,*) "Remapping landfrac"
     write (6,*) "Min/MaX before remap:", minval(landfrac), maxval(landfrac)
     landfrac_target = remap_field(landfrac,area_target,weights_eul_index_all(1:jall,:),weights_lgr_index_all(1:jall),&       
          weights_all(1:jall,:),ncube,jall,nreconstruction,ntarget)
     write (6,*) "Min/MaX after remap:", minval(landfrac_target), maxval(landfrac_target)
  end if
  !---aRH`

  write (6,*) "Remapping terrain"
  terr_target = remap_field(terr,area_target,weights_eul_index_all(1:jall,:),weights_lgr_index_all(1:jall),&
       weights_all(1:jall,:),ncube,jall,nreconstruction,ntarget)
  write (6,*) "Min/MaX:", minval(terr_target), maxval(terr_target)
  terr_uf_target = remap_field(terr,area_target,weights_eul_index_all(1:jall,:),weights_lgr_index_all(1:jall),&
       weights_all(1:jall,:),ncube,jall,nreconstruction,ntarget)
  write (6,*) "Min/MaX:", minval(terr_target), maxval(terr_target)


  write (6,*) "Remapping landm_coslat"
  landm_coslat_target = remap_field(landm_coslat,area_target,weights_eul_index_all(1:jall,:),weights_lgr_index_all(1:jall),&
       weights_all(1:jall,:),ncube,jall,nreconstruction,ntarget)
  write (6,*) "Min/MaX:", minval(landm_coslat_target), maxval(landm_coslat_target)

  write (6,*) "Remapping SGH30"
  sgh30_target = remap_field(var30,area_target,weights_eul_index_all(1:jall,:),weights_lgr_index_all(1:jall),&
       weights_all(1:jall,:),ncube,jall,nreconstruction,ntarget)
  write (6,*) "Min/MaX:", minval((sgh30_target)), maxval(sqrt(sgh30_target))
  deallocate(var30)
  deallocate(landm_coslat)

  write (6,*) "max difference between target grid area and remapping software area",&
       maxval(target_area-area_target)

  !
  ! Consistency checks  
  !
  do counti=1,ntarget
     if (terr_target(counti)>8848.0) then
        !
        ! max height is higher than Mount Everest
        !
        write (6,*) "FaTaL error: max height is higher than Mount Everest!"
        write (6,*) "terr_target",counti,terr_target(counti)
        write (6,*) "(lon,lat) locations of vertices of cell with excessive max height::"
        do i=1,ncorner
           write (6,*) target_corner_lon(i,counti),target_corner_lat(i,counti)
        end do
        stop
     else if (terr_target(counti)<-423.0) then
        !
        ! min height is lower than dead Sea
        !
        write (6,*) "FaTaL error: min height is lower than dead Sea!"
        write (6,*) "terr_target",counti,terr_target(counti)
        write (6,*) "(lon,lat) locations of vertices of cell with excessive min height::"
        do i=1,ncorner
           write (6,*) target_corner_lon(i,counti),target_corner_lat(i,counti)
        end do
        stop
     else 

     end if
  end do
  write (6,*) "Elevation data passed min/max consistency check!"
  write (6,*) " "

  !
  ! compute mean height (globally) of topography about sea-level for target grid unfiltered elevation
  !
  vol_target_un     = 0.0d0
  area_target_total = 0.0d0
  dO i=1,ntarget
     area_target_total = area_target_total+area_target(i)
     !    write (6,*) i,vol_target_un,terr_target(i),area_target(i)
     vol_target_un     = vol_target_un+terr_target(i)*area_target(i)
  end dO
  write (6,*) "mean height (globally) of topography about sea-level for target grid unfiltered elevation",&
       vol_target_un/area_target_total,vol_target_un,area_target_total

  !
  ! diagnostics
  !
  vol_source     = 0.0d0
  mea_source     = 0.0d0
  area_source    = 0.0d0
  !++jtb
  !   The two lines below moved up before call
  !   to smooth_intermediate_topo
  !
  !allocate ( da(ncube,ncube),stat=alloc_error )
  !call Equiangularallareas(ncube, da)
  !--jtb
  dO jp=1,6
     dO jy=1,ncube
        dO jx=1,ncube
           ii = (jp-1)*ncube*ncube+(jy-1)*ncube+jx
           vol_source = vol_source+terr(ii)*da(jx,jy)
           !+++aRH
           !if (landfrac(ii)>0.0d0) then
           mea_source   = mea_source  + terr(ii)*da(jx,jy)
           area_source  = area_source +          da(jx,jy)
           !else
           !end if
           !---aRH
        end dO
     end dO
  end dO
  write (6,*) "volume of input cubed-sphere terrain           :",vol_source
  write (6,*) "average elevation of input cubed-sphere terrain:",vol_source/(4.0d0*pi)
  write (6,*) "average elevation of input cubed-sphere terrain over land:",vol_source/area_source

  deallocate (da)
  !+++aRH
  deallocate(landfrac)
  !---aRH
  !
  ! compute variance with respect to cubed-sphere data
  !
  write (6,*) "Compute variance with respect to 3km cubed-sphere data: SGH"
  !  
  ! compute mean height (globally) of topography about sea-level for target grid filtered elevation
  !
  vol_target = 0d0
  dO i=1,ntarget
     vol_target = vol_target+terr_target(i)*area_target(i)
  end dO
  write (6,*) "mean height (globally) of topography about sea-level for target grid filtered elevation",&
       vol_target/area_target_total
  write (6,*) "percentage change in mean height between filtered and unfiltered elevations",&
       100.0d0*(vol_target-vol_target_un)/vol_target_un
  write (6,*) "percentage change in mean height between input cubed-sphere and unfiltered elevations",&
       100.0d0*(vol_source-vol_target_un)/vol_source    
  !
  ! done internal smoothing
  !
  terr_target   = 0d0
  sgh_target    = 0d0
  sgh_uf_target = 0d0
  
  do counti =  1, jall
     i    = weights_lgr_index_all(counti)!!
     !
     ix  = weights_eul_index_all(counti,1)
     iy  = weights_eul_index_all(counti,2)
     ip  = weights_eul_index_all(counti,3)
     !
     ! convert to 1d indexing of cubed-sphere
     !
     ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!

     wt = weights_all(counti,1)

     sgh_target  (i) = sgh_target  (i) + wt*(terr_dev(ix,iy,ip))**2/area_target(i)
     terr_target (i) = terr_target (i) + wt*(terr_sm(ix,iy,ip))/area_target(i) 
     sgh_uf_target(i) = sgh_uf_target(i)+wt*((terr_uf_target(i)-terr(ii))**2)/area_target(i)
  end do

  if (linterp_phis) then
     write (6,*) "bilinear interpolation of PHIS from intermediate cubed-sphere grid to target grid"
     call bilinear_interp(ncube,ntarget,target_center_lon,target_center_lat,terr_sm(1:ncube,1:ncube,:),terr_target)
  end if

  if (lfind_ridges) then
     allocate( uniqiC( ncube*ncube*6 ), uniqwC( ncube*ncube*6 ), wedgoC( ncube*ncube*6 )  )
     allocate( anisoC( ncube*ncube*6 ), anglxC( ncube*ncube*6 ), mxdisC( ncube*ncube*6 )  )
     allocate( hwdthC( ncube*ncube*6 ), clngtC( ncube*ncube*6 )  )
     allocate( riseqC( ncube*ncube*6 ), fallqC( ncube*ncube*6 )  )
     allocate( mxvrxC( ncube*ncube*6 ), mxvryC( ncube*ncube*6 )  )
     allocate( nodesC( ncube*ncube*6 ), cwghtC( ncube*ncube*6 ) )
     allocate( itrgtC( ncube*ncube*6 )  )

     call remapridge2cube( ncube,nhalo,nsw, &
          ncube_sph_smooth_coarse,ncube_sph_smooth_fine,lzero_negative_peaks, &
          ldevelopment_diags,lregional_refinement,  &
          rrfac, & 
          uniqiC, uniqwC, anisoC, &
          anglxC,mxdisC,hwdthC,clngtC, &
          riseqC,fallqC,mxvrxC,mxvryC, & 
          nodesC,cwghtC,wedgoC  )

     call remapridge2target(area_target,target_center_lon,target_center_lat, & 
          weights_eul_index_all(1:jall,:), & 
          weights_lgr_index_all(1:jall),weights_all(1:jall,:),ncube,jall,&
          nreconstruction,ntarget, &
          output_grid, ldevelopment_diags,&
          terr_dev, uniqiC, uniqwC, anisoC, &
          anglxC,mxdisC,hwdthC,clngtC, &
          riseqC,fallqC,mxvrxC,mxvryC,nodesC,cwghtC, itrgtC  )

     if (lridgetiles) then 
        call remapridge2tiles ( ntarget,ncube,jall,nreconstruction,     &
             area_target,target_center_lon,target_center_lat,         &
             weights_eul_index_all(1:jall,:), &
             weights_lgr_index_all(1:jall),  &
             weights_all(1:jall,:), &
             uniqiC,uniqwC,itrgtC,wedgoC )
     end if

     deallocate( uniqiC,uniqwC,anisoC,anglxC,mxdisC,hwdthC,clngtC, &
          riseqC,fallqC,mxvrxC,mxvryC,nodesC,cwghtC   )
  endif

  if (lwrite_rrfac_to_topo_file) then
     rrfac_target = 0.0_r8
     do counti=1,jall

        i    = weights_lgr_index_all(counti)!!
        !
        ix  = weights_eul_index_all(counti,1)
        iy  = weights_eul_index_all(counti,2)
        ip  = weights_eul_index_all(counti,3)
        !
        ! convert to 1d indexing of cubed-sphere
        !
        ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!

        wt = weights_all(counti,1)

        rrfac_target(i) = rrfac_target(i) + wt * rrfac(ix,iy,ip)/area_target(i)
     end do
  end if
  deallocate (weights_all,weights_eul_index_all)

  write (6,*) " !!!!!!!!  ******* maxval terr_target " , maxval(terr_target)

  !
  ! zero out small values
  !
  if (llandfrac) then
     do i=1,ntarget
        !+++ARH
        if (landfrac_target(i)<.001_r8)  landfrac_target(i) = 0d0
        !---ARH
     end do
     write (6,*) "min/max of landfac_target                : ",minval(landfrac_target    ),maxval(landfrac_target    )
  end if

  dO i=1,ntarget
     if (sgh_target(i)  <  0.5d0)  sgh_target(i)     = 0d0
     if (sgh30_target(i)<  0.5d0) sgh30_target(i)    = 0d0
  end dO
  sgh30_target = sqrt (sgh30_target)
  sgh_target   = sqrt (sgh_target)

  write (6,*) "min/max of terr source                   : ",minval(terr),maxval(terr)
  write (6,*) "min/max of terr_target                   : ",minval(terr_target    ),maxval(terr_target    )
  if (lwrite_rrfac_to_topo_file) then
     write (6,*) "min/max of rrfac                       : ",minval(rrfac_target),maxval(rrfac_target)
  end if
  write (6,*) "min/max of landm_coslat_target           : ",&
       minval(landm_coslat_target),maxval(landm_coslat_target)
  write (6,*) "min/max of var30_target                  : ",minval(sgh30_target   ),maxval(sgh30_target   )
  write (6,*) "min/max of var_target                    : ",minval(sgh_target   ),maxval(sgh_target   )

  write (6,*) " Model topo output file ",trim(output_fname)
  
  call wrtncdf_unstructured (ntarget, dom, idx, terr_target, landfrac_target, sgh_target, sgh30_target, &
       landm_coslat_target, target_center_lon, target_center_lat, target_area, &
       output_fname, lfind_ridges, command_line_arguments, &
       lwrite_rrfac_to_topo_file, rrfac_target, str_creator, area_target, llandfrac)
  
  deallocate (terr_target,landfrac_target,sgh30_target,sgh_target,landm_coslat_target)
  deallocate (target_center_lon, target_center_lat, target_area,area_target)
  if (lwrite_rrfac_to_topo_file) deallocate (rrfac_target,target_rrfac)
  
  !**********************************************************************************************************************************
  !
  ! dual grid physics grid configuration
  !
  !**********************************************************************************************************************************
  if (lphis_gll) then
     
     call read_target_grid (grid_descriptor_fname_gll, lregional_refinement, ltarget_latlon, lpole, &
          nlat, nlon, ntarget, ncorner, nrank, dom, idx, &
          target_corner_lon, target_corner_lat, target_center_lon, target_center_lat, target_area, target_rrfac)
     
     allocate (terr_target(ntarget))
     if (linterp_phis) then
        call bilinear_interp(ncube,ntarget,target_center_lon,target_center_lat,terr_sm(1:ncube,1:ncube,:),terr_target)
     else
        allocate (weights_all(jall_anticipated,nreconstruction),stat=alloc_error )
        allocate (weights_eul_index_all(jall_anticipated,3),stat=alloc_error )
        allocate (weights_lgr_index_all(jall_anticipated),stat=alloc_error )
        weights_all           = 0.0_r8
        weights_eul_index_all = 0
        weights_lgr_index_all = 0
        jall = jall_anticipated

        write (6,*) "Compute overlap weights for GLL grid: "
        call overlap_weights (weights_lgr_index_all, weights_eul_index_all, weights_all,&
             jall, ncube, ngauss, ntarget, ncorner, jmax_segments, target_corner_lon, target_corner_lat, nreconstruction, ldbg)

        allocate (area_target(ntarget))

        area_target = 0.0
        do counti=1,jall
           i    = weights_lgr_index_all(counti)
           wt   = weights_all(counti,1)
           area_target(i) = area_target(i) + wt
        end do

        write (6,*) "Remapping terrain"

        terr_target=0.0
        do counti=1,jall        
           i    = weights_lgr_index_all(counti)!!
           !
           ix  = weights_eul_index_all(counti,1)
           iy  = weights_eul_index_all(counti,2)
           ip  = weights_eul_index_all(counti,3)
           !
           ! convert to 1d indexing of cubed-sphere
           !
           ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!

           wt = weights_all(counti,1)

           terr_target (i) = terr_target (i) + wt * (terr_sm(ix,iy,ip))/area_target(i) 
        end do
        deallocate (weights_all,weights_eul_index_all)
     end if
     call wrtncdf_unstructured_append_phis (ntarget, terr_target, target_center_lon, target_center_lat, output_fname)
  end if
end program convterr

subroutine print_help
  write (6,*) "THIS NEEDS TO BE UPDATED"
  write (6,*) "Usage: cube_to_target [options] ..."
  write (6,*) "Options:"
  write (6,*) " "
  write (6,*) "    STaNdaRd OPTIONS"
  write (6,*) " "
  write (6,*) "-u, --name_email_of_creator=<string> [required] -> name and Email address of creator"
  write (6,*) "-g, --grid_descriptor_file=<string>  [required] -> ESMF or SCRIP compliant grid descriptor file"
  write (6,*) "-i, --intermediate_cs_name=<string>  [required] -> intermediate cubed-sphere topo file (usually ncube3000)"
  write (6,*) "-o, --output_grid=<string>           [required] -> identifier for output grid (e.g. ne30np4, f09x1.25, ..)"    
  write (6,*) "-c, --smoothing_scale=<real> (in km)            -> standard 'climate' smoothing is -c=100 for 1 degree"
  write (6,*) "-q, --output_data_directory=<string>            -> data output directory (default is output)"
  write (6,*) " "
  write (6,*) "    REGIONaL REFinEMENT OPTIONS"
  write (6,*) " "
  write (6,*) "-y, --rrfac_max=<int>   [required for var res]  -> maximum refinement level"
  write (6,*) "-v, --rrfac_manipulation                        -> enable manipulation of rrfac (used for spectral-element grids)"
  write (6,*) " "
  write (6,*) "    dISTaNCE WEIGHTEd SMOOTHER OPTIONS (dEFaULT IS LaPLaCIaN)"
  write (6,*) " "
  write (6,*) "-b, --distance_weighted_smoother                -> enable distance weighted smoother"
  write (6,*) "-p, --use_prefilter                             -> -b smoother option"
  write (6,*) "-f, --fine_radius=<int>                         -> "
  write (6,*) " "
  write (6,*) "   LaPLaCIaN SMOOTHER OPTIONS"
  write (6,*) " "
  write (6,*) "-m, --smoothing_over_ocean                      -> do not restrict smoother to only smooth over land"
  write (6,*) "-l, --smooth_phis_numcycle                      -> number of subcycles for Laplacian smoother (for stability)"
  write (6,*) " "
  write (6,*) "   MISCELLaNEOUS OPTIONS"
  write (6,*) " "
  write (6,*) "-r, --no_ridges                                 -> do not compute sub-grid-scale ridges"
  write (6,*) "-x, --stop_after_smooth                         -> stop after smoothing"
  write (6,*) "-1, --ridge2tiles                               -> ??? "
  write (6,*) "-z, --development_diags                         -> enable development diagnostics (for developers)"
  write (6,*) " "
  write (6,*) "-t, --smooth_topo_file=<string>                 -> use pre-compured smooth topo file"
  write (6,*) "-d, --write_rrfac_to_topo_file                  -> write rrfac to final topo file (usually for debugging)"
  write (6,*) " "
  write (6,*) "-n, --source_data_identifier=<string>           -> source topo data identifier"
  write (6,*) "                                                   (default gmted2010_modis_bedmachine)"
  write (6,*) " "
  write (6,*) "-a, --grid_descriptor_file_gll=<string>         -> grid descriptor file for dual grid configurations"
  write (6,*) "-s, --interpolate_phis                          -> bilinear interpolate PHIS to target grid (instead of remapping)"
  
  
  stop
end subroutine print_help

subroutine wrtncdf_unstructured (n, dom, idx, terr, landfrac, sgh, sgh30, landm_coslat, lon, lat, area, &
     output_fname, lfind_ridges, command_line_arguments, &
     lwrite_rrfac_to_topo_file, rrfac_target, str_creator, area_target, llandfrac)
  use shared_vars,  only: rad2deg
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shared_vars,  only: terr_uf_target, sgh_uf_target
  use ridge_ana,    only: nsubr, mxdis_target, mxvrx_target, mxvry_target, ang22_target, &
       anglx_target, aniso_target, anixy_target, hwdth_target, wghts_target, & 
       clngt_target, cwght_target, count_target,riseq_target,grid_length_scale, &
       fallq_target, isovar_target
  implicit none
#     include         <netcdf.inc>
  integer, intent(in) :: n
  real(r8),dimension(n), intent(in) :: terr, landfrac, sgh, sgh30, lon, lat, landm_coslat, area, area_target 
  integer, dimension(n), intent(in) :: dom, idx
  character(len=1024),   intent(in) :: output_fname
  logical,               intent(in) :: lfind_ridges
  character(len=1024),   intent(in) :: command_line_arguments
  logical,               intent(in) :: lwrite_rrfac_to_topo_file
  real(r8),dimension(n), intent(in) :: rrfac_target
  character(len=1024),   intent(in) :: str_creator
  logical,               intent(in) :: llandfrac

  integer            :: foutid     ! output file id
  integer            :: lonvid
  integer            :: latvid
  integer            :: terrid, areaid!,nid
  integer            :: domid, idxid, landfracid, sghid, sgh30id, landm_coslatid
  integer            :: mxdisid, ang22id, anixyid, anisoid, mxvrxid, mxvryid, hwdthid, wghtsid, anglxid, gbxarid
  integer            :: sghufid, terrufid, clngtid, cwghtid, countid, riseqid, fallqid, rrfacid, isovarid
  integer            :: THISID
  integer            :: status    ! return value for error control of netcdf routin
  integer, dimension(2) :: nid
  real(r8), parameter :: fillvalue = 1d36
  character(len=1024) :: str
  !
  !  Create NetCdF file for output
  !
  print *,"Create NetCdF file for output"
  status = nf_create (trim(output_fname), NF_64BIT_OFFSET , foutid)
  if (status /= NF_NOERR) call handle_err(status)
  !
  ! Create dimensions for output
  !
  status = nf_def_dim (foutid, 'ncol', n, nid(1))
  if (status /= NF_NOERR) call handle_err(status)
  !

  if (lfind_ridges) then 
     status = nf_def_dim (foutid, 'nrdg', nsubr, nid(2))
     if (status /= NF_NOERR) call handle_err(status)
  endif

  !
  !
  ! Create variable for output
  !
  print*, "Create variable for output"

  status = nf_def_var (foutid, 'PHIS', NF_DOUBLE, 1, nid(1), terrid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "PHIS error"
  end if

  status = nf_def_var (foutid, 'dom', NF_INT, 1, nid(1), domid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "dom error"
  end if

  status = nf_def_var (foutid, 'idx', NF_INT, 1, nid(1), idxid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "idx error"
  end if

  if (lwrite_rrfac_to_topo_file) then
     status = nf_def_var (foutid,'rrfac', NF_DOUBLE, 1, nid(1), rrfacid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "RRFaC error"
     end if
  end if

  if (llandfrac) then
     status = nf_def_var (foutid,'LANDFRAC', NF_DOUBLE, 1, nid(1), landfracid)
     if (status /= NF_NOERR) call handle_err(status)
  end if
  
  status = nf_def_var (foutid,'SGH', NF_DOUBLE, 1, nid(1), sghid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "SGH error"
  end if

  status = nf_def_var (foutid,'SGH30', NF_DOUBLE, 1, nid(1), sgh30id)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "SGH30 error"
  end if

  status = nf_def_var (foutid,'LANDM_COSLAT', NF_DOUBLE, 1, nid, landm_coslatid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "LANDM_COSLAT error"
  end if
  !
  status = nf_def_var (foutid,'area', NF_DOUBLE, 1, nid(1), areaid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "area error"
  end if

  status = nf_def_var (foutid, 'lat', NF_DOUBLE, 1, nid(1), latvid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "lat error"
  end if

  status = nf_def_var (foutid, 'lon', NF_DOUBLE, 1, nid(1), lonvid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "lon error"
  end if

  if (Lfind_ridges) then 
     status = nf_def_var (foutid,'ISOVAR', NF_DOUBLE, 1, nid(1), isovarid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "ISOVAR error"
     end if
     status = nf_def_var (foutid,'GBXAR', NF_DOUBLE, 1, nid(1), gbxarid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "GBXAR error"
     end if
     status = nf_def_var (foutid,'MXDIS', NF_DOUBLE, 2, nid , mxdisid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "MXDIS error"
     end if
     status = nf_def_var (foutid,'RISEQ', NF_DOUBLE, 2, nid , riseqid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "RISEQ error"
     end if
     status = nf_def_var (foutid,'FALLQ', NF_DOUBLE, 2, nid , fallqid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "FALLQ error"
     endif
     status = nf_def_var (foutid,'ANGLL', NF_DOUBLE, 2, nid , ang22id)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "ANGLL error"
     end if
     status = nf_def_var (foutid,'ANGLX', NF_DOUBLE, 2, nid , anglxid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "ANGLX error"
     endif
     status = nf_def_var (foutid,'ANISO', NF_DOUBLE, 2, nid , anisoid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "ANISO error"
     end if
     status = nf_def_var (foutid,'ANIXY', NF_DOUBLE, 2, nid , anixyid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "ANIXY error"
     end if
     status = nf_def_var (foutid,'HWDTH', NF_DOUBLE, 2, nid , hwdthid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "HWDTH error"
     end if
     status = nf_def_var (foutid,'CLNGT', NF_DOUBLE, 2, nid , clngtid)
     if (status /= NF_NOERR) then
        call handle_err(status)
        write (6,*) "CLNGT error"
     end if

  endif

  !
  ! Create attributes for output variables
  !
  status = nf_put_att_text   (foutid, terrid, 'long_name', 21, 'surface geopotential')
  status = nf_put_att_text   (foutid, terrid, 'units',     5, 'm2/s2')
  status = nf_put_att_double (foutid, terrid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, terrid, '_FillValue'   , nf_double, 1, fillvalue)

  status = nf_put_att_text   (foutid, domid, 'long_name', 6, 'domain')
  status = nf_put_att_text   (foutid, idxid, 'long_name', 10, 'node index')

  if (lwrite_rrfac_to_topo_file) then
     status = nf_put_att_text (foutid,rrfacid,'long_name', 17, 'refinement factor')
     status = nf_put_att_text (foutid,rrfacid,'units', 0, '')
     status = nf_put_att_double (foutid, rrfacid, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, rrfacid, '_FillValue'   , nf_double, 1, fillvalue)
  end if

  status = nf_put_att_double (foutid, sghid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, sghid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, sghid, 'long_name' , 48, &
       'standard deviation of 3km cubed-sphere elevation and target grid elevation')
  status = nf_put_att_text   (foutid, sghid, 'units'     , 1, 'm')


  status = nf_put_att_double (foutid, sgh30id, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, sgh30id, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, sgh30id, 'long_name' , 49, &
       'standard deviation of 30s elevation from 3km cubed-sphere cell average height')
  status = nf_put_att_text   (foutid, sgh30id, 'units'     , 1, 'm')


  status = nf_put_att_double (foutid, landm_coslatid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, landm_coslatid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, landm_coslatid, 'long_name' , 23, 'smoothed land fraction')
  status = nf_put_att_text   (foutid, landm_coslatid, 'filter'    , 4, 'none')
  if (llandfrac) then
     !+++aRH  
     status = nf_put_att_double (foutid, landfracid, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, landfracid, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, landfracid, 'long_name', 21, 'gridbox land fraction')
     !---aRH  
  end if

  status = nf_put_att_double (foutid, areaid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, areaid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, areaid, 'long_name' , 24, &
       'area of target grid cell')
  status = nf_put_att_text   (foutid, areaid, 'units'     , 1, 'm+2')

  status = nf_put_att_double (foutid, isovarid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, isovarid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, isovarid, 'long_name' , 30, &
       'residual variance after ridges')
  status = nf_put_att_text   (foutid, isovarid, 'units'     , 1, 'm+2')

  status = nf_put_att_text (foutid,latvid,'long_name', 8, 'latitude')
  if (status /= NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,latvid,'units', 13, 'degrees_north')
  if (status /= NF_NOERR) call handle_err(status)

  status = nf_put_att_text (foutid,lonvid,'long_name', 9, 'longitude')
  if (status /= NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,lonvid,'units', 12, 'degrees_east')
  if (status /= NF_NOERR) call handle_err(status)

  if (Lfind_ridges) then 
     THISID = mxdisid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 48, &
          'Obstacle height diagnosed by ridge-finding alg. ')
     status = nf_put_att_text   (foutid, THISID, 'units'     , 1, 'm')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = riseqid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 38, &
          'Rise to peak from left (ridge_finding)')
     !12345678901234567890123456789012345678901234567890
     !         10        20        30        40        
     status = nf_put_att_text   (foutid, THISID, 'units'     , 1, 'm')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = fallqid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 43, &
          'Fall from peak toward right (ridge_finding)')
     !12345678901234567890123456789012345678901234567890
     !         10        20        30        40        
     status = nf_put_att_text   (foutid, THISID, 'units'     , 1, 'm')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = ang22id
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 48, &
          'Ridge orientation clockwise from true north     ')
     !12345678901234567890123456789012345678901234567890
     !         10        20        30        40        
     status = nf_put_att_text   (foutid, THISID, 'units'     , 7, 'degrees')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = anglxid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 61, &
          'Ridge orientation clockwise from b-axis in cubed sphere panel')
     !1234567890123456789012345678901234567890123456789012345678901
     !         10        20        30        40        50        60
     status = nf_put_att_text   (foutid, THISID, 'units'     , 7, 'degrees')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = hwdthid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 21, &
          'Estimated Ridge width')
     !12345678901234567890123456789012345678901234567890
     !         10        20        30        40        
     status = nf_put_att_text   (foutid, THISID, 'units'     , 2, 'km')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = clngtid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 34, &
          'Estimated Ridge length along crest')
     !12345678901234567890123456789012345678901234567890
     !         10        20        30        40        
     status = nf_put_att_text   (foutid, THISID, 'units'     , 2, 'km')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = anixyid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 42, &
          'Variance ratio: cross/(cross+length) -wise')
     !12345678901234567890123456789012345678901234567890
     !         10        20        30        40        
     status = nf_put_att_text   (foutid, THISID, 'units'     , 1, '1')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = anisoid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 36, &
          'Variance fraction explained by ridge')
     !12345678901234567890123456789012345678901234567890
     !         10        20        30        40        
     status = nf_put_att_text   (foutid, THISID, 'units'     , 1, '1')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID = isovarid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 50, &
          'SQRT(Variance) from topo NOT represented by ridges')
     !12345678901234567890123456789012345678901234567890
     !         10        20        30        40        
     status = nf_put_att_text   (foutid, THISID, 'units'     , 1, '1')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')

     THISID=gbxarid
     status = nf_put_att_double (foutid, THISID, 'missing_value', nf_double, 1, fillvalue)
     status = nf_put_att_double (foutid, THISID, '_FillValue'   , nf_double, 1, fillvalue)
     status = nf_put_att_text   (foutid, THISID, 'long_name' , 46, &
          'angular area of target grid cell from scheme')
     !12345678901234567890123456789012345678901234
     !              10        20        30        40
     status = nf_put_att_text   (foutid, THISID, 'units'     , 7, 'm+2 m-2')
     status = nf_put_att_text   (foutid, THISID, 'filter'    , 4, 'none')
  end if

  !
  ! End define mode for output file
  !
  status = nf_enddef (foutid)
  if (status /= NF_NOERR) call handle_err(status)
  !
  ! Write variable for output
  !
  print*, "writing domain data", minval(dom), maxval(dom)
  status = nf_put_var_int (foutid, domid, dom)
  print*,"done writing domain data"

  print*, "writing node index data", minval(idx), maxval(idx)
  status = nf_put_var_int (foutid, idxid, idx)
  print*,"done writing node index  data"
  
  print*,"writing terrain data", minval(terr), maxval(terr)
  status = nf_put_var_double (foutid, terrid, terr*9.80616)

  !
  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing terrain data"

  if (lwrite_rrfac_to_topo_file) then
     status = nf_put_var_double (foutid, rrfacid, rrfac_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing rrfac data"
  end if
  
  if (llandfrac) then
     !+++aRH  
     print*,"writing landfrac data",minval(landfrac),maxval(landfrac)
     status = nf_put_var_double (foutid, landfracid, landfrac)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing landfrac data"
     !---aRH  
  end if
  print*,"writing sgh data",minval(sgh),maxval(sgh)
  status = nf_put_var_double (foutid, sghid, sgh)
  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing sgh data"

  print*,"writing sgh30 data",minval(sgh30),maxval(sgh30)
  status = nf_put_var_double (foutid, sgh30id, sgh30)
  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing sgh30 data"

  print*,"writing landm_coslat data",minval(landm_coslat),maxval(landm_coslat)
  status = nf_put_var_double (foutid, landm_coslatid, landm_coslat)
  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing sgh30 data"
  !
  print*,"writing area data",minval(area ),maxval(area)
  status = nf_put_var_double (foutid, areaid, area)
  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing area data"

  print*,"writing lat data"
  if (maxval(lat)<45.0) then
     status = nf_put_var_double (foutid, latvid, lat*rad2deg)
  else
     status = nf_put_var_double (foutid, latvid, lat)
  endif
  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing lat data"

  print*,"writing lon data"
  if (maxval(lon)<100.0) then
     status = nf_put_var_double (foutid, lonvid, lon*rad2deg)    
  else
     status = nf_put_var_double (foutid, lonvid, lon)
  end if

  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing lon data"

  if (Lfind_ridges) then 
     print*,"writing MXdIS data",minval(mxdis_target),maxval(mxdis_target)
     status = nf_put_var_double (foutid, mxdisid, mxdis_target )
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing MXdIS data"

     print*,"writing RISEQ  data",minval(riseq_target),maxval(riseq_target)
     status = nf_put_var_double (foutid, riseqid, riseq_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing RISEQ data"

     print*,"writing FaLLQ  data",minval(fallq_target),maxval(fallq_target)
     status = nf_put_var_double (foutid, fallqid, fallq_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing FaLLQ data"

     print*,"writing aNGLL data",minval(ang22_target),maxval(ang22_target)
     status = nf_put_var_double (foutid, ang22id, ang22_target )
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing aNGLL data"

     print*,"writing aNGLX  data",minval(anglx_target),maxval(anglx_target)
     status = nf_put_var_double (foutid, anglxid, anglx_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing aNGLX data"

     print*,"writing aNISO  data",minval(aniso_target),maxval(aniso_target)
     status = nf_put_var_double (foutid, anisoid, aniso_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing aNISO data"

     print*,"writing aNIXY  data",minval(anixy_target),maxval(anixy_target)
     status = nf_put_var_double (foutid, anixyid, anixy_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing aNIXY data"

     print*,"writing HWdTH  data",minval(hwdth_target),maxval(hwdth_target)
     status = nf_put_var_double (foutid, hwdthid, hwdth_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing HWdTH data"

     print*,"writing CLNGT  data",minval(clngt_target),maxval(clngt_target)
     status = nf_put_var_double (foutid, clngtid, clngt_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing CLNGT data"

     print*,"writing GBXaR  data",minval(area_target),maxval(area_target)
     status = nf_put_var_double (foutid, gbxarid, area_target)
     if (status /= NF_NOERR) call handle_err(status)
     print*,"done writing GBXaR data"

  endif
  !
  ! Close output file
  !
  print *,"close file"
  status = nf_close (foutid)
  if (status /= NF_NOERR) call handle_err(status)
end subroutine wrtncdf_unstructured

subroutine wrtncdf_unstructured_append_phis (n, terr, lon, lat, output_fname)
  !---aRH
  use shared_vars,  only: rad2deg
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
#     include         <netcdf.inc>
  integer,               intent(in) :: n
  real(r8),dimension(n), intent(in) :: terr,lon,lat
  character(len=1024),   intent(in) :: output_fname

  integer               :: foutid    ! output file id
  integer               :: lonvid
  integer               :: latvid
  integer               :: terrid
  integer               :: status    ! return value for error control of netcdf routin
  integer, dimension(2) :: nid
  real(r8), parameter   :: fillvalue = 1d36
  character(len=1024)   :: str

  !
  !  Create NetCdF file for output
  !
  print *,"Opening ",trim(output_fname)
  status = nf_open (trim(output_fname), nf_write, foutid)
  if (status /= NF_NOERR) call handle_err(status)

  status = nf_redef(foutid)
  if (status /= NF_NOERR) call handle_err(status)
  !
  ! Create dimensions for output
  !
  status = nf_def_dim (foutid, 'ncol_gll', n, nid(1))
  if (status /= NF_NOERR) call handle_err(status)
  !
  print *,"Create variable for output"
  status = nf_def_var (foutid,'PHIS_gll', NF_DOUBLE, 1, nid(1), terrid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "PHIS_gll error"
  end if

  status = nf_def_var (foutid,'lat_gll', NF_DOUBLE, 1, nid(1), latvid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "lat error"
  end if

  status = nf_def_var (foutid,'lon_gll', NF_DOUBLE, 1, nid(1), lonvid)
  if (status /= NF_NOERR) then
     call handle_err(status)
     write (6,*) "lon error"
  end if

  !
  ! Create attributes for output variables
  !
  status = nf_put_att_text (foutid,terrid,'long_name', 33, 'surface geopotential on GLL grid')
  status = nf_put_att_text (foutid,terrid,'units', 5, 'm2/s2')
  status = nf_put_att_double (foutid, terrid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, terrid, '_FillValue'   , nf_double, 1, fillvalue)

  status = nf_put_att_text (foutid,latvid,'long_name', 8, 'latitude')
  if (status /= NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,latvid,'units', 13, 'degrees_north')
  if (status /= NF_NOERR) call handle_err(status)

  status = nf_put_att_text (foutid,lonvid,'long_name', 9, 'longitude')
  if (status /= NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,lonvid,'units', 12, 'degrees_east')
  if (status /= NF_NOERR) call handle_err(status)
  !
  ! End define mode for output file
  !
  status = nf_enddef (foutid)
  if (status /= NF_NOERR) call handle_err(status)
  !
  ! Write variable for output
  !
  print*,"writing terrain data",minval(terr),maxval(terr)
  status = nf_put_var_double (foutid, terrid, terr*9.80616)
  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing terrain data"

  print*,"writing lat data"
  if (maxval(lat)<45.0) then
     status = nf_put_var_double (foutid, latvid, lat*rad2deg)
  else
     status = nf_put_var_double (foutid, latvid, lat)
  endif

  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing lat data"

  print*,"writing lon data"
  if (maxval(lon)<100.0) then
     status = nf_put_var_double (foutid, lonvid, lon*rad2deg)    
  else
     status = nf_put_var_double (foutid, lonvid, lon)
  end if

  if (status /= NF_NOERR) call handle_err(status)
  print*,"done writing lon data"

  !
  ! Close output file
  !
  print*, "close file"
  status = nf_close (foutid)
  if (status /= NF_NOERR) call handle_err(status)
end subroutine wrtncdf_unstructured_append_phis

!************************************************************************
!!handle_err
!************************************************************************
!
!!RoutinE:      handle_err
!!dESCRIPTION:  error handler
!--------------------------------------------------------------------------
subroutine handle_err(status)

  implicit         none

#     include          <netcdf.inc>

  integer          status

  if (status /= nf_noerr) then
     print *, nf_strerror(status)
     stop 'Stopped'
  endif

end subroutine handle_err


!*******************************************************************************
!  at this point mapping arrays are calculated
!
!      weights_lgr_index_all: dimension(JaLL). Index of target grid cell that contains
!                             current exchange grid cell
! 
!      weights_eul_index_all: dimension(JaLL,3). 3 indices of cubed-sphere grid cell that 
!                             contains current exchange grid cell:
!
!                                weights_eul_index_all(:,1) = x-index
!                                weights_eul_index_all(:,2) = y-index
!                                weights_eul_index_all(:,3) = panel/face number 1-6
!
!                             These are then converted to one-dimensional indices 
!                             for cubed sphere variables terr(n), ... etc. 
!
!      weights_all:           dimension(JaLL,nreconstrunction). Spherical area of
!                             exchange grid cell (steradians)
!
!********************************************************************************
subroutine overlap_weights (weights_lgr_index_all, weights_eul_index_all, weights_all, &
     jall, ncube, ngauss, ntarget, ncorner, jmax_segments, target_corner_lon, target_corner_lat, nreconstruction, ldbg)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use remap
  use shared_vars, only: progress_bar
  implicit NONE
  integer,                                     intent(inout):: jall !anticipated number of weights
  integer,                                     intent(in)   :: ncube, ngauss, ntarget, jmax_segments, ncorner, nreconstruction

  integer, dimension(jall,3),                  intent(out)  :: weights_eul_index_all
  real(r8), dimension(jall,nreconstruction)  , intent(out)  :: weights_all
  integer, dimension(jall)  ,                  intent(out)  :: weights_lgr_index_all

  real(r8), dimension(ncorner,ntarget),        intent(inout):: target_corner_lon, target_corner_lat
  LOGICaL,                                     intent(in)   :: ldbg

  integer,  dimension(9*(ncorner+1)) :: ipanel_tmp,ipanel_array
  real(r8), dimension(ncorner)  :: lat, lon
  real(r8), dimension(0:ncube+2):: xgno, ygno
  real(r8), dimension(0:ncorner+1) :: xcell, ycell

  real(r8), dimension(ngauss) :: gauss_weights, abscissae

  real(r8) :: da, tmp, alpha, beta
  real    (r8):: pi,piq,pih
  integer :: i, j,ncorner_this_cell,k,ip,ipanel,ii,jx,jy,jcollect
  integer :: alloc_error, ilon,ilat

  real    (r8) :: rad2deg
  real    (r8) :: deps

  real(r8), allocatable, dimension(:,:) :: weights
  integer , allocatable, dimension(:,:) :: weights_eul_index

  integer :: jall_anticipated, count

  pi = 4d0*datan (1d0)
  piq = pi/4d0
  pih = pi*0.5d0
  rad2deg = 180d0/pi

  deps = 10.0d0*pi/180.0_r8

  jall_anticipated = jall

  ipanel_array = -99
  !
  da = pih/dBLE(ncube)
  xgno(0) = -bignum
  dO i=1,ncube+1
     xgno(i) = tan (-piq+(i-1)*da)
  end dO
  xgno(ncube+2) = bignum
  ygno = xgno

  call glwp (ngauss, gauss_weights, abscissae)

  allocate (weights(jmax_segments,nreconstruction),stat=alloc_error )
  allocate (weights_eul_index(jmax_segments,2),stat=alloc_error )

  tmp = 0.0
  jall = 1
  dO i=1,ntarget
     if (MOd(i,10)==0)call progress_bar("# ", i, dBLE(100*i)/dBLE(ntarget))
     !
     !---------------------------------------------------          
     !
     ! determine how many vertices the cell has
     !
     !---------------------------------------------------
     !
     call remove_duplicates_latlon(ncorner,target_corner_lon(:,i),target_corner_lat(:,i),&
          ncorner_this_cell,lon,lat,1.0E-10)

     if (ldbg) THEN
        write (6,*) "number of vertices ",ncorner_this_cell
        write (6,*) "vertices locations lon,",lon(1:ncorner_this_cell)*rad2deg
        write (6,*) "vertices locations lat,",lat(1:ncorner_this_cell)*rad2deg
        dO j=1,ncorner_this_cell
           write (6,*) lon(j)*rad2deg, lat(j)*rad2deg
        end dO
        write (6,*) "  "
     end if
     !
     !---------------------------------------------------
     !
     ! determine how many and which panels the cell spans
     !
     !---------------------------------------------------          
     !
#ifdef old    
     dO j=1,ncorner_this_cell
        call CubedSphereaBPFromRLL(lon(j), lat(j), alpha, beta, ipanel_tmp(j), .TRUE.)
        if (ldbg) write (6,*) "ipanel for corner ",j," is ",ipanel_tmp(j)
     end dO
     ipanel_tmp(ncorner_this_cell+1) = ipanel_tmp(1)
     ! make sure to include possible overlap areas not on the face the vertices are located
     if (minval(lat(1:ncorner_this_cell))<-pi/6.0) THEN
        ! include South-pole panel in search
        ipanel_tmp(ncorner_this_cell+1) = 5
        if (ldbg) write (6,*)  "add panel 5 to search"
     end if
     if (maxval(lat(1:ncorner_this_cell))>pi/6.0) THEN
        ! include North-pole panel in search
        ipanel_tmp(ncorner_this_cell+1) = 6
        if (ldbg) write (6,*)  "add panel 6 to search"
     end if
     call remove_duplicates_integer(ncorner_this_cell+1,ipanel_tmp(1:ncorner_this_cell+1),&
          k,ipanel_array(1:ncorner_this_cell+1))
#endif
     !
     ! make sure to include possible overlap areas not on the face the vertices are located
     ! For example, a cell could be on panel 3 and 5 but have overlap area on panel 2
     count = 0
     do ilat=-1,1
        do ilon=-1,1
           dO j=1,ncorner_this_cell
              count=count+1
              call CubedSphereaBPFromRLL(lon(j)+ilon*deps, lat(j)+ilat*deps, alpha, beta, ipanel_tmp(count), .TRUE.)
           end dO
        end do
     end do

     !
     ! remove duplicates in ipanel_tmp
     !
     call remove_duplicates_integer(count,ipanel_tmp(1:count),&
          k,ipanel_array(1:count))
     !
     !---------------------------------------------------
     !
     ! loop over panels with possible overlap areas
     !
     !---------------------------------------------------          
     !
     dO ip = 1,k
        ipanel = ipanel_array(ip)
        dO j=1,ncorner_this_cell
           ii = ipanel
           call CubedSphereaBPFromRLL (lon(j), lat(j), alpha, beta, ii, .false.)            
           if (j==1) THEN
              jx = CEILinG((alpha + piq) / da)
              jy = CEILinG((beta  + piq) / da)
           end if
           xcell(ncorner_this_cell+1-j) = tan (alpha)
           ycell(ncorner_this_cell+1-j) = tan (beta)
        end dO
        xcell(0) = xcell(ncorner_this_cell)
        ycell(0) = ycell(ncorner_this_cell)
        xcell(ncorner_this_cell+1) = xcell(1)
        ycell(ncorner_this_cell+1) = ycell(1)

        jx = MaX(Min(jx,ncube+1),0)
        jy = MaX(Min(jy,ncube+1),0)

        call compute_weights_cell(xcell(0:ncorner_this_cell+1),ycell(0:ncorner_this_cell+1),&
             jx,jy,nreconstruction,xgno,ygno,&
             1, ncube+1, 1,ncube+1, tmp,&
             ngauss,gauss_weights,abscissae,weights,weights_eul_index,jcollect,jmax_segments,&
             ncube,0,ncorner_this_cell,ldbg,i)

        weights_all(jall:jall+jcollect-1,1:nreconstruction)  = weights(1:jcollect,1:nreconstruction)


        !weights_eul_index_all(jall:jall+jcollect-1,1:2) = weights_eul_index(1:jcollect,:)
        weights_eul_index_all(jall:jall+jcollect-1,  1) = weights_eul_index(1:jcollect,1)
        weights_eul_index_all(jall:jall+jcollect-1,  2) = weights_eul_index(1:jcollect,2)
        weights_eul_index_all(jall:jall+jcollect-1,  3) = ipanel
        weights_lgr_index_all(jall:jall+jcollect-1    ) = i

        jall = jall+jcollect
        if (jall>jall_anticipated) THEN
           write (6,*) "more weights than anticipated"
           write (6,*) "increase jall"
           stop
        end if
        if (ldbg) write (6,*) "jcollect",jcollect
     end dO
  end dO
  jall = jall-1
  write (6,*) "sum of all weights divided by surface area of sphere  =",tmp/(4.0*pi)
  write (6,*) "actual number of weights",jall
  write (6,*) "anticipated number of weights",jall_anticipated
  if (jall>jall_anticipated) THEN
     write (6,*) "anticipated number of weights < actual number of weights"
     write (6,*) "increase jall!"
     stop
  end if

  if (abs(tmp/(4.0*pi))-1.0>0.001) THEN
     write (6,*) "sum of all weights does not match the surface area of the sphere"
     write (6,*) "sum of all weights is : ",tmp
     write (6,*) "surface area of sphere: ",4.0*pi
     stop
  end if
end subroutine overlap_weights

subroutine bilinear_interp (ncube, ntarget, target_center_lon, target_center_lat, terr_cube, terr_target)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shared_vars,  only: progress_bar
  implicit NONE
  integer,                            intent(in) :: ncube, ntarget
  real(r8), dimension(ntarget),       intent(in) :: target_center_lon, target_center_lat
  real(r8), dimension(ncube,ncube,6), intent(in) :: terr_cube
  real(r8), dimension(ntarget),       intent(out):: terr_target

  real(r8)                           :: lat, lon
  real(r8), dimension(1:ncube+1)     :: xgno, ygno
  real(r8) :: da, alpha, beta, piq
  integer  :: i,ip,jx,jy,nhalo
  real(r8) :: x,y,x1,x2,y1,y2,w11,w12,w21,w22 !variables for bi-linear interpolation

  piq = datan (1d0)
  da = 2.0_r8*piq/dBLE(ncube)

  dO i=1,ncube+1
     xgno(i) = tan (-piq+(i-1)*da)
  end dO
  ygno = xgno
  dO i=1,ntarget
     if (MOd(i,10)==0)call progress_bar("# ", i, dBLE(100*i)/dBLE(ntarget))
     call CubedSphereaBPFromRLL(target_center_lon(i), target_center_lat(i), alpha, beta, ip, .TRUE.)
     jx = CEILinG((alpha + piq) / da)
     jy = CEILinG((beta  + piq) / da)
     jx = Min(MaX(1,jx),ncube-1); jy = Min(MaX(1,jy),ncube-1)

     x = tan(alpha);y = tan(beta)

     x1 = xgno(jx); x2 = xgno(jx+1); y1 = ygno(jy); y2 = ygno(jy+1)

     w11 = (x2-x )*(y2-y )/((x2-x1)*(y2-y1))
     w12 = (x2-x )*(y -y1)/((x2-x1)*(y2-y1))
     w21 = (x -x1)*(y2-y )/((x2-x1)*(y2-y1))
     w22 = (x -x1)*(y -y1)/((x2-x1)*(y2-y1))
     terr_target(i) = w11*terr_cube(jx  ,jy,ip)+w12*terr_cube(jx  ,jy+1,ip)+&
          w21*terr_cube(jx+1,jy,ip)+w22*terr_cube(jx+1,jy+1,ip)
  end dO
end subroutine bilinear_interp

!------------------------------------------------------------------------------
! subroutine CubedSphereaBPFromRLL
!
! description:
!   determine the (alpha,beta,panel) coordinate of a point on the sphere from
!   a given regular lat lon coordinate.
!
! Parameters:
!   lon - Coordinate longitude
!   lat - Coordinate latitude
!   alpha (out) - alpha coordinate
!   beta (out) - Beta coordinate
!   ipanel (out) - Face panel
!------------------------------------------------------------------------------
subroutine CubedSphereaBPFromRLL (lon, lat, alpha, beta, ipanel, ldetermine_panel)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shared_vars,  only: rotate_cube, pi, piq
  implicit NONE
  real(r8), intent(in)  :: lon, lat
  real(r8), intent(out) :: alpha, beta
  integer :: ipanel
  LOGICaL, intent(in) :: ldetermine_panel

  ! Local variables
  real(r8) :: xx, yy, zz, pm
  real(r8) :: sx, sy, sz
  integer  :: ix, iy, iz

  ! Translate to (x,y,z) space
  xx = COS(lon-rotate_cube) * COS(lat)
  yy = Sin(lon-rotate_cube) * COS(lat)
  zz = Sin(lat)

  pm = MaX(abs(xx), abs(yy), abs(zz))
  
  ! Check maximality of the x coordinate
  if (pm == abs(xx)) THEN
     if (xx > 0) THEN
        ix = 1
     ELSE
        ix = -1
     endif
  ELSE
     ix = 0
  endif

  ! Check maximality of the y coordinate
  if (pm == abs(yy)) THEN
     if (yy > 0) THEN
        iy = 1
     ELSE
        iy = -1
     endif
  ELSE
     iy = 0
  endif

  ! Check maximality of the z coordinate
  if (pm == abs(zz)) THEN
     if (zz > 0) THEN
        iz = 1
     ELSE
        iz = -1
     endif
  ELSE
     iz = 0
  endif
  
  ! Panel assignments
  if (ldetermine_panel) THEN
     if (iz  ==  1) THEN
        ipanel = 6; sx = yy; sy = -xx; sz = zz

     ELSEif (iz  == -1) THEN
        ipanel = 5; sx = yy; sy = xx; sz = -zz

     ELSEif ((ix == 1) .and. (iy /= 1)) THEN
        ipanel = 1; sx = yy; sy = zz; sz = xx

     ELSEif ((ix == -1) .and. (iy /= -1)) THEN
        ipanel = 3; sx = -yy; sy = zz; sz = -xx

     ELSEif ((iy == 1) .and. (ix /= -1)) THEN
        ipanel = 2; sx = -xx; sy = zz; sz = yy

     ELSEif ((iy == -1) .and. (ix /=  1)) THEN
        ipanel = 4; sx = xx; sy = zz; sz = -yy

     ELSE
        write (6,*) 'Fatal Error: CubedSphereaBPFromRLL failed'
        write (6,*) '(xx, yy, zz) = (', xx, ',', yy, ',', zz, ')'
        write (6,*) 'pm =', pm, ' (ix, iy, iz) = (', ix, ',', iy, ',', iz, ')'
        stop
     endif
  ELSE
     if (ipanel  ==  6) THEN
        sx = yy; sy = -xx; sz = zz
     ELSEif (ipanel  == 5) THEN
        sx = yy; sy = xx; sz = -zz
     ELSEif (ipanel == 1) THEN
        sx = yy; sy = zz; sz = xx        
     ELSEif (ipanel == 3) THEN
        sx = -yy; sy = zz; sz = -xx
     ELSEif (ipanel == 2) THEN
        sx = -xx; sy = zz; sz = yy
     ELSEif (ipanel == 4) THEN
        sx = xx; sy = zz; sz = -yy
     ELSE
        write (6,*) "ipanel out of range",ipanel
        stop
     end if
  end if

  ! Use panel information to calculate (alpha, beta) coords
  alpha = atan (sx / sz)
  beta = atan (sy / sz)
end subroutine CubedSphereaBPFromRLL

!------------------------------------------------------------------------------
! subroutine Equiangularallareas
!
! description:
!   Compute the area of all cubed sphere grid cells, storing the results in
!   a two dimensional array.
!
! Parameters: 
!   icube - Resolution of the cubed sphere
!   da (out) - Output array containing the area of all cubed sphere grid cells
!------------------------------------------------------------------------------
subroutine Equiangularallareas(icube, da)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  implicit NONE
  integer, intent(in)                           :: icube
  real (r8), dimension(icube,icube), intent(out) :: da

  ! Local variables
  integer                               :: k, k1, k2
  real (r8)                             :: a1, a2, a3, a4
  real (r8), dimension(icube+1,icube+1) :: ang
  real (r8), dimension(icube+1)         :: gp

  real    (r8):: pi, piq

  !#ifdef dBG 
  real (r8)   :: dbg1 !dBG
  !#endif


  pi = 4d0*datan (1d0)
  piq = pi/4d0
  ! Recall that we are using equi-angular spherical gridding
  !   Compute the angle between equiangular cubed sphere projection grid lines.
  dO k = 1, icube+1
     gp(k) = -piq + (pi/dBLE(2*(icube))) * dBLE(k-1)
  enddO

  dO k2=1,icube+1
     dO k1=1,icube+1
        ang(k1,k2) =aCOS(-Sin(gp(k1)) * Sin(gp(k2)))
     enddO
  enddO

  dO k2=1,icube
     dO k1=1,icube
        a1 =      ang(k1  , k2  )
        a2 = pi - ang(k1+1, k2  )
        a3 = pi - ang(k1  , k2+1)
        a4 =      ang(k1+1, k2+1)      
        ! area = r*r*(-2*pi+sum(interior angles))
        da(k1,k2) = -2d0*pi+a1+a2+a3+a4
     enddO
  enddO
end subroutine Equiangularallareas

!------------------------------------------------------------------------------
! subroutine CubedSphereRLLFromaBP
!
! description:
!   determine the lat lon coordinate of a point on a sphere given its
!   (alpha,beta,panel) coordinate.
!
! Parameters:
!   alpha - alpha coordinate
!   beta - Beta coordinate
!   panel - Cubed sphere panel id
!   lon (out) - Calculated longitude
!   lat (out) - Calculated latitude
!------------------------------------------------------------------------------
subroutine CubedSphereRLLFromaBP(alpha, beta, ipanel, lon, lat)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  use shared_vars, only: rotate_cube, pi, piq
  implicit NONE        
  real    (r8), intent(in)  :: alpha, beta
  integer     , intent(in)  :: ipanel
  real    (r8), intent(out) :: lon, lat        
  ! Local variables
  real    (r8) :: xx, yy, zz


  ! Convert to cartesian coordinates
  call CubedSphereXYZFromaBP(alpha, beta, ipanel, xx, yy, zz)        
  ! Convert back to lat lon
  lat = aSin(zz)
  if (xx==0.0.and.yy==0.0) THEN
     lon = 0.0
  else
     lon = aTaN2(yy, xx) +rotate_cube 
     if (lon<0.0) lon=lon+2d0*pi
     if (lon>2d0*pi) lon=lon-2d0*pi
  end if
end subroutine CubedSphereRLLFromaBP

!------------------------------------------------------------------------------
! subroutine CubedSphereXYZFromaBP
!
! description:
!   determine the Cartesian coordinate of a point on a sphere given its
!   (alpha,beta,panel) coordinate.
!
! Parameters:
!   alpha - alpha coordinate
!   beta - Beta coordinate
!   panel - Cubed sphere panel id
!   xx (out) - Calculated x coordinate
!   yy (out) - Calculated y coordinate
!   zz (out) - Calculated z coordinate
!------------------------------------------------------------------------------
subroutine CubedSphereXYZFromaBP(alpha, beta, ipanel, xx, yy, zz)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  implicit NONE
  real    (r8), intent(in)  :: alpha, beta
  integer     , intent(in)  :: ipanel
  real    (r8), intent(out) :: xx, yy, zz        
  ! Local variables
  real    (r8) :: a1, b1
  real    (r8) :: sx, sy, sz       

  ! Convert to Cartesian coordinates
  a1 = tan (alpha)
  b1 = tan (beta)

  sz = (1.0 + a1 * a1 + b1 * b1)**(-0.5)
  sx = sz * a1
  sy = sz * b1        
  ! Panel assignments
  if (ipanel == 6) THEN
     yy = sx; xx = -sy; zz = sz          
  ELSEif (ipanel == 5) THEN
     yy = sx; xx = sy; zz = -sz          
  ELSEif (ipanel == 1) THEN
     yy = sx; zz = sy; xx = sz          
  ELSEif (ipanel == 3) THEN
     yy = -sx; zz = sy; xx = -sz          
  ELSEif (ipanel == 2) THEN
     xx = -sx; zz = sy; yy = sz          
  ELSEif (ipanel == 4) THEN
     xx = sx; zz = sy; yy = -sz          
  ELSE
     write (6,*) 'Fatal Error: Panel out of range in CubedSphereXYZFromaBP'
     write (6,*) '(alpha, beta, panel) = (', alpha, ',', beta, ',', ipanel, ')'
     stop
  endif
end subroutine CubedSphereXYZFromaBP

subroutine remove_duplicates_integer(n_in,f_in,n_out,f_out)
  use shr_kind_mod, only: r8 => shr_kind_r8
  integer, intent(in) :: n_in
  integer,dimension(n_in), intent(in) :: f_in
  integer, intent(out) :: n_out
  integer,dimension(n_in), intent(out) :: f_out
  !
  ! local work space
  !
  integer :: k,i,j
  !
  ! remove duplicates in ipanel_tmp
  !
  k = 1
  f_out(1) = f_in(1)
  outer: do i=2,n_in
     do j=1,k
        !            if (f_out(j) == f_in(i)) then
        if (abs(f_out(j)-f_in(i))<1d-10) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     f_out(k) = f_in(i)
  end do outer
  n_out = k

  ! Hacke for Andes file
  do j=1,k
     if (f_out(j)>900) then
        n_out=k-1
     end if
  end do
end subroutine remove_duplicates_integer

subroutine remove_duplicates_latlon(n_in,lon_in,lat_in,n_out,lon_out,lat_out,tiny)
  use shr_kind_mod, only: r8 => shr_kind_r8
  integer, intent(in) :: n_in
  real(r8),dimension(n_in), intent(inout) :: lon_in,lat_in
  real, intent(in) :: tiny
  integer, intent(out) :: n_out
  real(r8),dimension(n_in), intent(out) :: lon_out,lat_out
  !
  ! local work space
  !
  integer  :: k,i,j
  real(r8) :: pi, pih

  pi = 4d0*datan (1d0)
  pih = pi*0.5d0
  !
  ! for pole points: make sure the longitudes are identical so that algorithm below works properly
  !
  do i = 2,n_in
     if (abs(lat_in(i)-pih)<tiny.or.abs(lat_in(i)+pih)<tiny) then 
        lon_in(i) = lon_in(i-1)    
        write (6,*) "pole fix"
     end if
  end do

  lon_out = -9999999.9
  lat_out = -9999999.9
  !
  k = 1
  lon_out(1) = lon_in(1)
  lat_out(1) = lat_in(1)
  outer: do i=2,n_in
     do j=1,k
        if (abs(lon_out(j)-lon_in(i))<tiny .and. abs(lat_out(j)-lat_in(i))<tiny) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     lon_out(k) = lon_in(i)
     lat_out(k) = lat_in(i)
  end do outer
  n_out = k
end subroutine remove_duplicates_latlon

