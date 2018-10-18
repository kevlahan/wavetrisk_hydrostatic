module test_case_mod
  ! Module file for flat_projection_data
  use shared_mod
  use domain_mod
  use comm_mpi_mod
  implicit none
  integer                              :: check_end, check_start, cp_2d, N, save_zlev
  real(8)                              :: initotalmass, mass_error, totalmass
  real(8), allocatable, dimension(:,:) :: threshold_def
  
  ! DCMIP2012c4
  real(8) :: eta_0, u_0 
  ! DCMIP2008c5
  real(8) :: d2, h_0, lat_c, lon_c
contains
  real(8) function surf_geopot (x_i)
    ! Surface geopotential
    implicit none
    Type(Coord) :: x_i
    real(8)     :: c1, cs2, sn2, lon, lat, rgrc

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)
    cs2 = cos(lat)**2
    sn2 = sin(lat)**2
    
    if (trim (test_case) == "DCMIP2012c4") then
       c1 = u_0*cos((1.0_8-eta_0)*MATH_PI/2)**1.5

       surf_geopot = c1 * (c1 * (-2*sn2**3*(cs2 + 1/3.0_8) + 10/63.0_8)  + &
        radius*omega*(8/5.0_8*cs2**1.5*(sn2 + 2/3.0_8) - MATH_PI/4))
    elseif (trim (test_case) == "DCMIP2008c5") then
       rgrc = radius*acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))

       surf_geopot = grav_accel*h_0*exp__flush (-rgrc**2/d2)
    elseif (trim (test_case) == "Held_Suarez") then
       surf_geopot = 0.0_8
    else
       write(6,'(A)') "Test case not supported"
       stop
    end if
  end function surf_geopot

 subroutine initialize_a_b_vert
    implicit none
    integer :: k

    ! Allocate vertical grid parameters
    allocate (a_vert(1:zlevels+1), b_vert(1:zlevels+1))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (uniform) then
       do k = 1, zlevels+1
          a_vert(k) = dble(k-1)/dble(zlevels) * press_infty/ref_press
          b_vert(k) = 1.0_8 - dble(k-1)/dble(zlevels)
       end do
    else
       if (zlevels==18) then
          a_vert=(/0.00251499_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
               0.07869805_8, 0.07463175_8, 0.06955308_8, 0.06339061_8, 0.05621774_8, 0.04815296_8, &
               0.03949230_8, 0.03058456_8, 0.02193336_8, 0.01403670_8, 0.007458598_8, 0.002646866_8, &
               0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.03756984_8, 0.08652625_8, 0.1476709_8, 0.221864_8, &
               0.308222_8, 0.4053179_8, 0.509588_8, 0.6168328_8, 0.7209891_8, 0.816061_8, 0.8952581_8, &
               0.953189_8, 0.985056_8, 1.0_8 /)
       elseif (zlevels==26) then
          a_vert=(/0.002194067_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
               0.07851243_8, 0.07731271_8, 0.07590131_8, 0.07424086_8, 0.07228744_8, 0.06998933_8, 0.06728574_8, 0.06410509_8, &
               0.06036322_8, 0.05596111_8, 0.05078225_8, 0.04468960_8, 0.03752191_8, 0.02908949_8, 0.02084739_8, 0.01334443_8, &
               0.00708499_8, 0.00252136_8, 0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.01505309_8, 0.03276228_8, 0.05359622_8, &
               0.07810627_8, 0.1069411_8, 0.1408637_8, 0.1807720_8, 0.2277220_8, 0.2829562_8, 0.3479364_8, 0.4243822_8, &
               0.5143168_8, 0.6201202_8, 0.7235355_8, 0.8176768_8, 0.8962153_8, 0.9534761_8, 0.9851122_8, 1.0_8 /)
       elseif (zlevels==30) then
          a_vert = (/ 0.00225523952394724, 0.00503169186413288, 0.0101579474285245, 0.0185553170740604, 0.0306691229343414, &
               0.0458674766123295, 0.0633234828710556, 0.0807014182209969, 0.0949410423636436, 0.11169321089983, & 
               0.131401270627975, 0.154586806893349, 0.181863352656364, 0.17459799349308, 0.166050657629967, &
               0.155995160341263, 0.14416541159153, 0.130248308181763, 0.113875567913055, 0.0946138575673103, &
               0.0753444507718086, 0.0576589405536652, 0.0427346378564835, 0.0316426791250706, 0.0252212174236774, &
               0.0191967375576496, 0.0136180268600583, 0.00853108894079924, 0.00397881818935275, 0.0, 0.0 /)
          b_vert = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0393548272550106, &
               0.0856537595391273, 0.140122056007385, 0.204201176762581, 0.279586911201477, 0.368274360895157,  &
               0.47261056303978, 0.576988518238068, 0.672786951065063, 0.753628432750702, 0.813710987567902, &
               0.848494648933411, 0.881127893924713, 0.911346435546875, 0.938901245594025, 0.963559806346893, &
               0.985112190246582, 1.0 /)
       elseif (zlevels==49) then
          a_vert=(/0.002251865_8, 0.003983890_8, 0.006704364_8, 0.01073231_8, 0.01634233_8, 0.02367119_8, &
               0.03261456_8, 0.04274527_8, 0.05382610_8, 0.06512175_8, 0.07569850_8, 0.08454283_8, &
               0.08396310_8, 0.08334103_8, 0.08267352_8, 0.08195725_8, 0.08118866_8, 0.08036393_8, &
               0.07947895_8, 0.07852934_8, 0.07751036_8, 0.07641695_8, 0.07524368_8, 0.07398470_8, &
               0.07263375_8, 0.07118414_8, 0.06962863_8, 0.06795950_8, 0.06616846_8, 0.06424658_8, &
               0.06218433_8, 0.05997144_8, 0.05759690_8, 0.05504892_8, 0.05231483_8, 0.04938102_8, &
               0.04623292_8, 0.04285487_8, 0.03923006_8, 0.03534049_8, 0.03116681_8, 0.02668825_8, &
               0.02188257_8, 0.01676371_8, 0.01208171_8, 0.007959612_8, 0.004510297_8, 0.001831215_8, &
               0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
               0.006755112_8, 0.01400364_8, 0.02178164_8, 0.03012778_8, 0.03908356_8, 0.04869352_8, &
               0.05900542_8, 0.07007056_8, 0.08194394_8, 0.09468459_8, 0.1083559_8, 0.1230258_8, &
               0.1387673_8, 0.1556586_8, 0.1737837_8, 0.1932327_8, 0.2141024_8, 0.2364965_8, &
               0.2605264_8, 0.2863115_8, 0.3139801_8, 0.3436697_8, 0.3755280_8, 0.4097133_8, &
               0.4463958_8, 0.4857576_8, 0.5279946_8, 0.5733168_8, 0.6219495_8, 0.6741346_8, &
               0.7301315_8, 0.7897776_8, 0.8443334_8, 0.8923650_8, 0.9325572_8, 0.9637744_8, &
               0.9851122_8, 1.0_8/)
       else
          write(0,*) "For this number of zlevels, no rule has been defined for a_vert and b_vert"
          stop
       end if

       ! DCMIP order is opposite to ours
       a_vert = a_vert(zlevels+1:1:-1)
       b_vert = b_vert(zlevels+1:1:-1)
    end if
    
    ! Set pressure at infinity
    press_infty = a_vert(zlevels+1)*ref_press ! note that b_vert at top level is 0, a_vert is small but non-zero

    ! Set mass coefficients
    b_vert_mass = b_vert(1:zlevels)-b_vert(2:zlevels+1)
    a_vert_mass = ((a_vert(1:zlevels)-a_vert(2:zlevels+1))*ref_press + b_vert_mass*press_infty)/grav_accel
  end subroutine initialize_a_b_vert

  subroutine read_test_case_parameters (filename)
    implicit none
    character(*)   :: filename
    
    integer        :: fid = 500
    real(8)        :: press_save
    character(255) :: varname

    open (unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, run_id
    read (fid,*) varname, check_start
    read (fid,*) varname, check_end
    read (fid,*) varname, cp_2d
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, uniform
    read (fid,*) varname, level_save
    read (fid,*) varname, N
    read (fid,*) varname, press_save
    close(fid)

    if (N > 2000) then
       write (6,'(A,i5,A)') "N = ", N, " is too large. Maximum allowed value of N = 2000."
       write (6,'(A)') "Modify export_2d in io.f90 if necessary."
       stop
    end if

    if (rank==0) then
       write (6,'(A,A)')      "test_case           = ", trim (test_case)
       write (6,'(A,A)')      "run_id              = ", trim (run_id)
       write (6,'(A,i4)')     "first file          = ", check_start
       write (6,'(A,i4)')     "last file           = ", check_end
       write (6,'(A,i4)')     "file to save as 2d  = ", cp_2d
       write (6,'(A,i3)')     "min_level           = ", min_level
       write (6,'(A,i3)')     "max_level           = ", max_level
       write (6,'(A,i3)')     "zlevels             = ", zlevels
       write(6,'(A,L1)')      "uniform             = ", uniform
       write (6,'(A,i3)')     "level_save          = ", level_save
       write (6,'(A,i5)')     "N                   = ", N
       write (6,'(A,es10.4)') "pressure_save (hPa) = ", press_save
       write (6,*) ' '
    end if
    
    allocate (pressure_save(1))
    pressure_save(1) = 100*press_save
  end subroutine read_test_case_parameters

  subroutine apply_initial_conditions
    implicit none
    integer :: k, l

    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (init_sol, l, k, -1, 1)
       end do
    end do
  end subroutine apply_initial_conditions

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Dummy routine
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
  end subroutine init_sol

  subroutine set_thresholds
    ! Dummy routine
  end subroutine set_thresholds

  subroutine initialize_thresholds
    ! Set default thresholds based on dimensional scalings of norms
    implicit none

    allocate (threshold(S_MASS:S_VELO,1:zlevels));     threshold     = 0.0_8
    allocate (threshold_def(S_MASS:S_VELO,1:zlevels)); threshold_def = 0.0_8
  end subroutine initialize_thresholds

  subroutine initialize_dt_viscosity 
    implicit none
    allocate (viscosity_divu(1:zlevels)); viscosity_divu = 0.0_8
  end subroutine initialize_dt_viscosity

  subroutine set_save_level
    implicit none
    save_zlev = zlevels
  end subroutine set_save_level

  subroutine dump (fid)
    implicit none
    integer :: fid

    write (fid) itime
    write (fid) iwrite
    write (fid) threshold
  end subroutine dump

  subroutine load (fid)
    implicit none
    integer :: fid

    read (fid) itime
    read (fid) iwrite
    read (fid) threshold
  end subroutine load
end module test_case_mod
