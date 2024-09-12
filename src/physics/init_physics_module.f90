!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: init_physics_module.f90
! Author: Gabrielle Ching-Johnson
! Date Revised: Nov 5 2023
! Description: Module contianing all initialization subroutines to initialize the Simple Physics
!               package or checkpoint and physics package plugins to be set during initialization.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module init_physics_mod
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use, intrinsic :: iso_c_binding, only : C_BOOL
  implicit none

  ! Physics test case arguments (set in test case main program)
  real(8) :: gas_molarmass, perihelion, aphelion, perihelion_day, obliquity, sea_surf, soil_surf, sea_inertia
  real(8) :: soil_inertia, sea_emissive, soil_emmisive, min_turbmix
  real(8) :: sw_atten, lw_atten
  
  real    :: sea_albedo, soil_albedo, emin_turb
  
  integer :: Nsoil                                                             ! number of soil layers
  
  logical :: radiation_mod, turbulence_mod, convecAdj_mod, seasons, diurnal

  logical(KIND=C_BOOL) :: physics_firstcall_flag = .true. ! flag for the physics package, true if call physics for 1st time
contains
  subroutine init_physics
    !-----------------------------------------------------------------------------------
    !
    !   Description: Sets necessary function pointers for the physics and initializes the
    !                main physics parameters.
    !
    !   Assumptions: Assumes that physics grid parameters has been initialized (ie called
    !                init_soild_grid(_default)) and the wavetrisk parameters have been
    !                initialized (ie called initialize(run_id)).
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    use iniphyparam_mod,   only : iniphyparam
    use single_column_mod, only : initialize_extra_levels
    use main_mod,          only : dt_new
    use read_param_mod
    use logging
    implicit none
    character(255) :: command, param_file

    ! Set physics function pointers (if using) ! it is here where I use soil_mod flag to set soil
    read_paramr_plugin => read_paramr
    read_parami_plugin => read_parami
    read_paramb_plugin => read_paramb
    flush_plugin       => flush_log_phys

    ! Write physics read in parameters to file (specific for each rank)
    write (param_file, '(A,I4.4)') trim(run_id)//'.physics_params.', rank
    call write_physics_params (9*rank, param_file)

    !Initialization of physics parameters
    open (unit=9*rank, file=trim (param_file), form="FORMATTED", action='READ')
    
    call iniphyparam (dt_new,  DAY, radius, grav_accel, R_d, c_p)
    
    close (9 * rank)

    ! Delete extra files
    write (command, '(a,a)') '\rm ', trim (param_file)
    call system (trim (command))

    ! Physics single column module extra levels initialization (as needs soil flag set in iniphyparam)
    call initialize_extra_levels (Nsoil+1)
  end subroutine init_physics

  subroutine write_physics_params (file_unit, file_params)
    !-----------------------------------------------------------------------------------
    !
    !   Description: Write desired physics parameters to a file, to be read by physics
    !                package during initialization by each rank.
    !
    !   Notes: Takes into account if mpi is being used, so file_params contains
    !           the file name that differs. Note that `mugaz` variable does not matter
    !           as set R manually to 278 in the initialization.
    !
    !   Used in: init_physics subroutine
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    integer      :: file_unit
    character(*) :: file_params
    logical      :: physics_write

    if (rank == 0) then
       physics_write = .true.
    else
       physics_write = .false.
    end if

    open (unit=file_unit, file=trim(file_params), form="FORMATTED", action='WRITE', status='REPLACE')

    ! Write physics parameters for reading by physics package
    write (file_unit,*) "planet_rat = ", radius
    write (file_unit,*) "g = ", grav_accel
    write (file_unit,*) "cpp = ", c_p
    write (file_unit,*) "mugaz = ", gas_molarmass 
    write (file_unit,*) "unjours = ", DAY
    write (file_unit,*) "year_day = ", int(YEAR / DAY)
    write (file_unit,*) "periheli = ", perihelion
    write (file_unit,*) "aphelie = ", aphelion
    write (file_unit,*) "peri_day = ", perihelion_day
    write (file_unit,*) "obliquit = ", obliquity
    write (file_unit,*) "Cd_mer = ", sea_surf
    write (file_unit,*) "Cd_ter = ", soil_surf
    write (file_unit,*) "I_mer = ", sea_inertia
    write (file_unit,*) "I_ter = ", soil_inertia
    write (file_unit,*) "alb_ter = ", sea_albedo
    write (file_unit,*) "alb_mer = ", soil_albedo
    write (file_unit,*) "emi_mer = ", sea_emissive
    write (file_unit,*) "emi_ter = ", soil_emmisive
    write (file_unit,*) "emin_turb = ", emin_turb
    write (file_unit,*) "lmixmin = ", min_turbmix
    write (file_unit,*) "coefvis = ", sw_atten
    write (file_unit,*) "coefir = ", lw_atten
    write (file_unit,*) "callrad = ", radiation_mod
    write (file_unit,*) "calldifv = ", turbulence_mod
    write (file_unit,*) "calladj = ", convecAdj_mod
    write (file_unit,*) "callsoil = ", soil_mod
    write (file_unit,*) "season = ", seasons
    write (file_unit,*) "diurnal = ", diurnal
    write (file_unit,*) "lverbose = ", physics_write
    write (file_unit,*) "period_sort = ", 1.
    close (file_unit)
  end subroutine write_physics_params

  subroutine init_soil_grid_default
    !-----------------------------------------------------------------------------------
    !
    !   Description: Initialize physics package with dummy longitude & latitude values
    !                and grid parameters. Also set number of soil levels and zmin to
    !                default value of the physics package (nsoilmx=10 currently)
    !
    !   Notes: To be used when number of soil levels not set in Simple_Physics.f90.
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    use comgeomfi, only : init_comgeomfi, nsoilmx

    real(8) :: lat(1), long(1)

    if (.not. soil_mod) then
       print*, 'STOP!! Cannot use default (init_soil_grid_default) to set zmin when soil_mod flag indicates false'
       stop
    end if

    ! Dummy latitude and longitude for initialization
    lat(1)  = 0d0
    long(1) = 0d0
    
    call init_comgeomfi(1, zlevels, long, lat) ! grid initialization for the physics

    ! Set the number of soil levels and zmin (lowest vertical level index)
    Nsoil = nsoilmx
    zmin  = - Nsoil
  end subroutine init_soil_grid_default

  subroutine init_soil_grid
    !-----------------------------------------------------------------------------------
    !
    !   Description: Initialize physics package with dummy longitude & latitude values
    !                and grid parameters. Also set number of soil levels and zmin.
    !
    !   Assumptions: Nsoil is set in Simple_Physics.f90 under the test case parameters.
    !                for case where soil model is turned on (ie soil_mod flag = true)
    !
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    use comgeomfi, only : init_comgeomfi

    real(8) :: lat(1), long(1)

    ! Dummy latitude and longitude for initialization
    lat(1)  = 0d0
    long(1) = 0d0

    ! Set the zmin (lowest vertical level index) (! See Assumptions)
    if (soil_mod) then
       zmin = - Nsoil
       call init_comgeomfi (1, zlevels, long, lat, Nsoil) ! grid initialization for the physics
    else
       Nsoil = 0
       zmin  = 0
       call init_comgeomfi (1, zlevels, long, lat)        ! grid initialization for the physics
    end if
  end subroutine init_soil_grid

  subroutine physics_checkpoint_restart
    !-----------------------------------------------------------------------------------
    !
    !   Description: Initialize physics cases for when checkpointing is used. Only to be
    !                called when cp_idx > 0.
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    use phyparam_mod, only : alloc, precompute, zday_last, icount

    real(8) :: day_fraction, nth_day

    ! Set flag for first call to physics false (thus the soil levels will get updated)
    physics_firstcall_flag = .false.

    ! Call allocation for physics call (usually done on the physics first call)
    call alloc(1, zlevels)

    ! Call precompute for physics call (usually done on the physics first call)
    call precompute ()

    ! Set the previous day in physics
    day_fraction = (time - dt) / DAY
    nth_day      = floor (day_fraction)
    day_fraction = day_fraction - nth_day

    zday_last = nth_day + day_fraction - dt/DAY
  end subroutine physics_checkpoint_restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Plugins !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_paramr (name, defval, val, comment)
    !-----------------------------------------------------------------------------------
    !
    !   Description: Physics plugin to read REALs from a file
    !
    !   Assumption: File is already open with unit number: 9*rank
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    character(*), intent(in)  :: name, comment
    real(8),      intent(in)  :: defval
    real(8),      intent(out) :: val

    character(15) :: param_name
    character     :: equal_sign

    ! Read line from input file
    read (9*rank,*) param_name, equal_sign, val

    if (trim (param_name) /= trim (name)) val = defval
  end subroutine read_paramr

  subroutine read_parami (name, defval, val, comment)
    !------------------- ----------------------------------------------------------------
    !
    !   Description: Physics plugin to read integers from a file
    !
    !   Assumption: File is already open with unit number: 9*rank
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    character(*), intent(IN)  :: name, comment
    integer,      intent(IN)  :: defval
    integer,      intent(OUT) :: val

    character(15) :: param_name
    character     :: equal_sign

    ! Read line from input file
    read (9*rank,*) param_name, equal_sign, val

    if (trim (param_name) /= trim (name)) val = defval
  end subroutine read_parami

  subroutine read_paramb(name, defval, val, comment)
    !-----------------------------------------------------------------------------------
    !
    !   Description: Physics plugin to read logicals from a file
    !
    !   Assumption: File is already open with unit number: 9*rank
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    character(*), intent(in)  :: name, comment
    logical,      intent(in)  :: defval
    logical,      intent(out) :: val

    character(15) :: param_name
    character     :: equal_sign

    ! Read line from input file
    read (9*rank,*) param_name, equal_sign, val

    if (trim (param_name) /= trim (name)) val = defval
  end subroutine read_paramb

  subroutine flush_log_phys (lev, taglen, tag, buflen, bufsize, buf) BIND(C)
    !-----------------------------------------------------------------------------------
    !
    !   Description: Physics plugin to flush the log buffer.
    !
    !   Notes: The goal was to only have the master (rank 0) flush (print to terminal)
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    use logging,                     only : dbtag
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int
    
    integer(c_int),         intent(in), value :: lev, taglen, buflen, bufsize
    character(KIND=c_char), intent(in)        :: tag(taglen), buf(buflen, bufsize)

    character(buflen+1) :: line
    character(100)      :: prefix
    integer             :: i

    if (rank == 0) then
       write (prefix,*) '[', dbtag(lev), ' ', tag, ']'
       do i = 1, bufsize
          write (line,*) buf(:,i)
          write (*,*) trim (prefix), trim (line)
       end do
    end if
  end subroutine flush_log_phys
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module init_physics_mod
