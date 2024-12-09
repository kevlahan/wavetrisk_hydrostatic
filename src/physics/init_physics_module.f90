!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: init_physics_module.f90
! Author: Gabrielle Ching-Johnson and Nicholas Kevlahan
! Date Revised: 2024-12-06
! Description: all initialization subroutines for Simple Physics package 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module init_physics_mod
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use, intrinsic :: iso_c_binding, only : C_BOOL
  implicit none
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Default values values physics Model Parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Sub-model activation
  logical :: convecAdj_model  = .true.                          ! convective adjustment module is on
  logical :: diurnal          = .true.                          ! diurnal cycle flag
  logical :: radiation_model  = .true.                          ! radiation module is on
  logical :: soil_model       = .true.                          ! soil model
  logical :: turbulence_model = .true.                          ! vertical diffusion module is on

  ! Physics Package planet test case parameters: set to Earth values by default
  real(8) :: gas_molarmass   = 28.9702532d0                     ! molar mass of main gain (used to set ideal gas const in pacakage)
  real(8) :: perihelion      = 150d0                            ! planet perihelion distance [MMkm]
  real(8) :: aphelion        = 150d0                            ! planet aphelion distance   [MMkm]
  real(8) :: perihelion_day  = 0d0                              ! perihelion day
  real(8) :: obliquity       = 23.5d0                           ! planet axial tilt/obliquity
  real(8) :: sea_surf        = 0.01d0                           ! sea surface roughness length scale  [m]
  real(8) :: soil_surf       = 0.01d0                           ! soil surface roughness length scale [m]
  real(8) :: sea_inertia     = 3000d0                           ! sea thermal  inertia [J/(m^3 K)]
  real(8) :: soil_inertia    = 3000d0                           ! soil thermal inertia [J/(m^3 K)]
  real(8) :: sea_emissive    = 1d0                              ! sea emissivity
  real(8) :: soil_emmisive   = 1d0                              ! soil emissivity
  real(8) :: min_turbmix     = 100d0                            ! minimum turbulent mixing length [m]
  real(8) :: sw_atten        = 0.99d0                           ! attenuation of shortwave radiation coefficient
  real(8) :: lw_atten        = 0.08d0                           ! attenuation of longwave radiation coefficient

  ! Single precision parameters
  real(8) :: sea_albedo      = 0.112e0                          ! sea albedo
  real(8) :: soil_albedo     = 0.112e0                          ! soil albedo
  real(8) :: Emin_turb       = 1e-16                            ! minimum turbulent kinetic energy

  logical(KIND=C_BOOL) :: physics_firstcall_flag = .true. ! flag for the physics package, true if call physics for 1st time
contains
  subroutine init_physics 
    !-----------------------------------------------------------------------------------
    !
    !   Description: Sets necessary function pointers for the physics and initializes the
    !                main physics parameters.
    !
    !   Assumptions: Assumes that physics grid parameters has been initialized (ie called
    !                init_soil_grid(_default)) and the wavetrisk parameters have been
    !                initialized (ie called initialize(run_id)).
    !
    !   if change dt_Phys and iradia in iniphyparam_mod.F90 to not compute every time step
    !   need to use fixed time step in test case
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    use iniphyparam_mod,   only : iniphyparam
    use single_column_mod, only : initialize_extra_levels
    use shared_mod,        only : dt
    use phyparam_mod,      only : alloc, precompute, zday_last, icount
    use read_param_mod
    use logging
    implicit none
    real(8) :: day_fraction, nth_day
    character(255) :: command, param_file

    ! Set physics function pointers (if using) ! it is here where I use soil_mod flag to set soil
    read_paramr_plugin => read_paramr
    read_parami_plugin => read_parami
    read_paramb_plugin => read_paramb
    flush_plugin       => flush_log_phys

    ! Write physics read in parameters to file (specific for each rank)
    write (param_file, '(a,i4.4)') trim(run_id)//'.physics_params.', rank
    call write_physics_params (9*rank, param_file)

    ! Initialization of physics parameters
    open (unit=9*rank, file=trim(param_file), form="FORMATTED", action='READ')
    
    call iniphyparam (dt, DAY, radius, grav_accel, R_d, c_p)
    
    close (9*rank)

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
    write (file_unit,*) "planet_rat  = ", radius
    write (file_unit,*) "g           = ", grav_accel
    write (file_unit,*) "cpp         = ", c_p
    write (file_unit,*) "mugaz       = ", gas_molarmass 
    write (file_unit,*) "unjours     = ", DAY
    write (file_unit,*) "year_day    = ", int (YEAR / DAY)
    write (file_unit,*) "periheli    = ", perihelion
    write (file_unit,*) "aphelie     = ", aphelion
    write (file_unit,*) "peri_day    = ", perihelion_day
    write (file_unit,*) "obliquit    = ", obliquity
    write (file_unit,*) "Cd_mer      = ", sea_surf
    write (file_unit,*) "Cd_ter      = ", soil_surf
    write (file_unit,*) "I_mer       = ", sea_inertia
    write (file_unit,*) "I_ter       = ", soil_inertia
    write (file_unit,*) "alb_ter     = ", sea_albedo
    write (file_unit,*) "alb_mer     = ", soil_albedo
    write (file_unit,*) "emi_mer     = ", sea_emissive
    write (file_unit,*) "emi_ter     = ", soil_emmisive
    write (file_unit,*) "emin_turb   = ", emin_turb
    write (file_unit,*) "lmixmin     = ", min_turbmix
    write (file_unit,*) "coefvis     = ", sw_atten
    write (file_unit,*) "coefir      = ", lw_atten
    write (file_unit,*) "callrad     = ", radiation_model
    write (file_unit,*) "calldifv    = ", turbulence_model
    write (file_unit,*) "calladj     = ", convecAdj_model
    write (file_unit,*) "callsoil    = ", soil_model
    write (file_unit,*) "diurnal     = ", diurnal
    write (file_unit,*) "lverbose    = ", physics_write
    write (file_unit,*) "period_sort = ", 1.0
    close (file_unit)
  end subroutine write_physics_params

   subroutine init_soil_grid
    !-----------------------------------------------------------------------------------
    !
    !   Description: Initialize physics package with dummy longitude & latitude values
    !                and grid parameters. Also set number of soil levels and zmin.
    !
    !   Assumptions: Nsoil is set in test case (e.g. climate.f90).
    !
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    use comgeomfi, only : init_comgeomfi, nsoilmx

    real(8) :: lat(1), long(1)

    ! Dummy latitude and longitude for initialization
    lat(1)  = 0d0
    long(1) = 0d0

    call init_comgeomfi (1, zlevels, long, lat)
    if (Nsoil /= 0) then
       soil_model = .true.
       nsoilmx    =   Nsoil
       zmin       = - Nsoil
    else ! Nsoil = 0 means there is no soil model
       soil_model = .false.
       nsoilmx    = 1
       zmin       = 0
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
    use phyparam_mod, only : alloc, precompute

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

  subroutine read_paramb (name, defval, val, comment)
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
