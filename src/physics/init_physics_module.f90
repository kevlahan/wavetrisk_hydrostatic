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
   use, intrinsic :: ISO_C_BINDING, only : C_BOOL ! require c booleans for physics
   implicit none

   ! Physics test case arguments (set in test case main program)
   real(8) :: e_thick
   real(8) :: gas_molarmass, perihelion, aphelion, perihelion_day, obliquity, sea_surf, soil_surf, sea_interia
   real(8) :: soil_interia, sea_emissive, soil_emmisive, min_turbmix
   real(8) :: sw_atten, lw_atten
   real    :: sea_albedo, soil_albedo, emin_turb
   integer :: Nsoil     ! Number of soil layers
   logical :: radiation_mod, turbulence_mod, convecAdj_mod, seasons, diurnal
   
   logical(kind=C_BOOL) :: physics_firstcall_flag = .true. ! flag for the physics package, true if call physics for 1st time
   
contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialization Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      use read_param_mod
      use logging
      use main_mod,          only : dt_new ! dynamics 
      implicit none
      
      character(255) :: command, param_file

      ! Default physics package parameters
      radiation_mod  = .true.                         ! (T) radiation module is on
      turbulence_mod = .true.                         ! (T) vertical diffusion module is on
      convecAdj_mod  = .true.                         ! (T) convective adjustment module is on
      soil_mod       = .true.                         ! (T) soil module is on

      ! Velocity (vertical mixing)
      emin_turb      = 1d-16                          ! minimum turbulent kinetic energy
      gas_molarmass  = 28.9702532d0                   ! molar mass of main gain (used to set ideal gas const in pacakage)
      min_turbmix    = 100d0                          ! minimum turbulent mixing length [m]

      ! Astronomical (orbit)
      aphelion       = 150d0                          ! planet aphelion distance [MMkm]
      diurnal        = .true.                         ! diurnal cycle flag
      obliquity      = 0d0                            ! planet axial tilt/obliquity
      perihelion     = 150d0                          ! planet perihelion distance [MMkm]
      perihelion_day = 0d0                            ! perihelion day
      seasons        = .false.                        ! seasons flag !! NOT USED !!

      ! Ocean
      sea_emissive   = 1d0                            ! sea emissivity
      sea_surf       = 0.01d0                         ! sea surface roughness length scale [m]
      sea_interia    = 3000d0                         ! sea thermal inertia [J/m^3K]

      ! Soil
      Nsoil          = 10                             ! number of soil layers
      soil_interia   = 3000.                          ! soil thermal inertial [J/m^3K]
      sea_albedo     = 0.112d0                        ! sea albedo
      soil_albedo    = 0.112d0                        ! soil albedo
      soil_emmisive  = 1d0                            ! soil emissivity
      soil_surf      = 0.01d0                         ! soil surface roughness length scale [m]

      ! Radiation
      lw_atten       = 0.08d0                         ! attenuation of longwave radiation coefficient
      sw_atten       = 0.99d0                         ! attenuation of shortwave radiation coefficient

      ! Set physics function pointers (if using) ! it is here where I use soil_mod flag to set soil
      read_paramr_plugin => read_paramr
      read_parami_plugin => read_parami
      read_paramb_plugin => read_paramb
      flush_plugin       => flush_log_phys

      ! Write physics read in parameters to file (specific for each rank)
      write(param_file, '(a,i4.4)') trim(run_id)//'.physics_params.', rank
      call write_physics_params (9*rank,param_file)

      ! Call initialization of physics parameters
      open (unit=9*rank, file=trim (param_file), form="FORMATTED", action='READ')
      call iniphyparam (dt_new,  DAY, radius, grav_accel, R_d, c_p)
      close (9*rank)

      ! Delete extra files
      write (command, '(a,a)') '\rm ', trim (param_file)
      call system (trim (command))

      ! Physics single column module extra levels initialization (as needs soil flag set in iniphyparam)
      call initialize_extra_levels (Nsoil + 1)
   end subroutine init_physics

   subroutine write_physics_params (file_unit,file_params)
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

      open(unit=file_unit, file=trim(file_params), form="FORMATTED", action='WRITE', status='REPLACE')

      ! Print physics parameters
      write (file_unit,*) "planet_rat  = ", radius
      write (file_unit,*) "g           = ", grav_accel
      write (file_unit,*) "cpp         = ", c_p
      write (file_unit,*) "mugaz       = ", gas_molarmass !28.9702532_8  !8314.46261815324/R_d
      write (file_unit,*) "unjours     = ", DAY
      write (file_unit,*) "year_day    = ", int (YEAR / DAY)
      write (file_unit,*) "periheli    = ", perihelion
      write (file_unit,*) "aphelie     = ", aphelion
      write (file_unit,*) "peri_day    = ", perihelion_day
      write (file_unit,*) "obliquit    = ", obliquity
      write (file_unit,*) "Cd_mer      = ", sea_surf
      write (file_unit,*) "Cd_ter      = ", soil_surf
      write (file_unit,*) "I_mer       = ", sea_interia
      write (file_unit,*) "I_ter       = ", soil_interia
      write (file_unit,*) "alb_ter     = ", sea_albedo
      write (file_unit,*) "alb_mer     = ", soil_albedo
      write (file_unit,*) "emi_mer     = ", sea_emissive
      write (file_unit,*) "emi_ter     = ", soil_emmisive
      write (file_unit,*) "emin_turb   = ", emin_turb
      write (file_unit,*) "lmixmin     = ", min_turbmix
      write (file_unit,*) "coefvis     = ", sw_atten
      write (file_unit,*) "coefir      = ", lw_atten
      write (file_unit,*) "callrad     = ", radiation_mod
      write (file_unit,*) "calldifv    = ", turbulence_mod
      write (file_unit,*) "calladj     = ", convecAdj_mod
      write (file_unit,*) "callsoil    = ", soil_mod
      write (file_unit,*) "season      = ", seasons
      write (file_unit,*) "diurnal     = ", diurnal
      write (file_unit,*) "lverbose    = ", physics_write
      write (file_unit,*) "period_sort = ", 1d0

      close(file_unit)
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
      ! Physics use cases
      use comgeomfi, only : init_comgeomfi, nsoilmx

      ! Local variables
      real(8) :: lat(1), long(1)

      if (.not. soil_mod) then
         print*, 'STOP!! Cannot use default (init_soil_grid_default) to set zmin when soil_mod flag indicates false'
         stop
      end if

      ! Dummy latatiude and longitude for initialization
      lat(1)  = 0
      long(1) = 0

      ! Call grid initialization for the physics
      call init_comgeomfi (1, zlevels, long, lat)

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
      ! Physics use cases
      use comgeomfi, only : init_comgeomfi

      ! Local variables
      real(8) :: lat(1), long(1)

      ! Dummy latatiude and longitude for initialization
      lat(1)  = 0
      long(1) = 0

      ! Set the zmin (lowest vertical level index) (! See Assumptions)
      if (soil_mod) then
         zmin = - Nsoil
         call init_comgeomfi (1, zlevels, long, lat, Nsoil) ! physics grid initialization 
      else
         Nsoil = 0
         zmin  = 0
         call init_comgeomfi (1, zlevels, long, lat)        ! physics grid initialization 
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
      ! Physics use cases
      use phyparam_mod, only : alloc, precompute, zday_last, icount

      ! Local variables
      real(8) :: day_fraction, nth_day

      ! Set flag for first call to physics false (thus the soil levels will get updated)
      physics_firstcall_flag = .false.

      ! Call allocation for physics call (usually done on the physics first call)
      call alloc (1, zlevels)

      ! Call precompute for physics call(usually done on the physics first call)
      call precompute ()

      ! Set the previous day in physics
      day_fraction = (time - dt) / DAY
      nth_day      = floor (day_fraction)
      day_fraction = day_fraction - nth_day
      zday_last    = nth_day + day_fraction - dt / DAY
      
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

      if (trim(param_name) /= trim(name)) val = defval
   end subroutine

   subroutine read_parami(name, defval, val, comment)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Physics plugin to read integers from a file
      !
      !   Assumption: File is already open with unit number: 9*rank
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      character(*), intent(in)  :: name, comment
      integer,      intent(in)  :: defval
      integer,      intent(out) :: val

      character(15) :: param_name
      character     :: equal_sign

      ! Read line from input file
      read (9*rank,*) param_name, equal_sign, val

      if (trim (param_name) /= trim (name)) val = defval
   end subroutine

   subroutine read_paramb (name, defval, val, comment)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Physics plugin to read LOGICALs from a file
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

      if (trim(param_name) /= trim(name)) val = defval
   end subroutine

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
      character(kind=c_char), intent(in)        :: tag(taglen), buf(buflen, bufsize)

      integer             :: i
      character(buflen+1) :: line
      character(100)      :: prefix

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
