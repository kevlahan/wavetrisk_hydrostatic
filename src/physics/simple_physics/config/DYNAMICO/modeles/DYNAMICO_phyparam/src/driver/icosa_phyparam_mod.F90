MODULE icosa_phyparam_mod
#include "use_logging.h"

  ! FCM gets confused when external modules are USEd at module level
  ! => USE statements to DYNAMICO modules go into subroutines

  USE, INTRINSIC :: ISO_C_BINDING
  USE icosa, ONLY : t_field

  IMPLICIT NONE
  PRIVATE
  SAVE

  LOGICAL(KIND=C_BOOL)            :: firstcall = .TRUE.
  LOGICAL(KIND=C_BOOL), PARAMETER :: lastcall = .FALSE.

  REAL, PARAMETER :: oneday = 86400. ! hard-coded
  INTEGER :: log_unit

  TYPE(t_field),POINTER :: f_write2d(:), f_write_llm(:), f_write_llmp1(:)

  PUBLIC :: init_physics, physics

CONTAINS

  SUBROUTINE init_physics
    ! DYNAMICO
    USE mpipara,  ONLY : is_mpi_master, mpi_rank
    USE icosa,    ONLY : llm

    USE icosa,     ONLY : g, radius, cpp, kappa
    USE getin_mod, ONLY : getin
    USE profiling_mod, ONLY : register_id, enter_profile, exit_profile
    USE physics_interface_mod, ONLY : inout => physics_inout
    ! phyparam
    USE logging, ONLY   : flush_plugin, dbtag, max_log_level
    USE profiling, ONLY : profile_register_plugin, profile_enter_plugin, profile_exit_plugin
    USE read_param_mod
    USE comgeomfi
    USE iniphyparam_mod
    INTEGER, PARAMETER :: dayref=0
    CHARACTER(10) :: physics_log_level
    INTEGER :: ngrid, lev
    REAL    :: timestep
    REAL    :: unjours, & ! solar day in seconds
         single_lon, single_lat ! longitude and lattitude if single column
    LOGICAL :: single_column

    log_unit = 10 + mpi_rank
    PRINT *, 'init_physics : MPI Rank => log_unit', mpi_rank, log_unit

    flush_plugin => flush_log_
    profile_register_plugin => register_id
    profile_enter_plugin    => enter_profile
    profile_exit_plugin     => exit_profile
    CALL init_plugin_writefield

    physics_log_level='INF'
    CALL getin('phyparam_log_level', physics_log_level)
    DO lev=1, SIZE(dbtag)
       IF(dbtag(lev)==TRIM(physics_log_level)) THEN
          max_log_level = lev
          EXIT
       END IF
    END DO

    read_paramr_plugin => read_paramr
    read_parami_plugin => read_parami
    read_paramb_plugin => read_paramb

    WRITELOG(*,*) 'init_physics called'
    WRITELOG(*,*) 'physics log level set to ', dbtag(max_log_level)
    LOG_INFO('phyparam')

    ngrid = inout%ngrid
    timestep = inout%dt_phys

    unjours = 86400.
    CALL getin('unjours', unjours)

    single_column=.FALSE.
    CALL getin('phyparam_single_column', single_column)
    IF(single_column) THEN
       single_lon=0.
       CALL getin('phyparam_single_lon', single_lon)
       single_lat=0.
       CALL getin('phyparam_single_lat', single_lat)
       inout%lon(:)=single_lon
       inout%lat(:)=single_lat
    END IF

    CALL init_comgeomfi(ngrid, llm, inout%lon, inout%lat)
    CALL iniphyparam(timestep, unjours, radius, g, cpp*kappa, cpp)

    WRITELOG(*,*) 'init_physics done'
    LOG_INFO('phyparam')

  END SUBROUTINE init_physics

  SUBROUTINE physics
    USE mpipara,  ONLY : is_mpi_master
    USE icosa,    ONLY : llm
    USE physics_interface_mod, ONLY : inout => physics_inout
    USE phyparam_mod
    USE error_mod
    REAL :: dps(inout%ngrid), play(inout%ngrid, llm), pphi(inout%ngrid, llm)
    REAL :: timestep, time, jourvrai, gmtime
    INTEGER :: l
    IF(is_mpi_master) WRITE(log_unit,*) 'phyparam/physics called', SHAPE(inout%p), SHAPE(inout%pk)

    timestep = inout%dt_phys
    time = timestep * inout%it
    gmtime = time/oneday
    jourvrai = FLOOR(gmtime)
    gmtime   = gmtime - jourvrai

    WRITELOG(*,*) 'it, timestep, time, jourvrai, gmtime', inout%it, timestep, time, jourvrai, gmtime
    LOG_DBG('physics')

    ! compute pressure and geopotential at full levels
    CALL compute_play(inout%ngrid, llm, inout%p, play)
    CALL compute_play(inout%ngrid, llm, inout%geopot, pphi)

    ! substract surface geopotential
    DO l=1,llm
       pphi(:,l) = pphi(:,l) - inout%geopot(:,1)
    END DO

    !      IF(is_mpi_master) PRINT *, 'phyparam phi :', pphi(inout%ngrid/2+1, :)

    CALL check_NaN('physics', 'ulon', inout%ulon)
    CALL check_NaN('physics', 'ulat', inout%ulat)
    CALL check_NaN('physics', 'temp', inout%temp)

    ! go
    CALL phyparam(inout%ngrid,llm,                       &
         &        firstcall,lastcall,                    &
         &        jourvrai, gmtime, timestep,            &
         &        inout%p, play, pphi,                   &
         &        inout%ulon,  inout%ulat,  inout%temp,  &
         &        inout%dulon, inout%dulat, inout%dtemp, dps)

    !      IF(is_mpi_master) PRINT *, 'phyparam dT :', inout%dtemp(inout%ngrid/2+1, :)

    CALL check_NaN('physics', 'dulon', inout%dulon)
    CALL check_NaN('physics', 'dulat', inout%dulat)
    CALL check_NaN('physics', 'dtemp', inout%dtemp)

    firstcall = .FALSE.
  END SUBROUTINE physics

  SUBROUTINE compute_play(ngrid, llm, plev, play)
    INTEGER, INTENT(IN) :: ngrid, llm
    REAL, INTENT(IN)    :: plev(ngrid, llm+1) ! pressure at interfaces (half-levels)
    REAL, INTENT(OUT)   :: play(ngrid, llm)   ! pressure in layers (full levels)
    INTEGER :: ij, l
    DO l = 1,llm
       DO ij = 1,ngrid
          play(ij,l) = .5*(plev(ij,l)+plev(ij,l+1))
       END DO
    END DO
  END SUBROUTINE compute_play

  !------------------------------------------------------------------------------------
  !------------------------------- Infrastructure plugins -----------------------------

  !--------------------------------------- Logging ------------------------------------

  SUBROUTINE flush_log_(lev, taglen, tag, buflen, bufsize, buf) BIND(C)
    USE mpipara, ONLY : is_mpi_master
    USE logging, ONLY : dbtag
    USE, INTRINSIC :: iso_c_binding, ONLY : c_char, c_null_char, c_int
    INTEGER(c_int), INTENT(IN), VALUE :: lev, taglen, buflen, bufsize
    CHARACTER(KIND=c_char), INTENT(IN) :: tag(taglen), buf(buflen, bufsize)

    CHARACTER(buflen+1) :: line
    CHARACTER(100) :: prefix
    LOGICAL :: verbose
    INTEGER :: i

    verbose = is_mpi_master .OR. (lev<log_level_info)

    IF(verbose) THEN
       WRITE(prefix,*) '[', dbtag(lev), ' ', tag, ']'
       DO i=1, bufsize
          WRITE(line,*) buf(:,i)
          WRITE(log_unit,*) TRIM(prefix) // TRIM(line)
       END DO
       WRITE(log_unit, *) ''
       FLUSH(log_unit)
    END IF
  END SUBROUTINE flush_log_

  !--------------------------------------- read_param ------------------------------------

  SUBROUTINE read_paramr(name, defval, val, comment)
    USE getin_mod, ONLY : getin
    CHARACTER(*), INTENT(IN) :: name, comment
    REAL, INTENT(IN)         :: defval
    REAL, INTENT(OUT)        :: val
    val = defval
    CALL getin('phyparam_'//name, val)
  END SUBROUTINE read_paramr

  SUBROUTINE read_parami(name, defval, val, comment)
    USE getin_mod, ONLY : getin
    CHARACTER(*), INTENT(IN) :: name, comment
    INTEGER, INTENT(IN)      :: defval
    INTEGER, INTENT(OUT)     :: val
    val = defval
    CALL getin('phyparam_'//name, val)
  END SUBROUTINE read_parami

  SUBROUTINE read_paramb(name, defval, val, comment)
    USE getin_mod, ONLY : getin
    CHARACTER(*), INTENT(IN) :: name, comment
    LOGICAL, INTENT(IN)      :: defval
    LOGICAL, INTENT(OUT)     :: val
    val = defval
    CALL getin('phyparam_'//name, val)
  END SUBROUTINE read_paramb

  !--------------------------------------- writefield ------------------------------------

  SUBROUTINE init_plugin_writefield
    USE icosa, ONLY : t_field, field_t, type_real, allocate_field, llm
    USE writefield_mod, ONLY : writefield1_plugin, writefield2_plugin
    CALL allocate_field(f_write2d,     field_t, type_real,        name='phyparam_write2d')
    CALL allocate_field(f_write_llm,   field_t, type_real, llm,   name='phyparam_write_llm')
    CALL allocate_field(f_write_llmp1, field_t, type_real, llm+1, name='phyparam_write_llmp1')
    writefield1_plugin => plugin_writefield1
    writefield2_plugin => plugin_writefield2
  END SUBROUTINE init_plugin_writefield

  SUBROUTINE plugin_writefield1(name,longname,unit, var)
    USE physics_interface_mod, ONLY : unpack_field, inout => physics_inout
    USE output_field_mod, ONLY : output_field, xios_output
    CHARACTER(*), INTENT(IN) :: name, longname, unit
    REAL, INTENT(IN)         :: var(:)
    IF(xios_output) THEN
       WRITELOG(*,*) TRIM(name), ' : ', TRIM(longname), SHAPE(var), inout%it
       WRITELOG(*,*) TRIM(name), ' : ', MINVAL(var), MAXVAL(var)
       LOG_INFO('writefield1')
       CALL unpack_field(f_write2d, var)
       CALL output_field('phyparam_'//TRIM(name), f_write2d)
    END IF
  END SUBROUTINE plugin_writefield1

  SUBROUTINE plugin_writefield2(name,longname,unit, var)
    USE physics_interface_mod, ONLY : unpack_field, inout => physics_inout
    USE output_field_mod, ONLY : output_field, xios_output
    USE icosa,    ONLY : llm
    CHARACTER(*), INTENT(IN) :: name, longname, unit
    REAL, INTENT(IN)         :: var(:,:)
    INTEGER :: nlev
    IF(xios_output) THEN
       WRITELOG(*,*) TRIM(name), ' : ', TRIM(longname), SHAPE(var), inout%it
       WRITELOG(*,*) TRIM(name), ' : ', MINVAL(var), MAXVAL(var)
       LOG_INFO('writefield2')
       nlev = SIZE(var, 2)
       IF(nlev==llm) THEN
          CALL unpack_field(f_write_llm, var)
          CALL output_field('phyparam_'//TRIM(name), f_write_llm)
       ELSEIF(nlev==llm+1) THEN
          CALL unpack_field(f_write_llmp1, var)
          CALL output_field('phyparam_'//TRIM(name), f_write_llmp1)
       END IF
    END IF
  END SUBROUTINE plugin_writefield2

END MODULE icosa_phyparam_mod
