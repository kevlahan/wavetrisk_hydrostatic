MODULE restart
  IMPLICIT NONE
  PRIVATE
  SAVE

  ! right after allocating arrays for persistent data (internal state), the physics
  ! code calls restart_register for each array
  ! it is then left to the driver code to either call cold_start or to read data from disk into these arrays
  ! it is also left to the driver code to write internal state to disk at the end of the time loop

  INTEGER, PARAMETER :: max_restart_field=100, max_len=100
  INTEGER :: nb_restart_field=0

  TYPE restart_record_t
     CHARACTER(max_len) :: name
     REAL, POINTER :: data1(:), data2(:,:)
  END TYPE restart_record_t

  TYPE(restart_record_t) :: restart_record(max_restart_field)

  INTERFACE restart_register
     MODULE PROCEDURE restart_register_1D, restart_register_2D
  END INTERFACE restart_register

  PUBLIC :: restart_register

CONTAINS

  SUBROUTINE check_name(name)
    CHARACTER(*), INTENT(IN) :: name ! unique name of the restart field
    INTEGER :: id
    DO id=1, nb_restart_field
       IF(TRIM(restart_record(id)%name) == name) THEN
          PRINT *, 'FATAL : Restart field ', name, ' already registered'
          STOP
       END IF
    END DO
    nb_restart_field = nb_restart_field+1
  END SUBROUTINE check_name

  SUBROUTINE restart_register_1D(name, field)
    CHARACTER(*), INTENT(IN) :: name ! unique name of the restart field
    REAL, INTENT(IN), TARGET :: field(:) ! 1D field to be read/written at checkpoint/restart
    CALL check_name(name)
    restart_record(nb_restart_field)%name = name
    restart_record(nb_restart_field)%data1 => field
  END SUBROUTINE restart_register_1D

  SUBROUTINE restart_register_2D(name, field)
    CHARACTER(*), INTENT(IN) :: name ! unique name of the restart field
    REAL, INTENT(IN), TARGET :: field(:,:) ! 1D field to be read/written at checkpoint/restart
    CALL check_name(name)
    restart_record(nb_restart_field)%name = name
    restart_record(nb_restart_field)%data2 => field
  END SUBROUTINE restart_register_2D

END MODULE restart
