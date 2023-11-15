MODULE read_param_mod
#include "use_logging.h"
  IMPLICIT NONE
  PRIVATE
  SAVE

  INTERFACE

     ! each of these plugins reads a parameter of a certain type and stores the value in 'val'
     ! a typical implementation reads from a configuration file containing 'name = value' statements
     ! if 'name' is not present in the config file, defval is used as default value
     ! 'comment' can be used for logging purposes

     SUBROUTINE plugin_read_paramr(name, defval, val, comment)
       CHARACTER(*), INTENT(IN) :: name, comment
       REAL, INTENT(IN)         :: defval
       REAL, INTENT(OUT)        :: val
     END SUBROUTINE plugin_read_paramr

     SUBROUTINE plugin_read_parami(name, defval, val, comment)
       CHARACTER(*), INTENT(IN) :: name, comment
       INTEGER, INTENT(IN)      :: defval
       INTEGER, INTENT(OUT)     :: val
     END SUBROUTINE plugin_read_parami

     SUBROUTINE plugin_read_paramb(name, defval, val, comment)
       CHARACTER(*), INTENT(IN) :: name, comment
       LOGICAL, INTENT(IN)      :: defval
       LOGICAL, INTENT(OUT)     :: val
     END SUBROUTINE plugin_read_paramb

  END INTERFACE

  PROCEDURE(plugin_read_paramr), POINTER, PUBLIC :: read_paramr_plugin => NULL()
  PROCEDURE(plugin_read_parami), POINTER, PUBLIC :: read_parami_plugin => NULL()
  PROCEDURE(plugin_read_paramb), POINTER, PUBLIC :: read_paramb_plugin => NULL()

  INTERFACE read_param
     PROCEDURE read_paramr, read_parami, read_paramb
  END INTERFACE read_param

  PUBLIC :: read_param

CONTAINS

  SUBROUTINE read_paramr(name, defval, val, comment)
    CHARACTER(*), INTENT(IN) :: name, comment
    REAL, INTENT(IN)         :: defval
    REAL, INTENT(OUT)        :: val

    IF(.NOT.ASSOCIATED(read_paramr_plugin)) THEN
       CALL missing_plugin('read_paramr','read_param_mod')
       read_paramr_plugin => default_read_paramr
    END IF

    CALL read_paramr_plugin(name, defval, val, comment)

    WRITELOG(*,*) name, ' = ', val
    LOG_INFO('read_param')
  END SUBROUTINE read_paramr

  SUBROUTINE default_read_paramr(name, defval, val, comment)
    CHARACTER(*), INTENT(IN) :: name, comment
    REAL, INTENT(IN)         :: defval
    REAL, INTENT(OUT)        :: val
    val = defval
  END SUBROUTINE default_read_paramr

  SUBROUTINE read_parami(name, defval, val, comment)
    CHARACTER(*), INTENT(IN) :: name, comment
    INTEGER, INTENT(IN)      :: defval
    INTEGER, INTENT(OUT)     :: val

    IF(.NOT.ASSOCIATED(read_parami_plugin)) THEN
       CALL missing_plugin('read_parami','read_param_mod')
       val = defval
    ELSE
       CALL read_parami_plugin(name, defval, val, comment)
    END IF

    WRITELOG(*,*) name, ' = ', val
    LOG_INFO('read_param')
  END SUBROUTINE read_parami

  SUBROUTINE read_paramb(name, defval, val, comment)
    CHARACTER(*), INTENT(IN) :: name, comment
    LOGICAL, INTENT(IN)      :: defval
    LOGICAL, INTENT(OUT)     :: val

    IF(.NOT.ASSOCIATED(read_paramb_plugin)) THEN
       CALL missing_plugin('read_paramb','read_param_mod')
       read_paramb_plugin => default_read_paramb
    END IF

    CALL read_paramb_plugin(name, defval, val, comment)

    WRITELOG(*,*) name, ' = ', val
    LOG_INFO('read_param')
  END SUBROUTINE read_paramb

  SUBROUTINE default_read_paramb(name, defval, val, comment)
    CHARACTER(*), INTENT(IN) :: name, comment
    LOGICAL, INTENT(IN)      :: defval
    LOGICAL, INTENT(OUT)     :: val
    val = defval
  END SUBROUTINE default_read_paramb

END MODULE read_param_mod
