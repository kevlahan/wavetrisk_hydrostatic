MODULE writefield_mod
  USE logging, ONLY : missing_plugin
  IMPLICIT NONE
  PRIVATE
  SAVE

  INTERFACE
     SUBROUTINE plugin_writefield1(name,longname,unit, var)
       CHARACTER(*), INTENT(IN) :: name, longname, unit
       REAL, INTENT(IN)         :: var(:)
     END SUBROUTINE plugin_writefield1
     SUBROUTINE plugin_writefield2(name,longname,unit, var)
       CHARACTER(*), INTENT(IN) :: name, longname, unit
       REAL, INTENT(IN)         :: var(:,:)
     END SUBROUTINE plugin_writefield2
  END INTERFACE

  PROCEDURE(plugin_writefield1), POINTER, PUBLIC :: writefield1_plugin => NULL()
  PROCEDURE(plugin_writefield2), POINTER, PUBLIC :: writefield2_plugin => NULL()

  INTERFACE writefield
     PROCEDURE writefield1, writefield2
  END INTERFACE writefield

  PUBLIC :: writefield

CONTAINS

  SUBROUTINE writefield2(name, longname, unit, var)
    CHARACTER(*), INTENT(IN) :: name, longname, unit
    REAL, INTENT(IN)         :: var(:,:)

    IF(ASSOCIATED(writefield2_plugin)) THEN
       CALL writefield2_plugin(name, longname, unit, var)
    ELSE
       CALL missing_plugin('writefield2','writefield_mod')
    END IF

  END SUBROUTINE writefield2

  SUBROUTINE writefield1(name, longname, unit, var)
    CHARACTER(*), INTENT(IN) :: name, longname, unit
    REAL, INTENT(IN)         :: var(:)

    IF(ASSOCIATED(writefield1_plugin)) THEN
       CALL writefield1_plugin(name, longname, unit, var)
    ELSE
       CALL missing_plugin('writefield1','writefield_mod')
    END IF

  END SUBROUTINE writefield1

END MODULE writefield_mod
