MODULE phyparam_plugins_lmdz
  USE getparam
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: read_paramr, read_parami, read_paramb

#include "iniprint.h"

CONTAINS
  
  SUBROUTINE read_paramr(name, defval, val, comment)
    CHARACTER(*), INTENT(IN) :: name, comment
    REAL, INTENT(IN)         :: defval
    REAL, INTENT(OUT)        :: val
    CALL getpar(name, defval, val, comment)
    WRITE(lunout, *) TRIM(name),'=',val
  END SUBROUTINE read_paramr
  
  SUBROUTINE read_parami(name, defval, val, comment)
    CHARACTER(*), INTENT(IN) :: name, comment
    INTEGER, INTENT(IN)      :: defval
    INTEGER, INTENT(OUT)     :: val
    CALL getpar(name, defval, val, comment)
    WRITE(lunout, *) TRIM(name),'=',val
  END SUBROUTINE read_parami

  SUBROUTINE read_paramb(name, defval, val, comment)
    CHARACTER(*), INTENT(IN) :: name, comment
    LOGICAL, INTENT(IN)      :: defval
    LOGICAL, INTENT(OUT)     :: val
    CALL getpar(name, defval, val, comment)
    WRITE(lunout, *) TRIM(name),'=',val
  END SUBROUTINE read_paramb  
  
END MODULE phyparam_plugins_lmdz
