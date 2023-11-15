MODULE error_mod
#include "use_logging.h"

#ifndef XCODEML
  USE, INTRINSIC :: IEEE_ARITHMETIC
#endif

  IMPLICIT NONE
  PRIVATE

  INTERFACE check_NaN
     MODULE PROCEDURE check_NaN1, check_NaN2
  END INTERFACE check_NaN

  PUBLIC :: check_NaN

CONTAINS


  SUBROUTINE check_NaN1(caller, name, data)
    CHARACTER(*), INTENT(IN) ::  caller, name
    REAL, INTENT(IN) :: data(:)
    LOGICAL :: isnan(SIZE(data,1))
    INTEGER :: i
#ifndef XCODEML
    isnan = IEEE_IS_NAN(data)
#else
    STOP
#endif
    IF(ANY(isnan)) THEN
       WRITELOG(*,*) 'In subroutine ', caller, ' array ', name, ' has NaN . Offending indices :'
       DO i=1, SIZE(isnan,1)
          IF(isnan(i)) THEN
             WRITELOG(*,*) i, data(i)
          END IF
       END DO
       LOG_DBG('check_NaN')
    END IF
  END SUBROUTINE check_NaN1

  SUBROUTINE check_NaN2(caller, name, data)
    CHARACTER(*), INTENT(IN) ::  caller, name
    REAL, INTENT(IN) :: data(:,:)
    LOGICAL :: isnan(SIZE(data,1), SIZE(data,2))
    INTEGER :: i,j
#ifndef XCODEML
    isnan = IEEE_IS_NAN(data)
#else
    STOP
#endif
    IF(ANY(isnan)) THEN
       WRITELOG(*,*) 'In subroutine ', caller, ' array ', name, ' has NaN . Offending indices :'
       DO i=1, SIZE(isnan,1)
          DO j=1, SIZE(isnan,2)
             IF(isnan(i,j)) THEN
                WRITELOG(*,*) i,j, data(i,j)
             END IF
          END DO
       END DO
       LOG_DBG('check_NaN')
    END IF
  END SUBROUTINE check_NaN2

END MODULE error_mod
