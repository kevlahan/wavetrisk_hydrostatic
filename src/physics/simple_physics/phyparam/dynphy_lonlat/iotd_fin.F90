      SUBROUTINE iotd_fin
      IMPLICIT NONE

!=======================================================================
!
!   Auteur:  F. Hourdin
!   -------
!
!   Objet:
!   ------
!   Light interface for netcdf outputs. can be used outside LMDZ
!
!=======================================================================


#include "netcdf.inc"
#include "iotd.h"
      integer ierr

!   Arguments:
!   ----------

      ierr=NF_close(nid)

      END
