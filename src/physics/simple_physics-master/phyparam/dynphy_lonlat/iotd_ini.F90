      SUBROUTINE iotd_ini(fichnom,iim,jjm,llm,prlonv,prlatu,pcoordv)
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
!-----------------------------------------------------------------------
!   Declarations:
!   -------------

#include "netcdf.inc"
#include "iotd.h"

!   Arguments:
!   ----------

      integer iim,jjm,llm
      real prlonv(iim),prlatu(jjm),pcoordv(llm),timestep
      INTEGER id_FOCE

      integer corner(4),edges(4),ndim
      real  px(1000)
      character (len=10) :: nom

!   Local:
!   ------
      INTEGER ierr

      integer :: nvarid
      integer, dimension(2) :: id  
      integer :: varid

      character*10 fichnom
      real*4 rlonv(iim),rlatu(jjm),coordv(llm)

      real pi

      print*,'INIIO prlonv ',prlonv
      imax=iim
      jmax=jjm
      lmax=llm

      rlonv=prlonv
      rlatu=prlatu
      coordv=pcoordv

!-----------------------------------------------------------------------


      pi=2.*asin(1.)

! Define dimensions
    
         ! Create the NetCDF file
         ierr=NF_CREATE(fichnom, NF_CLOBBER, nid)
         ! Define the 'Time' dimension
         ierr=nf_def_dim(nid,"Time",NF_UNLIMITED,dim_coord(4))
         ! Define the 'Time' variable
         ierr=NF_DEF_VAR(nid, "Time", NF_FLOAT, 1, dim_coord(4),varid)
!        ! Add a long_name attribute
!        ierr=NF_PUT_ATT_TEXT(nid, varid, "long_name",4,"Time")
!        ! Add a units attribute
         ierr=NF_PUT_ATT_TEXT(nid, varid,'units',29,"days since 0000-00-0 00:00:00")
         ! Switch out of NetCDF Define mode

      ierr=NF_DEF_DIM(nid, "longitude", iim, dim_coord(1))
      ierr=NF_DEF_DIM(nid, "latitude", jjm, dim_coord(2))
      ierr=NF_DEF_DIM(nid, "altitude", llm, dim_coord(3))


      ierr=NF_ENDDEF(nid)
!
!  Contol parameters for this run
! --------------------------

      ierr=NF_REDEF(nid)
      ierr=NF_DEF_VAR(nid,"longitude", NF_FLOAT, 1, dim_coord(1),nvarid)
!     ierr=NF_PUT_ATT_TEXT(nid,nvarid,"long_name", 14,
!    .      "East longitude")
!     ierr=NF_PUT_ATT_TEXT(nid,nvarid,'units',12,"degrees_east")
      ierr=NF_ENDDEF(nid)
      ierr=NF_PUT_VAR_REAL(nid,nvarid,rlonv)
       print*,ierr

! --------------------------
      ierr=NF_REDEF(nid)
      ierr=NF_DEF_VAR(nid, "latitude", NF_FLOAT, 1, dim_coord(2),nvarid)
!     ierr=NF_PUT_ATT_TEXT(nid,nvarid,'units',13,"degrees_north")
!     ierr=NF_PUT_ATT_TEXT(nid,nvarid,"long_name", 14,"North latitude")
      ierr=NF_ENDDEF(nid)
      ierr=NF_PUT_VAR_REAL(nid,nvarid,rlatu)
!
! --------------------------
      ierr=NF_REDEF(nid)
      ierr=NF_DEF_VAR(nid, "altitude", NF_FLOAT, 1,dim_coord(3),nvarid)
      ierr=NF_PUT_ATT_TEXT(nid,nvarid,"long_name",10,"pseudo-alt")
!     ierr=NF_PUT_ATT_TEXT(nid,nvarid,'units',2,"km")
      if ( pcoordv(2)>pcoordv(1) ) then
         ierr=NF_PUT_ATT_TEXT(nid,nvarid,"long_name",10,"pseudo-alt")
         ierr=NF_PUT_ATT_TEXT(nid,nvarid,'positive',2,"up")
      else
         ierr=NF_PUT_ATT_TEXT(nid,nvarid,"long_name",8,"pressure")
         ierr = NF_PUT_ATT_TEXT (nid,nvarid,'positive',4,"down")
      endif
      ierr=NF_ENDDEF(nid)

      ierr=NF_PUT_VAR_REAL(nid,nvarid,coordv)
!
      END
