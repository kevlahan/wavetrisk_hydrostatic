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

      integer imax,jmax,lmax,nid
      INTEGER dim_coord(4)
      real iotd_ts

      common/ecritd_c/imax,jmax,lmax,nid,dim_coord,iotd_ts
