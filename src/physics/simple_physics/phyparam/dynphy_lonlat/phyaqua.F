! Routines complementaires pour la physique planetaire.


      subroutine iniaqua(nlon,latfi,lonfi,iflag_phys)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Creation d'un etat initial et de conditions aux limites 
!  (resp startphy.nc et limit.nc) pour des configurations idealisees 
! du modele LMDZ dans sa version terrestre.
!  iflag_phys est un parametre qui controle
!  iflag_phys = N  
!    de 100 a 199 : aqua planetes avec SST forcees
!                 N-100 determine le type de SSTs
!    de 200 a 299 : terra planetes avec Ts calcule
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      integer nlon,iflag_phys
cIM ajout latfi, lonfi
      REAL, DIMENSION (nlon) :: lonfi, latfi


      return
      end

