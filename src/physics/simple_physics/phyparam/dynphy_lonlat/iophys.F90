MODULE iophys
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: iophys_ini
  
  CONTAINS

      subroutine iophys_ecrit(nom,lllm,titre,unite,px)
      USE dimphy
      USE mod_phys_lmdz_para
      USE mod_grid_phy_lmdz
      IMPLICIT NONE



!  Ecriture de variables diagnostiques au choix dans la physique 
!  dans un fichier NetCDF nomme  'diagfi'. Ces variables peuvent etre
!  3d (ex : temperature), 2d (ex : temperature de surface), ou
!  0d (pour un scalaire qui ne depend que du temps : ex : la longitude
!  solaire)
!  (ou encore 1d, dans le cas de testphys1d, pour sortir une colonne)
!  La periode d'ecriture est donnee par 
!  "ecritphy " regle dans le fichier de controle de run :  run.def
!
!    writediagfi peut etre appele de n'importe quelle subroutine
!    de la physique, plusieurs fois. L'initialisation et la creation du
!    fichier se fait au tout premier appel.
!
! WARNING : les variables dynamique (u,v,t,q,ps)
!  sauvees par writediagfi avec une
! date donnee sont legerement differentes que dans le fichier histoire car 
! on ne leur a pas encore ajoute de la dissipation et de la physique !!!
! IL est  RECOMMANDE d'ajouter les tendance physique a ces variables
! avant l'ecriture dans diagfi (cf. physiq.F)
!  
! Modifs: Aug.2010 Ehouarn: enforce outputs to be real*4
!
!  parametres (input) :
!  ----------
!      unit : unite logique du fichier de sortie (toujours la meme)
!      nom  : nom de la variable a sortir (chaine de caracteres)
!      titre: titre de la variable (chaine de caracteres)
!      unite : unite de la variable (chaine de caracteres)
!      px : variable a sortir (real 0, 1, 2, ou 3d)
!
!=================================================================

#include "dimensions.h"
#include "paramet.h"
#include "netcdf.inc"
#include "iotd.h"


! Arguments on input:
      integer lllm
      character (len=*) :: nom,titre,unite
      integer imjmax
      parameter (imjmax=100000)
      real px(klon_omp,lllm)
      real xglo(klon_glo,lllm)
      real zx(iim,jjp1,lllm)


      CALL Gather(px,xglo)
!$OMP MASTER
      IF (is_mpi_root) THEN       
        CALL Grid1Dto2D_glo(xglo,zx)
        call iotd_ecrit(nom,lllm,titre,unite,zx)
      ENDIF
!$OMP END MASTER

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE writefield1(name, longname, unit, var)
        CHARACTER(*), INTENT(IN) :: name, longname, unit
        REAL, INTENT(IN)         :: var(:)
        CALL iophys_ecrit(name, 1, longname, unit, var)
      END SUBROUTINE writefield1

      SUBROUTINE writefield2(name, longname, unit, var)
        CHARACTER(*), INTENT(IN) :: name, longname, unit
        REAL, INTENT(IN)         :: var(:,:)
        PRINT *, 'writefield2', name, SHAPE(var)
        CALL iophys_ecrit(name, SIZE(var,2), longname, unit, var)
      END SUBROUTINE writefield2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE iophys_ini(fichnom,presnivs)
      USE mod_phys_lmdz_para
      USE writefield_mod
      IMPLICIT NONE

!=======================================================================
!
!   Auteur:  L. Fairhead  ,  P. Le Van, Y. Wanherdrick, F. Forget
!   -------
!
!   Objet:
!   ------
!
!   'Initialize' the diagfi.nc file: write down dimensions as well
!   as time-independent fields (e.g: geopotential, mesh area, ...)
!
!=======================================================================
!-----------------------------------------------------------------------
!   Declarations:
!   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
! #include "comvert.h"
REAL presnivs(llm)

character*10 fichnom
real pi

!   Arguments:
!   ----------

writefield1_plugin => writefield1
writefield2_plugin => writefield2

!$OMP MASTER
    IF (is_mpi_root) THEN       
pi=2.*asin(1.)
call iotd_ini('phys.nc     ', &
iim,jjp1,llm,rlonv(1:iim)*180./pi,rlatu*180./pi,presnivs)
    ENDIF
!$OMP END MASTER

      END

END MODULE
