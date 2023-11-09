MODULE iniphyparam_mod
#include "use_logging.h"
  IMPLICIT NONE
  PRIVATE

  REAL, PARAMETER :: perfect_gas_const = 8314.46261815324 ! NB using g instead of kg for mass

  PUBLIC :: iniphyparam

CONTAINS

  SUBROUTINE read_params() BIND(C, name='phyparam_setup')
    !$cython header void phyparam_setup();
    !$cython wrapper def setup() : phy.phyparam_setup()
    USE read_param_mod

    USE phys_const, ONLY : planet_rad,g,r,cpp,rcp,dtphys,unjours,mugaz
    USE astronomy
    USE planet, ONLY : coefir, coefvis
    USE turbulence, ONLY : lmixmin, emin_turb
    USE soil_mod
    USE callkeys

    CALL read_param('planet_rad',6.4e6 ,planet_rad,'planet_rad')
    CALL read_param('g',9.8            ,g,'g')
    CALL read_param('cpp',1004.        ,cpp,'cpp')
    CALL read_param('mugaz',28.        ,mugaz,'mugaz')
    r=perfect_gas_const/mugaz
    rcp=r/cpp

    CALL read_param('unjours', 86400.,  unjours,'unjours')
    CALL read_param('year_day',360.    ,year_day,'year_day')
    CALL read_param('periheli',150.    ,periheli,'periheli')
    CALL read_param('aphelie',150.     ,aphelie,'aphelie')
    CALL read_param('peri_day',0.      ,peri_day,'peri_day')
    CALL read_param('obliquit',23.     ,obliquit,'obliquit')

    CALL read_param('Cd_mer',.01       ,Cd_mer,'Cd_mer')
    CALL read_param('Cd_ter',.01       ,Cd_ter,'Cd_ter')
    CALL read_param('I_mer',3000.      ,I_mer,'I_mer')
    CALL read_param('I_ter',3000.      ,I_ter,'I_ter')
    CALL read_param('alb_ter',.112     ,alb_ter,'alb_ter')
    CALL read_param('alb_mer',.112     ,alb_mer,'alb_mer')
    CALL read_param('emi_mer',1.       ,emi_mer,'emi_mer')
    CALL read_param('emi_ter',1.       ,emi_ter,'emi_ter')
    CALL read_param('emin_turb',1.e-16 ,emin_turb,'emin_turb')
    CALL read_param('lmixmin',100.     ,lmixmin,'lmixmin')

    CALL read_param('coefvis',.99      ,coefvis,'coefvis')
    CALL read_param('coefir',.08       ,coefir,'coefir')

    CALL read_param('callrad',  .true.,  callrad,   'appel rayonnement')
    CALL read_param('calldifv', .true.,  calldifv,  'appel difv')
    CALL read_param('calladj',  .true.,  calladj,   'appel adj')
    CALL read_param('callsoil', .true.,  callsoil,  'appel soil')
    CALL read_param('season',   .true.,  season,    'with seasonal cycle')
    CALL read_param('diurnal',  .true., diurnal,   'with diurnal cycle')
    CALL read_param('lverbose', .true.,  lverbose,  'lverbose')
    CALL read_param('period_sort', 1., period_sort, 'period sorties en jour')

  END SUBROUTINE read_params

  SUBROUTINE iniphyparam(ptimestep, punjours, prad, pg, pr, pcpp)
    USE profiling, ONLY : profile_register, id_phyparam
    USE comgeomfi, ONLY : nsoilmx
    USE soil_mod, ONLY : init_soil
    USE phys_const, ONLY : planet_rad,g,r,cpp,rcp,dtphys,unjours
    USE callkeys
    REAL, INTENT(IN)  :: ptimestep, punjours, prad, pg, pr, pcpp

    CALL profile_register('phyparam', id_phyparam)

    CALL read_params
    !   choice of the frequency of the computation of radiations
    IF(diurnal) THEN
       iradia=NINT(unjours/(20.*ptimestep))
    ELSE
       iradia=NINT(unjours/(4.*ptimestep))
    ENDIF
    iradia=1
    dtphys=ptimestep

    CALL check_mismatch('day length (s)', punjours, unjours)
    CALL check_mismatch('planetary radius (km)', prad/1000., planet_rad/1000.)
    CALL check_mismatch('gravity', pg, g)
    CALL check_mismatch('specific R', pr, r)
    CALL check_mismatch('specific heat capacity', pcpp, cpp)
    LOG_WARN('iniphyparam')

    WRITELOG(*,*) 'Activation de la physique:'
    WRITELOG(*,*) ' R=',r
    WRITELOG(*,*) ' Rayonnement ',callrad
    WRITELOG(*,*) ' Diffusion verticale turbulente ', calldifv
    WRITELOG(*,*) ' Ajustement convectif ',calladj
    WRITELOG(*,*) ' Sol ',callsoil
    WRITELOG(*,*) ' Cycle diurne ',diurnal

    WRITELOG(*,*) 'unjours',unjours
    WRITELOG(*,*) 'The radiative transfer is computed each ', &
         &   iradia,' physical time-step or each ', &
         &   iradia*ptimestep,' seconds'

    LOG_INFO('iniphyparam')

    CALL init_soil(nsoilmx)
  END SUBROUTINE iniphyparam

  SUBROUTINE check_mismatch(name, a,b)
    CHARACTER(*), INTENT(IN) :: name
    REAL, INTENT(IN) :: a,b
    IF(a /= b) THEN
       WRITELOG(*,*) 'Phys/dyn mismatch for ', name, ' : ',a,b
    END IF
  END SUBROUTINE check_mismatch

END MODULE iniphyparam_mod
