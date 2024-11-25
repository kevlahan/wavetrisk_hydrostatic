module iniphyparam_mod
#include "use_logging.h"
  implicit none
  private
  real, parameter :: perfect_gas_const = 8314.46261815324 ! nb using g instead of kg for mass
  public          :: iniphyparam
contains
  subroutine read_params () bind(c, name='phyparam_setup')
    USE read_param_mod
    USE phys_const, only : planet_rad, g, r, Cpp, rCp, dtPhys, Emin_turb, Unjours, mugaz
    USE astronomy
    USE planet,     only : CoefIR, CoefVis
    USE turbulence, only : LmixMin
    USE soil_mod
    USE callkeys

    call read_param ('planet_rad', 6.4e6, planet_rad,'planet_rad')
    call read_param ('g',           9.8,  g,          'g')
    call read_param ('cpp',      1004.0,  Cpp,        'Cpp')
    call read_param ('mugaz',      28.0,  mugaz,      'mugaz')

    r = perfect_gas_const / mugaz
    rCp = r / Cpp

    call read_param ('unjours', 86400.0,  unjours,   'Unjours')
    call read_param ('year_day',  360.0,  year_day,  'year_day')
    call read_param ('periheli',  150.0,  periheli,  'periheli')
    call read_param ('aphelie',   150.0,  aphelie,   'aphelie')
    call read_param ('peri_day',  0.0,    peri_day,  'peri_day')
    call read_param ('obliquit',  23.0,   obliquit,  'obliquit')

    call read_param ('cd_mer',    0.01,   Cd_mer,    'Cd_mer')
    call read_param ('cd_ter',    0.01,   Cd_ter,    'Cd_ter')
    call read_param ('i_mer',      3e3,   i_mer,     'i_mer')
    call read_param ('i_ter',      3e3,   i_ter,     'i_ter')
    call read_param ('alb_ter',  0.112,   alb_ter,   'alb_ter')
    call read_param ('alb_mer',  0.112,   alb_mer,   'alb_mer')
    call read_param ('emi_mer',   1.0,    emi_mer,   'emi_mer')
    call read_param ('emi_ter',   1.0,    emi_ter,   'emi_ter')
    call read_param ('emin_turb', 1e-16,  Emin_turb, 'Emin_turb')
    call read_param ('lmixmin',    1e2,   lmixmin,   'LmixMin')

    call read_param ('coefvis',   0.99,   CoefVis,   'CoefVis')
    call read_param ('coefir',    0.08,   CoefIR,    'CoefIR')

    call read_param ('callrad',  .true.,  callrad,    'with rayonnement')
    call read_param ('calldifv', .true.,  calldifv,   'with vertical turbulent diffusion')
    call read_param ('calladj',  .true.,  calladj,    'with adj')
    call read_param ('callsoil', .true.,  callsoil,   'with soil')
    call read_param ('diurnal',  .true.,  diurnal,    'with diurnal cycle')
    call read_param ('lverbose', .true.,  lverbose,   'lverbose')
    call read_param ('period_sort', 1.0, period_sort, 'period sorties en jour')
  end subroutine read_params

  subroutine iniphyparam (pTimestep, pUnjours, pRad, pg, pR, pCpp)
    USE profiling,  only : profile_register, id_phyparam
    use comgeomfi,  only : nsoilmx
    use soil_mod,   only : init_soil
    use phys_const, only : planet_rad, g, r, Cpp, rCp, dtPhys, Unjours
    use callkeys
    real, intent(in)  :: pTimestep, pUnjours, pRad, pg, pr, pCpp

    call profile_register ('phyparam', id_phyparam)

    call read_params

    ! Frequency of computation of radiation
    ! if (diurnal) then
    !    iradia = nint (Unjours / (20.0 * pTimestep))
    ! else
    !    iradia = nint (Unjours / ( 4.0 * pTimestep))
    ! end if

    ! Compute radiative transfer and physics every time step
    dtPhys = pTimestep
    iradia = 1

    call check_mismatch ('day length (s)',         pUnjours,    Unjours)
    call check_mismatch ('planetary radius (km)',  pRad/1000.0, planet_rad/1000.0)
    call check_mismatch ('gravity',                pg,          g)
    call check_mismatch ('specific r',             pR,          r)
    call check_mismatch ('specific heat capacity', pCpp,        Cpp)
    LOG_WARN ('iniphyparam')

    WRITELOG (*,'(a)')        'Activation de la physique:'
    WRITELOG (*,'(a,l)')      ' rayonnement                    ', callrad
    WRITELOG (*,'(a,l)')      ' diffusion verticale turbulente ', calldifv
    WRITELOG (*,'(a,l)')      ' ajustement convectif           ', calladj
    WRITELOG (*,'(a,l)')      ' sol                            ', callsoil
    WRITELOG (*,'(a,l)')      ' cycle diurne                   ', diurnal
    WRITELOG (*,'(a,es10.4)') ' r                              ', r
    WRITELOG (*,'(a,es10.4)') ' unjours                        ', Unjours
    WRITELOG (*,'(a,i3,a,es10.4,a)') 'Radiative transfer is computed each ', iradia,' time step, or each ', iradia * pTimestep,' s'

    LOG_INFO ('iniphyparam')

    if (callsoil) call init_soil (nsoilmx)
  end subroutine iniphyparam

  subroutine check_mismatch (name, a, b)
    character(*), intent(in) :: name
    real,         intent(in) :: a, b

    if (a /= b) then
       WRITELOG (*,*) 'phys/dyn mismatch for ', name, ' : ', a, b
    end if
  end subroutine check_mismatch
end module iniphyparam_mod
