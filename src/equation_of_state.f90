module equation_of_state_mod
  ! Equation of state for seawater based on simplified NEMO/TEOS10 simplified model that gives density in terms of temperature and salinity
  ! WAVETRISK-OCEAN uses density and temperature, and not potential density and potential temperature
  ! (provides a simple linear representation of thermal expansion and thermobaricity effects and nonlinear representation of cabbeling)
  ! (set lambda_1 = lambda_2 = 0 to remove cabbeling, set mu_1 = mu_2 = 0 to remove thermobaric effects)
  use shared_mod
  implicit none
contains
  real(8) function buoyancy_eos (salinity, temperature, z)
    implicit none
    real(8) :: temperature, salinity, z

    buoyancy_eos = 1 - density_eos (salinity, temperature, z) / ref_density
  end function buoyancy_eos

  real(8) function density_eos (salinity, temperature,  z)
    ! Equation of state: returns density as a function of temperature, salinity and depth z
    ! (includes cabbeling)
    implicit none
    real(8) :: salinity, temperature, z
    
    real(8) :: S_a, T_a

    S_a = salinity - S_ref
    T_a = temperature - T_ref

    density_eos = ref_density - a_0 * (1 + 0.5*lambda_1*T_a + mu_1*z) * T_a + b_0 * (1 - 0.5*lambda_2*S_a - mu_2*z) * S_a &
         - nu_0 * S_a * T_a
  end function density_eos

  real(8) function temperature_eos (density, salinity, z)
    ! Inverse equation of state: returns temperature from density, salinity Sal and depth z
    !!!! does not include (nonlinear) cabbeling effects: use only for visualization !!!!
    implicit none
    real(8) :: density, salinity, z

    real(8) :: rho_a, S_a

    rho_a = ref_density - density
    S_a   = salinity - S_ref

    temperature_eos = T_ref + (rho_a + b_0 * (1 - mu_2*z) * S_a) / (a_0 * (1 + mu_1*z))       
  end function temperature_eos

  real(8) function dk_buoyancy_eos (salinity, temperature, dk_salinity, dk_temperature, z)
    ! Vertical difference of buoyancy given salinity and temperature and vertical differences in salinity and temperature
    ! neglects any change in z
    ! (used in flux computations)
    implicit none
    real(8) :: dk_salinity, dk_temperature, temperature, salinity, z

    real(8) :: S_a, T_a

    S_a = salinity - S_ref
    T_a = temperature - T_ref

    dk_buoyancy_eos = (a_0 * (1 + 0.5*lambda_1*T_a + mu_1*z) * dk_temperature &
         - b_0 * (1 - 0.5*lambda_2*S_a - mu_2*z) * dk_salinity &
         + nu_0 * (S_a * dk_temperature + T_a * dk_salinity + dk_salinity * dk_temperature)) / ref_density
  end function dk_buoyancy_eos
end module equation_of_state_mod
