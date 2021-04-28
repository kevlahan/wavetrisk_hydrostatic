module equation_of_state_mod
  ! Equation of state for seawater based on simplified NEMO/TEOS10 simplified linear model (nemo_eos = .true.)
  ! (provides a simple linear representation of both cabbeling and thermobaricity effects)
  ! or a basic linear model (nemo_eos = .false.)
  use shared_mod
  implicit none
contains
  real(8) function buoyancy (salinity, temperature, z)
    implicit none
    real(8) :: temperature, salinity, z

    buoyancy = 1 - density (salinity, temperature, z) / ref_density
  end function buoyancy

  real(8) function density (salinity, temperature,  z)
    ! Equation of state: returns density as a function of temperature, salinity and depth z
    implicit none
    real(8) :: salinity, temperature, z

    density = ref_density - a_0 * (temperature - T_ref) * (1 + mu_1*z) + b_0 * (salinity - S_ref) * (1 - mu_2*z)       
  end function density

  real(8) function temperature (density, salinity, z)
    ! Equation of state: returns temperature from density, salinity Sal and depth z
    implicit none
    real(8) :: density, salinity, z

    temperature = T_ref + ((ref_density - density) - b_0 * (1 - mu_2*z) * (salinity - S_ref)) / (a_0 * (1 + mu_1*z))       
  end function temperature
end module equation_of_state_mod
