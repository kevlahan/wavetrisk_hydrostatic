module equation_of_state_mod
  ! Equation of state for seawater based on simplified NEMO/TEOS10 simplified model that gives density in terms of temperature and salinity
  ! WAVETRISK-OCEAN uses density and temperature, and not potential density and potential temperature
  ! (provides a simple linear representation of thermal expansion and thermobaricity effects and nonlinear representation of cabbeling)
  ! (set lambda_1 = lambda_2 = 0 to remove cabbeling, set mu_1 = mu_2 = 0 to remove thermobaric effects)

  use kind_mod,   only : dp
  use shared_mod, only : a_0, b_0, eos_nl, lambda_1, lambda_2, mu_1, mu_2, nu_0, S_ref, T_ref, ref_density

  implicit none


contains


  function buoyancy_eos (salinity, temperature, z) result(val)
    implicit none
    real(dp), intent(in) :: temperature, salinity, z
    real(dp)             :: val

    val = 1.0_dp - density_eos (salinity, temperature, z) / ref_density
  end function buoyancy_eos


  function density_eos (salinity, temperature,  z) result(val)
    ! Equation of state: returns density as a function of temperature, salinity and depth z
    ! (includes cabbeling)
    implicit none
    real(dp), intent(in) :: salinity, temperature, z
    real(dp)             :: val

    real(dp) :: S_a, T_a

    S_a = salinity - S_ref
    T_a = temperature - T_ref

    if (eos_nl) then  ! nonlinear equation of state
       val = ref_density &
            - a_0 * (1.0_dp + 0.5_dp * lambda_1 * T_a + mu_1 * z) * T_a &
            + b_0 * (1.0_dp - 0.5_dp * lambda_2 * S_a - mu_2 * z) * S_a &
            - nu_0 * S_a * T_a
    else
       val = ref_density - a_0 * (1.0_dp +  mu_1 * z) * T_a + b_0 * (1.0_dp - mu_2 * z) * S_a
    end if
  end function density_eos

  
  function temperature_eos (density, salinity, z) result(val)
    ! Inverse equation of state: returns temperature from density, salinity Sal and depth z
!!!! does not include (nonlinear) cabbeling effects: use only for visualization !!!!
    implicit none
    real(dp), intent(in) :: density, salinity, z
    real(dp)             :: val

    real(dp) :: rho_a, S_a

    rho_a = density - ref_density 
    S_a   = salinity - S_ref

    val = T_ref + (- rho_a + b_0 * (1.0_dp - mu_2 * z) * S_a) / (a_0 * (1.0_dp + mu_1 * z))       
  end function temperature_eos

  function dk_buoyancy_eos (salinity, temperature, dk_salinity, dk_temperature, z) result(val)
    ! Vertical difference of buoyancy given salinity and temperature and vertical differences in salinity and temperature
    ! neglects any change in z
    ! (used in flux computations)
    implicit none
    real(dp), intent(in) :: dk_salinity, dk_temperature, temperature, salinity, z
    real(dp)             :: val

    real(dp) :: S_a, T_a

    S_a = salinity - S_ref
    T_a = temperature - T_ref

    if (eos_nl) then ! nonlinear equation of state
       val = &
            (a_0 * (1.0_dp + 0.5_dp * lambda_1 * T_a + mu_1 * z) * dk_temperature &
            - b_0 * (1.0_dp - 0.5_dp * lambda_2 * S_a - mu_2 * z) * dk_salinity &
            + nu_0 * (S_a * dk_temperature + T_a * dk_salinity + dk_salinity * dk_temperature)) / ref_density
    else
       val = (a_0 * (1.0_dp +  mu_1 * z) * dk_temperature  - b_0 * (1.0_dp - mu_2 * z) * dk_salinity) / ref_density
    end if
  end function dk_buoyancy_eos

  
end module equation_of_state_mod
