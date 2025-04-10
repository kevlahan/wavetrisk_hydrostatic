module time_integr_mod
  use adapt_mod
  use barotropic_2d_mod
  implicit none
  type(Float_Field), dimension(:,:), allocatable :: q1, q2, q3, q4, dq1
  
  interface
     subroutine trend_fun (q, dq)
       use domain_mod
       implicit none
       type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq
     end subroutine trend_fun
  end interface
  
  abstract interface
     subroutine dt_integrator (q, wav, trend_fun, h)
       use domain_mod
       implicit none
       real(8)                                              :: h      ! time step
       type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q, wav ! solution and wavelet coefficients
       
       interface
          subroutine trend_fun (q, dq)
            use domain_mod
            implicit none
            type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q, dq
          end subroutine trend_fun
       end interface
     end subroutine dt_integrator
  end interface
  procedure (dt_integrator), pointer :: dt_step  => null ()

  abstract interface
     subroutine dt_integrator_split (h)
       ! Euler time step
       ! Stable for CFL<1, first order
       implicit none
       real(8) :: h ! time step
     end subroutine dt_integrator_split
  end interface
  procedure (dt_integrator_split), pointer :: dt_step_split  => null ()
contains
  subroutine Euler (q, wav, trend_fun, h)
    ! Euler time step
    ! Stable for CFL<1, first order
    implicit none
    real(8)                                              :: h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q, wav

    call trend_fun (q, trend)
    call RK_sub_step (q, trend, h, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine Euler

  subroutine RK3 (q, wav, trend_fun, h)
    ! Optimal third order, three stage strong stability preserving Runge-Kutta method
    ! Stable for hyperbolic equations for CFL<2
    ! Does not require extra solution variables.
    implicit none
    real(8)                                              :: h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q, wav

    call manage_q1_mem

    call trend_fun (q, trend) 
    call RK_sub_step (q, trend, h/3d0, q1)
    call WT_after_step (q1, wav)

    call trend_fun (q1, trend) 
    call RK_sub_step (q, trend, h/2d0, q1)
    call WT_after_step (q1, wav)

    call trend_fun (q1, trend) 
    call RK_sub_step (q, trend, h, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK3
  
  subroutine RK4 (q, wav, trend_fun, h)
    ! Low storage four stage second order accurate Runge-Kutta scheme used in Dubos et al (2015) Geosci. Model Dev., 8, 3131–3150, 2015.
    ! Fourth order accurate for linear equations, stable for CFL <= 2*sqrt(2) ~ 2.83.
    ! Does not require extra solution variables.
    implicit none
    real(8)                                              :: h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q, wav

    call manage_q1_mem

    call trend_fun (q, trend) 
    call RK_sub_step (q, trend, h/4d0, q1)
    call WT_after_step (q1, wav)

    call trend_fun (q1, trend) 
    call RK_sub_step (q, trend, h/3d0, q1)
    call WT_after_step (q1, wav)

    call trend_fun (q1, trend) 
    call RK_sub_step (q, trend, h/2d0, q1)
    call WT_after_step (q1, wav)

    call trend_fun (q1, trend) 
    call RK_sub_step (q, trend, h, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK4

  subroutine RK33_opt (q, wav, trend_fun, h)
    ! Optimal third order, three stage strong stability preserving Runge-Kutta method
    ! Stable for hyperbolic equations for CFL<2
    ! Spiteri and Ruuth (SIAM J. Numer. Anal., 40(2): 469-491, 2002) Appendix A.1
    implicit none
    real(8)                                               :: h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)  :: q, wav

    call manage_RK_mem

    call trend_fun (q, trend) 
    call RK_sub_step (q, trend, h, q1)
    call WT_after_step (q1, wav)

    call trend_fun (q1, trend) 
    call RK_sub_step2 (q, q1, trend, (/ 0.75d0, 0.25d0 /), h/4d0, q2)
    call WT_after_step (q2, wav)

    call trend_fun (q2, trend)
    call RK_sub_step2 (q, q2, trend, (/ 1d0/3d0, 2d0/3d0 /), h*2d0/3d0, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK33_opt

  subroutine RK34_opt (q, wav, trend_fun, h)
    ! Optimal third order, four stage strong stability preserving Runge-Kutta method
    ! Stable for hyperbolic equations for CFL<2
    ! Spiteri and Ruuth (SIAM J. Numer. Anal., 40(2): 469-491, 2002) Appendix A.1
    implicit none
    real(8)                                              :: h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q, wav

    call manage_RK_mem

    call trend_fun (q, trend) 
    call RK_sub_step (q, trend, h/2d0, q1)
    call WT_after_step (q1, wav)

    call trend_fun (q1, trend) 
    call RK_sub_step (q1, trend, h/2d0, q2)
    call WT_after_step (q2, wav)

    call trend_fun (q2, trend)
    call RK_sub_step2 (q, q2, trend, (/ 2d0/3d0, 1d0/3d0 /), h/6d0, q3)
    call WT_after_step (q3, wav)

    call trend_fun (q3, trend) 
    call RK_sub_step (q3, trend, h/2d0, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK34_opt

  subroutine RK45_opt (q, wav, trend_fun, h)
    ! Optimal fourth order, five stage strong stability preserving Runge-Kutta method stable with optimal maximum CFL coefficient of 1.51
    ! See A. Balan, G. May and J. Schoberl: "A Stable Spectral Difference Method for Triangles", 2011, Spiteri and Ruuth 2002
    implicit none
    real(8)                                              :: h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q, wav

    real(8), dimension(5,5) :: alpha, beta

    alpha = reshape ((/1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.444370494067_8, 0.555629505933_8, 0.0_8, 0.0_8, 0.0_8,  &
         0.620101851385_8, 0.0_8, 0.379898148615_8, 0.0_8, 0.0_8, 0.178079954108_8, 0.0_8, 0.0_8, 0.821920045892_8, 0.0_8,  &
         0.006833258840_8, 0.0_8, 0.517231672090_8, 0.127598311333_8, 0.348336757737_8/), (/5, 5/))

    beta = reshape ((/0.391752227004_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.36841059263_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
         0.251891774247_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.544974750212_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0846041633821_8, &
         0.226007483194_8/), (/5, 5/))

    call manage_RK_mem

    call trend_fun (q, trend) 
    call RK_sub_step1 (q, trend, alpha(1,1), h * beta(1,1), q1)
    call WT_after_step (q1, wav)

    call trend_fun (q1, trend)
    call RK_sub_step2 (q, q1, trend, alpha(1:2,2), h * beta(2,2), q2)
    call WT_after_step (q2, wav)

    call trend_fun (q2, trend)
    call RK_sub_step2 (q, q2, trend, (/alpha(1,3), alpha(3,3)/), h * beta(3,3), q3)
    call WT_after_step (q3, wav)

    call trend_fun (q3, trend)
    call RK_sub_step2 (q, q3, trend, (/alpha(1,4), alpha(4,4)/), h * beta(4,4), q4)
    call WT_after_step (q4, wav)

    call trend_fun (q4, dq1)
    call RK_sub_step4 (q, q2, q3, q4, trend, dq1, (/alpha(1,5), alpha(3:5,5)/), h * beta(4:5,5), q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK45_opt

  subroutine init_time_integr_mod
    logical :: initialized = .false.

    if (initialized) return ! initialize only once

    call init_comm_mod
    call init_ops_mod
    call init_multi_level_mod
    initialized = .true.
  end subroutine init_time_integr_mod

  subroutine RK_sub_step (sols, trends, h, dest)
    implicit none
    real(8)                                                             :: h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)                :: sols
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)                :: trends
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), intent(inout) :: dest
    
    integer :: d, ibeg, iend, k, v

    do v = 1, N_VARIABLE
       do d = 1, size(grid)
          ibeg = (1+2*(POSIT(v)-1)) * grid(d)%patch%elts(2+1)%elts_start + 1
          iend = dest(v,1)%data(d)%length
          do k = 1, zlevels
             dest(v,k)%data(d)%elts(ibeg:iend) = sols(v,k)%data(d)%elts(ibeg:iend) + h * trends(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do
    end do
    dest%bdry_uptodate = .false.
  end subroutine RK_sub_step

  subroutine RK_sub_step1 (sols, trends, alpha, h, dest)
    implicit none
    real(8)                                                             :: alpha, h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)                :: sols
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)                :: trends
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), intent(inout) :: dest

    integer :: k, v, d, ibeg, iend

    do v = 1, N_VARIABLE
       do d = 1, size(grid)
       ibeg = (1+2*(POSIT(v)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
       iend = sols(v,1)%data(d)%length
       do k = 1, zlevels
             dest(v,k)%data(d)%elts(ibeg:iend) = alpha * sols(v,k)%data(d)%elts(ibeg:iend) &
                  + h * trends(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do
       dest%bdry_uptodate = .False.
    end do
  end subroutine RK_sub_step1

  subroutine RK_sub_step2 (sol1, sol2, trends, alpha, h, dest)
    implicit none
    real(8)                                                              :: h
    real(8), dimension(2)                                                :: alpha
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)                :: sol1, sol2
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)                :: trends
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), intent(inout) :: dest
    
    integer :: k, v, d, ibeg, iend

    do v = 1, N_VARIABLE
       do d = 1, size(grid)
          ibeg = (1+2*(POSIT(v)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
          iend = dest(v,1)%data(d)%length
          do k = 1, zlevels
             dest(v,k)%data(d)%elts(ibeg:iend) = alpha(1) * sol1(v,k)%data(d)%elts(ibeg:iend) &
                                               + alpha(2) * sol2(v,k)%data(d)%elts(ibeg:iend) &
                  + h * trends(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do
       dest%bdry_uptodate = .false.
    end do
  end subroutine RK_sub_step2

  subroutine RK_sub_step4 (sol1, sol2, sol3, sol4, trend1, trend2, alpha, h, dest)
    implicit none
    real(8), dimension(2)                                               :: h
    real(8), dimension(4)                                               :: alpha
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)                :: sol1, sol2, sol3, sol4
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels)                :: trend1, trend2
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), intent(inout) :: dest

    integer :: k, v, d, ibeg, iend
    
    do v = 1, N_VARIABLE
       do d = 1, size(grid)
          ibeg = (1+2*(POSIT(v)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
          iend = dest(v,1)%data(d)%length
          do k = 1, zlevels
             dest(v,k)%data(d)%elts(ibeg:iend) = &
                    alpha(1) * sol1(v,k)%data(d)%elts(ibeg:iend) + alpha(2) * sol2(v,k)%data(d)%elts(ibeg:iend) &
                  + alpha(3) * sol3(v,k)%data(d)%elts(ibeg:iend) + alpha(4) * sol4(v,k)%data(d)%elts(ibeg:iend) &
                  + h(1) * trend1(v,k)%data(d)%elts(ibeg:iend) + h(2) * trend2(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do
       dest%bdry_uptodate = .False.
    end do
  end subroutine RK_sub_step4

  subroutine init_RK_mem
    implicit none
    integer :: d, k, v

    allocate (q1(1:N_VARIABLE,1:zmax), q2(1:N_VARIABLE,1:zmax), q3(1:N_VARIABLE,1:zmax), &
         q4(1:N_VARIABLE,1:zmax), dq1(1:N_VARIABLE,1:zmax))

    do k = 1, zmax
       do v = 1, N_VARIABLE
          call init_Float_Field (q1(v,k),  POSIT(v))
          call init_Float_Field (q2(v,k),  POSIT(v))
          call init_Float_Field (q3(v,k),  POSIT(v))
          call init_Float_Field (q4(v,k),  POSIT(v))
          call init_Float_Field (dq1(v,k), POSIT(v))
       end do

       do d = 1, size(grid)
          do v = 1, N_VARIABLE
             call init (q1(v,k)%data(d),  sol(v,k)%data(d)%length);  q1(v,k)%data(d)%elts = dble (N_VARIABLE-v)
             call init (q2(v,k)%data(d),  sol(v,k)%data(d)%length);  q2(v,k)%data(d)%elts = dble (N_VARIABLE-v)
             call init (q3(v,k)%data(d),  sol(v,k)%data(d)%length);  q3(v,k)%data(d)%elts = dble (N_VARIABLE-v)
             call init (q4(v,k)%data(d),  sol(v,k)%data(d)%length);  q4(v,k)%data(d)%elts = dble (N_VARIABLE-v)
             call init (dq1(v,k)%data(d), sol(v,k)%data(d)%length); dq1(v,k)%data(d)%elts = dble (N_VARIABLE-v) 
          end do
       end do
    end do
  end subroutine init_RK_mem

  subroutine manage_q1_mem
    implicit none
    integer :: d, k, v, n_new

    do k = 1, zmax
       do d = 1, size(grid)
          do v = 1, N_VARIABLE
             n_new = sol(v,k)%data(d)%length - q1(v,k)%data(d)%length
             if (n_new > 0) call extend (q1(v,k)%data(d), n_new, dble (N_VARIABLE-v))
          end do
       end do
    end do
  end subroutine manage_q1_mem

  subroutine manage_RK_mem
    implicit none
    integer :: d, k, v, n_new

    do k = 1, zmax
       do d = 1, size(grid)
          do v = 1, N_VARIABLE
             n_new = sol(v,k)%data(d)%length - q1(v,k)%data(d)%length
             if (n_new > 0) then
                call extend ( q1(v,k)%data(d), n_new, dble (N_VARIABLE-v))
                call extend ( q2(v,k)%data(d), n_new, dble (N_VARIABLE-v))
                call extend ( q3(v,k)%data(d), n_new, dble (N_VARIABLE-v))
                call extend ( q4(v,k)%data(d), n_new, dble (N_VARIABLE-v))
                call extend (dq1(v,k)%data(d), n_new, dble (N_VARIABLE-v))
             end if
          end do
       end do
    end do
  end subroutine manage_RK_mem

  subroutine RK4_split (h)
    ! Low storage four stage Runge-Kutta scheme used in Dubos et al (2015) Geosci. Model Dev., 8, 3131–3150, 2015.
    ! Fourth order accurate for linear equations, second order accurate for nonlinear equations.
    ! Stable for CFL <= 2*sqrt(2) ~ 2.83.
    ! Does not require extra solution variables.
    !
    ! This version implements the explicit-implicit free surface method used in the MITgcm.
    implicit none
    real(8)  :: h
    
    call manage_q1_mem

    call update_bdry (sol(:,1:zlevels+1), NONE, 968)

    ! Compute flux divergence of vertically integrated velocity at previous time step
    if (theta2 /= 1d0) call flux_divergence (sol, trend(S_TEMP,zlevels+1))

    call RK_split (h/4d0, sol, q1)
    call RK_split (h/3d0, q1,  q1)
    call RK_split (h/2d0, q1,  q1)
    call RK_split (h,     q1, sol)
    call free_surface_update 
  end subroutine RK4_split
  
  subroutine RK3_split (h)
    ! Low storage three stage Runge-Kutta from Kinnmark and Gray (Math Computers Simul 26 1984, 181-188)
    ! Third order accurate for linear equations, second order accurate for nonlinear equations.
    ! Stable for CFL <= sqrt(3) ~ 1.7321.
    ! Does not require extra solution variables.
    !
    ! This version implements the explicit-implicit free surface method used in the MITgcm.
    implicit none
    real(8)  :: h
    
    call manage_q1_mem

    call update_bdry (sol(:,1:zlevels+1), NONE, 969)

    ! Compute flux divergence of vertically integrated velocity at previous time step
    if (theta2 /= 1d0) call flux_divergence (sol, trend(S_TEMP,zlevels+1))

    call RK_split (h/3d0, sol,  q1)
    call RK_split (h/2d0,  q1,  q1)
    call RK_split (h,      q1, sol)
    call free_surface_update 
  end subroutine RK3_split

  subroutine RK2_split (h)
    ! Low storage two stage Runge-Kutta from Kinnmark and Gray (Math Computers Simul 26 1984, 181-188)
    ! Second order accurate.
    ! Stable for CFL <= sqrt(3) ~ 1.7321.
    ! Does not require extra solution variables.
    !
    ! This version implements the explicit-implicit free surface method used in the MITgcm.
    implicit none
    real(8)  :: h
    
    call manage_q1_mem

    call update_bdry (sol(:,1:zlevels+1), NONE, 970)

    ! Compute flux divergence of vertically integrated velocity at previous time step
    if (theta2 /= 1d0) call flux_divergence (sol, trend(S_TEMP,zlevels+1))

    call RK_split (h/2d0, sol, q1)
    call RK_split (h,     q1, sol)
    call free_surface_update 
  end subroutine RK2_split

  subroutine Euler_split (h)
    ! Euler time step for barotropic mode splitting
    ! Stable for CFL<1, first order
    implicit none        
    real(8) :: h

    call update_bdry (sol(:,1:zlevels+1), NONE, 971)

    ! Compute flux divergence of vertically integrated velocity at previous time step
    if (theta2 /= 1d0) call flux_divergence (sol, trend(S_TEMP,zlevels+1))

    call RK_split (h, sol, sol)
    call free_surface_update
  end subroutine Euler_split

  subroutine RK_split (h, sol1, sol2)
    ! Explicit Euler integration of velocity and scalars used in RK4_split
    implicit none
    real(8)                                                :: h
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels+1) :: sol1
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels+1) :: sol2

    ! Compute explicit trends
    call barotropic_correction (sol1(1:N_VARIABLE,1:zlevels+1))
    call trend_ml (sol1(1:N_VARIABLE,1:zlevels), trend)
    
    ! Explicit Euler step for scalars
    call scalar_star (h, sol2(1:N_VARIABLE,1:zlevels))
    
    ! Explicit Euler step for intermediate 3D baroclinic velocities u_star
    call u_star (h, sol2(1:N_VARIABLE,1:zlevels))
        
    ! Inverse wavelet transform of solution onto adaptive grid
    call WT_after_step (sol2(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels))
  end subroutine RK_split

  subroutine free_surface_update
    ! Backwards Euler implicit calculation of new free surface and correction of velocity and scalars
    implicit none

    ! Backwards Euler step for new free surface, updates sol(S_MASS,zlevels+1)
    call eta_update
    call barotropic_correction (sol(1:N_VARIABLE,1:zlevels+1))
    
    ! Explicit Euler step to update 3D baroclinic velocities with new external pressure gradient
    call u_update

    ! Inverse wavelet transform of solution onto adaptive grid
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine free_surface_update
end module time_integr_mod
