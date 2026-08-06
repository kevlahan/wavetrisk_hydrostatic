module time_integr_mod

  use kind_mod,   only : dp
  use shared_mod, only : N_VARIABLE, NONE, POSIT, S_TEMP, eps, level_start, theta2, zlevels, zmax
  
  use adapt_mod,         only : init_multi_level_mod, WT_after_step
  use barotropic_2d_mod, only : barotropic_correction, eta_update, flux_divergence, scalar_star, u_star, u_update
  use comm_mod,          only : init_comm_mod
  use comm_mpi_mod,      only : update_bdry
  use dyn_arrays,        only : extend, init
  use domain_mod,        only : Float_Field, init_Field, grid, sol, trend, wav_coeff
  use multi_level_mod,   only : trend_ml
  use ops_mod,           only : init_ops_mod

  implicit none

  private
  public :: dt_step, dt_step_split
  public :: init_RK_mem, init_time_integr_mod
  public :: Euler, Euler_split, RK2_split, RK3, RK3_split, RK33_opt, RK34_opt, RK4, RK45_opt, RK4_split
  public :: q1, q2, q3, q4, dq1
  
  type(Float_Field), dimension(:,:), allocatable :: q1, q2, q3, q4, dq1
  
  interface
     subroutine trend_sub (q, dq)
       use domain_mod, only : Float_Field
       use shared_mod, only : N_VARIABLE, zlevels
       implicit none
       type(Float_Field), intent(inout), target ::  q(1:N_VARIABLE,1:zlevels)
       type(Float_Field), intent(inout), target :: dq(1:N_VARIABLE,1:zlevels)
     end subroutine trend_sub
  end interface
  
  abstract interface
     
     subroutine dt_integrator (q, wav, routine, h)
       use kind_mod,   only : dp
       use domain_mod, only : Float_Field
       use shared_mod, only : N_VARIABLE, zlevels
       implicit none
       real(dp),          intent(in)    :: h      
       type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
       type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
       procedure (trend_sub)            :: routine
     end subroutine dt_integrator

     subroutine dt_integrator_split (h)
       use kind_mod, only : dp
       implicit none
       real(dp), intent(in) :: h 
     end subroutine dt_integrator_split
     
  end interface
  
  procedure (dt_integrator),       pointer :: dt_step        => null ()
  procedure (dt_integrator_split), pointer :: dt_step_split  => null ()
  
contains

  
  subroutine Euler (q, wav, routine, h)
    ! Euler time step
    ! Stable for CFL<1, first order
    
    implicit none
    
    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
    procedure (trend_sub)            :: routine

    call routine (q, trend)
    call RK_sub_step (q, trend, h, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine Euler
  

  subroutine RK3 (q, wav, routine, h)
    ! Optimal third order, three stage strong stability preserving Runge-Kutta method
    ! Stable for hyperbolic equations for CFL<2
    ! Does not require extra solution variables.
    
    implicit none

    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
    procedure (trend_sub)            :: routine

    call manage_q1_mem

    call routine (q, trend) 
    call RK_sub_step (q, trend, h/3, q1)
    call WT_after_step (q1, wav)

    call routine (q1, trend) 
    call RK_sub_step (q, trend, h/2, q1)
    call WT_after_step (q1, wav)

    call routine (q1, trend) 
    call RK_sub_step (q, trend, h, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK3
  
  
  subroutine RK4 (q, wav, routine, h)
    ! Low storage four stage second order accurate Runge-Kutta scheme used in Dubos et al (2015) Geosci. Model Dev., 8, 3131–3150, 2015.
    ! Fourth order accurate for linear equations, stable for CFL <= 2*sqrt(2) ~ 2.83.
    ! Does not require extra solution variables.
    
    implicit none

    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
    procedure (trend_sub)            :: routine

    call manage_q1_mem

    call routine (q, trend) 
    call RK_sub_step (q, trend, h/4, q1)
    call WT_after_step (q1, wav)

    call routine (q1, trend) 
    call RK_sub_step (q, trend, h/3, q1)
    call WT_after_step (q1, wav)

    call routine (q1, trend) 
    call RK_sub_step (q, trend, h/2, q1)
    call WT_after_step (q1, wav)

    call routine (q1, trend) 
    call RK_sub_step (q, trend, h, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK4
  

  subroutine RK33_opt (q, wav, routine, h)
    ! Optimal third order, three stage strong stability preserving Runge-Kutta method
    ! Stable for hyperbolic equations for CFL<2
    ! Spiteri and Ruuth (SIAM J. Numer. Anal., 40(2): 469-491, 2002) Appendix A.1
    
    implicit none

    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
    procedure (trend_sub)            :: routine

    call manage_RK_mem

    call routine (q, trend) 
    call RK_sub_step (q, trend, h, q1)
    call WT_after_step (q1, wav)

    call routine (q1, trend) 
    call RK_sub_step2 (q, q1, trend, [ 0.75_dp, 0.25_dp ], h/4, q2)
    call WT_after_step (q2, wav)

    call routine (q2, trend)
    call RK_sub_step2 (q, q2, trend, [ 1.0_dp/3, 2.0_dp/3 ], h * 2.0_dp/3, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK33_opt
  

  subroutine RK34_opt (q, wav, routine, h)
    ! Optimal third order, four stage strong stability preserving Runge-Kutta method
    ! Stable for hyperbolic equations for CFL<2.65
    ! Spiteri and Ruuth (SIAM J. Numer. Anal., 40(2): 469-491, 2002) Appendix A.1
    
    implicit none

    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
    procedure (trend_sub)            :: routine

    call manage_RK_mem

    call routine (q, trend) 
    call RK_sub_step (q, trend, h/2, q1)
    call WT_after_step (q1, wav)

    call routine (q1, trend) 
    call RK_sub_step (q1, trend, h/2, q2)
    call WT_after_step (q2, wav)

    call routine (q2, trend)
    call RK_sub_step2 (q, q2, trend, [ 2.0_dp/3, 1.0_dp/3 ], h/6, q3)
    call WT_after_step (q3, wav)

    call routine (q3, trend) 
    call RK_sub_step (q3, trend, h/2, q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK34_opt
  

  subroutine RK45_opt (q, wav, routine, h)
    ! Optimal fourth order, five stage strong stability preserving Runge-Kutta method stable with optimal maximum CFL coefficient of 2
    ! Spiteri and Ruuth (SIAM J. Numer. Anal., 40(2): 469-491, 2002) Appendix A.1
    
    implicit none

    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
    procedure (trend_sub)            :: routine

    real(dp), dimension(5,5) :: alpha, beta

    alpha = reshape ([ &
         1.0_dp,            0.0_dp,            0.0_dp,            0.0_dp,            0.0_dp, &
         0.444370494067_dp, 0.555629505933_dp, 0.0_dp,            0.0_dp,            0.0_dp, &
         0.620101851385_dp, 0.0_dp,            0.379898148615_dp, 0.0_dp,            0.0_dp, &
         0.178079954108_dp, 0.0_dp,            0.0_dp,            0.821920045892_dp, 0.0_dp, &
         0.006833258840_dp, 0.0_dp,            0.517231672090_dp, 0.127598311333_dp, 0.348336757737_dp ], [5, 5])

    beta = reshape ([&
         0.391752227004_dp, 0.0_dp,           0.0_dp,            0.0_dp,             0.0_dp, &
         0.0_dp,            0.36841059263_dp, 0.0_dp,            0.0_dp,             0.0_dp, &
         0.0_dp,            0.0_dp,           0.251891774247_dp, 0.0_dp,             0.0_dp, &
         0.0_dp,            0.0_dp,           0.0_dp,            0.544974750212_dp,  0.0_dp, &
         0.0_dp,            0.0_dp,           0.0_dp,            0.0846041633821_dp, 0.226007483194_dp ], [5, 5])

    call manage_RK_mem

    call routine (q, trend) 
    call RK_sub_step1 (q, trend, alpha(1,1), h * beta(1,1), q1)
    call WT_after_step (q1, wav)

    call routine (q1, trend)
    call RK_sub_step2 (q, q1, trend, alpha(1:2,2), h * beta(2,2), q2)
    call WT_after_step (q2, wav)

    call routine (q2, trend)
    call RK_sub_step2 (q, q2, trend, [alpha(1,3), alpha(3,3)], h * beta(3,3), q3)
    call WT_after_step (q3, wav)

    call routine (q3, trend)
    call RK_sub_step2 (q, q3, trend, [alpha(1,4), alpha(4,4)], h * beta(4,4), q4)
    call WT_after_step (q4, wav)

    call routine (q4, dq1)
    call RK_sub_step4 (q, q2, q3, q4, trend, dq1, [alpha(1,5), alpha(3:5,5)], h * beta(4:5,5), q)
    call WT_after_step (q, wav, level_start-1)
  end subroutine RK45_opt
  

  subroutine init_time_integr_mod
    implicit none
    
    logical :: initialized = .false.

    if (initialized) return ! initialize only once

    call init_comm_mod
    call init_ops_mod
    call init_multi_level_mod
    initialized = .true.
  end subroutine init_time_integr_mod
  

  subroutine RK_sub_step (sols, trends, h, dest)
    
    implicit none
    
    real(dp),          intent(in)    :: h
    type(Float_Field), intent(in)    :: sols(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: trends(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: dest(1:N_VARIABLE,1:zlevels)
    
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

    real(dp),          intent(in)    :: alpha, h
    type(Float_Field), intent(in)    :: sols(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: trends(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: dest(1:N_VARIABLE,1:zlevels)

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
    
    real(dp),          intent(in)    :: h
    real(dp),          intent(in)    :: alpha(2)
    type(Float_Field), intent(in)    :: sol1(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: sol2(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: trends(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: dest(1:N_VARIABLE,1:zlevels)
    
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
    
    real(dp),          intent(in)    :: h(2)
    real(dp),          intent(in)    :: alpha(4)
    type(Float_Field), intent(in)    :: sol1(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: sol2(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: sol3(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: sol4(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: trend1(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: trend2(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: dest(1:N_VARIABLE,1:zlevels)

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
       dest%bdry_uptodate = .false.
    end do
  end subroutine RK_sub_step4
  

  subroutine init_RK_mem
    
    implicit none
    
    integer :: d, k, v

    allocate (q1(1:N_VARIABLE,1:zmax), q2(1:N_VARIABLE,1:zmax), q3(1:N_VARIABLE,1:zmax), &
         q4(1:N_VARIABLE,1:zmax), dq1(1:N_VARIABLE,1:zmax))

    do k = 1, zmax
       do v = 1, N_VARIABLE
          call init_Field (q1(v,k),  POSIT(v))
          call init_Field (q2(v,k),  POSIT(v))
          call init_Field (q3(v,k),  POSIT(v))
          call init_Field (q4(v,k),  POSIT(v))
          call init_Field (dq1(v,k), POSIT(v))
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
    
    real(dp), intent(in)  :: h
    
    call manage_q1_mem

    call update_bdry (sol(:,1:zlevels+1), NONE, 968)

    ! Compute flux divergence of vertically integrated velocity at previous time step
    if (abs (theta2 - 1.0_dp) > eps(1.0_dp)) call flux_divergence (sol, trend(S_TEMP,zlevels+1))

    call RK_split (h/4, sol, q1)
    call RK_split (h/3, q1,  q1)
    call RK_split (h/2, q1,  q1)
    call RK_split (h,   q1, sol)
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
    
    real(dp), intent(in)  :: h
    
    call manage_q1_mem

    call update_bdry (sol(:,1:zlevels+1), NONE, 969)

    ! Compute flux divergence of vertically integrated velocity at previous time step
    if (abs (theta2 - 1.0_dp) > eps(1.0_dp)) call flux_divergence (sol, trend(S_TEMP,zlevels+1))

    call RK_split (h/3, sol,q1)
    call RK_split (h/2, q1, q1)
    call RK_split (h,   q1, sol)
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
    
    real(dp), intent(in)  :: h
    
    call manage_q1_mem

    call update_bdry (sol(:,1:zlevels+1), NONE, 970)

    ! Compute flux divergence of vertically integrated velocity at previous time step
    if (abs (theta2 - 1.0_dp) > eps(1.0_dp)) call flux_divergence (sol, trend(S_TEMP,zlevels+1))

    call RK_split (h/2, sol, q1)
    call RK_split (h,   q1, sol)
    call free_surface_update 
  end subroutine RK2_split
  

  subroutine Euler_split (h)
    ! Euler time step for barotropic mode splitting
    ! Stable for CFL<1, first order
    
    implicit none
    
    real(dp), intent(in) :: h

    call update_bdry (sol(:,1:zlevels+1), NONE, 971)

    ! Compute flux divergence of vertically integrated velocity at previous time step
    if (abs (theta2 - 1.0_dp) > eps(1.0_dp)) call flux_divergence (sol, trend(S_TEMP,zlevels+1))

    call RK_split (h, sol, sol)
    call free_surface_update
  end subroutine Euler_split
  

  subroutine RK_split (h, sol1, sol2)
    ! Explicit Euler integration of velocity and scalars used in RK4_split
    
    implicit none
    
    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: sol1(1:N_VARIABLE,1:zlevels+1)
    type(Float_Field), intent(inout) :: sol2(1:N_VARIABLE,1:zlevels+1)

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
