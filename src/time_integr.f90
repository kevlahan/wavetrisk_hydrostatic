module time_integr_mod

  use kind_mod,   only : dp
  use shared_mod, only : N_VARIABLE, NONE, POSIT, S_TEMP, eps, level_start, theta2, zlevels, zmax
  
  use adapt_mod,         only : WT_after_step
  use barotropic_2d_mod, only : barotropic_correction, eta_update, flux_divergence, scalar_star, u_star, u_update
  use comm_mpi_mod,      only : update_bdry
  use dyn_arrays,        only : extend, init
  use domain_mod,        only : Float_Field, init_Field, grid, sol, trend, wav_coeff
  use multi_level_mod,   only : trend_ml
  use parallel_block_mpi_mod, only : &
       begin_block_domain_multistage_candidate_stage, &
       begin_block_scalar_divergence_capture, &
       capture_block_domain_multistage_candidate_tendency, &
       finalize_block_scalar_divergence_capture, &
       accept_block_domain_multistage_candidate, &
       parallel_block_state_is_ready, &
       refresh_parallel_block_candidate_boundary_state, &
       refresh_parallel_block_domain_prognostic_state

  implicit none

  private
  public :: dt_step, dt_step_split
  public :: init_RK_mem
  public :: set_multistage_block_candidate_enabled
  public :: Euler, Euler_split, RK3, RK3_split, RK4, RK4_split
  public :: q1
  
  type(Float_Field), allocatable :: q1(:,:)
  
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
  logical :: multistage_block_candidate_enabled = .false.

  
contains


  subroutine set_multistage_block_candidate_enabled (enabled)
    ! Guard the production multistage candidate so other RK3 and RK4 callers
    ! remain on the unchanged Domain-only pathway.

    implicit none

    logical, intent(in) :: enabled

    multistage_block_candidate_enabled = enabled

  end subroutine set_multistage_block_candidate_enabled

  
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

    logical :: block_candidate
    logical :: block_state_ready

    call manage_q1_mem

    block_candidate = multistage_block_candidate_enabled
    if (block_candidate) then
       block_state_ready = parallel_block_state_is_ready()
       if (.not. block_state_ready) then
          error stop "guarded RK3 block candidate state is not ready"
       end if
    end if

    if (block_candidate) call begin_block_scalar_divergence_capture
    call routine (q, trend)
    if (block_candidate) then
       call finalize_block_scalar_divergence_capture
       call capture_block_domain_multistage_candidate_tendency(1,3,q)
       call begin_block_domain_multistage_candidate_stage(h/3,1,3)
    end if
    call RK_sub_step (q, trend, h/3, q1)
    call WT_after_step (q1, wav)
    if (block_candidate) then
       call update_bdry(q1,NONE,980)
       call refresh_parallel_block_candidate_boundary_state(q1,1,3)
    end if

    if (block_candidate) call begin_block_scalar_divergence_capture
    call routine (q1, trend)
    if (block_candidate) then
       call finalize_block_scalar_divergence_capture
       call capture_block_domain_multistage_candidate_tendency(2,3,q1)
       call begin_block_domain_multistage_candidate_stage(h/2,2,3)
    end if
    call RK_sub_step (q, trend, h/2, q1)
    call WT_after_step (q1, wav)
    if (block_candidate) then
       call update_bdry(q1,NONE,980)
       call refresh_parallel_block_candidate_boundary_state(q1,2,3)
    end if

    if (block_candidate) call begin_block_scalar_divergence_capture
    call routine (q1, trend)
    if (block_candidate) then
       call finalize_block_scalar_divergence_capture
       call capture_block_domain_multistage_candidate_tendency(3,3,q1)
       call begin_block_domain_multistage_candidate_stage(h,3,3)
    end if
    call RK_sub_step (q, trend, h, q)
    if (block_candidate) then
       call accept_block_domain_multistage_candidate(3)
    end if
    call WT_after_step (q, wav, level_start-1)
    if (block_candidate) then
       call refresh_parallel_block_domain_prognostic_state
    end if
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

    logical :: block_candidate
    logical :: block_state_ready

    call manage_q1_mem

    block_candidate = multistage_block_candidate_enabled
    if (block_candidate) then
       block_state_ready = parallel_block_state_is_ready()
       if (.not. block_state_ready) then
          error stop "guarded RK4 block candidate state is not ready"
       end if
    end if

    if (block_candidate) call begin_block_scalar_divergence_capture
    call routine (q, trend)
    if (block_candidate) then
       call finalize_block_scalar_divergence_capture
       call capture_block_domain_multistage_candidate_tendency(1,4,q)
       call begin_block_domain_multistage_candidate_stage(h/4,1,4)
    end if
    call RK_sub_step (q, trend, h/4, q1)
    call WT_after_step (q1, wav)
    if (block_candidate) then
       call update_bdry(q1,NONE,980)
       call refresh_parallel_block_candidate_boundary_state(q1,1,4)
    end if

    if (block_candidate) call begin_block_scalar_divergence_capture
    call routine (q1, trend)
    if (block_candidate) then
       call finalize_block_scalar_divergence_capture
       call capture_block_domain_multistage_candidate_tendency(2,4,q1)
       call begin_block_domain_multistage_candidate_stage(h/3,2,4)
    end if
    call RK_sub_step (q, trend, h/3, q1)
    call WT_after_step (q1, wav)
    if (block_candidate) then
       call update_bdry(q1,NONE,980)
       call refresh_parallel_block_candidate_boundary_state(q1,2,4)
    end if

    if (block_candidate) call begin_block_scalar_divergence_capture
    call routine (q1, trend)
    if (block_candidate) then
       call finalize_block_scalar_divergence_capture
       call capture_block_domain_multistage_candidate_tendency(3,4,q1)
       call begin_block_domain_multistage_candidate_stage(h/2,3,4)
    end if
    call RK_sub_step (q, trend, h/2, q1)
    call WT_after_step (q1, wav)
    if (block_candidate) then
       call update_bdry(q1,NONE,980)
       call refresh_parallel_block_candidate_boundary_state(q1,3,4)
    end if

    if (block_candidate) call begin_block_scalar_divergence_capture
    call routine (q1, trend)
    if (block_candidate) then
       call finalize_block_scalar_divergence_capture
       call capture_block_domain_multistage_candidate_tendency(4,4,q1)
       call begin_block_domain_multistage_candidate_stage(h,4,4)
    end if
    call RK_sub_step (q, trend, h, q)
    if (block_candidate) then
       call accept_block_domain_multistage_candidate(4)
    end if
    call WT_after_step (q, wav, level_start-1)
    if (block_candidate) then
       call refresh_parallel_block_domain_prognostic_state
    end if
  end subroutine RK4
  

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
  

  subroutine init_RK_mem
    
    implicit none
    
    integer :: d, k, v

    allocate (q1(1:N_VARIABLE,1:zmax))

    do k = 1, zmax
       do v = 1, N_VARIABLE
          call init_Field (q1(v,k), POSIT(v))
       end do

       do d = 1, size(grid)
          do v = 1, N_VARIABLE
             call init (q1(v,k)%data(d), sol(v,k)%data(d)%length)
             q1(v,k)%data(d)%elts = dble (N_VARIABLE-v)
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
