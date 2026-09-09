module time_integr_mod

  use, intrinsic :: iso_fortran_env, only : int64

  use kind_mod,   only : dp
  use shared_mod, only : N_VARIABLE, NONE, POSIT, S_TEMP, S_VELO, eps, level_start, theta2, zlevels, zmax
  use parallel_block_velocity_mod, only : native_velocity_rk, native_velocity
  
  use adapt_mod,         only : WT_after_step
  use barotropic_2d_mod, only : barotropic_correction, eta_update, flux_divergence, scalar_star, u_star, u_update
  use comm_mpi_mod,      only : update_bdry
  use dyn_arrays,        only : extend, init
  use domain_mod,        only : Float_Field, init_Field, grid, sol, trend, wav_coeff
  use multi_level_mod,   only : trend_ml, block_tendency_compatibility_ml
  use parallel_block_mpi_mod, only : &
       prepare_temperature_boundary_stage, apply_temperature_boundary_stage, &
       BLOCK_PROFILE_DOMAIN_TENDENCY, &
       BLOCK_PROFILE_DOMAIN_TENDENCY_COMPATIBILITY, &
       BLOCK_PROFILE_DOMAIN_RK_COMPATIBILITY, &
       begin_block_domain_multistage_candidate_stage, &
       begin_block_scalar_divergence_capture, &
       activate_block_native_inverse_transform, &
       capture_block_domain_multistage_candidate_tendency, &
       finalize_block_scalar_divergence_capture, &
       block_domain_production_writeback_count, &
       block_dynamics_validation_enabled, &
       parallel_block_profile_begin, parallel_block_profile_end, &
       parallel_block_grid_change_is_pending, &
       parallel_block_state_is_ready, &
       prepare_block_native_wavelet_compression, &
       prepare_block_native_multistage_wavelet_acceptance, &
       prepare_block_native_multistage_wavelet_stage, &
       prepare_block_native_velocity_remainder, &
       retain_block_native_multistage_candidate, &
       refresh_parallel_block_candidate_boundary_state, &
       refresh_parallel_block_domain_prognostic_state, &
       activate_block_native_wavelet_compression, &
       validate_candidate_block_outer_vector_wavelets, &
       validate_candidate_block_scalar_wavelets, &
       validate_candidate_block_velocity_restriction

  implicit none

  private
  public :: dt_step, dt_step_split
  public :: init_RK_mem
  public :: set_multistage_block_candidate_enabled
  public :: call_domain_tendency_consumer
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
  integer(int64) :: legacy_domain_tendency_call_count = 0_int64

  
contains


  subroutine call_domain_tendency_consumer (q,routine)
    ! Audit the complete legacy Domain tendency transaction rather than its
    ! per-cell callbacks.  This covers physics_scalar_flux and
    ! physics_velo_source without adding checks inside hot grid-point loops.

    implicit none

    type(Float_Field), intent(inout), target :: &
         q(1:N_VARIABLE,1:zlevels)
    procedure(trend_sub) :: routine

    integer(int64) :: writeback_before
    logical :: pending_before
    logical :: ready_before
    real(dp) :: profile_start

    ready_before = parallel_block_state_is_ready()
    pending_before = parallel_block_grid_change_is_pending()
    if (ready_before .and. pending_before) &
         error stop "tendency callback entered an invalid block phase"
    writeback_before = block_domain_production_writeback_count()

    profile_start = &
         parallel_block_profile_begin(BLOCK_PROFILE_DOMAIN_TENDENCY)
    call routine(q,trend)
    legacy_domain_tendency_call_count = &
         legacy_domain_tendency_call_count + 1_int64
    call parallel_block_profile_end( &
         BLOCK_PROFILE_DOMAIN_TENDENCY,profile_start)

    if (parallel_block_state_is_ready() .neqv. ready_before) &
         error stop "tendency callback changed authoritative block state"
    if (parallel_block_grid_change_is_pending() .neqv. &
         pending_before) &
         error stop "tendency callback changed the grid-change phase"
    if (block_domain_production_writeback_count() /= writeback_before) &
         error stop "tendency callback performed an unaccounted writeback"
  end subroutine call_domain_tendency_consumer


  subroutine prepare_block_multistage_tendency ( &
       domain_stage,routine,scale,stage,stage_count,validate_oracle)
    ! Production evaluates only the retained compatibility inputs.  Oracle
    ! builds execute the complete legacy tendency and keep the existing exact
    ! RK-stage comparison path.

    implicit none

    type(Float_Field), intent(inout), target :: &
         domain_stage(1:N_VARIABLE,1:zlevels)
    procedure(trend_sub) :: routine
    real(dp), intent(in) :: scale
    integer, intent(in) :: stage
    integer, intent(in) :: stage_count
    logical, intent(in) :: validate_oracle

    real(dp) :: profile_start

    call begin_block_scalar_divergence_capture
    profile_start = parallel_block_profile_begin( &
         BLOCK_PROFILE_DOMAIN_TENDENCY_COMPATIBILITY)
    call block_tendency_compatibility_ml(domain_stage,trend)
    call parallel_block_profile_end( &
         BLOCK_PROFILE_DOMAIN_TENDENCY_COMPATIBILITY,profile_start)
    call finalize_block_scalar_divergence_capture
    call prepare_block_native_velocity_remainder(domain_stage)
    if (validate_oracle) &
         call call_domain_tendency_consumer(domain_stage,routine)
    call capture_block_domain_multistage_candidate_tendency( &
         stage,stage_count,domain_stage, &
         compatibility_remainder=.true.)
    call begin_block_domain_multistage_candidate_stage( &
         scale,stage,stage_count)
  end subroutine prepare_block_multistage_tendency


  subroutine assert_multistage_legacy_tendency_calls ( &
       count_before,stage_count,validate_oracle)

    implicit none

    integer(int64), intent(in) :: count_before
    integer, intent(in) :: stage_count
    logical, intent(in) :: validate_oracle

    if (validate_oracle) then
       if (legacy_domain_tendency_call_count-count_before /= &
            int(stage_count,int64)) then
          error stop "oracle RK path did not execute every Domain tendency"
       end if
    else if (legacy_domain_tendency_call_count /= count_before) then
       error stop "optimized RK path executed a legacy Domain tendency"
    end if
  end subroutine assert_multistage_legacy_tendency_calls


  subroutine set_multistage_block_candidate_enabled (enabled)
    ! Guard the production multistage candidate so other RK3 and RK4 callers
    ! remain on the unchanged Domain-only pathway.

    implicit none

    logical, intent(in) :: enabled

    multistage_block_candidate_enabled = enabled

  end subroutine set_multistage_block_candidate_enabled


  subroutine native_provisional_WT (scaling,wavelet,stage,stage_count)
    ! Run one non-final RK transform over level_start:level_end without the
    ! fixed-coarse restriction used only by the accepted final stage.

    implicit none

    type(Float_Field), intent(inout) :: &
         scaling(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: &
         wavelet(1:N_VARIABLE,1:zlevels)
    integer, intent(in) :: stage
    integer, intent(in) :: stage_count

    call prepare_block_native_multistage_wavelet_stage(stage,stage_count)
    call WT_after_step( &
         scaling,wavelet,level_start, &
         validate_scalar_wavelets= &
         validate_candidate_block_scalar_wavelets, &
         validate_outer_vector_wavelets= &
         validate_candidate_block_outer_vector_wavelets, &
         native_wavelet_output=.true., &
         prepare_native_compression= &
         prepare_block_native_wavelet_compression, &
         activate_native_compression= &
         activate_block_native_wavelet_compression, &
         activate_native_inverse= &
         activate_block_native_inverse_transform)
    call refresh_parallel_block_candidate_boundary_state( &
         scaling,stage,stage_count,native_inverse=.true.)

  end subroutine native_provisional_WT

  
  subroutine Euler (q, wav, routine, h)
    ! Euler time step
    ! Stable for CFL<1, first order
    
    implicit none
    
    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
    procedure (trend_sub)            :: routine

    call call_domain_tendency_consumer(q,routine)
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
    logical :: validate_oracle
    integer(int64) :: legacy_call_count_before

    call manage_q1_mem

    block_candidate = multistage_block_candidate_enabled
    validate_oracle = .false.
    legacy_call_count_before = legacy_domain_tendency_call_count
    if (block_candidate) then
       block_state_ready = parallel_block_state_is_ready()
       if (.not. block_state_ready) then
          error stop "guarded RK3 block candidate state is not ready"
       end if
       validate_oracle = block_dynamics_validation_enabled()
    end if

    if (block_candidate) then
       call prepare_block_multistage_tendency( &
            q,routine,h/3,1,3,validate_oracle)
    else
       call call_domain_tendency_consumer(q,routine)
    end if
    ! Retain mass/velocity compatibility and initialize temperature scaffolding.
    ! Native temperature publication then supplies the wavelet/halo input;
    ! its block-owned RK candidate must not be left only in compact storage.
    if (block_candidate) then
       call RK_sub_step_compatibility(q,trend,h/3,q1)
    else
       call RK_sub_step(q,trend,h/3,q1)
    end if
    if (block_candidate) then
       call retain_block_native_multistage_candidate(q1,1,3)
    end if
    if (block_candidate) then
       call native_provisional_WT(q1,wav,1,3)
    else
       call WT_after_step(q1,wav)
    end if

    if (block_candidate) then
       call prepare_block_multistage_tendency( &
            q1,routine,h/2,2,3,validate_oracle)
    else
       call call_domain_tendency_consumer(q1,routine)
    end if
    if (block_candidate) then
       call RK_sub_step_compatibility(q,trend,h/2,q1)
    else
       call RK_sub_step(q,trend,h/2,q1)
    end if
    if (block_candidate) then
       call retain_block_native_multistage_candidate(q1,2,3)
    end if
    if (block_candidate) then
       call native_provisional_WT(q1,wav,2,3)
    else
       call WT_after_step(q1,wav)
    end if

    if (block_candidate) then
       call prepare_block_multistage_tendency( &
            q1,routine,h,3,3,validate_oracle)
    else
       call call_domain_tendency_consumer(q1,routine)
    end if
    if (block_candidate) then
       call RK_sub_step_compatibility(q,trend,h,q)
    else
       call RK_sub_step(q,trend,h,q)
    end if
    if (block_candidate) then
       call retain_block_native_multistage_candidate(q,3,3)
       call prepare_block_native_multistage_wavelet_acceptance(3)
    end if
    call WT_after_step( &
         q,wav,level_start-1, &
         validate_scalar_wavelets= &
         validate_candidate_block_scalar_wavelets, &
         validate_outer_vector_wavelets= &
         validate_candidate_block_outer_vector_wavelets, &
         validate_velocity_restriction= &
         validate_candidate_block_velocity_restriction, &
         native_wavelet_output=block_candidate, &
         prepare_native_compression= &
         prepare_block_native_wavelet_compression, &
         activate_native_compression= &
         activate_block_native_wavelet_compression, &
         activate_native_inverse= &
         activate_block_native_inverse_transform)
    if (block_candidate) then
       call refresh_parallel_block_domain_prognostic_state(native_inverse=.true.)
       call assert_multistage_legacy_tendency_calls( &
            legacy_call_count_before,3,validate_oracle)
    end if
  end subroutine RK3


  subroutine RK4 (q, wav, routine, h)
    ! Low-storage four-stage second-order Runge-Kutta scheme used in
    ! Dubos et al. (2015), Geosci. Model Dev., 8, 3131-3150.
    ! Fourth order accurate for linear equations, stable for CFL <= 2*sqrt(2) ~ 2.83.
    ! Does not require extra solution variables.
    
    implicit none

    real(dp),          intent(in)    :: h
    type(Float_Field), intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: wav(1:N_VARIABLE,1:zlevels)
    procedure (trend_sub)            :: routine

    logical :: block_candidate
    logical :: block_state_ready
    logical :: validate_oracle
    integer(int64) :: legacy_call_count_before

    call manage_q1_mem

    block_candidate = multistage_block_candidate_enabled
    validate_oracle = .false.
    legacy_call_count_before = legacy_domain_tendency_call_count
    if (block_candidate) then
       block_state_ready = parallel_block_state_is_ready()
       if (.not. block_state_ready) then
          error stop "guarded RK4 block candidate state is not ready"
       end if
       validate_oracle = block_dynamics_validation_enabled()
    end if

    if (block_candidate) then
       call prepare_block_multistage_tendency( &
            q,routine,h/4,1,4,validate_oracle)
    else
       call call_domain_tendency_consumer(q,routine)
    end if
    ! Retain mass/velocity compatibility and initialize temperature scaffolding.
    ! Native temperature publication then supplies the wavelet/halo input;
    ! its block-owned RK candidate must not be left only in compact storage.
    if (block_candidate) then
       call RK_sub_step_compatibility(q,trend,h/4,q1)
    else
       call RK_sub_step(q,trend,h/4,q1)
    end if
    if (block_candidate) then
       call retain_block_native_multistage_candidate(q1,1,4)
    end if
    if (block_candidate) then
       call native_provisional_WT(q1,wav,1,4)
    else
       call WT_after_step(q1,wav)
    end if

    if (block_candidate) then
       call prepare_block_multistage_tendency( &
            q1,routine,h/3,2,4,validate_oracle)
    else
       call call_domain_tendency_consumer(q1,routine)
    end if
    if (block_candidate) then
       call RK_sub_step_compatibility(q,trend,h/3,q1)
    else
       call RK_sub_step(q,trend,h/3,q1)
    end if
    if (block_candidate) then
       call retain_block_native_multistage_candidate(q1,2,4)
    end if
    if (block_candidate) then
       call native_provisional_WT(q1,wav,2,4)
    else
       call WT_after_step(q1,wav)
    end if

    if (block_candidate) then
       call prepare_block_multistage_tendency( &
            q1,routine,h/2,3,4,validate_oracle)
    else
       call call_domain_tendency_consumer(q1,routine)
    end if
    if (block_candidate) then
       call RK_sub_step_compatibility(q,trend,h/2,q1)
    else
       call RK_sub_step(q,trend,h/2,q1)
    end if
    if (block_candidate) then
       call retain_block_native_multistage_candidate(q1,3,4)
    end if
    if (block_candidate) then
       call native_provisional_WT(q1,wav,3,4)
    else
       call WT_after_step(q1,wav)
    end if

    if (block_candidate) then
       call prepare_block_multistage_tendency( &
            q1,routine,h,4,4,validate_oracle)
    else
       call call_domain_tendency_consumer(q1,routine)
    end if
    if (block_candidate) then
       call RK_sub_step_compatibility(q,trend,h,q)
    else
       call RK_sub_step(q,trend,h,q)
    end if
    if (block_candidate) then
       call retain_block_native_multistage_candidate(q,4,4)
       call prepare_block_native_multistage_wavelet_acceptance(4)
    end if
    call WT_after_step( &
         q,wav,level_start-1, &
         validate_scalar_wavelets= &
         validate_candidate_block_scalar_wavelets, &
         validate_outer_vector_wavelets= &
         validate_candidate_block_outer_vector_wavelets, &
         validate_velocity_restriction= &
         validate_candidate_block_velocity_restriction, &
         native_wavelet_output=block_candidate, &
         prepare_native_compression= &
         prepare_block_native_wavelet_compression, &
         activate_native_compression= &
         activate_block_native_wavelet_compression, &
         activate_native_inverse= &
         activate_block_native_inverse_transform)
    if (block_candidate) then
       call refresh_parallel_block_domain_prognostic_state(native_inverse=.true.)
       call assert_multistage_legacy_tendency_calls( &
            legacy_call_count_before,4,validate_oracle)
    end if
  end subroutine RK4


  subroutine RK_sub_step_compatibility (sols,trends,h,dest)
    ! Materialize mass compatibility and native velocity closure; preserve
    ! temperature scaffold values for native boundary/stage publication.

    implicit none

    real(dp), intent(in) :: h
    type(Float_Field), intent(in) :: &
         sols(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in) :: &
         trends(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: &
         dest(1:N_VARIABLE,1:zlevels)

    real(dp) :: profile_start
    type(Float_Field), allocatable :: temperature_reference(:,:)
    integer :: d,k,ibeg,iend

    profile_start = parallel_block_profile_begin( &
         BLOCK_PROFILE_DOMAIN_RK_COMPATIBILITY)
    call prepare_temperature_boundary_stage(sols,h)
    ! Compute before the legacy oracle: on the last RK substage sols and
    ! dest are the same field. Computing afterwards would advance it twice.
    do d=1,size(grid)
       ibeg=3*grid(d)%patch%elts(3)%elts_start+1
       iend=dest(S_VELO,1)%data(d)%length
       do k=1,zlevels
          call native_velocity_rk(d,k,ibeg,sols(S_VELO,k)%data(d)%elts(ibeg:iend),h, &
               native_velocity(d)%rk(ibeg:iend,k))
       end do
    end do
    call RK_sub_step(sols,trends,h,dest, &
         native_temperature=.not. block_dynamics_validation_enabled(), &
         native_velocity=.not. block_dynamics_validation_enabled())
    ! All velocity values needed by the compatibility consumer, including
    ! boundary/scaffold slots, come from the native workspace. The final-owner
    ! stage publication still supplies the integrated owned block state.
    do d=1,size(grid)
       ibeg=3*grid(d)%patch%elts(3)%elts_start+1
       iend=dest(S_VELO,1)%data(d)%length
       do k=1,zlevels
          if (block_dynamics_validation_enabled()) then
             if (any(transfer(native_velocity(d)%rk(ibeg:iend,k),[0_int64],iend-ibeg+1) /= &
                  transfer(dest(S_VELO,k)%data(d)%elts(ibeg:iend),[0_int64],iend-ibeg+1))) &
                  error stop "native velocity RK compatibility differs bit-for-bit"
          end if
          dest(S_VELO,k)%data(d)%elts(ibeg:iend)=native_velocity(d)%rk(ibeg:iend,k)
       end do
    end do
    if (block_dynamics_validation_enabled()) temperature_reference=dest(S_TEMP:S_TEMP,:)
    call apply_temperature_boundary_stage(dest)
    if (allocated(temperature_reference)) call assert_temperature_boundary_stage(dest,temperature_reference)
    call parallel_block_profile_end( &
         BLOCK_PROFILE_DOMAIN_RK_COMPATIBILITY,profile_start)
  end subroutine RK_sub_step_compatibility

  subroutine assert_temperature_boundary_stage(candidate,reference)
    ! Compare at the actual halo-completed consumer boundary. Individual
    ! storage aliases can be overwritten by communication before consumption.
    use arch_mod, only : rank
    type(Float_Field),intent(inout)::candidate(1:N_VARIABLE,1:zlevels),reference(:,:)
    integer :: d,k,i
    real(dp)::difference,allowed
    call update_bdry(reference,NONE,1260)
    call update_bdry(candidate(S_TEMP:S_TEMP,:),NONE,1261)
    do d=1,size(grid)
       do k=1,zlevels
          do i=grid(d)%patch%elts(3)%elts_start+1,candidate(S_TEMP,k)%data(d)%length
             difference=abs(candidate(S_TEMP,k)%data(d)%elts(i)-reference(1,k)%data(d)%elts(i))
             allowed=64*epsilon(1.0_dp)*max(1.0_dp,abs(reference(1,k)%data(d)%elts(i)))
             if (difference <= allowed) cycle
             write(6,*) 'Native temperature RK boundary mismatch: rank, Domain, level, node = ',rank,d,k,i-1
             write(6,*) 'native, reference = ',candidate(S_TEMP,k)%data(d)%elts(i),reference(1,k)%data(d)%elts(i)
             error stop "native temperature RK boundary differs"
          end do
       end do
    end do
  end subroutine assert_temperature_boundary_stage



  subroutine RK_sub_step (sols, trends, h, dest, native_temperature, native_velocity)
    
    implicit none
    
    real(dp),          intent(in)    :: h
    type(Float_Field), intent(in)    :: sols(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(in)    :: trends(1:N_VARIABLE,1:zlevels)
    type(Float_Field), intent(inout) :: dest(1:N_VARIABLE,1:zlevels)
    logical, optional, intent(in) :: native_temperature,native_velocity
    
    integer :: d, ibeg, iend, k, v
    logical :: copy_temperature,copy_velocity

    copy_temperature=.false.
    if (present(native_temperature)) copy_temperature=native_temperature
    copy_velocity=.false.
    if (present(native_velocity)) copy_velocity=native_velocity

    do v = 1, N_VARIABLE
       do d = 1, size(grid)
          ibeg = (1+2*(POSIT(v)-1)) * grid(d)%patch%elts(2+1)%elts_start + 1
          iend = dest(v,1)%data(d)%length
          do k = 1, zlevels
             if (v==S_VELO .and. copy_velocity) cycle
             if (v == S_TEMP .and. copy_temperature) then
                ! Seed untouched scaffolding. The native boundary adapter and
                ! integrated stage publication supply all evolved values.
                dest(v,k)%data(d)%elts(ibeg:iend)=sols(v,k)%data(d)%elts(ibeg:iend)
                cycle
             end if
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
