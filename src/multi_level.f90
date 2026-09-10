module multi_level_mod
  use parallel_block_profile_mod
  use mpi_f08
  use arch_mod, only : comm, MPI_DP, rank, n_process, glo_id
  use shared_mod, only : n_domain
  use parallel_block_mass_mod
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_is_nan
  
  use kind_mod,   only : dp
  use shared_mod, only : bfly_no2, nghb_pt, hex_sides, hex_s_offs, N_VARIABLE, zlevels, N_BDRY,  N_CHDRN, &
       LORT, UPLT, TRIAG, &
       ADJZONE, EDGE, MM, MP, PM, PP, UMZ, UPZ, UZM, UZP, VMM, VMPP, VMP, VPM, VPMM, VPP, WMM, WMP, WPM, WPP, WMMM, WPPP, &
       level_end, level_start, NONE,RT, DG, UP, z_null, RESTRCT, TRSK, &
       S_MASS, S_TEMP, S_VELO, S_DIVU, S_ROTU, &
       Laplace_divu, Laplace_rotu, Laplace_sclr, AT_EDGE, AT_NODE, scalars, eps, radius

  use comm_mpi_mod,    only : update_bdry, update_bdry__start, update_bdry__finish
  use diagnostics_mod, only : cal_div, cal_surf_press, gradi_e, &
       integrate_pressure_up, post_vort
  use domain_ops_mod,  only : apply_interscale_to_patch, apply_interscale_to_patch3, apply_onescale_to_patch, apply_to_penta_d 
  use init_mod,        only : physics_scalar_flux, &
       physics_velo_source, u_source
  use ops_mod,         only : du_grad, du_source, &
       post_step1, Qperp, scalar_trend, scalar_trend_pair, step1
  use patch_mod,       only : PATCH_SIZE
  use parallel_block_velocity_mod, only : velocity_dx, velocity_dy, &
       velocity_weights, velocity_qperp, velocity_source, &
       velocity_restrict_source, velocity_gradient, &
       Native_Velocity_Workspace, native_velocity, native_velocity_work, &
       native_velocity_transaction, native_velocity_oracle, &
       VELOCITY_DIRECT, VELOCITY_RESTRICT, validate_velocity_program, &
       execute_velocity_sources, execute_velocity_gradients
  use parallel_block_mpi_mod, only : &
       BLOCK_PROFILE_NATIVE_MASS_PLAN, BLOCK_PROFILE_NATIVE_MASS, &
       BLOCK_PROFILE_DOMAIN_MASS_COMPATIBILITY, &
       BLOCK_PROFILE_DOMAIN_OPERATOR_COMPATIBILITY, &
       BLOCK_PROFILE_DOMAIN_VELOCITY_COMPATIBILITY, &
       BLOCK_PROFILE_NATIVE_VELOCITY_PLAN, BLOCK_PROFILE_NATIVE_VELOCITY_SOURCE, BLOCK_PROFILE_NATIVE_VELOCITY_GRADIENT, &
       begin_block_velocity_source_transport, &
       block_dynamics_validation_enabled, &
       native_velocity_plan_generation, &
       block_scalar_capture_active, capture_block_scalar_physics_patch, &
       capture_block_scalar_divergence_level, &
       capture_block_velocity_source_level, &
       finalize_block_velocity_source_transport, &
       parallel_block_profile_begin, parallel_block_profile_end
  use utils_mod,       only : zero_float
  
  use domain_mod, only : Domain, Float_Field, get_offs_Domain, grid, &
       bernoulli, exner, exner_fun, ke, qe, mean_m, mean_t, sol, sol_mean, &
       horiz_flux, h_mflux, h_flux, divu, dvelo, mass, temp, velo, scalar, &
       vort, &
       Laplacian_scalar, Laplacian_vector, dscalar, Laplacian, &
       ed_idx, idx, idx2

    
  implicit none

  private
  public :: cpt_or_restr_flux, trend_ml, &
       block_tendency_compatibility_ml, cal_divu_ml

  type :: Velocity_Source_Measurement_Type
     real(dp), allocatable :: qperp(:)
     real(dp), allocatable :: physics(:)
     real(dp), allocatable :: edge_length(:)
     real(dp), allocatable :: integrated_source(:)
     ! Independent Stage 175 oracle workspace. Never initialized from dvelo.
     real(dp), allocatable :: native_source(:)
     logical, allocatable :: active(:)
     logical, allocatable :: covered(:)
     logical, allocatable :: direct(:)
  end type Velocity_Source_Measurement_Type

  type(Velocity_Source_Measurement_Type), allocatable, save :: &
       velocity_source_measurement(:)
  ! Serial oracle traversal context. Module callbacks avoid executable-stack
  ! trampolines on cluster builds; production does not enter this traversal.
  integer, save :: native_velocity_measurement_domain = 0
  integer, save :: velocity_plan_domain=0, velocity_plan_level=0, velocity_plan_patch=0
  integer, save :: velocity_plan_pass=0, velocity_plan_action=0, velocity_plan_direct=0, velocity_plan_restriction=0

  type :: Mass_Boundary_Plan
     integer, allocatable :: send_count(:),recv_count(:),send_displ(:),recv_displ(:)
     integer, allocatable :: source(:,:),destination(:,:),local(:,:)
     real(dp), allocatable :: send_value(:),recv_value(:)
  end type
  type(Mass_Boundary_Plan), allocatable, save :: mass_boundary(:,:)
  integer(int64), save :: mass_boundary_generation=-1_int64
  integer, save :: mass_plan_domain=0,mass_plan_level=0,mass_plan_pass=0,mass_plan_count=0

  
contains


  subroutine begin_velocity_source_measurement

    implicit none

    integer :: d
    integer :: n

    if (allocated(velocity_source_measurement)) then
       if (size(velocity_source_measurement) /= size(grid)) &
            deallocate(velocity_source_measurement)
    end if
    if (.not. allocated(velocity_source_measurement)) &
         allocate(velocity_source_measurement(size(grid)))
    do d = 1,size(grid)
       n = EDGE*grid(d)%node%length
       if (allocated(velocity_source_measurement(d)%qperp)) then
          if (size(velocity_source_measurement(d)%qperp) /= n) then
             deallocate(velocity_source_measurement(d)%qperp)
             deallocate(velocity_source_measurement(d)%physics)
             deallocate(velocity_source_measurement(d)%edge_length)
             deallocate(velocity_source_measurement(d)%integrated_source)
             deallocate(velocity_source_measurement(d)%native_source)
             deallocate(velocity_source_measurement(d)%active)
             deallocate(velocity_source_measurement(d)%covered)
             deallocate(velocity_source_measurement(d)%direct)
          end if
       end if
       if (.not. allocated(velocity_source_measurement(d)%qperp)) then
          allocate(velocity_source_measurement(d)%qperp(n))
          allocate(velocity_source_measurement(d)%physics(n))
          allocate(velocity_source_measurement(d)%edge_length(n))
          allocate(velocity_source_measurement(d)%integrated_source(n))
          allocate(velocity_source_measurement(d)%native_source(n))
          allocate(velocity_source_measurement(d)%active(n))
          allocate(velocity_source_measurement(d)%covered(n))
          allocate(velocity_source_measurement(d)%direct(n))
       end if
       velocity_source_measurement(d)%qperp = 0.0_dp
       velocity_source_measurement(d)%physics = 0.0_dp
       velocity_source_measurement(d)%edge_length = 0.0_dp
       velocity_source_measurement(d)%integrated_source = 0.0_dp
       velocity_source_measurement(d)%native_source = 0.0_dp
       velocity_source_measurement(d)%active = .false.
       velocity_source_measurement(d)%covered = .false.
       velocity_source_measurement(d)%direct = .false.
    end do
    call begin_block_velocity_source_transport

  end subroutine begin_velocity_source_measurement


  subroutine finish_velocity_source_measurement (field_level)

    implicit none

    integer, intent(in) :: field_level

    integer :: d

    do d = 1,size(grid)
       if (.not. allocated(velocity_source_measurement(d)%covered)) &
            error stop "velocity-source measurement storage is absent"
       if (.not. any(velocity_source_measurement(d)%covered)) &
            error stop "velocity-source measurement coverage is empty"
       call capture_block_velocity_source_level( &
            d,field_level,velocity_source_measurement(d)%qperp, &
            velocity_source_measurement(d)%physics, &
            velocity_source_measurement(d)%edge_length, &
            velocity_source_measurement(d)%integrated_source, &
            velocity_source_measurement(d)%active, &
            velocity_source_measurement(d)%covered, &
            velocity_source_measurement(d)%direct)
    end do
    call finalize_block_velocity_source_transport

  end subroutine finish_velocity_source_measurement

  
  subroutine trend_ml (q, dq)
    ! Compute trends of prognostic variables assuming Lagrangian vertical coordinates
    
    implicit none
    
    type(Float_Field), intent(inout), target :: q(1:N_VARIABLE,1:zlevels), dq(1:N_VARIABLE,1:zlevels)

    integer :: k, l

    call update_bdry (q, NONE, 967)

    ! Initialize trends
    call zero_float (dq)

    ! Compute surface pressure on all grids
    call cal_surf_press (q(1:N_VARIABLE,1:zlevels))

    ! Compute each vertical level starting from surface
    do k = 1, zlevels
       if (Laplace_divu /= 0) call cal_divu_ml (q(S_VELO,k))
       if (Laplace_sclr == 2) call cal_Laplacian_scalars (q, k)
       if (Laplace_divu == 2) call cal_Laplacian_divu ! requires divu

       ! Calculate trend on all scales, from fine to coarse
       do l = level_end, level_start, -1
          ! Finish non-blocking communication of dq from level (l+1)
          if (l < level_end) then
             call update_bdry__finish( &
                  dq(scalars(1):scalars(2),k),l+1)
             call capture_block_scalar_divergence_level( &
                  q,physics_scalar_flux,0,k,l+1, &
                  domain_tendency=dq,dscalar_only=.true.)
          end if

          call basic_operators  (q, dq, k, l, .false.)
          call cal_scalar_trend (q, dq, k, l, .false.)

          ! Start non-blocking communication of dq for use at next level (l-1)
          if (level_start /= level_end .and. l > level_start) call update_bdry__start (dq(scalars(1):scalars(2),k),l) 

          call velocity_trend_source (q, dq, k, l)
       end do
       call velocity_trend_grad (q, dq, k)
    end do
    dq%bdry_uptodate = .false.
  end subroutine trend_ml


  subroutine block_tendency_compatibility_ml (q, dq)
    ! Shared primitive/physics producer with native mass and velocity chains.
    ! Domain scalar restriction/divergence runs only for the oracle. Native
    ! mass executes on geometry owners at the velocity-consumption phase;
    ! temperature restriction remains on final owners. trend_ml is the
    ! independent complete validation oracle.

    implicit none

    type(Float_Field), intent(inout), target :: &
         q(1:N_VARIABLE,1:zlevels), dq(1:N_VARIABLE,1:zlevels)

    integer :: k, l, compatibility_last, velocity_domain
#ifdef WAVETRISK_TEST_MASS_CUT
    type :: Mass_Scratch_Backup
       real(dp), allocatable :: value(:)
    end type
    type(Mass_Scratch_Backup), allocatable :: mass_scratch(:)
#endif
#ifdef WAVETRISK_TEST_TEMPERATURE_CUT
    integer :: poison_domain
    type :: Temperature_Scratch_Backup
       real(dp), allocatable :: value(:)
    end type
    type(Temperature_Scratch_Backup), allocatable :: temperature_scratch(:)
#endif

    logical :: validate_velocity_source

    real(dp) :: profile_start

    call detail_enter(DP_SHARED)
    call detail_enter(DP_PRIMITIVE)
    call update_bdry(q,NONE,1067)
    call zero_float(dq)
    call cal_surf_press(q(1:N_VARIABLE,1:zlevels))
    call detail_leave(DP_PRIMITIVE)
    validate_velocity_source = block_dynamics_validation_enabled()
    profile_start=parallel_block_profile_begin(BLOCK_PROFILE_NATIVE_VELOCITY_PLAN)
    call prepare_native_velocity_programs
    call parallel_block_profile_end(BLOCK_PROFILE_NATIVE_VELOCITY_PLAN,profile_start)
    profile_start=parallel_block_profile_begin(BLOCK_PROFILE_NATIVE_MASS_PLAN)
    call prepare_native_mass
    call parallel_block_profile_end(BLOCK_PROFILE_NATIVE_MASS_PLAN,profile_start)
    native_velocity_transaction = .true.
    native_velocity_oracle = validate_velocity_source
    mass_transaction=.true.
    mass_oracle=validate_velocity_source
#ifdef WAVETRISK_TEST_MASS_CUT
    if (.not.mass_oracle) then
       allocate(mass_scratch(size(grid)))
       do velocity_domain=1,size(grid)
          do k=1,zlevels
             dq(S_MASS,k)%data(velocity_domain)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
          end do
       end do
    end if
#endif
    do velocity_domain=1,size(native_mass)
       native_mass(velocity_domain)%tendency=0.0_dp
       native_mass(velocity_domain)%ready=.false.
    end do
    do velocity_domain=1,size(native_velocity)
       native_velocity(velocity_domain)%ready=.false.
       native_velocity(velocity_domain)%tendency=0.0_dp
    end do
#ifdef WAVETRISK_TEST_VELOCITY_CUT
    if (.not. validate_velocity_source) then
       do velocity_domain=1,size(grid)
          do k=1,zlevels
             dq(S_VELO,k)%data(velocity_domain)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
          end do
       end do
    end if
#endif
    compatibility_last=S_MASS
    if (validate_velocity_source) compatibility_last=scalars(2)
#ifdef WAVETRISK_TEST_TEMPERATURE_CUT
    if (.not. validate_velocity_source) then
       allocate(temperature_scratch(size(grid)))
       do poison_domain=1,size(grid)
          do k=1,zlevels
             dq(S_TEMP,k)%data(poison_domain)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
          end do
       end do
    end if
#endif

    do k = 1,zlevels
       do velocity_domain=1,size(native_velocity)
          native_velocity(velocity_domain)%source=0.0_dp
       end do
       if (validate_velocity_source) &
            call begin_velocity_source_measurement
       if (Laplace_divu /= 0) call cal_divu_ml(q(S_VELO,k))
       if (Laplace_sclr == 2) call cal_Laplacian_scalars(q,k)
       if (Laplace_divu == 2) call cal_Laplacian_divu
#ifdef WAVETRISK_TEST_MASS_CUT
       if (.not.mass_oracle) then
          do velocity_domain=1,size(grid)
             mass_scratch(velocity_domain)%value=horiz_flux(S_MASS)%data(velocity_domain)%elts
             horiz_flux(S_MASS)%data(velocity_domain)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
          end do
       end if
#endif
#ifdef WAVETRISK_TEST_TEMPERATURE_CUT
       ! The retained physics interface uses this scratch field to form the
       ! biharmonic diffusion input. Poison after that single physics prepass:
       ! no advective flux, restriction or tendency may subsequently use it.
       if (.not. validate_velocity_source) then
          do poison_domain=1,size(grid)
             temperature_scratch(poison_domain)%value=horiz_flux(S_TEMP)%data(poison_domain)%elts
             horiz_flux(S_TEMP)%data(poison_domain)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
          end do
       end if
#endif

       do l = level_end,level_start,-1
          ! Complete the fine-level RHS before coarse restriction. The native
          ! graph preserves the phase of both local copies and remote delivery.
          if (l < level_end.and.validate_velocity_source) then
             call update_bdry__finish( &
                  dq(scalars(1):compatibility_last,k),l+1)
             call capture_block_scalar_divergence_level( &
                  q,physics_scalar_flux,0,k,l+1, &
                  domain_tendency=dq,dscalar_only=.true.)
          end if
          if (l<level_end) then
             profile_start=parallel_block_profile_begin(BLOCK_PROFILE_NATIVE_MASS)
             call exchange_native_mass(AT_NODE,l+1,k)
             call parallel_block_profile_end(BLOCK_PROFILE_NATIVE_MASS,profile_start)
          end if
          profile_start = parallel_block_profile_begin( &
               BLOCK_PROFILE_DOMAIN_OPERATOR_COMPATIBILITY)
          call basic_operators( &
               q,dq,k,l,.not. validate_velocity_source)
          call parallel_block_profile_end( &
               BLOCK_PROFILE_DOMAIN_OPERATOR_COMPATIBILITY,profile_start)
          profile_start = parallel_block_profile_begin( &
               BLOCK_PROFILE_NATIVE_MASS)
          call exchange_native_mass(AT_EDGE,l,k)
          do velocity_domain=1,size(native_mass)
             call execute_mass_divergence(native_mass(velocity_domain)%level(l)%divergence, &
                  native_mass(velocity_domain)%flux,native_mass(velocity_domain)%tendency(:,k))
          end do
          call parallel_block_profile_end(BLOCK_PROFILE_NATIVE_MASS,profile_start)
          if (validate_velocity_source) then
             profile_start=parallel_block_profile_begin(BLOCK_PROFILE_DOMAIN_MASS_COMPATIBILITY)
             call cal_scalar_trend_compatibility(q,dq,k,l)
             call assert_native_mass_level(dq,k,l)
             call parallel_block_profile_end(BLOCK_PROFILE_DOMAIN_MASS_COMPATIBILITY,profile_start)
          else
             call capture_block_scalar_divergence_level(q,physics_scalar_flux,0,k,l)
          end if
          if (level_start /= level_end .and. l > level_start.and.validate_velocity_source) then
             call update_bdry__start( &
                  dq(scalars(1):compatibility_last,k),l)
          end if
          profile_start = parallel_block_profile_begin( &
               BLOCK_PROFILE_NATIVE_VELOCITY_SOURCE)
          call compute_native_velocity_source(q,k,l)
          call parallel_block_profile_end( &
               BLOCK_PROFILE_NATIVE_VELOCITY_SOURCE,profile_start)
          if (validate_velocity_source) then
             profile_start=parallel_block_profile_begin(BLOCK_PROFILE_DOMAIN_VELOCITY_COMPATIBILITY)
             call velocity_trend_source(q,dq,k,l,.true.)
             call parallel_block_profile_end(BLOCK_PROFILE_DOMAIN_VELOCITY_COMPATIBILITY,profile_start)
          end if
       end do

       do velocity_domain=1,size(native_mass)
          native_mass(velocity_domain)%ready(k)=.true.
       end do

       ! Native gradient consumes FINAL restricted B/Exner from the shared
       ! primitive pass. No completed Domain source or gradient is an input.
       profile_start = parallel_block_profile_begin( &
            BLOCK_PROFILE_NATIVE_VELOCITY_GRADIENT)
       call compute_native_velocity_gradient(q,k)
       call parallel_block_profile_end( &
            BLOCK_PROFILE_NATIVE_VELOCITY_GRADIENT,profile_start)
       if (validate_velocity_source) then
          profile_start=parallel_block_profile_begin(BLOCK_PROFILE_DOMAIN_VELOCITY_COMPATIBILITY)
          call velocity_trend_grad(q,dq,k,.true.)
          call parallel_block_profile_end(BLOCK_PROFILE_DOMAIN_VELOCITY_COMPATIBILITY,profile_start)
       end if
       if (validate_velocity_source) &
            call finish_velocity_source_measurement(k)
#ifdef WAVETRISK_TEST_MASS_CUT
       if (.not.mass_oracle) then
          do velocity_domain=1,size(grid)
             if (.not.all(ieee_is_nan(dq(S_MASS,k)%data(velocity_domain)%elts))) &
                  error stop "native mass cut wrote Domain mass tendency"
             if (.not.all(ieee_is_nan(horiz_flux(S_MASS)%data(velocity_domain)%elts))) &
                  error stop "native mass cut wrote Domain mass flux"
             horiz_flux(S_MASS)%data(velocity_domain)%elts=mass_scratch(velocity_domain)%value
          end do
       end if
#endif
#ifdef WAVETRISK_TEST_TEMPERATURE_CUT
    if (.not. validate_velocity_source) then
       do poison_domain=1,size(grid)
          if (.not. all(ieee_is_nan(horiz_flux(S_TEMP)%data(poison_domain)%elts))) &
               error stop "temperature cut wrote the Domain flux workspace"
          if (.not. all(ieee_is_nan(dq(S_TEMP,k)%data(poison_domain)%elts))) &
               error stop "temperature cut wrote the Domain trend workspace"
          ! The next physical layer's diffusion prepass owns this scratch.
          horiz_flux(S_TEMP)%data(poison_domain)%elts=temperature_scratch(poison_domain)%value
       end do
    end if
#endif
    end do
#ifdef WAVETRISK_TEST_VELOCITY_CUT
    if (.not. validate_velocity_source) then
       do velocity_domain=1,size(grid)
          do k=1,zlevels
             if (.not. all(ieee_is_nan(dq(S_VELO,k)%data(velocity_domain)%elts))) &
                  error stop "native velocity cut wrote the Domain tendency"
          end do
       end do
    end if
#endif
    native_velocity_transaction=.false.
    mass_transaction=.false.
    dq%bdry_uptodate = .false.
    call detail_leave(DP_SHARED)

  end subroutine block_tendency_compatibility_ml


  subroutine prepare_native_mass
    integer :: d,l,p,c,s,pass,n
    integer(int64) :: generation
    generation=native_velocity_plan_generation()
    if (allocated(native_mass)) then
       if (size(native_mass)/=size(grid)) deallocate(native_mass)
    end if
    if (.not. allocated(native_mass)) allocate(native_mass(size(grid)))
    do d=1,size(grid)
       n=grid(d)%node%length
       if (native_mass(d)%generation==generation) then
          if (size(native_mass(d)%active)/=n) error stop "native mass layout changed without generation"
          if (all(native_mass(d)%active .eqv. (grid(d)%mask_n%elts(1:n)>=TRSK)) .and. &
               all(native_mass(d)%restricted .eqv. (grid(d)%mask_e%elts(1:EDGE*n)>=RESTRCT))) cycle
       end if
       native_mass(d)=Native_Mass_Workspace()
       associate(w=>native_mass(d))
       allocate(w%level(level_start:level_end),w%flux(EDGE*n),w%tendency(n,zlevels),w%rk(n,zlevels))
       allocate(w%active(n),w%restricted(EDGE*n),w%ready(zlevels))
       w%active=grid(d)%mask_n%elts(1:n)>=TRSK
       w%restricted=grid(d)%mask_e%elts(1:EDGE*n)>=RESTRCT
       w%flux=0.0_dp
       w%tendency=0.0_dp
       w%ready=.false.
       mass_plan_domain=d
       do l=level_end,level_start,-1
          mass_plan_level=l
          do pass=1,2
             mass_plan_pass=pass
             mass_plan_count=0
             if (l<level_end) then
                do s=1,grid(d)%lev(l)%length
                   p=grid(d)%lev(l)%elts(s)
                   do c=1,N_CHDRN
                      if (grid(d)%patch%elts(p+1)%children(c)<=0) cycle
                      call apply_interscale_to_patch3(compile_mass_restriction,grid(d),p,c,z_null,0,1)
                   end do
                end do
             end if
             if (pass==1) allocate(w%level(l)%restriction(mass_plan_count))
          end do
          do pass=1,2
             mass_plan_pass=pass
             mass_plan_count=0
             do s=1,grid(d)%lev(l)%length
                call apply_onescale_to_patch(compile_mass_divergence,grid(d),grid(d)%lev(l)%elts(s),z_null,0,1)
             end do
             if (pass==1) allocate(w%level(l)%divergence(mass_plan_count))
          end do
       end do
       w%generation=generation
       mass_work(1)=mass_work(1)+1_int64
       end associate
    end do
    if (mass_boundary_generation/=generation) then
       if (allocated(mass_boundary)) deallocate(mass_boundary)
       allocate(mass_boundary(AT_NODE:AT_EDGE,level_start:level_end))
       do l=level_start,level_end
          call compile_mass_boundary(AT_NODE,l)
          call compile_mass_boundary(AT_EDGE,l)
       end do
       mass_boundary_generation=generation
    end if
  end subroutine

  subroutine compile_mass_divergence(dom,i,j,k,offs,dims)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer :: id,iw,is,isw
    mass_plan_count=mass_plan_count+1
    if (mass_plan_pass==1) return
    id=idx(i,j,offs,dims)
    associate(s=>native_mass(mass_plan_domain)%level(mass_plan_level)%divergence(mass_plan_count))
    s%target=id+1
    s%active=dom%mask_n%elts(id+1)>=TRSK
    if (.not. s%active) return
    iw=idx(i-1,j,offs,dims)
    is=idx(i,j-1,offs,dims)
    isw=idx(i-1,j-1,offs,dims)
    s%edge=[EDGE*id+RT,EDGE*iw+RT,EDGE*isw+DG,EDGE*id+DG,EDGE*id+UP,EDGE*is+UP]+1
    s%inverse_area=dom%areas%elts(id+1)%hex_inv
    if (any(s%edge<1).or.any(s%edge>size(native_mass(mass_plan_domain)%flux))) &
         error stop "native mass divergence address is invalid"
    end associate
  end subroutine

  subroutine compile_mass_restriction(dom,p_chd,ip,jp,ic,jc,k,op,dp_,oc,dc)
    type(Domain), intent(inout) :: dom
    integer,intent(in)::p_chd,ip,jp,ic,jc,k,op(N_BDRY+1),dp_(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
    integer :: id,t(20),e,x,y,center,mz,pz,id_mp,id_pp,id_pm,id_mm,weights(2,4),n,m,nnode
    real(dp) :: a(2),o(4)
    if (ic>=PATCH_SIZE.or.jc>=PATCH_SIZE) return
    id=idx(ip,jp,op,dp_)+1
    if (.not. any(dom%mask_e%elts(EDGE*(id-1)+1:EDGE*id)>=RESTRCT)) return
    mass_plan_count=mass_plan_count+1
    if (mass_plan_pass==1) return
    nnode=grid(mass_plan_domain)%node%length
    associate(s=>native_mass(mass_plan_domain)%level(mass_plan_level)%restriction(mass_plan_count))
    s%target=id
    s%edge=dom%mask_e%elts(EDGE*(id-1)+1:EDGE*id)>=RESTRCT
    weights(:,1)=[idx(ic+1,jc-2,oc,dc),idx(ic+1,jc-1,oc,dc)]+1
    weights(:,2)=[idx(ic,jc,oc,dc),idx(ic,jc+1,oc,dc)]+1
    weights(:,3)=[idx(ic+1,jc,oc,dc),idx(ic+1,jc+1,oc,dc)]+1
    weights(:,4)=[idx(ic-2,jc,oc,dc),idx(ic-2,jc+1,oc,dc)]+1
    do n=1,4
       do m=1,2
          s%weight(:,m,n)=dom%R_F_wgt%elts(weights(m,n))%enc
       end do
    end do
    call get_indices(dom,ic+1,jc,RT,oc,dc,t)
    s%small(:,1,1)=t([WPM,UZM,VMM]+1)+1
    s%small(:,2,1)=t([VPM,WMMM,UMZ]+1)+1
    s%small(:,3,1)=t([UPZ,VPMM,WMM]+1)+1
    s%small(:,1,2)=t([WMP,UZP,VPP]+1)+1
    s%small(:,2,2)=t([VMP,WPPP,UPZ]+1)+1
    s%small(:,3,2)=t([UMZ,VMPP,WPP]+1)+1
    call get_indices(dom,ic,jc+1,UP,oc,dc,t)
    s%small(:,1,3)=t([UZM,VMM,WPM]+1)+1
    s%small(:,2,3)=t([WMMM,UMZ,VPM]+1)+1
    s%small(:,3,3)=t([VPMM,WMM,UPZ]+1)+1
    s%small(:,1,4)=t([UZP,VPP,WMP]+1)+1
    s%small(:,2,4)=t([WPPP,UPZ,VMP]+1)+1
    s%small(:,3,4)=t([VMPP,WPP,UMZ]+1)+1
    if (any(s%small<1).or.any(s%small>EDGE*nnode)) error stop "native mass small-flux address is invalid"
    do e=RT,UP
       if (.not.s%edge(e+1)) cycle
       x=ic
       y=jc
       if (e/=UP) x=x+1
       if (e/=RT) y=y+1
       center=idx(x,y,oc,dc)+1
       call get_indices(dom,x,y,e,oc,dc,t)
       s%partial_flux(:,e+1)=t([UPZ,UMZ,VMM,WMP,WPM,VPP]+1)+1
       s%partial_node(:,e+1)=t([PP,MM,MP,PM]+1)+1
       a=dom%overl_areas%elts(center)%a(1:2)
       o(1:2)=dom%overl_areas%elts(center)%split
       o(3:4)=dom%overl_areas%elts(center)%a(3:4)-o(1:2)
       a(1)=a(1)+o(1)+o(4)
       a(2)=a(2)+o(2)+o(3)
       s%area(:,e+1)=a/sum(a)
       o(1)=dom%overl_areas%elts(t(PP+1)+1)%split(1)
       o(2)=dom%overl_areas%elts(t(MM+1)+1)%split(2)
       o(3)=dom%overl_areas%elts(t(MP+1)+1)%a(3)-dom%overl_areas%elts(t(MP+1)+1)%split(1)
       o(4)=dom%overl_areas%elts(t(PM+1)+1)%a(4)-dom%overl_areas%elts(t(PM+1)+1)%split(2)
       s%overlap(:,e+1)=o
       mz=idx2(x,y,nghb_pt(:,hex_s_offs(e+1)+2),oc,dc)+1
       pz=idx2(x,y,nghb_pt(:,hex_s_offs(e+1)+5),oc,dc)+1
       id_mp=idx2(x,y,nghb_pt(:,hex_s_offs(e+1)+1),oc,dc)+1
       id_pp=idx2(x,y,nghb_pt(:,hex_s_offs(e+1)+6),oc,dc)+1
       id_pm=idx2(x,y,nghb_pt(:,hex_s_offs(e+1)+4),oc,dc)+1
       id_mm=idx2(x,y,nghb_pt(:,hex_s_offs(e+1)+3),oc,dc)+1
       s%coarse_node(:,e+1)=[mz,pz,idx2(x,y,bfly_no2(:,3,e+1),oc,dc)+1, &
            idx2(x,y,bfly_no2(:,2,e+1),oc,dc)+1,idx2(x,y,bfly_no2(:,4,e+1),oc,dc)+1, &
            idx2(x,y,bfly_no2(:,1,e+1),oc,dc)+1]
       s%coarse(1,e+1)= &
            dom%overl_areas%elts(center)%a(1)*dom%overl_areas%elts(center)%a(2)*dom%areas%elts(center)%hex_inv &
            +dom%overl_areas%elts(id_mp)%a(2)*dom%overl_areas%elts(id_mp)%a(3)*dom%areas%elts(id_mp)%hex_inv &
            +dom%overl_areas%elts(id_pp)%a(1)*dom%overl_areas%elts(id_pp)%a(3)*dom%areas%elts(id_pp)%hex_inv &
            +dom%overl_areas%elts(id_pm)%a(1)*dom%overl_areas%elts(id_pm)%a(4)*dom%areas%elts(id_pm)%hex_inv &
            +dom%overl_areas%elts(id_mm)%a(2)*dom%overl_areas%elts(id_mm)%a(4)*dom%areas%elts(id_mm)%hex_inv
       t(1:4)=[id_pp,id_pm,id_mp,id_mm]
       do n=1,4
          m=t(n)
          s%coarse(n+1,e+1)=dom%overl_areas%elts(m)%a(3)*dom%overl_areas%elts(m)%a(4)*dom%areas%elts(m)%hex_inv
       end do
       if (any(s%partial_flux(:,e+1)<1).or.any(s%partial_flux(:,e+1)>EDGE*nnode) .or. &
            any(s%partial_node(:,e+1)<1).or.any(s%partial_node(:,e+1)>nnode) .or. &
            any(s%coarse_node(:,e+1)<1).or.any(s%coarse_node(:,e+1)>nnode)) &
            error stop "native mass restriction address is invalid"
    end do
    end associate
  end subroutine

  subroutine compile_mass_boundary(pos,l)
    ! Compile the existing geometry-owner graph, including signed edges.
    ! Same-rank copies intentionally cover ALL levels in their original order:
    ! cp_bdry_inside has that contract even during a single-level exchange.
    integer,intent(in)::pos,l
    integer :: r,ds,dd,g,id,i,pass,ns,nr,nlocal,mult
    mult=1
    if (pos==AT_EDGE) mult=EDGE
    associate(p=>mass_boundary(pos,l))
    allocate(p%send_count(n_process),p%recv_count(n_process),p%send_displ(n_process),p%recv_displ(n_process))
    do pass=1,2
       ns=0
       nr=0
       nlocal=0
       do r=1,n_process
          p%send_displ(r)=ns
          p%recv_displ(r)=nr
          if (r==rank+1) cycle
          do ds=1,size(grid)
             do dd=1,n_domain(r)
                g=glo_id(r,dd)+1
                do i=1,grid(ds)%pack(pos,g)%length
                   id=grid(ds)%pack(pos,g)%elts(i)
                   if (grid(ds)%level%elts(id/mult+1)/=l) cycle
                   ns=ns+1
                   if(pass==2) p%source(:,ns)=[ds,id+1]
                end do
             end do
          end do
          do ds=1,n_domain(r)
             g=glo_id(r,ds)+1
             do dd=1,size(grid)
                do i=1,grid(dd)%unpk(pos,g)%length
                   id=grid(dd)%unpk(pos,g)%elts(i)
                   if (grid(dd)%level%elts(abs(id)/mult+1)/=l) cycle
                   nr=nr+1
                   if(pass==2) p%destination(:,nr)=[dd,abs(id)+1,merge(-1,1,id<0.and.pos==AT_EDGE)]
                end do
             end do
          end do
       end do
       do ds=1,size(grid)
          do dd=1,size(grid)
             g=glo_id(rank+1,dd)+1
             do i=1,grid(ds)%pack(pos,g)%length
                nlocal=nlocal+1
                if (pass==1) cycle
                id=grid(dd)%unpk(pos,glo_id(rank+1,ds)+1)%elts(i)
                p%local(:,nlocal)=[ds,grid(ds)%pack(pos,g)%elts(i)+1,dd,abs(id)+1, &
                     merge(-1,1,id<0.and.pos==AT_EDGE)]
             end do
          end do
       end do
       if(pass==1) then
          allocate(p%source(2,ns),p%destination(3,nr),p%local(5,nlocal),p%send_value(ns),p%recv_value(nr))
       end if
    end do
    do r=1,n_process-1
       p%send_count(r)=p%send_displ(r+1)-p%send_displ(r)
       p%recv_count(r)=p%recv_displ(r+1)-p%recv_displ(r)
    end do
    p%send_count(n_process)=ns-p%send_displ(n_process)
    p%recv_count(n_process)=nr-p%recv_displ(n_process)
    end associate
  end subroutine

  subroutine exchange_native_mass(pos,l,k)
    integer,intent(in)::pos,l,k
    integer::i,r,nreq,ierr,ds,dd,si,di
    type(MPI_Request)::requests(2*n_process)
    associate(p=>mass_boundary(pos,l))
    do i=1,size(p%source,2)
       ds=p%source(1,i)
       si=p%source(2,i)
       if(pos==AT_EDGE) then
          p%send_value(i)=native_mass(ds)%flux(si)
       else
          p%send_value(i)=native_mass(ds)%tendency(si,k)
       end if
    end do
    nreq=0
    do r=1,n_process
       if(p%recv_count(r)==0) cycle
       nreq=nreq+1
       call MPI_Irecv(p%recv_value(p%recv_displ(r)+1:),p%recv_count(r),MPI_DP,r-1,28471,comm,requests(nreq),ierr)
       if(ierr/=MPI_SUCCESS) error stop "native mass boundary receive failed"
    end do
    do r=1,n_process
       if(p%send_count(r)==0) cycle
       nreq=nreq+1
       call MPI_Isend(p%send_value(p%send_displ(r)+1:),p%send_count(r),MPI_DP,r-1,28471,comm,requests(nreq),ierr)
       if(ierr/=MPI_SUCCESS) error stop "native mass boundary send failed"
    end do
    do i=1,size(p%local,2)
       ds=p%local(1,i)
       si=p%local(2,i)
       dd=p%local(3,i)
       di=p%local(4,i)
       if(pos==AT_EDGE) then
          native_mass(dd)%flux(di)=p%local(5,i)*native_mass(ds)%flux(si)
       else
          native_mass(dd)%tendency(di,k)=native_mass(ds)%tendency(si,k)
       end if
    end do
    if(nreq>0) then
       call detail_enter(DP_MASS_WAIT)
       call MPI_Waitall(nreq,requests(1:nreq),MPI_STATUSES_IGNORE,ierr)
       call detail_leave(DP_MASS_WAIT)
       if(ierr/=MPI_SUCCESS) error stop "native mass boundary completion failed"
    end if
    do i=1,size(p%destination,2)
       dd=p%destination(1,i)
       di=p%destination(2,i)
       if(pos==AT_EDGE) then
          native_mass(dd)%flux(di)=p%destination(3,i)*p%recv_value(i)
       else
          native_mass(dd)%tendency(di,k)=p%recv_value(i)
       end if
    end do
    mass_work(4)=mass_work(4)+size(p%recv_value)+size(p%local,2)
    end associate
  end subroutine

  subroutine assert_native_mass_level(dq,k,l)
    type(Float_Field),intent(in)::dq(1:N_VARIABLE,1:zlevels)
    integer,intent(in)::k,l
    integer::d,s,id,e
    do d=1,size(grid)
       do s=1,size(native_mass(d)%level(l)%divergence)
          id=native_mass(d)%level(l)%divergence(s)%target
          if(transfer(native_mass(d)%tendency(id,k),0_int64)/=transfer(dq(S_MASS,k)%data(d)%elts(id),0_int64)) then
             write(6,*) "Native mass divergence differs: rank, owner, layer, level, node",rank,d,k,l,id
             write(6,*) native_mass(d)%tendency(id,k),dq(S_MASS,k)%data(d)%elts(id)
             error stop "native mass phase divergence differs bit-for-bit"
          end if
          do e=1,6
             if(.not.native_mass(d)%level(l)%divergence(s)%active) cycle
             id=native_mass(d)%level(l)%divergence(s)%edge(e)
             if(transfer(native_mass(d)%flux(id),0_int64)/=transfer(horiz_flux(S_MASS)%data(d)%elts(id),0_int64)) &
                  error stop "native mass phase flux differs bit-for-bit"
          end do
       end do
    end do
  end subroutine

  subroutine prepare_native_velocity_programs
    integer :: d,l,p,c,slot,pass,nnode,ngrad,i,j,s,id,n,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer(int64) :: generation
    generation=native_velocity_plan_generation()
    if (allocated(native_velocity)) then
       if (size(native_velocity)/=size(grid)) deallocate(native_velocity)
    end if
    if (.not. allocated(native_velocity)) allocate(native_velocity(size(grid)))
    do d=1,size(grid)
       nnode=grid(d)%node%length
       if (native_velocity(d)%generation==generation) then
          if (size(native_velocity(d)%source,2)/=nnode .or. size(native_velocity(d)%ready)/=zlevels) &
               error stop "native velocity layout changed without a new generation"
          ! A mask-only adaptation can change which source edges are direct
          ! without changing the address/ownership generation. Compare the
          ! actual predicates, not a collision-prone checksum or dimensions.
          if (all(native_velocity(d)%active_mask .eqv. (grid(d)%mask_n%elts(1:nnode)>=TRSK)) .and. &
               all(native_velocity(d)%restriction_mask .eqv. (grid(d)%mask_e%elts(1:EDGE*nnode)>=ADJZONE))) cycle
       end if
       native_velocity(d)=Native_Velocity_Workspace()
       associate(work=>native_velocity(d))
       allocate(work%level(level_start:level_end),work%physics_address(3,nnode))
       allocate(work%source(3,nnode),work%physics(3,nnode),work%physics_ready(nnode))
       allocate(work%tendency(EDGE*nnode,zlevels),work%rk(EDGE*nnode,zlevels),work%ready(zlevels))
       allocate(work%rho(nnode),work%rho_theta(nnode))
       allocate(work%active_mask(nnode),work%restriction_mask(EDGE*nnode))
       work%active_mask=grid(d)%mask_n%elts(1:nnode)>=TRSK
       work%restriction_mask=grid(d)%mask_e%elts(1:EDGE*nnode)>=ADJZONE
       work%physics_address=0
       work%physics=0.0_dp
       work%ready=.false.
       velocity_plan_domain=d
       do l=level_end,level_start,-1
          velocity_plan_level=l
          do pass=1,2
             velocity_plan_pass=pass
             velocity_plan_action=0
             velocity_plan_direct=0
             velocity_plan_restriction=0
             do slot=1,grid(d)%lev(l)%length
                p=grid(d)%lev(l)%elts(slot)
                velocity_plan_patch=p
                if (l==level_end) then
                   call apply_onescale_to_patch(compile_velocity_direct,grid(d),p,0,0,0)
                else
                   ! Repeated absent-child evaluations have identical inputs.
                   ! Compile one initial direct pass, then preserve every
                   ! mixed-mask direct/restrict operation in traversal order.
                   do c=1,N_CHDRN
                      if (grid(d)%patch%elts(p+1)%children(c)/=0) cycle
                      call apply_onescale_to_patch(compile_velocity_direct,grid(d),p,0,0,0)
                      exit
                   end do
                   call apply_interscale_to_patch(compile_velocity_restriction,grid(d),p,0,0,0)
                end if
             end do
             if (pass==1) then
                allocate(work%level(l)%action(velocity_plan_action),work%level(l)%operand(velocity_plan_action))
                allocate(work%level(l)%direct(velocity_plan_direct),work%level(l)%restriction(velocity_plan_restriction))
             end if
          end do
          call validate_velocity_program(work%level(l),nnode)
       end do
       ngrad=0
       do p=3,grid(d)%patch%length
          if (.not. grid(d)%patch%elts(p)%deleted) ngrad=ngrad+PATCH_SIZE**2
       end do
       allocate(work%gradient(ngrad))
       s=0
       do p=3,grid(d)%patch%length
          if (grid(d)%patch%elts(p)%deleted) cycle
          call get_offs_Domain(grid(d),p-1,offs,dims)
          do j=0,PATCH_SIZE-1
             do i=0,PATCH_SIZE-1
                s=s+1
                id=idx(i,j,offs,dims)+1
                work%gradient(s)%node(1)=id
                work%gradient(s)%active=grid(d)%mask_n%elts(id)>=TRSK
                if (.not. work%gradient(s)%active) cycle
                do n=2,4
                   work%gradient(s)%node(n)=idx(i+velocity_dx(n),j+velocity_dy(n),offs,dims)+1
                end do
                work%gradient(s)%length=grid(d)%len%elts(EDGE*(id-1)+1:EDGE*id)
                if (any(work%gradient(s)%node<1) .or. any(work%gradient(s)%node>nnode)) &
                     error stop "native velocity gradient dependency is invalid"
                if (any(work%gradient(s)%length<=0.0_dp)) error stop "native velocity gradient metric is invalid"
             end do
          end do
       end do
       work%generation=generation
       native_velocity_work(1)=native_velocity_work(1)+1_int64
       end associate
    end do
    velocity_plan_domain=0
  end subroutine prepare_native_velocity_programs

  subroutine compile_velocity_direct(dom,i,j,k,offs,dims)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer :: id,n,d,l,s
    real(dp) :: parts(6,4),inverse_area(4)
    velocity_plan_action=velocity_plan_action+1
    velocity_plan_direct=velocity_plan_direct+1
    if (velocity_plan_pass==1) return
    d=velocity_plan_domain
    l=velocity_plan_level
    s=velocity_plan_direct
    associate(program=>native_velocity(d)%level(l))
    program%action(velocity_plan_action)=VELOCITY_DIRECT
    program%operand(velocity_plan_action)=s
    id=idx(i,j,offs,dims)+1
    program%direct(s)%node(1)=id
    program%direct(s)%active=dom%mask_n%elts(id)>=TRSK
    if (.not. program%direct(s)%active) return
    native_velocity(d)%physics_address(:,id)=[velocity_plan_patch,i,j]
    do n=1,9
       program%direct(s)%node(n)=idx(i+velocity_dx(n),j+velocity_dy(n),offs,dims)+1
    end do
    do n=1,4
       parts(:,n)=dom%areas%elts(program%direct(s)%node(n))%part
       inverse_area(n)=dom%areas%elts(program%direct(s)%node(n))%hex_inv
    end do
    program%direct(s)%weights=velocity_weights(parts,inverse_area)
    program%direct(s)%length=dom%len%elts(EDGE*(id-1)+1:EDGE*id)
    end associate
  end subroutine compile_velocity_direct

  subroutine compile_velocity_restriction(dom,ip,jp,ic,jc,k,op,dp_,oc,dc)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: ip,jp,ic,jc,k,op(N_BDRY+1),dp_(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
    integer :: child,d,l,s
    logical :: restricted(EDGE)
    child=idx(ic,jc,oc,dc)+1
    restricted=dom%mask_e%elts(EDGE*(child-1)+1:EDGE*child)>=ADJZONE
    if (.not. all(restricted)) call compile_velocity_direct(dom,ip,jp,k,op,dp_)
    velocity_plan_action=velocity_plan_action+1
    velocity_plan_restriction=velocity_plan_restriction+1
    if (velocity_plan_pass==1) return
    d=velocity_plan_domain
    l=velocity_plan_level
    s=velocity_plan_restriction
    associate(program=>native_velocity(d)%level(l))
    program%action(velocity_plan_action)=VELOCITY_RESTRICT
    program%operand(velocity_plan_action)=s
    program%restriction(s)%target=idx(ip,jp,op,dp_)+1
    program%restriction(s)%child=child
    program%restriction(s)%neighbor=[idx(ic+1,jc,oc,dc),idx(ic+1,jc+1,oc,dc),idx(ic,jc+1,oc,dc)]+1
    program%restriction(s)%edge=restricted
    end associate
  end subroutine compile_velocity_restriction

  subroutine compute_native_velocity_source(q,k,l)
    type(Float_Field), target, intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    integer, intent(in) :: k,l
    integer :: d,s,id,p,last_p,offs(N_BDRY+1),dims(2,N_BDRY+1),nnode
    real(dp), pointer :: flux_view(:,:),pv_view(:,:)
    do d=1,size(grid)
       mass=>q(S_MASS,k)%data(d)%elts
       velo=>q(S_VELO,k)%data(d)%elts
       mean_m=>sol_mean(S_MASS,k)%data(d)%elts
       h_mflux=>native_mass(d)%flux
       qe=>grid(d)%qe%elts
       ke=>grid(d)%ke%elts
       vort=>grid(d)%vort%elts
       if (Laplace_divu==2) then
          divu=>Laplacian_vector(S_DIVU)%data(d)%elts
       else
          divu=>grid(d)%divu%elts
       end if
       associate(work=>native_velocity(d))
       nnode=size(work%source,2)
       work%physics_ready=.false.
       last_p=-1
       do s=1,size(work%level(l)%direct)
          if (.not. work%level(l)%direct(s)%active) cycle
          id=work%level(l)%direct(s)%node(1)
          if (work%physics_ready(id)) cycle
          p=work%physics_address(1,id)
          if (p/=last_p) call get_offs_Domain(grid(d),p,offs,dims)
          last_p=p
          work%physics(:,id)=physics_velo_source(grid(d),work%physics_address(2,id),work%physics_address(3,id),k,offs,dims)
          work%physics_ready(id)=.true.
          native_velocity_work(5)=native_velocity_work(5)+1_int64
       end do
       ! Consume the native mass flux at this exact restriction phase. No
       ! Domain mass flux, velocity source or gradient is a production input.
       flux_view(1:EDGE,1:nnode)=>native_mass(d)%flux
       pv_view(1:EDGE,1:nnode)=>grid(d)%qe%elts
       call execute_velocity_sources(work%level(l),flux_view,pv_view,work%physics,work%source)
       native_velocity_work(2)=native_velocity_work(2)+int(size(work%level(l)%direct),int64)
       do s=1,size(work%level(l)%restriction)
          native_velocity_work(3)=native_velocity_work(3)+int(count(work%level(l)%restriction(s)%edge),int64)
       end do
       end associate
       nullify(mass,velo,mean_m,h_mflux,qe,ke,vort,divu,flux_view,pv_view)
    end do
  end subroutine compute_native_velocity_source

  subroutine compute_native_velocity_gradient(q,k)
    type(Float_Field), target, intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    integer, intent(in) :: k
    integer :: d,nnode
    real(dp), pointer :: tendency_view(:,:)
    do d=1,size(grid)
       associate(work=>native_velocity(d))
       nnode=size(work%source,2)
       work%rho=sol_mean(S_MASS,k)%data(d)%elts(1:nnode)+q(S_MASS,k)%data(d)%elts(1:nnode)
       work%rho_theta=sol_mean(S_TEMP,k)%data(d)%elts(1:nnode)+q(S_TEMP,k)%data(d)%elts(1:nnode)
       tendency_view(1:EDGE,1:nnode)=>native_velocity(d)%tendency(:,k)
       call execute_velocity_gradients(work%gradient,grid(d)%bernoulli%elts(1:nnode),exner_fun(k)%data(d)%elts(1:nnode), &
            work%rho,work%rho_theta,work%source,tendency_view)
       work%ready(k)=.true.
       native_velocity_work(4)=native_velocity_work(4)+int(size(work%gradient),int64)
       end associate
    end do
    nullify(tendency_view)
  end subroutine compute_native_velocity_gradient

  
  subroutine basic_operators (q, dq, k, l, mass_only_compatibility)
    ! Evaluates basic operators on grid level l and computes/restricts Bernoulli, Exner and fluxes
    
    implicit none
    
    type(Float_Field), target, intent(inout) :: q(1:N_VARIABLE,1:zlevels), dq(1:N_VARIABLE,1:zlevels)
    integer,                   intent(in)    :: k, l
    logical,                   intent(in)    :: mass_only_compatibility

    integer :: d, j, v, scalar_last
    logical :: capture_scalar_physics
    real(dp) :: scalar_physics(EDGE,PATCH_SIZE**2,scalars(1):scalars(2))
    real(dp) :: mass_start

    call detail_enter(DP_BASIC)
    capture_scalar_physics = block_scalar_capture_active()

    do d = 1, size(grid)
       mass      => q(S_MASS,k)%data(d)%elts
       temp      => q(S_TEMP,k)%data(d)%elts
       velo      => q(S_VELO,k)%data(d)%elts
       mean_m    => sol_mean(S_MASS,k)%data(d)%elts
       mean_t    => sol_mean(S_TEMP,k)%data(d)%elts
       exner     => exner_fun(k)%data(d)%elts
       bernoulli => grid(d)%bernoulli%elts
       ke        => grid(d)%ke%elts
       vort      => grid(d)%vort%elts
       qe        => grid(d)%qe%elts
       
       ! Compute horizontal fluxes, potential vorticity (qe), Bernoulli, Exner (incompressible case) etc
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
          if (capture_scalar_physics) then
             call step1(dq,q,grid(d),grid(d)%lev(l)%elts(j),k,0, &
                  scalar_physics,mass_only_compatibility)
             call capture_block_scalar_physics_patch( &
                  d,grid(d)%lev(l)%elts(j),k,scalar_physics)
          else
             call step1(dq,q,grid(d),grid(d)%lev(l)%elts(j),k,0)
          end if
       end do
       call apply_to_penta_d (post_step1, grid(d), l, z_null)
       nullify (mass, velo, temp, mean_m, mean_t, ke, qe, vort)

       nullify (bernoulli, exner)
    end do

    ! Retain the direct step1 flux before a coarse level is overwritten by
    ! fine-to-coarse flux restriction.
    call capture_block_scalar_divergence_level( &
         q,physics_scalar_flux,0,k,l,.true.)

    ! Temperature restriction and its edge/node exchanges are block-native.

    ! Compute or restrict Bernoulli, Exner and fluxes only after every local
    ! Domain direct-flux shadow has been captured.
    if (l < level_end) then
       do d = 1,size(grid)
          scalar => grid(d)%bernoulli%elts
          call cpt_or_restr_scalar (grid(d), l)
          nullify (scalar)

          scalar => exner_fun(k)%data(d)%elts
          call cpt_or_restr_scalar (grid(d), l)
          nullify (scalar)

          scalar_last = scalars(2)
          if (mass_only_compatibility) scalar_last = S_MASS
          do v = scalars(1),scalar_last
             if (mass_transaction.and..not.mass_oracle) cycle
             dscalar => dq(v,k)%data(d)%elts
             h_flux  => horiz_flux(v)%data(d)%elts
             call cpt_or_restr_flux (grid(d), l)
             nullify (dscalar, h_flux)
          end do
       end do
    end if
    if (mass_transaction.and..not.mass_oracle) then
       ! Production mass storage and boundary processing are entirely native.
    else if (mass_only_compatibility) then
       horiz_flux(S_MASS)%bdry_uptodate = .false.
       if (level_start /= level_end) &
            call update_bdry(horiz_flux(S_MASS:S_MASS),l,1068)
    else
       horiz_flux%bdry_uptodate = .false.
       if (level_start /= level_end) call update_bdry(horiz_flux,l,968)
    end if

    if (mass_transaction) then
       mass_start=parallel_block_profile_begin(BLOCK_PROFILE_NATIVE_MASS)
       do d=1,size(native_mass)
          call execute_mass_restriction(native_mass(d)%level(l)%restriction, &
               native_mass(d)%flux,native_mass(d)%tendency(:,k))
       end do
       if (level_start/=level_end) call exchange_native_mass(AT_EDGE,l,k)
       call parallel_block_profile_end(BLOCK_PROFILE_NATIVE_MASS,mass_start)
    end if

    if (Laplace_rotu == 2) call cal_Laplacian_vector_rot (l) ! requires vorticity
    call detail_leave(DP_BASIC)

  end subroutine basic_operators


  subroutine cal_scalar_trend_compatibility (q,dq,k,l)
    ! Compute both compatibility scalar trends in one topology traversal.
    ! The complete trend_ml oracle continues to use cal_scalar_trend and its
    ! independent scalar_trend passes.

    implicit none

    type(Float_Field), target, intent(inout) :: &
         q(1:N_VARIABLE,1:zlevels),dq(1:N_VARIABLE,1:zlevels)
    integer, intent(in) :: k
    integer, intent(in) :: l

    integer :: d
    integer :: j

    call update_bdry(horiz_flux,l,1169)
    do d = 1,size(grid)
       do j = 1,grid(d)%lev(l)%length
          call scalar_trend_pair( &
               grid(d),grid(d)%lev(l)%elts(j), &
               horiz_flux(S_MASS)%data(d)%elts, &
               horiz_flux(S_TEMP)%data(d)%elts, &
               dq(S_MASS,k)%data(d)%elts, &
               dq(S_TEMP,k)%data(d)%elts)
       end do
    end do
    call capture_block_scalar_divergence_level( &
         q,physics_scalar_flux,0,k,l)
    call capture_block_scalar_divergence_level( &
         q,physics_scalar_flux,0,k,l, &
         domain_tendency=dq,dscalar_only=.true.)
    dq(S_MASS:S_TEMP,k)%bdry_uptodate = .false.

  end subroutine cal_scalar_trend_compatibility

  
  subroutine cal_scalar_trend (q, dq, k, l, mass_only_compatibility)
    ! Evaluate scalar trends at level l
    
    implicit none
    
    type(Float_Field), target, intent(inout) :: q(1:N_VARIABLE,1:zlevels), dq(1:N_VARIABLE,1:zlevels)
    integer,                   intent(in)    :: k, l
    logical,                   intent(in)    :: mass_only_compatibility

    integer :: d, j, v, scalar_last

    if (mass_only_compatibility) then
       call update_bdry(horiz_flux(S_MASS:S_MASS),l,1069)
       scalar_last = S_MASS
    else
       call update_bdry(horiz_flux,l,969)
       scalar_last = scalars(2)
    end if
    
    do d = 1, size(grid)
       do v = scalars(1),scalar_last
          dscalar => dq(v,k)%data(d)%elts
          h_flux  => horiz_flux(v)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (scalar_trend, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
    end do
    if (.not. mass_only_compatibility) then
       call capture_block_scalar_divergence_level( &
            q,physics_scalar_flux,0,k,l)
       call capture_block_scalar_divergence_level( &
            q,physics_scalar_flux,0,k,l, &
            domain_tendency=dq,dscalar_only=.true.)
       dq(S_MASS:S_TEMP,k)%bdry_uptodate = .false.
    else
       dq(S_MASS,k)%bdry_uptodate = .false.
    end if
  end subroutine cal_scalar_trend

  
  subroutine velocity_trend_source (q, dq, k, l, measure_source)
    ! Evaluate source part of velocity trends at level l
    
    implicit none
    
    type(Float_Field), target, intent(inout) :: q(1:N_VARIABLE,1:zlevels), dq(1:N_VARIABLE,1:zlevels)
    integer,                    intent(in)   :: k, l
    logical, optional,          intent(in)   :: measure_source

    integer :: d, j

    logical :: retain_measurement

    retain_measurement = .false.
    if (present(measure_source)) retain_measurement = measure_source

    if (native_velocity_transaction) then
       if (.not. native_velocity_oracle) error stop "production called Domain velocity source"
       native_velocity_work(8)=native_velocity_work(8)+1_int64
    end if

    u_source => du_source

    do d = 1, size(grid)
       mass    => q(S_MASS,k)%data(d)%elts
       velo    => q(S_VELO,k)%data(d)%elts
       mean_m  => sol_mean(S_MASS,k)%data(d)%elts
       dvelo   => dq(S_VELO,k)%data(d)%elts
       h_mflux => horiz_flux(S_MASS)%data(d)%elts
       ke      => grid(d)%ke%elts
       qe      => grid(d)%qe%elts
       vort    => grid(d)%vort%elts
       
       if (Laplace_divu == 2) then
          divu => Laplacian_vector(S_DIVU)%data(d)%elts
       else
          divu => grid(d)%divu%elts
       end if

       if (l < level_end) then
          call cpt_or_restr_u_source (grid(d), k, l)
       else
          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (u_source, grid(d), grid(d)%lev(level_end)%elts(j), k, 0, 0)
          end do
       end if

       if (retain_measurement) then
          call capture_velocity_source_measurement(d,k,l)
          call validate_native_velocity_source(d,k,l)
       end if

       nullify (mass, velo, mean_m, dvelo, h_mflux, divu, ke, qe, vort)
    end do
    dq(S_VELO,k)%bdry_uptodate = .false.

    nullify (u_source)
  end subroutine velocity_trend_source

  
  subroutine velocity_trend_grad (q, dq, k, measure_source)
    ! Evaluate complete velocity trend by adding gradient terms to previously calculated source terms on entire grid
    
    implicit none
    
    type(Float_Field), target, intent(inout) :: q(1:N_VARIABLE,1:zlevels), dq(1:N_VARIABLE,1:zlevels)
    integer,                   intent(in)    :: k
    logical, optional,         intent(in)    :: measure_source

    integer :: d, p

    logical :: validate_measurement

    validate_measurement = .false.
    if (present(measure_source)) validate_measurement = measure_source
    if (native_velocity_transaction) then
       if (.not. native_velocity_oracle) error stop "production called Domain velocity gradient"
       native_velocity_work(9)=native_velocity_work(9)+1_int64
    end if

    do d = 1, size(grid)
       mass      => q(S_MASS,k)%data(d)%elts
       temp      => q(S_TEMP,k)%data(d)%elts
       mean_m    => sol_mean(S_MASS,k)%data(d)%elts
       mean_t    => sol_mean(S_TEMP,k)%data(d)%elts
       dvelo     => dq(S_VELO,k)%data(d)%elts
       exner     => exner_fun(k)%data(d)%elts
       bernoulli => grid(d)%bernoulli%elts
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (du_grad, grid(d), p-1, k, 0, 0)
       end do
       if (validate_measurement) &
            call validate_velocity_source_measurement(d,k)
       nullify (mass, temp, mean_m, mean_t, dvelo, exner, bernoulli)
    end do
    dq(S_VELO,k)%bdry_uptodate = .false.
  end subroutine velocity_trend_grad

  subroutine capture_velocity_source_measurement (d,k,l)

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: k
    integer, intent(in) :: l

    integer :: dims(2,N_BDRY+1)
    integer :: e
    integer :: i
    integer :: id
    integer :: j
    integer :: p
    integer :: patch_slot
    integer :: pos
    integer :: offs(N_BDRY+1)
    real(dp) :: physics_value(EDGE)
    real(dp) :: qperp_value(EDGE)

    do patch_slot = 1,grid(d)%lev(l)%length
       p = grid(d)%lev(l)%elts(patch_slot)
       call get_offs_Domain(grid(d),p,offs,dims)
       do j = 0,PATCH_SIZE-1
          do i = 0,PATCH_SIZE-1
             id = idx(i,j,offs,dims)
             qperp_value = Qperp(grid(d),i,j,z_null,offs,dims)
             physics_value = physics_velo_source( &
                  grid(d),i,j,k,offs,dims)
             do e = RT,UP
                pos = EDGE*id+e+1
                velocity_source_measurement(d)%qperp(pos) = &
                     qperp_value(e+1)
                velocity_source_measurement(d)%physics(pos) = &
                     physics_value(e+1)
                velocity_source_measurement(d)%edge_length(pos) = &
                     grid(d)%len%elts(pos)
                velocity_source_measurement(d)%integrated_source(pos) = &
                     dvelo(pos)
                velocity_source_measurement(d)%active(pos) = &
                     grid(d)%mask_n%elts(id+1) >= TRSK
                velocity_source_measurement(d)%covered(pos) = .true.
                velocity_source_measurement(d)%direct(pos) = &
                     l == level_end
             end do
          end do
       end do
    end do

  end subroutine capture_velocity_source_measurement


  function native_qperp_sample(dom,i,j,offs,dims) result(value)
    ! Oracle adapter only. The numeric module cannot access Domain fields.
    ! These primitives are sampled at source time, not after later levels
    ! have overwritten/restricted them.
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,offs(N_BDRY+1),dims(2,N_BDRY+1)
    real(dp) :: value(EDGE), flux(3,9), pv(3,9), parts(6,4), inverse_area(4)
    integer :: n,id
    do n = 1,9
       id = idx(i+velocity_dx(n),j+velocity_dy(n),offs,dims)
       flux(:,n) = h_mflux(EDGE*id+1:EDGE*id+EDGE)
       pv(:,n) = qe(EDGE*id+1:EDGE*id+EDGE)
    end do
    do n = 1,4
       id = idx(i+velocity_dx(n),j+velocity_dy(n),offs,dims)
       parts(:,n) = dom%areas%elts(id+1)%part
       inverse_area(n) = dom%areas%elts(id+1)%hex_inv
    end do
    value = velocity_qperp(flux,pv,velocity_weights(parts,inverse_area))
  end function native_qperp_sample


  subroutine native_velocity_direct_sample(dom,i,j,k,offs,dims)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer :: id,first,last,d
    real(dp) :: qperp_value(EDGE)
    logical :: active
    d = native_velocity_measurement_domain
    if (d < 1) error stop "native velocity oracle has no traversal context"
    id = idx(i,j,offs,dims)
    first = EDGE*id+1
    last = first+EDGE-1
    if (.not. all(velocity_source_measurement(d)%covered(first:last))) &
         error stop "native velocity oracle direct primitive coverage is incomplete"
    active = dom%mask_n%elts(id+1) >= TRSK
    qperp_value = 0.0_dp
    if (active) then
       qperp_value = native_qperp_sample(dom,i,j,offs,dims)
       call assert_native_velocity_components(qperp_value,velocity_source_measurement(d)%qperp(first:last), &
            "Qperp",d,k,i,j,id)
    end if
    velocity_source_measurement(d)%native_source(first:last) = velocity_source(qperp_value, &
         velocity_source_measurement(d)%physics(first:last),dom%len%elts(first:last),active)
  end subroutine native_velocity_direct_sample


  subroutine native_velocity_restriction_sample(dom,ip,jp,ic,jc,k,op,dp_,oc,dc)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: ip,jp,ic,jc,k,op(N_BDRY+1),dp_(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
    integer :: parent,child,neighbor(EDGE),e,d
    real(dp) :: direct(EDGE), child_source(EDGE), neighbor_source(EDGE)
    logical :: restricted(EDGE)
    d = native_velocity_measurement_domain
    if (d < 1) error stop "native velocity oracle has no restriction context"
    parent = EDGE*idx(ip,jp,op,dp_)
    child = EDGE*idx(ic,jc,oc,dc)
    restricted = dom%mask_e%elts(child+1:child+EDGE) >= ADJZONE
    if (.not. all(restricted)) call native_velocity_direct_sample(dom,ip,jp,k,op,dp_)
    neighbor = EDGE*[idx(ic+1,jc,oc,dc),idx(ic+1,jc+1,oc,dc),idx(ic,jc+1,oc,dc)]
    direct = 0.0_dp
    child_source = 0.0_dp
    neighbor_source = 0.0_dp
    do e = 1,EDGE
       if (restricted(e)) then
          child_source(e) = velocity_source_measurement(d)%native_source(child+e)
          neighbor_source(e) = velocity_source_measurement(d)%native_source(neighbor(e)+e)
       else
          direct(e) = velocity_source_measurement(d)%native_source(parent+e)
       end if
    end do
    velocity_source_measurement(d)%native_source(parent+1:parent+EDGE) = &
         velocity_restrict_source(direct,child_source,neighbor_source,restricted)
  end subroutine native_velocity_restriction_sample


  subroutine validate_native_velocity_source(d,k,l)
    ! Reproduce the source schedule in separate storage, including absent
    ! children and mixed masks. The only reads of Domain dvelo are assertions.
    integer, intent(in) :: d,k,l
    integer :: p,slot,child,i,j,id,offs(N_BDRY+1),dims(2,N_BDRY+1)
    if (native_velocity_measurement_domain /= 0) error stop "nested native velocity oracle traversal"
    native_velocity_measurement_domain = d
    do slot = 1,grid(d)%lev(l)%length
       p = grid(d)%lev(l)%elts(slot)
       if (l == level_end) then
          call apply_onescale_to_patch(native_velocity_direct_sample,grid(d),p,k,0,0)
       else
          do child = 1,N_CHDRN
             if (grid(d)%patch%elts(p+1)%children(child) == 0) &
                  call apply_onescale_to_patch(native_velocity_direct_sample,grid(d),p,k,0,0)
          end do
          call apply_interscale_to_patch(native_velocity_restriction_sample,grid(d),p,k,0,0)
       end if
       call get_offs_Domain(grid(d),p,offs,dims)
       do j = 0,PATCH_SIZE-1
          do i = 0,PATCH_SIZE-1
             id = EDGE*idx(i,j,offs,dims)
             call assert_native_velocity_components(velocity_source_measurement(d)%native_source(id+1:id+EDGE), &
                  dvelo(id+1:id+EDGE),"restricted source",d,k,i,j,id/EDGE)
             call assert_native_velocity_components(native_velocity(d)%source(:,id/EDGE+1), &
                  dvelo(id+1:id+EDGE),"compiled production source",d,k,i,j,id/EDGE)
          end do
       end do
    end do
    native_velocity_measurement_domain = 0
  end subroutine validate_native_velocity_source


  subroutine assert_native_velocity_components(value,reference,label,d,k,i,j,node)
    real(dp), intent(in) :: value(EDGE),reference(EDGE)
    character(*), intent(in) :: label
    integer, intent(in) :: d,k,i,j,node
    if (all(ieee_is_finite(value)) .and. all(ieee_is_finite(reference))) then
       if (all(transfer(value,[0_int64],EDGE) == transfer(reference,[0_int64],EDGE))) return
    end if
    write(*,'(a,a)') "Stage 175 native velocity oracle mismatch: ",label
    write(*,'(a,4(i0,1x))') "local Domain, physical layer, i, j = ",d,k,i,j
    write(*,'(a,2(i0,1x))') "global Domain id, zero-based node = ",grid(d)%id,node
    write(*,'(a,3es25.16)') "native = ",value
    write(*,'(a,3es25.16)') "reference = ",reference
    error stop "native velocity component differs bit-for-bit"
  end subroutine assert_native_velocity_components


  subroutine validate_velocity_source_measurement (d,k)

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: k

    integer :: dims(2,N_BDRY+1)
    integer :: e
    integer :: i
    integer :: id
    integer :: id_e
    integer :: id_n
    integer :: id_ne
    integer :: j
    integer :: p
    integer :: pos
    integer :: offs(N_BDRY+1)
    real(dp) :: component_value
    real(dp) :: expected_value
    real(dp) :: grad_b(EDGE)
    real(dp) :: grad_e(EDGE)
    real(dp) :: rho(4)
    real(dp) :: rho_theta(4)
    real(dp) :: theta_edge(EDGE)
    real(dp) :: native_value(EDGE)
    integer :: node_ids(4)

    do p = 3,grid(d)%patch%length
       if (grid(d)%patch%elts(p)%deleted) cycle
       call get_offs_Domain(grid(d),p-1,offs,dims)
       do j = 0,PATCH_SIZE-1
          do i = 0,PATCH_SIZE-1
             id = idx(i,j,offs,dims)
             id_e = idx(i+1,j,offs,dims)
             id_ne = idx(i+1,j+1,offs,dims)
             id_n = idx(i,j+1,offs,dims)
             grad_b = gradi_e(bernoulli,grid(d),i,j,offs,dims)
             grad_e = gradi_e(exner,grid(d),i,j,offs,dims)
             rho = mean_m([id,id_e,id_ne,id_n]+1) + &
                  mass([id,id_e,id_ne,id_n]+1)
             rho_theta = mean_t([id,id_e,id_ne,id_n]+1) + &
                  temp([id,id_e,id_ne,id_n]+1)
             node_ids = [id,id_e,id_ne,id_n]+1
             if (all(velocity_source_measurement(d)%covered(EDGE*id+1:EDGE*id+EDGE))) then
                native_value = velocity_gradient(velocity_source_measurement(d)%native_source(EDGE*id+1:EDGE*id+EDGE), &
                     grid(d)%len%elts(EDGE*id+1:EDGE*id+EDGE),bernoulli(node_ids),exner(node_ids),rho,rho_theta, &
                     grid(d)%mask_n%elts(id+1) >= TRSK)
                call assert_native_velocity_components(native_value,dvelo(EDGE*id+1:EDGE*id+EDGE), &
                     "complete gradient",d,k,i,j,id)
                call assert_native_velocity_components(native_velocity(d)%tendency(EDGE*id+1:EDGE*id+EDGE,k), &
                     dvelo(EDGE*id+1:EDGE*id+EDGE),"compiled production gradient",d,k,i,j,id)
             end if
             theta_edge = [ &
                  0.5_dp*(rho_theta(1)/rho(1)+rho_theta(2)/rho(2)), &
                  0.5_dp*(rho_theta(1)/rho(1)+rho_theta(3)/rho(3)), &
                  0.5_dp*(rho_theta(1)/rho(1)+rho_theta(4)/rho(4))]
             do e = RT,UP
                pos = EDGE*id+e+1
                if (.not. velocity_source_measurement(d)%covered(pos)) &
                     cycle
                if (velocity_source_measurement(d)%active(pos)) then
                   if (velocity_source_measurement(d)% &
                        edge_length(pos) <= 0.0_dp) &
                        error stop "measured velocity edge length is nonpositive"
                   expected_value = velocity_source_measurement(d)% &
                        integrated_source(pos)/ &
                        velocity_source_measurement(d)%edge_length(pos) - &
                        grad_b(e+1) - theta_edge(e+1)*grad_e(e+1)
                   if (abs(dvelo(pos)-expected_value) > &
                        64.0_dp*epsilon(1.0_dp)*max( &
                        1.0_dp,abs(dvelo(pos)),abs(expected_value))) &
                        error stop "measured non-Exner velocity source differs"
                   if (velocity_source_measurement(d)%direct(pos)) then
                      component_value = &
                           -velocity_source_measurement(d)%qperp(pos) + &
                           velocity_source_measurement(d)%physics(pos)* &
                           velocity_source_measurement(d)%edge_length(pos)
                      if (abs(component_value- &
                           velocity_source_measurement(d)% &
                           integrated_source(pos)) > &
                           64.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
                           abs(component_value),abs( &
                           velocity_source_measurement(d)% &
                           integrated_source(pos)))) &
                           error stop "measured direct velocity source differs"
                   end if
                else if (abs(dvelo(pos)) > tiny(1.0_dp)) then
                   error stop "inactive measured velocity source is nonzero"
                end if
             end do
          end do
       end do
    end do

  end subroutine validate_velocity_source_measurement

  
  subroutine cal_Laplacian_scalars (q, k)
    ! Computes Laplacian of scalars q, div(grad q)
    
    implicit none
    
    type(Float_Field), target, intent(inout) :: q(1:N_VARIABLE,1:zlevels)
    integer,                   intent(in)    :: k
    
    integer :: d, j, l, v

    call update_bdry (q(scalars(1):scalars(2),k), NONE, 970)
    
    do l = level_end, level_start, -1
       ! Compute scalar fluxes
       do d = 1, size(grid)
          do v = scalars(1), scalars(2)
             scalar => q(v,k)%data(d)%elts
             h_flux => horiz_flux(v)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=1)
             end do
             nullify (scalar, h_flux)
          end do

          ! Compute or restrict fluxes
          if (l < level_end) then
             do v = scalars(1), scalars(2)
                dscalar => Laplacian_scalar(v)%data(d)%elts
                h_flux  => horiz_flux(v)%data(d)%elts
                call cpt_or_restr_flux (grid(d), l)  ! <= compute flux(l) using dscalar (l+1)
                nullify (dscalar, h_flux)
             end do
          end if
       end do
       horiz_flux%bdry_uptodate = .false.
       call update_bdry (horiz_flux, l, 971)

       do d = 1, size(grid)
          do v = scalars(1), scalars(2)
             dscalar => Laplacian_scalar(v)%data(d)%elts
             h_flux  => horiz_flux(v)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (dscalar, h_flux)
          end do
       end do
       Laplacian_scalar%bdry_uptodate = .false.
       call update_bdry (Laplacian_scalar, l, 972)
    end do
  end subroutine cal_Laplacian_scalars

  
  subroutine cal_Laplacian_vector_rot (l)
    ! Computes rot(rot(vorticity)) needed for second-order vector Laplacian
    
    implicit none
    
    integer, intent(in) :: l
    
    integer :: d, j

    ! Compute rot(vorticity)
    do d = 1, size(grid)
       vort      => grid(d)%vort%elts
       Laplacian => Laplacian_vector(S_ROTU)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_Laplacian_rotu, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
       end do
       nullify (vort, Laplacian)
    end do
    Laplacian_vector(S_ROTU)%bdry_uptodate = .false.
    call update_bdry (Laplacian_vector(S_ROTU), l, 973)
        
    ! Compute rot(rot(vorticity)) using previous result for rot(vorticity)
    !!! grid(d)%vort is now rot(rot(vorticity)), not vorticity !!!
    do d = 1, size(grid)
       velo => Laplacian_vector(S_ROTU)%data(d)%elts
       vort => grid(d)%vort%elts
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=8)
       end do
       call apply_to_penta_d (post_vort, grid(d), l, z_null)
       nullify (velo, vort)
    end do
  end subroutine cal_Laplacian_vector_rot

  
  subroutine cal_Laplacian_divu
    ! Computes Laplacian of divu, div(grad divu)
    
    implicit none
    
    integer :: d, j, l

    do l = level_end, level_start, -1
       ! Compute scalar fluxes
       do d = 1, size(grid)
          scalar => grid(d)%divu%elts
          h_flux => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=1)
          end do
          nullify (scalar, h_flux)

          ! Compute or restrict fluxes
          if (l < level_end) then
             dscalar => Laplacian_vector(S_DIVU)%data(d)%elts
             h_flux  => horiz_flux(S_MASS)%data(d)%elts
             call cpt_or_restr_flux (grid(d), l)  ! <= compute flux(l) using dscalar (l+1)
             nullify (dscalar, h_flux)
          end if
       end do
       horiz_flux(S_MASS)%bdry_uptodate = .false.
       call update_bdry (horiz_flux(S_MASS), l, 974)

       do d = 1, size(grid)
          dscalar => Laplacian_vector(S_DIVU)%data(d)%elts
          h_flux  => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
       Laplacian_vector(S_DIVU)%bdry_uptodate = .false.
       call update_bdry (Laplacian_vector(S_DIVU), l, 975)
    end do
  end subroutine cal_Laplacian_divu


  subroutine cal_Laplacian_rotu (dom, i, j, zlev, offs, dims)
    ! Curl of vorticity given at triangle circumcentres x_v, i.e. rotational part of vector Laplacian
    ! output is at edges x_e

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id, idS, idW

    id   = idx (i,   j,   offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)

    Laplacian(EDGE*id+RT+1) = - (vort(TRIAG*id+LORT+1) - vort(TRIAG*idS+UPLT+1)) / dom%pedlen%elts(EDGE*id+RT+1)

    if (dom%pedlen%elts(EDGE*id+DG+1) > eps (radius)) then
       Laplacian(EDGE*id+DG+1) = - (vort(TRIAG*id+LORT+1) - vort(TRIAG*id+UPLT+1)) / dom%pedlen%elts(EDGE*id+DG+1)
    else
       Laplacian(EDGE*id+DG+1) = 0.0_dp
    end if

    Laplacian(EDGE*id+UP+1) = - (vort(TRIAG*idW+LORT+1) - vort(TRIAG*id+UPLT+1)) / dom%pedlen%elts(EDGE*id+UP+1)
  end subroutine cal_Laplacian_rotu


  subroutine cpt_or_restr_scalar (dom, l)
    ! Restrict scalar if possible for grad(scalar) computation

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: l

    integer                     :: j, p_par, c, p_chd
    logical, dimension(N_CHDRN) :: restrict

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .false.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd > 0) restrict(c) = .true.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) then
             call apply_interscale_to_patch3 (scalar_cpt_restr, dom, p_par, c, z_null, 0, 1)
          end if
       end do
    end do
  end subroutine cpt_or_restr_scalar

  
  subroutine scalar_cpt_restr (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)

    integer :: id_par, id_chd

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx (i_par, j_par, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) >= RESTRCT) scalar(id_par+1) = scalar(id_chd+1)
  end subroutine scalar_cpt_restr

  
  subroutine cpt_or_restr_u_source (dom, zlev, l)
    ! Restrict velocity source if possible term u_source(velo)
    ! u_source is a pointer function
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: zlev, l

    integer :: c, j, p_par, p_chd

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd == 0) call apply_onescale_to_patch (u_source, dom, p_par, zlev, 0, 0)
       end do
       call apply_interscale_to_patch (u_source_cpt_restr, dom, dom%lev(l)%elts(j), zlev, 0, 0)
    end do
  end subroutine cpt_or_restr_u_source

  
  subroutine u_source_cpt_restr (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)
  
    integer :: id_par, id_chd, idE_chd, idNE_chd, idN_chd

    id_par = idx (i_par, j_par, offs_par, dims_par)

    id_chd   = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (minval(dom%mask_e%elts(EDGE*id_chd+RT+1:EDGE*id_chd+UP+1)) < ADJZONE) &
         call u_source (dom, i_par, j_par, zlev, offs_par, dims_par)

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) &
         dvelo(EDGE*id_par+RT+1) = dvelo(EDGE*id_chd+RT+1) + dvelo(EDGE*idE_chd+RT+1)

    if (dom%mask_e%elts(EDGE*id_chd+DG+1) >= ADJZONE) &
         dvelo(EDGE*id_par+DG+1) = dvelo(EDGE*id_chd+DG+1) + dvelo(EDGE*idNE_chd+DG+1)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) &
         dvelo(EDGE*id_par+UP+1) = dvelo(EDGE*id_chd+UP+1) + dvelo(EDGE*idN_chd+UP+1)
  end subroutine u_source_cpt_restr

  
  subroutine cpt_or_restr_flux (dom, l)
    ! Restrict flux if possible for dscalar = div(h_flux) computation
    ! requires dscalar = div(h_flux) in addition to h_flux
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: l

    integer                     :: j, p_par, c, p_chd
    logical, dimension(N_CHDRN) :: restrict

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .false.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd > 0) restrict(c) = .true.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) call apply_interscale_to_patch3 (flux_cpt_restr, dom, p_par, c, z_null, 0, 1)
       end do
    end do
  end subroutine cpt_or_restr_flux

  
  subroutine flux_cpt_restr (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute flux restriction by summing coarse, corrective and small fluxes

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)

    integer                :: id_par
    real(dp), dimension(4) :: sm_flux

    if (i_chd >= PATCH_SIZE .or. j_chd >= PATCH_SIZE) return

    id_par = idx (i_par, j_par, offs_par, dims_par)
    
    if (maxval(dom%mask_e%elts(EDGE*id_par+RT+1:EDGE*id_par+UP+1)) >= RESTRCT) &
         sm_flux = interp_flux (dom, i_chd, j_chd, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_par+RT+1) >= RESTRCT) h_flux(EDGE*id_par+RT+1) = &
         complete_coarse_flux (sm_flux, dom, i_par, j_par, i_chd, j_chd, RT, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_par+DG+1) >= RESTRCT) h_flux(EDGE*id_par+DG+1) = &
         complete_coarse_flux (sm_flux, dom, i_par, j_par, i_chd, j_chd, DG, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_par+UP+1) >= RESTRCT) h_flux(EDGE*id_par+UP+1) = &
         complete_coarse_flux (sm_flux, dom, i_par, j_par, i_chd, j_chd, UP, offs_chd, dims_chd)
  end subroutine flux_cpt_restr

  
  function interp_flux (dom, i, j, offs, dims) result(val)
    
    implicit none
    
    type(Domain), intent(in) :: dom
    integer,      intent(in) :: i, j
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)
    real(dp)                 :: val(4)
    
    integer                :: id, idE, idNE, idN, idSE, id2SE, id2W, idN2W
    integer, dimension(20) :: id_edges
    
    call get_indices (dom, i+1, j, RT, offs, dims, id_edges)

    id = idx (i, j, offs, dims)
    
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    
    idSE  = idx (i+1, j-1, offs, dims)
    id2SE = idx (i+1, j-2, offs, dims)
    id2W  = idx (i-2, j,   offs, dims)
    idN2W = idx (i-2, j+1, offs, dims)

    val(1) = - sum (h_flux(id_edges([WPM,UZM,VMM]+1)+1) * dom%R_F_wgt%elts(id2SE+1)%enc) &
         - sum ((h_flux(id_edges([VPM,WMMM,UMZ]+1)+1) - h_flux(id_edges([UPZ,VPMM,WMM]+1)+1)) * &
         dom%R_F_wgt%elts(idSE+1)%enc) ! UPLT S

    val(2) = sum (h_flux(id_edges([WMP,UZP,VPP]+1)+1)* dom%R_F_wgt%elts(id+1)%enc) &
         + sum ((h_flux(id_edges([VMP,WPPP,UPZ]+1)+1) - h_flux(id_edges([UMZ,VMPP,WPP]+1)+1))* &
         dom%R_F_wgt%elts(idN+1)%enc) ! LORT

    call get_indices (dom, i, j+1, UP, offs, dims, id_edges)

    val(3) = - sum (h_flux(id_edges([UZM,VMM,WPM]+1)+1) * dom%R_F_wgt%elts(idE+1)%enc) &
         - sum ((h_flux(id_edges([WMMM,UMZ,VPM]+1)+1) - h_flux(id_edges([VPMM,WMM,UPZ]+1)+1))* &
         dom%R_F_wgt%elts(idNE+1)%enc) ! UPLT

    val(4) = sum (h_flux(id_edges([UZP,VPP,WMP]+1)+1) * dom%R_F_wgt%elts(id2W+1)%enc) &
         + sum ((h_flux(id_edges([WPPP,UPZ,VMP]+1)+1) - h_flux(id_edges([VMPP,WPP,UMZ]+1)+1))* &
         dom%R_F_wgt%elts(idN2W+1)%enc) ! LORT W
  end function interp_flux

  function complete_coarse_flux (sm_flux, dom, i_par, j_par, i_chd, j_chd, e, offs_chd, dims_chd) result(val)
    
    implicit none

    real(dp),     intent(in) :: sm_flux(4)
    type(Domain), intent(in) :: dom
    integer,      intent(in) :: i_chd, j_chd, i_par, j_par
    integer,      intent(in) :: offs_chd(N_BDRY+1)
    integer,      intent(in) :: dims_chd(2,N_BDRY+1)
    real(dp)                 :: val
    
    integer  :: e
    real(dp) :: p_flux, c_flux

    val = 0.0_dp
    if (e == RT) then
       p_flux = part_coarse_flux (dom, i_chd+1, j_chd, RT, offs_chd, dims_chd)
       c_flux = coarse_flux (dom, i_par, j_par, i_chd+1, j_chd, offs_chd, dims_chd, RT)
       
       val = p_flux + c_flux + sm_flux(1) + sm_flux(2)
    elseif (e == DG) then
       p_flux = part_coarse_flux (dom, i_chd+1, j_chd+1, DG, offs_chd, dims_chd)
       c_flux = coarse_flux (dom, i_par, j_par, i_chd+1, j_chd+1, offs_chd, dims_chd, DG)
       
       val = p_flux + c_flux + sm_flux(2) + sm_flux(3)
    elseif (e == UP) then
       p_flux = part_coarse_flux (dom, i_chd, j_chd+1, UP, offs_chd, dims_chd)
       c_flux = coarse_flux (dom, i_par, j_par, i_chd, j_chd+1, offs_chd, dims_chd, UP)
       
       val = p_flux + c_flux + sm_flux(3) + sm_flux(4)
    end if
  end function complete_coarse_flux

  
  function part_coarse_flux (dom, i, j, e, offs, dims) result(val)
    
    implicit none
    
    type(Domain), intent(in) :: dom
    integer,      intent(in) :: i, j, e
    integer,      intent(in) :: offs(N_BDRY+1) 
    integer,      intent(in) :: dims(2,N_BDRY+1)
    real(dp)                 :: val

    integer                :: id
    integer, dimension(20) :: id_edges
    real(dp), dimension(2) :: area
    real(dp), dimension(4) :: ol_area

    id = idx (i, j, offs, dims)

    call get_indices (dom, i, j, e, offs, dims, id_edges)

    area         = dom%overl_areas%elts(id+1)%a(1:2)
    ol_area(1:2) = dom%overl_areas%elts(id+1)%split
    ol_area(3:4) = dom%overl_areas%elts(id+1)%a(3:4) - ol_area(1:2)

    area(1) = area(1) + ol_area(1) + ol_area(4)
    area(2) = area(2) + ol_area(2) + ol_area(3)
    area = area / sum(area)

    ol_area(1) = dom%overl_areas%elts(id_edges(PP+1)+1)%split(1)
    ol_area(2) = dom%overl_areas%elts(id_edges(MM+1)+1)%split(2)
    ol_area(3) = dom%overl_areas%elts(id_edges(MP+1)+1)%a(3) - dom%overl_areas%elts(id_edges(MP+1)+1)%split(1)
    ol_area(4) = dom%overl_areas%elts(id_edges(PM+1)+1)%a(4) - dom%overl_areas%elts(id_edges(PM+1)+1)%split(2)

    val = sum (h_flux(id_edges([UPZ,UMZ]+1)+1)*area) - sum (h_flux(id_edges([VMM,WMP]+1)+1))*area(2) &
         - sum (h_flux(id_edges([WPM,VPP]+1)+1))*area(1) &
         + ol_area(3) * dscalar(id_edges(MP+1)+1) - ol_area(4) * dscalar(id_edges(PM+1)+1)  &
         - ol_area(1) * dscalar(id_edges(PP+1)+1) + ol_area(2) * dscalar(id_edges(MM+1)+1)
  end function part_coarse_flux

  
  function coarse_flux (dom, i_par, j_par, i_chd, j_chd, offs_chd, dims_chd, e) result(val)
    
    implicit none
    
    type(Domain), intent(in) :: dom
    integer,      intent(in) :: i_par, j_par, i_chd, j_chd, e
    integer,      intent(in) :: offs_chd(N_BDRY+1) 
    integer,      intent(in) :: dims_chd(2,N_BDRY+1)
    real(dp)                 :: val

    integer :: id, id_mz, id_pz, id_mp,  id_pp, id_pm, id_mm, id_mm2, id_pm2, id_pp2,  id_mp2

    id_mz = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 1 + 1), offs_chd, dims_chd)
    id_pz = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 4 + 1), offs_chd, dims_chd)
    id_mp = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 0 + 1), offs_chd, dims_chd)
    id_pp = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 5 + 1), offs_chd, dims_chd)
    id_pm = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 3 + 1), offs_chd, dims_chd)
    id_mm = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 2 + 1), offs_chd, dims_chd)

    id_mm2 = idx2 (i_chd, j_chd, bfly_no2(:,1,e+1), offs_chd, dims_chd)
    id_pm2 = idx2 (i_chd, j_chd, bfly_no2(:,2,e+1), offs_chd, dims_chd)
    id_pp2 = idx2 (i_chd, j_chd, bfly_no2(:,3,e+1), offs_chd, dims_chd)
    id_mp2 = idx2 (i_chd, j_chd, bfly_no2(:,4,e+1), offs_chd, dims_chd)

    id = idx (i_chd, j_chd, offs_chd, dims_chd)

    val = (dom%overl_areas%elts(id+1)%a(1)*dom%overl_areas%elts(id+1)%a(2)*dom%areas%elts(id+1)%hex_inv &
         + dom%overl_areas%elts(id_mp+1)%a(2)*dom%overl_areas%elts(id_mp+1)%a(3)*dom%areas%elts(id_mp+1)%hex_inv &
         + dom%overl_areas%elts(id_pp+1)%a(1)*dom%overl_areas%elts(id_pp+1)%a(3)*dom%areas%elts(id_pp+1)%hex_inv &
         + dom%overl_areas%elts(id_pm+1)%a(1)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
         + dom%overl_areas%elts(id_mm+1)%a(2)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv) &
         * (dscalar(id_pz+1) - dscalar(id_mz+1)) + &
         dom%overl_areas%elts(id_pp+1)%a(3)*dom%overl_areas%elts(id_pp+1)%a(4)*dom%areas%elts(id_pp+1)%hex_inv &
         * 0.5_dp * (dscalar(id_pp2+1) - dscalar(id_mz+1)) + &
         dom%overl_areas%elts(id_pm+1)%a(3)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
         * 0.5_dp * (dscalar(id_pm2+1) - dscalar(id_mz+1)) + &
         dom%overl_areas%elts(id_mp+1)%a(3)*dom%overl_areas%elts(id_mp+1)%a(4)*dom%areas%elts(id_mp+1)%hex_inv &
         * 0.5_dp * (dscalar(id_pz+1) - dscalar(id_mp2+1)) + &
         dom%overl_areas%elts(id_mm+1)%a(3)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv &
         * 0.5_dp * (dscalar(id_pz+1) - dscalar(id_mm2+1))
  end function coarse_flux

  
  subroutine cal_divu_ml (q)
    ! Returns flux divergence of velocity in divF using solution q, stored in dom%divu
    
    implicit none
    
    type(Float_Field), target, intent(inout) :: q
    
    integer :: d, j, l

    call update_bdry (q, NONE, 976)

    do l = level_end, level_start, -1
       ! Calculate velocity flux
       do d = 1, size(grid)
          h_flux => horiz_flux(S_MASS)%data(d)%elts
          velo   => q%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=5)
          end do
          nullify (velo)
          if (l < level_end) then
             dscalar => grid(d)%divu%elts
             call cpt_or_restr_flux (grid(d), l) ! restrict flux if possible
             nullify (dscalar)
          end if
          nullify (h_flux)
       end do
       horiz_flux(S_MASS)%bdry_uptodate = .false.
       call update_bdry (horiz_flux(S_MASS), l, 977)

       ! Calculate divergence of velocity flux
       do d = 1, size(grid)
          dscalar => grid(d)%divu%elts
          h_flux  => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
    end do
  end subroutine cal_divu_ml

  
  subroutine get_indices (dom, i, j, e, offs, dims, id_edges)
    
    implicit none
    
    type(Domain), intent(in)  :: dom
    integer,      intent(in)  :: i, j, e
    integer,      intent(in)  :: offs(N_BDRY+1) 
    integer,      intent(in)  :: dims(2,N_BDRY+1)
    integer,      intent(out) :: id_edges(20)

    integer, dimension(2) :: ij_mp, ij_pp, ij_pm, ij_mm

    id_edges(UMZ+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 1 + 1), offs, dims)
    id_edges(UPZ+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 4 + 1), offs, dims)
    id_edges(WMP+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 0 + 1), offs, dims)
    id_edges(VPP+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 5 + 1), offs, dims)
    id_edges(WPM+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 3 + 1), offs, dims)
    id_edges(VMM+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 2 + 1), offs, dims)

    ij_mp = [i, j] + nghb_pt(:,hex_s_offs(e+1) + 0 + 1)
    id_edges(MP+1) = idx (ij_mp(1), ij_mp(2), offs, dims)

    ij_pp = [i, j] + nghb_pt(:,hex_s_offs(e+1) + 5 + 1)
    id_edges(PP+1) = idx (ij_pp(1), ij_pp(2), offs, dims)

    ij_pm = [i, j] + nghb_pt(:,hex_s_offs(e+1) + 3 + 1)
    id_edges(PM+1) = idx (ij_pm(1), ij_pm(2), offs, dims)

    ij_mm = [i, j] + nghb_pt(:,hex_s_offs(e+1) + 2 + 1)

    id_edges(MM+1)   = idx (ij_mm(1), ij_mm(2), offs, dims)

    id_edges(VMP+1)  = ed_idx (ij_mp(1), ij_mp(2), hex_sides (:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
    id_edges(VMPP+1) = ed_idx (ij_mp(1), ij_mp(2), hex_sides (:, hex_s_offs(e+1) + 1  + 4 + 1), offs, dims)
    id_edges(UZP+1)  = ed_idx (ij_mp(1), ij_mp(2), hex_sides (:, hex_s_offs(e+1) + 0  + 4 + 1), offs, dims)
    id_edges(WPPP+1) = ed_idx (ij_pp(1), ij_pp(2), hex_sides (:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
    id_edges(WPP+1)  = ed_idx (ij_pp(1), ij_pp(2), hex_sides (:, hex_s_offs(e+1) + 1  + 2 + 1), offs, dims)
    id_edges(VPM+1)  = ed_idx (ij_pm(1), ij_pm(2), hex_sides (:, hex_s_offs(e+1) + 1  + 4 + 1), offs, dims)
    id_edges(VPMM+1) = ed_idx (ij_pm(1), ij_pm(2), hex_sides (:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
    id_edges(UZM+1)  = ed_idx (ij_pm(1), ij_pm(2), hex_sides (:,(hex_s_offs(e+1) + 3) - 2 + 1), offs, dims)
    id_edges(WMMM+1) = ed_idx (ij_mm(1), ij_mm(2), hex_sides (:, hex_s_offs(e+1) + 1  + 2 + 1), offs, dims)
    id_edges(WMM+1)  = ed_idx (ij_mm(1), ij_mm(2), hex_sides (:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
  end subroutine get_indices

  
end module multi_level_mod
