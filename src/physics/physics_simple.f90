!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: physics_simple.f90
! Author: Gabrielle Ching-Johnson, Nicholas Kevlahan
!
! Date Revised: 2025-05
!
! Contains all subroutines needed to compute a Backwards Euler step using simple physics package
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module physics_simple_mod
  
  use, intrinsic :: iso_c_binding, only : C_BOOL
  
  use kind_mod,   only : dp, sp
  use shared_mod, only : DAY, EDGE, N_BDRY, N_VARIABLE, zlevels, S_MASS, S_TEMP, S_VELO, RT, DG, UP, &
       Nsoil, dt, grav_accel, p_top, time, zlevels, zmin, z_null 

  use diagnostics_mod,    only : cal_surf_press, kinetic_energy, integrate_pressure_up
  use domain_mod,         only : Domain, exner, exner_fun, mass, temp, mean_m, mean_t, grid, sol, sol_mean, topography, id_edge, idx
  use domain_ops_mod,     only : apply_no_bdry2
  use geom_mod,           only : cart2sph
  use init_physics_mod,   only : physics_firstcall_flag
  use single_column_mod,  only : change_latitude_longitude, physics_call_single_col
  use utils_mod,          only : interp

#ifdef PHYSICS_COLUMN_LOCAL
  use iso_fortran_env, only: int64
  use shared_mod, only: compressible,c_p,p_0,kappa
  use arch_mod, only: rank
  use domain_ops_mod, only: apply_onescale_to_patch
  use init_mod, only: surf_geopot
  use parallel_block_primitive_kernels_mod, only: pressure_layer,physics_kinetic_energy
#endif
#ifdef PHYSICS_OWNER
  use mpi_f08
  use arch_mod, only: comm,n_process,glo_id,loc_id,owner,block_catalog
  use domain_mod, only: get_offs_Domain
  use patch_mod, only: PATCH_SIZE
#endif
  implicit none

  private
  public :: physics_simple_step, Rayleigh_friction
  
#ifdef PHYSICS_COLUMN_LOCAL
  ! Integer visitation metadata only. Numerical hydrostatic state belongs to
  ! the active column; no replacement Domain-shaped pressure image is kept.
  type :: Physics_Visits
     integer,allocatable :: surface(:),seen(:)
  end type
  type(Physics_Visits),allocatable :: column_visits(:)
  integer(int64) :: column_count=0,layer_checks=0
#endif
#ifdef PHYSICS_VELOCITY_PLAN
  ! Local storage identities, not canonical block-owner identities. Stamp zero
  ! means the epoch-entry value; a positive stamp names the preceding writer.
  ! No numerical velocity image is retained. The walk also conservatively keeps
  ! read-before-write dependencies, including boundary aliases and pole reads.
  type :: Velocity_Plan
     integer,allocatable :: edges(:,:),sources(:,:),expected(:,:),version(:),ordinal(:)
#ifdef PHYSICS_OWNER
     integer,allocatable :: patch_catalog(:),node_catalog(:)
#endif
  end type
  type(Velocity_Plan),allocatable :: velocity_plan_by_domain(:)
  integer :: planned_columns=0,consumed_columns=0
  integer(int64) :: entry_reads=0,updated_reads=0,velocity_checks=0
#endif
#ifdef PHYSICS_OWNER
  integer,parameter :: OWNER_COMMAND=28101,OWNER_INPUT=28102,OWNER_OUTPUT=28103,OWNER_TOKEN=28104
  integer,parameter :: OWNER_HEADER=26
  integer,allocatable :: last_owner_column(:)
  logical :: owner_package_started=.false.
  integer(int64) :: owner_executed=0,owner_remote=0,owner_scaffold=0,owner_checks=0
#endif
  real(sp) :: phys_dt
contains
  subroutine physics_simple_step (h)
    ! Uses simple physics modules to take a Backwards Euler step for physics using time step dt set by dynamics
    implicit none
    real(dp) :: h ! time step
#ifdef PHYSICS_COLUMN_LOCAL
    integer :: d,p
#ifdef PHYSICS_OWNER
    integer :: source,destination,ierr,header(OWNER_HEADER),b
#endif
#endif
    
    phys_dt = real (h, kind=sp)
    
#ifdef PHYSICS_COLUMN_LOCAL
    if(.not.compressible) error stop "local Simple-physics columns require compressible dynamics"
    if(allocated(column_visits)) error stop "unfinished Simple-physics column epoch"
    allocate(column_visits(size(grid)))
    column_count=0;layer_checks=0
    do d=1,size(grid)
       allocate(column_visits(d)%surface(grid(d)%node%length),column_visits(d)%seen(grid(d)%node%length))
       column_visits(d)%surface=0;column_visits(d)%seen=0
       do p=3,grid(d)%patch%length
          call apply_onescale_to_patch(register_surface_visit,grid(d),p-1,0,0,1)
       end do
    end do
#endif
#if !defined(PHYSICS_COLUMN_LOCAL) || defined(PHYSICS_COLUMN_ORACLE)
    call cal_surf_press (sol(1:N_VARIABLE,1:zlevels))
#endif

#ifdef PHYSICS_VELOCITY_PLAN
    if(allocated(velocity_plan_by_domain)) error stop "unfinished physics velocity epoch"
    allocate(velocity_plan_by_domain(size(grid)))
    planned_columns=0;consumed_columns=0
    entry_reads=0;updated_reads=0;velocity_checks=0
    do d=1,size(grid)
       associate(plan=>velocity_plan_by_domain(d))
         allocate(plan%edges(6,grid(d)%node%length),plan%sources(6,grid(d)%node%length), &
              plan%expected(6,grid(d)%node%length),plan%version(grid(d)%node%length), &
              plan%ordinal(grid(d)%node%length))
         plan%version=0;plan%ordinal=0
       end associate
    end do
#ifdef PHYSICS_OWNER
    if(.not.allocated(block_catalog)) error stop "physics owner catalog absent"
    allocate(last_owner_column(n_process));last_owner_column=0
    owner_executed=0;owner_remote=0;owner_scaffold=0;owner_checks=0
    do d=1,size(grid)
       allocate(velocity_plan_by_domain(d)%patch_catalog(grid(d)%patch%length), &
            velocity_plan_by_domain(d)%node_catalog(grid(d)%node%length))
       velocity_plan_by_domain(d)%patch_catalog=0
       velocity_plan_by_domain(d)%node_catalog=0
    end do
    do b=1,size(block_catalog)
       if(owner(block_catalog(b)%root_domain+1)/=rank) cycle
       d=loc_id(block_catalog(b)%root_domain+1)+1
       call map_physics_subtree(d,block_catalog(b)%root_patch,b)
    end do
#endif
    call apply_no_bdry2(plan_velocity_column,z_null)
    do d=1,size(grid)
       velocity_plan_by_domain(d)%version=0
    end do
#endif
    ! Compute Simple Physics split step on all columns
#ifdef PHYSICS_OWNER
    ! Only one geometry source advances at a time. Every other rank services
    ! requests, including ranks with no Domains. A reply precedes the next
    ! dependent gather; neither entry aliases nor pending wavelets are refreshed.
    do source=0,n_process-1
       if(rank==source) then
          call apply_no_bdry2 (physics_call, z_null)
          header=0
          do destination=0,n_process-1
             if(destination==rank) cycle
             call MPI_Send(header,OWNER_HEADER,MPI_INTEGER,destination,OWNER_COMMAND,comm,ierr)
             call check_owner_mpi(ierr)
          end do
       else
          call serve_physics_owner(source)
       end if
       call MPI_Barrier(comm,ierr)
       call check_owner_mpi(ierr)
    end do
    write(*,'(a,i0,a,4(1x,i0))') 'Physics owner rank=',rank, &
         ' executed/remote/scaffold/checks=',owner_executed,owner_remote,owner_scaffold,owner_checks
    deallocate(last_owner_column)
#else
    call apply_no_bdry2 (physics_call, z_null)
#endif

#ifdef PHYSICS_VELOCITY_PLAN
    if(consumed_columns/=planned_columns) error stop "incomplete physics velocity plan"
    write(*,'(a,i0,a,4(1x,i0))') 'Physics velocity plan rank=',rank, &
         ' columns/entry/updated/checks=',consumed_columns,entry_reads,updated_reads,velocity_checks
    deallocate(velocity_plan_by_domain)
#endif
    sol%bdry_uptodate      = .false.
    physics_firstcall_flag = .false.
#ifdef PHYSICS_COLUMN_LOCAL
    write(*,'(a,i0,a,2(1x,i0))') 'Local physics columns rank=',rank,' columns/layer checks=',column_count,layer_checks
    deallocate(column_visits)
#endif
  end subroutine physics_simple_step

#ifdef PHYSICS_COLUMN_LOCAL
  subroutine register_surface_visit(dom,i,j,k,offs,dims)
    type(Domain),intent(inout)::dom
    integer,intent(in)::i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer::id
    id=idx(i,j,offs,dims)+1
    column_visits(dom%id+1)%surface(id)=column_visits(dom%id+1)%surface(id)+1
  end subroutine
#ifdef PHYSICS_COLUMN_ORACLE
  subroutine check_physics_value(actual,expected)
    real(dp),intent(in)::actual,expected
    if(transfer(actual,0_int64)/=transfer(expected,0_int64)) then
       write(*,*) 'Physics hydrostatic mismatch rank, local, reference:',rank,actual,expected
       error stop "independent physics hydrostatic mismatch"
    end if
  end subroutine
#endif
#endif

#ifdef PHYSICS_VELOCITY_PLAN
  subroutine plan_velocity_column(dom,p_null,i,j,zlev,offs,dims,is_pole)
    type(Domain),intent(inout)::dom
    integer,intent(in)::p_null,i,j,zlev,offs(N_BDRY+1),dims(2,N_BDRY+1),is_pole
    integer::d,id,south,southwest,west
    d=dom%id+1;id=idx(i,j,offs,dims)
    south=idx(i,j-1,offs,dims);southwest=idx(i-1,j-1,offs,dims);west=idx(i-1,j,offs,dims)
    associate(plan=>velocity_plan_by_domain(d))
      if(plan%ordinal(id+1)/=0) error stop "duplicate physics velocity plan column"
      planned_columns=planned_columns+1
      plan%ordinal(id+1)=planned_columns
#ifdef PHYSICS_OWNER
      plan%node_catalog(id+1)=plan%patch_catalog(p_null+1)
#endif
      plan%edges(:,id+1)=[EDGE*id+UP+1,EDGE*id+DG+1,EDGE*id+RT+1, &
           EDGE*south+UP+1,EDGE*southwest+DG+1,EDGE*west+RT+1]
      plan%sources(:,id+1)=[id+1,id+1,id+1,south+1,southwest+1,west+1]
      plan%expected(:,id+1)=plan%version(plan%sources(:,id+1))
      entry_reads=entry_reads+count(plan%expected(:,id+1)==0)
      updated_reads=updated_reads+count(plan%expected(:,id+1)>0)
      ! Pole columns update temperature/soil but do not write velocity.
      if(is_pole/=1) plan%version(id+1)=planned_columns
    end associate
  end subroutine plan_velocity_column
#endif

#ifdef PHYSICS_OWNER
  recursive subroutine map_physics_subtree(d,p,b)
    integer,intent(in)::d,p,b
    integer::offs(N_BDRY+1),dims(2,N_BDRY+1),i,j,c,child,id
    if(p<0.or.p>=grid(d)%patch%length) error stop "physics owner patch outside Domain"
    if(velocity_plan_by_domain(d)%patch_catalog(p+1)/=0) error stop "overlapping physics owner subtrees"
    velocity_plan_by_domain(d)%patch_catalog(p+1)=b
    call get_offs_Domain(grid(d),p,offs,dims)
    do j=0,PATCH_SIZE-1
       do i=0,PATCH_SIZE-1
          id=idx(i,j,offs,dims)+1
          velocity_plan_by_domain(d)%node_catalog(id)=b
       end do
    end do
    do c=1,size(grid(d)%patch%elts(p+1)%children)
       child=grid(d)%patch%elts(p+1)%children(c)
       if(child>0) call map_physics_subtree(d,child,b)
    end do
  end subroutine map_physics_subtree

  subroutine check_owner_mpi(ierr)
    integer,intent(in)::ierr
    if(ierr/=MPI_SUCCESS) error stop "physics owner MPI failure"
  end subroutine

  subroutine receive_owner_header(source,header)
    integer,intent(in)::source
    integer,intent(out)::header(OWNER_HEADER)
    type(MPI_Status)::status
    integer::ierr,n
    call MPI_Recv(header,OWNER_HEADER,MPI_INTEGER,source,OWNER_COMMAND,comm,status,ierr)
    call check_owner_mpi(ierr)
    call MPI_Get_count(status,MPI_INTEGER,n,ierr)
    call check_owner_mpi(ierr)
    if(n/=OWNER_HEADER) error stop "physics owner header extent mismatch"
  end subroutine

  subroutine execute_owner_packet(header,packet)
    integer,intent(in)::header(OWNER_HEADER)
    real(sp),intent(inout)::packet(8*zlevels+Nsoil+8)
    real(sp)::bootstrap(size(packet))
    logical(kind=C_BOOL)::first
    integer::n,f
    n=zlevels;f=8*n+Nsoil+3
    first=header(7)==1
    ! A previously idle owner may receive its first column after adaptation.
    ! Initialize package-private allocation state using a disposable column,
    ! then reload the requested soil and atmospheric inputs for the real call.
    if(.not.owner_package_started.and..not.first) then
       bootstrap=packet
       call invoke(.true._C_BOOL,bootstrap)
    end if
    call invoke(first,packet)
    owner_package_started=.true.
  contains
    subroutine invoke(first_call,q)
      logical(kind=C_BOOL),intent(in)::first_call
      real(sp),intent(inout)::q(:)
      call change_latitude_longitude(q(f+1),q(f+2))
      call physics_call_single_col(1,n,header(8),first_call,.false._C_BOOL,q(f+3),q(f+4),q(f+5), &
           q(1:n),q(n+1:2*n+1),q(2*n+2:3*n+1),q(8*n+2:8*n+2),q(3*n+2:4*n+1), &
           q(4*n+2:5*n+1),q(5*n+2:6*n+1),q(6*n+2:7*n+1),q(7*n+2:8*n+1),q(8*n+3:f))
    end subroutine
  end subroutine execute_owner_packet

  subroutine validate_owner_request(source,header)
    integer,intent(in)::source,header(OWNER_HEADER)
    integer::b,s
    b=header(2)
    if(header(1)/=1.or.header(3)<0.or.header(3)>=size(owner)) error stop "invalid physics owner request"
    if(owner(header(3)+1)/=source) error stop "physics owner source rank mismatch"
    if(b<0.or.b>size(block_catalog)) error stop "invalid physics owner block"
    if(b==0) then
       if(source/=rank) error stop "physics scaffold routed remotely"
    else
       if(block_catalog(b)%owner/=rank.or.block_catalog(b)%root_domain/=header(3)) &
            error stop "physics owner destination mismatch"
    end if
    if(header(6)<=last_owner_column(source+1)) error stop "physics owner duplicate or reordered request"
    if(any(header(9:14)<0).or.any(header(9:14)>=header(6))) error stop "physics owner invalid source version"
    if(any(header(15:20)<1)) error stop "physics owner invalid source identity"
    do s=21,26
       if(header(s)<0.or.header(s)>size(block_catalog)) error stop "physics owner invalid source block"
       if(header(s)>0) then
          if(block_catalog(header(s))%root_domain/=header(3)) error stop "physics owner source domain mismatch"
       end if
    end do
    last_owner_column(source+1)=header(6)
  end subroutine validate_owner_request

  subroutine serve_physics_owner(source)
    integer,intent(in)::source
    integer::header(OWNER_HEADER),ierr,n
    real(sp)::packet(8*zlevels+Nsoil+8)
    type(MPI_Status)::status
    do
       call receive_owner_header(source,header)
       if(header(1)==0) exit
       call validate_owner_request(source,header)
       call MPI_Recv(packet,size(packet),MPI_REAL,source,OWNER_INPUT,comm,status,ierr)
       call check_owner_mpi(ierr)
       call MPI_Get_count(status,MPI_REAL,n,ierr)
       call check_owner_mpi(ierr)
       if(n/=size(packet)) error stop "physics owner packet extent mismatch"
       call execute_owner_packet(header,packet)
       owner_executed=owner_executed+1;owner_remote=owner_remote+1
       call MPI_Send(header,OWNER_HEADER,MPI_INTEGER,source,OWNER_TOKEN,comm,ierr)
       call check_owner_mpi(ierr)
       call MPI_Send(packet,size(packet),MPI_REAL,source,OWNER_OUTPUT,comm,ierr)
       call check_owner_mpi(ierr)
    end do
  end subroutine serve_physics_owner

  subroutine dispatch_physics_owner(d,p,id,mask,latitude,longitude,day,fraction, &
       play,pint,phi,phisurf,umag,u,v,w,theta,soil)
    integer,intent(in)::d,p,id,mask
    real(sp),intent(in)::latitude,longitude,day,fraction,play(:),pint(0:),phi(:),phisurf(:),umag(:)
    real(sp),intent(inout)::u(:),v(:),w(:),theta(:),soil(:)
    integer::header(OWNER_HEADER),reply(OWNER_HEADER),destination,b,ierr,n,count_received
    real(sp)::packet(8*zlevels+Nsoil+8)
#ifdef PHYSICS_OWNER_ORACLE
    real(sp)::reference(size(packet))
#endif
    type(MPI_Status)::status
    n=zlevels
    associate(plan=>velocity_plan_by_domain(d))
      b=plan%patch_catalog(p+1)
      header(1:8)=[1,b,glo_id(rank+1,d),p,id,consumed_columns,merge(1,0,physics_firstcall_flag),mask]
      header(9:14)=plan%expected(:,id)
      header(15:20)=plan%sources(:,id)
      header(21:26)=plan%node_catalog(plan%sources(:,id))
    end associate
    packet=[play,pint,phi,umag,u,v,w,theta,phisurf,soil,latitude,longitude,day,fraction,phys_dt]
#ifdef PHYSICS_OWNER_ORACLE
    reference=packet
#endif
    destination=rank
    if(b>0) destination=block_catalog(b)%owner
    if(destination==rank) then
       call validate_owner_request(rank,header)
       call execute_owner_packet(header,packet)
       if(b>0) then
          owner_executed=owner_executed+1
       else
          owner_scaffold=owner_scaffold+1
       end if
    else
       call MPI_Send(header,OWNER_HEADER,MPI_INTEGER,destination,OWNER_COMMAND,comm,ierr)
       call check_owner_mpi(ierr)
       call MPI_Send(packet,size(packet),MPI_REAL,destination,OWNER_INPUT,comm,ierr)
       call check_owner_mpi(ierr)
       call MPI_Recv(reply,OWNER_HEADER,MPI_INTEGER,destination,OWNER_TOKEN,comm,status,ierr)
       call check_owner_mpi(ierr)
       call MPI_Get_count(status,MPI_INTEGER,count_received,ierr)
       call check_owner_mpi(ierr)
       if(count_received/=OWNER_HEADER.or.any(reply/=header)) error stop "physics owner reply identity mismatch"
       call MPI_Recv(packet,size(packet),MPI_REAL,destination,OWNER_OUTPUT,comm,status,ierr)
       call check_owner_mpi(ierr)
       call MPI_Get_count(status,MPI_REAL,count_received,ierr)
       call check_owner_mpi(ierr)
       if(count_received/=size(packet)) error stop "physics owner reply extent mismatch"
    end if
#ifdef PHYSICS_OWNER_ORACLE
    call execute_owner_packet(header,reference)
    if(any(transfer(packet,[0],size(packet))/=transfer(reference,[0],size(reference)))) &
         error stop "independent physics owner output mismatch"
    owner_checks=owner_checks+1
#endif
    u=packet(4*n+2:5*n+1);v=packet(5*n+2:6*n+1);w=packet(6*n+2:7*n+1)
    theta=packet(7*n+2:8*n+1);soil=packet(8*n+3:8*n+Nsoil+3)
  end subroutine dispatch_physics_owner
#endif

  subroutine Rayleigh_friction (dom, i, j, zlev, offs, dims)
    ! Rayleigh friction in atmospheric boundary layer
    ! (based on Held-Suarez)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    character(5) :: smth = "tropo"          ! where to apply Rayleigh friction (tropo, strat)
    real(dp)     :: k_f     = 1.0_dp  / DAY ! Rayleigh friction
    real(dp)     :: sigma_b1 = 0.7_dp       ! normalized tropopause pressure height
    real(dp)     :: sigma_b2 = 0.1_dp       ! normalized stratosphere friction pressure height

    integer                        :: d, id, id_i, k, l
    integer,  dimension(1:EDGE)    :: id_e
    real(dp)                       :: k_v, rho_dz, sigma
    real(dp), dimension(0:zlevels) :: Pl
    real(dp), dimension(1:zlevels) :: Pk

    d = dom%id + 1
    
    id   = idx (i, j, offs, dims)
    id_i = id + 1
    id_e = id_edge (id)
    
    Pl(zlevels) = p_top
    do l = zlevels-1, 0, -1
       k = l + 1
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(id_i) 
       Pl(l)  = Pl(l+1) + grav_accel * rho_dz
       Pk(k) = interp (Pl(l), Pl(l+1))
    end do

    do k = 1, zlevels
       sigma = (Pk(k) - p_top) / (Pl(0) - p_top)

       select case (smth)
       case ("tropo")
          k_v = k_f * max (0.0_dp,   (sigma - sigma_b1) / (1.0_dp - sigma_b1))
       case ("strat")
          k_v = k_f * max (0.0_dp, - (sigma - sigma_b2) / sigma_b2)
       case default
          k_v = 0.0_dp
       end select

       sol(S_VELO,k)%data(d)%elts(id_e) = (1.0_dp - dt * k_v) * sol(S_VELO,k)%data(d)%elts(id_e)
    end do
  end subroutine Rayleigh_friction

  subroutine physics_call (dom, p_null, i, j, zlev, offs, dims, is_pole)
    !-----------------------------------------------------------------------------------
    !
    !   Backwards Euler physics step on a single element/column
    !
    !-----------------------------------------------------------------------------------

    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, is_pole, p_null, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, id_i, k, mask

    ! Variables needed for simple physics
    real(sp)                       :: nth_day, day_fraction   ! day in simulation, fraction of the day
    real(sp), dimension(1:zlevels) :: phys_Play               ! pressure at layer centres
    real(sp), dimension(0:zlevels) :: phys_Pint               ! pressure at layer interfaces
    real(sp), dimension(1:1)       :: phys_Phisurf            ! surface geopotential

    real(sp), dimension(1:zlevels) :: phys_Theta              ! potential temperature
    real(sp), dimension(1:zlevels) :: phys_Phi                ! geopotential
    real(sp), dimension(1:zlevels) :: phys_U, phys_V, phys_W  ! velocities at edges
    real(sp), dimension(1:zlevels) :: phys_Umag               ! speed at nodes
    real(sp), dimension(1:Nsoil+1) :: Tsoil                   ! surface temp and soil temperatures

    real(dp)                       :: latitude, longitude     ! coordinates of the column
    real(dp), dimension(1:zlevels) :: rho_dz
    
#ifndef PHYSICS_OWNER
    logical(kind=C_BOOL)           :: lastcall_flag = .false.
#endif

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    day_fraction = (real (time, kind=sp) - phys_dt) / real (DAY, kind=sp)
    nth_day      = floor (day_fraction)
    day_fraction = real (day_fraction - nth_day, kind=sp)

#ifdef PHYSICS_VELOCITY_PLAN
    if(.not.allocated(velocity_plan_by_domain)) error stop "physics velocity epoch absent"
    consumed_columns=consumed_columns+1
    associate(plan=>velocity_plan_by_domain(d))
      if(plan%ordinal(id_i)/=consumed_columns) error stop "physics velocity consumer order mismatch"
      if(any(plan%version(plan%sources(:,id_i))/=plan%expected(:,id_i))) &
           error stop "physics velocity source version mismatch"
    end associate
#endif
    ! Save prognostic variables in physics data structure
    call pack_physics_vars

    ! Get latitude and longitude of the column
    call cart2sph (dom%node%elts(id_i), longitude, latitude)

#ifdef PHYSICS_OWNER
    phys_Phisurf=real(grav_accel*topography%data(d)%elts(id_i),kind=sp)
    mask=grid(d)%mask_n%elts(id_i)
    call dispatch_physics_owner(d,p_null,id_i,mask,real(latitude,sp),real(longitude,sp),nth_day,day_fraction, &
         phys_Play,phys_Pint,phys_Phi,phys_Phisurf,phys_Umag,phys_U,phys_V,phys_W,phys_Theta,Tsoil)
#else
    ! Update physics latitude and longitude
    call change_latitude_longitude (real (latitude, kind=sp), real (longitude, kind=sp))

    ! Surface geopotential for this column
    phys_Phisurf = real (grav_accel * topography%data(d)%elts(id_i), kind=sp)

    ! Mask for adaptive grid
    mask = grid(d)%mask_n%elts(id_i)

    ! Backwards Euler step on current column
    call physics_call_single_col (1, zlevels, mask, &
         physics_firstcall_flag, lastcall_flag, nth_day, day_fraction, phys_dt, &
         phys_Play, phys_Pint, phys_Phi, phys_Phisurf, phys_Umag, &
         phys_U, phys_V, phys_W, phys_Theta, Tsoil) ! updated values

#endif

    ! Assign solution at t+h
    do k = 1, zlevels
       if (is_pole /= 1) sol(S_VELO,k)%data(d)%elts(id_edge(id)) = real ((/ phys_U(k), phys_V(k), phys_W(k) /), kind=dp)
       sol(S_TEMP,k)%data(d)%elts(id_i) = rho_dz(k) * real (phys_Theta(k), kind=dp) - sol_mean(S_TEMP,k)%data(d)%elts(id_i)
    end do

#ifdef PHYSICS_VELOCITY_PLAN
    if(is_pole/=1) velocity_plan_by_domain(d)%version(id_i)=consumed_columns
#endif
    ! Assign soil column solution at t+phys_dt to WAVETRISK data structure
    do k = zmin, 0
       sol(S_TEMP,k)%data(d)%elts(id_i) = real (Tsoil(abs(k)+1), kind=dp)
    end do
  contains
    subroutine pack_physics_vars
      ! Gathers all required variables from wavetrisk data structure for all layers of column into physics data structure 
      integer  :: k
      real(dp) :: rho_dz_theta
#ifdef PHYSICS_VELOCITY_PLAN
      integer :: velocity_edges(6)
      real(dp) :: velocity_values(6),column_kinetic
#endif
#ifdef PHYSICS_COLUMN_LOCAL
      integer :: visit
      real(dp) :: column_surface,local_lower,local_geopot,local_pressure,local_exner,local_geopot_lower
      if(.not.allocated(column_visits)) error stop "Simple-physics column epoch absent"
      if(column_visits(d)%seen(id_i)/=0) error stop "duplicate Simple-physics column consumer"
      column_visits(d)%seen(id_i)=1
      column_count=column_count+1
      column_surface=0.0_dp
      do k=1,zlevels
         do visit=1,column_visits(d)%surface(id_i)
            column_surface=column_surface+(sol_mean(S_MASS,k)%data(d)%elts(id_i)+sol(S_MASS,k)%data(d)%elts(id_i))
         end do
      end do
      column_surface=grav_accel*column_surface+p_top
      local_lower=column_surface
      local_geopot=surf_geopot(d,id_i)
#ifdef PHYSICS_COLUMN_ORACLE
      call check_physics_value(column_surface,dom%surf_press%elts(id_i))
#endif
      phys_Pint(0)=real(column_surface,kind=sp)
#else

      phys_Pint(0)  = real (dom%surf_press%elts(id_i), kind=4)
#endif

      do k = 1, zlevels
         mass   =>      sol(S_MASS,k)%data(d)%elts
         temp   =>      sol(S_TEMP,k)%data(d)%elts
         mean_m => sol_mean(S_MASS,k)%data(d)%elts
         mean_t => sol_mean(S_TEMP,k)%data(d)%elts
#if !defined(PHYSICS_COLUMN_LOCAL) || defined(PHYSICS_COLUMN_ORACLE)
         exner  =>       exner_fun(k)%data(d)%elts
#endif

         rho_dz(k)    = mass(id_i) + mean_m(id_i)
         rho_dz_theta = temp(id_i) + mean_t(id_i)

         ! Pressure at layers and interface and geopotential at next interface (set in dom%geopot)
#ifdef PHYSICS_COLUMN_LOCAL
         call pressure_layer(mass(id_i),mean_m(id_i),temp(id_i),mean_t(id_i),grav_accel,c_p,p_0,kappa, &
              local_lower,local_geopot,local_pressure,local_exner,local_geopot_lower)
#ifdef PHYSICS_COLUMN_ORACLE
         call integrate_pressure_up (dom, i, j, zlev, offs, dims)
         call check_physics_value(local_lower,dom%press_lower%elts(id_i))
         call check_physics_value(local_pressure,dom%press%elts(id_i))
         call check_physics_value(local_geopot,dom%geopot%elts(id_i))
         call check_physics_value(local_geopot_lower,dom%geopot_lower%elts(id_i))
         call check_physics_value(local_exner,exner(id_i))
         layer_checks=layer_checks+1
#endif
#else
         call integrate_pressure_up (dom, i, j, zlev, offs, dims)
#endif

         ! Input variables for simple physics module
#ifdef PHYSICS_COLUMN_LOCAL
         phys_Pint(k)=real(local_lower,kind=sp)
         phys_Play(k)=real(local_pressure,kind=sp)
         phys_Phi(k)=real(interp(local_geopot_lower,local_geopot),kind=sp)
#else
         phys_Pint(k)  = real (dom%press_lower%elts(id_i),                                  kind=sp)
         phys_Play(k)  = real (dom%press%elts(id_i),                                        kind=sp)
         phys_Phi(k)   = real (interp (dom%geopot_lower%elts(id_i), dom%geopot%elts(id_i)), kind=sp)
#endif
#ifdef PHYSICS_VELOCITY_PLAN
         velocity_edges=velocity_plan_by_domain(d)%edges(:,id_i)
         velocity_values=sol(S_VELO,k)%data(d)%elts(velocity_edges)
         column_kinetic=physics_kinetic_energy(velocity_values,dom%len%elts(velocity_edges), &
              dom%pedlen%elts(velocity_edges),dom%areas%elts(id_i)%hex_inv)
#ifdef PHYSICS_VELOCITY_ORACLE
         if(transfer(column_kinetic,0_int64)/=transfer(kinetic_energy(dom,i,j,k,offs,dims),0_int64)) &
              error stop "independent physics velocity mismatch"
         velocity_checks=velocity_checks+1
#endif
         phys_Umag(k)=real(sqrt(2*column_kinetic),kind=sp)
         phys_U(k)=real(velocity_values(3),kind=sp)
         phys_V(k)=real(velocity_values(2),kind=sp)
         phys_W(k)=real(velocity_values(1),kind=sp)
#else
         phys_Umag(k)  = real (sqrt (2 * kinetic_energy (dom, i, j, k, offs, dims)),        kind=sp)
         phys_U(k)     = real (sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1),                    kind=sp)
         phys_V(k)     = real (sol(S_VELO,k)%data(d)%elts(EDGE*id+DG+1),                    kind=sp)
         phys_W(k)     = real (sol(S_VELO,k)%data(d)%elts(EDGE*id+UP+1),                    kind=sp)
#endif
         phys_Theta(k) = real (rho_dz_theta / rho_dz(k),                                    kind=sp) 

         nullify (mass, temp, mean_m, mean_t, exner)
      end do

      ! Retrieve surface temperature and soil column temperatures from dynamics
      ! (Tsoil(1) is surface temperature)
      do k = zmin, 0
         Tsoil(abs(k)+1) = real (sol(S_TEMP,k)%data(d)%elts(id_i), kind=sp)
      end do
    end subroutine pack_physics_vars
  end subroutine physics_call
end module physics_simple_mod
