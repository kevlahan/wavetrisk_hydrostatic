module types
implicit none
integer,parameter::N_BDRY=0,PATCH_SIZE=4,FROZEN=1,S_VELO=2,z_null=0,scalars(2)=[1,1]
type ints
 integer::length=1
 integer,allocatable::elts(:)
end type
type patch
 integer::children(4)=[1,0,0,0]
end type
type patches
 type(patch),allocatable::elts(:)
end type
type Domain
 type(ints),allocatable::lev(:)
 type(patches)::patch
 type(ints)::mask_n
end type
type values
 real,allocatable::elts(:)
end type
type Float_Field
 type(values),allocatable::data(:)
end type
type(Domain),allocatable::grid(:)
real,pointer::scalar(:),wc_s(:),velo(:),wc_u(:)
contains
integer function idx(i,j,o,d)
 integer::i,j,o(:),d(:,:)
 idx=0
end function
end module
module walker
use types
contains
subroutine apply_interscale_d(f,dom,level,z,a,b)
 interface
  subroutine f(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
   use types
   type(Domain),intent(inout)::dom
   integer,intent(in)::ip,jp,ic,jc,zlev,op(N_BDRY+1),dp0(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
  end subroutine
 end interface
 type(Domain)::dom
 integer::level,z,a,b
 call f(dom,0,0,0,0,0,[0],reshape([1,1],[2,1]),[0],reshape([1,1],[2,1]))
end subroutine
end module
module oracle
use types
use walker
use iso_fortran_env,only:int64
  logical :: subdivision_oracle_active = .false.
  integer :: subdivision_oracle_phase = 0
  integer(int64) :: subdivision_oracle_counts(4) = 0_int64

contains
  subroutine subdivision_inverse_oracle(wavelet,scaling,level,phase,counts)
    ! Independent legacy phase oracle for storage cuts. No production caller.
    ! Full geometry preserves parent/child operations across block boundaries.
    type(Float_Field),target,intent(inout) :: wavelet(:,:),scaling(:,:)
    integer,intent(in) :: level,phase
    integer(int64),intent(out) :: counts(4)
    integer :: d,k,v,p,i
    if(subdivision_oracle_active) error stop 'subdivision oracle is not reentrant'
    if(phase<1.or.phase>3) error stop 'invalid subdivision oracle phase'
    subdivision_oracle_active=.true.
    subdivision_oracle_phase=phase
    subdivision_oracle_counts=0_int64
    counts=0_int64
    do d=1,size(grid)
       do i=1,grid(d)%lev(level)%length
          p=grid(d)%lev(level)%elts(i)
          if(any(grid(d)%patch%elts(p+1)%children>0)) subdivision_oracle_counts(1)=subdivision_oracle_counts(1)+1_int64
       end do
       do k=1,size(scaling,2)
          if(phase<=2) then
             do v=scalars(1),scalars(2)
                scalar=>scaling(v,k)%data(d)%elts;wc_s=>wavelet(v,k)%data(d)%elts
                if(phase==1) then
                   call apply_interscale_d(subdivision_scalar_phase,grid(d),level,z_null,0,0)
                else
                   call apply_interscale_d(subdivision_scalar_phase,grid(d),level,z_null,0,0)
                end if
                nullify(scalar,wc_s)
             end do
          else
             velo=>scaling(S_VELO,k)%data(d)%elts;wc_u=>wavelet(S_VELO,k)%data(d)%elts
             call apply_interscale_d(subdivision_vector_phase,grid(d),level,z_null,0,0)
             nullify(velo,wc_u)
          end if
       end do
    end do
    counts=subdivision_oracle_counts
    subdivision_oracle_counts=0_int64
    subdivision_oracle_phase=0
    subdivision_oracle_active=.false.
  end subroutine subdivision_inverse_oracle

  subroutine subdivision_scalar_phase(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
      type(Domain),intent(inout) :: dom
      integer,intent(in) :: ip,jp,ic,jc,zlev,op(N_BDRY+1),dp0(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
      integer :: id
      if(.not.subdivision_oracle_active.or.subdivision_oracle_phase<1.or. &
           subdivision_oracle_phase>2) error stop 'inactive scalar subdivision oracle'
      if(subdivision_oracle_phase==1) then
         if(ic<PATCH_SIZE.and.jc<PATCH_SIZE) then
            subdivision_oracle_counts(2)=subdivision_oracle_counts(2)+1_int64;id=idx(ic,jc,oc,dc)
            if(dom%mask_n%elts(id+1)==FROZEN) then
               subdivision_oracle_counts(4)=subdivision_oracle_counts(4)+1_int64
            else
               subdivision_oracle_counts(3)=subdivision_oracle_counts(3)+1_int64
            end if
         end if
         call Prolong_scalar(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
      else
         subdivision_oracle_counts(2)=subdivision_oracle_counts(2)+3_int64
         subdivision_oracle_counts(3)=subdivision_oracle_counts(3)+3_int64
         call Reconstruct_scalar(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
      end if
    end subroutine
  subroutine subdivision_vector_phase(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
      type(Domain),intent(inout) :: dom
      integer,intent(in) :: ip,jp,ic,jc,zlev,op(N_BDRY+1),dp0(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
      if(.not.subdivision_oracle_active.or.subdivision_oracle_phase/=3) &
           error stop 'inactive vector subdivision oracle'
      subdivision_oracle_counts(2)=subdivision_oracle_counts(2)+1_int64
      subdivision_oracle_counts(3)=subdivision_oracle_counts(3)+6_int64
      call Reconstruct_inner_velo(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
    end subroutine


subroutine Prolong_scalar(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
 type(Domain),intent(inout)::dom
 integer,intent(in)::ip,jp,ic,jc,zlev,op(N_BDRY+1),dp0(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
 scalar=scalar+1
end subroutine
subroutine Reconstruct_scalar(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
 type(Domain),intent(inout)::dom
 integer,intent(in)::ip,jp,ic,jc,zlev,op(N_BDRY+1),dp0(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
 scalar=scalar+2
end subroutine
subroutine Reconstruct_inner_velo(dom,ip,jp,ic,jc,zlev,op,dp0,oc,dc)
 type(Domain),intent(inout)::dom
 integer,intent(in)::ip,jp,ic,jc,zlev,op(N_BDRY+1),dp0(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
 velo=velo+3
end subroutine
end module
program test
use oracle
use iso_fortran_env,only:int64
implicit none
type(Float_Field),target,allocatable::wav(:,:),sol(:,:)
integer(int64)::count(4),q
integer::i,k,d,nd,nk,frozen_value,repeat,probe
character(32)::arg
call get_command_argument(1,arg);read(arg,*)nd
call get_command_argument(2,arg);read(arg,*)nk
call get_command_argument(3,arg);read(arg,*)frozen_value
call get_command_argument(4,arg);read(arg,*)probe
allocate(grid(nd),wav(2,nk),sol(2,nk))
do d=1,nd
 allocate(grid(d)%lev(1),grid(d)%patch%elts(1),grid(d)%mask_n%elts(1))
 allocate(grid(d)%lev(1)%elts(1));grid(d)%lev(1)%elts=0;grid(d)%mask_n%elts=frozen_value
end do
do k=1,nk
 do i=1,2
  allocate(wav(i,k)%data(nd),sol(i,k)%data(nd))
  do d=1,nd
   allocate(wav(i,k)%data(d)%elts(1),sol(i,k)%data(d)%elts(1))
   wav(i,k)%data(d)%elts=0;sol(i,k)%data(d)%elts=0
  end do
 end do
end do
if(probe==1) call subdivision_inverse_oracle(wav,sol,1,0,count)
if(probe==2) call subdivision_inverse_oracle(wav,sol,1,4,count)
if(probe==3) then
 subdivision_oracle_active=.true.
 call subdivision_inverse_oracle(wav,sol,1,1,count)
end if
if(probe==4) call apply_interscale_d(subdivision_scalar_phase,grid(1),1,0,0,0)
if(probe==5) call apply_interscale_d(subdivision_vector_phase,grid(1),1,0,0,0)
q=int(nd*nk,int64)
do repeat=1,3
 do i=1,3
  call subdivision_inverse_oracle(wav,sol,1,i,count)
  if(i==1.and.any(count/=[int(nd,int64),q,merge(q,0_int64,frozen_value==0),merge(0_int64,q,frozen_value==0)])) stop 1
  if(i==2.and.any(count/=[int(nd,int64),3*q,3*q,0_int64])) stop 2
  if(i==3.and.any(count/=[int(nd,int64),q,6*q,0_int64])) stop 3
  if(subdivision_oracle_active.or.subdivision_oracle_phase/=0.or.any(subdivision_oracle_counts/=0)) stop 5
 end do
end do
do d=1,nd
 do k=1,nk
  if(sol(1,k)%data(d)%elts(1)/=9.or.sol(2,k)%data(d)%elts(1)/=9) stop 4
 end do
end do
print *,'PASS',nd,nk,frozen_value
end program
