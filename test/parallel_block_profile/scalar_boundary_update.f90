module fixture
use kind_mod, only: dp
use iso_fortran_env, only: int64
use parallel_block_scalar_storage_mod
implicit none
integer,parameter::EDGE=3,TRIAG=2,RT=0,UP=2,LORT=0,UPLT=1,TRSK=1,S_MASS=2,S_TEMP=3
include 'constants.inc'
real(dp),parameter::BLOCK_BOUNDARY_POISON=1.0e290_dp
integer::field_level,v_scalar=2,n_field_level=3
logical::validate_oracle,capture_dscalar=.false.,block_profile=.true.
integer(int64)::scalar_producer_work(7)=0_int64
type::Patch
 integer::elts_start=0
end type
type::Patches
 type(Patch)::elts(1)
end type
type::Areas
 real(dp)::hex_inv=0.25_dp
end type
type::AreaList
 type(Areas)::elts(32)
end type
type::Ints
 integer::elts(96)=2
end type
type::Reals
 real(dp)::elts(96)=2.0_dp
end type
type::Weight9
 real(dp)::enc(9)=0.5_dp
end type
type::Weights9
 type(Weight9)::elts(32)
end type
type::Weight3
 real(dp)::enc(3)=0.75_dp
end type
type::Weights3
 type(Weight3)::elts(32)
end type
type::Overlap
 real(dp)::a(4)=0.25_dp,split(2)=0.5_dp
end type
type::Overlaps
 type(Overlap)::elts(32)
end type
type::Domain
 type(Patches)::bdry_patch
 type(AreaList)::areas
 type(Ints)::mask_n,mask_e
 type(Reals)::len,triarea
 type(Weights9)::I_u_wgt
 type(Weights3)::R_F_wgt
 type(Overlaps)::overl_areas
end type
type(Domain)::grid(1)
type::Mass
 real(dp)::flux(96)
end type
type(Mass)::native_mass(1)
type::Closure
 integer::slot(96)
 real(dp)::value(2,96)
end type
type(Closure)::temperature_closure(1)
type::Field
 type(Reals)::data(1)
end type
type(Field)::horiz_flux(3),domain_tendency(3,2)
type::Storage
 type(Scalar_Record_Storage)::bdry
end type
type(Storage)::block_scalar_tendency(2)
type::Plan
 logical::full_transport
end type
type(Plan)::block_scalar_divergence_plan
contains
include 'new_update.inc'
include 'scalar_boundary_reference.inc'
include 'fill_node.inc'
end module
program test
use fixture
implicit none
integer::mode,full,v,k,i
logical::rebuilt
real(dp)::left(2400),right(2400)
native_mass(1)%flux=[(real(i,dp),i=1,96)]
temperature_closure(1)%slot=[(i,i=1,96)]
temperature_closure(1)%slot(2)=0
foreach_k: do k=1,2
 temperature_closure(1)%value(k,:)=[(-real(k*1000+i,dp),i=1,96)]
end do foreach_k
horiz_flux(2)%data(1)%elts=17.0_dp
horiz_flux(3)%data(1)%elts=29.0_dp
do i=1,2
 call scalar_allocate(block_scalar_tendency(i)%bdry,48,8,3,2,0,2,.true.,rebuilt)
end do
do mode=0,1
 validate_oracle=mode==1
 do full=0,1
  block_scalar_divergence_plan%full_transport=full==1
  do i=1,2
   call scalar_fill(block_scalar_tendency(i)%bdry,-991.0_dp)
  end do
  do field_level=1,2
   do v=1,2
    call fill_retained_boundary_record(1,0,2,4,v,field_level+1,2,.false.,1,8,2)
    call reference_boundary_record(1,0,2,4,v,field_level+1,2,.false.,2,8,2)
    left=scalar_read_range(block_scalar_tendency(1)%bdry,1,2400)
    right=scalar_read_range(block_scalar_tendency(2)%bdry,1,2400)
    if(any(abs(left-right)>0.0_dp))error stop 'narrow update changed full-record semantics'
   end do
  end do
 end do
end do
if(scalar_producer_work(6)/=16)error stop 'narrow coverage differs'
if(scalar_producer_work(7)/=1528)error stop 'logical slot count differs'
print *, 'PASS boundary update'
end program
