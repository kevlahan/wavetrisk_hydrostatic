module fixture
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, parameter :: zlevels=4, EDGE=3, N_BDRY=8, ADJZONE=3, S_VELO=1, S_MASS=2, S_TEMP=3
  integer, parameter :: RT=0, DG=1, UP=2
  real(dp), parameter :: p_top=100.0_dp, grav_accel=9.81_dp
  real(dp) :: a_vert(0:zlevels), b_vert(0:zlevels)
  type :: Array
     real(dp) :: elts(120)
  end type
  type :: Float_Field
     type(Array) :: data(1)
  end type
  type :: Mask
     integer :: elts(40)=ADJZONE
  end type
  type :: Domain
     integer :: id=0
     type(Mask) :: mask_n
  end type
  type(Float_Field) :: sol(3,zlevels), sol_mean(3,zlevels), old_mass(zlevels)
contains
  integer function idx(i,j,offs,dims)
    integer, intent(in) :: i,j,offs(:),dims(:,:)
    idx=i+4*j
  end function
  subroutine interpolate(n,new,p_new,old,p_old)
    integer, intent(in) :: n
    real(dp), intent(in) :: p_new(0:n),old(n),p_old(0:n)
    real(dp), intent(out) :: new(n)
    ! Deterministic test interpolation sensitive to both coordinate arrays.
    ! This test checks the caller contract, not the production interpolator.
    new=old+(p_new(1:n)-p_old(1:n))*1.0e-5_dp
  end subroutine
  include 'polar_remap.inc'
end module

program check
  use fixture
  implicit none
  type(Domain) :: dom
  type(Float_Field) :: initial(3,zlevels), reference(3,zlevels)
  integer :: k,v,n,offs(N_BDRY+1),dims(2,N_BDRY+1),mode,pole,i,j
  do k=0,zlevels
     a_vert(k)=p_top*real(k,dp)/zlevels
     b_vert(k)=1.0_dp-real(k,dp)/zlevels
  end do
  do k=1,zlevels
     do v=1,3
        do n=1,120
           sol(v,k)%data(1)%elts(n)=real(11*k+3*v+n,dp)/7
           sol_mean(v,k)%data(1)%elts(n)=real(19*k+v+n,dp)*13
        end do
     end do
  end do
  offs=0
  dims=4
  initial=sol
  do pole=1,2
   i=merge(0,4,pole==1)
   j=merge(4,0,pole==1)
   do mode=1,2
     dom%mask_n%elts=ADJZONE
     if (mode==2) dom%mask_n%elts=ADJZONE-1
     sol=initial
     old_mass=sol(S_MASS,:)
     call remap_compressible(dom,0,i,j,0,offs,dims,1)
     reference=sol
     sol=initial
     call remap_compressible_pole(dom,0,i,j,0,offs,dims,1)
     do k=1,zlevels
        do v=1,3
           if (any(sol(v,k)%data(1)%elts/=reference(v,k)%data(1)%elts)) &
                error stop 'polar remap differs from legacy scalar contract'
        end do
     end do
   end do
  end do
  print *, 'polar remap contract passed'
end program
