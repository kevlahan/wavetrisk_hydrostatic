module parallel_block_mass_mod
  ! Numerical storage and ordered operators, independent of Domain and MPI.
  use kind_mod, only : dp
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  public :: Mass_Restriction, Mass_Divergence, Mass_Level, Native_Mass_Workspace
  public :: native_mass, mass_transaction, mass_oracle, execute_mass_restriction, execute_mass_divergence
  public :: mass_work

  type :: Mass_Restriction
     integer :: target=0
     logical :: edge(3)=.false.
     integer :: small(3,3,4)=0, partial_flux(6,3)=0, partial_node(4,3)=0, coarse_node(6,3)=0
     real(dp) :: weight(3,2,4)=0.0_dp, area(2,3)=0.0_dp, overlap(4,3)=0.0_dp, coarse(5,3)=0.0_dp
  end type
  type :: Mass_Divergence
     integer :: target=0, edge(6)=0
     real(dp) :: inverse_area=0.0_dp
     logical :: active=.false.
  end type
  type :: Mass_Level
     type(Mass_Restriction), allocatable :: restriction(:)
     type(Mass_Divergence), allocatable :: divergence(:)
  end type
  type :: Native_Mass_Workspace
     integer(int64) :: generation=-1_int64
     type(Mass_Level), allocatable :: level(:)
     real(dp), allocatable :: flux(:), tendency(:,:), rk(:,:)
     logical, allocatable :: active(:), restricted(:), ready(:)
  end type
  type(Native_Mass_Workspace), allocatable, target, save :: native_mass(:)
  logical, save :: mass_transaction=.false., mass_oracle=.false.
  ! plan builds, restricted edges, divergences, boundary values, published RHS, RK values
  integer(int64), save :: mass_work(6)=0_int64
contains
  subroutine execute_mass_restriction(program,flux,rhs)
    type(Mass_Restriction), intent(in) :: program(:)
    real(dp), intent(inout) :: flux(:)
    real(dp), intent(in) :: rhs(:)
    real(dp) :: small(4),partial,coarse
    integer :: s,e,n
    do s=1,size(program)
       associate(p=>program(s))
       do n=1,4
          small(n)=sum(flux(p%small(:,1,n))*p%weight(:,1,n)) + &
               sum((flux(p%small(:,2,n))-flux(p%small(:,3,n)))*p%weight(:,2,n))
          if (mod(n,2)==1) small(n)=-small(n)
       end do
       do e=1,3
          if (.not. p%edge(e)) cycle
          partial=sum(flux(p%partial_flux(1:2,e))*p%area(:,e)) &
               -sum(flux(p%partial_flux(3:4,e)))*p%area(2,e) &
               -sum(flux(p%partial_flux(5:6,e)))*p%area(1,e) &
               +p%overlap(3,e)*rhs(p%partial_node(3,e))-p%overlap(4,e)*rhs(p%partial_node(4,e)) &
               -p%overlap(1,e)*rhs(p%partial_node(1,e))+p%overlap(2,e)*rhs(p%partial_node(2,e))
          coarse=p%coarse(1,e)*(rhs(p%coarse_node(2,e))-rhs(p%coarse_node(1,e))) &
               +p%coarse(2,e)*0.5_dp*(rhs(p%coarse_node(3,e))-rhs(p%coarse_node(1,e))) &
               +p%coarse(3,e)*0.5_dp*(rhs(p%coarse_node(4,e))-rhs(p%coarse_node(1,e))) &
               +p%coarse(4,e)*0.5_dp*(rhs(p%coarse_node(2,e))-rhs(p%coarse_node(5,e))) &
               +p%coarse(5,e)*0.5_dp*(rhs(p%coarse_node(2,e))-rhs(p%coarse_node(6,e)))
          flux(3*(p%target-1)+e)=partial+coarse+small(e)+small(e+1)
       end do
       mass_work(2)=mass_work(2)+count(p%edge)
       end associate
    end do
  end subroutine

  subroutine execute_mass_divergence(program,flux,rhs)
    type(Mass_Divergence), intent(in) :: program(:)
    real(dp), intent(in) :: flux(:)
    real(dp), intent(inout) :: rhs(:)
    integer :: s
    do s=1,size(program)
       associate(p=>program(s))
       rhs(p%target)=0.0_dp
       if (p%active) rhs(p%target)=-(flux(p%edge(1))-flux(p%edge(2))+flux(p%edge(3))-flux(p%edge(4)) &
            +flux(p%edge(5))-flux(p%edge(6)))*p%inverse_area
       end associate
    end do
    mass_work(3)=mass_work(3)+size(program)
  end subroutine
end module parallel_block_mass_mod
