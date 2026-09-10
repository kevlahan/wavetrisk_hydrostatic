program test_mass_kernels
  use kind_mod, only : dp
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
  use parallel_block_mass_mod
  implicit none
  type(Mass_Restriction) :: program(2)
  type(Mass_Divergence) :: divergence(2)
  real(dp) :: flux(36),reference(36),rhs(12),small(4),partial,coarse,expected(12)
  integer :: s,e,n,i,mask,impulse
  do s=1,2
     program(s)%target=s
     do n=1,4
        do i=1,3
           program(s)%small(:,i,n)=modulo([n+i,n+2*i,n+3*i],12)+1
           program(s)%weight(i,1,n)=real(3*n+i,dp)/13.0_dp
           program(s)%weight(i,2,n)=real(n+2*i,dp)/17.0_dp
        end do
     end do
     do e=1,3
        program(s)%partial_flux(:,e)=modulo([1,2,3,4,5,6]+e,12)+1
        program(s)%partial_node(:,e)=[1,4,8,12]
        program(s)%coarse_node(:,e)=[2,3,5,7,9,11]
        program(s)%area(:,e)=[0.25_dp,0.75_dp]
        program(s)%overlap(:,e)=[1.1_dp,2.2_dp,3.3_dp,4.4_dp]
        program(s)%coarse(:,e)=[0.13_dp,0.27_dp,0.39_dp,0.41_dp,0.57_dp]
     end do
  end do
  rhs=[(real((-1)**i*i,dp)/7.0_dp,i=1,12)]
  ! Exercise every edge mask, asymmetric impulse inputs, and ordered writes:
  ! the second action reads edges overwritten by the first action.
  do mask=1,7
     do s=1,2
        program(s)%edge=[btest(mask,0),btest(mask,1),btest(mask,2)]
     end do
     do impulse=0,36
        flux=0.0_dp
        if(impulse==0) then
           flux=[(real((-1)**i*i,dp)/11.0_dp,i=1,36)]
        else
           flux(max(1,impulse))=13.0_dp/3.0_dp
        end if
        reference=flux
        do s=1,2
           associate(p=>program(s))
           do n=1,4
              if(mod(n,2)==1) then
                 small(n)=-sum(reference(p%small(:,1,n))*p%weight(:,1,n)) &
                      -sum((reference(p%small(:,2,n))-reference(p%small(:,3,n)))*p%weight(:,2,n))
              else
                 small(n)=sum(reference(p%small(:,1,n))*p%weight(:,1,n)) &
                      +sum((reference(p%small(:,2,n))-reference(p%small(:,3,n)))*p%weight(:,2,n))
              end if
           end do
           do e=1,3
              if(.not.p%edge(e)) cycle
              partial=sum(reference(p%partial_flux(1:2,e))*p%area(:,e)) &
                   -sum(reference(p%partial_flux(3:4,e)))*p%area(2,e) &
                   -sum(reference(p%partial_flux(5:6,e)))*p%area(1,e) &
                   +p%overlap(3,e)*rhs(8)-p%overlap(4,e)*rhs(12)-p%overlap(1,e)*rhs(1)+p%overlap(2,e)*rhs(4)
              coarse=p%coarse(1,e)*(rhs(3)-rhs(2))+p%coarse(2,e)*0.5_dp*(rhs(5)-rhs(2)) &
                   +p%coarse(3,e)*0.5_dp*(rhs(7)-rhs(2))+p%coarse(4,e)*0.5_dp*(rhs(3)-rhs(9)) &
                   +p%coarse(5,e)*0.5_dp*(rhs(3)-rhs(11))
              reference(3*(s-1)+e)=partial+coarse+small(e)+small(e+1)
           end do
           end associate
        end do
        call execute_mass_restriction(program,flux,rhs)
        if(any(transfer(flux,[0_int64],36)/=transfer(reference,[0_int64],36))) &
             error stop "mass restriction mask/impulse/order differs"
     end do
  end do
  divergence(1)%target=1
  divergence(1)%active=.true.
  divergence(1)%edge=[1,5,9,13,17,21]
  divergence(1)%inverse_area=0.125_dp
  divergence(2)%target=2
  divergence(2)%inverse_area=ieee_value(0.0_dp,ieee_quiet_nan)
  expected=rhs
  expected(1)=-(flux(1)-flux(5)+flux(9)-flux(13)+flux(17)-flux(21))*0.125_dp
  expected(2)=0.0_dp
  call execute_mass_divergence(divergence,flux,rhs)
  if(any(transfer(rhs,[0_int64],12)/=transfer(expected,[0_int64],12))) error stop "mass divergence differs"
  print *,"Native mass kernel tests PASS"
end program
