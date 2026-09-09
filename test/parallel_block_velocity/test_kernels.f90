program test_velocity_kernels
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use kind_mod, only : dp
  use parallel_block_velocity_mod
  implicit none
  integer, parameter :: nodes(10,3) = reshape([ &
       1,1,5,6,7,7,8,2,2,2, &
       1,5,6,7,1,2,3,3,3,4, &
       5,6,7,1,1,4,4,4,9,5],[10,3])
  integer, parameter :: edges(10,3) = reshape([ &
       2,3,1,2,3,2,3,1,2,3, &
       3,1,2,3,1,3,1,2,3,1, &
       1,2,3,1,2,1,2,3,1,2],[10,3])
  real(dp) :: flux(3,9),pv(3,9),weights(5,2,3),parts(6,4),area(4)
  real(dp) :: result(3),expected(3),direct(3),child(3),neighbor(3),poison
  real(dp) :: p,upper,exner,phi,lower,lower_phi,mass,temperature,reference_p,reference_e
  real(dp) :: cumulative,reference_w
  logical :: selected(3)
  integer :: n,e,out,term,side,t,node,offset,mask,k
  type(Velocity_Level_Program) :: program
  type(Velocity_Gradient_Stencil) :: gradient_plan(2)
  real(dp) :: source_workspace(3,9),physics_workspace(3,9),gradient_workspace(3,9)
  real(dp) :: bfield(9),efield(9),mfield(9),tfield(9)

  do n = 1,9
     do e = 1,3
        pv(e,n) = real(3*n+e,dp)/7.0_dp
     end do
  end do
  do out = 1,3
     do side = 1,2
        do t = 1,5
           weights(t,side,out) = real(t+10*side+100*out,dp)/17.0_dp
        end do
     end do
  end do
  ! Every input edge in the 3x3 neighborhood, not only a symmetric stencil.
  do n = 1,9
     do e = 1,3
        flux = 0.0_dp
        flux(e,n) = 13.0_dp/3.0_dp
        expected = 0.0_dp
        do out = 1,3
           do term = 1,10
              if (nodes(term,out) /= n .or. edges(term,out) /= e) cycle
              side = (term-1)/5+1
              t = modulo(term-1,5)+1
              expected(out) = flux(e,n)*(0.5_dp*(pv(e,n)+pv(out,1)))*weights(t,side,out)
           end do
        end do
        result = velocity_qperp(flux,pv,weights)
        call exact(result,expected,"Qperp impulse stencil")
     end do
  end do
  do n = 1,9
     do e = 1,3
        flux(e,n) = real((-1)**n*(3*n+e),dp)/11.0_dp
     end do
  end do
  expected = 0.0_dp
  do out = 1,3
     do term = 1,10
        n = nodes(term,out)
        e = edges(term,out)
        side = (term-1)/5+1
        t = modulo(term-1,5)+1
        expected(out) = expected(out)+flux(e,n)*(0.5_dp*(pv(e,n)+pv(out,1)))*weights(t,side,out)
     end do
  end do
  result = velocity_qperp(flux,pv,weights)
  call exact(result,expected,"dense Qperp term order and cancellation")

  do n = 1,4
     do t = 1,6
        parts(t,n) = real(n*10+t,dp)/13.0_dp
     end do
     area(n) = 1.0_dp/sum(parts(:,n))
  end do
  weights = velocity_weights(parts,area)
  do out = 1,3
     do side = 1,2
        node = 1
        offset = out-1
        if (side == 2) then
           node = out+1
           offset = out+2
        end if
        cumulative = 0.0_dp
        do term = 1,5
           cumulative = cumulative+parts(modulo(offset+term-1,6)+1,node)
           reference_w = 0.5_dp-cumulative*area(node)
           if (modulo(term,2) == 0) reference_w = -reference_w
           call exact([weights(term,side,out),0.0_dp,0.0_dp],[reference_w,0.0_dp,0.0_dp],"Qperp geometry")
        end do
     end do
  end do

  result = velocity_source([1.0_dp,2.0_dp,3.0_dp],[4.0_dp,5.0_dp,6.0_dp],[7.0_dp,8.0_dp,9.0_dp],.true.)
  call exact(result,[27.0_dp,38.0_dp,51.0_dp],"direct source")
  poison = ieee_value(0.0_dp,ieee_quiet_nan)
  do mask = 0,7
     direct = poison
     child = poison
     neighbor = poison
     do e = 1,3
        selected(e) = btest(mask,e-1)
        if (selected(e)) then
           child(e) = real(e,dp)
           neighbor(e) = real(10*e,dp)
           expected(e) = real(11*e,dp)
        else
           direct(e) = real(100*e,dp)
           expected(e) = direct(e)
        end if
     end do
     result = velocity_restrict_source(direct,child,neighbor,selected)
     call exact(result,expected,"mixed child masks and unused poisoned operands")
  end do

  result = velocity_gradient([40.0_dp,80.0_dp,160.0_dp],[2.0_dp,4.0_dp,8.0_dp], &
       [1.0_dp,5.0_dp,9.0_dp,17.0_dp],[2.0_dp,6.0_dp,10.0_dp,18.0_dp], &
       [1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2.0_dp,4.0_dp,6.0_dp,8.0_dp],.true.)
  call exact(result,[14.0_dp,26.0_dp,14.0_dp],"gradient orientation and metric")
  result = velocity_gradient(spread(poison,1,3),spread(poison,1,3), &
       spread(poison,1,4),spread(poison,1,4),spread(poison,1,4),spread(poison,1,4),.false.)
  call exact(result,[0.0_dp,0.0_dp,0.0_dp],"inactive gradient must not evaluate inputs")
  result = velocity_source(spread(poison,1,3),spread(poison,1,3),spread(poison,1,3),.false.)
  call exact(result,[0.0_dp,0.0_dp,0.0_dp],"inactive source must not evaluate inputs")

  lower = 100000.0_dp
  lower_phi = 123.0_dp
  do k = 1,30
     mass = 270.0_dp+real(k,dp)/7.0_dp
     temperature = mass*(280.0_dp+real(k,dp))
     reference_p = 0.5_dp*(lower+(lower-9.81_dp*mass))
     reference_e = 1004.0_dp*(reference_p/100000.0_dp)**(2.0_dp/7.0_dp)
     call velocity_pressure_layer(mass,temperature,9.81_dp,1004.0_dp,100000.0_dp,2.0_dp/7.0_dp, &
          lower,lower_phi,p,upper,exner,phi)
     call exact([p,exner,phi], &
          [reference_p,reference_e,lower_phi+9.81_dp*(2.0_dp/7.0_dp)*temperature*reference_e/reference_p], &
          "ordered physical-layer integration")
     lower = upper
     lower_phi = phi
  end do
  ! The instruction tape must preserve mixed direct/restricted overwrite
  ! order. Regrouping all direct operations before restriction is incorrect.
  allocate(program%action(5),program%operand(5),program%direct(2),program%restriction(2))
  program%action=[VELOCITY_DIRECT,VELOCITY_RESTRICT,VELOCITY_DIRECT,VELOCITY_RESTRICT,VELOCITY_DIRECT]
  program%operand=[1,1,1,2,2]
  program%direct(1)%node=[1,2,3,4,5,6,7,8,9]
  program%direct(1)%active=.true.
  program%direct(1)%weights=weights
  program%direct(1)%length=[2.0_dp,3.0_dp,4.0_dp]
  program%direct(2)%node(1)=9
  ! Inactive stencil has deliberately invalid unused addresses.
  program%direct(2)%active=.false.
  program%restriction(1)%target=1
  program%restriction(1)%child=2
  program%restriction(1)%neighbor=[3,0,0]
  program%restriction(1)%edge=[.true.,.false.,.false.]
  program%restriction(2)%target=1
  program%restriction(2)%child=2
  program%restriction(2)%neighbor=[0,3,0]
  program%restriction(2)%edge=[.false.,.true.,.false.]
  call validate_velocity_program(program,9)
  source_workspace=13.0_dp
  physics_workspace=2.0_dp
  expected=velocity_source(velocity_qperp(flux,pv,weights),physics_workspace(:,1),program%direct(1)%length,.true.)
  expected(2)=26.0_dp
  call execute_velocity_sources(program,flux,pv,physics_workspace,source_workspace)
  call exact(source_workspace(:,1),expected,"ordered source/restriction program")
  call exact(source_workspace(:,9),[0.0_dp,0.0_dp,0.0_dp],"inactive source program")
  gradient_plan(1)%node=[1,2,3,4]
  gradient_plan(1)%active=.true.
  gradient_plan(1)%length=program%direct(1)%length
  gradient_plan(2)%node(1)=9
  bfield=[(real(n*n,dp),n=1,9)]
  efield=2.0_dp*bfield
  mfield=3.0_dp
  tfield=7.0_dp
  gradient_workspace=poison
  expected=velocity_gradient(source_workspace(:,1),gradient_plan(1)%length, &
       bfield(1:4),efield(1:4),mfield(1:4),tfield(1:4),.true.)
  call execute_velocity_gradients(gradient_plan,bfield,efield,mfield,tfield,source_workspace,gradient_workspace)
  call exact(gradient_workspace(:,1),expected,"gradient program")
  call exact(gradient_workspace(:,9),[0.0_dp,0.0_dp,0.0_dp],"inactive gradient program")
  print '(a)',"Stage 175 velocity numeric kernels and programs: PASS"
contains
  subroutine exact(value,reference,label)
    real(dp), intent(in) :: value(3),reference(3)
    character(*), intent(in) :: label
    if (all(transfer(value,[0_int64],3) == transfer(reference,[0_int64],3))) return
    print *,label
    print *,value
    print *,reference
    error stop "kernel check is not bitwise equal"
  end subroutine exact
end program test_velocity_kernels
