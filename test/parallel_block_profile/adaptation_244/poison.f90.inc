  subroutine poison_adaptation_fields(restore)
    logical,intent(in)::restore
    integer::d,k,v,n
    real(dp)::poison
    poison=ieee_value(0.0_dp,ieee_quiet_nan)
    if(.not.restore) allocate(adapt_backup(N_VARIABLE,zmin:zmax,size(grid)))
    do d=1,size(grid)
       do k=zmin,zmax
          do v=1,N_VARIABLE
             associate(a=>adapt_backup(v,k,d))
               if(restore) then
                  n=size(a%s)
                  if(size(sol(v,k)%data(d)%elts)<n.or.size(wav_coeff(v,k)%data(d)%elts)<n) &
                       error stop 'adaptation poison storage shrank'
                  sol(v,k)%data(d)%elts(:n)=a%s
                  wav_coeff(v,k)%data(d)%elts(:n)=a%w
               else
                  n=grid(d)%node%length
                  if(v==S_VELO) n=EDGE*n
                  allocate(a%s,source=sol(v,k)%data(d)%elts(:n))
                  allocate(a%w,source=wav_coeff(v,k)%data(d)%elts(:n))
                  sol(v,k)%data(d)%elts=poison
                  wav_coeff(v,k)%data(d)%elts=poison
               end if
             end associate
          end do
       end do
    end do
    if(restore) then
       deallocate(adapt_backup)
       write(*,'(a,i0)') 'Native adaptation poison restored rank=',rank
    end if
  end subroutine poison_adaptation_fields
