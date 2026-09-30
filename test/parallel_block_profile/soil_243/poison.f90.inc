  subroutine poison_soil_fields(restore)
    logical,intent(in)::restore
    integer::d,k,v
    real(dp)::poison
    poison=ieee_value(0.0_dp,ieee_quiet_nan)
    if(.not.restore) allocate(soil_backup(N_VARIABLE,zmin:0,size(grid)))
    do d=1,size(grid)
       do k=zmin,0
          do v=1,N_VARIABLE
             associate(a=>soil_backup(v,k,d))
               if(restore) then
                  if(size(sol(v,k)%data(d)%elts)/=size(a%s).or. &
                       size(wav_coeff(v,k)%data(d)%elts)/=size(a%w)) error stop 'soil poison allocation changed'
                  sol(v,k)%data(d)%elts(:)=a%s
                  wav_coeff(v,k)%data(d)%elts(:)=a%w
               else
                  allocate(a%s,source=sol(v,k)%data(d)%elts)
                  allocate(a%w,source=wav_coeff(v,k)%data(d)%elts)
                  if(size(a%s)/=size(a%w)) write(*,'(a,5i10)') &
                       'Native soil poison independent capacities rank/d/k/s/w=',rank,d,k,size(a%s),size(a%w)
                  sol(v,k)%data(d)%elts=poison
                  wav_coeff(v,k)%data(d)%elts=poison
               end if
             end associate
          end do
       end do
    end do
    if(restore) then
       deallocate(soil_backup)
       write(*,'(a,i0)') 'Native soil Domain poison restored after transform rank=',rank
    end if
  end subroutine poison_soil_fields
