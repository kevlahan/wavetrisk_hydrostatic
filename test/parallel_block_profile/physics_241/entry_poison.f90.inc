  subroutine poison_native_entry
    integer::d,k,id,key(4),nn
    allocate(entry_backup(size(grid)))
    do d=1,size(grid)
       nn=size(sol(S_MASS,1)%data(d)%elts)
       associate(a=>entry_backup(d))
         allocate(a%mass(nn,zlevels),a%temp(nn,zlevels),a%velo(EDGE*nn,zlevels), &
              a%mean_m(nn,zlevels),a%mean_t(nn,zlevels),a%soil(nn,zmin:0))
         do k=1,zlevels
            a%mass(:,k)=sol(S_MASS,k)%data(d)%elts
            a%temp(:,k)=sol(S_TEMP,k)%data(d)%elts
            a%velo(:,k)=sol(S_VELO,k)%data(d)%elts
            a%mean_m(:,k)=sol_mean(S_MASS,k)%data(d)%elts
            a%mean_t(:,k)=sol_mean(S_TEMP,k)%data(d)%elts
            sol(S_MASS,k)%data(d)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
            sol(S_TEMP,k)%data(d)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
            sol(S_VELO,k)%data(d)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
         end do
         do k=zmin,0
            a%soil(:,k)=sol(S_TEMP,k)%data(d)%elts
         end do
         do id=1,grid(d)%node%length
            call inverse_operand_key(d,id-1,key)
            if(key(1)==0) cycle
            do k=1,zlevels
               sol_mean(S_MASS,k)%data(d)%elts(id)=ieee_value(0.0_dp,ieee_quiet_nan)
               sol_mean(S_TEMP,k)%data(d)%elts(id)=ieee_value(0.0_dp,ieee_quiet_nan)
            end do
            do k=zmin,0
               sol(S_TEMP,k)%data(d)%elts(id)=ieee_value(0.0_dp,ieee_quiet_nan)
            end do
         end do
       end associate
    end do
  end subroutine poison_native_entry

  subroutine restore_native_entry
    integer::d,k
    do d=1,size(grid)
       associate(a=>entry_backup(d))
         do k=1,zlevels
            sol(S_MASS,k)%data(d)%elts=a%mass(:,k)
            sol(S_TEMP,k)%data(d)%elts=a%temp(:,k)
            sol(S_VELO,k)%data(d)%elts=a%velo(:,k)
            sol_mean(S_MASS,k)%data(d)%elts=a%mean_m(:,k)
            sol_mean(S_TEMP,k)%data(d)%elts=a%mean_t(:,k)
         end do
         do k=zmin,0
            sol(S_TEMP,k)%data(d)%elts=a%soil(:,k)
         end do
       end associate
    end do
    deallocate(entry_backup)
    write(*,'(a)') 'Native physics entry poison restored after gather'
    write(*,'(a,i0,a,i0)') 'Native physics primary aliases rank=',rank, &
         ' columns=',count(batch_entry_key(1,1,:)==0)
  end subroutine restore_native_entry
