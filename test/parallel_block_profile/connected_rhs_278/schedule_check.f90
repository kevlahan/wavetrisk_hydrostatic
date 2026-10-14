module rhs_schedule_check_mod
 use iso_fortran_env, only : int64
 use shared_mod, only : N_BDRY,NORTH,EAST,SOUTH,WEST
 use patch_mod, only : PATCH_SIZE
 use parallel_block_mod, only : Block_Data
 use parallel_block_rhs_layout_mod
 use parallel_block_rhs_schedule_mod
 use rhs_schedule_oracle_mod
 implicit none
 contains
 subroutine check_schedule(block,generation,counts)
 type(Block_Data), intent(in) :: block
 integer(int64), intent(in) :: generation
 integer(int64), intent(out) :: counts(5)
 type(Rhs_Input_Layout) :: layout
 type(Rhs_Producer_Schedule) :: program
 integer :: first(0:N_BDRY),dims(2,N_BDRY),base(3),p,s,l,r,n,q,at,t,pass
 integer :: trace(7,(PATCH_SIZE+1)**2),link,record
 integer, allocatable :: order(:),fw(:,:),pvw(:,:),vw(:,:),expect_f(:,:),expect_v(:,:)
 logical :: boundary(4),sw((PATCH_SIZE+1)**2),rebuilt
 counts=0
 call compile_rhs_input_layout(block,generation,layout,rebuilt)
 base=[0,size(block%node),size(block%node)+size(block%bdry_node)]
 do l=minval(block%patch%level),maxval(block%patch%level)
   ! Both input orders verify the compiler never silently sorts patches. Actual
   ! geometry-owner level order must still be distributed during integration.
   do pass=1,2
     order=pack([(p,p=1,size(block%patch))],block%patch%level==l)
     if (pass==2) order=order(size(order):1:-1)
     call compile_rhs_producer_schedule(block,layout,generation,l,order,[integer::],[integer::],program)
     at=0
     if (allocated(fw)) deallocate(fw,pvw,vw,expect_f,expect_v)
     allocate(fw(3,program%nslot),pvw(3,program%nslot),vw(2,program%nslot),expect_f(3,program%nslot),expect_v(2,program%nslot))
     expect_f=0;expect_v=0
     do q=1,size(order)
       p=order(q);first(0)=block%patch(p)%elts_start;boundary=.false.
       do s=1,N_BDRY
         r=block%patch(p)%neigh(s)
         if (r>0) then
           first(s)=block%patch(r+1)%elts_start;dims(:,s)=PATCH_SIZE
         else
           link=-r;record=block%block_bdry(link)%ghost_id
           if (record>0) then
             first(s)=base(3)+block%ghost_storage(record)%local_start;dims(:,s)=PATCH_SIZE
           else
             record=block%block_bdry(link)%storage_id
             first(s)=base(2)+block%bdry_storage(record)%local_start;dims(:,s)=block%block_bdry(link)%dims
           end if
           if (s<=4) boundary(s)=block%block_bdry(link)%bdry_side>0
         end if
       end do
       call legacy_trace(first,dims,boundary,trace,sw,n)
       if (at+n>size(program%patch)) error stop 'RHS schedule missing events'
       if (any(program%node(:,at+1:at+n)/=trace(:,:n))) error stop 'RHS schedule legacy addresses differ'
       if (any(program%southwest(at+1:at+n).neqv.sw(:n))) error stop 'RHS schedule SW sequence differs'
       if (any(program%patch(at+1:at+n)/=p)) error stop 'RHS schedule patch sequence differs'
       do t=1,n
         r=at+t
         expect_f(:,trace(1,t))=r
         expect_v(1,trace(1,t))=r;expect_v(2,trace(1,t))=r
         expect_v(1,trace(5,t))=r;expect_v(2,trace(4,t))=r
         if (sw(t)) then
           expect_f(1,trace(5,t))=-r;expect_f(2,trace(7,t))=-r;expect_f(3,trace(4,t))=-r
           expect_v(1,trace(7,t))=-r;expect_v(2,trace(7,t))=-r
         end if
       end do
       at=at+n
     end do
     if (at/=size(program%patch)) error stop 'RHS schedule extra events'
     call rhs_schedule_last_writes(program,fw,pvw,vw)
     if (any(fw/=expect_f).or.any(pvw/=expect_f).or.any(vw/=expect_v)) error stop 'RHS schedule final local writer differs'
     if (pass==1) counts=counts+[int(at,int64),int(count(program%southwest),int64), &
          int(at-size(order)*PATCH_SIZE**2,int64),rhs_schedule_bytes(program),1_int64]
   end do
 end do
 end subroutine
end module
