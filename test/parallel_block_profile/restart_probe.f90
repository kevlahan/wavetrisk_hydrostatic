module restart_probe_mod
  ! Isolated diagnostic builds only. Read-only Domain snapshots; no boundary
  ! refresh, reduction, mask update, pointer rebinding or production writeback.
  use iso_fortran_env, only : int64
  use shared_mod, only : N_VARIABLE, MULT, zmin, zmax, time, istep_cumul, iremap, threshold
  use arch_mod, only : rank, glo_id
  use domain_mod, only : grid, sol, wav_coeff
  use patch_mod, only : PATCH_SIZE
  implicit none
contains
  subroutine restart_probe(tag)
    character(*), intent(in) :: tag
    integer :: u,d,p,c,q,k,v,b,n,status,selected,nd
    integer(int64), allocatable :: key(:)
    character(256) :: filename
    character(16) :: enabled
    call get_environment_variable('WAVETRISK_RESTART_PROBE',enabled,status=status)
    if (status/=0.or.trim(enabled)/='1') return
    call get_environment_variable('WAVETRISK_PROBE_DOMAIN',enabled,status=status)
    selected=-1
    if (status==0.and.len_trim(enabled)>0) read(enabled,*) selected
    nd=size(grid)
    if (selected>=0) nd=count(glo_id(rank+1,1:size(grid))==selected)
    write(filename,'(a,a,a,i0,a,i0,a)') 'probe-',tag,'-step-',istep_cumul,'-rank-',rank,'.bin'
    open(newunit=u,file=trim(filename),status='new',access='stream',form='unformatted',action='write')
    write(u) 179401,nd,N_VARIABLE,zmin,zmax,PATCH_SIZE,istep_cumul,iremap
    write(u) time,threshold
    do d=1,size(grid)
       if (selected>=0.and.glo_id(rank+1,d)/=selected) cycle
       allocate(key(grid(d)%patch%length))
       key=0_int64
       key(2)=1_int64
       do p=2,size(key)
          if (key(p)==0_int64) error stop 'restart probe: patch ordering/key missing'
          do c=1,4
             q=grid(d)%patch%elts(p)%children(c)
             if (q<=0) cycle
             if (q+1>size(key)) error stop 'restart probe: invalid child'
             key(q+1)=4_int64*key(p)+int(c-1,int64)
          end do
       end do
       write(u) glo_id(rank+1,d),size(key)-1
       do p=2,size(key)
          b=grid(d)%patch%elts(p)%elts_start
          n=PATCH_SIZE**2
          write(u) key(p),grid(d)%patch%elts(p)%level
          write(u) grid(d)%mask_n%elts(b+1:b+n),grid(d)%mask_e%elts(3*b+1:3*(b+n))
          do k=zmin,zmax
             do v=1,N_VARIABLE
                write(u) sol(v,k)%data(d)%elts(MULT(v)*b+1:MULT(v)*(b+n))
                write(u) wav_coeff(v,k)%data(d)%elts(MULT(v)*b+1:MULT(v)*(b+n))
             end do
          end do
       end do
       deallocate(key)
    end do
    close(u)
    if (tag=='after-dynamics'.and.selected<0) call alias_probe(tag)
  end subroutine

  subroutine alias_probe(tag)
    ! First-step diagnosis only: full raw Domain solution, including aliases.
    ! Compare raw addresses only when the two Domain layouts are identical.
    character(*), intent(in) :: tag
    character(256) :: filename
    integer :: u,d,k,v,n
    write(filename,'(a,a,a,i0,a,i0,a)') 'aliases-',tag,'-step-',istep_cumul,'-rank-',rank,'.bin'
    open(newunit=u,file=trim(filename),status='new',access='stream',form='unformatted',action='write')
    write(u) 179402,size(grid),N_VARIABLE,zmin,zmax,istep_cumul
    do d=1,size(grid)
       n=size(sol(2,zmin)%data(d)%elts)
       write(u) glo_id(rank+1,d),n
       do k=zmin,zmax
          do v=1,N_VARIABLE
             write(u) sol(v,k)%data(d)%elts
          end do
       end do
    end do
    close(u)
  end subroutine alias_probe
end module
