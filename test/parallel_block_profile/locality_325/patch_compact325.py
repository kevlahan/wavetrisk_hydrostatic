"""Exact bounded real digests and physical-key hooks for Stage325 qualification."""
from pathlib import Path
from patch_owner325 import replace_once


def apply(tree):
    root=Path(tree)/'src'
    p=root/'audit_compress322.c'
    p.write_text(p.read_text().replace(' * The compressed stream contains the original Stage316 bytes without changes.\n * The caller enforces its uncompressed512MiB bound and lifecycle.', ' * Stage325 writes exact metadata and finite real-value digests.\n * The Fortran caller enforces the uncompressed128MiB bound and lifecycle.')+'\n'+Path(__file__).with_name('digest325.c').read_text())
    p=root/'audit316.f90';s=p.read_text()
    s=replace_once(s,'c_loc,c_null_char','c_loc,c_null_char,c_null_ptr')
    s=replace_once(s,'  ! Qualification-only component records, losslessly compressed in this process.',
                  '  ! Qualification records: exact integer metadata and finite real SHA-256 digests.')
    s=replace_once(s,'limit316=512_int64*1024_int64*1024_int64','limit316=128_int64*1024_int64*1024_int64')
    s=replace_once(s,'Stage316 audit file exceeds 512MiB bound','Stage325 compact audit file exceeds 128MiB bound')
    s=replace_once(s,"    magic='RK316V01'","    magic='RK325V01'")
    s=replace_once(s,'    if(.not.enabled316) return\n',
                  "    if(.not.enabled316) return\n"
                  "    option='';call get_environment_variable('WAVETRISK_COMPACT_AUDIT325',option,status=status)\n"
                  "    if(status/=0.or.trim(option)/='1') error stop 'Stage325 audit requires compact mode1'\n")
    s=replace_once(s,'  interface\n', '''  interface
     integer(c_int) function digest_write325(address,length) bind(C,name='audit_digest325_write')
       import :: c_int,c_ptr,c_int64_t
       type(c_ptr),value :: address
       integer(c_int64_t),value :: length
     end function
''')
    s=replace_once(s,'    added=16_int64+k+n*width\n',
                  '    added=16_int64+k+n*width\n    if(kind==13) added=16_int64+k+32_int64\n')
    start=s.index('  subroutine audit316_r(');end=s.index('  subroutine audit316_l(',start)
    s=s[:start]+'''  subroutine audit316_r(name,value)
    character(*),intent(in) :: name
    real(real64),intent(in),target :: value(:)
    real(real64),allocatable,target :: packed(:)
    integer(c_int) :: status
    if(.not.audit316_enabled()) return
    ! Kind13 preserves the original logical element count; payload is32 bytes.
    call header316(name,13,int(size(value),int64),0)
    if(size(value)==0) then
       status=digest_write325(c_null_ptr,0_c_int64_t)
    else if(is_contiguous(value)) then
       status=digest_write325(c_loc(value(1)),8_c_int64_t*size(value))
    else
       allocate(packed(size(value)));packed=value
       status=digest_write325(c_loc(packed(1)),8_c_int64_t*size(packed))
    end if
    if(status/=0_c_int) error stop 'Stage325 nonfinite real or digest write failed'
  end subroutine

'''+s[end:]
    p.write_text(s)
    p=root/'parallel_block_rhs_backend.f90';s=p.read_text()
    s=replace_once(s,"       call audit316_i('rk.adapter.shape',shape(adapter_value))\n",'''       call audit316_i('rk.adapter.shape',shape(adapter_value))
       call audit316_i('rk.adapter.numeric.shape',shape(adapter_value(:n,:)))
       call audit316_i('rk.adapter.domain.shape',shape(adapter_value(n+1:,:)))
       call audit316_r('rk.adapter.numeric',reshape(adapter_value(:n,:),[size(adapter_value(:n,:))]))
       call audit316_r('rk.adapter.domain',reshape(adapter_value(n+1:,:),[size(adapter_value(n+1:,:))]))
''')
    p.write_text(s)
    p=root/'parallel_block_transform_backend.f90';s=p.read_text()
    s=replace_once(s,'  integer,save :: nnode=0,nrow=0,nfixed=0,npublication=0,nboundary=0\n',
                  '  integer,save :: nnode=0,nrow=0,nfixed=0,npublication=0,nboundary=0\n'
                  '  integer,allocatable,save :: audit_publication_key325(:,:)\n')
    s=replace_once(s,'    placement=Transform_Placement();generation=-1\n',
                  '    if(allocated(audit_publication_key325)) deallocate(audit_publication_key325)\n'
                  '    placement=Transform_Placement();generation=-1\n')
    s=replace_once(s,'    pkeys(:,:np)=physical(:,selected);pkeys(:,np+1:)=inventory%physical\n',
                  '    pkeys(:,:np)=physical(:,selected);pkeys(:,np+1:)=inventory%physical\n'
                  '    if(audit316_enabled()) audit_publication_key325=pkeys\n')
    s=replace_once(s,"    call audit316_i('extent',[nnode,nrow,nfixed,npublication,nboundary,level_start,level_end,zlevels])\n",
                  "    call audit316_i('extent',[nnode,nrow,nfixed,npublication,nboundary,level_start,level_end,zlevels])\n"
                  "    call audit316_i('publication.key325.shape',shape(audit_publication_key325))\n"
                  "    call audit316_i('publication.key325',reshape(audit_publication_key325,[size(audit_publication_key325)]))\n")
    s=replace_once(s,"    call audit316_i('transform.scalar.shape',shape(scalar))\n",'''    call audit316_i('transform.scalar.shape',shape(scalar))
    call audit316_i('transform.numeric.scalar.shape',shape(scalar(:nnode,:)))
    call audit316_i('transform.numeric.vector.shape',shape(vector(:3*nnode,:)))
    call audit316_r('transform.numeric.scalar',reshape(scalar(:nnode,:),[size(scalar(:nnode,:))]))
    call audit316_r('transform.numeric.vector',reshape(vector(:3*nnode,:),[size(vector(:3*nnode,:))]))
    call audit316_i('transform.numeric.saved_scalar.shape',shape(scalar(nnode+1:nnode+size(save_node)+size(fixed_node),:)))
    call audit316_i('transform.numeric.saved_vector.shape', &
         shape(vector(3*nnode+1:3*(nnode+size(save_node)+size(fixed_node)),:)))
    call audit316_r('transform.numeric.saved_scalar', &
         reshape(scalar(nnode+1:nnode+size(save_node)+size(fixed_node),:), &
         [(size(save_node)+size(fixed_node))*size(scalar,2)]))
    call audit316_r('transform.numeric.saved_vector', &
         reshape(vector(3*nnode+1:3*(nnode+size(save_node)+size(fixed_node)),:), &
         [3*(size(save_node)+size(fixed_node))*size(vector,2)]))
''')
    p.write_text(s)
    p=root/'parallel_block_rhs_state.f90';s=p.read_text()
    s=replace_once(s,'  use audit316_mod, only : audit316_i,audit316_i64,audit316_l',
                  '  use audit316_mod, only : audit316_i,audit316_i64,audit316_l,audit316_enabled')
    s=replace_once(s,'     integer :: nprovider=0,nslot=0\n',
                  '     integer :: nprovider=0,nslot=0\n'
                  '     integer,allocatable :: audit_provider_key325(:,:),audit_physical_key325(:,:)\n')
    s=replace_once(s,"    call audit316_i('input.extent',[input%nprovider,input%nslot])\n",'''    call audit316_i('input.extent',[input%nprovider,input%nslot])
    call audit316_i('input.provider_key325.shape',shape(input%audit_provider_key325))
    call audit316_i('input.provider_key325',reshape(input%audit_provider_key325,[size(input%audit_provider_key325)]))
    call audit316_i('input.physical_key325.shape',shape(input%audit_physical_key325))
    call audit316_i('input.physical_key325',reshape(input%audit_physical_key325,[size(input%audit_physical_key325)]))
''')
    s=replace_once(s,'    input%nprovider=np;input%nslot=nd;input%topology=topology;input%masks=masks\n',
                  '    input%nprovider=np;input%nslot=nd;input%topology=topology;input%masks=masks\n'
                  '    if(audit316_enabled()) then\n'
                  '       input%audit_provider_key325=provider_key;input%audit_physical_key325=physical_key\n'
                  '    end if\n')
    p.write_text(s)


if __name__=='__main__':
    import sys
    apply(sys.argv[1])
