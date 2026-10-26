"""Compile the actual selector/codec and verify exact maps and byte contracts."""
from pathlib import Path
import ctypes
import gzip
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from patch_owner325 import maps


def run_gate(source, out):
    source=Path(source).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    fc=shutil.which('gfortran');cc=shutil.which('cc')
    if not fc or not cc:raise RuntimeError('C and Fortran compilers required')
    a,b,p=maps(Path(__file__).parent)
    text=(source/'src/parallel_block_mpi.f90').read_text()
    start=text.index('  integer function owner_bridge322_mode()')
    end=text.index('  subroutine check_owner_bridge322_catalog(',start)
    functions=text[start:end]
    harness='module selector_gate325\nimplicit none\ncontains\n'+functions+'''end module
program gate
use selector_gate325
implicit none
integer :: i,m
m=owner_bridge322_mode()
if(m==4) call check_locality325_permutation()
print *,m
if(m>0) then
do i=1,640
print *,selected_owner322(i,m)
end do
end if
end program
'''
    (out/'selector.f90').write_text(harness)
    subprocess.run([fc,'-std=f2008','-Wall','-Wextra','-Werror','-fcheck=all','-ffree-line-length-132','selector.f90','-o','selector'],cwd=out,check=True)
    for selector,mode,expected in [('off',0,[]),('A160',1,[i//4 for i in range(640)]),('B256',2,[i//4 for i in range(640)]),('C256',3,a),('L256',4,b)]:
        env={**os.environ,'WAVETRISK_OWNER_MAP322':selector}
        result=subprocess.run([str(out/'selector')],env=env,capture_output=True,text=True,check=True)
        assert list(map(int,result.stdout.split()))==[mode]+expected,selector
    result=subprocess.run([str(out/'selector')],env={**os.environ,'WAVETRISK_OWNER_MAP322':'unsealed'},capture_output=True)
    assert result.returncode!=0
    # No existing fixed-case guard is relaxed by the new selector.
    for marker in ['N_GLO_DOMAIN/=160','min_level/=5','max_level/=5','resume<0','zlevels/=30',
                   "timeint_type)/='RK4'",'NCAR_topo.or.sso','size(block_catalog)/=640',
                   'packed(k-3:k)/=0','selected(i)/=selected_owner322(i,native_mode)',
                   'mode==2.or.mode==3.or.mode==4']:
        assert marker in text,marker
    lib=out/'codec.so'
    subprocess.run([cc,'-std=c99','-Wall','-Wextra','-Werror','-O2','-shared','-fPIC',str(source/'src/audit_compress322.c'),'-lz','-o',str(lib)],check=True)
    codec=ctypes.CDLL(str(lib));codec.audit_sha325.argtypes=[ctypes.c_void_p,ctypes.c_int64,ctypes.c_void_p]
    vectors=[b'',b'abc',b'a'*1000000]+[bytes((i*13+7)%256 for i in range(n)) for n in [1,7,8,31,32,55,56,57,63,64,65,119,120,127,128,129,1024]]
    for raw in vectors:
        data=ctypes.create_string_buffer(raw);digest=ctypes.create_string_buffer(32)
        assert codec.audit_sha325(data,len(raw),digest)==0
        assert digest.raw==hashlib.sha256(raw).digest()
    codec.audit_compress322_open.argtypes=[ctypes.c_char_p]
    codec.audit_digest325_write.argtypes=[ctypes.c_void_p,ctypes.c_int64]
    assert codec.audit_compress322_open(str(out/'digest.bin.gz').encode())==0
    for x in [float('nan'),float('inf'),-float('inf')]:
        data=ctypes.create_string_buffer(struct.pack('<d',x));assert codec.audit_digest325_write(data,8)!=0
    assert codec.audit_digest325_write(None,0)==0
    raw=struct.pack('<4d',0.,-0.,1.25,-1e100);data=ctypes.create_string_buffer(raw)
    assert codec.audit_digest325_write(data,len(raw))==0
    assert codec.audit_compress322_close()==0
    assert gzip.decompress((out/'digest.bin.gz').read_bytes())==hashlib.sha256(b'').digest()+hashlib.sha256(raw).digest()
    # Exercise actual Fortran noncontiguous/empty wrappers and footer lengths.
    (out/'emit.f90').write_text('''program emit
use audit316_mod
use iso_fortran_env, only: real64,int64
implicit none
real(real64) :: x(8)
x=[1._real64,2._real64,3._real64,4._real64,5._real64,6._real64,7._real64,8._real64]
call audit316_init(0)
call audit316_begin('test',1_int64)
call audit316_i('shape',[4])
call audit316_r('strided',x(1:8:2))
call audit316_r('empty',x(:0))
call audit316_end()
call audit316_finish()
end program
''')
    subprocess.run([cc,'-c','-O2',str(source/'src/audit_compress322.c'),'-o',str(out/'codec.o')],check=True)
    subprocess.run([fc,'-fcheck=all','-ffree-line-length-132',str(source/'src/audit316.f90'),'emit.f90','codec.o','-lz','-o','emit'],cwd=out,check=True)
    env={**os.environ,'WAVETRISK_EXECUTION_AUDIT316':'1','WAVETRISK_COMPACT_AUDIT325':'1','WAVETRISK_WORK_LEDGER325':'0'}
    subprocess.run([str(out/'emit')],cwd=out,env=env,check=True)
    raw=gzip.decompress((out/'execution316-rank-0.bin.gz').read_bytes())
    assert struct.unpack('<8s5i',raw[:28])==(b'RK325V01',1,0,4,8,0x01020304)
    offset=28;index=0;entries=[]
    while offset<len(raw):
        begin=offset;k,kind,count=struct.unpack_from('<iiq',raw,offset);offset+=16
        name=raw[offset:offset+k].decode();offset+=k
        size=32 if kind==13 else count*{1:4,2:8,3:8,4:1,10:8,11:8,12:8}[kind]
        data=raw[offset:offset+size];offset+=size
        if kind==12:assert struct.unpack('<2q',data)==(index,begin)
        entries.append((name,kind,count,data));index+=1
    assert entries[2][0:3]==('strided',13,4)
    assert entries[2][3]==hashlib.sha256(struct.pack('<4d',1.,3.,5.,7.)).digest()
    assert entries[3]==('empty',13,0,hashlib.sha256(b'').digest())
    assert offset==len(raw) and entries[-1][0]=='complete'
    report=dict(passed=True,selectors=5,invalid_selector_rejected=True,sha256_vectors=len(vectors),
                nonfinite_cases=3,fortran_codec_noncontiguous_and_empty=True,physical_footer_verified=True,
                current_home_blocks=sum(i//4==r for i,r in enumerate(a)),candidate_home_blocks=400,
                compiled_source_sha256={name:hashlib.sha256((source/'src'/name).read_bytes()).hexdigest()
                  for name in ['parallel_block_mpi.f90','audit316.f90','audit_compress322.c']},
                native_qualification=False)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    print(json.dumps(run_gate(*sys.argv[1:]),indent=2))
