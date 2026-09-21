"""Bounded serial I/O reproducer for the Stage223 checkpoint prefix failure.

Extracts the actual distribution routine, stubs its MPI broadcasts and tests
one versus two concatenated Domain records. This is not an MPI/solver gate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

PRELUDE = '''module probe_mod
use iso_fortran_env, only: int8,int64
implicit none
integer :: rank=0,owner(1)=[0],n_process=2,comm=0
integer,parameter :: MPI_INTEGER8=1,MPI_BYTE=2,MPI_SUCCESS=0
integer(int64) :: cp_nbytes(1)=[16_int64]
interface MPI_Bcast
module procedure bcast_extent,bcast_bytes
end interface
contains
subroutine bcast_extent(data,n,kind,source,communicator,ierr)
integer(int64) :: data
integer :: n,kind,source,communicator,ierr
ierr=0
end subroutine
subroutine bcast_bytes(data,n,kind,source,communicator,ierr)
integer :: n,kind,source,communicator,ierr
integer(int8) :: data(n)
if(transfer(data,0_int64)/=123456789_int64) error stop 'Wrong global prefix bytes'
ierr=0
end subroutine
subroutine read_checkpoint_global_state(fid,position)
integer :: fid
integer(int64) :: position
error stop 'Source-only reproducer must not call receiver'
end subroutine
'''
DRIVER = '''
end module
program probe
use probe_mod
implicit none
integer :: fid,i,records
integer(int64) :: positions(1),prefix_bytes
character(8) :: argument
call get_command_argument(1,argument)
read(argument,*) records
open(newunit=fid,status='scratch',form='unformatted',access='stream',action='readwrite')
do i=1,records
write(fid) 123456789_int64,987654321_int64
end do
flush(fid)
! Reconstruction consumes the first record and advances its cursor to 17.
positions=[17_int64]
prefix_bytes=8_int64
call distribute_checkpoint_global_state(fid,positions,prefix_bytes)
close(fid)
print *, 'PASS: exact global metadata prefix'
end program
'''


def run(source, output):
    text=source.read_text()
    start=text.index('  subroutine distribute_checkpoint_global_state (')
    end=text.index('  end subroutine distribute_checkpoint_global_state',start)+len('  end subroutine distribute_checkpoint_global_state')
    original=text[start:end]
    # The first Domain record starts at position 1 in this scratch stream.
    # Test this seek correction in isolation; do not edit production sources.
    anchor='read(fid,pos=domain_pos(1)) prefix'
    assert original.count(anchor)==1
    fixed=original.replace(anchor,'read(fid,pos=1_int64) prefix')
    rows=[]
    with tempfile.TemporaryDirectory(prefix='stage223-prefix-probe-') as tmp:
        root=Path(tmp)
        for name,routine in [('original',original),('corrected-seek',fixed)]:
            src=root/(name+'.f90');src.write_text(PRELUDE+routine+DRIVER)
            subprocess.run(['gfortran','-g','-fcheck=all','-o',str(root/name),str(src)],cwd=root,check=True,capture_output=True,text=True)
            for records in (1,2):
                r=subprocess.run([str(root/name),str(records)],cwd=root,capture_output=True,text=True)
                expected=0 if name=='corrected-seek' or records==2 else 2
                assert r.returncode==expected,(name,records,r.stderr)
                if expected:assert 'End of file' in r.stderr
                rows.append(dict(variant=name,domain_records=records,returncode=r.returncode,stdout=r.stdout,stderr=r.stderr))
    report=dict(status='passed',scope='Serial source-routine I/O reproduction; MPI broadcasts stubbed, not a full solver validation',
                source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cases=rows)
    output.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: original single-record EOF; multi-record masking; corrected seek succeeds for both.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();run(a.source,a.output)
