#!/usr/bin/env python3
"""Compile a small MPI test of the actual checkpoint output helpers."""
import argparse
from pathlib import Path
import re
import subprocess

DRIVER = '''program probe
  use checkpoint_output_test_mod
  implicit none
  type(MPI_File) :: fh
  type(MPI_Status) :: status
  integer :: ierr,i,r,fid,ios
  integer(int8) :: payload(8192),values(8192)
  character(30) :: mode
  call MPI_Init(ierr)
  comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr)
  call MPI_Comm_size(comm,n_process,ierr)
  call get_command_argument(1,mode)
  allocate(cp_load(N_GLO_DOMAIN),cp_offset(N_GLO_DOMAIN),cp_nbytes(N_GLO_DOMAIN))
  cp_load=1
  cp_offset=[(CP_DATA_POS+16_int64*(i-1),i=1,N_GLO_DOMAIN)]
  cp_nbytes=16_int64
  payload=int(rank+1,int8)
  call open_checkpoint_output(fh,'output.bin')
  if(rank==0)call write_checkpoint_directory(fh)
  call sync_checkpoint_output(fh,'metadata')
  call MPI_File_write_at_all(fh,int(CP_DATA_POS+8192*rank,MPI_OFFSET_KIND),payload,8192,MPI_BYTE,status,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'payload write failed'
  call require_mpi_transfer_count(status,MPI_BYTE,8192,'test payload')
  call sync_checkpoint_output(fh,'payload')
  call MPI_File_close(fh,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'close failed'
  call MPI_Barrier(comm,ierr)
  if(rank==0)then
    call verify_checkpoint_output_directory('output.bin')
    open(newunit=fid,file='output.bin',status='old',access='stream',form='unformatted')
    do r=0,n_process-1
      read(fid,pos=CP_DATA_POS+8192*r+1,iostat=ios)values
      if(ios/=0)error stop 'short payload'
      if(any(values/=int(r+1,int8)))error stop 'wrong payload'
    end do
    select case(trim(mode))
    case('bad-header')
      write(fid,pos=1)0_int64
    case('bad-directory')
      write(fid,pos=CP_LOAD_POS+1)999
    end select
    close(fid)
    call verify_checkpoint_output_directory('output.bin')
    write(*,'(a)')'PASS actual checkpoint writer metadata and payload'
  end if
  call MPI_Barrier(comm,ierr)
  call MPI_Finalize(ierr)
end program probe
'''


def prepare(source, out, compiler='mpif90', opt='-O2'):
    source = source.read_text()
    constants = source[source.index('  ! Magic number'):source.index('contains')]
    names = ['require_mpi_transfer_count', 'open_checkpoint_output', 'sync_checkpoint_output',
             'verify_checkpoint_output_directory', 'write_checkpoint_directory']
    helpers = []
    for name in names:
        match = re.search(r'^  subroutine ' + name + r'\b.*?^  end subroutine ' + name + r'\b', source, re.M | re.S)
        if not match:
            raise ValueError('Missing actual checkpoint helper: ' + name)
        helpers.append(match[0])
    module = '''module checkpoint_output_test_mod
use mpi_f08
use iso_fortran_env, only: int8,int64,error_unit
implicit none
integer,parameter :: N_GLO_DOMAIN=160
integer :: rank,n_process
integer,allocatable :: cp_load(:)
type(MPI_Comm) :: comm
''' + constants + '\ncontains\n' + '\n\n'.join(helpers) + '\nend module\n'
    out.mkdir(parents=True, exist_ok=False)
    (out / 'probe.f90').write_text(module + DRIVER)
    subprocess.run([compiler, opt, '-g', '-Wall', '-Wextra', '-Werror', '-fcheck=all',
                    '-ffree-line-length-132', 'probe.f90', '-o', 'probe'], cwd=out, check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[2] / 'src/checkpoint.f90')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--compiler', default='mpif90')
    parser.add_argument('--optimization', choices=('O0', 'O2'), default='O2')
    args = parser.parse_args()
    prepare(args.source, args.out, args.compiler, '-' + args.optimization)
