"""Index Stage251 owner operations by depth without changing their order."""
import hashlib
from pathlib import Path
import sys

BASELINE_SHA = 'eff68d917b06eebea30aa81f12ec7bb1c8163de3c502e27b087025da4a66a8e8'
OPTIMIZED_SHA = 'a65f465f483e789f68f0f454bb91088fe81309d4cf51d62950a0d5bdd299f2dc'


def replace_once(source, before, after):
    if source.count(before) != 1:
        raise ValueError('Missing or ambiguous Stage251 anchor: '+before[:65])
    return source.replace(before, after, 1)


def transform(program):
    program = Path(program)
    original = program.read_bytes()
    if hashlib.sha256(original).hexdigest() != BASELINE_SHA:
        raise ValueError('Require the exact Stage251 c5 owner-program source')
    code = original.decode()
    code = replace_once(code,
        '     integer, allocatable :: table(:)\n     type(Owner_Read), allocatable :: read(:)\n',
        '     integer, allocatable :: table(:)\n'
        '     integer, allocatable :: depth_start(:),depth_order(:)\n'
        '     type(Owner_Read), allocatable :: read(:)\n')
    code = replace_once(code,
        '    real(dp),allocatable::sg(:),rg(:)\n'
        '    integer::i,r,a,b,n,h,t,depth,ierr,ni,ng,ns,cursor\n',
        '    real(dp),allocatable::sg(:),rg(:)\n'
        '    integer,allocatable::depth_cursor(:)\n'
        '    integer::i,r,a,b,n,h,t,depth,ierr,ni,ng,ns,cursor\n')
    code = replace_once(code,
        "    if(ierr/=MPI_SUCCESS) call die('program depth')\n"
        '    allocate(p%read(0:p%maxdepth))\n',
        "    if(ierr/=MPI_SUCCESS) call die('program depth')\n"
        '    ! Stable buckets retain the original owned-operation order within\n'
        '    ! each dependency depth. The execution loop then touches only\n'
        '    ! operations that are ready at that depth.\n'
        '    allocate(p%depth_start(0:p%maxdepth+1),p%depth_order(n), &\n'
        '         depth_cursor(0:p%maxdepth))\n'
        '    p%depth_start=0\n'
        '    do i=1,n\n'
        '       depth=p%owned(i)%depth\n'
        '       p%depth_start(depth+1)=p%depth_start(depth+1)+1\n'
        '    end do\n'
        '    p%depth_start(0)=1\n'
        '    do depth=1,p%maxdepth+1\n'
        '       p%depth_start(depth)=p%depth_start(depth-1)+p%depth_start(depth)\n'
        '    end do\n'
        '    depth_cursor=p%depth_start(0:p%maxdepth)\n'
        '    do i=1,n\n'
        '       depth=p%owned(i)%depth\n'
        '       p%depth_order(depth_cursor(depth))=i\n'
        '       depth_cursor(depth)=depth_cursor(depth)+1\n'
        '    end do\n'
        '    allocate(p%read(0:p%maxdepth))\n')
    code = replace_once(code,
        '    integer::depth,i,a,n\n'
        "    if(.not.p%compiled) call die('execute before compile')\n",
        '    integer::depth,i,a,n,slot\n'
        "    if(.not.p%compiled) call die('execute before compile')\n")
    code = replace_once(code,
        '       do i=1,size(p%owned)\n'
        '          associate(o=>p%owned(i))\n'
        '          if(o%depth/=depth) cycle\n'
        '          a=o%input_start;n=size(o%input,2)\n',
        '       do slot=p%depth_start(depth),p%depth_start(depth+1)-1\n'
        '          i=p%depth_order(slot)\n'
        '          associate(o=>p%owned(i))\n'
        '          a=o%input_start;n=size(o%input,2)\n')
    digest = hashlib.sha256(code.encode()).hexdigest()
    if digest != OPTIMIZED_SHA:
        raise ValueError('Stage255 transformed source hash differs')
    program.write_text(code)
    return digest


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: transform.py STAGE251_SOURCE/src/parallel_block_program.f90')
    print(transform(sys.argv[1]))
