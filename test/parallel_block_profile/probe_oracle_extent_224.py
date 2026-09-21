"""Reproduce the oracle extent predicate using the exact cluster source.

Serial checked Fortran predicate test only; no solver changes or MPI claims.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def probe(source,output):
    text=source.read_text()
    start=text.index('    if (production_vector_transform_displ(1) /= 0 .or. &')
    end=text.index('    end if',start)+len('    end if')
    original=text[start:end]
    assert 'complete vector transform stage extent mismatch' in original
    anchor='production_vector_transform_displ(size(grid)+1) /= &'
    replacement='max(1,production_vector_transform_displ(size(grid)+1)) /= &'
    input_corrected=replacement in original
    if input_corrected:
        fixed=original
        original=original.replace(replacement,anchor)
    else:
        assert original.count(anchor)==1
        fixed=original.replace(anchor,replacement)
    rows=[]
    with tempfile.TemporaryDirectory(prefix='stage224-oracle-extent-') as tmp:
        root=Path(tmp)
        for name,condition in [('original',original),('allocation-consistent',fixed)]:
            body=condition.replace('call fail("complete vector transform stage extent mismatch")','rejected = .true.')
            program='''program probe
implicit none
integer,allocatable :: grid(:),production_vector_transform_displ(:)
real,allocatable :: production_vector_transform_stage(:),production_vector_wavelet_stage(:)
logical,allocatable :: production_vector_transform_covered(:)
logical :: rejected
integer :: domains,logical_values,capacity,covered_capacity
character(16) :: argument
call get_command_argument(1,argument)
read(argument,*) domains
call get_command_argument(2,argument)
read(argument,*) logical_values
call get_command_argument(3,argument)
read(argument,*) capacity
call get_command_argument(4,argument)
read(argument,*) covered_capacity
allocate(grid(domains),production_vector_transform_displ(domains+1))
production_vector_transform_displ=0
production_vector_transform_displ(domains+1)=logical_values
allocate(production_vector_transform_stage(capacity),production_vector_wavelet_stage(capacity))
allocate(production_vector_transform_covered(covered_capacity))
rejected=.false.
'''+body+'''
if(rejected) then
 print *, 'REJECT'
else
 print *, 'ACCEPT'
end if
end program
'''
            path=root/(name+'.f90');path.write_text(program)
            subprocess.run(['gfortran','-fcheck=all','-o',str(root/name),str(path)],cwd=root,check=True,capture_output=True,text=True)
            cases=[('empty-valid',0,0,1,1),('empty-missing-dummy',0,0,0,0),('nonempty-valid',1,8,8,8),
                   ('nonempty-wrong-capacity',1,8,9,9),('nonempty-wrong-covered',1,8,8,9)]
            for case,*args in cases:
                result=subprocess.check_output([str(root/name),*map(str,args)],text=True).strip()
                expected='ACCEPT' if case in (('empty-missing-dummy','nonempty-valid') if name=='original' else ('empty-valid','nonempty-valid')) else 'REJECT'
                assert result==expected,(name,case,result)
                rows.append(dict(variant=name,case=case,domains=args[0],logical_values=args[1],capacity=args[2],covered_capacity=args[3],result=result))
    output.write_text(json.dumps(dict(status='passed',scope='Exact source predicate reproduced; candidate correction tested only in generated Fortran',
                                     input_predicate_corrected=input_corrected,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cases=rows),indent=2)+'\n')
    print('PASS: original rejects valid empty-rank storage; allocation-consistent predicate accepts it and rejects tested malformed extents.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();probe(a.source,a.output)
