"""Review the numerical extraction evidence without claiming owner migration."""
import argparse,json,re,sys,hashlib
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import save,digest

def read_case(path):
    data=json.loads((path/'results.json').read_text())
    assert data['status']=='passed'
    for name,file in [('binary_sha256','climate'),('log_sha256','run.log'),('input_sha256','simple.in')]:
        if name in data:assert digest(path/file)==data[name],(path,name)
    return data

def review(root):
    root=root.resolve()
    long=read_case(root/'checked-oracles-4')
    short=read_case(root/'all-kernels-one-step')
    negative=read_case(root/'corrupt-pressure')
    assert long['steps']==12 and long['coverage']['remaps']==2
    assert long['numerical']['exact_observed'] and long['step_states_exact'] and long['vtk_exact_members']==62
    raw=(root/'checked-oracles-4/run.log').read_text()
    counts=[[int(v) for v in s.split()] for s in re.findall(r'Primitive232 exact counts:([^\n]*)',raw)]
    assert len(counts)==4 and all(len(c)==7 for c in counts),counts
    totals=[sum(c[i] for c in counts) for i in range(7)]
    assert all(c>0 for c in totals)
    assert all(c>0 for c in short['counts'].values())
    assert negative['expected_failure'] and negative['returncode']!=0
    kernel=digest(BASE/'src/parallel_block_primitive_kernels.f90')
    builds={}
    for name,case in [('candidate',long),('complete-kernel-check',short)]:
        tree=root/name
        assert digest(tree/'src/parallel_block_primitive_kernels.f90')==kernel
        assert digest(tree/'bin-230-checked/climate')==case['binary_sha256']
        names=['Makefile','src/parallel_block_primitive_kernels.f90','src/primitive_oracle_232.f90',
               'src/ops.f90','src/diagnostics.f90','test/climate/climate.f90','test/climate/test_case_module.f90',
               'src/parallel_block.f90','src/parallel_block_mpi.f90','src/time_integr.f90']
        builds[name]={p:hashlib.sha256((tree/p).read_bytes()).hexdigest() for p in names}
    assert negative['original_source_sha256']==kernel
    result=dict(status='arithmetic-extraction-validated',step1_complete=False,production_owner_migration=False,
        kernel_sha256=kernel,test_source_hashes=builds,long_trajectory=long,long_kernel_counts=totals,
        all_kernel_short_check=short,corruption_control=negative,
        limitations=['Geometry-owner scheduling and Domain input adapters remain in production.',
                    'No all-inputs-unavailable production run or block-owner tendency path has been delivered.',
                    'These are local correctness checks, not timing, scaling or cluster qualification.'])
    save(BASE/'docs/parallel-blocks-stage232-local-results.json',result)
    print('PASS numerical extraction checks; step 1 is incomplete')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=BASE/'.stage232-validation')
    review(p.parse_args().root)
