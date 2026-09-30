"""Audit frozen sources, scientific outputs, storage absence and rejected invalid epochs."""
import argparse,json,re,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save

def review(tree,root,checked_case):
    tree,root,checked_case=(p.resolve() for p in (tree,root,checked_case))
    build=json.loads((tree/'build-247.json').read_text())
    assert set(build['binaries'])=={'storage','production','oracle','control'}
    for key in ['changed_sources','package_sources']:
        for name,sha in build[key].items():assert digest(BASE/name)==digest(tree/name)==sha,name
    for name,sha in build['resolved_sources'].items():assert digest((tree/name).resolve())==sha,name
    for name,sha in build['binaries'].items():assert digest(tree/('bin-'+name)/'climate')==sha,name
    cases={}
    for out,variant,ranks,steps,views in [
        (checked_case,'storage',4,12,False),(root/'production4-4','production',4,12,False),
        (root/'production4-41','production',41,12,False),(root/'oracle4-4','oracle',4,1,True),
        (root/'control4-4','control',4,1,True)]:
        p=out/'storage-results.json';r=json.loads(p.read_text())
        assert r['status']=='passed' and r['ranks']==ranks and r['steps']==steps and r['legacy_views']==views
        assert digest(out/'climate')==r['binary_sha256']==build['binaries'][variant]
        assert digest(out/'run.log')==r['log_sha256']
        assert digest(out/'lifecycle-results.json')==r['lifecycle_results_sha256']
        if not views:assert r['absence_checks']>=ranks*steps and r['bounded_output_checks']>=ranks
        if steps==12:
            f=json.loads((out/'results.json').read_text())
            assert f['step_states_exact'] and f['coverage']['remaps']==2
            assert f['vtk_exact_members']==62 and f['numerical']['passed']
        if variant=='oracle':assert (out/'run.log').read_text().count('Retained RK exact aliases rank=')==16
        cases[out.name]=dict(path=str(p.relative_to(BASE)),sha256=digest(p),steps=steps,ranks=ranks)
        if steps==12:cases[out.name]['trajectory_sha256']=digest(out/'results.json')
        if ranks==41:
            layouts=re.findall(r'Block columns plan rank=(\d+) geometry/owner/aliases/domains=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',(out/'run.log').read_text())
            empty=sorted({int(rank) for rank,geometry,owner,aliases,domains in layouts if int(geometry)==0 and int(owner)>0 and int(domains)==0})
            assert empty,'41-rank case must exercise compute owners without Domain geometry'
            cases[out.name]['empty_geometry_compute_ranks']=empty
    for mode in ['lifecycle','coefficients','mirror']:
        out=root/('reject-'+mode);p=out/'results.json';r=json.loads(p.read_text())
        source='src/domain.f90' if mode=='mirror' else 'src/parallel_block_inverse.f90'
        assert r['status']=='passed' and r['expected_failure'] and r['returncode']!=0
        assert r['original_source_sha256']==digest(tree/source)
        assert r['injected_source_sha256']==digest(out/'injected.f90')
        assert r['binary_sha256']==digest(out/'climate') and r['log_sha256']==digest(out/'run.log')
        assert r['diagnostic'] in (out/'run.log').read_text()
        cases['reject-'+mode]=dict(path=str(p.relative_to(BASE)),sha256=digest(p),expected_failure=True)
    for stage,path in [(244,'.stage244-validation/candidate3'),(245,'.stage245-validation/candidate2'),(246,'.stage246-validation/candidate')]:
        old=BASE/path;record=json.loads((old/f'build-{stage}.json').read_text())
        for n,sha in record['changed_sources'].items():assert digest(old/n)==sha,n
        for n,sha in record['binaries'].items():assert digest(old/('bin-'+n)/'climate')==sha,n
    audit=BASE/'docs/parallel-blocks-stage247-kernels.json'
    kernels=json.loads(audit.read_text());assert kernels['status']=='passed' and len(kernels['kernels'])==15
    for n,sha in kernels['sources'].items():assert digest(BASE/n)==sha,n
    main=subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip()
    assert main=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
    helpers={str(p.relative_to(BASE)):digest(p) for p in sorted(Path(__file__).parent.iterdir()) if p.is_file()}
    save(BASE/'docs/parallel-blocks-stage247-results.json',dict(stage=247,status='passed',build=build,cases=cases,
         helpers=helpers,kernel_audit_sha256=digest(audit),main_unchanged=main,frozen_candidate=str(tree.relative_to(BASE)),
         scope='Restarted compressible climate RK4 with Simple physics. Legacy SOL/WAV numerical views exist only inside bootstrap/output/checkpoint adapters; obsolete RK snapshots and persistent scalar writeback images are absent.',
         retained_workspace='Native vector restriction/wavelet scratch, native stencil/alias workspaces, fixed reference fields and communication buffers remain active production storage.',
         performance_claim=None,bbserv_status='not run'))
    print('PASS integrated source-bound storage review: eight cases',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','root','checked-case'):p.add_argument('--'+n,type=Path,required=True)
    review(**vars(p.parse_args()))
