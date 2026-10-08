import copy,json,shutil,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import recovery
from common import save,sha,VALIDATIONS

class Recovery(unittest.TestCase):
    def fixture(self,root):
        kit=root/'kit';prior=root/'prior';kit.mkdir();prior.mkdir()
        save(kit/'files-sha256.json',{})
        save(kit/'RECOVERY_SOURCE.json',dict(v1_manifest_sha256='v1'))
        manifests={}
        for kind in ('main','stage266','candidate'):
            tree=prior/'sources'/kind;tree.mkdir(parents=True)
            (tree/'source.f90').write_text(kind)
            manifests[kind]={'source.f90':sha(tree/'source.f90')}
        save(kit/'source-files.json',manifests)
        save(kit/'RECOVERY_SOURCE.json',dict(v1_manifest_sha256='v1',source_manifest_sha256=sha(kit/'source-files.json')))
        inputs=prior/'input';inputs.mkdir();(inputs/'simple.in').write_text('fixture')
        hashes={'simple.in':sha(inputs/'simple.in')}
        notes=dict(status='passed',level='J5',source_manifest_sha256=sha(kit/'source-files.json'),
                   physics_objects={str(i):str(i) for i in range(23)},builds={})
        for kind in ('main','stage266','candidate','checked'):
            (prior/kind).write_text('binary '+kind);(prior/('build-'+kind+'.log')).write_text('build '+kind)
            notes['builds'][kind]=dict(binary=kind,binary_sha256=sha(prior/kind),physics_objects=notes['physics_objects'],log_sha256=sha(prior/('build-'+kind+'.log')))
        save(prior/'builds.json',notes)
        result=dict(status='failed',error=recovery.EXPECTED_ERROR,kit_manifest_sha256='v1',builds=notes,schedule=[],runs={},job_id='111')
        for index,(n,kind,oracle) in enumerate(VALIDATIONS[:3],1):
            label=f'J5-n{n}-validation-{index}-{kind}';d=prior/label;d.mkdir()
            shutil.copy2(prior/kind,d/'climate');shutil.copy2(inputs/'simple.in',d/'simple.in');(d/'run.log').write_text('completed')
            result['schedule'].append(label);result['runs'][label]=dict(returncode=0,ranks=n,kind=kind,oracle=oracle)
        save(prior/'results.json',result)
        return kit,prior,hashes
    def test_reuse_preserves_original_build_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);kit,prior,hashes=self.fixture(root)
            with patch.object(recovery,'J5_HASHES',hashes):record=recovery.load_recovery(kit,prior)
            before={p:sha(p) for p in prior.rglob('*') if p.is_file()}
            out=root/'new';out.mkdir();notes=recovery.reuse_builds(kit,out,record)
            self.assertEqual(set(notes['builds']),{'main','stage266','candidate','checked'})
            for row in notes['builds'].values():self.assertEqual(sha(out/row['binary']),row['binary_sha256'])
            self.assertEqual({p:sha(p) for p in before},before)
            self.assertEqual(sha(out/'recovery-source-results.json'),record['results_sha256'])
    def test_changed_source_binary_and_log_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            kit,prior,hashes=self.fixture(Path(tmp))
            for file in (prior/'sources/candidate/source.f90',prior/'candidate',prior/'build-checked.log',prior/'input/simple.in'):
                old=file.read_bytes();file.write_bytes(b'changed')
                with patch.object(recovery,'J5_HASHES',hashes),self.assertRaises(ValueError):recovery.load_recovery(kit,prior)
                file.write_bytes(old)
    def test_wrong_failure_campaign_or_incomplete_solver_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            kit,prior,hashes=self.fixture(Path(tmp));p=prior/'results.json';original=json.loads(p.read_text())
            for field,value in (('status','complete'),('error','solver failed'),('kit_manifest_sha256','other'),('schedule',[])):
                result=copy.deepcopy(original);result[field]=value;save(p,result)
                with patch.object(recovery,'J5_HASHES',hashes),self.assertRaises(ValueError):recovery.load_recovery(kit,prior)
            result=copy.deepcopy(original);result['runs'][result['schedule'][-1]]['returncode']=143;save(p,result)
            with patch.object(recovery,'J5_HASHES',hashes),self.assertRaises(ValueError):recovery.load_recovery(kit,prior)
    def test_submit_prevents_duplicate_recovery(self):
        import control,subprocess
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);kit,prior,hashes=self.fixture(root)
            save(root/'rk274-latest.json',dict(job=111,phase='initial',out=str(prior),fixture=str(prior/'input')))
            with patch.object(recovery,'J5_HASHES',hashes),patch.object(Path,'home',return_value=root), \
                 patch.object(control,'preflight'),patch.object(control,'queued',return_value=False), \
                 patch.object(subprocess,'check_output',return_value='222\n') as submit:
                recovery.submit(kit)
                with self.assertRaises(ValueError):recovery.submit(kit)
                self.assertEqual(submit.call_count,1)
                state=json.loads((root/'rk274-latest.json').read_text())
                self.assertEqual(state['previous_job'],111)
                self.assertEqual(submit.call_args.args[0][-1],str(prior))

if __name__=='__main__':unittest.main()
