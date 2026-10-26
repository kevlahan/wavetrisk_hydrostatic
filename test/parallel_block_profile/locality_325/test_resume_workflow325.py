"""Composed offline checks; constructed data never qualify the native campaign."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import signal
import sqlite3
import subprocess
import sys
import tempfile
import unittest

import inverse_identity325 as inv
from resume_state325 import ResumeStateError325
import resume_workflow325 as flow
import verify325 as q
from test_provider_paths325 import inverses
from test_stream_compare325 import fixture, change


class FixtureBackend325:
    native = False
    expected_ranks = 256
    steps = 1
    expected_families = 1
    expected_family_keys = [('static', 'inverse.alias/fixture')]
    expected_table_rows = dict(route=256, domain=256, context=512, feature=8)
    arm_names = ('qual-C0', 'qual-L256')

    def __init__(self, inputs):
        self.inputs = Path(inputs)
        self.protected_roots = (self.inputs,)
        self.calls = dict(cache=0, work=0, preservation=0)
        self.interrupt_work = False
        self.interrupt_binding = False
        self.false_cache = False
        self.incomplete_reads = False
        self._binding = None
        if not self.inputs.exists():
            self.inputs.mkdir()
            for arm in self.arm_names:
                path = fixture(self.inputs/(arm+'.sqlite'), arm, ranks=256, steps=1)
                plans, contexts, _ = inverses('inverse.alias', unallocated=True)
                c = sqlite3.connect(path)
                for rank in range(256):
                    c.execute('insert into route values(?,?,?,?,?)',
                        (arm,rank,'static','inverse.alias/fixture',json.dumps(plans[rank])))
                    c.execute('insert into context values(?,?,?,?,?)',
                        (arm,rank,'static','inverse.alias/fixture',json.dumps(contexts[rank])))
                c.commit(); c.close()
            for key in ('source','checker','config','captures'):
                (self.inputs/key).write_text('original-'+key)

    def binding(self):
        if self.interrupt_binding:
            os.kill(os.getpid(), signal.SIGTERM)
        values = {key: flow.sha256325(self.inputs/key) for key in ('source','checker','config','captures')}
        values['captures'] = flow._sha_value(dict(original=values['captures'], databases={
            arm:flow.sha256325(self.inputs/(arm+'.sqlite')) for arm in self.arm_names}))
        self._binding = values
        return dict(values)

    def bind_sources(self):
        inv.bind_inverse_bridge_sources325({'candidate':dict(inv.EXPECTED_SOURCE_HASHES325)})

    def cache(self, stage, allowance):
        self.calls['cache'] += 1
        for arm in self.arm_names:
            (stage/(arm+'.sqlite')).write_bytes((self.inputs/(arm+'.sqlite')).read_bytes())
        return dict(artifacts=dict(left_db='qual-C0.sqlite',right_db='qual-L256.sqlite'),
                    census=dict(passed=not self.false_cache,streams=512,
                                capture_reads=511 if self.incomplete_reads else 512))

    def work(self, stage, allowance):
        self.calls['work'] += 1
        if self.interrupt_work:
            (stage/'incomplete').write_bytes(b'partial work')
            os.kill(os.getpid(), signal.SIGTERM)
        flow.save325(stage/'work.json',dict(passed=True,fixture=True))
        return dict(artifacts={'proof':'work.json'},census=dict(passed=True,fixture=True))

    def preservation(self, stage, allowance):
        self.calls['preservation'] += 1
        original = dict(self._binding)
        if self.binding() != original:
            raise flow.SemanticFailure325('fixture inputs changed')
        flow.save325(stage/'preservation.json',dict(passed=True,fixture=True))
        return dict(artifacts={'proof':'preservation.json'},census=dict(passed=True,fixture=True))


class ComposedWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name).resolve()
        self.backend = FixtureBackend325(self.base/'inputs')
        self.root = self.base/'state'
        self.limits = {name:2*flow.MiB for name in
            ('providers-C0','providers-L256','pair','work','preservation')}
        self.limits['cache'] = 80*flow.MiB

    def run_workflow(self, **kwargs):
        return flow.run325(self.backend,self.root,max_seconds=30,cap_bytes=128*flow.MiB,
                          phase_limits=self.limits,**kwargs)

    def test_real_provider_feature_bridge_and_certification_complete(self):
        before = {path:flow.sha256325(path) for path in self.backend.inputs.iterdir()}
        report = self.run_workflow()
        self.assertEqual(report['status'],'passed')
        self.assertFalse(report['native_qualification'])
        self.assertEqual(len(report['certificate_sha256']),6)
        self.assertEqual(report['progress_only']['completed'],1)
        self.assertEqual({path:flow.sha256325(path) for path in self.backend.inputs.iterdir()},before)
        self.assertEqual(self.backend.calls,dict(cache=1,work=1,preservation=1))
        pair = json.loads((self.root/'phases/pair/pair.json').read_text())
        self.assertEqual(pair['physical']['inverse_bridge']['proved'],256)
        self.assertEqual(pair['physical']['features']['equal'],8)
        self.assertLess((self.root.with_name('state.attempt.json')).stat().st_size,flow.REPORT_LIMIT325)
        self.assertTrue(all('commit_and_full_rehash_seconds' in item for item in report['phases']))
        used = sum(p.stat().st_size for p in self.base.rglob('*') if p.is_file() and not p.is_relative_to(self.backend.inputs))
        self.assertLess(used,128*flow.MiB)

    def test_real_sigterm_retains_four_completed_phases_and_resume_skips_them(self):
        self.backend.interrupt_work = True
        with self.assertRaises(flow.ReplayInterrupted325): self.run_workflow()
        report = json.loads((self.base/'state.attempt.json').read_text())
        self.assertEqual(report['status'],'interrupted')
        self.assertEqual(report['signal']['number'],signal.SIGTERM)
        self.assertEqual(report['signal']['name'],'SIGTERM')
        self.assertEqual(report['signal']['phase'],'work')
        self.assertEqual(len(list((self.root/'phases').iterdir())),4)
        proof = flow.sha256325(self.root/'phases/providers-C0/providers.json')
        self.backend.interrupt_work = False
        with self.assertRaisesRegex(ResumeStateError325,'explicit_resume_required'):
            with contextlib.redirect_stderr(io.StringIO()): self.run_workflow()
        with self.assertRaisesRegex(ResumeStateError325,'explicit_partial_discard_required'):
            self.run_workflow(explicit_resume=True)
        report = self.run_workflow(explicit_resume=True,discard_partial=True)
        self.assertEqual(self.backend.calls,dict(cache=1,work=2,preservation=1))
        self.assertEqual([x['reused'] for x in report['phases']],[True]*4+[False]*2)
        self.assertEqual(flow.sha256325(self.root/'phases/providers-C0/providers.json'),proof)
        self.assertEqual(report['attempt_history'][0]['signal']['number'],signal.SIGTERM)

    def test_equal_missing_feature_is_not_mistaken_for_equivalence(self):
        for arm in self.backend.arm_names:
            change(self.backend.inputs/(arm+'.sqlite'),"delete from feature where key='key-00000000'")
        with self.assertRaisesRegex(flow.SemanticFailure325,'exact table census'): self.run_workflow()
        self.assertTrue((self.root/'semantic-failure.json').is_file())
        self.assertFalse((self.root/'phases/cache').exists())
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaisesRegex(ResumeStateError325,'semantic_failed_latch'):
            self.run_workflow(explicit_resume=True,discard_partial=True)

    def test_false_backend_closure_and_missing_read_are_rejected(self):
        for name in ('false_cache','incomplete_reads'):
            with self.subTest(name=name):
                self.root = self.base/name
                setattr(self.backend,name,True)
                with self.assertRaisesRegex(flow.SemanticFailure325,'actual-read closure'): self.run_workflow()
                self.assertTrue((self.root/'semantic-failure.json').exists())
                setattr(self.backend,name,False)

    def test_physical_difference_latches_failure_before_resume(self):
        change(self.backend.inputs/'qual-L256.sqlite','update feature set digest=?',(hashlib.sha256(b'changed').hexdigest(),))
        with self.assertRaisesRegex(flow.SemanticFailure325,'cross-arm difference'): self.run_workflow()
        self.assertTrue((self.root/'semantic-failure.json').exists())
        self.assertFalse((self.root/'phases/pair').exists())
        self.assertTrue((self.root/'partial/pair/pair.json').exists())
        self.assertEqual(self.backend.calls['work'],0)

    def test_same_size_preserved_mtime_input_changes_reject_resume(self):
        self.run_workflow()
        for key in ('source','checker','config','captures'):
            path = self.backend.inputs/key
            original = path.read_bytes(); metadata=path.stat()
            path.write_bytes(b'X'+original[1:])
            os.utime(path,ns=(metadata.st_atime_ns,metadata.st_mtime_ns))
            with self.subTest(key=key), contextlib.redirect_stderr(io.StringIO()), self.assertRaisesRegex(ResumeStateError325,'binding_changed'):
                self.run_workflow(explicit_resume=True)
            path.write_bytes(original)

    def test_preownership_signal_has_bounded_exact_stderr_record(self):
        self.backend.interrupt_binding = True
        output = io.StringIO()
        with contextlib.redirect_stderr(output), self.assertRaises(flow.ReplayInterrupted325): self.run_workflow()
        event=json.loads(output.getvalue())['preownership_termination']
        self.assertEqual(event['signal']['number'],signal.SIGTERM)
        self.assertEqual(event['signal']['phase'],'fresh_input_binding')
        self.assertFalse((self.base/'state.attempt.json').exists())

    def test_full_resume_reuses_every_phase_but_freshly_hashes_inputs(self):
        self.run_workflow()
        report=self.run_workflow(explicit_resume=True)
        self.assertTrue(all(item['reused'] for item in report['phases']))
        self.assertEqual(self.backend.calls,dict(cache=1,work=1,preservation=1))

    def test_real_sigkill_has_no_fabricated_signal_and_certified_phases_survive(self):
        code = '''
import os,signal,sys
from pathlib import Path
from test_resume_workflow325 import FixtureBackend325
from resume_workflow325 import run325,MiB
base=Path(sys.argv[1]); backend=FixtureBackend325(base/'inputs')
def killed_work(stage,allowance):
    (stage/'partial').write_bytes(b'uncommitted')
    os.kill(os.getpid(),signal.SIGKILL)
backend.work=killed_work
limits={name:2*MiB for name in ('providers-C0','providers-L256','pair','work','preservation')}
limits['cache']=80*MiB
run325(backend,base/'state',max_seconds=30,cap_bytes=128*MiB,phase_limits=limits)
'''
        result=subprocess.run([sys.executable,'-c',code,str(self.base)],capture_output=True,timeout=30)
        self.assertEqual(result.returncode,-signal.SIGKILL,result.stderr.decode())
        prior=json.loads((self.base/'state.attempt.json').read_text())
        self.assertIsNone(prior['signal'])  # SIGKILL cannot run a logging handler.
        self.assertEqual(prior['status'],'running')
        report=self.run_workflow(explicit_resume=True,discard_partial=True)
        self.assertEqual(self.backend.calls,dict(cache=0,work=1,preservation=1))
        self.assertEqual([x['reused'] for x in report['phases']],[True]*4+[False]*2)
        self.assertEqual(report['attempt_history'][0]['status'],'running')

    def test_runtime_alarm_is_identified_separately_and_handler_restored(self):
        prior=signal.getsignal(signal.SIGALRM)
        self.backend.binding=lambda: signal.pause()
        error=io.StringIO()
        with contextlib.redirect_stderr(error),self.assertRaises(flow.ReplayInterrupted325):
            flow.run325(self.backend,self.root,max_seconds=0.02,cap_bytes=128*flow.MiB,
                        phase_limits=self.limits)
        self.assertEqual(signal.getsignal(signal.SIGALRM),prior)
        event=json.loads(error.getvalue())['preownership_termination']
        self.assertEqual(event['signal']['number'],signal.SIGALRM)
        self.assertEqual(event['signal']['name'],'SIGALRM')

    def test_real_accounting_valueerror_creates_semantic_latch(self):
        import account325
        from resume_native325 import semantic_phase325
        @semantic_phase325
        def bad_work(stage,allowance):
            path=stage/'bad.tsv'
            path.write_bytes(b'wrong work header\n')
            account325.parse(path,0,{})
        self.backend.work=bad_work
        with self.assertRaisesRegex(flow.SemanticFailure325,'account325.header'):
            self.run_workflow()
        latch=json.loads((self.root/'semantic-failure.json').read_text())
        self.assertEqual(latch['phase'],'work')
        self.assertEqual(latch['status'],'semantic_failed')
        self.assertEqual(self.backend.calls['preservation'],0)

    def test_fixture_native_flag_cannot_grant_native_gate(self):
        self.backend.native=True
        with self.assertRaisesRegex(ResumeStateError325,'native_backend_required'):
            self.run_workflow()
        self.assertFalse(self.root.exists())


if __name__=='__main__': unittest.main()
