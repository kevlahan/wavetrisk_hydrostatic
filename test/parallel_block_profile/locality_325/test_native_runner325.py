"""Allocation, resource and owned-process cleanup gates; no solver is invoked."""
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import native_runner325 as n


class NativeRunnerGuards(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name).resolve()

    def test_allocation_rejects_wrong_host_cpu_count_partition_and_limit(self):
        env=dict(SLURM_JOB_ID='1',SLURM_JOB_PARTITION='nk',SLURM_JOB_NUM_NODES='1',
                 SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='16')
        good=dict(env=env,hostname='bb13.rhpcs.internal',cpus=list(range(16)),limits=(n.CAP_MEMORY325,)*2)
        with patch.object(n.platform,'system',return_value='Linux'):
            self.assertEqual(n.allocation325(**good)['cpus'],list(range(16)))
            for key,value in [('hostname','bbserv'),('cpus',list(range(256))),('limits',(-1,-1)),
                              ('env',dict(env,SLURM_JOB_PARTITION='bb'))]:
                with self.subTest(key=key),self.assertRaises(ValueError):
                    n.allocation325(**dict(good,**{key:value}))

    def test_derived_logical_bytes_include_nested_and_temporary_files(self):
        (self.root/'state').mkdir();(self.root/'state/data').write_bytes(b'a'*19)
        (self.root/'state.attempt.tmp').write_bytes(b'b'*23)
        self.assertEqual(n.tree_bytes325(self.root),42)
        (self.root/'linked').symlink_to(self.root/'state/data')
        with self.assertRaisesRegex(ValueError,'derived_symlink'):n.tree_bytes325(self.root)

    def test_step_installs_exact_soft_and_hard_caps_even_if_inherited_hard_is_unlimited(self):
        before=(n.CAP_MEMORY325,n.resource.RLIM_INFINITY)
        with patch.object(n.resource,'getrlimit',side_effect=[before,(n.CAP_MEMORY325,)*2]),\
             patch.object(n.resource,'setrlimit') as setter:
            proof=n.install_address_limit325()
        setter.assert_called_once_with(n.resource.RLIMIT_AS,(n.CAP_MEMORY325,)*2)
        self.assertEqual(proof['before'],list(before));self.assertEqual(proof['after'],[n.CAP_MEMORY325]*2)
        with patch.object(n.resource,'getrlimit',return_value=(1024,1024)),\
             patch.object(n.resource,'setrlimit') as setter,\
             self.assertRaisesRegex(ValueError,'inherited_hard_limit_below_required_cap'):
            n.install_address_limit325()
        setter.assert_not_called()
        with patch.object(n.resource,'getrlimit',side_effect=[before,(n.CAP_MEMORY325,n.resource.RLIM_INFINITY)]),\
             patch.object(n.resource,'setrlimit'),\
             self.assertRaisesRegex(ValueError,'installed_address_limit'):
            n.install_address_limit325()

    @unittest.skipUnless(sys.platform.startswith('linux'),'native Linux address-space limits')
    def test_exact_caps_are_installed_and_inherited_in_real_subprocesses(self):
        script='''
import json,resource,subprocess,sys
from native_runner325 import install_address_limit325,CAP_MEMORY325
proof=install_address_limit325()
assert proof['after']==[CAP_MEMORY325]*2
child=subprocess.check_output([sys.executable,'-B','-c','import resource,json; print(json.dumps(resource.getrlimit(resource.RLIMIT_AS)))'],text=True)
assert json.loads(child)==proof['after']
print(json.dumps(proof))
'''
        result=subprocess.run([sys.executable,'-B','-c',script],capture_output=True,text=True,timeout=10)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(json.loads(result.stdout)['after'],[n.CAP_MEMORY325]*2)

    @unittest.skipUnless(sys.platform.startswith('linux'),'native Linux CLI bootstrap')
    def test_cli_logs_installed_limits_before_failing_allocation_without_reading_inputs(self):
        env=dict(os.environ);env.pop('SLURM_JOB_ID',None)
        result=subprocess.run([sys.executable,'-B',n.__file__,'supervise',str(self.root),str(self.root/'out')],
                              capture_output=True,text=True,env=env,timeout=10)
        self.assertNotEqual(result.returncode,0)
        lines=result.stderr.splitlines()
        self.assertIn('before',json.loads(lines[0])['address_limit_bootstrap'])
        proof=json.loads(lines[1])['address_limit_bootstrap_complete']
        self.assertEqual(proof['after'],[n.CAP_MEMORY325]*2)
        self.assertIn('job_identity',result.stderr)
        self.assertFalse((self.root/'out').exists())

    def test_process_sampler_includes_owned_group_after_reparenting(self):
        for pid,parent,group in [(10,1,10),(11,10,11),(12,1,11),(13,1,13)]:
            p=self.root/str(pid);p.mkdir()
            (p/'status').write_text(f'PPid:\t{parent}\nVmSize:\t10 kB\nVmRSS:\t3 kB\n')
            (p/'stat').write_text(f'{pid} (test name) S {parent} {group} 0 0')
        rows=n.processes325(10,group=11,proc=self.root)
        self.assertEqual(set(rows),{10,11,12})
        self.assertEqual(sum(row['virtual_bytes'] for row in rows.values()),30*1024)

    def test_unreapable_wait_remains_bounded_and_returns_unknown(self):
        class Child:
            pid=12345
            def poll(self):return None
            def wait(self,timeout):raise subprocess.TimeoutExpired('child',timeout)
        with patch.object(n,'kill_group325') as kill:
            self.assertIsNone(n.stop_group325(Child(),grace=0))
        self.assertEqual([x.args[1] for x in kill.call_args_list],[signal.SIGTERM,signal.SIGKILL])

    def test_package_inventory_and_byte_changes_fail_closed(self):
        kit=self.root/'kit';kit.mkdir();manifest={}
        for i in range(51):
            p=kit/f'f{i}';p.write_bytes(str(i).encode());manifest[p.name]=n.sha325(p)
        (kit/'native-files-sha256.json').write_text(json.dumps(manifest))
        self.assertEqual(len(n.package_check325(kit)),64)
        (kit/'f1').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'package_bytes'):n.package_check325(kit)
        (kit/'f1').write_bytes(b'1');(kit/'unknown').write_bytes(b'new')
        with self.assertRaisesRegex(ValueError,'package_inventory'):n.package_check325(kit)

    def test_preflight_predecessor_and_binding_paths_protected_before_writes(self):
        for protected in (self.root/'old/preflight',self.root/'replay/results',self.root/'kit/binding.json'):
            with self.subTest(protected=protected),self.assertRaisesRegex(ValueError,'protected_output_overlap'):
                n.protect_output325(protected,protected/'out',[protected])
            self.assertFalse(protected.exists())

    def test_completion_requires_exact_six_matching_real_certificates(self):
        from resume_state325 import PHASES325
        hashes={}
        for phase in PHASES325:
            p=self.root/'state/phases'/phase/'certificate.json';p.parent.mkdir(parents=True)
            p.write_text('{}');hashes[phase]=n.sha325(p)
        n.save325(self.root/'worker-complete.json',dict(passed=True,phase_certificates=hashes))
        n.save325(self.root/'state.attempt.json',dict(status='passed',native_gate_passed=True,certificate_sha256=hashes))
        self.assertTrue(n.completion_check325(self.root)['passed'])
        (self.root/'state/phases/cache/certificate.json').write_text('{"changed":true}')
        with self.assertRaisesRegex(ValueError,'phase_certificate_bytes'):n.completion_check325(self.root)
        n.save325(self.root/'worker-complete.json',dict(passed=True,phase_certificates={str(i):'a'*64 for i in range(6)}))
        with self.assertRaisesRegex(ValueError,'complete_gate_certificate'):n.completion_check325(self.root)


@unittest.skipUnless(sys.platform.startswith('linux'),'actual /proc cleanup runs in native preflight')
class NativeSupervisorProcesses(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name).resolve();self.kit=self.root/'kit';self.kit.mkdir()
        self.work=self.root/'work';self.work.mkdir(mode=0o700);self.out=self.work/'out'
        (self.kit/'recovery_v4p1').mkdir()
        (self.kit/'recovery_v4p1/binding325.json').write_text(json.dumps(dict(
            original_results=str(self.root/'protected'),original_kit=str(self.root/'oldkit'))))

    def fake_child(self,body):
        # This is a constructed process test. Package and allocation checks are
        # separately exercised above and never relaxed in the native entrypoint.
        script="import sys,os,signal,time,json,subprocess\nfrom pathlib import Path\nout=Path(sys.argv[3])\n"+body
        (self.kit/'native_runner325.py').write_text(script)

    def run_fake(self,**extra):
        with ExitStack() as stack:
            stack.enter_context(patch.object(n,'allocation325',return_value={'fixture':True}))
            stack.enter_context(patch.object(n,'install_address_limit325',side_effect=AssertionError(
                'library supervisor must not lower its caller hard limit')))
            stack.enter_context(patch.object(n,'package_check325',return_value='a'*64))
            stack.enter_context(patch.object(n,'protected_paths325',return_value=(self.kit,self.root/'protected')))
            stack.enter_context(patch.object(n,'SUPERVISOR_SECONDS325',.3))
            stack.enter_context(patch.object(n,'GRACE_SECONDS325',.2))
            for name,value in extra.items():stack.enter_context(patch.object(n,name,value))
            return n.supervise325(self.kit,self.out)

    def test_timeout_kills_term_ignoring_child_and_grandchild_preserves_certificate(self):
        self.fake_child("""
phase=out/'state/phases/cache';phase.mkdir(parents=True)
(phase/'certificate.json').write_text('{"fixture":"retain"}')
signal.signal(signal.SIGTERM,signal.SIG_IGN)
p=subprocess.Popen([sys.executable,'-c','import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)'])
(out/'grandchild.json').write_text(json.dumps({'pid':p.pid}))
while True:time.sleep(.05)
""")
        with self.assertRaises(TimeoutError):self.run_fake()
        report=json.loads((self.out/'supervisor-report.json').read_text())
        self.assertEqual(report['status'],'failed')
        self.assertEqual(report['termination_initiator'],'supervisor')
        self.assertEqual(report['remaining_live_child_group_pids'],[])
        self.assertEqual(report['child_returncode'],-signal.SIGKILL)
        self.assertTrue((self.out/'state/phases/cache/certificate.json').is_file())

    def test_sampled_memory_guard_stops_owned_group(self):
        self.fake_child('time.sleep(60)\n')
        original=n.processes325;calls=[]
        def oversized(*args,**kwargs):
            if not calls:
                calls.append(True)
                return {os.getpid():dict(group=os.getpgrp(),state='R',virtual_bytes=n.CAP_MEMORY325+1,rss_bytes=1)}
            return original(*args,**kwargs)
        with self.assertRaisesRegex(ValueError,'aggregate_virtual_memory_cap'):
            self.run_fake(processes325=oversized)
        self.assertEqual(json.loads((self.out/'supervisor-report.json').read_text())['remaining_live_child_group_pids'],[])

    def test_sampled_disk_guard_retains_failure_evidence(self):
        self.fake_child("(out/'large.tmp').write_bytes(b'x'*20000)\ntime.sleep(60)\n")
        with self.assertRaisesRegex(ValueError,'derived_byte_cap'):self.run_fake(CAP_BYTES325=10000)
        report=json.loads((self.out/'supervisor-report.json').read_text())
        self.assertEqual(report['status'],'failed');self.assertTrue((self.out/'large.tmp').exists())

    def test_child_signal_is_distinct_from_supervisor_deadline(self):
        self.fake_child('os.kill(os.getpid(),signal.SIGTERM)\n')
        with self.assertRaisesRegex(ValueError,'verification_child_failed'):self.run_fake()
        report=json.loads((self.out/'supervisor-report.json').read_text())
        self.assertEqual(report['child_returncode'],-signal.SIGTERM)
        self.assertEqual(report['signals'],[])
        self.assertEqual(report['termination_initiator'],'child_or_preflight')

    def test_unreapable_cleanup_still_records_failure(self):
        self.fake_child("""
import hashlib
hashes={}
for phase in ('cache','providers-C0','providers-L256','pair','work','preservation'):
 p=out/'state/phases'/phase/'certificate.json';p.parent.mkdir(parents=True);p.write_text('{}')
 hashes[phase]=hashlib.sha256(p.read_bytes()).hexdigest()
(out/'worker-complete.json').write_text(json.dumps({'passed':True,'phase_certificates':hashes}))
(out/'state.attempt.json').write_text(json.dumps({'status':'passed','native_gate_passed':True,'certificate_sha256':hashes}))
""")
        with self.assertRaisesRegex(ValueError,'supervisor_final_verdict'):
            self.run_fake(stop_group325=lambda child:None)
        report=json.loads((self.out/'supervisor-report.json').read_text())
        self.assertTrue(report['cleanup_wait_timed_out']);self.assertEqual(report['status'],'failed')

    def test_repeated_supervisor_term_preserves_first_signal_and_restores_handler(self):
        prior=signal.getsignal(signal.SIGTERM)
        self.fake_child("signal.signal(signal.SIGTERM,signal.SIG_IGN)\nos.kill(os.getppid(),signal.SIGTERM)\ntime.sleep(.01)\nos.kill(os.getppid(),signal.SIGTERM)\ntime.sleep(60)\n")
        with self.assertRaisesRegex(InterruptedError,'supervisor_received_SIGTERM'):self.run_fake()
        report=json.loads((self.out/'supervisor-report.json').read_text())
        self.assertEqual(len(report['signals']),1)
        self.assertEqual(report['signals'][0]['number'],signal.SIGTERM)
        self.assertEqual(signal.getsignal(signal.SIGTERM),prior)
        self.assertEqual(report['remaining_live_child_group_pids'],[])


if __name__=='__main__':unittest.main()
