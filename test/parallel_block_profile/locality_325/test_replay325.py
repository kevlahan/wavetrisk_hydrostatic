"""Verification-only recovery preserves partial evidence without solver calls."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import verify325 as q
import verify322 as v

class ReplayReportTests(unittest.TestCase):
    def completed_fixture(self):
        def cache(c,directory,arm,scratch,started,stats,summaries,progress):
            summaries.extend({'rank':rank,'storage_orders':[q.expected_storage_order(rank,arm)]*5}
                             for rank in range(q.N))
        return (patch.object(q,'cache_run',side_effect=cache),
                patch.object(q,'route_certificates',return_value=({},{})),
                patch.object(q,'compare_features',return_value={'passed':True}))
    def test_final_success_report_write_failure_resets_passed(self):
        original=q.write_report325
        def fail_final(path,report):
            if report['passed']:raise OSError('injected final report save failure')
            original(path,report)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.json';a,b,c=self.completed_fixture()
            with a,b,c,patch.object(q,'write_report325',side_effect=fail_final):
                with self.assertRaisesRegex(OSError,'injected final report save failure'):
                    q.qualify(tmp,report_path=path,scratch_root=tmp)
            saved=json.loads(path.read_text())
            self.assertFalse(saved['passed']);self.assertEqual(saved['status'],'failed')
            self.assertEqual(saved['failure']['message'],'injected final report save failure')
            self.assertEqual(list(Path(tmp).iterdir()),[path])
    def test_secondary_failure_report_write_preserves_original_error(self):
        original=q.write_report325
        def fail_final_and_recovery(path,report):
            if report['passed']:raise OSError('original final report save failure')
            if report['status']=='failed':raise RuntimeError('secondary failure report save failure')
            original(path,report)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.json';a,b,c=self.completed_fixture()
            with a,b,c,patch.object(q,'write_report325',side_effect=fail_final_and_recovery):
                with self.assertRaisesRegex(OSError,'original final report save failure') as caught:
                    q.qualify(tmp,report_path=path,scratch_root=tmp)
            self.assertIn('secondary failure report save failure',' '.join(caught.exception.__notes__))
            self.assertFalse(json.loads(path.read_text())['passed'])
            self.assertEqual(list(Path(tmp).iterdir()),[path])
    def test_parse_failure_retains_completed_streams_and_exact_location(self):
        def cache(c,directory,arm,scratch,started,stats,summaries,progress):
            summaries.append({'rank':0,'stream_sha256':'retained-proof'})
            stats['sampled_peak_scratch_bytes']=8192
            progress('streams',arm,1,256,{'rank':1,'file':'execution316-rank-1.bin.gz'})
            raise v.AuditError('physical.rhs_components: phase4 slots999')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);captures=root/'captures';captures.mkdir()
            sentinel=captures/'unchanged';sentinel.write_bytes(b'original capture')
            recovery=root/'recovery';recovery.mkdir();report=recovery/'semantic.json'
            with patch.object(q,'cache_run',side_effect=cache):
                with self.assertRaisesRegex(v.AuditError,'rhs_components'):
                    q.qualify(captures,report_path=report,scratch_root=recovery)
            saved=json.loads(report.read_text())
            self.assertFalse(saved['passed']);self.assertEqual(saved['status'],'failed')
            self.assertEqual(saved['runs']['checked-L256']['streams'][0]['stream_sha256'],'retained-proof')
            self.assertEqual(saved['progress']['detail']['rank'],1)
            self.assertIn('phase4 slots999',saved['failure']['message'])
            self.assertEqual(saved['sampled_peak_scratch_bytes'],8192)
            self.assertEqual(list(recovery.iterdir()),[report])
            self.assertEqual(sentinel.read_bytes(),b'original capture')
    def test_provider_failure_retains_arm_streams_and_family(self):
        def cache(c,directory,arm,scratch,started,stats,summaries,progress):
            summaries.append({'rank':0,'stream_sha256':'retained-proof'})
        def providers(c,arm,progress):
            progress('providers',arm,7,140,{'kind':'routes','key':'phase504'},True)
            raise v.AuditError('provider mismatch')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.json'
            with patch.object(q,'cache_run',side_effect=cache),patch.object(q,'route_certificates',side_effect=providers):
                with self.assertRaisesRegex(v.AuditError,'provider mismatch'):
                    q.qualify(tmp,report_path=path,scratch_root=tmp)
            saved=json.loads(path.read_text())
            self.assertEqual(saved['progress'],{'stage':'providers','arm':'checked-L256','completed':7,'total':140,'detail':{'kind':'routes','key':'phase504'}})
            self.assertEqual(saved['runs']['checked-L256']['status'],'providers')
    def test_interruption_keeps_failure_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.json'
            with patch.object(q,'cache_run',side_effect=KeyboardInterrupt('bounded stop')):
                with self.assertRaises(KeyboardInterrupt):q.qualify(tmp,report_path=path,scratch_root=tmp)
            self.assertEqual(json.loads(path.read_text())['failure']['type'],'KeyboardInterrupt')
    def test_existing_call_interface_without_reports(self):
        with patch.object(q,'cache_run',side_effect=v.AuditError('expected fixture stop')):
            with self.assertRaisesRegex(v.AuditError,'expected fixture stop'):q.qualify('/fixture',None,None)
    def test_report_bound_preserves_previous_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.json';q.write_report325(path,{'passed':False,'status':'running'})
            before=path.read_bytes()
            with self.assertRaisesRegex(v.AuditError,'report_bound'):q.write_report325(path,{'large':'x'*(16*1024*1024)})
            self.assertEqual(path.read_bytes(),before)
            self.assertFalse(path.with_name(path.name+'.tmp').exists())

if __name__=='__main__':unittest.main()
