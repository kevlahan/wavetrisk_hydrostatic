"""Offline read-only binding and recovery failure tests."""
from pathlib import Path
from unittest import TestCase,main,mock
import gzip,json,os,tempfile,types,shutil
import recover325 as r
import control_recovery325 as ctl
class RecoveryTests(TestCase):
 def capture(self,root):
  p=root/'checked-L256/execution316-rank-0.bin.gz';p.parent.mkdir();p.write_bytes(gzip.compress(b'logical capture bytes'*100,mtime=0));st=p.stat()
  return p,dict(relative=str(p.relative_to(root)),stored_bytes=st.st_size,mtime_ns=st.st_mtime_ns,sha256=r.sha(p),family='capture_inventory')
 def test_actual_read_matches_anchor(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root);seen={}
   with r.verified_stream(p,root,{row['relative']:row},seen) as f:self.assertEqual(f.read(),b'logical capture bytes'*100)
   self.assertEqual(seen[row['relative']]['compressed_sha256'],row['sha256'])
 def test_same_size_and_timestamp_mutation_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root);data=bytearray(p.read_bytes());data[9]=1 if data[9]!=1 else 2;p.write_bytes(data);os.utime(p,ns=(row['mtime_ns'],row['mtime_ns']))
   with self.assertRaisesRegex(ValueError,'changed between binding'):
    with r.verified_stream(p,root,{row['relative']:row},{}) as f:f.read()
 def test_historical_anchor_accepts_unchanged_bytes(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root)
   self.assertEqual(r.hash_bound_input(p,row,{row['relative']:row}),row['sha256'])
 def test_historical_anchor_rejects_changed_bytes_with_original_metadata(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root);data=bytearray(p.read_bytes());data[9]=1 if data[9]!=1 else 2;p.write_bytes(data);os.utime(p,ns=(row['mtime_ns'],row['mtime_ns']))
   with self.assertRaisesRegex(ValueError,'job55571 historical anchor'):r.hash_bound_input(p,row,{row['relative']:row})
 def test_historical_anchor_missing_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root)
   with self.assertRaisesRegex(ValueError,'job55571 historical anchor'):r.hash_bound_input(p,row,{})
 def test_historical_certificate_and_inventory(self):
  kit=Path(__file__).resolve().parent;binding=json.loads((kit/'binding325.json').read_text())
  history,anchors=r.load_history(kit,binding);self.assertEqual(len(anchors),1536);self.assertEqual(history['cumulative_maximum_seconds'],2681)
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)
   for name in ('recovery-history325.json','prior-native-input-anchors.json'):shutil.copy2(kit/name,p/name)
   (p/'prior-native-input-anchors.json').write_text('{}')
   with self.assertRaisesRegex(ValueError,'certificate changed'):r.load_history(p,binding)
 def test_history_mutations_rejected(self):
  kit=Path(__file__).resolve().parent;binding=json.loads((kit/'binding325.json').read_text())
  original=json.loads((kit/'recovery-history325.json').read_text())
  import copy
  for mutation in ('budget','earlier_cost','predecessor','scope','accounting','anchors'):
   with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp);h=copy.deepcopy(original)
    if mutation=='budget':h['recovery_walltime_cap_seconds']=2100
    elif mutation=='earlier_cost':h['earlier_replay']['predecessor_elapsed_seconds']=0
    elif mutation=='predecessor':h['predecessor_job']=55571
    elif mutation=='scope':h['predecessor_report']['solver_runs']=1
    elif mutation=='accounting':h['predecessor_accounting']['stdout']=h['predecessor_accounting']['stdout'].replace('|476|','|475|')
    elif mutation=='anchors':h['earlier_replay']['anchors_sha256']='0'*64
    (root/'recovery-history325.json').write_text(json.dumps(h))
    shutil.copy2(kit/'prior-native-input-anchors.json',root/'prior-native-input-anchors.json')
    with self.assertRaises(ValueError):r.load_history(root,binding)
 def test_native_limits_and_budget(self):
  kit=Path(__file__).resolve().parent
  self.assertIn('--time=00:30:00',(kit/'verification.sbatch').read_text())
  self.assertIn('signal.alarm(29*60)',(kit/'recover325.py').read_text())
  self.assertEqual(392+13+476+1800,2681)
 def test_parser_error_not_masked(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root);seen={}
   with self.assertRaisesRegex(ValueError,'original parser failure'):
    with r.verified_stream(p,root,{row['relative']:row},seen):raise ValueError('original parser failure')
   self.assertEqual(seen,{})
 def test_changed_mtime_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root);os.utime(p,ns=(row['mtime_ns']+1000000,row['mtime_ns']+1000000))
   with self.assertRaisesRegex(ValueError,'inventory changed'):r.stat_match(p,row)
 def test_unbound_capture_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root)
   with self.assertRaisesRegex(ValueError,'Unbound'):
    with r.verified_stream(p,root,{},{}):pass
 def test_symlink_capture_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p,row=self.capture(root);link=root/'link';link.symlink_to(p)
   with self.assertRaisesRegex(ValueError,'regular native'):r.stat_match(link,row)
 def test_shadow_work_new_output_only(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);native=root/'native';native.mkdir();out=root/'out';out.mkdir();d=native/'qual-C0';d.mkdir();p=d/'work325-rank-0.tsv';p.write_bytes(b'raw ledger');st=p.stat();row=dict(relative='qual-C0/work325-rank-0.tsv',stored_bytes=st.st_size,mtime_ns=st.st_mtime_ns)
   shadow=r.shadow_work(dict(original_results=str(native),work_inventory=[row]),{row['relative']:{**row,'sha256':r.sha(p)}},out)
   self.assertTrue((shadow/row['relative']).is_symlink());(shadow/'qual-C0/work325-rank-vectors.json').write_text('{}');self.assertFalse((native/'qual-C0/work325-rank-vectors.json').exists());self.assertEqual(p.read_bytes(),b'raw ledger')
 def test_package_extra_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'a').write_bytes(b'a');(p/'files-sha256.json').write_text(json.dumps({'a':r.sha(p/'a')}));r.package_check(p);(p/'extra').write_bytes(b'x')
   with self.assertRaisesRegex(ValueError,'inventory'):r.package_check(p)
 def test_escaping_paths(self):
  with tempfile.TemporaryDirectory() as tmp:
   for name in ('../else','/tmp/else'):
    with self.assertRaises(ValueError):r.safe(Path(tmp),name)
 def test_allocation_scope(self):
  with mock.patch.dict(os.environ,{},clear=True):
   with self.assertRaisesRegex(ValueError,'one-task'):r.require_allocation()
 def test_allocation_requires_both_memory_limits(self):
  env=dict(SLURM_JOB_ID='55572',SLURM_JOB_PARTITION='nk',SLURM_JOB_NUM_NODES='1',SLURM_NTASKS='1')
  with mock.patch.dict(os.environ,env,clear=True),mock.patch.object(r.socket,'gethostname',return_value='bb13.rhpcs.internal'):
   with mock.patch.object(r.resource,'getrlimit',return_value=(16*1024**3,16*1024**3)):r.require_allocation()
   for caps in ((-1,-1),(16*1024**3,-1)):
    with mock.patch.object(r.resource,'getrlimit',return_value=caps):
     with self.assertRaisesRegex(ValueError,'hard and soft'):r.require_allocation()
 def test_partial_inventory_fails_before_read(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);native=p/'native';native.mkdir();out=p/'out';out.mkdir()
   for name in r.LABELS:(native/name).mkdir()
   with self.assertRaisesRegex(ValueError,'inventory differs'):r.anchor_inputs(dict(original_results=str(native),capture_inventory=[],work_inventory=[]),out,{})
 def test_collector_excludes_work_links_and_vectors(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);kit=p/'kit';kit.mkdir();out=p/'out';out.mkdir();native=p/'original.tsv';native.write_bytes(b'native');(out/'input.tsv').symlink_to(native)
   with (out/'work325-rank-vectors.json').open('wb') as f:f.truncate(17*1024**2)
   (out/'report.json').write_text('{}');files,total=ctl.selected_files(dict(kit=str(kit),out=str(out)))
   self.assertEqual([x[1] for x in files],['results/report.json']);self.assertEqual(total,2);self.assertEqual(native.read_bytes(),b'native')
 def test_collector_cap(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'kit').mkdir();(p/'out').mkdir()
   with (p/'out/large.json').open('wb') as f:f.truncate(17*1024**2)
   with self.assertRaisesRegex(ValueError,'exceeds16MiB'):ctl.selected_files(dict(kit=str(p/'kit'),out=str(p/'out')))
 def test_scheduler_terminal_with_squeue_error(self):
  answers=[types.SimpleNamespace(returncode=1,stdout='',stderr='Invalid job id'),types.SimpleNamespace(returncode=0,stdout='55569|FAILED|\n',stderr='')]
  with mock.patch.object(ctl.subprocess,'run',side_effect=answers):self.assertFalse(ctl.queued(55569))
 def test_unknown_scheduler(self):
  with mock.patch.object(ctl.subprocess,'run',return_value=types.SimpleNamespace(returncode=1,stdout='',stderr='offline')):
   with self.assertRaises(RuntimeError):ctl.queued(55569)
if __name__=='__main__':main()
