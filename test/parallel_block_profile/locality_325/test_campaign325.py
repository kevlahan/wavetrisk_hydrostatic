"""Offline failure/preservation checks; no compilation or solver launch."""
from pathlib import Path
from unittest import TestCase,main,mock
import hashlib,io,json,os,sys,tarfile,tempfile,types
import common,config325,control325,build,run325
class GateTests(TestCase):
 def test_schedule(self):
  self.assertEqual([(s['label'],s['ranks'],s['steps']) for s in config325.SCHEDULE],[('checked-L256',256,4),('qual-C0',256,4),('qual-L256',256,4)])
 def test_private_copy(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);a=p/'a';b=p/'b';a.write_bytes(b'original');a.chmod(0o444)
   self.assertTrue(common.private_copy_checkpoint(a,b)['private_copy']);self.assertEqual(a.read_bytes(),b.read_bytes());self.assertNotEqual(a.stat().st_ino,b.stat().st_ino)
   with self.assertRaises(ValueError):common.private_copy_checkpoint(a,b)
 def test_checkpoint_symlink_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'a').write_bytes(b'yes');(p/'link').symlink_to(p/'a')
   with self.assertRaises(ValueError):common.private_copy_checkpoint(p/'link',p/'dest')
 def test_environment(self):
  verifier=types.SimpleNamespace(capture_environment=lambda *x:{'WAVETRISK_EXECUTION_AUDIT316':'1','WAVETRISK_COMPACT_AUDIT325':'1'})
  with mock.patch.dict(sys.modules,verify325=verifier),mock.patch.dict(os.environ,{'WAVETRISK_RHS_PAIR_PROBE':'1','WAVETRISK_WINDOW322':'1','OMP_NUM_THREADS':'16'}):env=config325.environment(Path(__file__).parent,'candidate','L256')
  self.assertEqual(env['WAVETRISK_WINDOW322'],'0');self.assertEqual(env['WAVETRISK_WORK_LEDGER325'],'1');self.assertEqual(env['WAVETRISK_RHS_PAIR_PROBE'],'0');self.assertEqual(env['OMP_NUM_THREADS'],'1')
 def test_unbounded_capture_rejected(self):
  with mock.patch.dict(sys.modules,verify325=types.SimpleNamespace(capture_environment=lambda *x:{'WAVETRISK_EXECUTION_AUDIT316':'1'})):
   with self.assertRaises(ValueError):config325.environment(Path(__file__).parent,'candidate','L256')
 def test_invalid_arms(self):
  for kind,owner in [('candidate','B256'),('main','off'),('checked','off')]:
   with self.assertRaises(ValueError):config325.environment('.',kind,owner)
 def test_completed_squeue_error(self):
  answers=[types.SimpleNamespace(returncode=1,stdout='',stderr='Invalid job id'),types.SimpleNamespace(returncode=0,stdout='55570|COMPLETED|\n',stderr='')]
  with mock.patch.object(control325.subprocess,'run',side_effect=answers):self.assertFalse(control325.queued(55570))
 def test_unknown_scheduler(self):
  with mock.patch.object(control325.subprocess,'run',return_value=types.SimpleNamespace(returncode=1,stdout='',stderr='offline')):
   with self.assertRaises(RuntimeError):control325.queued(55570)
 def test_boolean_job(self):
  with self.assertRaises(ValueError):control325.queued(True)
 def roots(self,p):
  for name in ('kit','preflight','results'):(p/name).mkdir()
  return dict(kit=str(p/'kit'),preflight=str(p/'preflight'),out=str(p/'results'))
 def test_excludes_large_evidence_without_reading(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);state=self.roots(p);r=p/'results';(r/'qual-L256').mkdir();(r/'qual-L256/execution316-rank-0.bin.gz').write_bytes(b'raw');(r/'qual-L256/test_checkpoint_9901.bin.zst').write_bytes(b'checkpoint');(r/'results.json').write_text('{}');(r/'sources').mkdir();(r/'sources/x.json').write_text('{}')
   entries,excluded,total=control325.export_files(state);self.assertEqual([e[1] for e in entries],['results/results.json']);self.assertEqual(len(excluded),2);self.assertEqual(total,2)
 def test_excludes_raw_work_and_full_peer_vectors_without_read(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);state=self.roots(p);d=p/'results/qual-L256';d.mkdir()
   excluded_names=('work325-rank-0.tsv','work325-rank-vectors.json')
   for name in excluded_names:
    with (d/name).open('wb') as f:f.truncate(control325.MAX_MEMBER+1)
   (d/'work325-review-rank-0.json').write_text('{}')
   original=Path.read_bytes
   def guarded(path):
    if path.name in excluded_names:raise AssertionError('Excluded native work was read')
    return original(path)
   with mock.patch.object(Path,'read_bytes',guarded):entries,excluded,total=control325.export_files(state)
   self.assertEqual([e[1] for e in entries],['results/qual-L256/work325-review-rank-0.json']);self.assertEqual(total,2)
   self.assertEqual({Path(e['path']).name for e in excluded},set(excluded_names));self.assertTrue(all(e['reason']=='native_artifact_retained_no_read' for e in excluded))
 def test_collection_member_cap(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);state=self.roots(p)
   with (p/'results/large.json').open('wb') as f:f.truncate(control325.MAX_MEMBER+1)
   with self.assertRaises(ValueError):control325.export_files(state)
 def test_collection_total_cap(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);state=self.roots(p);(p/'results/a.json').write_bytes(b'1234');(p/'results/b.json').write_bytes(b'5678')
   with mock.patch.object(control325,'MAX_TOTAL',7):
    with self.assertRaises(ValueError):control325.export_files(state)
 def test_no_symlink_follow(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);state=self.roots(p);(p/'secret').write_bytes(b'secret');(p/'results/link.json').symlink_to(p/'secret')
   entries,excluded,total=control325.export_files(state);self.assertEqual(entries,[]);self.assertEqual(total,0);self.assertEqual(excluded[0]['reason'],'symlink_not_followed')
 def test_compact_archive_roundtrip(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);state=self.roots(p);state.update(job=None,status='preflight_failed',run_id='test',work=str(p))
   (p/'kit/files-sha256.json').write_text('{}');(p/'results/results.json').write_text('{"status":"failed"}');(p/'results/test_checkpoint_9901.bin.zst').write_bytes(b'keep native')
   with mock.patch.object(control325.Path,'home',return_value=p),mock.patch.object(control325.subprocess,'run',return_value=types.SimpleNamespace(returncode=0)):
    archive=control325.collect(state)
   self.assertTrue(archive.is_file());self.assertEqual((p/'results/test_checkpoint_9901.bin.zst').read_bytes(),b'keep native')
   with tarfile.open(archive) as t:self.assertFalse(any('checkpoint' in m.name for m in t.getmembers()))
   self.assertFalse(any(p.glob('rk325-collect-*')));self.assertGreater(control325.verify_archive(archive),4)
 def test_archive_traversal(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'source-files.json').write_text('{"candidate":{}}')
   with tarfile.open(p/'sources.tar.gz','w:gz') as t:
    info=tarfile.TarInfo('candidate/../../escape');info.size=1;t.addfile(info,io.BytesIO(b'x'))
   with self.assertRaises(ValueError):build.extract_sources(p,p/'sources')
 def test_archive_safe_link(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);content=b'program';link='../test/program.f90';manifest={'candidate':{'test/program.f90':hashlib.sha256(content).hexdigest(),'src/test.f90':hashlib.sha256(link.encode()).hexdigest()}};(p/'source-files.json').write_text(json.dumps(manifest))
   with tarfile.open(p/'sources.tar.gz','w:gz') as t:
    info=tarfile.TarInfo('candidate/test/program.f90');info.size=len(content);t.addfile(info,io.BytesIO(content));info=tarfile.TarInfo('candidate/src/test.f90');info.type=tarfile.SYMTYPE;info.linkname=link;t.addfile(info)
   tree,_=build.extract_sources(p,p/'sources');self.assertEqual((tree/'src/test.f90').read_bytes(),content)
 def test_allocation_scope(self):
  with mock.patch.dict(os.environ,{},clear=True):
   with self.assertRaises(ValueError):run325.require_allocation()
 def test_capture_monitor_run_cap(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);d=p/'qual-C0';d.mkdir();(d/'execution316-rank-0.bin.gz').write_bytes(b'12345')
   with mock.patch.object(run325,'CAPTURE_RUN_BYTES',4):
    with self.assertRaises(ValueError):run325.capture_usage(p,d)
 def test_capture_monitor_campaign_cap(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)
   for name in ('qual-C0','qual-L256'):
    d=p/name;d.mkdir();(d/'execution316-rank-0.bin.gz').write_bytes(b'123')
   with mock.patch.object(run325,'CAPTURE_CAMPAIGN_BYTES',5):
    with self.assertRaises(ValueError):run325.capture_usage(p,p/'qual-L256')
 def test_capture_monitor_kills_srun_on_threshold(self):
  process=types.SimpleNamespace(pid=12345,poll=mock.Mock(return_value=None),wait=mock.Mock(return_value=-15))
  with mock.patch.object(run325.subprocess,'Popen',return_value=process),mock.patch.object(run325,'capture_usage',side_effect=[{},ValueError('capture cap')]),mock.patch.object(run325.os,'killpg') as kill:
   with self.assertRaisesRegex(ValueError,'capture cap'):run325.run_bounded(['srun'],'run',{},io.StringIO(),'results')
   kill.assert_called_once_with(12345,run325.signal.SIGTERM);process.wait.assert_called_once_with(timeout=10)
 def test_toolchain_mutation(self):
  keys=('wrapper_resolved','wrapper_sha256','compiler_resolved','compiler_sha256','compiler_version','c_compiler_resolved','c_compiler_sha256','c_compiler_version','target','wrapper_command','expected','library_version');old={k:k for k in keys};new=dict(old);new['compiler_version']='14.0'
  with self.assertRaises(ValueError):run325.bind_prior_toolchain(new,{'original_toolchain':old})
 def test_extra_package_file(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'good').write_bytes(b'ok');(p/'files-sha256.json').write_text(json.dumps({'good':common.sha(p/'good')}));common.verify(p);(p/'unexpected').write_bytes(b'bad')
   with self.assertRaises(ValueError):common.verify(p)
if __name__=='__main__':main()
