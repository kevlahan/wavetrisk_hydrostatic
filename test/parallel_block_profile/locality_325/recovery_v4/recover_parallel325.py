"""One bounded continuation of immutable job55569 capture verification.

The checked-arm native closure is reused from55573. Both optimized arms are
recomputed, including new physical inverse invariants and provider certificates.
No compiler, MPI launcher or solver path exists in this module.
"""
from pathlib import Path
import json,os,resource,signal,socket,sys,time,shutil,platform
from recover325 import require,save,sha,safe,package_check,metadata_check,load_history,anchor_inputs,shadow_work,stat_match,LABELS
import parallel_gate325 as gate

def history_check(kit,binding):
 earlier,anchors=load_history(kit,binding)
 h=json.loads((kit/'parallel-history325.json').read_text())
 require(h['source_job']==55569 and h['predecessor_job']==55573 and h['predecessor_elapsed_seconds']==1511,'Unexpected final predecessor')
 require(h['prior_elapsed_seconds']==2392 and h['walltime_cap_seconds']==300 and h['cumulative_maximum_seconds']==2692 and 392+13+476+1511+300==2692<=2700,'Cumulative budget mismatch')
 for name,hsh in h['certificate_sha256'].items():
  require(sha(safe(kit,name))==hsh,'Prior verification certificate changed')
 prior=json.loads((kit/'predecessor-v3/recovery-report.json').read_text())
 semantic=json.loads((kit/'predecessor-v3/semantic-replay.json').read_text())
 actual=json.loads((kit/'predecessor-v3/capture-read-checks-partial.json').read_text())
 require(prior['verification_job']=='55573' and prior['source_job']==55569 and prior['status']=='failed' and prior['builds']==prior['solver_runs']==prior['performance_runs']==0,'Prior scope mismatch')
 require(prior['historical_input_hashes_matched']==1536 and len(actual)==768,'Prior byte bindings incomplete')
 require(semantic['failure']['type']=='AuditError' and semantic['failure']['message'].startswith('physical.ordered_program_or_value_difference:') and semantic['progress']['stage']=='cross-arm','Prior failure not the reviewed cross-arm assertion')
 require(prior['package_manifest_sha256']=='386333bfa0bb686f4dbec9cdc8f940df50d8c020dc02c76ff5d2a8e7003d70ad','Unexpected prior verifier package')
 for label in LABELS:
  run=semantic['runs'][label]
  require(run['status']=='complete' and run['families']==254 and [x['rank'] for x in run['streams']]==list(range(256)),'Prior per-arm closure incomplete')
 for name,row in actual.items():
  require(name in anchors and row==dict(compressed_sha256=anchors[name]['sha256'],stored_bytes=anchors[name]['stored_bytes'],complete=True),'Prior actual-read binding differs')
 accounting=json.loads((kit/'predecessor-v3/accounting.json').read_text())
 require(accounting['returncode']==0 and '|FAILED|' in accounting['stdout'] and '|1511|' in accounting['stdout'],'Prior scheduler cost differs')
 # Byte-bind the remote original certificates as well as the packaged copies.
 for name,hsh in h['remote_predecessor_files'].items():require(sha(Path(name))==hsh,'Original55573certificate changed')
 return h,anchors,semantic

def allocation():
 require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_JOB_PARTITION')=='nk' and os.environ.get('SLURM_JOB_NUM_NODES')=='1' and os.environ.get('SLURM_NTASKS')=='1' and os.environ.get('SLURM_CPUS_PER_TASK')=='16' and socket.gethostname().split('.')[0]=='bb13','Require one-task16CPUbb13/nk allocation')
 require(len(os.sched_getaffinity(0))==16,'Require16allocated CPUs')
 require(resource.getrlimit(resource.RLIMIT_AS)==(16*1024**3,16*1024**3),'Require16GiBprocess hard/soft VM limits')

def run(kit,out):
 sys.dont_write_bytecode=True
 kit=Path(kit).resolve(strict=True);package_sha=package_check(kit)
 binding=json.loads((kit/'binding325.json').read_text());history,historical,prior_semantic=history_check(kit,binding);allocation()
 out=Path(out).resolve();require(not out.exists(),'Fresh continuation output required')
 original=Path(binding['original_results']).resolve(strict=True)
 for protected in (original,Path(binding['original_kit']).resolve(strict=True),kit):require(not out.is_relative_to(protected) and not protected.is_relative_to(out),'Output overlaps protected input')
 out.mkdir();scratch=out/'scratch';scratch.mkdir();started=time.monotonic();stats={}
 report=dict(verifier_environment=dict(python=platform.python_version(),executable=sys.executable,hostname=socket.gethostname(),allowed_cpus=sorted(os.sched_getaffinity(0)),process_virtual_memory_limits=list(resource.getrlimit(resource.RLIMIT_AS)),aggregate_virtual_memory_limit_bytes=gate.CAP_MEMORY),stage=325,revision='verification-recovery-v4',status='running',source_job=55569,verification_job=os.environ['SLURM_JOB_ID'],predecessor_job=55573,builds=0,solver_runs=0,performance_runs=0,performance_qualified=False,walltime_cap_seconds=300,runner_alarm_seconds=270,prior_elapsed_seconds=2392,cumulative_maximum_seconds=2692,package_manifest_sha256=package_sha,verifier_sha256=sha(kit/'verify325.py'),parallelism=16,checked_arm_reuse_policy='55573hash-bound full256stream/254provider checkedclosure retained; new optimized physical-pair proof is recomputed, not inferred from earlier counts.')
 def persist():
  report.update(stats);report['elapsed_seconds']=time.monotonic()-started;save(out/'recovery-report.json',report)
 def stop(signum,frame):raise TimeoutError('Bounded five-minute continuation ending; no retry or budget extension')
 signal.signal(signal.SIGALRM,stop);signal.signal(signal.SIGTERM,stop);signal.alarm(270);persist()
 try:
  report['phase']='bind_original_source_metadata_and_inputs';persist()
  report['binding_before']=metadata_check(binding)
  anchors=anchor_inputs(binding,out,historical);report['historical_input_hashes_matched']=len(anchors)
  # Bridge proof uses exact inspected source and the fixed5/5configured schedule.
  import inverse_identity325
  manifest=json.loads((Path(binding['original_kit'])/'source-files.json').read_text())
  report['inverse_bridge_source_binding']=inverse_identity325.bind_inverse_bridge_sources325(manifest,source_root=original/'sources/candidate')
  report['checked_arm_reuse']=dict(passed=True,source_job=55573,streams=256,families=254,report_sha256=sha(kit/'predecessor-v3/semantic-replay.json'),actual_capture_hashes_matched=256)
  report['phase']='parallel_optimized_streams';persist()
  command=[sys.executable,'-u',str(kit/'parallel_streams325.py'),str(kit),str(original),str(scratch),str(out/'parallel-streams.json')]
  gate.supervise([(command,out/'parallel-streams.log')],scratch,stats,persist)
  stream_report=json.loads((out/'parallel-streams.json').read_text());require(stream_report['passed'] is True,'Optimized stream/shard continuation incomplete')
  report['optimized_streams']=dict(passed=True,report_sha256=sha(out/'parallel-streams.json'),streams=512);persist()
  observed=stream_report['actual_capture_reads'];require(len(observed)==512,'Optimized actual-read census incomplete')
  for name,row in observed.items():
   require(anchors[name]['family']=='capture_inventory' and row==dict(compressed_sha256=anchors[name]['sha256'],stored_bytes=anchors[name]['stored_bytes'],complete=True),'Changed optimized actual-read input')
  report['phase']='parallel_optimized_providers';persist();commands=[]
  for arm in ('qual-C0','qual-L256'):
   commands.append(([sys.executable,'-u',str(kit/'parallel_gate325.py'),'providers',str(scratch/(arm+'.sqlite')),arm,str(scratch/(arm+'-proof.json')),str(out/(arm+'-provider-progress.json'))],out/(arm+'-providers.log')))
  gate.supervise(commands,scratch,stats,persist)
  report['phase']='complete_physical_cross_arm_comparison';persist()
  for arm in ('qual-C0','qual-L256'):
   obj=json.loads((scratch/(arm+'-proof.json')).read_text());save(out/(arm+'-provider-certificates.json'),dict(arm=arm,proofs=obj['proofs'],proofs_sha256=obj['proofs_sha256'],counts=obj['counts'],features_sha256=obj['features_sha256'],feature_count=obj['feature_count']))
  try:
   pair=gate.compare_payloads(scratch/'qual-C0-proof.json',scratch/'qual-L256-proof.json',out/'cross-arm.json')
  except ValueError as exc:
   if not str(exc).startswith('Cross-arm physical feature/provider mismatch'):raise
   pair=json.loads((out/'cross-arm.json').read_text())
  report['physical_pair']=dict(passed=pair['passed'],report_sha256=sha(out/'cross-arm.json'));persist()
  # Semantic files are no longer needed; retain bounded certificates and summaries.
  shutil.rmtree(scratch);scratch.mkdir()
  report['phase']='work_replay_new_output_only';persist();shadow=shadow_work(binding,anchors,out)
  from account325 import qualify
  work=qualify(shadow,kit);require(work['passed'] is True,'Work replay failed');save(out/'work-replay.json',work)
  count=0
  for arm in LABELS:
   vectors=json.loads((shadow/arm/'work325-rank-vectors.json').read_text())['per_rank']
   require([x['rank'] for x in vectors]==list(range(256)),'Work replay rank census incomplete')
   for rank,row in enumerate(vectors):
    anchor=anchors[f'{arm}/work325-rank-{rank}.tsv'];require(row['sha256']==anchor['sha256'] and row['bytes']==anchor['stored_bytes'],'Work actual read differs');count+=1
  require(count==768,'Work actual-read census incomplete');report['work_actual_read_hashes_matched']=count;report['work_replay']=dict(passed=True,report_sha256=sha(out/'work-replay.json'));persist()
  report['phase']='final_preservation';persist();report['binding_after']=metadata_check(binding);history_check(kit,binding)
  for row in binding['capture_inventory']+binding['work_inventory']:stat_match(safe(original,row['relative']),row)
  for row in binding['work_inventory']:require(sha(safe(original,row['relative']))==anchors[row['relative']]['sha256'],'Work changed after verification')
  require(pair['passed'],'Cross-arm physical pair failed; work and preservation reports retained')
  report.update(status='passed',phase='complete',decision='CORRECTED_PHYSICAL_PROVIDER_AND_WORK_GATE_PASSED_REQUIRES_INDEPENDENT_REVIEW_NO_PERFORMANCE',original_ledger_preserved=True);persist()
 except BaseException as exc:
  # Preserve independent post-failure evidence even when the semantic gate fails.
  # The270s alarm leaves30s inside the300s allocation for shutdown/reporting.
  try:
   report['failure_preservation']=metadata_check(binding)
   for row in binding['capture_inventory']+binding['work_inventory']:
    from recover325 import hash_bound_input
    hash_bound_input(safe(original,row['relative']),row,historical)
   report['failure_input_hashes_matched']=1536
  except BaseException as preserve_error:report['failure_preservation_error']=repr(preserve_error)
  report.update(status='failed',error=repr(exc),decision='STOP_NO_AUTOMATIC_RETRY_OR_BUDGET_EXTENSION');persist();raise
 finally:
  signal.alarm(0);shutil.rmtree(scratch,ignore_errors=True)
if __name__=='__main__':
 require(len(sys.argv)==3,'Usage: recover_parallel325.py KIT NEW_OUTPUT');run(*sys.argv[1:])
