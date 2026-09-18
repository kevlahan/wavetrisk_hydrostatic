#!/usr/bin/env python3
"""Read-only Stage215 archive audit; reselect the first full step after reload."""
import argparse,hashlib,importlib.util,json,re,statistics as st,tarfile,tempfile
from pathlib import Path
from collect_cluster_logs_215 import verify_archive
from build_restriction_208 import verify_helpers,digest
from analyze import analyze
from checkpoint_gate import accept
from cluster_199 import workload_gate
from cluster_215 import campaign,cases_for
from placement_200 import verify
from rank_metrics_190 import rank_resources
from protocol_215 import ordinary,POLICY

def read(p):return json.loads(p.read_text())
def require(test,message):
    if not test:raise ValueError(message)

def review(archive,kit):
    inventory=verify_archive(archive);verify_helpers(kit)
    require(inventory['stage']==215 and inventory['job']=='54714','Unexpected campaign identity')
    manifest=read(kit/'stage215-sources.json')
    spec=importlib.util.spec_from_file_location('frozen_protocol215',kit/'protocol_215.py')
    frozen=importlib.util.module_from_spec(spec);spec.loader.exec_module(frozen)
    report=dict(stage=215,job=54714,archive=str(archive),archive_sha256=digest(archive),
                inventoried_files_verified=len(inventory['files']),runs={},status='reviewing')
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(archive) as t:t.extractall(tmp,filter='data')
        root=Path(tmp);base=root/'results';r=read(base/'results.json');notes=r['build_notes']
        require(r['status']=='passed' and r['policy']==POLICY,'Campaign status/policy differs')
        require(notes==read(root/'builds/wt-stage215/restriction-215-builds.json'),'Build records differ')
        require(notes['manifest_sha256']==digest(kit/'stage215-sources.json'),'Frozen source manifest differs')
        for kind,b in notes['builds'].items():
            require(b['sources']==manifest['reference_sources' if kind=='baseline' else 'candidate_sources'],'Wrong build sources')
            log=root/f'builds/wt-stage215/build-215-{kind}.log'
            require(digest(log)==b['log_sha256'],'Build log hash differs')
            text=log.read_text();require('-ffp-contract=off' in text and '-DSTAGE202_RECORDS' in text and '-DSTAGE202_INVERSE' in text,'Build flags missing')
            require(('-DSTAGE215_AUDIT' in text)==(kind=='checked'),'Audit flag differs')
            require(not re.search(r'-DSTAGE203\w*',text),'Stage203 enabled')
            require('STAGE202=on' in b['command'] and 'STAGE203=off' in b['command'],'Build selectors differ')
        require(len(notes['physics_objects'])==23,'Physics inventory differs')
        require(r['cases']==cases_for(read(base/'topology.json')),'Placement map differs from topology')
        analyses={};medians={};anchors={};original_medians={};total_ranks=0
        for name,run in r['runs'].items():
            directory=base/name;a=analyze(directory/'run.log');logged=read(directory/'analysis.json')
            log=(directory/'run.log').read_text();require(a['completed'] and not a['failures'],'Solver log failed')
            require(a['step_records']==logged['step_records'] and len(a['step_records'])==12,'Raw steps differ')
            require(re.findall(r'Restarting from checkpoint\s+(\d+)',log)==['3','4'],'Reload coverage differs')
            require('Remapping vertical coordinates' in log.split('Restarting from checkpoint    4')[-1],'Post-reload remap absent')
            require(run['binary_sha256']==notes['builds'][run['kind']]['sha256']==r['executables'][run['kind']]['sha256'],'Binary identities differ')
            env=run['environment']
            for key in ['WAVETRISK_PROFILE_PARALLEL_BLOCKS','WAVETRISK_PROFILE_BLOCK_DETAIL','WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS']:
                require(env[key]=='0','Profiling active')
            require((env['WAVETRISK_VALIDATE_BLOCK_DYNAMICS']=='1')==run['oracle'],'Oracle setting differs')
            old=frozen.ordinary(log,a);require(old==run['ordinary'],'Frozen timing reconstruction differs')
            new=ordinary(log,a)
            require([v['index'] for v in new['excluded']]==[0,7,8],'Unexpected strict exclusions')
            require(len(new['included'])==9,'Wrong ordinary step count')
            medians[name]=new['median'];original_medians[name]=old['median']
            metrics=rank_resources(directory);require(metrics==run['rank_metrics'],'Rank resources differ')
            verify(metrics,r['cases'][run['case']]);total_ranks+=len(metrics)
            require(all(v['slurm']['SLURM_JOB_ID']=='54714' for v in metrics),'Rank job identity differs')
            case=run['case']
            if case not in anchors:anchors[case]=name
            analyses[name]=a
            require(workload_gate(analyses[anchors[case]],a)['passed'],'Within-count workload differs')
            require(workload_gate(analyses['n16-s1-validation-validate-A'],a)['passed'],'Cross-count workload differs')
            gate=r['gates'][name]
            require(gate['passed'] and not gate['failures'],'Recorded checkpoint gate failed')
            for cp in gate['checkpoints'].values():
                require(not accept(cp,exact=True),'Checkpoint report fails exact gate')
                require(all(v['different']==0 for v in cp['fields'].values()),'Nonidentical field report')
            require(r['gates']['work-'+name]['passed'],'Recorded workload failed')
            vm=run['memory_after']['vmstat'];before=run['memory_before']['vmstat']
            report['runs'][name]=dict(kind=run['kind'],phase=run['phase'],ranks=run['ranks'],
                archived_median=old['median'],corrected_ordinary=new,
                max_rank_rss_kib=max(v['max_rss_kib'] for v in metrics),
                sum_rank_peak_rss_kib=sum(v['max_rss_kib'] for v in metrics),
                swapin_pages=vm.get('pswpin',0)-before.get('pswpin',0),
                swapout_pages=vm.get('pswpout',0)-before.get('pswpout',0))
        require(r['gates']['cross-count']['passed'],'Cross-count tolerance report failed')
        for cp in r['gates']['cross-count']['checkpoints'].values():require(not accept(cp,exact=False),'Cross-count fixed tolerance exceeded')
        def replay(values):
            seen=[]
            def execute(case,kind,phase,suffix,oracle):
                name=f'{case}-{phase}-{suffix}'
                seen.append(dict(label=name,kind=kind,oracle=oracle));return values[name]
            result=campaign(execute,lambda _:None)
            require(seen==r['schedule'],'Scheduler decision changes executed sequence')
            return result
        recorded=replay(original_medians);require(recorded==r['decisions'],'Archived decisions do not reproduce')
        corrected=replay(medians);require(corrected['confirmation']['decision']=='meets_target','Corrected result not accepted')
        gains=corrected['confirmation']['paired_gains_seconds']
        report.update(status='passed',decision='meets_target',archived_decisions=r['decisions'],corrected_decisions=corrected,
            rank_run_records_verified=total_ranks,checkpoint_report_gates=len(r['runs']),workload_gates=len(r['runs']),
            cross_count_tolerance='passed',build_sources_and_logs='match frozen kit',physics_objects_recorded=23,
            protocol_sha256=digest(Path(__file__).with_name('protocol_215.py')),
            stdout_rounding_gain_intervals_seconds=[[g-.001,g+.001] for g in gains],
            correction='Frozen parser consumes restart flag on the checkpoint transaction. Corrected parser retains it for the next full integration step. Exclude 1,8,9 (one-based); nine retained steps/run. Original kit/results preserved.',
            limitations=['No checkpoint or binary payloads in review archive: checked stored field reports and identities, not absent payloads.',
                         'One ABBA batch, two independent runs/version; timesteps are not independent repetitions.',
                         'No contemporaneous legacy run and no mechanism ablation; total gain cannot be divided into packing, allocation and MPI waiting savings.'])
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('archive','kit','out'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();r=review(a.archive,a.kit);a.out.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ('status','decision','rank_run_records_verified','corrected_decisions')},indent=2))
