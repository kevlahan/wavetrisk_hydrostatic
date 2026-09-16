#!/usr/bin/env python3
"""Local restart/correctness gate for the frozen focused observers, never timing."""
import argparse,json
from pathlib import Path
from analyze_lifecycle_205 import collect,attach_sites
from analyze_hardware_201 import collect as hardware
from analyze_scaling_200 import collect as scaling
from local_j4 import compare_windows
from cluster_199 import workload_gate
from compare_restart_checkpoints import compare
from checkpoint_gate import accept
from build_restriction_205 import verify,digest

def validate(a):
    manifest=json.loads((a.kit/'stage205-sources.json').read_text());sites=json.loads((a.kit/'operation-sites-205.json').read_text())
    out={'stage':205,'status':'passed','scope':'Four-rank J4/J6 numerical and observer validation only. No local speed claim; cluster calibration remains mandatory.','runs':{}}
    for side,key in [('baseline','baseline_profile_sources'),('legacy','legacy_profile_sources')]:
        ref=getattr(a,'reference_'+side);rr=json.loads((ref/'results.json').read_text())
        for mode in (0,1):
            root=getattr(a,side+('_off' if mode==0 else '_on'));r=json.loads((root/'results.json').read_text())
            if r['status']!='passed' or r['ranks']!=4 or r['environment'].get('WAVETRISK_LIFECYCLE205')!=str(mode):raise ValueError('Local mode/configuration differs')
            verify(Path(r['build']),manifest[key])
            if digest(root/'climateJ4')!=r['binary_sha256']:raise ValueError('Binary identity differs')
            equality=compare_windows(ref,root,'localJ4J6_checkpoint_0002.bin.zst',expected_new=[3])
            if not equality['passed'] or not workload_gate(rr['analysis'],r['analysis'])['passed']:raise ValueError('Observer changed numerical workload')
            q=collect(root,4,mode);attach_sites(q,sites,side);h=hardware(root,4,0);sc=scaling(root,4,0)
            if h['quality_flags'] or len(sc['workflow']['steps'])!=12:raise ValueError('Clock or step coverage differs')
            out['runs'][side+'-'+str(mode)]={'path':str(root),'mode':mode,'steps':12,'ordinary_steps':len(q['ordinary_steps']),
                'exact_to_unmodified_reference':str(ref),'reloads':r['reloads'],'binary_sha256':r['binary_sha256'],
                'source_inventory':manifest[key],'event_counts':q['lifecycle']['event_counts'],'plans':q['lifecycle']['plans'],
                'collective_sites':len(q['lifecycle']['collective_sites']),
                'max_rank_events':max(map(len,q['lifecycle']['rank_events'].values())),
                'observer_checks':'exclusive/inclusive bounds, common-context conservation, MPI-owner and collective-site duration conservation; restart excluded from ordinary accounting',
                'build_log_sha256':digest(Path('/private/tmp/stage205-build-'+side+'.log'))}
    report=compare(a.baseline_on/'localJ4J6_checkpoint_0003.bin.zst',a.legacy_on/'localJ4J6_checkpoint_0003.bin.zst',expected_domains=40)
    if accept(report,exact=False):raise ValueError('Legacy fixed-tolerance screen failed')
    out['legacy_screen']='passed unchanged fixed tolerance (not bitwise)'
    out['tests']=['Actual Fortran O0/fcheck=all observer: nested scopes, conditional call injection, disabled mode, lifecycle events, site timing conservation and corruption rejection.',
      'Synthetic adjacent/reversed calibration and signed comparison; missing controls rejected.',
      'Real six-build frozen patch/hash/submit protocol with compiler stand-in; modified binary rejected.']
    out['excluded_attempts']=['Initial sandboxed MPI launch failed before solver initialization; successful runs used local socket permission.',
      'Initial parallel local builds exposed existing incomplete Makefile dependencies; final builds use the frozen serial -j1 protocol.',
      'Initial generated event lines exceeded 132 columns; corrected before successful solver tests.']
    a.out.write_text(json.dumps(out,indent=2)+'\n');print('PASS: four local observer on/off runs, exact restart equality, lifecycle coverage, clock conservation and fixed legacy screen')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('kit','baseline-on','baseline-off','legacy-on','legacy-off','reference-baseline','reference-legacy','out'):p.add_argument('--'+n,type=Path,required=True)
    validate(p.parse_args())
