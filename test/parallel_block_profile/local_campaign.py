#!/usr/bin/env python3
"""Gate a lighter local timing experiment against unchanged legacy checkpoints.

Validation includes two write/reload cycles and remapping, separately from
alternating timing pairs. No 83-rank forecast and no bitwise-equality claim.
"""
import argparse
import json
from pathlib import Path
import platform
import statistics
from build_legacy import tracked_entries, verify
from checkpoint_gate import compare_runs, POLICY
from experiment import LEGACY_REF, digest
from profile_protocol import checkpoint_name, environment, execute, fixture_identity


def save(path, value):
    path.write_text(json.dumps(value,indent=2)+'\n')


def binding(binaries, fixture, ranks):
    return {'binaries':{k:digest(v) for k,v in binaries.items()},
            'fixture':fixture_identity(fixture),'ranks':ranks,
            'threads':{k:v for k,v in environment().items()
                       if k.endswith('NUM_THREADS') or k=='VECLIB_MAXIMUM_THREADS'}}


def require_certificate(certificate, current, time_end):
    if certificate.get('passed') is not True or certificate.get('policy')!=POLICY:
        raise ValueError('Missing or incompatible full-field validation gate')
    if certificate.get('binding')!=current:
        raise ValueError('Validation identity differs: binary, fixture, rank count or threads changed')
    if not float(current['fixture']['parameters']['resume'])>=0:
        raise ValueError('An evolved fixture is required')
    if float(time_end)>float(certificate['validation_time_end']):
        raise ValueError('Timing window exceeds validated interval')


def ordinary(record):
    return [s['seconds'] for s in record['analysis']['step_records']
            if not s['checkpoint'] and not s['remap']]


def verify_instrumented(repo, root, entries):
    identity=json.loads((root/'experimental-identity.json').read_text())
    if identity['revision']!=LEGACY_REF or identity['authoritative'] is not False:
        raise ValueError('Unexpected experimental legacy identity')
    changed=identity['instrumented_files']
    if not set(changed).issubset({name for _,_,name in entries}):
        raise ValueError('Unknown tracked instrumentation file')
    verify(root,[entry for entry in entries if entry[2] not in changed])
    for name,expected in changed.items():
        if digest(root/name)!=expected: raise ValueError('Instrumented source changed: '+name)
    additions={}
    for name,source in [('parallel_block_profile.f90',repo/'src/parallel_block_profile.f90'),
                        ('legacy_profile_io.f90',Path(__file__).parent/'legacy_profile_io.f90')]:
        actual=digest(root/'src'/name)
        if actual!=digest(source): raise ValueError('Instrumentation module differs: '+name)
        additions[name]=actual
    return {**identity,'additional_module_sha256':additions}


def run(args):
    repo=args.repo.resolve(strict=True);fixture=args.fixture.resolve(strict=True)
    out=args.out.resolve();legacy=args.legacy_source.resolve(strict=True)
    instrumented=args.instrumented_source.resolve(strict=True);block=args.block_source.resolve(strict=True)
    if any(out.is_relative_to(p) or p.is_relative_to(out) for p in (repo,fixture,legacy,instrumented,block)):
        raise ValueError('Use a new external output directory, outside all sources and inputs')
    entries=tracked_entries(repo,LEGACY_REF)
    verify(legacy,entries)
    instrumented_identity=verify_instrumented(repo,instrumented,entries)
    block_identity=json.loads((block/'build-identity.json').read_text())
    verify(block,tracked_entries(repo,block_identity['revision']))
    if digest(block/'bin/climate')!=block_identity['binary_sha256']:
        raise ValueError('Block binary differs from its clean build identity')
    binaries={'legacy':legacy/'bin/climate','block':block/'bin/climate',
              'instrumented':instrumented/'bin/climate'}
    initial=binding(binaries,fixture,args.ranks)
    out.mkdir(parents=True,exist_ok=False)
    result={'host':platform.uname()._asdict(),'placement_verified':False,
            'legacy_ref':LEGACY_REF,'block_build':block_identity,
            'legacy_instrumentation':instrumented_identity,
            'binding':initial,'validation':{},'timing':[],
            'scope':'Local J6 diagnosis; not a cluster speedup forecast'}
    result['helper_sha256']={p.name:digest(p) for p in Path(__file__).parent.glob('*.py')}
    save(out/'results.json',result)
    seed=checkpoint_name((fixture/'simple.in').read_text())
    schedule=[('legacy','legacy',False,False),('block','block',False,False),
              ('block-oracle','block',True,False),('legacy-off','instrumented',False,False),
              ('legacy-detail','instrumented',False,True),('block-detail','block',False,True)]
    for label,key,oracle,detail in schedule:
        if binding(binaries,fixture,args.ranks)!=initial: raise ValueError('Experiment inputs changed')
        record=execute(binaries[key],fixture,out/label,args.ranks,args.validation_end,
                       oracles=oracle,detail=detail,monitor=True)
        if binding(binaries,fixture,args.ranks)!=initial: raise ValueError('Experiment inputs changed during run')
        try:
            record['gate']=compare_runs(out/'legacy',out/label,seed,exact=key in ('legacy','instrumented'))
        except (OSError,ValueError) as error:
            record['gate']={'passed':False,'failures':[str(error)],'policy':POLICY,'checkpoints':{}}
        result['validation'][label]=record
        save(out/'results.json',result)
        if not record['gate']['passed']:
            save(out/'validation.json',{'passed':False,'policy':POLICY,'binding':initial,
                                        'failed_run':label,'failures':record['gate']['failures']})
            raise ValueError('Numerical validation failed; no timing pairs launched. See '+str(out/'results.json'))
    # Direct oracle-on/off comparison additionally guards oracle-dependent repair.
    result['oracle_transparency']=compare_runs(out/'block',out/'block-oracle',seed)
    if not result['oracle_transparency']['passed']:
        save(out/'results.json',result);raise ValueError('Oracle changes production results')
    verify(legacy,entries)
    verify_instrumented(repo,instrumented,entries)
    verify(block,tracked_entries(repo,block_identity['revision']))
    certificate={'passed':True,'policy':POLICY,'binding':initial,
                 'validation_time_end':args.validation_end,
                 'evidence':{label:{'gate':r['gate'],
                             'log_sha256':digest(Path(r['directory'])/'run.log')}
                             for label,r in result['validation'].items()}}
    save(out/'validation.json',certificate)
    pressure={label:result['validation'][label]['memory']['reasons']
              for label in ('legacy','block')
              if not result['validation'][label]['memory']['usable_for_local_timing']}
    if pressure:
        result['timing_blocked_reasons']=pressure
        save(out/'results.json',result)
        print('Numerical gates passed, but memory preflight rejected timing:',pressure,flush=True)
        return
    for pair in range(args.pairs):
        for label in (('legacy','block') if pair%2==0 else ('block','legacy')):
            require_certificate(certificate,binding(binaries,fixture,args.ranks),args.time_end)
            record=execute(binaries[label],fixture,out/f'timing-{pair}-{label}',args.ranks,args.time_end,monitor=True)
            record.update(pair=pair,label=label)
            if not ordinary(record): raise ValueError('No ordinary timing steps')
            result['timing'].append(record);save(out/'results.json',result)
            if not record['memory']['usable_for_local_timing']:
                result['timing_blocked_reasons']={record['directory']:record['memory']['reasons']}
                save(out/'results.json',result)
                print('Memory quality failed; remaining timing pairs cancelled',flush=True)
                return
    ratios=[]
    for pair in range(args.pairs):
        rows={r['label']:r for r in result['timing'] if r['pair']==pair}
        usable=all(r['memory']['usable_for_local_timing'] for r in rows.values())
        states_equal=rows['legacy']['analysis']['printed_states']==rows['block']['analysis']['printed_states']
        ratios.append({'pair':pair,'usable':usable and states_equal,
                       'printed_states_equal':states_equal,
                       'ordinary_mean_ratio':statistics.mean(ordinary(rows['block']))/
                                             statistics.mean(ordinary(rows['legacy'])) if usable and states_equal else None})
    result['paired_ratios']=ratios
    verify(legacy,entries)
    verify_instrumented(repo,instrumented,entries)
    if binding(binaries,fixture,args.ranks)!=initial: raise ValueError('Experiment inputs changed')
    save(out/'results.json',result)
    print(json.dumps(ratios,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','legacy-source','instrumented-source','block-source','fixture','out'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    p.add_argument('--pairs',type=int,default=3)
    p.add_argument('--validation-end',default='0.3350')
    p.add_argument('--time-end',default='0.3190')
    a=p.parse_args()
    if a.ranks<1 or a.pairs<0: p.error('Require positive ranks and nonnegative pairs')
    run(a)
