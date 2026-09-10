#!/usr/bin/env python3
"""Run an isolated allocation census and require checkpoint transparency.

No timing comparisons. Rank peaks are kept separate; aligned samples are not
claimed to be simultaneous wall-clock observations because no barrier is added.
"""
import argparse
import json
from pathlib import Path
from checkpoint_gate import compare_runs
from experiment import digest
from analyze import analyze
from profile_protocol import execute,fixture_identity,checkpoint_name


def read_census(path,expected):
    samples=[];current=None
    for line in path.read_text().splitlines():
        words=line.split()
        if words[0]=='sample':
            current={'sequence':int(words[1]),'phase':words[2],'allocations':{}}
            samples.append(current)
        else:
            if current is None or len(words)!=2:raise ValueError('Malformed allocation record')
            value=int(words[0]);key=words[1]
            if key in current['allocations'] or value<0:raise ValueError('Duplicate/negative allocation')
            current['allocations'][key]=value
    if not samples:raise ValueError('No allocation samples')
    for i,sample in enumerate(samples,1):
        if sample['sequence']!=i or set(sample['allocations'])!=expected:
            raise ValueError('Incomplete census sample or unknown allocation path')
        sample['total_bytes']=sum(sample['allocations'].values())
    return samples


def summarize(directory,manifest,ranks):
    expected={key for coverage in manifest['coverage'].values() for key in coverage['allocations']}
    records=[];rank_samples=[]
    for rank in range(ranks):
        samples=read_census(directory/f'allocation-rank-{rank}.txt',expected);rank_samples.append(samples)
        peak=max(samples,key=lambda s:s['total_bytes'])
        roots={}
        for key,value in peak['allocations'].items():
            root='.'.join(key.split('.')[:2]);roots[root]=roots.get(root,0)+value
        records.append({'rank':rank,'samples':len(samples),'peak_sample':peak['sequence'],'peak_phase':peak['phase'],
                        'peak_counted_bytes':peak['total_bytes'],
                        'owners_at_peak':{owner:sum(v for k,v in peak['allocations'].items() if k.startswith(owner+'.'))
                                          for owner in manifest['coverage']},
                        'largest_roots_at_peak':dict(sorted(roots.items(),key=lambda kv:kv[1],reverse=True)[:15])})
    aligned=all([(s['sequence'],s['phase']) for s in r]==
                [(s['sequence'],s['phase']) for s in rank_samples[0]] for r in rank_samples)
    result={'scope':manifest['scope'],'excluded':manifest['excluded'],'ranks':records,'phase_sequences_align':aligned}
    if aligned:
        sums=[sum(r[i]['total_bytes'] for r in rank_samples) for i in range(len(rank_samples[0]))]
        peak=max(range(len(sums)),key=sums.__getitem__)
        roots={}
        for r in rank_samples:
            for key,value in r[peak]['allocations'].items():
                root='.'.join(key.split('.')[:2]);roots[root]=roots.get(root,0)+value
        result['largest_aligned_phase_sum']={'sequence':peak+1,'phase':rank_samples[0][peak]['phase'],
                                            'bytes':sums[peak],
                                            'caution':'Same logical sample, not a synchronized instantaneous global peak',
                                            'roots':dict(sorted(roots.items(),key=lambda kv:kv[1],reverse=True))}
    return result


def run(args):
    build=args.build.resolve(strict=True);reference=args.reference.resolve(strict=True)
    fixture=args.fixture.resolve(strict=True);out=args.out.resolve()
    manifest=json.loads((build/'allocation-census-identity.json').read_text())
    for path,expected in manifest['changed_files'].items():
        if digest(build/path)!=expected:raise ValueError('Census source changed: '+path)
    identity=json.loads((reference/'identity.json').read_text())
    if identity['sha256']!=manifest['baseline_binary_sha256']:raise ValueError('Wrong numerical reference binary')
    if identity['ranks']!=args.ranks or float(identity['time_end'])!=float(args.time_end):
        raise ValueError('Reference rank count/window differ')
    for name,value in identity['environment'].items():
        if name.startswith('WAVETRISK_') and value!='0':raise ValueError('Reference must be production, oracles/profiling off')
    if fixture_identity(reference)!=fixture_identity(fixture):raise ValueError('Reference input/seed/grid assets differ')
    if not analyze(reference/'run.log')['completed']:raise ValueError('Reference run did not complete')
    initial=fixture_identity(fixture)
    record=execute(build/'bin/climate',fixture,out,args.ranks,args.time_end,census=True,monitor=True)
    gate=compare_runs(reference,out,checkpoint_name((fixture/'simple.in').read_text()),exact=True)
    result={'identity':record['identity'],'numerical_transparency':gate,
            'census':summarize(out,manifest,args.ranks),'memory':record['memory']}
    (out/'census-results.json').write_text(json.dumps(result,indent=2)+'\n')
    if fixture_identity(fixture)!=initial:raise ValueError('Fixture changed')
    print('Checkpoint transparency:',gate['passed'],gate['failures'],flush=True)
    print(json.dumps(result['census'],indent=2),flush=True)
    if not gate['passed']:raise ValueError('Census changes numerical results; do not accept it')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','reference','fixture','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4);p.add_argument('--time-end',default='0.3350')
    run(p.parse_args())
