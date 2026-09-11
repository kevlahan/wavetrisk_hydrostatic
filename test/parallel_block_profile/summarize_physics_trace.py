"""Compare completed, isolated phase/physics traces and checkpoint transparency."""
import argparse
import json
from pathlib import Path
from compare_physics_probe import compare as physics_compare
from compare_restart_probe import compare as phase_compare
from compare_restart_checkpoints import compare as checkpoint_compare


def summarize(root,reference):
    result={'root':str(root),'reference':str(reference),'transparency':{},'steps':{}}
    for label in ('legacy','block'):
        identity=json.loads((root/label/'identity.json').read_text())
        result.setdefault('identities',{})[label]=identity
        result['transparency'][label]={}
        for cp in (5,6):
            name=f'test_checkpoint_{cp:04d}.bin.zst'
            c=checkpoint_compare(reference/label/name,root/label/name)
            result['transparency'][label][str(cp)]=c
            if (c['header_differences'] or c['topology_differences'] or c['threshold_max_abs'] or
                    c['threshold_nonfinite'] or any(s['different'] or s['nonfinite'] for s in c['fields'].values())):
                raise ValueError(f'Probe altered {label} checkpoint {cp}')
    steps=sorted({int(p.name.split('-step-')[1].split('-rank-')[0])
                  for p in (root/'legacy').glob('physics-input-step-*.bin')})
    for step in steps:
        entry={}
        for tag in ('input','output'):
            pattern=f'physics-{tag}-step-{step}-rank-*.bin'
            entry['physics-'+tag]=physics_compare(str(root/'legacy'/pattern),str(root/'block'/pattern))
        for tag in ('after-dynamics','after-physics','after-adaptation','step-end'):
            pattern=f'probe-{tag}-step-{step}-rank-*.bin'
            entry[tag]=phase_compare(str(root/'legacy'/pattern),str(root/'block'/pattern))
        result['steps'][str(step)]=entry
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    result=summarize(args.root,args.reference)
    with args.out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
