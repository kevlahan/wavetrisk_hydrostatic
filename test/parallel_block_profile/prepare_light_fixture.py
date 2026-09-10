"""Derive an evolved J6 seed using unchanged legacy's supported capped restart."""
import argparse
import json
from pathlib import Path
import shutil
from build_legacy import tracked_entries,verify
from experiment import LEGACY_REF,digest
from profile_protocol import execute,replace_parameter,fixture_identity,parameters
from checkpoint_gate import topology


def prepare(repo,legacy_source,fixture,out,ranks):
    repo=repo.resolve(strict=True);legacy_source=legacy_source.resolve(strict=True);out=out.resolve()
    if out.is_relative_to(repo) or out.is_relative_to(fixture.resolve()) or out.is_relative_to(legacy_source):
        raise ValueError('Use a new external output directory')
    entries=tracked_entries(repo,LEGACY_REF);verify(legacy_source,entries)
    original=fixture_identity(fixture)
    p=original['parameters']
    if (p.get('zlevels'),p.get('Nsoil'),p.get('resume'),p.get('max_level'))!=('30','10','3','7'):
        raise ValueError('Expected the original J7/30+soil checkpoint-3 fixture')
    out.mkdir(parents=True,exist_ok=False)
    binary=legacy_source/'bin/climate'
    record=execute(binary,fixture,out/'generation',ranks,'0.3105',overrides={'max_level':6})
    checkpoint=out/'generation'/f"{p['run_id']}_checkpoint_0004.bin.zst"
    shape=topology(checkpoint)
    if shape['maximum_level']!=6: raise ValueError('Expected an evolved checkpoint with level 6 but no level 7')
    seed=out/'seed';seed.mkdir()
    shutil.copy2(checkpoint,seed/checkpoint.name)
    (seed/'grids').symlink_to((fixture/'grids').resolve(strict=True),target_is_directory=True)
    text=(out/'generation/simple.in').read_text()
    text=replace_parameter(text,'resume',4);text=replace_parameter(text,'time_end','0.3190')
    (seed/'simple.in').write_text(text)
    verify(legacy_source,entries)
    if fixture_identity(fixture)!=original: raise ValueError('Original fixture changed')
    manifest={'role':'derived-local-benchmark','legacy_ref':LEGACY_REF,'legacy_binary_sha256':digest(binary),
              'source_fixture':original,'fixture':fixture_identity(seed),'topology':shape,
              'generation':record,'validation_time_end':'0.3350','timing_time_end':'0.3190'}
    (out/'fixture.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'seed':str(seed),'topology':shape},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','legacy-source','fixture','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    a=p.parse_args();prepare(a.repo,a.legacy_source,a.fixture,a.out,a.ranks)
