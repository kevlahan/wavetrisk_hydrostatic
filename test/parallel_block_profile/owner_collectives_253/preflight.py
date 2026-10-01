"""Verify sealed Stage251 evidence and the one-file Stage253 replacement."""
import hashlib
import json
from pathlib import Path
import sys

PRIOR_JOB = '55256'
PRIOR_RESULTS_SHA = '138ec179faea1cea0ed4c8afb1f6b126f5982b5d64e49cf39de8c67014565b51'
PRIOR_PLACEMENT_SHA = 'f5badbd905365bd459266bbeab5cb9215fc3ad22f8baedd0725fdcbf40ec3246'
BASELINE_SHA = 'e3feffb3eaddfbf49c2745cac84807f02eceabecba42eb2967133781e258fa91'
ORIGINAL_SHA = 'eff68d917b06eebea30aa81f12ec7bb1c8163de3c502e27b087025da4a66a8e8'
OPTIMIZED_SHA = 'c26cf0da11547145a6fccef2524d81b3b60563c6b2bbc230c1d19c258205672f'
SOURCE = 'src/parallel_block_program.f90'
GENERATED_PHYSICS = 'source/src/physics/simple_physics/phyparam/'
INPUTS = {'simple.in': 'ed173774b8400630657243fa67b4065faff1faf409ba93c1712bb80174d6ae8d',
          'test_checkpoint_0003.bin.zst': '8c7b651c2d6bc6f481838163591565771e1f9a95fd54ec9adf0ce9b74a117e2d'}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def preflight(kit, prior_kit, candidate, prior_results):
    kit, prior_kit, candidate, prior_results = map(Path, (kit, prior_kit, candidate, prior_results))
    for root, label in ((kit, 'Stage253'), (prior_kit, 'Stage251')):
        for name, expected in json.loads((root/'files-sha256.json').read_text()).items():
            if sha(root/name) != expected:
                raise ValueError(f'{label} kit file differs: {name}')
    if sha(kit/'optimized-parallel_block_program.f90') != OPTIMIZED_SHA:
        raise ValueError('Optimized source differs')
    if sha(prior_results/'results.json') != PRIOR_RESULTS_SHA or sha(prior_results/'placement-83.json') != PRIOR_PLACEMENT_SHA:
        raise ValueError('Stage251 result or placement differs')
    prior = json.loads((prior_results/'results.json').read_text())
    if prior['status'] != 'complete' or prior['slurm_job_id'] != PRIOR_JOB:
        raise ValueError('Wrong Stage251 job')
    manifest = json.loads((candidate.parent/'files-sha256.json').read_text())
    source_key = 'source/'+SOURCE
    if manifest[source_key] != ORIGINAL_SHA or sha(candidate/SOURCE) != ORIGINAL_SHA:
        raise ValueError('Original Stage251 program source differs')
    # Stage251 prepares Makefile.inc and may regenerate the phyparam/tmp
    # Fortran copies. Those are build products; the already-built physics
    # objects are pinned below. All other sealed source remains immutable.
    generated_changes = []
    for name, expected in manifest.items():
        if not name.startswith('source/'):
            continue
        path = candidate/name.removeprefix('source/')
        generated = (name.startswith(GENERATED_PHYSICS+'tmp/') or
                     name == GENERATED_PHYSICS+'Makefile.inc')
        if generated:
            if not path.is_file() or sha(path) != expected:
                generated_changes.append(name)
            continue
        if sha(path) != expected:
            raise ValueError('Stage251 source differs: '+name)
    objects = prior['build']['builds']['main']['physics_objects_sha256']
    directory = candidate/'src/physics/simple_physics/phyparam/obj'
    actual_objects = {path.name: sha(path) for path in sorted(directory.glob('*.o'))}
    if len(objects) != 23 or actual_objects != objects:
        raise ValueError('Stage251 compiled physics objects differ')
    baseline = prior_results/'J5-n83-validation-c5-c5/climate'
    if sha(baseline) != BASELINE_SHA or prior['build']['builds']['c5']['sha256'] != BASELINE_SHA:
        raise ValueError('Stage251 baseline binary differs')
    for name, expected in INPUTS.items():
        if sha(prior_results/'input'/name) != expected:
            raise ValueError('J5 input differs: '+name)
    print('PASS: Stage253 replacement, fixed source, 23 physics objects, binary, input and core map; '
          f'{len(generated_changes)} regenerated physics files', flush=True)
    return generated_changes


if __name__ == '__main__':
    if len(sys.argv) != 5:
        raise SystemExit('Usage: preflight.py KIT STAGE251_KIT CANDIDATE_SOURCE STAGE251_RESULTS')
    preflight(*sys.argv[1:])
