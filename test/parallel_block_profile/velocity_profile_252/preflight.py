"""Pin the verified Stage251 result and binaries before bb13 submission."""
import hashlib
import json
from pathlib import Path
import sys


PRIOR_JOB = '55256'
PRIOR_RESULTS_SHA = '138ec179faea1cea0ed4c8afb1f6b126f5982b5d64e49cf39de8c67014565b51'
PRIOR_PLACEMENT_SHA = 'f5badbd905365bd459266bbeab5cb9215fc3ad22f8baedd0725fdcbf40ec3246'
BINARIES = {'c4': '02ce866f8c82763eb5b600ea5563e7ce484a2e7f88fb63a1f3b8f161778712ec',
            'c5': 'e3feffb3eaddfbf49c2745cac84807f02eceabecba42eb2967133781e258fa91'}
INPUTS = {'simple.in': 'ed173774b8400630657243fa67b4065faff1faf409ba93c1712bb80174d6ae8d',
          'test_checkpoint_0003.bin.zst': '8c7b651c2d6bc6f481838163591565771e1f9a95fd54ec9adf0ce9b74a117e2d'}


def sha(path):
    with Path(path).open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def preflight(perf, prior_kit, candidate, prior_results):
    perf, prior_kit, candidate, prior_results = map(Path, (perf, prior_kit, candidate, prior_results))
    for name, expected in json.loads((perf/'files-sha256.json').read_text()).items():
        if sha(perf/name) != expected:
            raise ValueError('Profile kit file differs: '+name)
    for name, expected in json.loads((prior_kit/'files-sha256.json').read_text()).items():
        if sha(prior_kit/name) != expected:
            raise ValueError('Stage251 kit file differs: '+name)
    if sha(prior_results/'results.json') != PRIOR_RESULTS_SHA:
        raise ValueError('Stage251 result differs from the verified review archive')
    if sha(prior_results/'placement-83.json') != PRIOR_PLACEMENT_SHA:
        raise ValueError('Stage251 core map differs')
    prior = json.loads((prior_results/'results.json').read_text())
    if prior['status'] != 'complete' or prior['slurm_job_id'] != PRIOR_JOB:
        raise ValueError('Wrong Stage251 job')
    if not candidate.is_dir():
        raise ValueError('Stage251 source tree is absent')
    manifest = json.loads((candidate.parent/'files-sha256.json').read_text())
    launcher = 'source/test/parallel_block_profile/rank_launch_223.py'
    if sha(candidate.parent/launcher) != manifest[launcher]:
        raise ValueError('Rank launcher differs from sealed candidate source')
    for kind, expected in BINARIES.items():
        binary = prior_results/f'J5-n83-validation-{kind}-{kind}/climate'
        if sha(binary) != expected or prior['build']['builds'][kind]['sha256'] != expected:
            raise ValueError('Stage251 binary differs: '+kind)
    for name, expected in INPUTS.items():
        if sha(prior_results/'input'/name) != expected:
            raise ValueError('J5 input differs: '+name)
    print('PASS: Stage252 kit, Stage251 result, binaries, inputs and core map', flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 5:
        raise SystemExit('Usage: preflight.py PROFILE_KIT STAGE251_KIT CANDIDATE_SOURCE STAGE251_RESULTS')
    preflight(*sys.argv[1:])
