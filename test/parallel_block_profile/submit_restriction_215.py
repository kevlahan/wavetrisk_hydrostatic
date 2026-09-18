#!/usr/bin/env python3
"""Verify kit, local validation, sources and executable identities before running."""
import json,subprocess,sys
from pathlib import Path
from build_restriction_208 import digest,verify,verify_helpers
from build_restriction_215 import KINDS

def submit(kit,tree,fixture,out):
    kit=kit.resolve(strict=True);tree=tree.resolve(strict=True);verify_helpers(kit)
    m=json.loads((kit/'stage215-sources.json').read_text());local=json.loads((kit/'LOCAL_VALIDATION.json').read_text())
    if local.get('stage')!=215 or local.get('status')!='passed' or local.get('candidate_sources')!=m['candidate_sources']:
        raise ValueError('Exact candidate local validation required')
    if set(local.get('runs',{}))!={'optimized','checked-off','oracles-on'}:raise ValueError('Incomplete local numerical gates')
    path=tree/'restriction-215-builds.json';r=json.loads(path.read_text())
    if r.get('status')!='passed' or r.get('fp_contract')!='off' or r.get('manifest_sha256')!=digest(kit/'stage215-sources.json'):
        raise ValueError('Build identity differs')
    verify(tree,m['candidate_sources']);command=[sys.executable,str(kit/'cluster_215.py')]
    for kind in KINDS:
        record=r['builds'][kind];binary=tree/record['binary']
        if record['sources']!=m['reference_sources' if kind=='baseline' else 'candidate_sources'] or digest(binary)!=record['sha256']:
            raise ValueError('Source/binary mismatch '+kind)
        command+=['--'+kind,str(binary)]
    command+=['--fixture',str(fixture.resolve(strict=True)),'--out',str(out.resolve()),'--build-notes-file',str(path)]
    subprocess.run(command,check=True)
if __name__=='__main__':submit(*(Path(p) for p in sys.argv[1:]))
