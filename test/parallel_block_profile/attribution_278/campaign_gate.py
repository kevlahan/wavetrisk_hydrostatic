"""Validate the connected RHS campaign entry gate; never submits a Slurm job.

Fail closed until the active candidate is qualified. The evidence must bind full local
production qualification to exact source bytes before a kit may be sealed.
"""
import argparse, hashlib, json
from pathlib import Path

REQUIRED = ('production_active', 'ownership_credit', 'full_step_science_qualified',
            'connected_retirement_qualified', 'replacement_cost_case_qualified')
BASE = 'ff740744608731a05d53771c32a3055309f3a586'

def check(root, evaluation):
    errors=[]
    if evaluation.get('accepted_source') != BASE:
        errors.append('accepted source is not Stage266')
    for key in REQUIRED:
        if evaluation.get(key) is not True:
            errors.append(key+' is not qualified')
    manifest=evaluation.get('qualified_source_sha256',{})
    if not manifest:
        errors.append('source-bound production qualification manifest is missing')
    for relative,wanted in manifest.items():
        target=(root/relative).resolve()
        if not target.is_relative_to(root.resolve()) or not target.is_file():
            errors.append('invalid qualified source path: '+relative)
        elif hashlib.sha256(target.read_bytes()).hexdigest()!=wanted:
            errors.append('qualified source changed: '+relative)
    for mode in ('optimized','checked','oracle'):
        record=evaluation.get('full_step_local_runs',{}).get(mode,{})
        if record.get('passed') is not True or record.get('ranks',0) not in range(1,5):
            errors.append('complete local '+mode+' replay is missing')
    if not evaluation.get('retired_production_work'):
        errors.append('no retired production work is recorded')
    if evaluation.get('validation_only_reference_packets_in_production') is not False:
        errors.append('absence of ordinary-stage reference packets is not established')
    return {'status':'blocked' if errors else 'entry_gate_passed', 'launch_authorized':False,
            'sealing_ready':not errors, 'errors':errors,
            'scope':'Local production entry gate only; cluster performance/promotion remain separate'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--evaluation',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=check(a.root,json.loads(a.evaluation.read_text()))
    a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
    raise SystemExit(0 if r['sealing_ready'] else 2)
