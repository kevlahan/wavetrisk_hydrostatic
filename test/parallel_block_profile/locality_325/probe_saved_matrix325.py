"""Offline preflight only: exact per-rank gates on a bounded saved-capture subset."""
import argparse,ast,json,tempfile,time
from pathlib import Path
import verify325 as q

EXPECTED={'checked-L256':[1],'qual-C0':[1,2,69,255],'qual-L256':[1,52,160,255]}

def capture_inventory(base):
    base=Path(base);result=[]
    for arm,ranks in EXPECTED.items():
        expected=[f'execution316-rank-{rank}.bin.gz' for rank in ranks]
        found=sorted(p.name for p in (base/arm).glob('execution316-rank-*'))
        q.require(found==sorted(expected),'probe.exact_capture_inventory',{'arm':arm,'expected':expected,'found':found})
        result.extend((arm,rank,base/arm/name) for rank,name in zip(ranks,expected))
    return result

def common_feature_keys(a,b):
    common=a.keys()&b.keys()
    q.require(bool(common),'probe.empty_common_features')
    return common

def resolve_probe_paths(capture_root,output):
    base=Path(capture_root).resolve();output=Path(output).resolve()
    q.require(output!=base and base not in output.parents,'probe.output_in_capture_root')
    return base,output

def main(capture_root,output):
    base,output=resolve_probe_paths(capture_root,output)
    started=time.monotonic();report={'passed':False,'scope':'nine selected saved native captures; not world qualification','ranks':[]}
    captures=capture_inventory(base)
    with tempfile.TemporaryDirectory(prefix='rk325-saved-matrix-') as scratch:
        c=q._create_cache(Path(scratch)/'routes.sqlite')
        for arm,rank,path in captures:
            print('Exact per-rank gates',arm,rank,flush=True)
            summary=q.cache_rank325(c,path,arm,rank)
            report['ranks'].append({'arm':arm,'rank':rank,'reader':summary})
        a=dict(c.execute('select key,digest from feature where arm=?',('qual-C0',)))
        b=dict(c.execute('select key,digest from feature where arm=?',('qual-L256',)))
        common=common_feature_keys(a,b);bad=[k for k in common if a[k]!=b[k]]
        roles={}
        for key in common:
            role=ast.literal_eval(key)[0][0];roles[role]=roles.get(role,0)+1
        report['common_features']={'count':len(common),'roles':roles,'differences':bad[:30],
                                  'control_only':len(a.keys()-b.keys()),'candidate_only':len(b.keys()-a.keys())}
        report['database_bytes']=sum(p.stat().st_size for p in Path(scratch).iterdir() if p.is_file())
        report['elapsed_seconds']=time.monotonic()-started
        c.close()
        report['passed']=not bad
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='ranks'},indent=2))
    q.require(not bad,'probe.common_feature_difference',bad[:3])

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('capture_root',type=Path,help='directory containing the exact checked/control/locality nine-capture subset')
    parser.add_argument('output_report',type=Path,help='JSON report destination outside the capture root')
    args=parser.parse_args()
    main(args.capture_root,args.output_report)
