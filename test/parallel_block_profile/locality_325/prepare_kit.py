"""Build a fresh, sealed qualification kit from this checkout's exact sources."""
from pathlib import Path
import argparse
import hashlib
import io
import json
import shutil
import tarfile

NAME='rk325-locality-qualification-kit-v1'


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def prepare(checkout,destination):
    checkout=Path(checkout).resolve(strict=True);destination=Path(destination).resolve()
    template=Path(__file__).resolve().parent
    kit=destination/NAME
    if kit.exists() or (destination/(NAME+'.tar.gz')).exists():raise ValueError('Fresh kit output required')
    destination.mkdir(parents=True,exist_ok=True)
    shutil.copytree(template,kit,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    names=json.loads((kit/'baseline-source-files.json').read_text())
    manifest={}
    with tarfile.open(kit/'sources.tar.gz','w:gz',compresslevel=6) as archive:
        for name in sorted(names):
            rel=Path(name);path=checkout/rel
            if rel.is_absolute() or '..' in rel.parts or not path.resolve().is_relative_to(checkout):raise ValueError('Unsafe source path')
            raw=path.read_bytes();manifest[name]=hashlib.sha256(raw).hexdigest()
            entry=tarfile.TarInfo('candidate/'+name);entry.size=len(raw);entry.mode=0o644;entry.mtime=0
            archive.addfile(entry,io.BytesIO(raw))
    (kit/'source-files.json').write_text(json.dumps({'candidate':manifest},indent=2,sort_keys=True)+'\n')
    old=json.loads((kit/'baseline-source-files.json').read_text())
    changes={n:{'baseline':old[n],'candidate':h} for n,h in manifest.items() if h!=old[n]}
    (kit/'source-provenance.json').write_text(json.dumps({'stage':325,'baseline':json.loads((kit/'baseline-provenance.json').read_text()),
        'changed_sources':changes,'performance_qualified':False,'native_qualification_required':True},indent=2)+'\n')
    (kit/'IDENTITY.json').write_text(json.dumps({'stage':325,'revision':'locality-qualification-v1','builds':2,'solver_launches':3,'performance_launches':0,
        'source_manifest_sha256':digest(kit/'source-files.json')},indent=2)+'\n')
    files={str(p.relative_to(kit)):digest(p) for p in sorted(kit.rglob('*')) if p.is_file() and p.name!='files-sha256.json' and '__pycache__' not in p.parts}
    (kit/'files-sha256.json').write_text(json.dumps(files,indent=2,sort_keys=True)+'\n')
    target=destination/(NAME+'.tar.gz')
    with tarfile.open(target,'w:gz',compresslevel=6) as archive:
        for p in sorted(kit.rglob('*')):
            if p.is_file():archive.add(p,arcname=NAME+'/'+str(p.relative_to(kit)),recursive=False)
    checksum=destination/(NAME+'.tar.gz.sha256')
    checksum.write_text(digest(target)+'  '+target.name+'\n')
    return {'kit':str(kit),'archive':str(target),'sha256':digest(target),'source_files':len(manifest),'changed_sources':len(changes),'files':len(files)}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('checkout');parser.add_argument('destination');args=parser.parse_args()
    print(json.dumps(prepare(args.checkout,args.destination),indent=2))
