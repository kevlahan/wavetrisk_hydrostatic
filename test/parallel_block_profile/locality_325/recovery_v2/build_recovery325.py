#!/usr/bin/env python3
"""Reproduce sealed recovery payload from unique files and shared Stage325 helpers.

No solver, compilation, scheduler, or remote operation is available here.
The tar/gzip metadata is normalized; its checksum differs from the deployed
archive, while every packaged file and its sealed manifest remain identical.
"""
from pathlib import Path
import argparse,gzip,hashlib,io,json,subprocess,sys,tarfile

def digest(data):
    return hashlib.sha256(data).hexdigest()

def require(ok,message):
    if not ok:
        raise ValueError(message)

def build(recovery,shared,output,run_tests=False):
    recovery=Path(recovery).resolve(strict=True)
    shared=Path(shared).resolve(strict=True)
    output=Path(output).resolve()
    require(not output.exists(),'Fresh output directory required')
    recipe=json.loads((recovery/'package_recipe325.json').read_text())
    manifest_bytes=(recovery/'files-sha256.json').read_bytes()
    require(digest(manifest_bytes)==recipe['sealed_manifest_sha256'],'Sealed manifest changed')
    manifest=json.loads(manifest_bytes)
    require(set(manifest)==set(recipe['sources']),'Recipe and sealed inventory differ')
    payload={}
    for name,row in recipe['sources'].items():
        rel=Path(name)
        require(not rel.is_absolute() and '..' not in rel.parts and str(rel)==name,'Unsafe package path')
        kind=row['source']
        require(kind in ('recovery','shared','shared-normalize-final-newlines'),'Unknown source kind')
        source=(recovery if kind=='recovery' else shared)/rel
        require(source.is_file() and not source.is_symlink(),'Expected regular source file: '+name)
        data=source.read_bytes()
        if kind=='shared-normalize-final-newlines':
            require(name=='verify322.py','Unexpected newline-normalized helper')
            data=data.rstrip(b'\n')+b'\n\n'
        require(digest(data)==row['sha256']==manifest[name],'Source differs from sealed payload: '+name)
        payload[name]=data
    payload['files-sha256.json']=manifest_bytes
    output.mkdir()
    kit=output/recipe['artifact'];kit.mkdir()
    for name,data in payload.items():
        dest=kit/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
        dest.chmod(0o755 if name.endswith('.sh') else 0o644)
    subprocess.run([sys.executable,'-c','from pathlib import Path; from recover325 import package_check; print(package_check(Path.cwd()))'],cwd=kit,check=True)
    if run_tests:
        with (output/'reproduction-tests.log').open('w') as log:
            subprocess.run([sys.executable,'-m','unittest','discover','-v'],cwd=kit,stdout=log,stderr=subprocess.STDOUT,check=True)
    archive=output/(recipe['artifact']+'.tar.gz')
    with archive.open('wb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0,compresslevel=9) as gz,tarfile.open(fileobj=gz,mode='w',format=tarfile.PAX_FORMAT) as tar:
        for name,data in sorted(payload.items()):
            info=tarfile.TarInfo(recipe['artifact']+'/'+name)
            info.size=len(data);info.mode=0o755 if name.endswith('.sh') else 0o644
            info.mtime=0;info.uid=info.gid=0;info.uname=info.gname=''
            tar.addfile(info,io.BytesIO(data))
    h=digest(archive.read_bytes())
    Path(str(archive)+'.sha256').write_text(h+'  '+archive.name+'\n')
    result={'archive':str(archive),'sha256':h,'sealed_manifest_sha256':recipe['sealed_manifest_sha256'],'deployed_archive_sha256':recipe['deployed_archive_sha256'],'packaged_files':len(payload),'payload_identical_to_deployed':True,'normalized_archive_metadata':True,'tests_requested':run_tests}
    (output/'reproduction.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recovery',type=Path,default=here)
    parser.add_argument('--shared',type=Path,default=here.parent)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--run-tests',action='store_true')
    args=parser.parse_args()
    print(json.dumps(build(args.recovery,args.shared,args.output,args.run_tests),indent=2))
