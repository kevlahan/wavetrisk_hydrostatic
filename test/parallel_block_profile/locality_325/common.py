"""Stage325 immutable-package, atomic metadata and private-copy helpers."""
from pathlib import Path
import hashlib,json,os,shutil
PLACEMENT_RANKS=(256,)

def sha(path):
 p=Path(path)
 if p.is_symlink():return hashlib.sha256(str(p.readlink()).encode()).hexdigest()
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def save(path,value):
 p=Path(path);tmp=p.with_suffix(p.suffix+'.tmp')
 with tmp.open('w') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 tmp.replace(p)

def verify(root,filename='files-sha256.json'):
 root=Path(root).resolve(strict=True);manifest=json.loads((root/filename).read_text())
 if not isinstance(manifest,dict):raise ValueError('Invalid sealed kit manifest')
 actual={str(p.relative_to(root)) for p in root.rglob('*') if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.relative_to(root).parts and str(p.relative_to(root))!=filename}
 if actual!=set(manifest):raise ValueError('Sealed kit inventory differs')
 for name,want in manifest.items():
  p=Path(name)
  if p.is_absolute() or '..' in p.parts or sha(root/p)!=want:raise ValueError('Sealed file differs: '+name)
 return manifest

def private_copy_checkpoint(source,destination):
 source,destination=Path(source),Path(destination)
 if source.is_symlink() or not source.is_file():raise ValueError('Checkpoint source must be regular')
 if destination.exists() or destination.is_symlink():raise ValueError('Checkpoint destination already exists')
 before=sha(source);ss=source.stat()
 with source.open('rb') as src,destination.open('xb') as dst:shutil.copyfileobj(src,dst,1024*1024);dst.flush();os.fsync(dst.fileno())
 ds=destination.stat()
 if ds.st_nlink!=1 or (ss.st_dev,ss.st_ino)==(ds.st_dev,ds.st_ino):raise ValueError('Checkpoint is not a private copy')
 if sha(source)!=before or sha(destination)!=before:raise ValueError('Checkpoint changed during copy')
 return dict(source=str(source),destination=str(destination),sha256=before,private_copy=True,destination_links=ds.st_nlink)

def preflight(kit):
 from check_kit import check
 return check(Path(kit))
