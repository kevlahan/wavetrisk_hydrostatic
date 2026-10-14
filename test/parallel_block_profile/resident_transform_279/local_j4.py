#!/usr/bin/env python3
"""Run the explicit J4/J6 local fixture and measure refinement and VM pressure.

Existing J5 acceptance policy/defaults remain unchanged. Above-threshold wavelet
locations measure refinement demand; they are not a reconstructed node mask.
"""
import argparse
from collections import Counter,deque
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import tarfile
import numpy as np
from analyze import analyze
from compare_restart_checkpoints import read,compare
from checkpoint_gate import accept
from experiment import digest
from profile_protocol import MemoryMonitor,environment,parameters,grid_identity,checkpoint_name


def refinement(path):
    patches=Counter();nodes=Counter();edges=Counter();nonzero=Counter()
    for header,threshold,poles,coarse,waves,required in read(path,-10,30,expected_domains=40):
        if not all(np.isfinite(v).all() for v in (threshold,poles,coarse,waves)):
            raise ValueError('Nonfinite J4 checkpoint')
        levels=deque([4]*4)
        limits=threshold.reshape(41,3)[11:,:]
        for wave,flags in zip(waves,required):
            if not levels:raise ValueError('Trailing wavelet records')
            level=levels.popleft()
            if level>6 or any(f not in (0,1) for f in flags):raise ValueError('Invalid J4/J6 topology')
            patches[level]+=1
            w=np.abs(wave[11:,:]);nonzero[level]+=int(np.count_nonzero(w))
            edges[level]+=int(np.count_nonzero(np.any(w[:,:48]>limits[:,0,None],axis=0)))
            active=np.any(w[:,48:64]>limits[:,1,None],axis=0)|np.any(w[:,64:80]>limits[:,2,None],axis=0)
            nodes[level]+=int(np.count_nonzero(active))
            levels.extend([level+1]*int(np.count_nonzero(flags)))
        if levels:raise ValueError('Missing wavelet records')
    total=sum(nodes.values())+sum(edges.values());finest=nodes[6]+edges[6]
    return {'checkpoint':str(path),'sha256':digest(path),'patches_by_level':dict(patches),
            'above_threshold_node_wavelet_locations':dict(nodes),'above_threshold_edge_wavelet_locations':dict(edges),
            'nonzero_atmospheric_coefficients':dict(nonzero),'finest_above_threshold_locations':finest,
            'finest_fraction_of_above_threshold_locations':finest/total if total else 0.,
            'scope':'Atmospheric wavelet locations exceeding stored adaptation thresholds in any physical layer; excludes extra mask/adjacency nodes'}


def vtk_refinement(path,layer=1):
    """Count the actually exported adaptive triangles and their vertices."""
    with tarfile.open(path,'r:gz') as archive:
        members=[m for m in archive.getmembers() if m.isfile() and re.search(rf'_tri_{layer:03d}_\d{{4}}\.vtk$',m.name)]
        if len(members)!=1:raise ValueError('Expected one atmospheric VTK layer')
        raw=archive.extractfile(members[0]).read()
    pos=0
    def line():
        nonlocal pos
        end=raw.index(b'\n',pos);value=raw[pos:end];pos=end+1;return value
    if [line() for _ in range(4)]!=[b'# vtk DataFile Version 2.0',b'WAVETRISK adaptive data',b'BINARY',b'DATASET POLYDATA']:
        raise ValueError('Unsupported VTK format')
    parts=line().split()
    if parts[0]!=b'POINTS' or parts[2]!=b'float':raise ValueError('Unsupported VTK points')
    points=int(parts[1]);pos+=12*points
    parts=line().split()
    if parts[0]!=b'POLYGONS':raise ValueError('Missing VTK triangles')
    cells,nints=map(int,parts[1:])
    if nints!=4*cells:raise ValueError('Expected triangles')
    conn=np.frombuffer(raw,dtype='>i4',count=nints,offset=pos).reshape(cells,4);pos+=4*nints
    if np.any(conn[:,0]!=3) or np.any(conn[:,1:]<0) or np.any(conn[:,1:]>=points):raise ValueError('Invalid connectivity')
    if line().split()!=[b'CELL_DATA',str(cells).encode()] or line()!=b'SCALARS Level int' or line()!=b'LOOKUP_TABLE default':
        raise ValueError('Missing VTK level data')
    levels=np.frombuffer(raw,dtype='>i4',count=cells,offset=pos)
    if np.any((levels<4)|(levels>6)):raise ValueError('VTK levels outside J4/J6')
    counts={int(l):int(np.count_nonzero(levels==l)) for l in np.unique(levels)}
    finest_vertices=len(np.unique(conn[levels==6,1:]))
    return {'archive':str(path),'layer':layer,'active_triangles_by_level':counts,
            'vertices':points,'vertices_on_level6_triangles':finest_vertices,
            'level6_triangle_fraction':counts.get(6,0)/cells if cells else 0.,
            'scope':'Exported adaptive mesh; unique vertices touching level-6 triangles, not raw checkpoint patch capacity'}


def compare_windows(left,right,seed_name,*,expected_new):
    """Exact J4 comparison with the explicitly selected single-restart window."""
    reports={};failures=[]
    identities=[json.loads((root/'results.json').read_text()) for root in (left,right)]
    for key in ('parameters','input_sha256','seed_sha256','grids','ranks'):
        if identities[0][key]!=identities[1][key]:raise ValueError('Unmatched J4 inputs: '+key)
    run_id=identities[0]['parameters']['run_id']
    expected={f'{run_id}_checkpoint_{index:04d}.bin.zst' for index in expected_new}
    for root in (left,right):
        found={p.name for p in root.glob('*_checkpoint_*.bin.zst') if p.name!=seed_name}
        if found!=expected:raise ValueError('Unexpected J4 checkpoint set')
        log=(root/'run.log').read_text()
        if not analyze(root/'run.log')['completed']:failures.append(str(root)+': run incomplete')
        for index in expected_new:
            matches=list(re.finditer(r'Restarting from checkpoint\s+'+str(index)+r'\b',log))
            if len(matches)!=1:failures.append(str(root)+': missing/duplicate reload '+str(index));continue
            if 'Remapping vertical coordinates' not in log[matches[0].end():]:
                failures.append(str(root)+': no remap following reload '+str(index))
    for name in sorted(expected):
        report=compare(left/name,right/name,expected_domains=40);reports[name]=report
        failures.extend(name+': '+reason for reason in accept(report,exact=True))
    return {'passed':not failures,'failures':failures,'exact_fields_required':True,
            'domains':40,'min_level':4,'expected_new_checkpoint_indices':list(expected_new),'checkpoints':reports}


def run(args):
    fixture=args.fixture.resolve(strict=True);out=args.out.resolve();build=args.build.resolve(strict=True)
    if any(out.is_relative_to(p) or p.is_relative_to(out) for p in (fixture,build)):
        raise ValueError('Use a fresh sibling output directory outside fixture/build')
    text=(fixture/'simple.in').read_text();p=parameters(text)
    if not re.fullmatch(r'[A-Za-z0-9_-]+',p['run_id']):raise ValueError('Unsafe run_id')
    if p['max_level']!='6' or p['zlevels']!='30' or p['Nsoil']!='10':raise ValueError('Expected J4/J6 Z30 soil10 fixture')
    if args.ranks<1 or args.ranks>40:raise ValueError('J4 rank count invalid')
    binary=build/args.binary
    if not binary.is_file():raise ValueError('Missing requested binary')
    out.mkdir(parents=True,exist_ok=False)
    (out/'simple.in').write_text(text)
    (out/'grids').symlink_to((fixture/'grids').resolve(strict=True),target_is_directory=True)
    seed=None
    if int(p['resume'])>=0:
        seed=checkpoint_name(text);shutil.copy2(fixture/seed,out/seed)
    shutil.copy2(binary,out/'climateJ4')
    env=environment(oracles=args.oracles,detail=args.profile)
    result={'status':'running','binary':str(binary),'binary_sha256':digest(binary),
            'build':str(build),'param':'param_J4','parameters':p,'ranks':args.ranks,
            'input_sha256':digest(out/'simple.in'),'grids':grid_identity(fixture/'grids'),
            'seed_sha256':digest(out/seed) if seed else None,
            'environment':{k:v for k,v in env.items() if k.startswith('WAVETRISK_') or k.endswith('NUM_THREADS') or k=='VECLIB_MAXIMUM_THREADS'},
            'source_sha256':{str(f.relative_to(build)):digest(f) for f in sorted((build/'src').rglob('*')) if f.is_file() and f.suffix.lower()=='.f90'}}
    def save():(out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    save();monitor=MemoryMonitor(out);start=time.monotonic()
    with (out/'run.log').open('x') as log:
        process=subprocess.Popen(['mpirun','-n',str(args.ranks),'./climateJ4','simple.in'],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT)
        monitor.start(process.pid);code=process.wait()
    result['memory']=monitor.finish();result['elapsed_seconds']=time.monotonic()-start;result['returncode']=code
    a=analyze(out/'run.log');result['analysis']=a
    log=(out/'run.log').read_text()
    result['status']='passed' if code==0 and a['completed'] else 'failed';save()
    if result['status']!='passed':raise RuntimeError('Run failed: '+str(out/'run.log'))
    if not re.search(r'min_level\s*=\s*4\b',log) or not re.search(r'DOMAIN_LEVEL\s*=\s*1\b',log):
        result['status']='wrong_build';save();raise ValueError('Executable is not param_J4')
    result['refinement']=[refinement(f) for f in sorted(out.glob('*_checkpoint_*.bin.zst')) if f.name!=seed]
    result['vtk_refinement']=[vtk_refinement(f) for f in sorted(out.glob('*_tri_*.vtk.tgz'))]
    result['reloads']=list(map(int,re.findall(r'Restarting from checkpoint\s+(\d+)',log)))
    result['remaps']=log.count('Remapping vertical coordinates');save()
    print(json.dumps({k:result[k] for k in ['status','elapsed_seconds','reloads','remaps','refinement']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--binary',default='bin/climate')
    parser.add_argument('--fixture',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--ranks',type=int,default=4)
    parser.add_argument('--oracles',action='store_true')
    parser.add_argument('--profile',action='store_true')
    run(parser.parse_args())
