"""Read-only, isolated physics bridge probes. No numerical acceptance thresholds."""
import argparse
import glob
import json
from pathlib import Path
import struct
import numpy as np


def load(pattern):
    records={}
    files=glob.glob(pattern)
    for name in files:
        data=Path(name).read_bytes();pos=0
        while pos<len(data):
            h=struct.unpack_from('<11i',data,pos);pos+=44
            magic,step,gid,node,patch,i,j,pole,mask,nz,ns=h
            if magic!=181402 or nz<1 or ns<0: raise ValueError('Invalid physics probe')
            count=3+8*nz
            raw=np.frombuffer(data,'<f8',count,pos);pos+=8*count
            count=6+8*nz+ns
            packed=np.frombuffer(data,'<f4',count,pos);pos+=4*count
            fields={};at=0
            for field,n in [('dt',1),('day',1),('day_fraction',1),('Phisurf',1),
                            ('Play',nz),('Pint',nz+1),('Phi',nz),('Umag',nz),
                            ('U',nz),('V',nz),('W',nz),('Theta',nz),('Tsoil',ns+1)]:
                fields[field]=packed[at:at+n];at+=n
            if at!=count: raise ValueError('Unexpected physics layout')
            key=(step,gid,node)
            if key in records: raise ValueError('Duplicate physics column')
            records[key]={'header':h,'raw':raw,'packed':packed,'fields':fields}
        if pos!=len(data): raise ValueError('Truncated physics probe')
    if not records: raise ValueError('No physics columns found')
    return records


def compare(left,right):
    a=load(left);b=load(right)
    if a.keys()!=b.keys(): raise ValueError('Physics column keys differ')
    result={'columns':len(a),'header_differences':0,'raw_max_abs':0.,'fields':{}}
    for name in next(iter(a.values()))['fields']:
        result['fields'][name]={'different':0,'max_abs':0.,'first':None}
    for key,av in sorted(a.items()):
        bv=b[key]
        result['header_differences']+=int(av['header']!=bv['header'])
        if not np.isfinite(av['raw']).all() or not np.isfinite(bv['raw']).all():
            raise ValueError('Nonfinite raw physics input')
        result['raw_max_abs']=max(result['raw_max_abs'],float(np.max(np.abs(av['raw']-bv['raw']))))
        if not np.isfinite(av['packed']).all() or not np.isfinite(bv['packed']).all():
            raise ValueError('Nonfinite physics array')
        if np.array_equal(av['packed'],bv['packed']):continue
        for name,aa in av['fields'].items():
            bb=bv['fields'][name]
            if not np.isfinite(aa).all() or not np.isfinite(bb).all(): raise ValueError('Nonfinite physics array')
            delta=np.abs(aa.astype('f8')-bb.astype('f8'))
            s=result['fields'].setdefault(name,{'different':0,'max_abs':0.,'first':None})
            s['different']+=int(np.count_nonzero(delta))
            s['max_abs']=max(s['max_abs'],float(delta.max()))
            if np.any(delta) and s['first'] is None:
                loc=int(np.flatnonzero(delta)[0])
                s['first']={'step':key[0],'domain':key[1],'node':key[2],'index':loc,
                            'header':av['header'],'left':float(aa[loc]),'right':float(bb[loc]),
                            'left_bits':int(aa.view('u4')[loc]),'right_bits':int(bb.view('u4')[loc])}
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--left',required=True);p.add_argument('--right',required=True)
    args=p.parse_args();print(json.dumps(compare(args.left,args.right),indent=2))
