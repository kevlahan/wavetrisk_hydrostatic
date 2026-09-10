#!/usr/bin/env python3
"""Compare read-only Domain interior snapshots by global domain/tree key.

Boundary aliases are not dumped. Pre-writeback Domain data may be stale by
design; such a difference alone is NOT a failure of authoritative block state.
The reported 1e-12 scaled count is diagnostic, not an acceptance tolerance.
"""
import argparse,glob,json,struct
from pathlib import Path
import numpy as np


def load(pattern):
    patches={}
    headers=[]
    for name in sorted(glob.glob(pattern)):
        raw=Path(name).read_bytes()
        magic,nd,nv,zlo,zhi,ps,step,remap=struct.unpack_from('<8i',raw)
        if magic!=179401 or nv!=3 or ps!=4:raise ValueError('Unsupported probe layout')
        nk=zhi-zlo+1
        now=struct.unpack_from('<d',raw,32)[0]
        headers.append({'step':step,'remap':remap,'time':now})
        pos=40+8*nv*nk
        for d in range(nd):
            gid,count=struct.unpack_from('<2i',raw,pos);pos+=8
            for p in range(count):
                key,level=struct.unpack_from('<qi',raw,pos);pos+=12
                mask=np.frombuffer(raw,'<i4',64,pos);pos+=256
                values=np.frombuffer(raw,'<f8',nk*160,pos).reshape(nk,160);pos+=nk*1280
                if (gid,key) in patches:raise ValueError('Duplicate patch key')
                patches[gid,key]=(level,zlo,mask,values)
        if pos!=len(raw):raise ValueError('Probe trailing/truncated bytes')
    if not headers:raise ValueError('No probe files match')
    return patches,headers


def compare(left,right):
    a,ah=load(left);b,bh=load(right)
    stats={}
    mask_differences=0
    for key in sorted(a.keys() & b.keys()):
        level,zlo,mask,av=a[key];bl,bz,bmask,bv=b[key]
        if (level,zlo,av.shape)!=(bl,bz,bv.shape):raise ValueError('Patch layout differs')
        mask_differences+=int(np.count_nonzero(mask!=bmask))
        for field,start,n,mask_slice in [('velocity',0,48,slice(16,64)),('mass',96,16,slice(0,16)),
                                       ('temperature',128,16,slice(0,16))]:
            for family,offset in [('sol',0),('wav',n)]:
                aa=av[:,start+offset:start+offset+n];bb=bv[:,start+offset:start+offset+n]
                finite=np.isfinite(aa)&np.isfinite(bb)
                delta=np.zeros_like(aa)
                np.subtract(aa,bb,out=delta,where=finite);np.abs(delta,out=delta)
                scale=np.maximum(1,np.maximum(np.abs(aa),np.abs(bb)))
                for zone in ('soil','atmosphere'):
                    vertical=np.arange(zlo,zlo+aa.shape[0])
                    layers=(vertical<=0) if zone=='soil' else (vertical>0)
                    for active_only in (False,True):
                        select=np.broadcast_to(layers[:,None],aa.shape).copy()
                        if active_only:select &= (mask[mask_slice]>=2)[None,:] & (bmask[mask_slice]>=2)[None,:]
                        name='/'.join([family,field,zone,'active' if active_only else 'all'])
                        s=stats.setdefault(name,{'samples':0,'different':0,'scaled_above_1e12':0,
                                                'nonfinite':0,'max_abs':0.,'max_scaled':0.,'worst':None})
                        s['samples']+=int(np.count_nonzero(select))
                        s['nonfinite']+=int(np.count_nonzero(select&~finite))
                        s['different']+=int(np.count_nonzero(select&finite&(delta!=0)))
                        s['scaled_above_1e12']+=int(np.count_nonzero(select&finite&(delta>1e-12*scale)))
                        vals=np.where(select&finite,delta,0)
                        if vals.max()>s['max_abs']:
                            k,i=np.unravel_index(vals.argmax(),vals.shape)
                            s['max_abs']=float(vals[k,i]);s['worst']={'domain':key[0],'tree_key':key[1],'level':level,
                                'k':int(k+zlo),'index':int(i),'left':float(aa[k,i]),'right':float(bb[k,i])}
                        s['max_scaled']=max(s['max_scaled'],float(np.max(np.where(select&finite,delta/scale,0))))
    return {'left':left,'right':right,'left_headers':ah,'right_headers':bh,
            'left_only_patches':len(a.keys()-b.keys()),'right_only_patches':len(b.keys()-a.keys()),
            'common_patches':len(a.keys()&b.keys()),'mask_differences':mask_differences,'fields':stats}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--left',required=True);p.add_argument('--right',required=True)
    args=p.parse_args();print(json.dumps(compare(args.left,args.right),indent=2))
