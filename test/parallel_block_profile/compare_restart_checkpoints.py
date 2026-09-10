"""Compare semantic payloads of the J5 climate checkpoint fixture, including poles.

This deliberately supports only the three-field, PATCH_SIZE=4, no-TKE climate
layout with an empty testcase checkpoint extension. Reject other layouts.
Directory load weights are metadata, not numerical solution equivalence.
"""
import argparse
import json
from pathlib import Path
import struct
import subprocess
import numpy as np


def read(path, zmin, zmax):
    raw = subprocess.check_output(['zstd', '-dc', str(path)])
    magic, version, nd = struct.unpack_from('<3q', raw)
    if magic != 0x5741564554524953 or version != 1 or nd != 160:
        raise ValueError('Expected J5 climate checkpoint v1 with 160 Domains')
    offset = np.frombuffer(raw, '<i8', nd, 24+4*nd)
    size = np.frombuffer(raw, '<i8', nd, 24+12*nd)
    expected = 24+20*nd
    for start,length in zip(offset,size):
        if int(start)!=expected or int(length)<=0:
            raise ValueError('Noncontiguous or invalid checkpoint directory')
        expected += int(length)
    if expected!=len(raw):
        raise ValueError('Truncated checkpoint or unsupported trailing extension')
    nk = zmax-zmin+1
    metadata = 24+3*nk*8
    coarse_bytes = nk*80*8
    record_bytes = coarse_bytes+16
    pole_bytes = nk*2*8
    result = []
    for start, length in zip(offset, size):
        payload = raw[int(start):int(start+length)]
        if len(payload) != length:
            raise ValueError('Truncated Domain')
        remainder = len(payload)-metadata-coarse_bytes
        candidates = [n for n in range(3) if remainder>=n*pole_bytes and (remainder-n*pole_bytes)%record_bytes==0]
        if len(candidates) != 1:
            raise ValueError('Unsupported Domain payload layout')
        npole = candidates[0]
        header = struct.unpack_from('<idqi', payload)
        thresholds = np.frombuffer(payload, '<f8', 3*nk, 24)
        pos = metadata
        poles = np.frombuffer(payload, '<f8', npole*nk*2, pos).reshape(npole,nk,2)
        pos += npole*pole_bytes
        coarse = np.frombuffer(payload, '<f8', nk*80, pos).reshape(1,nk,80)
        pos += coarse_bytes
        wavelets, topology = [], []
        while pos<len(payload):
            wavelets.append(np.frombuffer(payload, '<f8', nk*80, pos).reshape(nk,80))
            pos += coarse_bytes
            topology.append(np.frombuffer(payload, '<i4', 4, pos))
            pos += 16
        result.append((header,thresholds,poles,coarse,np.array(wavelets),np.array(topology)))
    return result


def compare(left, right, zmin=-10, zmax=30):
    a, b = read(left,zmin,zmax), read(right,zmin,zmax)
    result = {'left':str(left),'right':str(right),'header_differences':0,
              'topology_differences':0,'threshold_max_abs':0.,'threshold_nonfinite':0,'fields':{}}
    for gid,(aa,bb) in enumerate(zip(a,b)):
        result['header_differences'] += int(aa[0]!=bb[0])
        result['threshold_nonfinite'] += int(np.count_nonzero(~np.isfinite(aa[1]) | ~np.isfinite(bb[1])))
        result['threshold_max_abs'] = max(result['threshold_max_abs'],float(np.max(np.abs(aa[1]-bb[1]))))
        if aa[5].shape!=bb[5].shape:
            raise ValueError(f'Domain {gid}: checkpoint tree shapes differ')
        result['topology_differences'] += int(np.count_nonzero(aa[5]!=bb[5]))
        for family,index,fields in [('pole',2,[('mass',0,1),('temperature',1,1)]),
                                    ('coarse',3,[('velocity',0,48),('mass',48,16),('temperature',64,16)]),
                                    ('wavelet',4,[('velocity',0,48),('mass',48,16),('temperature',64,16)])]:
            if aa[index].shape!=bb[index].shape: raise ValueError('Field layouts differ')
            if not aa[index].size: continue
            for field,first,count in fields:
                for zone,levels in [('soil',slice(0,1-zmin)),('atmosphere',slice(1-zmin,None))]:
                    av=aa[index][:,levels,first:first+count];bv=bb[index][:,levels,first:first+count]
                    finite=np.isfinite(av)&np.isfinite(bv)
                    delta=np.where(finite,np.abs(av-bv),0)
                    name='/'.join((family,field,zone))
                    s=result['fields'].setdefault(name,{'different':0,'nonfinite':0,'max_abs':0.,'worst':None})
                    s['different']+=int(np.count_nonzero(delta))
                    s['nonfinite']+=int(np.count_nonzero(~finite))
                    if delta.size and delta.max()>s['max_abs']:
                        loc=np.unravel_index(delta.argmax(),delta.shape)
                        s['max_abs']=float(delta[loc])
                        s['worst']={'domain':gid,'record':int(loc[0]),'k':int(loc[1]+(zmin if zone=='soil' else 1)),
                                    'index':int(loc[2]),'left':float(av[loc]),'right':float(bv[loc])}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--left',type=Path,required=True)
    parser.add_argument('--right',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(compare(args.left,args.right),indent=2))
