"""Raw-address diagnosis BEFORE adaptation; requires identical Domain layouts.

Unlike the canonical patch comparator, this includes aliases/poles/scaffolds.
It is not suitable for comparing different topology or patch storage orders.
"""
import argparse
import glob
import json
from pathlib import Path
import struct
import numpy as np


def load(pattern):
    domains = {}
    for name in sorted(glob.glob(pattern)):
        raw = Path(name).read_bytes()
        magic, nd, nv, zlo, zhi, step = struct.unpack_from('<6i', raw)
        if magic != 179402 or nv != 3:
            raise ValueError('Unsupported alias snapshot')
        pos = 24
        for _ in range(nd):
            gid, n = struct.unpack_from('<2i', raw, pos)
            pos += 8
            values = np.frombuffer(raw, '<f8', (zhi-zlo+1)*5*n, pos).reshape(zhi-zlo+1, 5*n)
            pos += values.nbytes
            if gid in domains:
                raise ValueError('Duplicate Domain')
            domains[gid] = zlo, step, n, values
        if pos != len(raw):
            raise ValueError('Alias layout mismatch')
    if not domains:
        raise ValueError('No alias snapshots')
    return domains


def compare(left, right):
    a, b = load(left), load(right)
    if a.keys() != b.keys():
        raise ValueError('Domain sets differ')
    result = {}
    for gid, (zlo, step, n, av) in a.items():
        bz, bs, bn, bv = b[gid]
        if (zlo,step,n,av.shape) != (bz,bs,bn,bv.shape):
            raise ValueError('Raw layout differs')
        for name, start, size in [('velocity',0,3*n),('mass',3*n,n),('temperature',4*n,n)]:
            aa, bb = av[:,start:start+size], bv[:,start:start+size]
            if not (np.isfinite(aa).all() and np.isfinite(bb).all()):
                raise ValueError('Nonfinite alias snapshot')
            delta = np.abs(aa-bb)
            s = result.setdefault(name, {'max_abs':0.,'different':0,'worst':None})
            s['different'] += int(np.count_nonzero(delta))
            if delta.max() > s['max_abs']:
                k, i = np.unravel_index(delta.argmax(), delta.shape)
                s['max_abs'] = float(delta[k,i])
                s['worst'] = {'domain':gid,'k':int(k+zlo),'raw_index':int(i),
                              'left':float(aa[k,i]),'right':float(bb[k,i])}
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--left',required=True)
    parser.add_argument('--right',required=True)
    args=parser.parse_args()
    print(json.dumps(compare(args.left,args.right),indent=2))
