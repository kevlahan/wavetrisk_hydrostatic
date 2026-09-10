import shutil
import struct
import subprocess
import tempfile
import unittest
from pathlib import Path
import numpy as np
from compare_restart_checkpoints import compare


class CheckpointComparisonTests(unittest.TestCase):
    def write_fixture(self, path, perturb=False):
        chunks=[]
        for d in range(160):
            data=struct.pack('<idqi',252,100.,10000,3)+np.zeros(6,dtype='<f8').tobytes()
            if d==0:
                pole=np.zeros((1,2,2),dtype='<f8')
                if perturb: pole[0,1,1]=0.25
                data+=pole.tobytes()
            coarse=np.zeros((2,80),dtype='<f8')
            if perturb and d==1: coarse[1,48]=1.5
            data+=coarse.tobytes()
            wave=np.zeros((2,80),dtype='<f8')
            if perturb and d==2: wave[1,64]=2.5
            data+=wave.tobytes()+np.zeros(4,dtype='<i4').tobytes()
            chunks.append(data)
        sizes=np.array([len(p) for p in chunks],dtype='<i8')
        offsets=np.concatenate(([24+20*160],24+20*160+np.cumsum(sizes)[:-1])).astype('<i8')
        raw=struct.pack('<3q',0x5741564554524953,1,160)+np.zeros(160,dtype='<i4').tobytes()
        raw+=offsets.tobytes()+sizes.tobytes()+b''.join(chunks)
        subprocess.run(['zstd','-q','-o',str(path)],input=raw,check=True)

    @unittest.skipUnless(shutil.which('zstd'),'zstd required')
    def test_identity_poles_and_field_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/'a.zst',Path(tmp)/'b.zst'
            self.write_fixture(a)
            self.write_fixture(b,True)
            same=compare(a,a,0,1)
            self.assertEqual(same['header_differences'],0)
            self.assertEqual(same['topology_differences'],0)
            self.assertTrue(all(s['different']==0 for s in same['fields'].values()))
            r=compare(a,b,0,1)
            self.assertEqual(r['fields']['pole/temperature/atmosphere']['max_abs'],0.25)
            self.assertEqual(r['fields']['coarse/mass/atmosphere']['max_abs'],1.5)
            self.assertEqual(r['fields']['wavelet/temperature/atmosphere']['max_abs'],2.5)
            self.assertEqual(r['fields']['wavelet/velocity/atmosphere']['max_abs'],0.)


if __name__=='__main__': unittest.main()
