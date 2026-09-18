#!/usr/bin/env python3
"""Interrupted collections never publish a final archive; readable is insufficient."""
import io,json,tarfile,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from collect_cluster_logs_215 import collect,verify_archive

class Collector(unittest.TestCase):
    def setup(self,root):
        results=root/'results';results.mkdir();(results/'results.json').write_text(json.dumps({'status':'passed','slurm_environment':{'SLURM_JOB_ID':'123'}}))
        (results/'run.log').write_text('complete run\n');return results

    def test_complete_verifies_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            home=Path(d);results=self.setup(home);archive=collect(215,'123',home,results)
            before=archive.read_bytes();m=verify_archive(archive);self.assertEqual(len(m['files']),2)
            with self.assertRaises(FileExistsError):collect(215,'123',home,results)
            self.assertEqual(before,archive.read_bytes());self.assertFalse(list(home.glob('*.partial-*')))

    def test_interrupt_does_not_publish(self):
        with tempfile.TemporaryDirectory() as d:
            home=Path(d);results=self.setup(home);original=tarfile.TarFile.addfile;calls=[]
            def stop(t,info,*args,**kwargs):
                calls.append(info.name)
                if len(calls)==2:raise KeyboardInterrupt()
                return original(t,info,*args,**kwargs)
            with patch.object(tarfile.TarFile,'addfile',stop):
                with self.assertRaises(KeyboardInterrupt):collect(215,'123',home,results)
            self.assertFalse((home/'stage215-review-123.tar.gz').exists())
            self.assertFalse(list(home.glob('.stage215-review-*.partial-*')))
            self.assertEqual(len(verify_archive(collect(215,'123',home,results))['files']),2)

    def test_readable_partial_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'partial.tar.gz'
            with tarfile.open(p,'w:gz') as t:
                data=b'{}';i=tarfile.TarInfo('results/results.json');i.size=len(data);t.addfile(i,io.BytesIO(data))
            with tarfile.open(p) as t:self.assertEqual(len(t.getmembers()),1)
            with self.assertRaisesRegex(ValueError,'Incomplete review'):verify_archive(p)

    def test_member_hash_checked(self):
        with tempfile.TemporaryDirectory() as d:
            home=Path(d);archive=collect(215,'123',home,self.setup(home));bad=home/'bad.tar.gz'
            with tarfile.open(archive) as src,tarfile.open(bad,'w:gz') as dst:
                for member in src.getmembers():
                    data=src.extractfile(member).read()
                    if member.name=='results/run.log':data=b'corrupted\n';member.size=len(data)
                    dst.addfile(member,io.BytesIO(data))
            with self.assertRaisesRegex(ValueError,'Corrupt member'):verify_archive(bad)

if __name__=='__main__':unittest.main()
