#!/usr/bin/env python3
"""Parser adversarial tests and a source-derived 256-rank conservation fixture."""
import copy,json,tempfile,unittest
from pathlib import Path
import account325 as a

HERE=Path(__file__).resolve().parent

def lines_for(rank,records):
    lines=[f'RK325WORK1 {rank}']; rows=0
    for tag,phase,calls,peak,peers in records:
        n=len(peak);lines.append(' '.join(map(str,['C',tag,phase,calls,n,*[x*calls for x in peak],*peak])));rows+=1
        for peer,(sent,recv) in sorted(peers.items()):
            if sent==recv==0:continue
            lines.append(' '.join(map(str,['P',tag,phase,peer,sent*calls,recv*calls,sent,recv])));rows+=1
    lines.append(f'COMPLETE {len(records)} {rows}')
    return '\n'.join(lines)+'\n'

class ParserTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)/'work325-rank-0.tsv'
        self.schema={'example':{'fields':['a','b'],'peer_unit':'bytes'}}
        self.text=lines_for(0,[('example',3,7,[4,5],{0:(2,2),1:(3,4)})])
    def tearDown(self):self.tmp.cleanup()
    def parse(self,text):self.p.write_text(text);return a.parse(self.p,0,self.schema)
    def test_valid(self):
        r=self.parse(self.text);self.assertEqual(r['entries']['example',3]['sum'],[28,35])
    def test_mutations(self):
        changes={
          'wrong_rank':self.text.replace('RK325WORK1 0','RK325WORK1 1'),
          'wrong_magic':self.text.replace('RK325WORK1','RK325WORK2'),
          'unknown_tag':self.text.replace('example','unknown'),
          'negative':self.text.replace('C example 3 7 2','C example -3 7 2'),
          'overflow':self.text.replace('28 35','9223372036854775808 35'),
          'zero_calls':self.text.replace('C example 3 7 2','C example 3 0 2'),
          'wrong_width':self.text.replace('C example 3 7 2','C example 3 7 3'),
          'counter_sum':self.text.replace('28 35 4 5','29 35 4 5'),
          'counter_peak':self.text.replace('28 35 4 5','28 35 29 5'),
          'peer_sum':self.text.replace('P example 3 1 21 28 3 4','P example 3 1 22 28 3 4'),
          'peer_peak':self.text.replace('P example 3 1 21 28 3 4','P example 3 1 21 28 30 4'),
          'zero_peer':self.text.replace('P example 3 1 21 28 3 4','P example 3 1 0 0 0 0'),
          'missing_footer':self.text.rsplit('COMPLETE',1)[0],
          'wrong_footer':self.text.replace('COMPLETE 1 3','COMPLETE 1 2'),
          'trailing':self.text+'C example 3 1 2 1 1 1 1\n',
          'duplicate_counter':self.text.replace('COMPLETE','C example 3 7 2 28 35 4 5\nCOMPLETE'),
          'duplicate_peer':self.text.replace('COMPLETE','P example 3 1 21 28 3 4\nCOMPLETE'),
          'orphan_peer':self.text.replace('P example 3 1','P example 4 1'),
          'noninteger':self.text.replace('C example 3 7','C example 3 seven'),
        }
        for name,text in changes.items():
            with self.subTest(name=name),self.assertRaises(ValueError):self.parse(text)

def make_run(path,owners,schema):
    path.mkdir(parents=True,exist_ok=True);nr=256;epochs=5
    source=[{} for _ in range(nr)]; incoming=[{} for _ in range(nr)]
    outblocks=[0]*nr;inblocks=[0]*nr;retained=[0]*nr
    geometry=[{} for _ in range(nr)];ingeometry=[{} for _ in range(nr)]
    for b,dest in enumerate(owners,1):
        home=(b-1)//4
        geometry[home][dest]=10;ingeometry[dest][home]=10
        if home==dest:retained[home]+=1;continue
        source[home][dest]=source[home].get(dest,0)+100+b
        incoming[dest][home]=incoming[dest].get(home,0)+100+b
        outblocks[home]+=1;inblocks[dest]+=1
    for rank in range(nr):
        peers={p:(source[rank].get(p,0),incoming[rank].get(p,0)) for p in source[rank].keys()|incoming[rank].keys()}
        sb=sum(source[rank].values());rb=sum(incoming[rank].values())
        rec=[('migration',0,epochs,[4 if rank<160 else 0,retained[rank],outblocks[rank],inblocks[rank],sb,rb,
                                   owners.count(rank),17,max(1,sb)+max(1,rb)],peers)]
        for b in range(rank*4+1,min(rank*4+5,641)):
            rec.append(('migration.source',b,epochs,[100+b,int(owners[b-1]!=rank)],{}))
        g=geometry[rank];i=ingeometry[rank]
        rec.append(('rhs.geometry',0,epochs,[len(g),sum(g.values()),sum(i.values()),g.get(rank,0),i.get(rank,0)],
                    {p:(g.get(p,0),i.get(p,0)) for p in g.keys()|i.keys()}))
        for tag in ('transform.geometry','route.exchange','transform.workspace','selector.call'):
            rec.append((tag,0,epochs,[0]*len(schema[tag]['fields']),{}))
        (path/f'work325-rank-{rank}.tsv').write_text(lines_for(rank,rec))

class ConservationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema=json.loads((HERE/'fixture/work325-schema.json').read_text())['records']
        cls.maps={k:json.loads((HERE/(k+'.json')).read_text()) for k in ('C0','L256')}
    def test_both_maps_and_byte_peaks(self):
        with tempfile.TemporaryDirectory() as tmp:
            reports={}
            for arm in self.maps:
                path=Path(tmp)/arm;make_run(path,self.maps[arm],self.schema)
                reports[arm]=a.validate_run(path,schema=self.schema,arm=arm)
            self.assertEqual(reports['C0']['summary']['migration:0']['sum'][2],634*5)
            self.assertEqual(reports['L256']['summary']['migration:0']['sum'][2],240*5)
            # Every block differs in size: do not multiply a block-count peak by average bytes.
            self.assertGreater(reports['C0']['summary']['migration:0']['max_single_call'][4],400)
    def test_three_run_qualification_and_source_extent_mutation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)
            (path/'work325-schema.json').write_text(json.dumps({'records':self.schema}))
            for label in a.LABELS:
                make_run(path/label,self.maps['C0' if label=='qual-C0' else 'L256'],self.schema)
            report=a.qualify(path,path)
            self.assertTrue(report['passed']);self.assertFalse(report['performance_qualified'])
            p=path/'checked-L256/work325-rank-0.tsv';text=p.read_text()
            old='C migration.source 1 5 2 505 0 101 0';new='C migration.source 1 5 2 510 0 102 0'
            self.assertIn(old,text);p.write_text(text.replace(old,new))
            with self.assertRaisesRegex(ValueError,'paired_source_serialized_extents'):a.qualify(path,path)
    def test_cross_rank_mutations(self):
        for mutation in ('drop_rank','wrong_peer','payload_omission','wrong_rebuilds','wrong_arm'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp);make_run(path,self.maps['L256'],self.schema)
                p=path/'work325-rank-0.tsv';text=p.read_text();kwargs={}
                if mutation=='drop_rank':p.unlink()
                elif mutation=='wrong_peer':
                    text=text.replace('P migration 0 1 ', 'P migration 0 255 ');p.write_text(text)
                elif mutation=='payload_omission':
                    text=text.replace('C migration.source 1 5 2 505 0 101 0','C migration.source 1 5 2 500 0 100 0');p.write_text(text)
                    # Block1 is retained, so mutate the first remotely sent source instead.
                    p=path/'work325-rank-1.tsv';text=p.read_text().replace('C migration.source 7 5 2 535 5 107 1',
                        'C migration.source 7 5 2 530 5 106 1');p.write_text(text)
                elif mutation=='wrong_rebuilds':kwargs['expected_rebuilds']=6
                elif mutation=='wrong_arm':kwargs['arm']='C0'
                with self.assertRaises(ValueError):a.validate_run(path,schema=self.schema,arm=kwargs.pop('arm','L256'),**kwargs)

if __name__=='__main__':unittest.main(verbosity=2)
