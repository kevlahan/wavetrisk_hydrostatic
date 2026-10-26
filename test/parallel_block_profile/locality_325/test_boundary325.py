import copy,json,unittest
from pathlib import Path
from reconcile_boundary325 import reconcile_boundary325,prefix

HERE=Path(__file__).resolve().parent;P='block_writeback_plan%'

def fixture(owners):
    contexts=[{'blocks':[]} for _ in range(256)];blocks={}
    for cat,rank in enumerate(owners,1):
        # Distinct source boundary IDs and variable sizes; cannot pass by equal widths.
        bb=[[cat%7,5,100+cat,2,2,4,0],[cat%7+9,5,200+cat,1,2+cat%3,2+cat%3,4]]
        identity=[cat-1,(cat-1)//4,2+(cat-1)%4,5,2,2,1,-10,41,1,3,-10,41]
        block={'catalog':cat,'fields':{'identity':identity,'bdry_storage':[x for b in bb for x in b]}}
        contexts[rank]['blocks'].append(block);blocks[cat]=bb
    plans=[]
    for rank in range(256):
        p={};man={}
        for direction in ('send','recv'):
            cats=sorted([cat for cat in range(1,641)
                  if (owners[cat-1]==rank and (cat-1)//4!=rank if direction=='send'
                  else (cat-1)//4==rank and owners[cat-1]!=rank)],
                  key=lambda cat:((cat-1)//4 if direction=='send' else owners[cat-1],cat))
            p['n_'+direction]=[len(cats)];p[P+direction+'_block']=cats or [0]
            counts=[len(blocks[c]) for c in cats];nn=[sum(b[5] for b in blocks[c]) for c in cats]
            values={'count':counts,'displ':prefix(counts),'scalar_nvalue':[82*x for x in nn],
                    'vector_nvalue':[123*x for x in nn],
                    'source':[b[0] for c in cats for b in blocks[c]],
                    'level':[b[1] for c in cats for b in blocks[c]],
                    'nnode':[b[5] for c in cats for b in blocks[c]]}
            for name,val in values.items():p[P+direction+'_boundary_'+name]=val or [-1 if name=='level' else 0]
            p['n_'+direction+'_boundary']=[sum(counts)];man[direction]=(cats,nn)
        for comp,width in (('scalar',82),('vector',123)):
            for flow,direction in (('domain_send','recv'),('block_recv','send')):
                cats,nn=man[direction];counts=[0]*256
                for c,number in zip(cats,nn):counts[owners[c-1] if flow=='domain_send' else (c-1)//4]+=width*number
                field=P+'boundary_'+comp+'_'+flow
                p[field+'_count']=counts;p[field+'_displ']=prefix(counts);p[field+'_buffer.shape']=[max(1,sum(counts))]
        p['n_retained']=[sum(r==rank==(c-1)//4 for c,r in enumerate(owners,1))]
        plans.append(p)
    return plans,contexts

class BoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.maps={arm:json.loads((HERE/(arm+'.json')).read_text()) for arm in ('C0','L256')}
        cls.fixtures={arm:fixture(owners) for arm,owners in cls.maps.items()}
    def validate(self,arm,fixtures=None):
        p,c=fixtures or self.fixtures[arm]
        return reconcile_boundary325(p,c,lambda cat:self.maps[arm][cat-1])
    def test_complete_permutation_invariant_physical_providers(self):
        a=self.validate('C0');b=self.validate('L256');self.assertEqual(a,b);self.assertEqual(len(a),1280)
    def test_corruptions(self):
        for mutation in ('source','level','nnode','offset','scalar_width','vector_count','prefix','retained','capacity','duplicate_block','missing_block'):
            with self.subTest(mutation=mutation):
                plans,ctx=copy.deepcopy(self.fixtures['L256'])
                rank=next(r for r,p in enumerate(plans) if p['n_send'][0]>0);p=plans[rank]
                field={'source':'send_boundary_source','level':'send_boundary_level','nnode':'send_boundary_nnode',
                       'offset':'send_boundary_displ','scalar_width':'send_boundary_scalar_nvalue',
                       'vector_count':'boundary_vector_block_recv_count','prefix':'boundary_scalar_domain_send_displ',
                       'capacity':'boundary_scalar_block_recv_buffer.shape'}
                if mutation in field:p[P+field[mutation]][0]+=1
                elif mutation=='retained':p['n_retained'][0]+=1
                elif mutation=='duplicate_block':ctx[rank]['blocks'].append(copy.deepcopy(ctx[rank]['blocks'][0]))
                elif mutation=='missing_block':ctx[rank]['blocks'].pop()
                with self.assertRaises(ValueError):self.validate('L256',(plans,ctx))

if __name__=='__main__':unittest.main(verbosity=2)
