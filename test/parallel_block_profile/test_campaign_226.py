"""Campaign coverage and full-node physical placement for production integration."""
import unittest
from capability_226 import schedule
from placement_226 import plan

class Campaign(unittest.TestCase):
    def test_higher_ratios_have_checked_and_plain_runs(self):
        s=list(schedule())
        self.assertEqual(len(s),21)
        for level,ranks in [('J4',41),('J4',80),('J4',120),('J4',160),('J5',167),('J5',256)]:
            for kind in ('checked','candidate'):self.assertIn((level,ranks,kind,False,'native'),s)
        self.assertIn(('J4',160,'candidate',True,'native'),s)
        self.assertIn(('J5',256,'candidate',True,'native'),s)
        self.assertIn(('J4',41,'checked',False,'empty-root'),s)
        for level in ('J4','J5'):
            self.assertEqual(next(x for x in s if x[0]==level)[2],'baseline')
    def test_full_node_uses_every_core_once(self):
        top=dict(host='test',rows=[dict(cpu=i,core=i%128,socket=i//128,node=i//128,llc=[i//8]) for i in range(256)])
        cases=plan(top)
        self.assertEqual(set(cases),{40,41,80,120,160,167,256})
        self.assertEqual(cases[256]['socket_rank_counts'],[128,128])
        for n,c in cases.items():
            self.assertEqual(len({r['cpu'] for r in c['mapping']}),n)
            self.assertEqual(len({(r['socket'],r['core']) for r in c['mapping']}),n)

if __name__=='__main__':unittest.main()
