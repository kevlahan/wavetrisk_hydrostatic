"""Regression tests for campaign pairing and the substantial-gain decision."""
import copy
import unittest
from cluster_202 import schedule,KINDS
from analyze_batching_202 import summarize,PHASES

def evidence(recovery=.15):
    runs={}
    for side,time in [('baseline',1.06),('candidate',1.06-recovery),('legacy',.54)]:
        for repeat in (1,2):
            label=f'{side}-bare-{repeat}'
            runs[label]=dict(label=label,side=side,phase='bare',round=repeat,ordinary_step_median=time)
        if side=='legacy':continue
        for phase,repeat in [('coarse',0),('core',1),('core',2)]:
            label=f'{side}-{phase}-{repeat}'
            counts={'instructions':1000*time,'cycles':2000*time}
            phases={p:dict(wall_seconds_rank_mean=time/5,counts_rank_mean_per_step=counts) for p in ('whole step',*PHASES)}
            runs[label]=dict(label=label,side=side,phase=phase,round=repeat,ordinary_step_median=time,
                             hardware={'quality_flags':[],'phases':phases})
    return {'runs':runs}

class Campaign202(unittest.TestCase):
    def test_frozen_versions_reversed_pairs_and_no_sampling(self):
        jobs=schedule({})
        self.assertEqual(len(jobs),17)
        self.assertEqual(set(k for _,k,_,_,_ in jobs),set(KINDS))
        self.assertEqual(len(set(jobs)),17)
        self.assertEqual([k for _,k,p,r,_ in jobs if p=='bare' and r==2],['legacy','candidate','baseline'])
        self.assertEqual([k for _,k,p,r,_ in jobs if p=='core' and r==2],['candidate_profile','baseline_profile'])
        self.assertFalse(any(p in ('sample','memory') for _,_,p,_,_ in jobs))
    def test_one_percent_gain_is_not_promoted(self):
        self.assertEqual(summarize(evidence(.011))['decision'],'do_not_promote_insufficient_gain')
    def test_each_pair_must_meet_target(self):
        r=evidence(.2);r['runs']['candidate-bare-2']['ordinary_step_median']=1.05
        self.assertGreater(summarize(r)['recovery_seconds'],.1)
        self.assertFalse(summarize(r)['substantial_gain_screen_passed'])
    def test_substantial_gain_requires_calibration_review(self):
        r=evidence();self.assertEqual(summarize(r)['decision'],'review_for_promotion')
        r['runs']['candidate-core-1']['ordinary_step_median']=1.2
        self.assertEqual(summarize(r)['decision'],'review_observer_effects')
    def test_missing_repetition_is_an_error(self):
        r=evidence();del r['runs']['baseline-bare-1']
        with self.assertRaisesRegex(ValueError,'Two bare'):summarize(r)
if __name__=='__main__':unittest.main()
