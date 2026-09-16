"""Protocol and substantial-gain decisions using independent synthetic timings."""
import unittest
from cluster_203 import schedule,KINDS
from analyze_batching_203 import summarize,PHASES,inverse_ledger,phase_balance
from test_campaign_202 import evidence

def records(gain=.15):
 r=evidence(gain)
 phases={p:dict(wall_seconds_rank_mean=.1,counts_rank_mean_per_step={}) for p in ('whole step',*PHASES)}
 r['runs']['legacy-coarse']=dict(label='legacy-coarse',side='legacy',phase='coarse',round=0,ordinary_step_median=.54,hardware={'quality_flags':[],'phases':phases})
 return r
class Campaign203(unittest.TestCase):
 def test_schedule_contains_fresh_legacy_profile_and_reversed_pairs(self):
  s=schedule({});self.assertEqual(len(s),18);self.assertEqual({x[1] for x in s},set(KINDS))
  self.assertEqual(sum(x[1]=='legacy_profile' for x in s),1)
  self.assertEqual([x[1] for x in s if x[2]=='bare' and x[3]==2],['legacy','candidate','baseline'])
 def test_phase_gap_uses_same_campaign_legacy(self):
  r=summarize(records());p=r['matched_coarse_phases']['whole step']
  self.assertAlmostEqual(p['remaining_excess'],(.91/5)-.1)
  self.assertEqual(r['decision'],'review_for_promotion')
 def test_one_percent_not_promoted(self):
  self.assertEqual(summarize(records(.01))['decision'],'do_not_promote_insufficient_gain')
 def test_absolute_rank_excess_is_not_relative_ratio(self):
  cells=[dict(step=k,rank=r,phase=p,wall_seconds=v,thread_cpu_seconds=v/2) for p in ('whole step',*PHASES) for k,values in [(1,[1.,3.]),(2,[3.,5.])] for r,v in enumerate(values)]
  b=phase_balance({'cells':cells})['phases']['whole step']['wall_seconds']
  self.assertEqual(b['rank_mean'],3.)
  self.assertEqual(b['mean_step_rank_max'],4.)
  self.assertEqual(b['mean_step_max_minus_mean'],1.)
  self.assertEqual(b['mean_step_max_over_mean'],1.375)
 def test_ledger_coverage_and_sums(self):
  text='Stage203 inverse rank=0 work= 5 6 7 8 9 10\nStage203 inverse rank=1 work= 4 3 2 1 0 0\n'
  self.assertEqual(inverse_ledger(text,2)['per_rank_totals']['0'],[5,6,7,8,9,10])
  with self.assertRaisesRegex(ValueError,'Missing'):inverse_ledger(text,3)
  with self.assertRaisesRegex(ValueError,'transaction'):inverse_ledger(text+text.splitlines()[0]+'\n',2)
 def test_legacy_observer_bias_requires_review(self):
  r=records();r['runs']['legacy-coarse']['ordinary_step_median']=.65
  self.assertEqual(summarize(r)['decision'],'review_observer_effects')
 def test_missing_legacy_profile_rejected(self):
  r=records();del r['runs']['legacy-coarse']
  with self.assertRaises(StopIteration):summarize(r)
if __name__=='__main__':unittest.main()
