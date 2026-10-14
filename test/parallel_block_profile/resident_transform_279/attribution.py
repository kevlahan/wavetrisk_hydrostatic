from pathlib import Path
import re,json,statistics,sys,math,hashlib

def analyze(directory):
 """Attribute ordinary steps on the same slowest rank; never sum rank maxima.

 Input is an exact-science-validated run with the existing bounded detail profile.
 This is diagnostic attribution, not a transport-only or production speedup claim.
 """
 d=Path(directory);text=(d/'run.log').read_text();v=json.loads((d/'validation.json').read_text());keep={x['index']+1 for x in v['ordinary']['included']}; ranks=v['ranks']
 if v['status']!='passed' or v.get('oracles'):raise ValueError('Require successful oracle-free science validation')
 if re.search(r'detail dropped step records = [1-9]',text):raise ValueError('Profile overflow')
 names={};steps={}
 for line in text.splitlines():
  m=re.match(r'^\s*(\d+)\s+(.{32})\s+[-+\d.E]+',line)
  if m:names[int(m[1])]=m[2].strip()
  m=re.match(r'\s*detail step (\d+) rank (\d+) restarts (\d+) remaps (\d+) wall\s+(\S+) CPU\s+(\S+)',line)
  if m:
   key=(int(m[1]),int(m[2]))
   if key in steps:raise ValueError('Duplicate rank/step')
   steps[key]=dict(wall=float(m[5]),cpu=float(m[6]),remaps=int(m[4]))
  m=re.match(r'\s*detail step (\d+) rank (\d+) self-wall = (.*)',line)
  if m:steps[(int(m[1]),int(m[2]))]['self']=[float(x) for x in m[3].split()]
 if len(steps)!=12*ranks or len(keep)!=9:raise ValueError('Incomplete twelve-step profile')
 if set(steps)!={(s,r) for s in range(1,13) for r in range(ranks)}:raise ValueError('Missing rank/step')
 widths={len(row.get('self',[])) for row in steps.values()}
 if len(widths)!=1 or min(widths)<67:raise ValueError('Missing or inconsistent self-time rows')
 for row in steps.values():
  if not math.isfinite(row['wall']) or row['wall']<=0 or any(not math.isfinite(x) or x<0 for x in row['self']):raise ValueError('Invalid profile time')
  if abs(sum(row['self'])-row['wall'])>0.002*row['wall']:raise ValueError('Self times do not partition the step')
 slow=[]
 for s in sorted(keep):
  rows=[(r,row) for (step,r),row in steps.items() if step==s];r,row=max(rows,key=lambda x:x[1]['wall']);slow.append((r,row))
 n=len(slow[0][1]['self']);means={}
 for i in range(n):
  values=[row['self'][i] for r,row in slow]
  means[i+1]=dict(name=names.get(i+1,str(i+1)),mean_on_slowest_step_rank_seconds=statistics.mean(values),mean_all_ranks_seconds=statistics.mean(row['self'][i] for (s,r),row in steps.items() if s in keep),max_rank_ordinary_mean_seconds=max(statistics.mean(steps[s,r]['self'][i] for s in keep) for r in range(ranks)))
 groups={}
 for label,event in [('non_remap',False),('remap',True)]:
  selected=[(step,rank,row) for step,(rank,row) in zip(sorted(keep),slow) if bool(row['remaps'])==event]
  if selected:
   groups[label]=dict(step_indices=[step-1 for step,rank,row in selected],samples=len(selected),
    mean_slowest_step_seconds=statistics.mean(row['wall'] for step,rank,row in selected),
    mean_region_seconds={str(i+1):statistics.mean(row['self'][i] for step,rank,row in selected) for i in range(n)})
 ordered=sorted(means.items(),key=lambda x:-x[1]['mean_on_slowest_step_rank_seconds'])
 rank_means={str(r):dict(whole_step_seconds=statistics.mean(steps[s,r]['wall'] for s in keep),
  region_seconds={str(i+1):statistics.mean(steps[s,r]['self'][i] for s in keep) for i in range(n)}) for r in range(ranks)}
 result=dict(status='passed',ranks=ranks,rank_ordinary_means=rank_means,remap_groups=groups,fine_exchange_profile=v.get('environment',{}).get('WAVETRISK_PROFILE_RHS_FINE','legacy-full'),binary_sha256=v['binary_sha256'],log_sha256=hashlib.sha256((d/'run.log').read_bytes()).hexdigest(),ordinary_indices=sorted(keep),ordinary_bare_log_median=v['ordinary']['median'],mean_slowest_step_seconds=statistics.mean(row['wall'] for r,row in slow),slowest_ranks=[r for r,row in slow],regions=means,scope='Profiled ordinary steps; hardware/ranks belong to the input run. Self times from the SAME slowest whole-step rank, not sums of region maxima. MPI wait includes producer lateness; process CPU includes MPI spinning. No isolated transport attribution or bare speedup forecast.')
 (d/'attribution.json').write_text(json.dumps(result,indent=2)+'\n')
 print(d.name,'whole',result['mean_slowest_step_seconds'])
 for i,r in ordered[:22]:print(i,r['name'],round(r['mean_on_slowest_step_rank_seconds'],5),round(r['mean_all_ranks_seconds'],5))
if __name__=='__main__':analyze(sys.argv[1])
