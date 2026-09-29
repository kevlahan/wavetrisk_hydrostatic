"""Stage229 run-level timing decisions. Timesteps are never repetitions."""
import math, statistics, re
from analyze import STEP

POLICY = dict(screen_gain_fraction=0.01, clear_regression_fraction=0.03,
              acceptance_seconds=0.025, near_threshold_seconds=0.005,
              maximum_confirmation_batches=2)


def precise_seconds(text, analysis):
    rows=[line.split() for line in text.splitlines() if line.strip()]
    if len(rows)!=len(analysis['step_records']):raise ValueError('test_log row count differs')
    values=[]
    for fields,record in zip(rows,analysis['step_records']):
        if len(fields)!=6:raise ValueError('Invalid test_log row')
        date,t,dt,level,dof,cpu=map(float,fields)
        if not all(math.isfinite(v) for v in (date,t,dt,level,dof,cpu)) or cpu<=0:raise ValueError('Invalid precise timing')
        state=re.fullmatch(r'\s*(\d+)\s+([\d.eE+-]+) d dt =\s*([\d.eE+-]+) s Jmax =\s*(\d+) dof =\s*(\d+)',record['state'])
        if state is None:raise ValueError('Unknown stdout state')
        sd,st,sdt,sl,sn=map(float,state.groups())
        if (date,level,dof)!=(sd,sl,sn) or abs(t/86400-st)>0.000050001 or abs(dt-sdt)>0.050001:
            raise ValueError('test_log state differs from stdout')
        rounding=0.50001*10**(math.floor(math.log10(cpu))-2)
        if abs(cpu-record['seconds'])>rounding:raise ValueError('test_log time differs from stdout rounding')
        values.append(cpu)
    return values


def ordinary(log, analysis, precise=None):
    # Exclude the first timed step after every restart, including initial load.
    # Checkpoint-write steps are excluded independently. Keep ordinary remaps.
    values=precise_seconds(precise,analysis) if precise is not None else [r['seconds'] for r in analysis['step_records']]
    startup=True; saving=False; index=0; included=[]; excluded=[]
    for line in log.splitlines():
        if 'Restarting from checkpoint' in line: startup=True
        if 'Saving field' in line or 'Saving data' in line: saving=True
        if STEP.match(line):
            row=analysis['step_records'][index]; reasons=[]
            if startup and not row['checkpoint']: reasons.append('first step after startup/reload')
            if row['checkpoint']: reasons.append('checkpoint write')
            if saving: reasons.append('data saving')
            value=values[index]
            if not math.isfinite(value) or value<=0: raise ValueError('Invalid timestep duration')
            entry=dict(index=index,state=row['state'],seconds=value)
            if reasons: excluded.append(dict(entry,reasons=reasons))
            else: included.append(entry)
            # write_checkpoint reloads BEFORE printing the checkpoint step's
            # timing. That row is still the old integration transaction;
            # retain the restart exclusion for the next full step.
            startup=startup and row['checkpoint']; saving=False; index+=1
    if index!=len(analysis['step_records']) or len(included)<5: raise ValueError('Insufficient ordinary steps')
    return dict(precision_source='test_log' if precise is not None else 'stdout',median=statistics.median(r['seconds'] for r in included),included=included,excluded=excluded,
                scope='One run-level median; timesteps share run conditions and are not independent repetitions.')


def screen(a,b):
    gain=a-b; fraction=gain/a
    decision=('advance' if fraction>=POLICY['screen_gain_fraction'] else
              'reject_regression' if fraction<=-POLICY['clear_regression_fraction'] else
              'stop_no_clear_gain')
    return dict(baseline=a,candidate=b,gain_seconds=gain,gain_fraction=fraction,decision=decision,
                scope='Exploratory single pair, not a confirmed speedup.')


def confirm(batches):
    # Each batch is strictly A1 B1 B2 A2. Pair early and late adjacent runs.
    gains=[]; a=[]; b=[]
    for batch in batches:
        if len(batch)!=4 or [r[0] for r in batch]!=['baseline','candidate','candidate','baseline']:
            raise ValueError('Confirmation must be A-B-B-A')
        values=[r[1] for r in batch]
        a.extend([values[0],values[3]]); b.extend(values[1:3])
        gains.extend([values[0]-values[1],values[3]-values[2]])
    mean=statistics.mean(a)-statistics.mean(b); target=POLICY['acceptance_seconds']
    near=abs(mean-target)<=POLICY['near_threshold_seconds']
    disagree=min(gains)<0<max(gains) or min(gains)<target<=max(gains)
    repeat=near or disagree
    if repeat: decision='repeat_abba' if len(batches)<POLICY['maximum_confirmation_batches'] else 'inconclusive'
    elif min(gains)>=target: decision='meets_target'
    elif max(gains)<=0: decision='reject_regression'
    else: decision='below_target'
    return dict(decision=decision,mean_baseline_run_median=statistics.mean(a),
                mean_candidate_run_median=statistics.mean(b),gain_seconds=mean,
                paired_gains_seconds=gains,near_threshold=near,disagreement=disagree,
                independent_runs_per_version=len(a),
                scope='Mean of separate run medians; screening pairs excluded. No timestep-level significance test.')
