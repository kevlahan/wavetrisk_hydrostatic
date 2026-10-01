"""Separate run medians and bounded paired comparison; timesteps are not replicates."""
import math
import re
import statistics
from analyze import STEP


def ordinary(log, analysis, precise):
    rows = [line.split() for line in precise.splitlines() if line.strip()]
    records = analysis['step_records']
    if len(rows) != len(records):
        raise ValueError('Solver high-precision timing row count differs')
    values = []
    for fields, record in zip(rows, records):
        if len(fields) != 6:
            raise ValueError('Malformed high-precision timing')
        date, t, dt, level, dof, cpu = map(float, fields)
        if not all(math.isfinite(v) for v in (date, t, dt, level, dof, cpu)) or cpu <= 0:
            raise ValueError('Nonfinite timing')
        state = re.fullmatch(r'\s*(\d+)\s+([\d.eE+-]+) d dt =\s*([\d.eE+-]+) s Jmax =\s*(\d+) dof =\s*(\d+)', record['state'])
        if state is None:
            raise ValueError('Unknown stdout state')
        sd, st, sdt, sl, sn = map(float, state.groups())
        if (date, level, dof) != (sd, sl, sn) or abs(t / 86400 - st) > 0.000050001 or abs(dt - sdt) > 0.050001:
            raise ValueError('Timing row state differs from stdout')
        rounding = 0.50001 * 10 ** (math.floor(math.log10(cpu)) - 2)
        if abs(cpu - record['seconds']) > rounding:
            raise ValueError('Timing row duration differs from stdout')
        values.append(cpu)
    startup, saving, index = True, False, 0
    included, excluded = [], []
    for line in log.splitlines():
        if 'Restarting from checkpoint' in line:
            startup = True
        if 'Saving field' in line or 'Saving data' in line:
            saving = True
        if STEP.match(line):
            record = records[index]
            reasons = []
            if startup and not record['checkpoint']:
                reasons.append('first step after startup/reload')
            if record['checkpoint']:
                reasons.append('checkpoint write')
            if saving:
                reasons.append('data saving')
            row = dict(index=index, seconds=values[index], state=record['state'])
            (excluded if reasons else included).append(dict(row, reasons=reasons) if reasons else row)
            startup = startup and record['checkpoint']
            saving = False
            index += 1
    if index != len(records) or len(included) < 5 or [r['index'] for r in excluded] != [0, 7, 8]:
        raise ValueError('Unexpected J4/J5 timing window')
    return dict(median=statistics.median(r['seconds'] for r in included), included=included,
                excluded=excluded, scope='One run median; steps within a run are not independent.')


def screen(a, b):
    change = b/a - 1
    return dict(baseline=a, candidate=b, regression_fraction=change,
                decision='review_regression' if change >= 0.03 else 'advance',
                scope='Exploratory single pair, not a confirmed performance result.')


def confirm(batches):
    baseline, candidate, pairs = [], [], []
    for batch in batches:
        if [row[0] for row in batch] != ['baseline', 'candidate', 'candidate', 'baseline']:
            raise ValueError('Confirmation order must be A-B-B-A')
        a1, b1, b2, a2 = [row[1] for row in batch]
        baseline += [a1, a2]
        candidate += [b1, b2]
        pairs += [b1/a1 - 1, b2/a2 - 1]
    change = statistics.mean(candidate)/statistics.mean(baseline) - 1
    uncertain = abs(change-0.03) <= 0.005 or min(pairs) < 0.03 <= max(pairs)
    decision = ('repeat_abba' if len(batches) < 2 else 'inconclusive') if uncertain else (
        'review_regression' if change >= 0.03 else 'within_regression_guard')
    return dict(decision=decision, baseline_run_median_mean=statistics.mean(baseline),
                candidate_run_median_mean=statistics.mean(candidate), gain_seconds=statistics.mean(baseline)-statistics.mean(candidate),
                speedup_ratio=statistics.mean(baseline)/statistics.mean(candidate), regression_fraction=change,
                paired_regression_fractions=pairs, independent_runs_per_version=len(baseline),
                scope='Separate full-run medians; screens excluded; no timestep-level significance claim.')
