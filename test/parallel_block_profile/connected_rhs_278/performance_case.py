"""Offline performance decision screen from sealed Stage270/272 evidence.

No fitted runtime forecast: all scenarios expose their assumptions. Does not run
solvers, change accepted source, or use uncalibrated 160/256 phase timings.
"""
import argparse
import hashlib
import json
from pathlib import Path
import runpy
import statistics


def analyse(history):
    provenance = {}
    def read_json(relative):
        path = history / relative
        provenance[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
        return json.loads(path.read_text())
    r270 = read_json('docs/parallel-blocks-stage270-results.json')
    r272 = read_json('docs/parallel-blocks-stage272-results.json')
    r275 = read_json('docs/parallel-blocks-stage275-results.json')
    r277 = read_json('docs/parallel-blocks-stage277-results.json')
    ordinary = (2, 3, 4, 5, 6, 7, 10, 11, 12)
    checked = {}
    def rows(stage, job, label, ranks):
        review = history / f'.stage{stage}-validation/review-{job}/review'
        manifest = json.loads((review / 'files-sha256.json').read_text())
        phases_path = history / f'test/parallel_block_profile/{"cost_audit" if stage == 270 else "scaling"}_{stage}/phases.py'
        provenance[str(phases_path.relative_to(history))] = hashlib.sha256(phases_path.read_bytes()).hexdigest()
        names = runpy.run_path(str(phases_path))['NAMES']
        folder = review / 'results' / label
        files = sorted(folder.glob('cost-rank-*.txt'))
        assert len(files) == ranks, label
        result = {}
        for rank in range(ranks):
            path = folder / f'cost-rank-{rank:06d}.txt'
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            assert digest == manifest[str(path.relative_to(review))], str(path)
            checked[str(path.relative_to(history))] = digest
            lines = path.read_text().splitlines()
            assert lines[0].split() == [f'STAGE{stage}_COST_V1', str(rank), '12', str(len(names))]
            assert len(lines) == 1 + 12 * len(names)
            result[rank] = {}
            for step in range(1, 13):
                phases = {}
                for index, name in enumerate(names):
                    v = lines[1 + (step - 1) * len(names) + index].split()
                    assert len(v) == 10 and list(map(int, v[:2])) == [step, index + 1]
                    phases[name] = dict(self_s=float(v[3]), inclusive_s=float(v[5]), work=list(map(int, v[7:])))
                assert abs(sum(p['self_s'] for p in phases.values()) - phases['step']['inclusive_s']) < 1e-6
                if step in ordinary:
                    result[rank][step] = phases
        return result

    phase_names = ('shared compatibility producer', 'scalar capture preparation',
                   'scalar capture publication', 'velocity remainder preparation',
                   'scalar restriction transaction', 'RK state assembly',
                   'canonical publication', 'native scalar forward', 'migration')
    costs = {name: [] for name in phase_names}
    for label in r270['candidate_only']['per_run']:
        raw = rows(270, 55467, label, 83)
        for name in phase_names:
            metric = 'self_s' if name == 'shared compatibility producer' else 'inclusive_s'
            value = statistics.mean(raw[rank][s][name][metric] for rank in raw for s in ordinary)
            expected = r270['candidate_only']['per_run'][label]['phases'][name][metric.replace('_s', '_wall')]
            assert abs(value - expected) < 1e-12, name
            costs[name].append(value * 1000)
    phase_ms = {name: statistics.mean(values) for name, values in costs.items()}
    producer = phase_ms['shared compatibility producer']
    pool = sum(phase_ms[name] for name in ('scalar capture preparation', 'scalar capture publication', 'velocity remainder preparation'))
    baseline = 1000 * r270['comparison']['means_seconds']['stage266']
    gap = 1000 * r270['comparison']['candidate_minus_main_seconds']
    scenarios = []
    for fraction, slowdown, overhead in ((.25, 0, 10), (.5, 0, 10), (.5, .05, 10), (.5, .10, 10), (.5, .05, 20), (.75, 0, 10)):
        saved = fraction * pool - slowdown * producer - overhead
        scenarios.append(dict(removed_fraction_of_transfer_pool=fraction, producer_cost_increase_fraction=slowdown,
            new_overhead_ms=overhead, conditional_saved_ms=saved, conditional_time_reduction_percent=100 * saved / baseline,
            break_even_new_overhead_ms=fraction * pool - slowdown * producer))

    traffic = []
    for n, index in ((83, 6), (160, 7), (256, 8)):
        for b in r275['cases'][f'J5-n{n}-validation-{index}-candidate']['budgets']:
            f = b['fields']; export = b['current_tendency_export_remote_bytes_per_layer']
            resident_rho = sum(f[k]['hypothetical_upstream_remote_bytes'] for k in ('rho', 'rho_theta'))
            revised = b['hypothetical_upstream_remote_bytes_per_layer'] - resident_rho
            physics = f['physics']['Domain_remote_bytes']
            traffic.append(dict(ranks=n, ordinal=b['rank_local_ordinal'], export_bytes=export,
                residual_upstream_bytes_with_resident_rho=revised, retained_rhs_physics_bytes=physics,
                local_rhs_physics_ratio=revised/export, retained_rhs_physics_ratio=(revised+physics)/export,
                remaining_partial_byte_headroom=export-revised))

    inventories = {}
    for n, indices in ((160, (2, 7)), (256, (4, 5))):
        runs = []
        for index in indices:
            raw = rows(272, 55471, f'J5-n{n}-timing-{index}-candidate', n)
            inventory = []
            for step in ordinary:
                d = [raw[rank][step]['patch inventory']['work'][0] for rank in raw]
                b = [raw[rank][step]['patch inventory']['work'][1] for rank in raw]
                assert sum(d) == sum(b) == r272['ordinary_global_work'][str(step)]['domain_patch_samples']
                inventory.append(dict(step=step, sum_patch_stage_visits=sum(b), max_domain_patch_stage_visits=max(d),
                    max_block_patch_stage_visits=max(b), block_to_domain_max_ratio=max(b)/max(d),
                    max_to_mean_block_ratio=max(b)/statistics.mean(b)))
            runs.append(inventory)
        assert runs[0] == runs[1], n
        inventories[str(n)] = runs[0]
    count_unchanged = sum(a['max_block_patch_stage_visits'] == b['max_block_patch_stage_visits']
        for a, b in zip(inventories['160'], inventories['256']))
    inventory_ratio = (statistics.mean(x['max_block_patch_stage_visits'] for x in inventories['256']) /
                       statistics.mean(x['max_block_patch_stage_visits'] for x in inventories['160']))
    latest = r277['comparison']['control_step_seconds']
    return dict(status='conditional_design_screen_not_runtime_forecast', production_qualified=False,
        source_reports_sha256=provenance, verified_raw_files=len(checked), raw_file_hashes=checked,
        calibrated_scope='Stage270 accepted266 candidate-only observer at83; mean-rank per-step costs; not synchronized critical-path savings',
        phase_ms=phase_ms, transfer_pool_ms=pool, stage270_bare_baseline_ms=baseline, stage270_legacy_gap_ms=gap,
        zero_replacement_pool_fraction_of_baseline_percent=100*pool/baseline,
        zero_replacement_pool_fraction_of_legacy_gap_percent=100*pool/gap,
        model='saved_ms = removed_fraction * transfer_pool_ms - producer_cost_increase_fraction * producer_ms - new_overhead_ms',
        scenarios=scenarios, sampled_partial_traffic=traffic,
        traffic_scope='Rank-local ordinal observations per layer; density reuse and upstream ownership conditional; excludes producer halos and scalar pipeline, not a global epoch or whole-step total',
        inventory= inventories, identical_160_256_max_steps=count_unchanged,
        mean_max_block_256_to_160_ratio=inventory_ratio,
        inventory_scope='Exact stage-summed integer inventories; equal patch cost hypothetical; no calibrated timing or per-stage critical path',
        latest_256_extra_saved_ms_needed_to_beat_candidate160_offset=1000*(latest['256-stage266']-latest['160-stage266']),
        decision='Revise production layout and transfer retirement contract; no unconditional speedup justification; require removal/replacement measurements before launch')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = analyse(args.history.resolve())
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    summary = {k: result[k] for k in ('verified_raw_files', 'phase_ms', 'transfer_pool_ms', 'scenarios',
                                    'identical_160_256_max_steps', 'mean_max_block_256_to_160_ratio',
                                    'latest_256_extra_saved_ms_needed_to_beat_candidate160_offset')}
    print(json.dumps(summary, indent=2))
