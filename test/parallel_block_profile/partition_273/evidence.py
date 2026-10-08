"""Validate publication proofs and the actual installation-call counter contract."""
import re


def publication_evidence(log, kind, oracle):
    headers = log.count('Stage265 direct inverse SOL return = on')
    tiled = log.count('Stage265 tiled canonical publication = on')
    deferred = log.count('Stage266 provisional first installation = deferred; oracle = retained')
    rows = [(final == 'T', int(patches), int(values)) for final, patches, values in
            re.findall(r'Stage265 exact inverse return: final=([TF]) patches/values=\s*(\d+)\s+(\d+)', log)]
    installs = [(final == 'T', int(count)) for final, count in
                re.findall(r'Stage266 inverse installations: final=([TF]) count=(\d+)', log)]
    if kind == 'main':
        if headers or tiled or deferred or rows or installs:
            raise ValueError('Block markers unexpectedly present in main')
        return dict(mode='main_reference')
    if headers != 1 or tiled != 1:
        raise ValueError('Missing/duplicate tiled publication marker')
    if kind == 'stage265t':
        if oracle or deferred or rows or installs:
            raise ValueError('Unexpected optimization/oracle markers in accepted baseline')
        return dict(mode='tiled_baseline', provisional_installations=2, accepted_installations=1)
    if kind not in ('candidate', 'checked') or deferred != 1:
        raise ValueError('Missing/duplicate deferred installation marker')
    order = [False, False, False, True] * 12
    if oracle:
        if len(rows) != 48 or any(p <= 0 or v <= 0 for _, p, v in rows):
            raise ValueError('Incomplete direct-return transaction proof')
        if [final for final, _, _ in rows] != order:
            raise ValueError('Direct-return RK4 transaction order differs')
        if installs != [(final, 1 if final else 2) for final in order]:
            raise ValueError('Oracle must retain both provisional installations')
    elif rows or installs:
        raise ValueError('Unexpected validation cost in production run')
    return dict(mode='deferred_first_installation', exact_transactions=len(rows),
                accepted_transactions=sum(f for f, _, _ in rows),
                coverage=[dict(patches=p, values=v) for _, p, v in rows],
                oracle_installations=sum(count for _, count in installs),
                production_provisional_installations=1, production_accepted_installations=1,
                production_full_SOL_returns_per_inverse=0,
                limitation='Actual installation calls counted on every rank and asserted at each Stage203 exit; per-transaction logging occurs only in the oracle.')
