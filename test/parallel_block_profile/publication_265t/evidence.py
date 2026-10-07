"""Validate direct-return evidence without relaxing any science gate."""
import re


def publication_evidence(log, kind, oracle):
    headers = log.count('Stage265 direct inverse SOL return = on')
    tiled = log.count('Stage265 tiled canonical publication = on')
    rows = [(final == 'T', int(patches), int(values)) for final, patches, values in
            re.findall(r'Stage265 exact inverse return: final=([TF]) patches/values=\s*(\d+)\s+(\d+)', log)]
    if kind not in ('candidate', 'checked'):
        if headers or tiled or rows:
            raise ValueError('Publication markers unexpectedly present in reference')
        return dict(mode='reference')
    if headers != 1 or tiled != 1:
        raise ValueError('Missing/duplicate direct publication marker')
    if oracle:
        if len(rows) != 48 or sum(final for final, _, _ in rows) != 12 or any(p <= 0 or v <= 0 for _, p, v in rows):
            raise ValueError('Incomplete direct-return transaction proof')
        if [final for final, _, _ in rows] != [False, False, False, True] * 12:
            raise ValueError('Direct-return RK4 transaction order differs')
    elif rows:
        raise ValueError('Unexpected validation cost in production run')
    return dict(mode='tiled_direct', exact_transactions=len(rows), accepted_transactions=sum(f for f, _, _ in rows),
                coverage=[dict(patches=p, values=v) for _, p, v in rows],
                production_full_SOL_returns_per_inverse=0,
                limitation='Removed-return count enforced by the production runtime assertion; oracle proofs are logged.')
