"""Collective proof of absent Domain trend storage on the production path."""
import re


def storage_evidence(log,ranks,kind,oracles):
    counts=[tuple(map(int,m)) for m in re.findall(r'Domain tendency storage ranks = (\d+) / (\d+)',log)]
    absent=log.count('Domain tendency storage: absent (resident block dynamics)')
    retained=log.count('Domain tendency storage: legacy/oracle')
    if kind=='baseline':
        if counts or absent or retained:raise ValueError('Baseline unexpectedly contains Stage230 storage census')
        return dict(instrumentation='not present in committed baseline')
    expected=[(ranks if oracles else 0,ranks)]*2
    if counts!=expected or (absent,retained)!=((0,2) if oracles else (2,0)):
        raise ValueError(f'Wrong storage policy: {counts}, absent={absent}, retained={retained}')
    return dict(storage_ranks=counts,absent=absent,retained=retained)
