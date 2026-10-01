"""Verified, non-overlapping physical-core maps for the nk node."""
from placement_200 import topology, verify


def plan(top, ranks):
    if ranks not in (16, 41, 83, 167):
        raise ValueError('Unqualified rank count')
    cores = {}
    for row in sorted(top['rows'], key=lambda r: r['cpu']):
        cores.setdefault((row['socket'], row['core']), row)
    sockets = sorted({s for s, _ in cores})
    if len(sockets) != 2 or any(sum(s == sock for s, _ in cores) != 128 for sock in sockets):
        raise ValueError('Require two complete 128-core nk sockets')
    orders = {}
    for sock in sockets:
        groups = {}
        for (s, _), row in cores.items():
            if s == sock:
                groups.setdefault((row['node'], tuple(row['llc'])), []).append(row)
        groups = [sorted(g, key=lambda r: r['cpu']) for _, g in sorted(groups.items())]
        orders[sock] = [g[i] for i in range(max(map(len, groups))) for g in groups if i < len(g)]
    chosen = ([orders[sockets[0]][i] for i in range(ranks)] if ranks == 16 else
              [orders[sockets[i % 2]][i // 2] for i in range(ranks)])
    return dict(ranks=ranks, host=top['host'], mapping=[dict(row, rank=i) for i, row in enumerate(chosen)],
                socket_rank_counts=[ranks, 0] if ranks == 16 else [(ranks+1)//2, ranks//2],
                policy='One rank per physical core; 16 uses one socket, larger cases alternate sockets.')
