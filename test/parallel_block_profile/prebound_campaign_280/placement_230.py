"""nk: two 128-core sockets; pin one rank per physical core across both sockets."""
from placement_200 import topology, verify


def plan(top):
    cores = {}
    for row in sorted(top['rows'], key=lambda r: r['cpu']):
        cores.setdefault((row['socket'], row['core']), row)
    sockets = sorted({s for s, _ in cores})
    if len(sockets) != 2 or any(sum(s == sock for s, _ in cores) != 128 for sock in sockets):
        raise ValueError('Require nk allocation with two complete 128-physical-core sockets')
    orders = {}
    for sock in sockets:
        groups = {}
        for (s, _), row in cores.items():
            if s == sock:
                groups.setdefault((row['node'], tuple(row['llc'])), []).append(row)
        groups = [sorted(g, key=lambda r: r['cpu']) for _, g in sorted(groups.items())]
        orders[sock] = [g[i] for i in range(max(map(len, groups))) for g in groups if i < len(g)]
    result = {}
    for n in (16, 83, 167):
        mapping = ([dict(orders[sockets[0]][i],rank=i) for i in range(n)] if n==16 else
                   [dict(orders[sockets[i % 2]][i // 2], rank=i) for i in range(n)])
        result[n] = dict(ranks=n, sockets=1 if n==16 else 2, host=top['host'], mapping=mapping,
                         socket_rank_counts=[16,0] if n==16 else [(n+1)//2,n//2],
                         policy='One hardware thread per physical core; 16 on one socket, otherwise alternate sockets; spread within LLC/NUMA groups.')
    return result
