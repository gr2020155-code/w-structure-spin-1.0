"""A moving particle core: its centre shifts one step every tick (P1) and in T ticks
follows a simple closed path (no repeated places). The core makes its exchanges relative to the
current centre. At every tick all points are partitioned into pairs of neighbours; far from the path — the medium.
Question: for which lengths T of the circulation does a solution exist for the electron and for the proton."""
import sys, itertools
from pysat.formula import WCNF, CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
L = 10
MED = "H0 H1 V0 H1 V0 V1 H0 V1".split()
sites = [(x, y) for x in range(L) for y in range(L)]
def add(s, d): return ((s[0] + d[0]) % L, (s[1] + d[1]) % L)
edges = [(s, add(s, d)) for s in sites for d in ((1, 0), (0, 1))]
eid = {e: i for i, e in enumerate(edges)}
def ek(a, b): return (a, b) if (a, b) in eid else (b, a)
def med(kind):
    out = set()
    for (x, y) in sites:
        if kind == 'H0' and x % 2 == 0: out.add(ek((x, y), add((x, y), (1, 0))))
        if kind == 'H1' and x % 2 == 1: out.add(ek((x, y), add((x, y), (1, 0))))
        if kind == 'V0' and y % 2 == 0: out.add(ek((x, y), add((x, y), (0, 1))))
        if kind == 'V1' and y % 2 == 1: out.add(ek((x, y), add((x, y), (0, 1))))
    return out
DIRS = ((1, 0), (-1, 0), (0, 1), (0, -1))
CORES = {
    'electron': [((0, 0), (-1, 0)), ((0, 0), (-1, 0)), ((0, 0), (1, 0)), ((0, 0), (1, 0))],
    'proton':   [((-1, 0), (0, 0)), ((-1, 0), (0, 0)), ((0, 0), (1, 0)), ((1, 0), (2, 0)), ((1, 0), (2, 0)), ((0, 0), (1, 0))],
}
def feasible(core, T, medium_far=True, far_r=3):
    cnf = CNF(); top = [0]
    def nv(): top[0] += 1; return top[0]
    E = [[nv() for _ in edges] for _ in range(T)]
    P = [[nv() for _ in sites] for _ in range(T)]
    sid = {s: i for i, s in enumerate(sites)}
    for t in range(T):
        for s in sites:   # exactly one pair for every point
            vs = [E[t][eid[e]] for e in edges if s in e]
            cnf.append(vs)
            for a, b in itertools.combinations(vs, 2): cnf.append([-a, -b])
        # exactly one position of the centre
        enc = CardEnc.equals(lits=P[t], bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); cnf.extend(enc.clauses)
        a_rel, b_rel = core[t % len(core)]
        for s in sites:
            p = P[t][sid[s]]
            cnf.append([-p, E[t][eid[ek(add(s, a_rel), add(s, b_rel))]]])          # core exchange
            nxt = [P[(t + 1) % T][sid[add(s, d)]] for d in DIRS]                   # shift by one step
            cnf.append([-p] + nxt)
    # simple path: every place is visited at most once per circulation
    for s in sites:
        vs = [P[t][sid[s]] for t in range(T)]
        enc = CardEnc.atmost(lits=vs, bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); cnf.extend(enc.clauses)
    # fix the start (symmetry)
    cnf.append([P[0][sid[(L // 2, L // 2)]]])
    # far from any position of the centre — the medium (soft constraints are not available in SAT: require the medium for pairs whose two points
    # are farther than far_r from the start in the Chebyshev metric — the path is kept within a window of radius far_r - 2)
    if medium_far:
        c = (L // 2, L // 2)
        def dist(p): return max(abs((p[0] - c[0] + L // 2) % L - L // 2), abs((p[1] - c[1] + L // 2) % L - L // 2))
        for t in range(T):
            M = med(MED[t % 8])
            for e in edges:
                if dist(e[0]) > far_r and dist(e[1]) > far_r:
                    cnf.append([E[t][eid[e]]] if e in M else [-E[t][eid[e]]])
    with Cadical153(bootstrap_with=cnf.clauses) as s:
        ok = s.solve()
        if not ok: return None
        m = set(x for x in s.get_model() if x > 0)
        path = [next(sites[i] for i in range(len(sites)) if P[t][i] in m) for t in range(T)]
        return path
if __name__ == "__main__":
    for name, core in CORES.items():
        for T in (4, 6, 8, 12, 16, 24):
            for r in (3, 4):
                path = feasible(core, T, far_r=r)
                if path:
                    print(f"{name}: circulation of length {T} exists (medium beyond radius {r}); path of the centre: {path}")
                    break
            else:
                print(f"{name}: circulation of length {T} — no solution")
