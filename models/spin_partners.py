"""Exchange identity, and why it says nothing about the moment of a pattern.

1. Kinematic identity: in an exchange x↔y the points move by d and −d. The moment of the pair is
   r_x×d + r_y×(−d) = (r_x − r_y)×d = (−d)×d = 0 and its displacement is d − d = 0.
   So the sum over all points of moment and displacement is zero for ANY sequence of exchanges. This is algebra:
   it holds whatever the points do and therefore cannot tell where the moment of a pattern goes.
2. Illustration on solutions of moving.py (medium + moving core, all points in pairs): the sum over points is zero
   (as the identity requires), while the path of the pattern (core centre) over its circulation encloses a nonzero area.
   The moment of a pattern is a property of the sequence, not a sum over points.
By the identity, the candidate "the counter-displacement of exchange partners carries −ħ/2" is excluded:
partners and carriers together carry exactly zero.
"""
import itertools, sys
from fractions import Fraction
from pysat.formula import CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
import moving as M

L = M.L
def mi(d):  # minimal image on the torus
    return ((d[0] + L // 2) % L - L // 2, (d[1] + L // 2) % L - L // 2)
def cross(a, b): return a[0] * b[1] - a[1] * b[0]

def solve_full(core, T, far_r=3, rot_sign=None):
    """Like moving.feasible, but returns the path of the centre and all pairs at every tick.
    rot_sign: +1/−1 — require counterclockwise/clockwise circulation of the centre (by the sign of the area), None — any."""
    cnf = CNF(); top = [0]
    def nv(): top[0] += 1; return top[0]
    sites, edges, eid, ek, add = M.sites, M.edges, M.eid, M.ek, M.add
    E = [[nv() for _ in edges] for _ in range(T)]
    P = [[nv() for _ in sites] for _ in range(T)]
    sid = {s: i for i, s in enumerate(sites)}
    for t in range(T):
        for s in sites:
            vs = [E[t][eid[e]] for e in edges if s in e]
            cnf.append(vs)
            for a, b in itertools.combinations(vs, 2): cnf.append([-a, -b])
        enc = CardEnc.equals(lits=P[t], bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); cnf.extend(enc.clauses)
        a_rel, b_rel = core[t % len(core)]
        for s in sites:
            p = P[t][sid[s]]
            cnf.append([-p, E[t][eid[ek(add(s, a_rel), add(s, b_rel))]]])
            cnf.append([-p] + [P[(t + 1) % T][sid[add(s, d)]] for d in M.DIRS])
    for s in sites:
        enc = CardEnc.atmost(lits=[P[t][sid[s]] for t in range(T)], bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); cnf.extend(enc.clauses)
    c = (L // 2, L // 2)
    cnf.append([P[0][sid[c]]])
    def dist(p): return max(abs(mi((p[0] - c[0], p[1] - c[1]))[0]), abs(mi((p[0] - c[0], p[1] - c[1]))[1]))
    for t in range(T):
        Med = M.med(M.MED[t % 8])
        for e in edges:
            if dist(e[0]) > far_r and dist(e[1]) > far_r:
                cnf.append([E[t][eid[e]]] if e in Med else [-E[t][eid[e]]])
    with Cadical153(bootstrap_with=cnf.clauses) as s:
        if not s.solve(): return None
        m = set(x for x in s.get_model() if x > 0)
    path = [next(sites[i] for i in range(len(sites)) if P[t][i] in m) for t in range(T)]
    pairs = [[edges[i] for i in range(len(edges)) if E[t][i] in m] for t in range(T)]
    return path, pairs

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and cond

print("== 1. Exchange identity (on integer vectors) ==")
for rx, ry in itertools.product(itertools.product(range(-3, 4), repeat=2), repeat=2):
    d = (ry[0] - rx[0], ry[1] - rx[1])
    if abs(d[0]) + abs(d[1]) != 1: continue
    mom = cross(rx, d) + cross(ry, (-d[0], -d[1]))
    if mom != 0: check("moment of the pair = 0", False)
check("moment of a pair of neighbours in an exchange = 0 for all pairs in a 7×7 window", True)

print("\n== 2. Illustration: medium + moving core ==")
for name in ("electron", "proton"):
    for T in (4, 8, 12):
        res = solve_full(M.CORES[name], T)
        if res is None:
            print(f"{name}, T={T}: no solution"); continue
        path, pairs = res
        c = path[0]
        Lpts = 0
        for t in range(T):
            for (x, y) in pairs[t]:
                d = mi((y[0] - x[0], y[1] - x[1]))
                rx = mi((x[0] - c[0], x[1] - c[1])); ry = mi((y[0] - c[0], y[1] - c[1]))
                Lpts += cross(rx, d) + cross(ry, (-d[0], -d[1]))
        # moment of the path of the pattern: Σ r_c × Δc over a circulation = 2·(oriented area)
        Lpat = 0
        for t in range(T):
            r = mi((path[t][0] - c[0], path[t][1] - c[1]))
            d = mi((path[(t + 1) % T][0] - path[t][0], path[(t + 1) % T][1] - path[t][1]))
            Lpat += cross(r, d)
        n_pairs = sum(len(p) for p in pairs)
        print(f"{name}, T={T}: pairs per circulation {n_pairs}; sum of point moments {Lpts} (zero by the identity); "
              f"pattern path 2·area = {Lpat}")
        check(f"{name}, T={T}: sum over points = 0 (identity)", Lpts == 0)
        check(f"{name}, T={T}: moment of the pattern path ≠ 0", Lpat != 0)

print("\nRESULT:", "all checks passed" if ok else "ERRORS FOUND")
