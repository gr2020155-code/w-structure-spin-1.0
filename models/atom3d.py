"""Hydrogen atom in a three-dimensional exchange model: the electron goes round the proton to the other side.

Model premises: at every tick every point takes part in exactly one exchange with a neighbour; far away — the medium (eights in
xy planes as in defect3d.py; this medium is not isotropic — a simplification); the proton core A–P₁–P₂–C stays at the centre and
makes its 6 exchanges (A↔P₁, A↔P₁, P₁↔P₂, P₂↔C, P₂↔C, P₁↔P₂); the electron core O with neighbours ±x makes its
4 exchanges (O↔A, O↔A, O↔B, O↔B) relative to its centre, and the centre shifts one step per tick (P1).
Points of the electron core never coincide with points of the proton core: the electron does not pass through the proton.
Window of 24 ticks (a multiple of 4, 6 and 8). At tick 0 the electron is on one side of the proton (beyond end C),
at tick 12 on the other side (beyond end A), at tick 24 back at the start.
Variants:
  A) free path;
  B) the way back repeats the way out in reverse order (the moment of the first half is cancelled by the second —
     a degenerate eight, there and back);
  C) forced: a straight path through the proton (a check that the model indeed forbids it).
"""
import itertools, sys, time
from pysat.formula import CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153

L = 12
T = 24
MED = "X0 X1 Y0 X1 Y0 Y1 X0 Y1".split()
sites = [(x, y, z) for x in range(L) for y in range(L) for z in range(L)]
sid = {s: i for i, s in enumerate(sites)}
def add(s, d): return tuple((s[k] + d[k]) % L for k in range(3))
AX = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
edges = [(s, add(s, d)) for s in sites for d in AX]
eid = {e: i for i, e in enumerate(edges)}
def ek(a, b): return (a, b) if (a, b) in eid else (b, a)
inc = {s: [] for s in sites}
for e in edges: inc[e[0]].append(e); inc[e[1]].append(e)
def med(kind):
    out = set()
    for s in sites:
        x, y, z = s
        if kind == 'X0' and x % 2 == 0: out.add(ek(s, add(s, (1, 0, 0))))
        if kind == 'X1' and x % 2 == 1: out.add(ek(s, add(s, (1, 0, 0))))
        if kind == 'Y0' and y % 2 == 0: out.add(ek(s, add(s, (0, 1, 0))))
        if kind == 'Y1' and y % 2 == 1: out.add(ek(s, add(s, (0, 1, 0))))
    return out
C0 = (L // 2, L // 2, L // 2)
def rel(d): return add(C0, d)
P_PTS = [rel((-1, 0, 0)), rel((0, 0, 0)), rel((1, 0, 0)), rel((2, 0, 0))]           # A, P₁, P₂, C
A_, P1, P2, C_ = P_PTS
P_CORE = [(A_, P1), (A_, P1), (P1, P2), (P2, C_), (P2, C_), (P1, P2)]
E_CORE = [((0, 0, 0), (-1, 0, 0)), ((0, 0, 0), (-1, 0, 0)), ((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (1, 0, 0))]
E_OFFS = [(0, 0, 0), (-1, 0, 0), (1, 0, 0)]
DIRS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
START = rel((4, 0, 0))     # electron points x = 3, 4, 5 — beyond end C (x = 2)
OTHER = rel((-4, 0, 0))    # electron points x = −5, −4, −3 — beyond end A (x = −1); displacement 8 is even, like the 12 steps (the lattice is bipartite)
def cheb(s):
    return max(abs((s[k] - C0[k] + L // 2) % L - L // 2) for k in range(3))
FAR = 6                    # beyond — the medium

def build(variant):
    cnf = CNF(); top = [0]
    def nv(): top[0] += 1; return top[0]
    E = [[nv() for _ in edges] for _ in range(T)]
    P = [[nv() for _ in sites] for _ in range(T)]
    for t in range(T):
        for s in sites:
            vs = [E[t][eid[e]] for e in inc[s]]
            cnf.append(vs)
            for a, b in itertools.combinations(vs, 2): cnf.append([-a, -b])
        enc = CardEnc.equals(lits=P[t], bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); cnf.extend(enc.clauses)
        a, b = P_CORE[t % 6]
        cnf.append([E[t][eid[ek(a, b)]]])                                  # proton core exchange
        Medt = med(MED[t % 8])
        for e in edges:
            if cheb(e[0]) >= FAR and cheb(e[1]) >= FAR:
                cnf.append([E[t][eid[e]]] if e in Medt else [-E[t][eid[e]]])
        ea, eb = E_CORE[t % 4]
        for s in sites:
            p = P[t][sid[s]]
            pts = [add(s, o) for o in E_OFFS]
            if any(q in P_PTS for q in pts) or any(cheb(q) >= FAR for q in pts):
                cnf.append([-p]); continue                                   # not through the proton, not in the medium
            cnf.append([-p, E[t][eid[ek(add(s, ea), add(s, eb))]]])           # electron core exchange
            cnf.append([-p] + [P[(t + 1) % T][sid[add(s, d)]] for d in DIRS])  # one step per tick
    cnf.append([P[0][sid[START]]]); cnf.append([P[T // 2][sid[OTHER]]])
    if variant == 'B':
        for t in range(1, T // 2):
            for s in sites:
                cnf.append([-P[t][sid[s]], P[T - t][sid[s]]])
    if variant == 'C':   # straight path through the proton: the centre on the x axis at every tick
        for t in range(T):
            for s in sites:
                if s[1] != C0[1] or s[2] != C0[2]: cnf.append([-P[t][sid[s]]])
    return cnf, E, P

def solve(variant):
    t0 = time.time()
    cnf, E, P = build(variant)
    with Cadical153(bootstrap_with=cnf.clauses) as s:
        ok = s.solve()
        m = set(x for x in s.get_model() if x > 0) if ok else None
    dt = time.time() - t0
    if not ok: return None, dt
    path = [next(sites[i] for i in range(len(sites)) if P[t][i] in m) for t in range(T)]
    pairs = [[edges[i] for i in range(len(edges)) if E[t][i] in m] for t in range(T)]
    return (path, pairs), dt

def relc(s): return tuple((s[k] - C0[k] + L // 2) % L - L // 2 for k in range(3))
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

ok_all = True
def check(name, cond):
    global ok_all
    print(("PASS " if cond else "FAIL ") + name); ok_all = ok_all and cond

for variant, title in (('A', "free path"), ('B', "the way back is the way out in reverse"), ('C', "straight through the proton")):
    res, dt = solve(variant)
    print(f"\n== Variant {variant}: {title} ({dt:.0f} s) ==")
    if res is None:
        print(" no solution")
        if variant == 'C': check("a path through the proton is impossible", True)
        else: check(f"variant {variant}: a solution exists", False)
        continue
    path, pairs = res
    check(f"variant {variant}: a solution exists", variant != 'C')
    # checks of the solution
    allpts = all(sum(1 for e in pairs[t] if s in e) == 1 for t in range(T) for s in sites)
    check(" every point in exactly one exchange at every tick", allpts)
    noover = all(not any(add(path[t], o) in P_PTS for o in E_OFFS) for t in range(T))
    check(" electron points never coincided with proton points", noover)
    print(" path of the electron centre (relative to P₁):", [relc(p) for p in path])
    Lm = [0, 0, 0]; half = [0, 0, 0]
    for t in range(T):
        r = relc(path[t]); d = relc(path[(t + 1) % T]); d = tuple(d[k] - r[k] for k in range(3))
        c = cross(r, d)
        Lm = [Lm[k] + c[k] for k in range(3)]
        if t < T // 2: half = [half[k] + c[k] for k in range(3)]
    print(f" moment of the path over the first half {tuple(half)}, over the whole cycle {tuple(Lm)}")
    changed = sum(1 for t in range(T) for e in pairs[t] if e not in med(MED[t % 8]))
    print(f" exchanges differing from the medium in 24 ticks: {changed}")

print("\nRESULT:", "all checks passed" if ok_all else "ERRORS FOUND")
