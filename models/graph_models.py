"""Exchange models without a lattice: places are random points equal in all directions in a ball, Delaunay neighbours.

Checks what was checked on lattices before:
  1) the medium: every point in exactly one exchange at every tick; does a medium exist
     in which every point follows an eight (two exchanges with one neighbour, two with another);
  2) static cores of the electron (A–O–B, exchanges O↔A, O↔A, O↔B, O↔B) and of the proton
     (A–P₁–P₂–C, exchanges A↔P₁, A↔P₁, P₁↔P₂, P₂↔C, P₂↔C, P₁↔P₂) next to each other in such a medium;
  3) the atom: the proton core at the centre, the electron core moves to a neighbouring point every tick,
     goes round the proton to the other side and back; points of the cores never coincide.
Beyond radius 6 from the centre (ball of radius 7.8) — the medium (eights); inside, the solver finds any consistent rearrangement.
Parity: a Delaunay graph is not bipartite (it has triangles), so there is no lattice parity restriction.
"""
import itertools, sys, time
import numpy as np
import networkx as nx
from scipy.spatial import Delaunay
from pysat.formula import CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and cond

# --- arrangement: N points in a ball, density 1, Delaunay neighbours (without long edges at the boundary)
N = 2000
Rb = (3 * N / (4 * np.pi)) ** (1/3)
p = rng.normal(size=(N, 3)); p /= np.linalg.norm(p, axis=1)[:, None]
P = p * Rb * rng.random(N)[:, None] ** (1/3)
s = Delaunay(P).simplices
E = {tuple(sorted((int(a), int(b)))) for t in s for a, b in itertools.combinations(t, 2)}
E = [e for e in E if np.linalg.norm(P[e[0]] - P[e[1]]) < 2.2]
G = nx.Graph(); G.add_nodes_from(range(N)); G.add_edges_from(E)
G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
nodes = sorted(G.nodes)
if len(nodes) % 2: nodes = nodes[:-1]; G = G.subgraph(nodes).copy()
edges = [tuple(sorted(e)) for e in G.edges]
eid = {e: i for i, e in enumerate(edges)}
def ek(a, b): return eid[(a, b) if a < b else (b, a)]
nbr = {v: list(G.neighbors(v)) for v in G}
print(f"points {G.number_of_nodes()}, edges {len(edges)}, mean number of neighbours {2*len(edges)/G.number_of_nodes():.1f}")
odd = any(len(c) % 2 for c in nx.cycle_basis(G)[:200])
print(f"odd cycles present: {odd} (the graph is not bipartite)")

print("\n== 1. Medium: two disjoint perfect partitions into pairs, eight M1 M1 M2 M2 ==")
M1 = nx.max_weight_matching(G, maxcardinality=True)
G2 = G.copy(); G2.remove_edges_from(M1)
M2 = nx.max_weight_matching(G2, maxcardinality=True)
n = G.number_of_nodes()
check("partition 1 is perfect (every point in a pair)", 2 * len(M1) == n)
check("partition 2 is perfect and shares no edges with partition 1", 2 * len(M2) == n)
MED = [M1, M1, M2, M2]
pos = {v: v for v in G}
for M in MED:
    sw = {}
    for a, b in M: sw[a] = b; sw[b] = a
    pos = {v: sw[pos[v]] for v in G}
check("after 4 ticks every point is back; each path is two loops (an eight)", all(pos[v] == v for v in G))
MEDk = [{tuple(sorted(e)) for e in M} for M in MED]

# --- centre, cores
c0 = min(G, key=lambda v: np.linalg.norm(P[v]))
def dir_nb(v, d, excl=()):
    """the neighbour of v closest to the direction d"""
    cand = [u for u in nbr[v] if u not in excl]
    return max(cand, key=lambda u: (P[u] - P[v]) @ d / np.linalg.norm(P[u] - P[v]))
X = np.array([1.0, 0, 0])
P1 = c0; P2 = dir_nb(P1, X); A_ = dir_nb(P1, -X, (P2,)); C_ = dir_nb(P2, X, (P1, A_))
PCH = [A_, P1, P2, C_]
PCORE = [(A_, P1), (A_, P1), (P1, P2), (P2, C_), (P2, C_), (P1, P2)]
def ecore(o):
    a = dir_nb(o, -X); b = dir_nb(o, X, (a,))
    return a, b
ang = lambda a, o, b: np.degrees(np.arccos((P[a]-P[o]) @ (P[b]-P[o]) / np.linalg.norm(P[a]-P[o]) / np.linalg.norm(P[b]-P[o])))
print(f"\nproton: chain angles A–P₁–P₂ {ang(A_,P1,P2):.0f}°, P₁–P₂–C {ang(P1,P2,C_):.0f}°")
dist0 = {v: np.linalg.norm(P[v] - P[c0]) for v in G}

def solve(T, far, extra):
    cnf = CNF(); top = [len(edges) * T]
    Ev = lambda t, i: t * len(edges) + i + 1
    for t in range(T):
        for v in G:
            lits = [Ev(t, ek(v, u)) for u in nbr[v]]
            enc = CardEnc.equals(lits=lits, bound=1, top_id=top[0], encoding=EncType.seqcounter)
            top[0] = max(top[0], enc.nv); cnf.extend(enc.clauses)
        for e in edges:
            if dist0[e[0]] > far and dist0[e[1]] > far:
                cnf.append([Ev(t, eid[e])] if e in MEDk[t % 4] else [-Ev(t, eid[e])])
    extra(cnf, Ev, top)
    with Cadical153(bootstrap_with=cnf.clauses) as sv:
        r = sv.solve(); m = set(x for x in sv.get_model() if x > 0) if r else None
    return m, Ev

print("\n== 2. Static cores of the electron and the proton next to each other, 12 ticks ==")
oe = min((v for v in G if 2.8 < dist0[v] < 4.0 and v not in PCH), key=lambda v: -P[v][1])
ae, be = ecore(oe)
print(f"electron: angle A–O–B {ang(ae, oe, be):.0f}°, distance to the proton {dist0[oe]:.1f} steps")
ECORE = [(oe, ae), (oe, ae), (oe, be), (oe, be)]
def static_cores(cnf, Ev, top):
    for t in range(12):
        cnf.append([Ev(t, ek(*PCORE[t % 6]))]); cnf.append([Ev(t, ek(*ECORE[t % 4]))])
t0 = time.time(); m, Ev = solve(12, 6.0, static_cores)
check(f"cores of the electron and the proton exist in the medium ({time.time()-t0:.0f} s)", m is not None)
if m:
    allok = all(sum(Ev(t, ek(v, u)) in m for u in nbr[v]) == 1 for t in range(12) for v in G)
    check(" every point in exactly one exchange at every tick", allok)
    diff = sum(1 for t in range(12) for e in edges if Ev(t, eid[e]) in m and e not in MEDk[t % 4])
    print(f" pairs differing from the medium: {diff} of {12*n//2} in 12 ticks; the medium is fixed beyond radius 6")

print("\n== 3. Atom: the electron goes round the proton to the other side and back, 24 ticks ==")
T = 24
side_p = min((v for v in G if 3.0 < dist0[v] < 4.5), key=lambda v: -P[v][0])
side_m = min((v for v in G if 3.0 < dist0[v] < 4.5), key=lambda v: P[v][0])
CAND = [v for v in G if dist0[v] < 5.5 and v not in PCH and len(nbr[v]) >= 3]
EC = {}
for v in CAND:
    a, b = ecore(v)
    if a in PCH or b in PCH or ang(a, v, b) < 120: continue
    EC[v] = (a, b)
CAND = [v for v in CAND if v in EC]
print(f"places for the electron centre: {len(CAND)}; start x = {P[side_p][0]:+.1f}, other side x = {P[side_m][0]:+.1f}")
def atom(cnf, Ev, top):
    Pv = {}
    for t in range(T):
        for v in CAND: top[0] += 1; Pv[t, v] = top[0]
    for t in range(T):
        enc = CardEnc.equals(lits=[Pv[t, v] for v in CAND], bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); cnf.extend(enc.clauses)
        cnf.append([Ev(t, ek(*PCORE[t % 6]))])
        for v in CAND:
            a, b = EC[v]
            e = (v, a) if t % 4 < 2 else (v, b)
            cnf.append([-Pv[t, v], Ev(t, ek(*e))])
            nxt = [Pv[(t + 1) % T, u] for u in nbr[v] if u in EC]
            cnf.append([-Pv[t, v]] + nxt)
    cnf.append([Pv[0, side_p]]); cnf.append([Pv[T // 2, side_m]])
    atom.Pv = Pv
t0 = time.time(); m, Ev = solve(T, 6.0, atom)
check(f"a path round the proton exists ({time.time()-t0:.0f} s)", m is not None)
if m:
    path = [next(v for v in CAND if atom.Pv[t, v] in m) for t in range(T)]
    allok = all(sum(Ev(t, ek(v, u)) in m for u in nbr[v]) == 1 for t in range(T) for v in G)
    check(" every point in exactly one exchange at every tick", allok)
    nov = all(not ({path[t], *EC[path[t]]} & set(PCH)) for t in range(T))
    check(" points of the electron core never coincided with points of the proton", nov)
    xs = [P[v][0] for v in path]
    print(" x of the electron centre by tick:", " ".join(f"{x:+.1f}" for x in xs))

print("\nRESULT:", "all checks passed" if ok else "ERRORS FOUND")
