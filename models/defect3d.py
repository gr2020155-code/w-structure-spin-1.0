"""Three-dimensional exchange model. Cube L×L×L on a torus. At every tick all points are partitioned
into pairs of neighbours. Medium: the sequence X0 X1 Y0 X1 Y0 Y1 X0 Y1 (every point follows an
eight in its xy plane; layers along z are identical). Cores of particles of different dimension:
e — a line (O with A, B), μ — a cross in a plane (O with A, B, C, D), τ — a cross + a perpendicular
eight (O with A, B, C, D, E, F), p — a chain of 4 points.
Find the smallest rearrangement of the medium (MaxSAT) for the given core exchanges.
EXPECTATION, written before the run: costs of one order (a few rearrangements per tick),
no ratios in the hundreds or thousands; the model checks the consistency of forms, not masses."""
import sys, itertools, time, json, os, tempfile
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2
L = 6
MED = "X0 X1 Y0 X1 Y0 Y1 X0 Y1".split()
sites = [(x, y, z) for x in range(L) for y in range(L) for z in range(L)]
U = {'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1)}
def add(s, d): return tuple((s[i] + d[i]) % L for i in range(3))
edges = [(s, add(s, U[a])) for s in sites for a in 'xyz']
eid = {e: i for i, e in enumerate(edges)}
def ek(a, b): return (a, b) if (a, b) in eid else (b, a)
inc = {s: [] for s in sites}
for e in edges: inc[e[0]].append(e); inc[e[1]].append(e)
def med(kind):
    ax, par = kind[0].lower(), int(kind[1]); i = 'xyz'.index(ax)
    return {ek(s, add(s, U[ax])) for s in sites if s[i] % 2 == par}
def dist(p, c):
    return max(abs((p[i] - c[i] + L // 2) % L - L // 2) for i in range(3))
def solve(core, T, radius, c):
    MEDT = [med(MED[t % len(MED)]) for t in range(T)]
    w = WCNF()
    V = lambda t, e: 1 + t * len(edges) + eid[e]
    for t in range(T):
        for s in sites:
            vs = [V(t, e) for e in inc[s]]
            w.append(vs)
            for a, b in itertools.combinations(vs, 2): w.append([-a, -b])
        for e in edges:
            lit = V(t, e) if e in MEDT[t] else -V(t, e)
            if dist(e[0], c) > radius and dist(e[1], c) > radius: w.append([lit])
            else: w.append([lit], weight=1)
    for t in range(T):
        a, b = core[t % len(core)]
        w.append([V(t, ek(add(c, a), add(c, b)))])
    with RC2(w) as rc:
        m = rc.compute()
        return None if m is None else rc.cost
O = (0, 0, 0)
def ax(i, s): d = [0, 0, 0]; d[i] = s; return tuple(d)
A, B, C, D, E, F = ax(0, -1), ax(0, 1), ax(1, -1), ax(1, 1), ax(2, -1), ax(2, 1)
CORES = {
    'e (1D)': [(O, A), (O, A), (O, B), (O, B)],
    'μ (2D)': [(O, A), (O, A), (O, B), (O, B), (O, C), (O, C), (O, D), (O, D)],
    'τ (3D)': [(O, A), (O, A), (O, B), (O, B), (O, C), (O, C), (O, D), (O, D), (O, E), (O, E), (O, F), (O, F)],
    'p (1D)': [(A, O), (A, O), (O, B), (B, ax(0, 2)), (B, ax(0, 2)), (O, B)],
}
if __name__ == "__main__":
    T = 24
    c = (L // 2, L // 2, L // 2)
    res = {}
    for name, core in CORES.items():
        for r in (1, 2):
            t0 = time.time(); k = solve(core, T, r, c)
            if k is not None:
                print(f"{name}: radius {r}, rearrangements in {T} ticks {k}, per tick {k/T:.2f} ({time.time()-t0:.0f} s)", flush=True)
                res[name] = (r, k); break
            print(f"{name}: radius {r} — no solution ({time.time()-t0:.0f} s)", flush=True)
    json.dump(res, open(os.path.join(tempfile.gettempdir(), 'defect3d.json'), 'w'), ensure_ascii=False)
