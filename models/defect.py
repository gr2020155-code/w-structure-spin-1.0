"""A particle as a local deviation from the medium of eights.
At every tick all points are partitioned into pairs of neighbours (a perfect matching).
Far from the particle — the medium (H0 H1 V0 H1 V0 V1 H0 V1). At the particle — the given core exchanges.
Find the solution with the fewest differences from the medium (MaxSAT), then count how many
points of the medium deviate from their eights over a cycle ("reach")."""
import sys, json, os, tempfile
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2
L = int(sys.argv[1]) if len(sys.argv) > 1 else 10
MED = "H0 H1 V0 H1 V0 V1 H0 V1".split()
T = 24
MEDT_PERIOD = 8
sites = [(x, y) for x in range(L) for y in range(L)]
def nb(s, d): return ((s[0] + (d == 'x')) % L, (s[1] + (d == 'y')) % L)
edges = [(s, nb(s, d)) for s in sites for d in 'xy']
eid = {}
for e in edges: eid[e] = len(eid)
def ekey(a, b):
    return (a, b) if (a, b) in eid else (b, a)
def medium_edges(kind):
    out = set()
    for (x, y) in sites:
        if kind == 'H0' and x % 2 == 0: out.add(ekey((x, y), ((x+1) % L, y)))
        if kind == 'H1' and x % 2 == 1: out.add(ekey((x, y), ((x+1) % L, y)))
        if kind == 'V0' and y % 2 == 0: out.add(ekey((x, y), (x, (y+1) % L)))
        if kind == 'V1' and y % 2 == 1: out.add(ekey((x, y), (x, (y+1) % L)))
    return out
MEDT = [medium_edges(MED[t % 8]) for t in range(T)]
def var(t, e): return 1 + t * len(edges) + eid[e]
def inc(s): return [e for e in edges if s in e]

def solve(core, radius, c0):
    """core: list of (t, (a, b)) — required core exchanges; radius: beyond this distance — the medium."""
    w = WCNF()
    for t in range(T):
        for s in sites:
            vs = [var(t, e) for e in inc(s)]
            w.append(vs)
            for i in range(len(vs)):
                for j in range(i + 1, len(vs)): w.append([-vs[i], -vs[j]])
        for e in edges:
            far = all(max(abs((p[0] - c0[0] + L//2) % L - L//2), abs((p[1] - c0[1] + L//2) % L - L//2)) > radius for p in e)
            v = var(t, e); m = e in MEDT[t]
            if far: w.append([v] if m else [-v])
            else: w.append([v] if m else [-v], weight=1)   # soft: penalty for a difference from the medium
    for t, (a, b) in core: w.append([var(t, ekey(a, b))])
    with RC2(w) as rc:
        m = rc.compute()
        if m is None: return None, None
        cost = rc.cost
    on = set(x for x in m if x > 0)
    tilings = [[e for e in edges if var(t, e) in on] for t in range(T)]
    return tilings, cost

def trajectories(tilings, periods=6):
    occ = {s: s for s in sites}; path = {s: [s] for s in sites}
    for k in range(periods * T):
        for a, b in tilings[k % T]:
            pa, pb = occ[a], occ[b]; occ[a], occ[b] = pb, pa
        for s, p in occ.items(): path[p].append(s)
    return path

def medium_paths(periods=6):
    return trajectories([sorted(m) for m in MEDT], periods)

if __name__ == "__main__":
    c = (L // 2, L // 2)
    O, A, B = c, ((c[0] - 1) % L, c[1]), ((c[0] + 1) % L, c[1])
    # electron: O↔A for two ticks, O↔B for two ticks, repeat (period 4; the window T = 24 is a multiple of 4, 6 and 8)
    pat_e = [(O, A), (O, A), (O, B), (O, B)]
    core_e = [(t, pat_e[t % 4]) for t in range(T)]
    P1, P2 = c, ((c[0] + 1) % L, c[1]); Ap, Cp = ((c[0] - 1) % L, c[1]), ((c[0] + 2) % L, c[1])
    # proton: A↔P1 ×2, P1↔P2 ×1, P2↔C ×2, P1↔P2 ×1 (period 6; the window T = 24 is a multiple of 4, 6 and 8)
    pat_p = [(Ap, P1), (Ap, P1), (P1, P2), (P2, Cp), (P2, Cp), (P1, P2)]
    core_p6 = [(t, pat_p[t % 6]) for t in range(T)]
    out = {}
    for name, core in (("electron", core_e), ("proton", core_p6)):
        for r in range(1, L // 2):
            til, cost = solve(core, r, c)
            if til is not None:
                med = medium_paths(); par = trajectories(til)
                moved = [s for s in sites if par[s][:T * 6 + 1] != med[s][:T * 6 + 1]]
                print(f"{name}: found at radius {r}; differences from the medium in 24 ticks: {cost}; points whose paths changed: {len(moved)}")
                out[name] = {"radius": r, "cost": cost, "moved": len(moved)}
                break
        else:
            print(f"{name}: no solution up to radius {L//2 - 1}")
    json.dump(out, open(os.path.join(tempfile.gettempdir(), f"defect_L{L}.json"), "w"), ensure_ascii=False)
