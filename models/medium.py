"""Two-dimensional exchange model. At every tick every point takes part in exactly one exchange
with a neighbour (a partition of the lattice into pairs = a domino tiling). Enumerate periodic
sequences of four "straight" tilings and classify the trajectories of points."""
import itertools, math
L = 8
def tiling(kind):
    pairs = []
    for x in range(L):
        for y in range(L):
            if kind == 'H0' and x % 2 == 0: pairs.append(((x, y), ((x+1) % L, y)))
            if kind == 'H1' and x % 2 == 1: pairs.append(((x, y), ((x+1) % L, y)))
            if kind == 'V0' and y % 2 == 0: pairs.append(((x, y), (x, (y+1) % L)))
            if kind == 'V1' and y % 2 == 1: pairs.append(((x, y), (x, (y+1) % L)))
    return pairs
T = {k: tiling(k) for k in ('H0', 'H1', 'V0', 'V1')}

def run(seq, periods=None):
    sites = [(x, y) for x in range(L) for y in range(L)]
    pos = {s: s for s in sites}            # point (by starting place) -> place
    occ = {s: s for s in sites}            # place -> point
    traj = {s: [(0, 0)] for s in sites}    # unwrapped displacements
    steps = len(seq) * (periods or 2 * L)
    for t in range(steps):
        for a, b in T[seq[t % len(seq)]]:
            pa, pb = occ[a], occ[b]
            occ[a], occ[b] = pb, pa
            for p, frm, to in ((pa, a, b), (pb, b, a)):
                dx = (to[0] - frm[0] + L // 2) % L - L // 2
                dy = (to[1] - frm[1] + L // 2) % L - L // 2
                x0, y0 = traj[p][-1]
                traj[p].append((x0 + dx, y0 + dy))
    return traj

def classify(path, k):
    # period of a point's motion: the smallest multiple of k after which the displacement is (0,0) — a closed path
    for P in range(k, len(path), k):
        if path[P] == (0, 0):
            loop = path[:P + 1]
            area = 0.5 * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(loop, loop[1:]))
            visits = {}
            for q in loop[:-1]: visits[q] = visits.get(q, 0) + 1
            rep = sum(1 for v in visits.values() if v > 1)
            xs = [q[0] for q in loop]; ys = [q[1] for q in loop]
            return ('closed', P, round(area, 2), rep, (max(xs)-min(xs), max(ys)-min(ys)))
    d = path[-1]
    return ('escapes', None, None, None, d)

if __name__ == "__main__":
    results = {}
    for k in range(2, 9):
        for seq in itertools.product(['H0', 'H1', 'V0', 'V1'], repeat=k):
            if any(seq[i] == seq[(i + 1) % k] for i in range(k)):   # no repetition in a row (otherwise the pair stands still)
                continue
            if seq[0] != 'H0':                                     # a cyclic shift is the same order
                continue
            traj = run(seq)
            kinds = {}
            for s, path in traj.items():
                c = classify(path, k)
                key = (c[0], c[1], c[2], c[3], c[4])
                kinds[key] = kinds.get(key, 0) + 1
            results[seq] = kinds
    
    # summary: which sequences give every point a closed path with zero area (an eight) and a self-crossing
    for seq, kinds in results.items():
        keys = list(kinds)
        if all(k[0] == 'closed' for k in keys):
            tag = []
            for kk, n in kinds.items():
                tag.append(f"{n} points: period {kk[1]}, area {kk[2]}, repeated places {kk[3]}, extent {kk[4]}")
            if any(kk[2] == 0 and kk[3] > 0 for kk in keys):
                print('EIGHT?', ' '.join(seq)); [print('   ', t) for t in tag]
    print('---- sequences checked:', len(results))
    closed = [s for s, k in results.items() if all(x[0] == 'closed' for x in k)]
    print('all points closed:', len(closed))
    for s in closed[:12]:
        print(' '.join(s), dict(list(results[s].items())[:3]))
    