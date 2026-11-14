"""Search for sequences of pairings in which EVERY point follows a true figure-eight:
a closed path, exactly one place visited twice (the crossing), two loops
of opposite direction (areas of opposite sign), no other repeats."""
import itertools
from medium import run
def is_eight(path):
    loop = path[:-1]
    seen = {}
    for i, q in enumerate(loop): seen.setdefault(q, []).append(i)
    reps = {q: v for q, v in seen.items() if len(v) > 1}
    if len(reps) != 1: return None
    (q, idx), = reps.items()
    if len(idx) != 2: return None
    i, j = idx
    def area(seg): return 0.5 * sum(x0*y1 - x1*y0 for (x0, y0), (x1, y1) in zip(seg, seg[1:]))
    a1 = area(path[i:j+1]); a2 = area(path[j:] + path[1:i+1])
    if a1 * a2 < 0: return (abs(a1), abs(a2), j - i, len(loop) - (j - i))
    return None
found = []
for k in range(2, 9):
    for seq in itertools.product(['H0', 'H1', 'V0', 'V1'], repeat=k):
        if seq[0] != 'H0' or any(seq[i] == seq[(i+1) % k] for i in range(k)): continue
        traj = run(seq, periods=4)
        ok = True; shapes = set()
        for s, path in traj.items():
            P = next((P for P in range(k, len(path), k) if path[P] == (0, 0)), None)
            if P is None: ok = False; break
            e = is_eight(path[:P+1])
            if e is None: ok = False; break
            shapes.add((P,) + e)
        if ok: found.append((seq, shapes))
for seq, shapes in found: print(' '.join(seq), '| path period, loop areas, steps in loops:', sorted(shapes))
print('found:', len(found))
