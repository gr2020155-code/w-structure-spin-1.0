"""Isotropy in one realization: how the random direction dependence of propagation falls with distance.

isotropy_check.py shows that the ensemble average has no preferred direction. A physical arrangement of places is one
realization. Here, for single clouds without a lattice (random points equal in all directions in a ball, Delaunay neighbours,
one move per tick), we measure at distance r (in mean spacings):
  s(r) — the scatter over end points of the speed d/h (Euclidean distance per hop), relative to its mean,
         for points in a thin shell around r, after removing the mean;
  D(r) — the per-realization ratio of the mean speed in 20° cones around the diagonals to that around the axes, minus 1;
         its mean (with standard error over realizations) and its spread over realizations.
A fit ∝ r^(−β) on r = 4 … 32 gives the exponent with which these fluctuations fall. s(r) is a scatter over end points:
it mixes local randomness at each end point and the whole-number hop count h (∝ 1/r) with any dependence on direction,
so it is not itself a direction dependence; D(r) is closer to one. For first-passage distances in three dimensions the
expected asymptotic value is β = 1 − χ ≈ 0.76 (χ ≈ 0.24). Premise of the measure: places are embedded in Euclidean space
and rods measure the Euclidean distance; the speed of light is the speed of the shortest chain of moves.
"""
import sys, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.spatial import Delaunay, cKDTree

rng = np.random.default_rng(11)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300_000
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 8
AX = np.vstack([np.eye(3), -np.eye(3)])
DG = np.array([[sx, sy, sz] for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]) / np.sqrt(3)
cosT = np.cos(np.radians(20))
shells = [(4, 5), (6, 7), (8, 10), (11, 14), (15, 19), (20, 25), (26, 32)]

def cloud():
    Rb = (3 * N / (4 * np.pi)) ** (1/3)
    p = rng.normal(size=(N, 3)); p /= np.linalg.norm(p, axis=1)[:, None]
    pts = p * Rb * rng.random(N)[:, None] ** (1/3)
    s = Delaunay(pts).simplices
    e = np.unique(np.sort(np.vstack([s[:, [i, j]] for i in range(4) for j in range(i+1, 4)]), axis=1), axis=0)
    G = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(N, N)).tocsr()
    src = cKDTree(pts).query([0, 0, 0])[1]
    h = shortest_path(G, method='D', directed=False, unweighted=True, indices=src)
    v = pts - pts[src]; d = np.linalg.norm(v, axis=1)
    out = []
    for lo, hi in shells:
        sel = (d >= lo) & (d < hi) & (h > 0)
        sp = d[sel] / h[sel]
        # remove the trend inside the shell (speed grows slowly with r): residual of a linear fit in d
        A = np.vstack([np.ones(sel.sum()), d[sel]]).T
        res = sp - A @ np.linalg.lstsq(A, sp, rcond=None)[0]
        u = v[sel] / d[sel, None]
        ina = (np.abs(u @ AX.T) > cosT).any(axis=1); ind = ((u @ DG.T) > cosT).any(axis=1)
        out.append((res.std() / sp.mean(), sp[ind].mean() / sp[ina].mean() - 1))
    return np.array(out), Rb

t0 = time.time()
res = []
for i in range(NC):
    o, Rb = cloud(); res.append(o)
res = np.array(res)                       # (clouds, shells, 2)
print(f"clouds {NC}, points in each {N}, ball radius {Rb:.0f} mean spacings ({time.time()-t0:.0f} s)")
print(" r (mean spacings) | scatter s(r) | per-realization diagonals/axes − 1: mean ± standard error, spread")
rm = np.array([(lo + hi) / 2 for lo, hi in shells])
for k, (lo, hi) in enumerate(shells):
    s = res[:, k, 0]; D = res[:, k, 1]
    print(f"  {lo:>2}–{hi:<2}            {s.mean():.4f}        {D.mean():+.4f} ± {D.std(ddof=1)/np.sqrt(NC):.4f}, spread {D.std(ddof=1):.4f}")
sm = res[:, :, 0].mean(axis=0)
beta = -np.polyfit(np.log(rm), np.log(sm), 1)[0]
Dsd = res[:, :, 1].std(axis=0, ddof=1)
betaD = -np.polyfit(np.log(rm), np.log(Dsd), 1)[0]
print(f"fit s(r) ∝ r^(−β): β = {beta:.2f}")
print(f"fit spread of the per-realization diagonals/axes ratio ∝ r^(−β): β = {betaD:.2f}")
print("the range r = 4 … 32 is short; the exponent is an estimate, not a limit value")
