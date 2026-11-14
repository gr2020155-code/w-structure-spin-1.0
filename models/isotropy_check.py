"""Isotropy of propagation: one move to a neighbour per tick; speed along an axis and along a diagonal.

Three arrangements of places:
  A1) a regular cubic lattice, all moves equal;
  A2) the same lattice of places, but random moments of exchange (random time of a move along a link):
      dynamics is present, the lattice of places remains;
  B)  places form no lattice: random points equal in all directions (a Poisson
      arrangement in a ball), neighbours by Delaunay tetrahedra; one move per tick.
Measure for B: ratio of speeds in equal 20° cones around the axes and around the diagonals,
sums over many independent clouds, error by jackknife; control — the same cones
rotated at random (must give 1 for any arrangement).
"""
import sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra, shortest_path
from scipy.spatial import Delaunay, cKDTree
from scipy.spatial.transform import Rotation

rng = np.random.default_rng(7)
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and cond

print("== A1. Regular lattice, all moves equal ==")
print(f" speed along the diagonal (1,1,1) / along an axis = 1/√3 = {1/np.sqrt(3):.3f} — anisotropy of order one")

print("\n== A2. Lattice of places, random moments of exchange ==")
n = 71; c = n // 2
idx = np.arange(n**3).reshape(n, n, n)
rows = np.concatenate([np.take(idx, range(0, n - 1), axis=ax).ravel() for ax in range(3)])
cols = np.concatenate([np.take(idx, range(1, n), axis=ax).ravel() for ax in range(3)])
AX = np.vstack([np.eye(3), -np.eye(3)])
DG = np.array([[sx, sy, sz] for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]) / np.sqrt(3)
cosT = np.cos(np.radians(20))
g = np.indices((n, n, n)).reshape(3, -1).T - c
dg = np.linalg.norm(g, axis=1); selg = (dg > 24) & (dg < 33); ug = g[selg] / dg[selg, None]
def cones(u, rot=None):
    uu = u if rot is None else rot.apply(u)
    return (np.abs(uu @ AX.T) > cosT).any(axis=1), ((uu @ DG.T) > cosT).any(axis=1)
def ratio(x): return (x[:, 0].sum() / x[:, 1].sum()) / (x[:, 2].sum() / x[:, 3].sum())
def jack(a):
    a = np.array(a); r = ratio(a); jk = np.array([ratio(np.delete(a, j, 0)) for j in range(len(a))])
    return r, np.sqrt((len(a) - 1) * np.mean((jk - jk.mean()) ** 2))
A2 = {"axes and diagonals": [], "control: random rotation": []}
for seed in range(30):
    G = coo_matrix((rng.exponential(1.0, size=rows.size), (rows, cols)), shape=(n**3, n**3)).tocsr()
    t = dijkstra(G, directed=False, indices=idx[c, c, c])[selg]
    for kk, rot in (("axes and diagonals", None), ("control: random rotation", Rotation.random(random_state=rng.integers(1 << 30)))):
        ia, idd = cones(ug, rot)
        A2[kk].append((t[ia].sum(), dg[selg][ia].sum(), t[idd].sum(), dg[selg][idd].sum()))
for kk, arr in A2.items():
    r, e = jack(arr); print(f" {kk}: speed along the diagonal / along an axis = {r:.4f} ± {e:.4f} (30 realizations, r from 24 to 33)")
    if kk == "axes and diagonals": mA2, eA2 = r, e
    else: check("A2 control: the measure is unbiased", abs(r - 1) < 3 * e)
check("A2: dynamics on a lattice of places does not remove anisotropy (difference above 3 errors)", abs(mA2 - 1) > 3 * eA2)

print("\n== B. Places without a lattice: random equal points (ball), Delaunay neighbours ==")
def cloud(N=70_000):
    Rball = (3 * N / (4 * np.pi)) ** (1/3)
    p = rng.normal(size=(N, 3)); p /= np.linalg.norm(p, axis=1)[:, None]
    pts = p * Rball * rng.random(N)[:, None] ** (1/3)
    s = Delaunay(pts).simplices
    e = np.unique(np.sort(np.vstack([s[:, [i, j]] for i in range(4) for j in range(i+1, 4)]), axis=1), axis=0)
    G = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(N, N)).tocsr()
    src = cKDTree(pts).query([0, 0, 0])[1]
    h = shortest_path(G, method='D', directed=False, unweighted=True, indices=src)
    v = pts - pts[src]; d = np.linalg.norm(v, axis=1)
    sel = (d > 8) & (d < 13); u = v[sel] / d[sel, None]
    out = {}
    for name, rot in (("axes and diagonals", None), ("control: random rotation", Rotation.random(random_state=rng.integers(1 << 30)))):
        uu = u if rot is None else rot.apply(u)
        ina = (np.abs(uu @ AX.T) > cosT).any(axis=1); ind = ((uu @ DG.T) > cosT).any(axis=1)
        out[name] = (h[sel][ina].sum(), d[sel][ina].sum(), h[sel][ind].sum(), d[sel][ind].sum())
    return out, 2 * len(e) / N
NC = int(sys.argv[1]) if len(sys.argv) > 1 else 16
acc = {}
for i in range(NC):
    o, deg = cloud()
    for kk, vv in o.items(): acc.setdefault(kk, []).append(vv)
print(f" clouds {NC}, points in each 70 000, mean number of neighbours {deg:.1f}")
resB = {}
for kk, arr in acc.items():
    a = np.array(arr); r = ratio(a)
    jk = np.array([ratio(np.delete(a, j, 0)) for j in range(len(a))])
    err = np.sqrt((len(a) - 1) * np.mean((jk - jk.mean()) ** 2))
    resB[kk] = (r, err)
    print(f" {kk}: speed along the diagonal / along an axis = {r:.4f} ± {err:.4f}")
r, err = resB["axes and diagonals"]
check("B: speeds along axes and diagonals agree within 3 errors", abs(r - 1) < 3 * err)
rc, ec = resB["control: random rotation"]
check("control: the measure is unbiased (randomly rotated cones give 1)", abs(rc - 1) < 3 * ec)
print("\nRESULT:", "all checks passed" if ok else "ERRORS FOUND")
