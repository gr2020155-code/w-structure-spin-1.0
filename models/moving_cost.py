"""Selection of the circulation: for each length T find the path with the fewest differences from the medium
(MaxSAT) and compare the cost per tick for the electron and the proton."""
import itertools, time
from pysat.formula import WCNF
from pysat.card import CardEnc, EncType
from pysat.examples.rc2 import RC2
import moving as M
def cost(core, T, far_r=3):
    w = WCNF(); top = [0]
    def nv(): top[0] += 1; return top[0]
    E = [[nv() for _ in M.edges] for _ in range(T)]
    P = [[nv() for _ in M.sites] for _ in range(T)]
    sid = {s: i for i, s in enumerate(M.sites)}
    c = (M.L // 2, M.L // 2)
    def dist(p): return max(abs((p[0]-c[0]+M.L//2) % M.L - M.L//2), abs((p[1]-c[1]+M.L//2) % M.L - M.L//2))
    for t in range(T):
        for s in M.sites:
            vs = [E[t][M.eid[e]] for e in M.edges if s in e]
            w.append(vs)
            for a, b in itertools.combinations(vs, 2): w.append([-a, -b])
        enc = CardEnc.equals(lits=P[t], bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); [w.append(cl) for cl in enc.clauses]
        a_rel, b_rel = core[t % len(core)]
        for s in M.sites:
            p = P[t][sid[s]]
            w.append([-p, E[t][M.eid[M.ek(M.add(s, a_rel), M.add(s, b_rel))]]])
            w.append([-p] + [P[(t+1) % T][sid[M.add(s, d)]] for d in M.DIRS])
        Med = M.med(M.MED[t % 8])
        for e in M.edges:
            v = E[t][M.eid[e]]; lit = v if e in Med else -v
            if dist(e[0]) > far_r and dist(e[1]) > far_r: w.append([lit])
            else: w.append([lit], weight=1)
    for s in M.sites:
        enc = CardEnc.atmost(lits=[P[t][sid[s]] for t in range(T)], bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); [w.append(cl) for cl in enc.clauses]
    w.append([P[0][sid[c]]])
    with RC2(w) as rc:
        m = rc.compute()
        return None if m is None else rc.cost
for name, core in M.CORES.items():
    for T in (4, 8, 12, 24):
        t0 = time.time(); k = cost(core, T)
        print(f"{name}: T={T}: differences from the medium {k}, per tick {k/T:.2f}  ({time.time()-t0:.0f} s)", flush=True)
