"""Which two-valued signs of a particle are invariant under rotations and reflection, and which of them are odd
under loop reversal (candidates for charge).

A particle in three-dimensional space (the 24 rotations of a cube serve as a finite check; the analytic argument for all
rotations is given at the end):
  v — direction of motion along the circulation (polar vector), the "front" end of the chain is along v;
  c — axis of the circulation (axial vector, c ⟂ v);
  n — normal of the loop at the front end (axial vector, n ⟂ v: an exchange loop lies in a plane containing the line).
Forms:
  electron  A–O–B: front loop +n, back loop −n;
  proton    A–P₁–P₂–C: front +n, middle −n, back +n;
  neutrino  two loops in one point, in time order: first n₁, then −n₁ (ordered pair, no front end, n₁ any direction);
  photon    a point whose side turns around v: sense h = ±1, no circulation.
Operations: rotations; reflection P (polar → −, axial unchanged); loop reversal C (n → −n, everything else unchanged);
reversal of motion T (v → −v, c → −c, sides of loops reversed, front and back ends swap).
A charge candidate is a function of the state with values ±1 that is invariant under rotations and P and odd under C.
"""
import itertools
import numpy as np

E = [np.array(e) for e in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))]
def key(a): return tuple(int(x) for x in a)
ROT = []
for perm in itertools.permutations(range(3)):
    for signs in itertools.product((1, -1), repeat=3):
        M = np.zeros((3, 3), int)
        for i, j in enumerate(perm): M[i, j] = signs[i]
        if round(np.linalg.det(M)) == 1: ROT.append(M)
assert len(ROT) == 24
perp = lambda v: [e for e in E if e @ v == 0]
rot3 = lambda R, s: tuple(key(R @ np.array(x)) for x in s)
P = lambda s: (key(-np.array(s[0])), s[1], s[2])
C = lambda s: (s[0], s[1], key(-np.array(s[2])))

def orbits(states, act):
    left, out = set(states), []
    while left:
        s = left.pop(); orb = {s}
        for R in ROT: orb.add(act(R, s))
        left -= orb; out.append(frozenset(orb))
    return out

def charge_like(states):
    """All ±1 functions on rotation classes, invariant under P, odd under C, taking both values."""
    orbs = orbits(states, rot3); idx = {s: i for i, o in enumerate(orbs) for s in o}
    rep = [next(iter(o)) for o in orbs]
    found = []
    for f in itertools.product((1, -1), repeat=len(orbs)):
        if all(f[idx[P(rep[i])]] == f[i] and f[idx[C(rep[i])]] == -f[i] for i in range(len(orbs))):
            found.append(f)
    return orbs, idx, rep, found

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and cond

def describe(s):
    v, c, n = map(np.array, s)
    if n @ c: return f"n·c = {int(n @ c):+d}"
    if n @ v: return f"n·v = {int(n @ v):+d}"
    return f"n·(c×v) = {int(n @ np.cross(c, v)):+d}"

print("== 1. Electron and proton (loops carried by the particle) ==")
chain = [(key(v), key(c), key(n)) for v in E for c in perp(v) for n in perp(v)]
orbs, idx, rep, found = charge_like(chain)
print(f"states {len(chain)}, classes under rotations {len(orbs)}:")
for i, s in enumerate(rep):
    print(f"  class {i}: {describe(s):<16} P → class {idx[P(s)]}, C → class {idx[C(s)]}")
print(f"charge-like functions on all states: {len(found)}")
check("on all states no ±1 function is invariant under rotations and P and odd under C "
      "(the n ⟂ c classes are swapped by P and by C together)", len(found) == 0)
sub = [s for s in chain if np.array(s[2]) @ np.array(s[1]) != 0]
orbs_s, idx_s, rep_s, found_s = charge_like(sub)
print(f"restricted to loops in the plane of circulation (n ∥ c): classes {len(orbs_s)}, charge-like functions {len(found_s)}")
check("with n ∥ c there is exactly one charge-like sign (up to overall sign): n·c",
      len(found_s) == 2 and all(f[i] * f[i] == 1 for f in found_s for i in range(len(orbs_s)))
      and all(f[i] == f[0] * int(np.array(rep_s[i][2]) @ np.array(rep_s[i][1])) * int(np.array(rep_s[0][2]) @ np.array(rep_s[0][1]))
              for f in found_s for i in range(len(orbs_s))))
nperp = [s for s in chain if np.array(s[2]) @ np.array(s[1]) == 0]
orbs_p, _, _, found_p = charge_like(nperp)
check("with n ⟂ c there is no charge-like sign (the same form would be neutral)", len(found_p) == 0)
binary = []
for f in itertools.product((1, -1), repeat=len(orbs)):
    if len(set(f)) == 2 and all(f[idx[P(rep[i])]] == f[i] for i in range(len(orbs))):
        binary.append(f)
print(f"two-valued signs invariant under rotations and P (any behaviour under C): {len(binary) // 2} up to overall sign")
print("  among them: n·c on the n ∥ c classes, and 'n ∥ c versus n ⟂ c' (even under C — not a charge)")

def T_of(form):
    def T(s):
        v, c, n = map(np.array, s)
        n2 = n if form == 'electron' else -n     # front ← back; side reversed
        return (key(-v), key(-c), key(n2))
    return T
for form in ('electron', 'proton'):
    T = T_of(form)
    flips = all(idx_s[T(rep_s[i])] != i for i in range(len(orbs_s)))
    print(f"{form}: under T the sign n·c {'changes' if flips else 'does not change'}")

print("\n== 2. Spin flip: the circulation axis reverses ==")
flip_c = lambda s: (s[0], key(-np.array(s[1])), s[2])
flip_cn = lambda s: (s[0], key(-np.array(s[1])), key(-np.array(s[2])))
nc = lambda s: int(np.array(s[2]) @ np.array(s[1]))
check("if only c reverses, n·c changes sign (the charge would flip with the spin)", all(nc(flip_c(s)) == -nc(s) for s in sub))
check("if the loops are carried by the circulation (n reverses with c), n·c is unchanged", all(nc(flip_cn(s)) == nc(s) for s in sub))

print("\n== 3. Neutrino: two loops in one point, in time order (n₁ first, then −n₁) ==")
nu = [(key(v), key(c), key(n)) for v in E for c in perp(v) for n in E]
orbs_n, idx_n, rep_n, found_n = charge_like(nu)
print(f"states {len(nu)}, classes {len(orbs_n)}: " + ", ".join(describe(s) for s in rep_n))
nu_par = [s for s in nu if np.array(s[2]) @ np.array(s[1]) != 0]
nu_v = [s for s in nu if np.array(s[2]) @ np.array(s[0]) != 0]
check("if the first loop lies in the plane of circulation (n₁ ∥ c), the loop order is a charge-like sign — "
      "the neutrino would be charged", len(charge_like(nu_par)[3]) == 2)
o_v, i_v, r_v, f_v = charge_like(nu_v)
pswap = all(i_v[P(r_v[i])] != i for i in range(len(o_v)))
check("if the loops lie across the direction of motion (n₁ ∥ v), there is no charge-like sign; the loop order is "
      "n₁·v = ±1, which P reverses (a helicity)", len(f_v) == 0 and pswap and len(o_v) == 2)

print("\n== 4. Photon: the side turns around v, sense h ==")
ph = [(key(v), h) for v in E for h in (1, -1)]
o = orbits(ph, lambda R, s: (key(R @ np.array(s[0])), s[1]))
check("photon: two classes, helicity h = ±1; no circulation, so no sign of the form n·c", len(o) == 2)

print("\n== 5. A side fixed by a rule n = f(v) alone ==")
inv = 0
for r in itertools.product(*[range(4) for _ in E]):
    nrule = {key(v): perp(v)[r[i]] for i, v in enumerate(E)}
    st = [(key(v), key(c)) for v in E for c in perp(v)]
    s_of = lambda s: int(nrule[s[0]] @ np.array(s[1]))
    inv += all(s_of((key(R @ np.array(s[0])), key(R @ np.array(s[1])))) == s_of(s) for s in st for R in ROT)
check(f"no rule n = f(v) gives a rotation-invariant n·c ({inv} of 4096; trivially: a rotation by 180° about v "
      "sends c to −c and leaves f(v))", inv == 0)

print("\nAnalytic argument for all rotations (SO(3)) and reflection: on states with n ∥ c write n = s·c, s = ±1.")
print("SO(3) acts transitively on pairs (v, c) of orthogonal unit vectors, so an invariant can depend only on s;")
print("P leaves the axial vectors n, c unchanged, so s is P-even; C sends s → −s. Hence s = n·c is the only charge-like sign there.")
print("\nRESULT:", "all checks passed" if ok else "ERRORS FOUND")
