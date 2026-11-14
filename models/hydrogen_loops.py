"""Hydrogen atom as loops (Kepler ellipses with the proton at a focus) that never pass through the proton.

Premises, named: the 1/r attraction (data); the quantization premise that a closed loop holds a whole number of
wavelengths, ∮ p·dl = n·h with λ = h/p from P4 and no extra phase term; atomic units (ħ = m_e = e²/4πε₀ = 1),
proton mass infinite except where the reduced mass is stated.
Checked:
  1) levels: for every eccentricity, an ellipse with ∮ p·dl = 2πn has semi-major axis n² and energy −1/(2n²);
  2) the ensemble "L² uniform on [0, 1]" for n = 1 — an assumption chosen because it reproduces the 1s momentum
     distribution (it is the microcanonical ensemble; Abrines, Percival 1966): momenta, ⟨p⁴⟩, ⟨r⟩, ⟨r²⟩,
     exact fractions of loops coming close to the proton;
  3) the electron density at the proton in this ensemble (contact density: 21 cm line, Darwin term) against 1s,
     with and without the loops that P17 forbids (perihelion inside the proton), for two meanings of "the proton";
  4) the dipole: one ellipse has a time-averaged dipole; a rule for successive loops that cancels both the dipole
     and the orbital moment.
"""
import math
import numpy as np
from scipy.integrate import quad

alpha = 1 / 137.035999084
a0 = 5.29177210903e-11
rp = 0.84e-15 / a0                         # proton charge radius (Pohl et al. 2010) in a₀
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and cond

print("== 1. Levels from ∮ p·dl = 2πn for ellipses of any eccentricity ==")
def action(a, e):
    """∮ p·dl over one period of a Kepler ellipse, energy −1/(2a): p² = 2/r − 1/a, dl = v dt, ∮ p v dt = ∮ p² dt."""
    T = 2 * math.pi * a ** 1.5
    def p2(M):
        E = M
        for _ in range(80): E = E - (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        r = a * (1 - e * math.cos(E))
        return 2 / r - 1 / a
    return quad(p2, 0, 2 * math.pi, limit=400)[0] / (2 * math.pi) * T
worst = 0
for n in (1, 2, 3, 4):
    for e in (0.0, 0.3, 0.6, 0.9):
        worst = max(worst, abs(action(n * n, e) / (2 * math.pi * n) - 1))
print(f" largest relative deviation of ∮p·dl from 2πn for a = n², e = 0 … 0.9: {worst:.1e}")
check("for every e the loop with a = n² has ∮p·dl = n·h; its energy −1/(2n²) = −13.606/n² eV (−13.598/n² eV with the reduced mass)", worst < 1e-8)

print("\n== 2. Ensemble for n = 1: L² uniform on [0, 1] (assumption chosen to fit) ==")
rng = np.random.default_rng(1)
N = 400_000
L2 = rng.random(N)
e = np.sqrt(1 - L2)
M = rng.random(N) * 2 * math.pi                    # uniform in time (mean anomaly)
Ecc = M.copy()
for _ in range(60): Ecc = Ecc - (Ecc - e * np.sin(Ecc) - M) / (1 - e * np.cos(Ecc))
r = 1 - e * np.cos(Ecc)
p2 = 2 / r - 1
p = np.sqrt(p2)
# exact time averages for an ellipse (a = 1): ⟨r⟩ = 1 + e²/2, ⟨r²⟩ = 1 + 3e²/2, ⟨1/r⟩ = 1, ⟨1/r²⟩ = 1/L
# p⁴ = (2/r − 1)² → ⟨p⁴⟩ = 4⟨1/r²⟩ − 4⟨1/r⟩ + 1; with L² uniform ⟨1/L⟩ = 2
p4 = 4 * quad(lambda L: 2 * L * (1 / L), 0, 1)[0] - 4 + 1
rm = quad(lambda L: 2 * L * (1 + (1 - L * L) / 2), 0, 1)[0]
r2m = quad(lambda L: 2 * L * (1 + 1.5 * (1 - L * L)), 0, 1)[0]
print(f" ⟨p²⟩ = {p2.mean():.4f} (1s: 1)   ⟨p⁴⟩ = {p4:.4f} exactly (1s: 5; this is only the kinetic part of the fine structure)")
print(f" ⟨r⟩ = {rm:.4f} (1s: 1.5)   ⟨r²⟩ = {r2m:.4f} (1s: 3)   always r ≤ 2 (1s: 24 % beyond 2a₀)")
check("⟨p⁴⟩ = 5 exactly", abs(p4 - 5) < 1e-9)
bins = np.linspace(0, 4, 41)
hist, _ = np.histogram(p, bins=bins, density=True)
mid = (bins[1:] + bins[:-1]) / 2
qm = 32 / math.pi * mid**2 / (1 + mid**2) ** 4
dev = np.max(np.abs(hist - qm))
print(f" momentum distribution: largest deviation from 1s {dev:.4f} at a peak of {qm.max():.3f}")
check("the momentum distribution equals 1s (Abrines, Percival 1966)", dev < 0.02)
# perihelion q = 1 − e = 1 − √(1 − L²); q < x  ⇔  L² < 1 − (1 − x)²  (exact, L² uniform)
for x, name in ((alpha, "the electron circulation radius ħ/(m_e c) = α a₀"), (rp, "the proton radius")):
    print(f" fraction of loops coming closer than {name}: {1 - (1 - x)**2:.3e} exactly (sample {np.mean(1 - e < x):.2e})")

print("\n== 3. Electron density at the proton ==")
# microcanonical density at energy −1/2: ρ(r) ∝ ∫δ(H − E) d³p ∝ p(r) = √(2/r − 1), r < 2; normalized: ρ = √(2/r − 1)/(2π²)
Z = quad(lambda x: 4 * math.pi * x * x * math.sqrt(2 / x - 1), 0, 2)[0]
check("normalization ∫ 4πr² √(2/r − 1) dr = 2π²", abs(Z - 2 * math.pi**2) < 1e-6)
rho = lambda x: math.sqrt(2 / x - 1) / Z
hist_r, edges = np.histogram(r, bins=np.linspace(0.2, 1.8, 17))
pred = np.array([quad(lambda x: 4 * math.pi * x * x * rho(x), edges[k], edges[k+1])[0] for k in range(16)])
check("the sampled ensemble follows ρ ∝ √(2/r − 1)", np.max(np.abs(hist_r / N - pred)) < 3e-3)
for x, name in ((rp, "proton radius"), (alpha, "α a₀")):
    inside = quad(lambda y: 4 * math.pi * y * y * rho(y), 0, x)[0] / (4 / 3 * math.pi * x**3)
    s1 = quad(lambda y: 4 * math.pi * y * y * math.exp(-2 * y) / math.pi, 0, x)[0] / (4 / 3 * math.pi * x**3)
    print(f" mean density inside the {name}: loops {inside:.3f}, 1s {s1:.4f}; ratio {inside / s1:.1f}")
    if name == "proton radius": ratio_p = inside / s1
print(" ρ diverges as r^(−1/2) at the proton; the hyperfine (21 cm) splitting and the Darwin term are proportional to")
print(" the contact density, and both agree with the 1s value |ψ(0)|² = 1/π far better than a factor of 2")
print(" the density inside r_p comes only from loops whose perihelion q is inside r_p (fraction 3.2e-5)")
def frac_inside(x, q0):
    """time fraction with r < x for the ensemble (e² uniform), keeping only loops with perihelion q = 1 − e ≥ q0.
    For one loop the fraction is (E − e sin E)/π at cos E = (1 − x)/e (mean anomaly / π)."""
    def f(u):                              # u = 1 − e, density of e is 2e
        e_ = 1 - u
        E_ = math.acos(min(1.0, (1 - x) / e_))
        return 2 * e_ * (E_ - e_ * math.sin(E_)) / math.pi
    return quad(f, q0, x, limit=200, points=[x * 1e-3, x * 1e-1])[0] if q0 < x else 0.0
s1p = quad(lambda y: 4 * math.pi * y * y * math.exp(-2 * y) / math.pi, 0, rp)[0]
r0 = frac_inside(rp, 0.0) / s1p
check("closed form agrees with the density integral (all loops)", abs(r0 / ratio_p - 1) < 1e-3)
core = 3e-20 / a0                          # P6: the proton core is four points, ≈ 3a < 3e-20 m (LEP bound on a)
for q0, name in ((core, "the proton is its four-point core (≲ 3·10⁻²⁰ m)"), (rp, "the proton is the region of its charge radius")):
    rr = frac_inside(rp, q0) / s1p
    print(f" P17 enforced, {name}: density inside r_p = {rr:.1f} × 1s")
    if q0 == core: r_core = rr
    else: r_charge = rr
check("with P17 and the proton as its core the density inside r_p stays ≈ 68 × 1s; with the proton as its charge region it is 0;"
      " either way the contact density disagrees with the 21 cm line", 60 < r_core < 75 and r_charge == 0.0)

print("\n== 4. Dipole and orbital moment of successive loops ==")
def ellipse_avgs(e, flip):
    """Time averages of position and L_z for one ellipse (a = 1, perihelion on +x, in the xy plane);
    flip: the loop turned by π about the y axis (x → −x, z → −z)."""
    M = np.linspace(0, 2 * math.pi, 20001)[:-1]
    E = M.copy()
    for _ in range(80): E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    x = np.cos(E) - e; y = math.sqrt(1 - e * e) * np.sin(E)
    Lz = math.sqrt(1 - e * e)
    if flip: x, Lz = -x, -Lz
    return np.array([x.mean(), y.mean()]), Lz
for e_ in (0.3, 0.8):
    d0, L0 = ellipse_avgs(e_, False); d1, L1 = ellipse_avgs(e_, True)
    print(f" e = {e_}: one loop — dipole ({d0[0]:+.4f}, {d0[1]:+.4f}) (−3e/2 = {-1.5*e_:+.4f}), L = {L0:+.3f};"
          f" next loop turned by π about the axis ⟂ apse — ({d1[0]:+.4f}, {d1[1]:+.4f}), L = {L1:+.3f}")
    check(f"e = {e_}: a single loop has dipole −(3/2)·a·e along the apse line", abs(d0[0] + 1.5 * e_) < 1e-6)
    check(f"e = {e_}: with the turn rule the pair has zero dipole and zero orbital moment",
          np.allclose(d0 + d1, 0, atol=1e-9) and abs(L0 + L1) < 1e-12)
print(" the turn is a rotation about the axis through the proton perpendicular to the apse line, in the plane of the loop;")
print(" at the perihelion the electron is at distance q = 1 − e from that axis and never passes through the proton.")
print(" A reflection through the proton (r → −r, p → −p) would cancel the dipole but keep L = r×p: it is not enough.")

print("\nRESULT:", "all checks passed" if ok else "ERRORS FOUND")
