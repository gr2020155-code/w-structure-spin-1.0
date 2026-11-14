"""Arithmetic checks of "The Electron in the W-Structure" (sections 1, 2, 3, 4, 8, 10, 11).
Most checks here are algebra once the premises are typed in (R = ħ/(mc), S = ħ/2); they guard against
arithmetic slips, they do not derive anything. Constants: CODATA 2018 (Tiesinga et al. 2021).
Run: python3 verify_spin.py  (needs numpy, scipy)."""
import math
import numpy as np
from scipy.integrate import quad

# CODATA 2018
h = 6.62607015e-34; hbar = h / (2 * math.pi); c = 299792458.0
e = 1.602176634e-19; me = 9.1093837015e-31
mu_B = e * hbar / (2 * me)
ok = True
def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)

print("== 1. Circulation at speed c, ω = c/R, E = ħω (algebra) ==")
R = hbar / (me * c)                       # from m c^2 = ħ ω, ω = c / R
print(f"R_e = ħ/(m_e c) = {R:.6e} m")
L = me * c * R
check("angular momentum of the circulation m c R = ħ (not ħ/2: the open contradiction of section 4)", abs(L / hbar - 1) < 1e-12)
check("action per cycle ∮p dl = m c · 2πR = h", abs(me * c * 2 * math.pi * R / h - 1) < 1e-12)
mu = e * c * R / 2
check("μ = e c R / 2 = μ_B", abs(mu / mu_B - 1) < 1e-12)
g_half = 2 * me * mu / (e * hbar / 2)
g_full = 2 * me * mu / (e * hbar)
print(f" g with S = ħ/2: {g_half:.6f};  g with S = m c R = ħ: {g_full:.6f}")
check("g = 2 if the spin is ħ/2, g = 1 if it is the circulation moment ħ", abs(g_half - 2) < 1e-12 and abs(g_full - 1) < 1e-12)
R2 = hbar / (2 * me * c)                  # Schrödinger's radius, for comparison (contradicts P8)
print(f" for comparison, R = ħ/(2mc) = {R2:.3e} m would give μ = {e*c*R2/2/mu_B:.2f} μ_B; it contradicts P8 and the data μ_e ≈ μ_B")

print("\n== 1a. A circulation slower than c ==")
# P8 fixes R = ħ/(mc); μ = e v R/2 = μ_B · v/c
for frac in (1.0, 0.99, 0.9):
    v = frac * c; muv = e * v * R / 2
    print(f" v = {frac} c: μ/μ_B = {muv/mu_B:.6f}")
print(" with R fixed by P8, μ/μ_B = v/c (algebra); data: μ_e/μ_B = 1.00116 > 1 — a slower circulation gives μ < μ_B, the wrong direction")

print("\n== 2. Internal clock in motion ==")
def period(v, contracted):
    """Period of the circulation (in units of the rest period) for an orbit in a plane containing v.
    The point always moves at c. contracted: orbit contracted by γ along v."""
    gam = 1 / math.sqrt(1 - v * v)
    k = gam if contracted else 1.0
    def integrand(th):
        tx, ty = -math.sin(th) / k, math.cos(th)
        dl = math.hypot(tx, ty)
        tx, ty = tx / dl, ty / dl
        b = v * tx
        w = -b + math.sqrt(b * b + 1 - v * v)   # |w t + v x| = 1
        return dl / w
    return quad(integrand, 0, 2 * math.pi, epsabs=1e-13, epsrel=1e-13)[0] / (2 * math.pi)
print(" across the motion the speed left for the circulation is √(c² − v²), so the period is γ (algebra)")
print(" v/c    γ         in plane      in plane+contraction")
for v in (0.1, 0.5, 0.9, 0.99):
    gam = 1 / math.sqrt(1 - v * v)
    t_in = period(v, False); t_con = period(v, True)
    print(f" {v:<5}  {gam:.6f}  {t_in:.6f}     {t_con:.6f}")
    check(f"v={v}: in plane with contraction = γ", abs(t_con / gam - 1) < 1e-9)
    check(f"v={v}: in plane without contraction ≠ γ", abs(t_in / gam - 1) > 1e-3)

print("\n== 2a. Exact formula and uniqueness of the contraction ==")
# T = γ · P(γ·p, q) / (2π): the orbit is a circle contracted by p along v and by q across;
# P(A, B) is the perimeter of an ellipse with semi-axes A, B.
from scipy.special import ellipe
def perim(A, B):
    A, B = max(A, B), min(A, B)
    return 4 * A * ellipe(1 - (B / A) ** 2)
for v in (0.5, 0.9):
    gam = 1 / math.sqrt(1 - v * v)
    check(f"v={v}: formula γ·P(γ,1)/2π equals the integral without contraction",
          abs(gam * perim(gam, 1) / (2 * math.pi) - period(v, False)) < 1e-8)
# the perimeter grows strictly with a semi-axis ⇒ P(A, q) = 2πq only for A = q ⇒ p = q/γ (unique)
As = np.linspace(0.5, 2.0, 301)
P = np.array([perim(A, 1.0) for A in As])
check("the perimeter of an ellipse grows strictly with a semi-axis (unique among ellipses with unchanged transverse size)", np.all(np.diff(P) > 0))

print("\n== 3. Field energy and the size of the core ==")
eps0 = 8.8541878128e-12
Rc = 2e-20
U = e**2 / (8 * math.pi * eps0 * Rc) / e / 1e9   # GeV, Coulomb field outside Rc
print(f" Coulomb field energy outside {Rc:.0e} m: {U:.1f} GeV = {U*1e3/0.51099895:.1e} m_e c^2")
check("field energy outside 2e-20 m is 36 GeV (e²/8πε₀R)", abs(U - 36.0) < 0.1)
print(f" ħ/(m_e c) / (2e-20 m) = {R/Rc:.1e}: the circulation radius is at least this many times the core")

print("\n== 4. Limits on the step a ==")
hc_MeV_m = 197.3269804e-15
R_lep = hc_MeV_m / 1e7                    # Λ ≈ 10 TeV
print(f"LEP: ħc/Λ = {R_lep:.2e} m [estimate]; electron ≈ 3a  ⇒  a < {R_lep/3:.1e} m")
print(f"resolution time at √s = 209 GeV: ħ/E ≈ {6.582119569e-22/2.09e5:.1e} s; circulation period / this time ≈ {2*math.pi*R/c/(6.582119569e-22/2.09e5):.1e}")
print(f"g−2 (quadratic estimate): R < {hc_MeV_m/(0.51099895/math.sqrt(0.7e-12)):.2e} m")
print(f"steps on the electron circulation: N_e = 2πR_e/a > {2*math.pi*R/(R_lep/3):.2e}")
# LHAASO GRB 221009A, quadratic case: E_QG,2 > 6e-8 E_Pl (LHAASO convention includes the factor (n+1)/2 = 3/2).
# The step bound depends on how a enters: corrections (a/ƛ)², ƛ = ħc/E, give a < ħc/E_QG,2; corrections (a/λ)², λ = hc/E, give a < hc/E_QG,2.
E_Pl = 1.22089e19                          # GeV
aL = hc_MeV_m * 1e-3 / (6e-8 * E_Pl)
print(f"LHAASO [estimate]: corrections (a/ƛ)² → a < {aL:.1e} m; corrections (a/λ)² → a < {2*math.pi*aL:.1e} m")
lam_C = h / (me * c)
check("N·a = h/(m c) — the Compton wavelength", abs(2 * math.pi * R / lam_C - 1) < 1e-12)
for name, M in [("Sgr A*", 4.3e6), ("3 M☉", 3.0)]:
    r = 1476.625 * M
    print(f"{name}: (a/r)^2 < {(R_lep/3/r)**2:.1e}")

print("\n== 5. Hydrogen ground state with reduced mass ==")
mp = 1.67262192369e-27; Ry_inf = 13.605693122994
print(f" E_1 = {Ry_inf:.4f} eV × μ/m_e = {Ry_inf * mp / (mp + me):.4f} eV")
check("reduced-mass ground state 13.598 eV", abs(Ry_inf * mp / (mp + me) - 13.598) < 1e-3)
# a charge circulating on R = ħ/(m_e c) = α a₀ smears the Coulomb potential: ΔV = (R²/6)∇²V, ∇²V = 4π δ(r) (atomic units)
# with the 1s density |ψ(0)|² = 1/π: ΔE = (2/3) α² hartree; the Dirac–Darwin term of 1s is α²/2 hartree
alpha = 1 / 137.035999084; Eh = 27.211386245988; h_eV = 4.135667696e-15
smear = 2 / 3 * alpha**2; darwin = alpha**2 / 2
print(f" circulation smearing / Darwin = {smear/darwin:.4f}; if it replaced the Darwin term, extra 1S shift = {(smear-darwin)*Eh/h_eV/1e9:.0f} GHz")
check("smearing on ħ/(mc) is 4/3 of the Darwin term", abs(smear / darwin - 4 / 3) < 1e-12)

print("\n== 6. Mass of a composite: acts (energy) against closures (clock) ==")
# total electron binding energy (fit of Lunney, Pearson, Thibault 2003), eV; kinetic energy of the shells = binding (virial, 1/r)
Bel = lambda Z: 14.4381 * Z**2.39 + 1.55468e-6 * Z**5.35
u_eV = 931.49410242e6
fAl = Bel(13) / (26.98 * u_eV); fFe = Bel(26) / (55.845 * u_eV)
print(f" electron shells: Al {Bel(13)/1e3:.1f} keV ({fAl:.2e} of the mass), Fe {Bel(26)/1e3:.1f} keV ({fFe:.2e})")
print(f" if the source counted closures (clock, m/γ), active/passive mass of Fe and Al would differ by ~{fFe - fAl:.1e}")
print(" lunar laser ranging: equal to 3.9e-14 (Singh et al. 2023)")
check("counting closures for the source is excluded by lunar laser ranging", fFe - fAl > 1e3 * 3.9e-14)

print("\nRESULT:", "all checks passed" if ok else "ERRORS FOUND")
