"""Hydrogen atom, excluded model: the electron follows radial "petals" that leave the proton and return to it.

This model is kept as a record of a test that failed: section 4 shows that radial petals put far too much electron
density at the proton (21 cm line), and hydrogen_loops.py section 1 shows that ∮ p·dl = n·h gives the same levels
for loops of every eccentricity, so the levels do not single out petals. The model of the paper is hydrogen_loops.py.

Premises (named explicitly):
  P4 (W): λ = h/p — a whole number of wavelengths fits a closed path: ∮ p dr = n·h.
  The 1/r law of attraction of charges — DATA, not derived in W (W has no mechanism for it yet).
  The closed path is a petal: from the proton to the turning point and back (the premise under test).
Atomic units: ħ = m_e = e²/(4πε₀) = 1; length a₀, energy 27.2 eV; α = 1/137.036.
Comparison with measurements and with the quantum description of the 1s ground state.
"""
import math
from scipy.integrate import quad
from scipy.optimize import brentq

alpha = 1 / 137.035999084
Eh = 27.211386245988  # eV
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and cond

def action_petal(R):      # ∮ p dr along the petal 0 → R → 0 at E = −1/R
    f = lambda r: math.sqrt(max(2 * (1 / r - 1 / R), 0.0))
    return 2 * quad(f, 0, R, limit=200)[0]

print("== 1. Levels: petal 0 → R → 0, ∮ p dr = 2πn ==")
for n in (1, 2, 3, 4):
    R = brentq(lambda R: action_petal(R) - 2 * math.pi * n, 0.1, 1000)
    E = -1 / R
    print(f" n={n}: turning point R = {R:.4f} a₀, E = {E*Eh:.4f} eV, Rydberg −13.6057/n² = {-Eh/2/n**2:.4f} eV")
    check(f"n={n}: E = −1/(2n²) (Bohr–Rydberg levels)", abs(E + 0.5 / n**2) < 1e-6)
R1 = brentq(lambda R: action_petal(R) - 2 * math.pi, 0.1, 100)
check("n=1: turning point at 2a₀ — where E = V for the 1s state", abs(R1 - 2) < 1e-6)

print("\n== 2. If the closed path is a passage through the proton (−R → R → −R) instead of a petal ==")
R = brentq(lambda R: 2 * action_petal(R) - 2 * math.pi, 0.01, 100)
print(f" n=1: E = {-1/R*Eh:.3f} eV — four times deeper than the data (−13.606 eV)")
check("within this model, a passage through the proton counted as one cycle gives four times the energy", abs(-1/R + 2.0) < 1e-6)

print("\n== 3. Time averages for n=1 (R = 2) against 1s ==")
R = 2.0
v = lambda r: math.sqrt(2 * (1 / r - 1 / R))
Thalf = quad(lambda r: 1 / v(r), 0, R, limit=200)[0]
avg = lambda g: quad(lambda r: g(r) / v(r), 0, R, limit=400)[0] / Thalf
rows = [("⟨r⟩", avg(lambda r: r), 1.5), ("⟨1/r⟩", avg(lambda r: 1 / r), 1.0), ("⟨r²⟩", avg(lambda r: r * r), 3.0),
        ("⟨p²⟩", avg(lambda r: v(r) ** 2), 1.0)]
for name, cl, qm in rows:
    print(f" {name:6s} petal {cl:.4f}   1s {qm:.4f}")
check("⟨r⟩ = 1.5 a₀ equals 1s", abs(rows[0][1] - 1.5) < 1e-4)
check("⟨1/r⟩ and ⟨p²⟩ agree (energy)", abs(rows[1][1] - 1) < 1e-4 and abs(rows[3][1] - 1) < 1e-4)
print(f" ⟨r²⟩: {rows[2][1]:.3f} against 3 — a difference of 1/6 (testable through the diamagnetic susceptibility and polarizability of atomic H)")

print("\n== 4. Density at the proton (hyperfine splitting, 21 cm line) ==")
qm0 = 1 / math.pi
for rc, name in ((alpha, "the electron circulation radius R_e = α·a₀"), (0.84e-15 / 5.29177e-11, "the proton radius")):
    frac = quad(lambda r: 1 / v(r), 0, rc, limit=200)[0] / Thalf
    rho = frac / (4 / 3 * math.pi * rc**3)
    print(f" mean density within {name}: {rho:.3g} / a₀³; for 1s |ψ(0)|² = {qm0:.3f} / a₀³; ratio {rho/qm0:.3g}")
check("a radial petal gives a density at the proton far above 1s — radial petals are excluded",
      quad(lambda r: 1 / v(r), 0, alpha, limit=200)[0] / Thalf / (4/3*math.pi*alpha**3) > 10 * qm0)

print("\n== 5. ⟨p⁴⟩ (relativistic correction of the 1s level) ==")
for rc in (alpha, 0.1):
    p4 = quad(lambda r: v(r) ** 4 / v(r), rc, R, limit=400)[0] / Thalf
    print(f" excluding r < {rc:.4f} a₀: ⟨p⁴⟩ = {p4:.1f}; for 1s ⟨p⁴⟩ = 5 (confirmed by the fine structure and the 1S–2S line)")
print(" as r → 0 the ⟨p⁴⟩ of a petal diverges: ∫ dr/(v r²) ~ ∫ r^(−3/2) dr")

print("\n== 6. The nucleus does not matter: binding energy = μα²c²/2, μ the reduced mass ==")
me = 0.51099895e6
for name, M in (("hydrogen (proton)", 938.27208816e6), ("deuterium", 1875.61294257e6), ("positronium (positron)", me), ("muonium (μ⁺)", 105.6583755e6)):
    mu = me * M / (me + M)
    print(f" {name:22s} {mu*alpha**2/2:.4f} eV")
print(" measured: H 13.598 eV; Ps 6.80 eV — binding depends only on the charges and μ, not on the structure of the nucleus")

print("\nRESULT:", "all checks passed" if ok else "ERRORS FOUND")
