# The Electron in the W-Structure: Spin, Magnetic Moment, Charge and Motion

**Author:** Hrachya Danielyan. **Version:** 1.0, November 2026. **License:** CC BY 4.0.

## Abstract

The W-Structure describes space as points that never rest: every point exchanges places with a neighbour at every tick, and one exchange is one step and one tick. A particle is a stable pattern of exchanges that moves along its own closed path, the particle circulation. This paper works out what the author's positions imply for the electron and its neighbours, and marks for every statement whether it follows, is assumed, is computed or fails.

A circulation at c with radius ħ/(mc) gives the magnetic moment μ_B. If the observable spin is the orientation moment ħ/2 of the rule "orientation returns after two cycles", then g = 2. But the circulation itself carries angular momentum ħ, and this contradiction is open. If the loops of a charged particle lie in the plane of its circulation, charge is the only two-valued sign that is invariant under rotations and reflection and odd under loop reversal. The loops of a neutrino cannot lie in the plane of its circulation, or the order of its two loops would be a charge. Places without a lattice show no preferred direction on average. In a single realization the fluctuations of the propagation speed fall as r^(−0.73…0.80) over 4–32 mean spacings, as expected for first passage in three dimensions; a bound on the step from isotropy needs a direction-correlated measure, which is open. Lorentz transformations follow from the causal order only under a stated premise. What is computed is that a clock moving in its own plane keeps the period γ only if it is contracted by γ. In hydrogen, Kepler loops with ∮p·dl = nh give the levels, and their microcanonical ensemble gives the 1s momentum distribution. Its electron density at the proton disagrees with the 21 cm line: 68 times the 1s value if the proton is its four-point core, zero if the electron may not enter the proton's charge region. Mass is defined by closures per own tick of a particle; the mass of a composite system is not yet defined. The cycle lengths that fix mass ratios are not derived. Programs in this repository reproduce every computed number. Open problems are listed in section 12.

**Labels.**
- [author] — a position of the author (an axiom of the theory).
- [definition] — a definition used in this paper.
- [assumption] — a working assumption, named, not derived.
- [derived] — follows from the items named; [derived, algebra] — follows by algebra once the premises are written down.
- [computed] — computed by a program in this repository.
- [estimate] — approximate, up to a factor of order one.
- [candidate] — a proposed step, not proven.
- [data] — a measurement.
- [open] — not solved.
- [contradiction] — two statements of the theory or of the theory and data that do not agree; not resolved.

## 1. Positions, definitions and assumptions

### Positions

- **P1** [author]. One indivisible change is one step and one tick. Length is the number of steps, time is the number of ticks. The speed c is one step per tick.
- **P2** [author]. A link is a stable sequence of acts of change (movements of points). The arrangement is dynamic: a state is a place in a sequence, not a snapshot. No point is ever at rest.
- **P3** [author]. Local self-closure is an exchange of places of neighbouring points. An exchange is an indivisible act: one step, one tick. Its side — the side on which the two points pass each other — is a property of the act. The path of an exchange (a loop) has length one step; there is no metric inside a step.
- **P4** [author]. E = ħω and λ = h/p.
- **P5** [author]. At every step the budget of change is shared between translation and internal change in quadrature. Hence E² = (mc²)² + (pc)², and E = mc² at rest.
- **P6** [author]. Points have nonzero size, one step. The electron is three points on a line, A–O–B, with two loops of opposite orientation (exchanges O↔A and O↔B). The proton is four points on a line, A–P₁–P₂–C, with three loops: A–P₁ (+), P₁–P₂ (−), P₂–C (+). Both are one-dimensional objects in three-dimensional space. The muon is a cross (two dimensions); the tau is a cross with a perpendicular eight (three dimensions).
- **P7** [author]. Rule 2:1: the orientation of a particle returns after two cycles.
- **P8** [author]. Three levels of motion: (1) the exchange of neighbours; (2) the particle circulation — the path of its pattern, of radius ħ/(mc); (3) the global eight that every point of the whole performs.
- **P9** [author]. The closure element is the element common to the eights of a particle: the point O for the electron, the exchange P₁↔P₂ for the proton.
- **P10** [author]. A point closes on itself by one loop. A neutrino closes on itself by an eight, the smallest possible self-closure. The central point O of the electron is such an eight-point, and the whole electron is an eight: one eight on two levels.
- **P11** [author]. A closed system as a whole does not change (conservation laws).
- **P12** [author]. Indistinguishable things cannot carry different weight.
- **P13** [author]. Charge is the direction of a loop. Everything is counted from the direction of the global eight.
- **P14** [author]. The photon is a point without particle self-closure, moving one step per tick; each of its steps is an exchange with a point of the medium.
- **P15** [author]. There are no probabilities: a transition to the next configuration happens with probability 1. Locally the previous state can be computed but not restored.
- **P16** [author]. The arrangement of places forms no lattice.
- **P17** [author]. The hydrogen atom is seven points and five loops (3 + 4 points, 2 + 3 loops). The electron cannot pass through the proton; it changes sides of the proton by exchange of places, with the directions of links taken into account. The bound atom is a system of its own, not a sum of parts. The element common to the electron and the proton in the atom is the link between them, not a point or a loop.
- **P18** [author]. Mass is the density of self-closed links. Its value is set by the large scale (levels, global eights), not by the local dynamics of a core.

### Definitions

- **Exchange loop** [definition, from P3]. The joint path of the two points of one exchange: one step, one tick, with a side.
- **Point loop** [definition]. A closed path of one point. Without a lattice (section 9) the shortest point loop is two exchanges with the same partner, there and back: zero area, and its orientation is the side of the exchange (P3). The programs do not represent sides. On the square lattice (section 9) a point loop is a unit square of four exchanges, and its orientation is the sign of its area.
- **Eight** [definition]. Two successive point loops of one point with opposite orientation. The two loops of the electron in P6 are the point loops of O: O↔A, O↔A, then O↔B, O↔B. This is why each core exchange appears twice in the models; the proton core is A↔P₁, A↔P₁, P₁↔P₂, P₂↔C, P₂↔C, P₁↔P₂. P6 names the loops by their exchanges, and by P3 each exchange is an exchange loop; reading them as point loops of O, each of two exchanges with one partner, is a reading of P6. The orientation (the side) is represented in no program, so `graph_models.py` checks the paths of an eight, not its orientations.
- **Eight-point** [author, P10; reading stated by the author]. A point closed on itself by an eight contracted into it: the neutrino, and the point O of the electron. The two levels of P10 are the eight through the three points of the electron (the point loops O↔A and O↔B) and the eight contracted into O. Neither is on the circulation level of P8.
- **Self-closed link** (P18) [definition]. The closure of a point on itself by an eight contracted into it (an eight-point). An eight through several points — the electron's three-point eight, or the eight that a point of the medium follows among its neighbours (section 9) — is a path, not a self-closed link. A closure by one loop (a point, the photon, P14) is not a self-closed link either. [open] How a contracted eight is represented in models where every point exchanges with a neighbour at every tick.
- **Pattern, circulation, carrier** [definition]. The pattern of a particle is its sequence of exchanges. The circulation is the closed path of the pattern (level 2 of P8). The carrier is the point that plays O at a given tick.
- **Rest frame of places** [definition]. The frame in which the arrangement of places is at rest.
- **Mass** [definition, reading of P18]. m = M₀·ν, where ν is the number of closures per tick and M₀ = h/(ca) is one unit for all particles. A particle closes once per cycle of N ticks, so ν = 1/N and m = M₀/N = h/(cNa). This reads "density" in P18 as density in time; [candidate] that this is what P18 means. ν is counted in the particle's own ticks [author]: a free particle has the same mass at any speed, and the binding energy of an atom is not counted twice (section 10). [open] The mass of a composite system: the parts in their own ticks, the shared link (P17) and, by P5, the energy of motion of the parts; a rule is not yet given. The same counting (acts per common tick) is used for the source of gravity in the author's work on gravity, in preparation.

### Working assumptions

- **W1** [assumption]. All processes of a particle, decay included, run in step with its circulation (its own clock).
- **W2** [assumption]. The calculations of section 8 use the rest frame of places.
- **W3** [assumption]. The pattern of a particle moves one step per tick along its circulation, so the circulation runs at c. P1 fixes the speed of a change, not of a pattern; W3 is selected by data (section 3).
- **W4** [derived from P8, W3]. The frequency ω in E = ħω (P4), for a particle at rest, is the frequency of one circulation cycle: with R = ħ/(mc) (P8) and speed c (W3), one cycle has ω = c/R = mc²/ħ = E/ħ.
- **W5** [assumption]. A closed path holds a whole number of wavelengths: ∮p·dl = nh, with λ = h/p from P4 and no extra phase term.
- **W6** [assumption]. Places are embedded in Euclidean three-dimensional space, and rods measure the Euclidean distance. The speed of light is the speed of the shortest chain of moves (section 7).

Constants are CODATA 2018 values (Tiesinga et al. 2021).

## 2. The electron

**Form** [author, P6]. Three points A–O–B one step apart; the exchanges O↔A and O↔B are two point loops of opposite orientation through O. The core measures about 3a, where a is the step.

**Electron and medium** [derived from P2, P3, P10, in the medium with repeated pairings, whose choice is open (section 9)]. In the medium without a lattice (section 9) every point exchanges with one neighbour twice and then with another twice. The sequence of exchanges of O is the same as that of a point of this medium. The difference is the kind of closure: a point of the medium follows its eight as a path; O in addition carries an eight contracted into itself (P10). On an arrangement without a lattice a straight line A–O–B is defined only approximately, but the kind of closure is always defined. (On the square lattice the medium path is two unit squares instead; section 9.)

**The electron is a pattern** [author]. The electron is not three particular points. The points it passes through change all the time (section 3).

**Limits on the core** [estimate, data]:

| Data | What it limits | Limit |
|---|---|---|
| Contact interactions at LEP (LEP EW 2013): Λ ≈ 10 TeV (8–20 TeV by model, coupling g²/4π = 1) | core ≈ 3a < ħc/Λ ≈ 2·10⁻²⁰ m | a < 7·10⁻²¹ m |
| g/2 = 1.001 159 652 180 59(13) (Fan et al. 2023) | core correction ~ (R/ƛ_C)², ƛ_C = ħ/(m_e c) the reduced Compton wavelength (Brodsky, Drell 1980) | core < 3·10⁻¹⁹ m |
| Proton charge radius 0.84·10⁻¹⁵ m (Pohl et al. 2010); the proton is made of the same points | step | a < 10⁻¹⁵ m |

The second row uses the quadratic case of Brodsky and Drell, which holds if the substructure is chirally protected; the linear case gives a much stronger limit. The level of agreement between measurement and QED (~10⁻¹²) depends on which value of α is used; the caesium and rubidium recoil values disagree. This row is an estimate.

**The magnetic moment comes from the circulation, not from the core** [derived, algebra; computed, `verify_spin.py`]. A charge moving at no more than c on a circle of radius r gives μ_B only for r ≥ ħ/(m_e c). This radius is at least 2·10⁷ times the size of the core. The two point loops of the core have opposite orientation (P6), so their contributions cancel.

**Circulation and scattering** [estimate]. The circulation is the path of a pattern, not an extended body. At LEP energies (√s ≤ 209 GeV) the resolution time, about 3·10⁻²⁷ s, is 2.6·10⁶ times shorter than the period of the circulation (8.1·10⁻²¹ s): scattering sees a point.

**Field energy** [computed, `verify_spin.py`; data]. A classical Coulomb field outside the core (R = 2·10⁻²⁰ m) would hold e²/(8πε₀R) = 36 GeV, 7·10⁴ times the electron mass (the classical self-energy problem: Weisskopf 1939). The field of a charge below ħ/(m_e c) therefore cannot be classical in W either [open: what it is]. Interaction energy must enter the mass: platinum and titanium, with different fractions of electrostatic energy, fall alike to ~10⁻¹⁵ (MICROSCOPE 2022).

## 3. The particle circulation

**Speed** [assumption W3; data]. The pattern moves one step per tick: the circulation runs at c. A slower circulation would give μ < μ_B (below); [data] μ_e = 1.00116 μ_B.

**Radius** [author, P8; W4]. For a circulation at c on radius R, ω = c/R. With E = ħω (W4: ω of one cycle) and E = mc²: R = ħ/(mc), cycle length 2πR = h/(mc), momentum p = mc, frequency mc²/ħ — de Broglie's internal clock. Schrödinger's radius ħ/(2mc) (frequency of two cycles) contradicts P8, and it would give μ = μ_B/2, which the data exclude.

**Magnetic moment** [derived, algebra]. A charge e going round a loop of radius R at speed v is a current with moment μ = e·v·R/2. At v = c and R = ħ/(mc): μ = eħ/(2m) = μ_B. With R fixed by P8, μ/μ_B = v/c.

**Clock in motion** [candidate]. For a moving particle P5 gives E = γmc², and E = ħω would give ω = γω₀. The circulation, however, slows to ω₀/γ (section 8). This is de Broglie's clock paradox. He resolved it by separating the internal clock from the phase wave, and that resolution is adopted here: for a moving particle, ω in P4 is the frequency of the phase wave; the circulation is the internal clock; the two coincide at rest.

**The internal exchange is the step of the circulation** [author]. If a point moved along the core loops and along the circulation at the same time, the speed c would be shared and the circulation would be slower than c. There are no points at rest and no gaps (P2, P6): every movement is an exchange. The internal exchange of the core is itself the step of the circulation: the carrier moves one step forward, its partner one step back. The pattern moves on to new points; the points stay in the medium. [open] The rule by which the pattern passes to new points is not formulated. The models of section 9 do not implement it: they let the centre of a core move to any neighbour.

**Length and π** [derived from P1, P3]. An exchange loop is one step and one tick. There is no metric inside a step, so such a loop has no "π diameters". π appears only on long paths, as the ratio of steps along a path to steps across it. For the circulation of a particle (more than 3.5·10⁸ steps) the relation "length = 2πR" holds approximately and requires isotropy (section 7).

**Mass** [definition]. m = h/(c·N·a) = M₀/N, with N ≡ 2πR/a = λ_C/a. Until N is determined this is a change of variables. Mass ratios are inverse ratios of N; a and ħ drop out. [computed] Local constraints admit closed paths of every even length tested (4–24) at costs of order one (section 9). The models have no transition rule, so they cannot say which N a dynamics would select. [author, P18] N is set by the large scale. Neither N nor the step a is derived; deriving one N would give a and the tick a/c at once.

## 4. Spin and g = 2

**Angular momentum of the circulation** [derived, algebra]. A pattern with momentum mc on a circle of radius ħ/(mc) has angular momentum L_c = m·c·R = ħ, for any mass. The action per cycle is ∮p·dl = mc·2πR = h.

**What 2:1 gives** [author, P7; candidate]. The orientation returns after two cycles. If states carry a sign that changes under a rotation by 2π and returns under 4π, the representation is half-integer: j = ½, 3/2, … Two premises are not in P1–P18: that states have such a sign (there are no amplitudes in P1–P18; the phase e^(iπ) is a candidate, section 12), and that the angular momentum is ħ·j with the smallest j. With both, the orientation moment is ħ/2.

**ħ against ħ/2** [contradiction, open]. The circulation carries L_c = ħ. The measured spin is ħ/2. With S = ħ/2 and μ = μ_B, g = 2, as measured. With S = L_c = ħ, g = 1, which is excluded. No rule in this paper makes L_c unobservable or halves it. m·c·R cannot be set aside as "the action of a cycle": it is by definition an angular momentum, and the action per cycle is h. Nor are spin ½ and g = 2 one condition. They are two inputs: μ comes from the circulation (P8, W3, W4), S from P7 with the premises above.

**Carrier of 2:1** [candidate]. The eight of the core — the two point loops of O, of opposite orientation — returns the orientation after both loops. The circulation cycle cannot itself be one loop of an eight. Successive circulation cycles of opposite sense would reverse the circulation axis every cycle. Then either the magnetic moment averages to zero or, if the loops reverse with the axis, the sign of charge alternates. Both contradict the data. [open] How the period of the core's eight is tied to the circulation cycle so that the orientation returns after two cycles. P10's two levels are the eight through the three points and the eight contracted into O (section 1); neither is on the circulation level, so P10 and this paragraph agree.

**g = 2** [derived, conditional]. If S = ħ/2 and μ = μ_B, then g = 2mμ/(eS) = 2. A single circulation of a charge with S = L_c gives g = 1 for any radius. [data] Gyromagnetic measurements on ferromagnets (the Einstein–de Haas effect; Einstein and de Haas 1915 reported g ≈ 1, later measurements gave g ≈ 2) show that magnetization comes from electron spin with g ≈ 2.

**Exchange identity** [derived, algebra; illustrated in `models/spin_partners.py`]. In an exchange x↔y the points move by d and −d. The moment of the pair is r_x×d + r_y×(−d) = (r_x − r_y)×d = (−d)×d = 0, and its displacement is d − d = 0. For any sequence of exchanges the sum over all points of moment and displacement is therefore zero. This is a kinematic fact that holds whatever the points do; it has no consequence for patterns. It excludes the candidate that the counter-displacement of exchange partners carries −ħ/2: partners and carriers together carry exactly zero. The program shows the identity on model solutions and shows that the path of the pattern encloses a nonzero area while the sum over points is zero. The moment of a pattern is a property of its sequence, not a sum over points.

**Conservation for patterns** [candidate]. P2 says a state is a place in a sequence. Read together with P11, what is conserved are properties of whole sequences (complete cycles), not of instants. [open] State this rule precisely and derive collision kinematics from it.

**Proton** [open]. g = 2 refers to a particle with one eight. The proton has three loops and two eights with a common exchange P₁↔P₂ (P6, P9). Its value g = 5.586 is not derived from its form.

**Anomaly** [open]. The measured (g − 2)/2 = 0.00116 ≈ α/2π must come from the interaction of the circulation with the medium; this calculation is not done.

**Spin and statistics** [candidate]. In three dimensions an exchange of two identical objects gives only the signs ±1 (Leinaas, Myrheim 1977); for extended objects tied to a medium an exchange is topologically a rotation by 2π (Finkelstein, Rubinstein 1968). If self-closure is an exchange (P3), the sign −1 under exchange and the sign −1 under rotation by 2π are one property. This rests on the same phase premise as spin ½.

**No massless spin ½** [derived, conditional]. A particle whose spin ½ is carried by an eight-point (the electron's O, the neutrino) has a self-closed link (definition) and therefore mass (P18): the lightest neutrino has nonzero mass. For the proton, whose closure element is an exchange (P9), the argument does not apply; it is massive anyway. "No massless spin ½" in general needs every carrier of spin ½ to be a self-closed link [open]. The photon closes by one loop, has no self-closed link, is massless and has integer spin, which agrees.

## 5. Photon, neutrino and lepton number

**Photon** [author, P14; derived]. Each step of the photon is an exchange with a point of the medium. The side of the exchange is transverse to the step. Polarization is this side. If the side turns from step to step around the direction of motion, the polarization is circular with helicity ±1. Consequences:
- integer spin: one loop, no eight;
- no helicity 0: the side is transverse. [data] Only ±1 is observed (spin 1: Beth 1936; Wigner 1939 gives the same for any massless spin-1 particle; W gives a mechanism);
- one step per tick at every frequency, by P14; whether discreteness adds corrections at the scale of a step is not derived (section 11);
- no mass: no self-closed link.

[open] Why the moment is ħ; why the side turns uniformly; how a source sets the frequency.

**Neutrino** [author, P10; candidate]. A neutrino is a point closed by an eight. After one loop its orientation is reversed, after two restored: spin ½. It has mass (section 4) and a circulation of radius ħ/(mc) ≈ 4·10⁻⁶ m for m ≈ 0.05 eV. The two loops come in a time order. [computed, `models/charge_invariants.py`, section 3] If the loops lay in the plane of the circulation, the order of the loops would be a sign invariant under rotations and reflection and odd under loop reversal — a charge — and the neutrino would be charged; so they do not lie there [derived, within the state model of section 6]. If they lie across the circulation tangent v, the order is n₁·v = ±1, which P and C each reverse. This sign is not the measured helicity, which is the spin along the translational momentum: with the spin along c, spin·v = 0. Taking ν versus ν̄ to be this sign would make ν̄ the mirror image of ν (a Majorana-like pair), would map ν to ν̄ under C as under P, and would leave CP acting trivially. The observed pattern is different: ν left-handed and ν̄ right-handed (Goldhaber, Grodzins, Sunyar 1958), C and P each maximally violated, CP nearly conserved. [open] Where the loops of a neutrino lie, and what distinguishes ν from ν̄.

**Series of closures** [candidate]. A point: one loop, integer spin, massless (photon). A neutrino: an eight in a point, spin ½ (with the phase premise of section 4), nonzero mass. An electron: an eight through three points around an eight-point. A point of the medium follows an eight among its neighbours as its path; no eight is contracted into it (section 1).

**Lepton number** [candidate]. A lepton is a form whose closure element is an eight-point; the proton's closure element is the exchange P₁↔P₂. Lepton number is then the number of eight-points, counted with a sign, and baryon number the number of closure exchanges. Checks:
- electron capture p + e⁻ → n + ν: the eight-point of the electron leaves as the neutrino;
- pair creation γ → e⁻ + e⁺ (in the field of a nucleus): eight-points appear in pairs of opposite order;
- muon decay μ⁻ → e⁻ + ν̄_e + ν_μ: the number is conserved.

If lepton number is conserved, there is no neutrinoless double beta decay [candidate consequence]. [data] Total lepton number is conserved but flavour is not: neutrinos change flavour in flight (Super-Kamiokande 1998, SNO 2001). [open] What the sign of lepton number is. For the electron, if it is n·c, it equals the charge and the observed correlation (no e⁻ with lepton number −1) is automatic; for the neutrino it cannot be n·c, since the neutrino has no charge-like sign. One definition for both is not found. If ν and ν̄ were helicity states of one object, the candidate "no neutrinoless double beta decay" would lose its basis. [open] What flavour is in W.

## 6. Charge and antiparticle

**The reference of the sign** [definition; reading of P13]. P13 is read as follows: the sign of charge is n·c, set by the particle's own loops and circulation (below), and the global eight fixes which of the two signs counts as matter. A line in three-dimensional space goes into itself under any rotation about itself, so the line itself selects no plane and is not the reference; why the loops lie in the plane of the circulation is open problem 3. Read literally ("counted from the direction of the global eight"), P13 would let the spin-flip rule below reverse the charge, since n reverses in space; the reading above avoids this.

**The sign of charge** [computed, `models/charge_invariants.py`; derived]. A state of a particle is v (direction of motion along the circulation, which fixes the front end of the chain), c (axis of the circulation) and n (normal of the front loop, n ⟂ v). The operations are rotations, reflection P and loop reversal C (n → −n). A charge is a two-valued sign that is invariant under rotations and P and odd under C. Results:
- on all states there is no such sign. The states with n ⟂ c fall into two classes that P and C exchange together;
- restricted to loops in the plane of the circulation (n ∥ c), there is exactly one such sign: n·c = ±1. Analytic argument for all rotations: write n = s·c with s = ±1. SO(3) acts transitively on pairs (v, c) of orthogonal unit vectors, so an invariant depends only on s. P leaves the axial vectors n and c unchanged, so s is P-even, and C sends s to −s;
- with n ⟂ c the same form would be neutral;
- "n ∥ c versus n ⟂ c" is another two-valued invariant, but it is even under C, so it is not a charge;
- a photon has no circulation, only helicity: neutral. A neutrino: section 5;
- no rule fixing the side by the direction of motion alone, n = f(v), gives a rotation-invariant n·c: a rotation by 180° about v sends c to −c and leaves f(v).

The charge result therefore holds under the condition that the loops of a charged particle lie in the plane of its circulation. That condition is open problem 3, so uniqueness is [candidate] until it is solved.

**Spin flip** [candidate]. A spin flip in a magnetic field, as observed routinely (Fan et al. 2023), reverses the magnetic moment and therefore the circulation axis c. If n stayed fixed, n·c, the charge, would change sign [derived, algebra; `models/charge_invariants.py`, section 2]. The rule needed is that the loops are carried by the circulation: n reverses with c, and n·c is unchanged. The spin direction is along c. The current of a circulating negative charge then gives a magnetic moment opposite to the spin for the electron and along it for the positron, as observed [derived, conditional on S ∥ c].

**World and antiworld** [author; derived]. Reversing the global direction reverses the sides of all loops relative to their circulations, so all signs reverse at once. This is the transition "world ↔ antiworld" (P13), which is not annihilation. Locally a particle can carry either sign (pair creation).

**Requirements from data** [data]. Charge is conserved: neighbouring points cannot change sign separately. No gamma rays from matter–antimatter boundaries are seen, so such regions are not smaller than the visible universe (Cohen, De Rújula, Glashow 1998). Particle ↔ antiparticle oscillations occur only for neutral particles (K⁰, B⁰, D⁰).

[open] (1) Why the electron and the proton carry opposite signs in our world. (2) Under reversal of motion the sign of the electron changes and the proton's does not (two loops versus three) [computed]; W has no reversal (P15), and which W operation corresponds to the T symmetry tested in experiments is to be shown. (3) Why the loops of a charged particle lie in the plane of its circulation. (4) CP violation (Christenson et al. 1964) as a trace of the present global direction [candidate].

## 7. Isotropy

**The problem** [data]. On a regular lattice the path along a diagonal in steps differs from the path along an axis (Weyl 1949); the speed along the body diagonal of a cube is 0.577 of the speed along an axis. Resonator experiments limit any dependence of the speed of light on direction to ~10⁻¹⁸ (Nagel et al. 2015).

**Premise** [assumption W6]. The measure below uses Euclidean distance per hop and Delaunay neighbours, which are defined by Euclidean geometry. P1 defines length as a number of steps. The calculation therefore presupposes an embedding of places in Euclidean space whose metric is what rods measure.

**Average over realizations** [computed, `models/isotropy_check.py`]. One move to a neighbour per tick; speed along the axes compared with speed along the diagonals, in 20° cones, at 8–13 mean spacings (lattice: 24–33):

| Arrangement of places | Diagonal / axis |
|---|---|
| regular cubic lattice, equal moves | 0.577 |
| lattice of places, random times of moves along links (30 realizations) | 0.982 ± 0.002 |
| no lattice: random points equal in all directions, Delaunay neighbours (16 clouds) | 0.995 ± 0.003 |

A control with randomly rotated cones gives 1 within errors (0.9994 for the lattice of places, 1.0007 without a lattice). The second row is first-passage percolation with independent random times; it is not W dynamics, in which every point exchanges at every tick. Its conclusion — that random times on a lattice of places do not remove the anisotropy — holds for that model only.

**One realization** [computed, `models/isotropy_scaling.py`]. A physical arrangement is one realization. In 24 single clouds of 300 000 points, two quantities fall with distance r: the scatter over end points of the speed (distance per hop), and the spread over clouds of the per-cloud ratio of mean speeds along the diagonals and along the axes:

| r (mean spacings) | 4–5 | 8–10 | 15–19 | 26–32 |
|---|---|---|---|---|
| scatter over end points | 0.134 | 0.082 | 0.052 | 0.034 |
| spread of the per-cloud ratio | 0.033 | 0.017 | 0.010 | 0.007 |

Fits give r^(−β) with β = 0.73 and β = 0.80. The range is short, and the scatter includes the whole-number hop count. For first-passage distances in three dimensions the expected asymptotic value is β = 1 − χ ≈ 0.76 (fluctuation exponent χ ≈ 0.24), consistent with both fits. [open] The mean of the per-cloud ratio is +1.0 % ± 0.35 % at 8–10 spacings and falls with r (+0.2 % ± 0.15 % at 26–32). The construction is invariant under rotations, so the expectation is zero; `isotropy_check.py`, with a ratio of sums at 8–13 spacings, gives −0.5 % ± 0.3 %. The difference is not explained.

**Bound on the step** [open]. Nagel et al. bound a modulation of the speed of light correlated with the orientation of a resonator; the beam averages over many paths and, in time, over the places the laboratory moves through. A bound on the step needs the direction dependence of such a beam-averaged passage time in one realization. The scatter s(r) is not that: it is a scatter over end points and includes the whole-number hop count. No bound on the step is drawn from isotropy here.

**Result** [author, P16]. Places form no lattice. Systematic anisotropy is then zero on average; in one realization fluctuations of the propagation speed remain and fall with distance as above. P16 cannot be derived from P12 ("a lattice would distinguish axes without basis"): that argument assumes isotropy. P16 is a position.

## 8. Motion and clocks

**Across the motion** [derived, algebra; W1, W2]. A particle moves at speed v while its points still make one step per tick (P1). For a circulation in the plane perpendicular to v, the speed left for the circulation is √(c² − v²), and the period grows by γ = 1/√(1 − v²/c²).

**In the plane of motion** [computed, `verify_spin.py`]. For a point moving at c along a closed curve in a plane containing v, the period is (γ/c)·L, where L is the length of the curve stretched by γ along v. A circle contracted by γ along v gives exactly γ; a circle not contracted gives more (1.246 instead of 1.155 at v/c = 0.5). Two premises are needed:
- the period does not depend on the orientation of the orbit;
- the deformation is affine and leaves the transverse size unchanged.

Under these premises the contraction by γ along v is unique, since the perimeter of an ellipse grows strictly with each semi-axis. Without them, many non-elliptic closed curves, and a one-parameter family of ellipses, give the same period.

**Lorentz transformations** [candidate]. The Alexandrov–Zeeman theorem (Alexandrov 1950; Zeeman 1964): the one-to-one maps of Minkowski space (space of dimension at least two) that preserve the causal order of light cones are Lorentz transformations, translations and dilations. Statement: if the large-scale dynamics of W is invariant under the automorphisms of its causal order, and if that order is at large scales the light-cone order of Minkowski space, then a moving structure is the Lorentz transform of the same structure at rest (contraction by γ, dilation by γ, local time −γvx/c², reciprocity), with dilations excluded because the step is the same everywhere. Neither condition is shown:
- W has a discrete random set of places. Its order automorphisms are generically trivial, and "at large scales the causal order is a light cone" is not shown (section 7 computes isotropy of hop speeds only).
- W has a rest frame of places on which its dynamics is defined (below). The dynamics is therefore not set by the causal order alone, and the invariance premise is what would have to be proved.
- A core of three points one step apart (P6) cannot be contracted by γ. At the scale of the core the corrections are of order one.

What is derived is the explicit calculation above: period γ for both orientations of the orbit, with contraction by γ.

**Preferred frame**. The rest frame of places exists by definition (section 1): W is a theory of the Lorentz–Poincaré type (Lorentz 1904; FitzGerald 1889; Bell 1976; Brown 2005). [candidate] It shows only through discreteness. The size of the corrections and their power of a/L are not derived. If light has such corrections, LHAASO bounds the step (section 11).

## 9. Exchange models

Programs are in `models/` (Python, python-sat, networkx, scipy). The output of each program is in the `*_output.txt` file of the same name. In every model each point takes part in exactly one exchange with a neighbour at every tick. This follows from the positions: no point is ever at rest (P2), and one act takes one tick (P1). The models are constraint searches: they show that a configuration exists, they do not contain a transition rule, and they do not minimize rearrangements.

**Two media** [definition]. On the square lattice, `medium.py` and `eights.py` exclude a repeated pairing on consecutive ticks; there a point loop is a unit square. Without a lattice, `graph_models.py` builds its medium on exactly such a repetition (M₁ M₁ M₂ M₂); there a point loop is "there and back". [open] Which of the two the positions require.

[computed]
- **Line.** There are only two partitions into pairs. If repetition is excluded, they alternate: half the points move right and half left, one step per tick — two counter-streams at c, the massless motion of Feynman's checkerboard. If repetition is allowed (M₁ M₁ M₂ M₂), there is no streaming.
- **Plane** (8×8; `medium.py`, `eights.py`). Periodic sequences of the four pairings H0, H1, V0, V1 (horizontal and vertical domino pairings, even and odd), up to 8 ticks, without consecutive repeats, starting with H0: 2460 sequences, 1130 of them with all points on closed paths. `eights.py` finds 8 sequences in which every point follows an eight of two unit squares of opposite sign, for example H0 H1 V0 H1 V0 V1 H0 V1, period 8.
- **Core in the medium** (10×10; `defect.py`). Static cores of the electron and of the proton exist. The medium rearranges within radius 2 (electron) and 3 (proton): 96 exchanges in 24 ticks differ from the medium.
- **Moving core** (`moving.py`, `moving_cost.py`). A core whose centre shifts one step every tick along a closed path of length 4, 6, 8, 12, 16 or 24 exists for both forms. Exchanges differing from the medium per tick: electron 1.00–2.00, proton 1.67–3.00.
- **Volume** (6×6×6; `defect3d.py`). Cores e (line), p (chain of four), μ (cross) and τ (cross with a perpendicular eight) exist. Rearrangements per tick: e 3.50, μ 2.50, τ 3.67, p 3.50. The muon core is cheaper than the electron core: local cost does not order the masses.
- **Without a lattice** (`graph_models.py`). 2000 random points in a ball, Delaunay neighbours (14 on average). The graph is not bipartite, so the parity restriction of lattices disappears. (1) A medium exists: two disjoint partitions into pairs M₁, M₂ in the order M₁ M₁ M₂ M₂, and every point is back after four ticks. (2) Static cores of the electron and the proton next to each other exist; 5542 of 12 000 exchanges (46 %) in 12 ticks differ from the medium, which is fixed beyond radius 6. (3) The electron core goes round the proton core to the other side and back in 24 ticks without their points ever coinciding.

**What the models do not do** [by inspection of the programs]. The centre of a moving core moves to any neighbour each tick, independently of which core exchange happens, so "the internal exchange is the step of the circulation" (section 3) is not implemented. The point that plays O is not tracked: because the centre changes every tick, "O↔A, O↔A" acts on different pairs. In `graph_models.py` (3) the electron at each tick is a single forced exchange at its centre, which any point of the medium satisfies anyway; the constraint is weak.

[derived] (1) A dynamics in which all points always exchange is consistent, on a lattice and without one; the forms of e, p, μ and τ and a moving core are consistent with it. (2) If repetition is excluded, local patterns need at least two dimensions, although e and p are one-dimensional. (3) Local constraints admit closed paths of every even length tested (4–24) at costs of order one and give no mass hierarchy (ratios of order one, not 1836, with μ cheaper than e), in agreement with P18.

## 10. The hydrogen atom

**Positions** [author, P17]. Seven points and five loops. The electron cannot pass through the proton; it changes sides by exchange of places. The atom is a system of its own. The element common to the electron and the proton is the link between them.

**Why the dipole and the orbital moment must cancel** [data; computed]. [data] Ground-state hydrogen has no permanent dipole: there is no linear Stark effect (Bethe, Salpeter 1957). Its magnetic moment is that of the electron spin alone: no orbital moment (Phipps, Taylor 1927). [computed, `models/hydrogen_loops.py`, section 4] A single Kepler ellipse with the proton at a focus has a time-averaged dipole (3/2)·a·e along its apse line and an orbital moment. So successive loops must cancel both. Alternation of sides (below) is one way [candidate]; a precession of the apse line would cancel the dipole too.

**Levels** [computed, `models/hydrogen_loops.py`, section 1; premises: W5 and the 1/r attraction as data]. For every eccentricity, an ellipse with ∮p·dl = nh has semi-major axis n²a₀ and energy −13.606/n² eV; with the reduced mass, −13.598/n² eV. W5 fixes only the total action of a loop, not its angular momentum. W5 is applied to each ellipse; the action of the turns between loops (below) is not counted [open]. These are the nonrelativistic levels only; the 1S–2S frequency is measured to 10⁻¹⁵ (Parthey et al. 2011). Radial petals that leave the proton and return to it give the same levels but put far too much density at the proton; that model is excluded (`models/hydrogen_radial.py`, section 4).

**The ensemble of loops for n = 1** [assumption, chosen to fit]. Loops with semi-major axis a₀ and L² uniform on [0, 1] — the microcanonical ensemble. It is chosen because it reproduces the 1s momentum distribution exactly (Abrines, Percival 1966) [computed]. It also gives ⟨p⁴⟩ = 5 as for 1s, which is only the kinetic part of the 1S fine structure. Exact fractions: 1.45·10⁻² of loops come closer to the proton than the circulation radius of the electron (α·a₀), and 3.2·10⁻⁵ closer than the proton radius.

**Density at the proton** [computed, `models/hydrogen_loops.py`, section 3; contradiction with data]. The ensemble has density ρ(r) = √(2/r − 1)/(2π²) (atomic units), which diverges as r^(−1/2) at the proton. The density inside a radius x comes only from loops whose perihelion lies inside x. What P17 forbids depends on what "the proton" is in W: by P6 it is four points (≲ 3·10⁻²⁰ m with the LEP bound); its measured charge radius is 0.84 fm. If the proton is its four-point core, P17 removes almost nothing, and the mean density inside the charge radius is 68 times the 1s value |ψ(0)|² = 1/π (3.2 times inside α·a₀). If the electron may not enter the charge region, the density there is zero. The hyperfine splitting (21 cm line) and the Darwin term of the fine structure follow the 1s contact density. **Either way the ensemble is excluded as the spatial structure of the ground state.** It reproduces momenta and levels but not the electron density at the nucleus. [open] Whether some distribution of loops gives both ⟨1/L⟩ = 2 (needed for ⟨p⁴⟩ = 5) and the 1s contact density has not been examined. Separately, a charge circulating on R = ħ/(mc) smears the Coulomb potential by (R²/6)∇²V. That is 4/3 of the Dirac–Darwin term; if it took the place of the Darwin term, the 1S level would shift by an extra 58 GHz [computed, `verify_spin.py`; open, together with the field below ħ/(mc), section 2]. Its spatial averages ⟨r⟩ = 1.25 a₀ (1s: 1.5) and ⟨r²⟩ = 1.75 a₀² (1s: 3) must likewise be confronted with existing data (diamagnetic susceptibility, elastic scattering form factors) before anything is called a prediction.

**The other side** [candidate; computed, `models/hydrogen_loops.py`, section 4]. Rule: each next loop is the previous one turned by π about the axis through the proton that is perpendicular to the apse line and lies in the plane of the loop. Over a pair of loops the dipole and the orbital moment are then both zero (by construction). At the perihelion the electron is at distance q = a(1 − e) from that axis, so the turn never takes it through the proton. For a circular loop (e = 0) the axis is undefined. A reflection through the proton (r → −r) would cancel the dipole but keep L, so it is not enough. A reflection of the relative position through the centre of mass cannot be the change of sides. It would move the electron by about 2r in one act, faster than c (P1), and it is not an exchange in the sense of P3, since the proton does not move to the electron's place. [open] The duration of the turn and its speed (≤ c); the angular momentum 2L that changes at each turn — whether, under the sequence reading of P11 (section 4), the conserved moment is that of the pair of loops, which is zero; how the turn relates to "exchange of places with the directions of links taken into account" in P17.

**No radiation in the ground state** [candidate]. Below the ground sequence there is no stable sequence to pass to; excited states radiate because a lower sequence exists. The argument "a closed system does not change, so it does not radiate" fails: it decides closure by the absence of change, and it would also forbid the radiation of excited states (2P lives 1.6 ns). [open] The mechanism of emission of a photon in a transition.

**Lighter than its parts** [author, P17; candidate]. With the mass definition of section 1, a bound system that shares closures counts them once: ν(A ∪ B) = ν(A) + ν(B) − ν(A ∩ B). P17 counts 3 + 4 points and 2 + 3 loops with none shared. The shared element is the link between the electron and the proton, a closure common to both. Closures are counted in each particle's own ticks (section 1), so the electron's count does not depend on its motion in the atom, and the whole binding energy — kinetic and potential together — is carried by the link: M₀c²·ν(link) = 13.6 eV. Counting in the ticks of the rest frame of places is not used. [derived, algebra; W1] There the time dilation of the moving electron alone (⟨v²⟩ = α²c², virial theorem for 1/r) would lower its count by m⟨v²⟩/2 = 13.6 eV, so the link would count the binding a second time (27.2 eV against 13.6 eV measured). That count would also give a moving body the mass m/γ and a hot body less mass than a cold one, against P5. [data] The binding energy depends only on the charges and the reduced mass μ, as μα²c²/2: hydrogen 13.598 eV, deuterium 13.602 eV, positronium 6.80 eV. So the shared link is set by the pair of charges, not by the structure of the nucleus. [open] Derive ν(link) from the form of the link; a rule for the mass of a composite system (section 1).

**Model in three dimensions** [computed, `models/atom3d.py`, `models/graph_models.py`]. With every point in exactly one exchange per tick, a path of the electron core round the proton core to the other side and back exists, without their points ever coinciding, on a lattice and without one. In `atom3d.py` variant A, 15 574 of 20 736 exchanges (75 %) in the 12³ window differ from the medium, which is fixed only on the outer shell. Variant C ("through the proton is impossible") is excluded by a clause of the model and is a consistency check, not a result. The models implement neither loops nor the turn rule.

## 11. Limits from data, and determinism

- **Step** [estimate]. a < 7·10⁻²¹ m (LEP); independently a < 10⁻¹⁵ m (proton). If discreteness gives light quadratic corrections: a ≲ 3·10⁻²⁸ … 2·10⁻²⁷ m (below). Isotropy gives no bound yet (section 7). The tick a/c < 2.3·10⁻²⁹ s from LEP alone.
- **Steps in the electron circulation.** N_e = 2πR/a > 3.5·10⁸ (LEP bound).
- **Discreteness at black holes.** Corrections (a/r)² at r = GM/c² are ≈ 10⁻⁶⁰ or less for Sgr A* (with the LEP bound; Danielyan 2026).
- **Dispersion of light** [open]. P14 gives one step per tick at every frequency. Whether discreteness adds corrections of order (a/λ)ⁿ is not derived. The rest frame of places (section 8) suggests it might. If the corrections are quadratic, LHAASO (GRB 221009A) gives E_QG,2 > 6·10⁻⁸ E_Pl ≈ 7·10¹¹ GeV (its convention includes a factor 3/2). With corrections (a/ƛ)², E_a = ħc/a and a ≲ 3·10⁻²⁸ m; with (a/λ)², a ≲ 1.7·10⁻²⁷ m [estimate, `verify_spin.py`].
- **No probabilities** [author, P15] **and Bell** [candidate]. The violation of Bell inequalities (Bell 1964; loophole-free experiments: Hensen et al., Giustina et al., Shalm et al. 2015) excludes determinism with local rules and independent measurement settings. W keeps local exchanges (P1, P3). The option consistent with that is that the global eight — an ongoing motion of every point (P8), shared by sources and detectors from the start — correlates them. This is a violation of measurement independence, known as superdeterminism (Hossenfelder, Palmer 2020). Cosmic Bell tests (Rauch et al. 2018) push any such common cause back by billions of years, which the global eight allows. The other option, a dynamical non-local influence, would break the premise of local exchanges. This is a choice, not a derivation.
- **Linearity** [derived, mathematics]. A nonlinear rule for single configurations generates an exactly linear evolution of weights over all configurations (Koopman 1931). This does not give interference or the rule |Σ|². [open] Why frequencies of outcomes follow |Σ|²; why decays are exponential.
- **Particles as orbits** [candidate]. A particle is a periodic orbit of the transition rule; stability follows from discreteness (an orbit on a discrete set of states repeats exactly).

## 12. Open problems

1. **ħ against ħ/2** [contradiction]: the circulation carries ħ, the spin is ħ/2 (section 4).
2. How the core's eight is tied to the circulation cycle so that the orientation returns after two cycles; the phase e^(iπ) as a place in a cycle [candidate]; the rule 2:1 from a single exchange [candidate].
3. Why the loops of a charged particle lie in the plane of the circulation, and where those of a neutrino lie; what distinguishes ν from ν̄; why the electron and the proton have opposite signs; one definition of the sign of lepton number for the electron and the neutrino; which W operation corresponds to T symmetry; CP violation.
4. The cycle lengths N of particles (mass ratios), and with them the step and the tick without calibration (P18; the theory of levels); whether P18 means density in time; the mass of a composite system (closures in own ticks, the shared link and the energy of motion of the parts).
5. The value g = 5.586 of the proton from its three loops; the anomaly (g − 2)/2 ≈ α/2π; α.
6. The moment ħ of the photon; why its side turns uniformly; how a source sets the frequency; corrections to the propagation of light at the scale of a step.
7. What neutrino flavour is, given that it changes in flight.
8. The field of a charge below ħ/(mc); the mechanism of attraction of charges and the origin of the 1/r law.
9. Hydrogen: a distribution of loops with the right contact density, or another model; the turn between loops, its speed and its angular momentum; the shared link ν(link); the mechanism of emission.
10. The rule by which the pattern of a core passes to new points; which medium (with or without repeated pairings) the positions require; how a contracted eight is represented in models where every point exchanges at every tick.
11. Conservation laws as properties of whole sequences.
12. Lorentz invariance at large scales: whether the causal order of a random arrangement of places becomes a light-cone order, and whether the dynamics respects it; the direction dependence of a beam-averaged passage time in one realization, and with it a bound on the step from isotropy.
13. Frequencies |Σ|² and exponential decay.

## 13. What is known and what W adds

| Result | Known from | Here |
|---|---|---|
| Internal clock with frequency mc²/ħ; phase wave versus clock | de Broglie 1924 | radius ħ/(mc) from P8; de Broglie's resolution adopted |
| Circulation or trembling at ħ/(2mc), frequency 2mc²/ħ | Schrödinger 1930 | W uses ħ/(mc) (P8); ħ/(2mc) would give μ_B/2, excluded by data |
| Spin ½, μ_B and g = 2 from a charge circulating at c | Barut, Zanghi 1984; Hestenes 1990; Rivas 1994, 2001 | μ from the circulation, S from the rule 2:1: two inputs; ħ against ħ/2 is open |
| g = 2 without relativity | Lévy-Leblond 1967 | — |
| Dirac equation from one step per tick on a lattice | Feynman, Hibbs 1965; Bialynicki-Birula 1994; D'Ariano, Perinotti 2014 | on a line: two counter-streams (massless case), if repeated pairings are excluded; no derivation of the Dirac equation |
| Massless spin 1 has only helicities ±1 | Wigner 1939 | a mechanism: the side of exchange is transverse to the step |
| Neutrinos left-handed | Goldhaber et al. 1958 | the loops of a neutrino cannot lie in its circulation plane; ν versus ν̄ is open |
| Spin and statistics | Fierz 1939; Pauli 1940; Finkelstein, Rubinstein 1968 | one property of exchange [candidate] |
| Causal order determines Lorentz transformations | Alexandrov 1950; Zeeman 1964 | applies only under a premise not shown in W [candidate] |
| Contraction from internal dynamics | Lorentz 1904; FitzGerald 1889; Bell 1976; Brown 2005 | contraction unique among affine deformations with unchanged transverse size, for a clock in its plane |
| Momentum distribution of 1s from classical orbits | Abrines, Percival 1966 | the same ensemble gives the wrong electron density at the proton (68× or 0, by the meaning of "proton"): excluded |

**What W adds:**
1. Charge as the only sign invariant under rotations and reflection and odd under loop reversal, on condition that the loops lie in the plane of the circulation; the rule that keeps the charge under a spin flip (section 6).
2. The photon as a point whose side of exchange turns around the direction of motion: spin 1 without helicity 0 and without mass (section 5).
3. The neutrino as an eight in a point whose loops cannot lie in its circulation plane (section 5).
4. Places without a lattice: isotropy on average, and a measured fall of the propagation fluctuations in one realization (section 7).
5. A clock moving in its own plane keeps the period γ only if contracted by γ (section 8).
6. Exchange models in which every point exchanges at every tick, on lattices and without one, with cores of e, p, μ, τ and an atom path (section 9).
7. Tests that the theory's own models fail and that are recorded here: radial petals and the microcanonical loop ensemble in hydrogen (section 10); the contradiction ħ against ħ/2 (section 4); mass counted in the ticks of the rest frame of places (section 10); the neutrino loop order as a helicity (section 5).
8. Numerical limits on the step of space (sections 2, 11).

## References

- L. de Broglie. Recherches sur la théorie des quanta [Research on the theory of quanta]. Thèse, Paris, 1924.
- E. Schrödinger. Über die kräftefreie Bewegung in der relativistischen Quantenmechanik [On the free motion in relativistic quantum mechanics]. Sitzungsber. Preuss. Akad. Wiss., Phys.-Math. Kl. 24, 418 (1930).
- G. F. FitzGerald. The ether and the Earth's atmosphere. Science 13, 390 (1889).
- H. A. Lorentz. Electromagnetic phenomena in a system moving with any velocity smaller than that of light. Proc. R. Neth. Acad. Arts Sci. 6, 809 (1904).
- A. Einstein, W. J. de Haas. Experimenteller Nachweis der Ampèreschen Molekularströme [Experimental proof of Ampère's molecular currents]. Verh. Dtsch. Phys. Ges. 17, 152 (1915).
- P. Phipps, J. B. Taylor. The magnetic moment of the hydrogen atom. Phys. Rev. 29, 309 (1927).
- B. O. Koopman. Hamiltonian systems and transformations in Hilbert space. Proc. Natl. Acad. Sci. 17, 315 (1931).
- R. A. Beth. Mechanical detection and measurement of the angular momentum of light. Phys. Rev. 50, 115 (1936).
- E. Wigner. On unitary representations of the inhomogeneous Lorentz group. Ann. Math. 40, 149 (1939).
- V. F. Weisskopf. On the self-energy and the electromagnetic field of the electron. Phys. Rev. 56, 72 (1939).
- M. Fierz. Helv. Phys. Acta 12, 3 (1939); W. Pauli. The connection between spin and statistics. Phys. Rev. 58, 716 (1940).
- H. Weyl. Philosophy of Mathematics and Natural Science. Princeton, 1949.
- A. D. Alexandrov. On Lorentz transformations. Uspekhi Mat. Nauk 5, 187 (1950).
- H. A. Bethe, E. E. Salpeter. Quantum Mechanics of One- and Two-Electron Atoms. Springer, 1957.
- M. Goldhaber, L. Grodzins, A. W. Sunyar. Helicity of neutrinos. Phys. Rev. 109, 1015 (1958).
- J. S. Bell. On the Einstein Podolsky Rosen paradox. Physics 1, 195 (1964).
- E. C. Zeeman. Causality implies the Lorentz group. J. Math. Phys. 5, 490 (1964).
- J. H. Christenson, J. W. Cronin, V. L. Fitch, R. Turlay. Evidence for the 2π decay of the K₂⁰ meson. Phys. Rev. Lett. 13, 138 (1964).
- R. P. Feynman, A. R. Hibbs. Quantum Mechanics and Path Integrals. McGraw-Hill, 1965.
- R. Abrines, I. C. Percival. Classical theory of charge transfer and ionization of hydrogen atoms by protons. Proc. Phys. Soc. 88, 861 (1966).
- J.-M. Lévy-Leblond. Nonrelativistic particles and wave equations. Commun. Math. Phys. 6, 286 (1967).
- D. Finkelstein, J. Rubinstein. Connection between spin, statistics, and kinks. J. Math. Phys. 9, 1762 (1968).
- J. S. Bell. How to teach special relativity. Progress in Scientific Culture 1 (1976).
- J. M. Leinaas, J. Myrheim. On the theory of identical particles. Nuovo Cimento B 37, 1 (1977).
- S. J. Brodsky, S. D. Drell. Anomalous magnetic moment and limits on fermion substructure. Phys. Rev. D 22, 2236 (1980).
- A. O. Barut, N. Zanghi. Classical model of the Dirac electron. Phys. Rev. Lett. 52, 2009 (1984).
- D. Hestenes. The zitterbewegung interpretation of quantum mechanics. Found. Phys. 20, 1213 (1990).
- M. Rivas. J. Math. Phys. 35, 3380 (1994); Kinematical Theory of Spinning Particles. Kluwer, 2001.
- I. Bialynicki-Birula. Weyl, Dirac, and Maxwell equations on a lattice as unitary cellular automata. Phys. Rev. D 49, 6920 (1994).
- A. G. Cohen, A. De Rújula, S. L. Glashow. A matter–antimatter universe? Astrophys. J. 495, 539 (1998).
- Super-Kamiokande Collaboration, Y. Fukuda et al. Evidence for oscillation of atmospheric neutrinos. Phys. Rev. Lett. 81, 1562 (1998).
- SNO Collaboration, Q. R. Ahmad et al. Measurement of the rate of νe + d → p + p + e⁻ interactions produced by ⁸B solar neutrinos at the Sudbury Neutrino Observatory. Phys. Rev. Lett. 87, 071301 (2001).
- H. R. Brown. Physical Relativity. Oxford University Press, 2005.
- R. Pohl et al. The size of the proton. Nature 466, 213 (2010).
- C. G. Parthey et al. Improved measurement of the hydrogen 1S–2S transition frequency. Phys. Rev. Lett. 107, 203001 (2011).
- ALEPH, DELPHI, L3, OPAL, LEP Electroweak Working Group. Electroweak measurements in electron–positron collisions at W-boson-pair energies at LEP. Phys. Rept. 532, 119 (2013).
- G. M. D'Ariano, P. Perinotti. Derivation of the Dirac equation from principles of information processing. Phys. Rev. A 90, 062106 (2014).
- M. Nagel et al. Direct terrestrial test of Lorentz symmetry in electrodynamics to 10⁻¹⁸. Nat. Commun. 6, 8174 (2015).
- B. Hensen et al., Nature 526, 682 (2015); M. Giustina et al., Phys. Rev. Lett. 115, 250401 (2015); L. K. Shalm et al., Phys. Rev. Lett. 115, 250402 (2015).
- D. Rauch et al. Cosmic Bell test using random measurement settings from high-redshift quasars. Phys. Rev. Lett. 121, 080403 (2018).
- S. Hossenfelder, T. Palmer. Rethinking superdeterminism. Front. Phys. 8, 139 (2020).
- E. Tiesinga, P. J. Mohr, D. B. Newell, B. N. Taylor. CODATA recommended values of the fundamental physical constants: 2018. Rev. Mod. Phys. 93, 025010 (2021).
- MICROSCOPE Collaboration, P. Touboul et al. MICROSCOPE mission: final results of the test of the equivalence principle. Phys. Rev. Lett. 129, 121102 (2022).
- X. Fan, T. G. Myers, B. A. D. Sukra, G. Gabrielse. Measurement of the electron magnetic moment. Phys. Rev. Lett. 130, 071801 (2023).
- LHAASO Collaboration, Z. Cao et al. Stringent tests of Lorentz invariance violation from LHAASO observations of GRB 221009A. Phys. Rev. Lett. 133, 071501 (2024).
- H. Danielyan. Black Holes in the W-Structure. Zenodo, 2026. https://doi.org/10.5281/zenodo.24016340

Prepared with the assistance of AI models (Claude, ChatGPT, DeepSeek).
