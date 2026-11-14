# The Electron in the W-Structure: Spin, Magnetic Moment, Charge and Motion

**Author:** Hrachya Danielyan · **Version:** 1.1 · **License:** CC BY 4.0 · **DOI of version 1.0:** [10.5281/zenodo.24220735](https://doi.org/10.5281/zenodo.24220735)

In the W-Structure, space is points that never rest: every point exchanges places with a neighbour at every tick, and one exchange is one step and one tick. A particle is a stable pattern of exchanges that moves along its own closed path, the particle circulation. This work traces what the author's positions imply for the electron, and labels every statement as a position, definition, assumption, derivation, computation, candidate, measurement or open problem.

Main results, with their status:
- a circulation at c with radius ħ/(mc) gives the magnetic moment μ_B; with spin ħ/2 from the rule "orientation returns after two cycles", g = 2. The circulation itself carries ħ: the contradiction ħ against ħ/2 is open;
- if the loops of a charged particle lie in the plane of its circulation, charge is the only sign invariant under rotations and reflection and odd under loop reversal; the loops of a neutrino cannot lie in its circulation plane, or it would be charged;
- places without a lattice are isotropic on average; in one realization the fluctuations of the propagation speed fall as r^(−0.73…0.80), as expected for first passage in three dimensions;
- a clock moving in its own plane keeps the period γ only if contracted by γ; full Lorentz invariance is a candidate under a stated premise;
- hydrogen as Kepler loops with ∮p·dl = nh gives the levels and the 1s momentum distribution, but the electron density at the proton is wrong (68 times the 1s value, or zero, depending on what "the proton" is), so the ensemble is excluded by the 21 cm line;
- the mass of a particle is its closures per own tick; the mass of a composite system is its energy at rest (P5), in which the motion of its parts counts — as lunar laser ranging confirms.

Open problems are listed in section 12 of the paper.

## Contents

- `electron.md` — the paper.
- `verify_spin.py` — arithmetic checks for sections 1, 2, 3, 4, 8, 10 and 11; output in `verify_spin_output.txt`.
- `models/` — exchange models and checks. The output of each program is in the `*_output.txt` file of the same name:
  - `medium.py`, `eights.py` — a medium of figure-eights on a square lattice;
  - `defect.py`, `defect3d.py` — particle cores in the medium (2D, 3D);
  - `moving.py`, `moving_cost.py` — moving cores;
  - `spin_partners.py` — the exchange identity and why it says nothing about the moment of a pattern;
  - `charge_invariants.py` — two-valued signs under rotations, reflection and loop reversal; neutrino loop order;
  - `isotropy_check.py`, `isotropy_scaling.py` — isotropy on average and in one realization;
  - `graph_models.py` — the medium, cores and the atom path without a lattice;
  - `hydrogen_loops.py` — the hydrogen model of the paper; `hydrogen_radial.py` — an excluded model kept as a record; `atom3d.py` — the atom path on a lattice.

## Changes in 1.1

- Mass of a composite system defined as its energy at rest (P5): the motion of the parts counts, acts common to several parts are counted once.
- Hydrogen: the electron's motion adds +13.6 eV, the link carries −27.2 eV; version 1.0 assigned the whole −13.6 eV to the link.
- `verify_spin.py` section 6: counting closures (instead of acts) for the source of gravity is excluded by lunar laser ranging.

## Running the checks

Python 3 with `numpy`, `scipy`, `python-sat` and `networkx`:

```
pip install numpy scipy python-sat networkx
python3 verify_spin.py
cd models && python3 charge_invariants.py
```

Most programs finish in under a minute. `isotropy_check.py` takes about 1.5 minutes, `atom3d.py` under a minute, `isotropy_scaling.py` (24 clouds) about 10 minutes and under 1 GB of memory.

## Related work

- H. Danielyan. Black Holes in the W-Structure. Zenodo, 2026. https://doi.org/10.5281/zenodo.24016340

## How to cite

H. Danielyan. The Electron in the W-Structure: Spin, Magnetic Moment, Charge and Motion. Version 1.0. Zenodo, 2026. https://doi.org/10.5281/zenodo.24220735 (the DOI of version 1.1 will be added after its release)

See also `CITATION.cff`.

Prepared with the assistance of AI models (Claude, ChatGPT, DeepSeek).
