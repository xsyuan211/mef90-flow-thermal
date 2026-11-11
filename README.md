# This branch: petroleum-engineering Flow extension (2026-11)

Upstream `mef90/vDef` is the reference implementation of variational fracture.
This branch adds a configurable Darcy/Biot Flow side, time-dependent material
laws, conservation/energy diagnostics, and a CI-enforced hard-verdict system:

* `FLOW_FIX_CHECKLIST.md` — gate status, root causes and CI acceptance records.
* `RESEARCH_POSITIONING.md` — comparison with prior work, capability matrix,
  honest limitations and citation leads.
* `RATE_FRACABILITY_ROADMAP.md` — injection-rate / shale-fracability matching:
  literature survey, code entry-point audit, similarity mapping, the Gate 4.15
  route, Stage B0 acceptance (CI #125), Stage B1 acceptance (CI #127),
  Stage B2 acceptance (CI #128; 19 verdict groups green in one run) and the
  B1/B2 identifiability record (sec. 12-13). Two corrections: propagation-
  regime scalings do NOT apply to the pre-initiation storage-dominated window
  (Lecampion et al. 2017, doi:10.1002/2016JB013183), and bedding anisotropy is
  likewise NOT identifiable in the surface energy before initiation (CI #127
  measured only a 0.41% shift). Sec. 14 recomputes tables A-C as A'-C' from
  published China continental shale-oil parameters (S01-S20) and documents why
  a field rate-optimization table CANNOT be delivered from a pre-initiation
  storage-dominated window (CI #130 accepted this recomputation as the 20th
  verdict group). Sec. 15 documents Stage B3: the AT2 toughness ladder DID
  cross the stiffness-half threshold, but CI #132 classified the result as
  UNIFORM_SOFTENING, i.e. B3_INITIATION_ACCESS = NOT_ESTABLISHED, NOT a
  localised crack (sec. 15.9). No Fortran solver changes at any stage.
* `vDef/data/terzaghi415_scale.yaml` / `terzaghi415_inject.yaml` — Stage B0
  similarity-calibrated zero-source control and constant-rate injection decks
  (T = 1.57e-7, storage-dominated); 4-line / 1-line code-diff contracts.
* `vDef/data/terzaghi415_b1_*.yaml` / `terzaghi415_b2_*.yaml` — Stage B1
  Q x S x Gc x k x mu one-factor identifiability matrix, and Stage B2
  variable-rate / bedding-anisotropy decks (each one code line from B0).
* `vDef/data/check_gate415_rate_fracability.py` — Stage B0 V1-V8 hard verdict,
  including the `[1/s]` vs `[Pa/s]` source-unit falsification test.
* `vDef/data/check_gate415_b1_identifiability.py` (A1-A8, CI #127 PASS) /
  `check_gate415_b2_extension.py` (L1-L8, CI #128 PASS) — Stage B1/B2 verdicts: beta = f/S
  sees only (f,S), Gc is visible in the phase-field surface energy; f(t)=f0*t
  has a discrete t^2 anchor (CI #127: dp10 = 5.5e4 Pa = 0.55 x B0); 3D tetra4
  consistency. L4 is now a falsifiable NULL proposition for pre-initiation
  anisotropy: wired but not identifiable (two-sided band + ds < damage_dev).
* `vDef/data/shale_oil_params.json` / `check_shale_oil_params.py` — shale-oil
  literature parameter table and the M1-M8 recomputation verdict. Reproduces
  the first-edition tables A/B/C at printed precision as a calibration anchor,
  then recomputes A'-C' with three independent formula implementations
  (reference / literal / expanded power law) so any exponent drift is
  falsifiable. Keeps k_core (core gas permeability) strictly separate from
  k_eff (pressure-decline inverted permeability, ~37.5x smaller at rep) and
  derives S from phi*(c_pore+c_fluid). Pure post-processing: no Fortran, no
  deck, no solver run; `--md` regenerates the published tables, `--selftest`
  runs 1 positive + 33 mutations + 3 invariances. CI #130 `success`:
  SHALE_OIL_PARAMS_RESULT: PASS (M1-M8), 20 verdict groups green in one run.
  Acceptance is inference-level (run/job/step success + deferred exit-1 on
  gso_rc), NOT byte-level log replay: the logs/artifact endpoints returned 401
  in that run; see roadmap sec. 14.9.
* `vDef/data/terzaghi415_b3_*.yaml` / `check_gate415_b3_initiation.py` — Stage
  B3 fracture-INITIATION ACCESS TEST: four one-line-per-decade Gc decks
  (27000 -> 2.7 N/m, only `cs0001.fracturetoughness`) plus two reused B0/B1
  EXO anchors, and the N1-N8 hard verdict. AT2 analytic anchor
  alpha_eq = 2*Psi/(Gc/ell + 2*Psi), calibrated on the CI #125 real damage
  field (2*Psi = 1.998714e1 Pa, Gc* = 9.650647 N/m at a(alpha)=1/2) and
  independently checked against CI #127's gb1_surface_ratio_gc = 2.000410.
  CI #132 (first run) DID cross the stiffness-half threshold (alpha=1,
  measured bracket [2.7,27] N/m) but classified the last rung as
  UNIFORM_SOFTENING (dmax/dmean = 1.000004), NOT a localised crack. N5 was
  therefore rewritten from "must localise" into a falsifiable mechanism
  classification (LOCALISED > 1.5 / UNIFORM <= 1.05; the ambiguous middle
  band is a hard FAIL). The verdict now emits B3_REGIME = UNIFORM_SOFTENING
  and B3_INITIATION_ACCESS = NOT_ESTABLISHED: a PASS certifies the DIAGNOSIS,
  never the desired physical outcome. beta = f/S stays bitwise invariant;
  zero Fortran / zero deck / zero production changes. No propagation-regime
  scaling, fracture length or rate-optimization claim.
  `--selftest` 1 positive + 32 mutations + 5 invariances; a full replay of
  the #132 real artefacts gives N1-N8 PASS with all 20 gb3_* metrics
  bit-identical to the original job log. CI re-verification pending; see
  roadmap sec. 15.9.
* `vDef/data/check_gate*.py` / `check_shale_oil_params.py` — the CI hard-verdict
  suite: 18 gate verdicts plus the M1-M8 literature recomputation (20 verdict
  groups green in one run at CI #130; B3 is the 21st, pending re-verification),
  plus two auxiliary Exodus / displacement checks. Each verdict has a
  positive/negative offline selftest.
* `.github/workflows/build.yml` — build, diagnostic decks and deferred gate exit.

Every verdict runs in CI; auditing a change requires no local compilation or
solver run. Upstream must still be cited (Zenodo DOI 10.5281/zenodo.3242131 and
the four Francfort–Marigo / Bourdin papers listed below).

---

mef90 / vDef: A reference implementation of the variational approach to fracture,
as described in:

* Francfort, G. and Marigo, J.-J. (1998). Revisiting brittle fracture as an energy minimization problem. J. Mech. Phys. Solids., 46(8):1319–1342. [DOI:10.1016/S0022-5096(98)00034-9](http://doi.org/10.1016/S0022-5096(98)00034-9)
* Bourdin, B., Francfort, G., and Marigo, J.-J. (2000). Numerical experiments in revisited brittle fracture. J. Mech. Phys. Solids, 48(4):797–826. [DOI:10.1016/S0022-5096(99)00028-9](http://doi.org/10.1016/S0022-5096(99)00028-9)
* Bourdin, B., Francfort, G., and Marigo, J.-J. (2008). The variational approach to fracture. J. Elasticity, 91(1-3):1–148. [DOI:10.1007/s10659-007-9107-3](http://doi.org/10.1007/s10659-007-9107-3)
* Bourdin, B. (2007). Numerical implementation of a variational formulation of quasi-static brittle fracture. Interfaces Free Bound., 9:411–430. [DOI:10.4171/IFB/171](http://doi.org/10.4171/IFB/171)


This software is free software, distributed under the 2-clause BSD license. A copy of the license is included in the LICENSE file.
We cordially ask that any published work derived from this code, or utilizing it references the software as [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3242131.svg)](https://doi.org/10.5281/zenodo.3242131), and the above-mentioned published works:

```

@software{mef90,
	author = {Bourdin, B.},
	title = {{mef90/vDef}: variational models of defect mechanics},
	url = {https://github.com/bourdin/mef90},
	year = {2019},
	doi = {10.5281/zenodo.3242131}}

@article{Bourdin-Francfort-EtAl-2000a,
	Author = {Bourdin, B. and Francfort, G.A. and Marigo, J.-J.},
	Journal = {J. Mech. Phys. Solids},
	Number = {4},
	Pages = {797--826},
	Title = {Numerical experiments in revisited brittle fracture},
	Volume = {48},
	Year = {2000},
	DOI = {10.1016/S0022-5096(99)00028-9}}

@article{Bourdin-2007a,
	Author = {Bourdin, B.},
	Journal = {Interfaces Free Bound.},
	Month = {08},
	Pages = {411--430},
	Title = {Numerical implementation of a variational formulation of quasi-static brittle fracture},
	Volume = {9},
	Year = {2007},
	DOI = {10.4171/IFB/171}}

@article{Bourdin-Francfort-EtAl-2008b,
	Author = {Bourdin, B. and Francfort, G.A. and Marigo, J.-J.},
	Journal = {J. Elasticity},
	Number = {1-3},
	Pages = {1--148},
	Title = {The Variational Approach to Fracture},
	Volume = {91},
	Year = {2008},
	DOI = {10.1007/s10659-007-9107-3}}

@article{Francfort-Marigo-1998,
	Author = {Francfort, G.A. and Marigo, J.-J.},
	Journal = {J. Mech. Phys. Solids.},
	Number = 8,
	Pages = {1319--1342},
	Title = {Revisiting brittle fracture as an energy minimization problem},
	Volume = 46,
	Year = 1998,
	DOI = {10.1016/S0022-5096(98)00034-9}}
```
