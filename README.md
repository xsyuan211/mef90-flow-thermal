# This branch: petroleum-engineering Flow extension (2026-11)

Upstream `mef90/vDef` is the reference implementation of variational fracture.
This branch adds a configurable Darcy/Biot Flow side, time-dependent material
laws, conservation/energy diagnostics, and a CI-enforced hard-verdict system:

* `FLOW_FIX_CHECKLIST.md` — gate status, root causes and CI acceptance records.
* `RESEARCH_POSITIONING.md` — comparison with prior work, capability matrix,
  honest limitations and citation leads.
* `RATE_FRACABILITY_ROADMAP.md` — injection-rate / shale-fracability matching:
  literature survey, code entry-point audit, similarity mapping, the Gate 4.15
  route, and the Stage B0 similitude-calibration + source-unit calibration
  record (no Fortran solver changes at any stage).
* `vDef/data/terzaghi415_scale.yaml` / `terzaghi415_inject.yaml` — Stage B0
  similarity-calibrated zero-source control and constant-rate injection decks
  (T = 1.57e-7, storage-dominated); 4-line / 1-line code-diff contracts.
* `vDef/data/check_gate415_rate_fracability.py` — Stage B0 V1-V8 hard verdict,
  including the `[1/s]` vs `[Pa/s]` source-unit falsification test.
* `vDef/data/check_gate*.py` — 15+ falsifiable verdict scripts, each with a
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
