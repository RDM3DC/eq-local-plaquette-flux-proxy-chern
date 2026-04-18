# Local Plaquette-Flux Proxy Chern

**ID:** `eq-local-plaquette-flux-proxy-chern`  
**Tier:** derived  
**Score:** 92  
**Units:** OK  
**Theory:** PASS-WITH-ASSUMPTIONS  

## Equation

$$
\rho_p=\frac{1}{2\pi}\arg\!\big(g_{p,0}g_{p,1}\,\overline{g_{p,2}}\,\overline{g_{p,3}}\big),\quad C_{\mathrm{loc}}=\sum_p \rho_p
$$

## Description

Solver-backed local plaquette-flux density and summed proxy Chern for damaged HAFC/EGATL lattices. This makes the previously abstract local Chern signal computable from an oriented plaquette product, giving a scalable post-damage topological proxy when dense Bott-style postprocessing becomes impractical.

## Assumptions

- Plaquette bond ordering and orientation are fixed consistently so the oriented plaquette product is comparable across runs.
- The complex bond responses g_{p,k} remain defined on the monitored plaquettes, and the plaquette-product phase is interpreted on a principal branch before summation.
- C_{\mathrm{loc}} is used as a scalable proxy diagnostic on damaged finite lattices, not claimed to be an exactly quantized Bott invariant in every finite-size regime.
- The same local phase convention and damage protocol are used when comparing healthy, central-strip, and top-edge runs.

## Repository Structure

```
images/       # Visualizations, plots, diagrams
derivations/  # Step-by-step derivations and proofs
simulations/  # Computational models and code
data/         # Numerical data, experimental results
notes/        # Research notes and references
```

## Links

- [TopEquations Registry](https://rdm3dc.github.io/TopEquations/registry.html)
- [TopEquations Main Repo](https://github.com/RDM3DC/TopEquations)
- [Certificates](https://rdm3dc.github.io/TopEquations/certificates.html)

---
*Part of the [TopEquations](https://github.com/RDM3DC/TopEquations) project.*
