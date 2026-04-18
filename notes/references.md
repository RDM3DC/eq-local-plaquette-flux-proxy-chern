# References & related work

## Primary

- T. Fukui, Y. Hatsugai, and H. Suzuki,
  *Chern numbers in discretized Brillouin zone: Efficient method of
  computing (spin) Hall conductances*,
  J. Phys. Soc. Jpn. **74**, 1674 (2005).
  arXiv: [cond-mat/0503172](https://arxiv.org/abs/cond-mat/0503172).
  Establishes the lattice-link-variable construction that the registered
  equation specialises to when $g_{p,k}$ are unit-modulus link
  variables.

## Adjacent constructions

- T. A. Loring and M. B. Hastings,
  *Disordered topological insulators via $C^*$-algebras*,
  EPL **92**, 67004 (2010). The Bott index that this proxy is meant to
  replace in the damaged regime.
- E. Prodan,
  *Disordered topological insulators: a non-commutative geometry
  perspective*,
  J. Phys. A **44**, 113001 (2011).
- B. A. Bernevig and T. L. Hughes,
  *Topological Insulators and Topological Superconductors* (Princeton,
  2013), Ch. 8.

## Underlying band model

- X.-L. Qi, Y.-S. Wu, and S.-C. Zhang,
  *Topological quantization of the spin Hall effect in two-dimensional
  paramagnetic semiconductors*,
  Phys. Rev. B **74**, 085308 (2006).
  Source of the two-band Hamiltonian
  $H(\mathbf k) = \sin k_x\sigma_x + \sin k_y\sigma_y + (m + \cos k_x + \cos k_y)\sigma_z$
  that we use to validate the FHS quantisation.

## In-project links

- TopEquations registry entry:
  https://rdm3dc.github.io/TopEquations/registry.html#eq-local-plaquette-flux-proxy-chern
- HAFC / EGATL solver context:
  https://github.com/RDM3DC/eq-egatl-hlatn-adaptiveruler
- Sister observable for adaptive rulers in the same lattice family:
  https://github.com/RDM3DC/eq-flat-adaptive-annular-capacity-law
