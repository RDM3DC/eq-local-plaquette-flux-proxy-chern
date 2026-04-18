# Derivation: Local Plaquette-Flux Proxy Chern

**Equation ID:** `eq-local-plaquette-flux-proxy-chern`

$$
\rho_p \;=\; \frac{1}{2\pi}\,\arg\!\Big( g_{p,0}\, g_{p,1}\, \overline{g_{p,2}}\, \overline{g_{p,3}} \Big),
\qquad
C_{\mathrm{loc}} \;=\; \sum_p \rho_p .
$$

---

## 1. Setting

Consider a square plaquette lattice $\Lambda$ — either a real-space lattice
of bonds (HAFC / EGATL solver output) or a discretised Brillouin-zone
torus of a band model. To each oriented bond $b$ on $\Lambda$ we attach a
**complex bond response** $g_b \in \mathbb{C}$. Each plaquette $p$ is
bounded by four oriented bonds; ordering them counter-clockwise around
the plaquette anchored at site $n$ we write

$$
\big(g_{p,0},\, g_{p,1},\, g_{p,2},\, g_{p,3}\big)
\;=\;
\big(g_{n\to n+\hat x},\; g_{n+\hat x \to n+\hat x+\hat y},\;
g_{n+\hat y \to n+\hat x+\hat y},\; g_{n \to n+\hat y}\big).
$$

The two trailing bonds appear *reversed* in the loop traversal, so the
oriented loop product around the plaquette is

$$
W_p \;=\; g_{p,0}\, g_{p,1}\, g_{p,2}^{-1}\, g_{p,3}^{-1}.
$$

When $|g_b| = 1$ — i.e. when $g_b$ is taken as the link variable
$U_b = \langle\psi_a|\psi_b\rangle / |\langle\psi_a|\psi_b\rangle|$ on
the bond $b = (a\to b)$ — we have $g^{-1} = \bar g$, and $W_p$ becomes
exactly the form in the registered equation:

$$
W_p \;=\; g_{p,0}\, g_{p,1}\, \overline{g_{p,2}}\, \overline{g_{p,3}}.
$$

## 2. Local plaquette flux density $\rho_p$

Define the local flux density by lifting $W_p$ to the principal branch
of the argument:

$$
\rho_p \;=\; \frac{1}{2\pi}\,\arg(W_p) \;\in\; \Big(-\tfrac12,\,\tfrac12\Big].
$$

This is **bounded** and **gauge-invariant** under per-site phase changes
$\psi_n \to e^{i\alpha_n}\psi_n$, because the four phases entering each
$W_p$ telescope:

$$
W_p \;\to\; e^{i(\alpha_{n+\hat x}-\alpha_n)} e^{i(\alpha_{n+\hat x+\hat y}-\alpha_{n+\hat x})}
e^{-i(\alpha_{n+\hat x+\hat y}-\alpha_{n+\hat y})} e^{-i(\alpha_{n+\hat y}-\alpha_n)} \;=\; W_p.
$$

So $\rho_p$ is a meaningful local observable: a *flux density per
plaquette* with no gauge ambiguity.

## 3. Summed proxy Chern $C_{\mathrm{loc}}$

The summed observable

$$
C_{\mathrm{loc}} \;=\; \sum_{p \in \Lambda}\rho_p
$$

is the **Fukui–Hatsugai–Suzuki (FHS) proxy Chern number** when $\Lambda$
is a closed Brillouin-zone torus and $g$'s are the link variables of a
single (or projected) Bloch band:

$$
C \;=\; \frac{1}{2\pi}\!\int_{\mathrm{BZ}} F(\mathbf k)\, d^2k
\;\xleftarrow{L\to\infty}\;
\sum_p \rho_p .
$$

**Quantisation.** When the eigenvectors $\psi_n$ are chosen by an
arbitrary local diagonalisation (e.g. `numpy.linalg.eigh`), the gauge
varies wildly between neighbouring sites and individual $\rho_p$ values
spread across $(-\tfrac12, \tfrac12]$. Yet on a **closed** torus

$$
\sum_p \rho_p \;\in\; \mathbb{Z},
$$

*exactly*, for any mesh size — because the global product
$\prod_p W_p$ telescopes to $1$, so the sum of plaquette arguments is a
multiple of $2\pi$. This is the FHS theorem and is what the test
`test_chern_quantised_in_topological_phase` exercises in
[`tests/test_plaquette_chern.py`](../tests/test_plaquette_chern.py).

## 4. Why a *proxy* on damaged finite lattices

On HAFC / EGATL lattices in the wild the situation differs:

- $\Lambda$ is **not** a closed torus (open boundaries, broken bonds).
- Some plaquettes lose the link variables that bound them (the bond
  responses $g$ become undefined or noise-dominated on damaged bonds).
- Bott-style projector postprocessing is dense and scales poorly.

In this regime $C_{\mathrm{loc}}$ is **not** an exactly quantised
invariant. It is a **scalable, local, computable proxy** with three
properties that justify its use as a damage-resilient topological signal:

1. **Locality.** $\rho_p$ depends only on the four bonds bounding $p$,
   so the sum can be restricted to any subset
   $\Lambda' \subseteq \Lambda$ via a boolean mask:
   $C_{\mathrm{loc}}(\Lambda') = \sum_{p\in\Lambda'}\rho_p$.
2. **Boundedness.** $|\rho_p| \le \tfrac12$ uniformly, so a single bad
   plaquette can shift $C_{\mathrm{loc}}$ by at most $\tfrac12$.
3. **Continuity in the bulk.** When the BZ-torus reference value is
   $C \in \{0, \pm 1, \pm 2, \dots\}$, smooth bulk perturbations leave
   $C_{\mathrm{loc}}$ close to that integer; topological transitions
   show up as $\mathcal O(1)$ jumps.

## 5. Sanity scenarios run in this repo

[`simulations/run_scenarios.py`](../simulations/run_scenarios.py) executes
three scenarios on a $48\times 48$ BZ mesh of the Qi–Wu–Zhang lattice
Chern insulator $H(\mathbf k) = \sin k_x \sigma_x + \sin k_y \sigma_y + (m + \cos k_x + \cos k_y)\sigma_z$ at $m = -1.5$:

| Scenario        | Plaquettes kept | $C_{\mathrm{loc}}$ |
|-----------------|-----------------|---------------------|
| `healthy`       | 2304 / 2304     | $\mathbf{-1.000000}$ (FHS-quantised) |
| `central_strip` | 1920 / 2304     | $-0.452$            |
| `top_edge`      | 2208 / 2304     | $-0.995$            |

Reproduce with::

    python simulations/run_scenarios.py

The sign is gauge-arbitrary; the magnitude $|C| = 1$ is the topological
content. Damage that excises plaquettes near the Berry-curvature
concentration (the diamond near the band-touching point) costs much more
$C_{\mathrm{loc}}$ than damage to the same number of plaquettes near the
edge of the BZ — exactly what one expects from a proxy that retains
local sensitivity to where the topological weight lives.

## 6. Assumptions (see registry entry)

- Plaquette bond ordering and orientation are fixed consistently across
  runs.
- The $g_{p,k}$ remain defined on monitored plaquettes; the plaquette
  product is interpreted on the principal branch of $\arg$ before
  summation.
- $C_{\mathrm{loc}}$ is used as a *scalable proxy* on damaged finite
  lattices, not claimed to be an exactly quantised Bott invariant in
  every finite-size regime.
- The same local phase convention and damage protocol are used when
  comparing healthy / central-strip / top-edge runs.
