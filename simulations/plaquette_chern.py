"""Local plaquette-flux proxy Chern: core implementation.

Implements the equation registered as
``eq-local-plaquette-flux-proxy-chern``::

    rho_p     = (1 / 2*pi) * arg( g_{p,0} g_{p,1} conj(g_{p,2}) conj(g_{p,3}) )
    C_loc     = sum_p rho_p

The complex bond responses ``g_{p,k}`` are interpreted as oriented link
variables on a square plaquette lattice. With link variables built from
overlaps of a smooth Bloch eigenstate this reduces to the well-known
Fukui-Hatsugai-Suzuki (FHS) lattice Chern formula, which gives a finite
proxy for the first Chern number on a discretised Brillouin-zone torus.

We use this concrete construction so the equation can be exercised on a
real physical model (a two-band Chern insulator), then we expose generic
hooks so arbitrary user-supplied bond responses (HAFC / EGATL solver
output) plug into the same observable.

References
----------
Fukui, Hatsugai, Suzuki, J. Phys. Soc. Jpn. 74, 1674 (2005).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# Band model used to generate physically meaningful link variables.
# ---------------------------------------------------------------------------

def two_band_lower_eigenstate(kx: np.ndarray, ky: np.ndarray, m: float) -> np.ndarray:
    """Lower-band eigenstate of a minimal Chern-insulator Hamiltonian.

    H(k) = sin(kx) sigma_x + sin(ky) sigma_y + (m + cos(kx) + cos(ky)) sigma_z

    Returns an array of shape ``(*kx.shape, 2)`` containing a unit-norm
    eigenvector of the lower band at each (kx, ky).

    The eigenvector is obtained per-site by ``numpy.linalg.eigh``. This
    yields a *non-smooth* gauge across the lattice (eigenvectors at
    neighbouring sites may carry arbitrary relative phase), which is
    exactly the regime in which the Fukui-Hatsugai-Suzuki construction
    delivers an integer Chern number on a finite mesh: large per-link
    phase jumps are absorbed by ``arg`` wrapping into (-pi, pi].

    For ``-2 < m < 0`` the lower band carries Chern number C = +1; for
    ``m < -2`` or ``m > 0`` the band is trivial.
    """
    if kx.shape != ky.shape:
        raise ValueError("kx and ky must have the same shape")

    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    sy = np.array([[0.0, -1j], [1j, 0.0]], dtype=np.complex128)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)

    dx = np.sin(kx)
    dy = np.sin(ky)
    dz = m + np.cos(kx) + np.cos(ky)

    H = (
        dx[..., None, None] * sx
        + dy[..., None, None] * sy
        + dz[..., None, None] * sz
    )
    # eigh returns eigenvalues in ascending order; column 0 is the lower band.
    _w, V = np.linalg.eigh(H)
    psi = V[..., :, 0]  # shape (..., 2)
    return psi


# ---------------------------------------------------------------------------
# Link variables and plaquette flux.
# ---------------------------------------------------------------------------

def link_variables(psi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute U-link variables on a periodic square mesh.

    Parameters
    ----------
    psi : ndarray, shape (Lx, Ly, n)
        Eigenstate at each mesh site, with periodic indexing in (x, y).

    Returns
    -------
    Ux, Uy : ndarray, shape (Lx, Ly), complex unit modulus
        Link variables in the +x and +y directions, ``U_mu(n) = <psi(n) | psi(n + mu)> / |...|``.
    """
    if psi.ndim != 3:
        raise ValueError("psi must have shape (Lx, Ly, n)")

    psi_xp = np.roll(psi, -1, axis=0)
    psi_yp = np.roll(psi, -1, axis=1)

    inner_x = np.einsum("ijk,ijk->ij", psi.conj(), psi_xp)
    inner_y = np.einsum("ijk,ijk->ij", psi.conj(), psi_yp)

    Ux = inner_x / (np.abs(inner_x) + 1.0e-30)
    Uy = inner_y / (np.abs(inner_y) + 1.0e-30)
    return Ux, Uy


def plaquette_g(Ux: np.ndarray, Uy: np.ndarray) -> np.ndarray:
    """Pack link variables into the four oriented bond responses per plaquette.

    Returns array of shape ``(Lx, Ly, 4)`` where the four bonds are taken
    in the canonical CCW order around the plaquette anchored at site n:

        g_{p,0} = U_x(n)             (n  -> n + x)
        g_{p,1} = U_y(n + x)         (n + x -> n + x + y)
        g_{p,2} = U_x(n + y)         (n + y -> n + x + y)   [reversed in product]
        g_{p,3} = U_y(n)             (n  -> n + y)          [reversed in product]
    """
    g0 = Ux
    g1 = np.roll(Uy, -1, axis=0)
    g2 = np.roll(Ux, -1, axis=1)
    g3 = Uy
    return np.stack([g0, g1, g2, g3], axis=-1)


def plaquette_flux(g: np.ndarray) -> np.ndarray:
    """Apply the registered equation: rho_p = arg(g0 g1 conj(g2) conj(g3)) / (2 pi).

    ``g`` has shape ``(..., 4)``. The result is a real array with the
    same leading shape, taking values in ``(-1/2, 1/2]``.
    """
    if g.shape[-1] != 4:
        raise ValueError("g must have last dimension 4")
    prod = g[..., 0] * g[..., 1] * np.conj(g[..., 2]) * np.conj(g[..., 3])
    return np.angle(prod) / (2.0 * np.pi)


def proxy_chern(rho_p: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Sum the local plaquette flux into the proxy Chern number.

    Parameters
    ----------
    rho_p : ndarray
        Local flux density per plaquette.
    mask : ndarray of bool, optional
        Boolean mask the same shape as ``rho_p``. ``True`` keeps the
        plaquette in the sum; ``False`` drops it (modelling damage that
        has rendered the local plaquette product undefined).
    """
    if mask is None:
        return float(rho_p.sum())
    if mask.shape != rho_p.shape:
        raise ValueError("mask shape mismatch")
    return float(rho_p[mask].sum())


# ---------------------------------------------------------------------------
# High-level scenarios.
# ---------------------------------------------------------------------------

@dataclass
class ChernResult:
    rho_p: np.ndarray            # (Lx, Ly) local flux density
    C_loc: float                 # summed proxy Chern
    mask: np.ndarray             # (Lx, Ly) bool plaquette mask used
    label: str                   # human-readable scenario label
    L: int                       # mesh size
    m: float                     # band-mass parameter


def build_chern_field(L: int = 48, m: float = -1.5) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Discretise the Brillouin zone and compute (rho_p, Ux, Uy)."""
    k1 = np.linspace(-np.pi, np.pi, L, endpoint=False)
    kx, ky = np.meshgrid(k1, k1, indexing="ij")
    psi = two_band_lower_eigenstate(kx, ky, m)
    Ux, Uy = link_variables(psi)
    g = plaquette_g(Ux, Uy)
    rho_p = plaquette_flux(g)
    return rho_p, Ux, Uy


def central_strip_mask(L: int, half_width: int = 4) -> np.ndarray:
    """Mask out a horizontal strip of plaquettes through the middle."""
    mask = np.ones((L, L), dtype=bool)
    cy = L // 2
    mask[:, cy - half_width : cy + half_width] = False
    return mask


def top_edge_mask(L: int, depth: int = 2) -> np.ndarray:
    """Mask out a thin band of plaquettes along the top edge."""
    mask = np.ones((L, L), dtype=bool)
    mask[:, -depth:] = False
    return mask


def run_scenario(label: str, mask: np.ndarray | None, L: int = 48, m: float = -1.5) -> ChernResult:
    rho_p, _Ux, _Uy = build_chern_field(L=L, m=m)
    if mask is None:
        mask = np.ones_like(rho_p, dtype=bool)
    C = proxy_chern(rho_p, mask=mask)
    return ChernResult(rho_p=rho_p, C_loc=C, mask=mask, label=label, L=L, m=m)


__all__ = [
    "ChernResult",
    "build_chern_field",
    "central_strip_mask",
    "link_variables",
    "plaquette_flux",
    "plaquette_g",
    "proxy_chern",
    "run_scenario",
    "top_edge_mask",
    "two_band_lower_eigenstate",
]
