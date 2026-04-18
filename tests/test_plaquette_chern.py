"""Tests for the local plaquette-flux proxy Chern observable."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulations"))

from plaquette_chern import (  # noqa: E402
    build_chern_field,
    central_strip_mask,
    plaquette_flux,
    plaquette_g,
    proxy_chern,
    top_edge_mask,
)


def test_plaquette_flux_zero_for_unit_links() -> None:
    L = 6
    Ux = np.ones((L, L), dtype=complex)
    Uy = np.ones((L, L), dtype=complex)
    g = plaquette_g(Ux, Uy)
    rho = plaquette_flux(g)
    assert np.allclose(rho, 0.0)
    assert abs(proxy_chern(rho)) < 1.0e-12


def test_plaquette_flux_branch_bounds() -> None:
    rng = np.random.default_rng(0)
    g = np.exp(1j * rng.uniform(-np.pi, np.pi, size=(20, 20, 4)))
    rho = plaquette_flux(g)
    assert rho.shape == (20, 20)
    assert (rho > -0.5).all() and (rho <= 0.5 + 1.0e-12).all()


def test_chern_quantised_in_topological_phase() -> None:
    rho_p, _Ux, _Uy = build_chern_field(L=64, m=-1.5)
    C = proxy_chern(rho_p)
    # Lower band of the |C|=1 Chern insulator phase. FHS rounds to an
    # integer on a closed BZ torus regardless of mesh size; the sign of
    # C depends on the (gauge-arbitrary) eigenvector convention chosen
    # at each site by eigh and on the loop orientation, so we check
    # quantised magnitude.
    assert abs(abs(C) - 1.0) < 1.0e-9


def test_chern_trivial_phase() -> None:
    rho_p, _Ux, _Uy = build_chern_field(L=64, m=-3.0)
    C = proxy_chern(rho_p)
    assert abs(C) < 1.0e-9


def test_damaged_strip_reduces_chern_magnitude() -> None:
    rho_p, _Ux, _Uy = build_chern_field(L=48, m=-1.5)
    full = proxy_chern(rho_p)
    strip_mask = central_strip_mask(48, half_width=4)
    damaged = proxy_chern(rho_p, mask=strip_mask)
    # Strip removes nontrivial flux density => proxy Chern shifts away from 1.
    assert abs(damaged - full) > 1.0e-3
    # And does not exceed the healthy reference in magnitude.
    assert abs(damaged) <= abs(full) + 1.0e-9


def test_top_edge_smaller_perturbation_than_central_strip() -> None:
    rho_p, _Ux, _Uy = build_chern_field(L=48, m=-1.5)
    full = proxy_chern(rho_p)
    strip = proxy_chern(rho_p, mask=central_strip_mask(48, half_width=4))
    edge = proxy_chern(rho_p, mask=top_edge_mask(48, depth=2))
    # A 2-row top-edge cut perturbs C_loc less than a 8-row central strip.
    assert abs(edge - full) < abs(strip - full)
