"""Run the three canonical scenarios for the proxy Chern observable.

Outputs
-------
data/proxy_chern_summary.csv   - one row per scenario (label, L, m, C_loc, n_plaquettes_kept)
data/<label>_rho_p.csv         - per-plaquette flux density on the lattice
images/<label>_rho_p.png       - heatmap of rho_p with damaged region marked
images/c_loc_bar.png           - bar chart comparing C_loc across scenarios

Run from the repository root::

    python simulations/run_scenarios.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from plaquette_chern import (
    central_strip_mask,
    proxy_chern,
    run_scenario,
    top_edge_mask,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
IMAGES = ROOT / "images"


def _ensure_dirs() -> None:
    DATA.mkdir(exist_ok=True)
    IMAGES.mkdir(exist_ok=True)


def _save_rho_csv(label: str, rho_p: np.ndarray) -> Path:
    path = DATA / f"{label}_rho_p.csv"
    L = rho_p.shape[0]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["px", "py", "rho_p"])
        for i in range(L):
            for j in range(L):
                w.writerow([i, j, f"{rho_p[i, j]:.10e}"])
    return path


def _plot_rho(label: str, rho_p: np.ndarray, mask: np.ndarray, c_loc: float) -> Path:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.4, 4.8))
    vmax = max(1.0e-12, float(np.max(np.abs(rho_p))))
    im = ax.imshow(
        rho_p.T,
        origin="lower",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    # overlay the damaged plaquettes (mask == False) as hatched
    damaged = ~mask
    if damaged.any():
        ax.contourf(
            np.arange(rho_p.shape[0]),
            np.arange(rho_p.shape[1]),
            damaged.T.astype(float),
            levels=[0.5, 1.5],
            colors="none",
            hatches=["////"],
        )
    ax.set_title(f"{label}\n$C_{{\\mathrm{{loc}}}} = {c_loc:+.4f}$")
    ax.set_xlabel("plaquette x")
    ax.set_ylabel("plaquette y")
    fig.colorbar(im, ax=ax, label=r"$\rho_p$")
    fig.tight_layout()
    out = IMAGES / f"{label}_rho_p.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _plot_bar(results: list) -> Path:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    labels = [r.label for r in results]
    values = [r.C_loc for r in results]
    bars = ax.bar(labels, values, color=["#2a7", "#d62", "#46c"])
    ax.axhline(1.0, color="black", lw=0.6, ls="--", alpha=0.7)
    ax.axhline(-1.0, color="black", lw=0.6, ls="--", alpha=0.7,
               label=r"$|C| = 1$ topological")
    ax.axhline(0.0, color="grey", lw=0.4)
    ax.set_ylabel(r"$C_{\mathrm{loc}} = \sum_p \rho_p$")
    ax.set_title("Proxy Chern under plaquette damage")
    ax.legend(loc="lower left")
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val + 0.02 * (1 if val >= 0 else -1),
            f"{val:+.3f}",
            ha="center",
            va="bottom" if val >= 0 else "top",
            fontsize=9,
        )
    fig.tight_layout()
    out = IMAGES / "c_loc_bar.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> None:
    _ensure_dirs()
    L = 48
    m = -1.5

    healthy = run_scenario("healthy", mask=None, L=L, m=m)
    strip = run_scenario(
        "central_strip", mask=central_strip_mask(L, half_width=4), L=L, m=m
    )
    edge = run_scenario(
        "top_edge", mask=top_edge_mask(L, depth=2), L=L, m=m
    )
    results = [healthy, strip, edge]

    summary = DATA / "proxy_chern_summary.csv"
    with summary.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "label",
                "L",
                "m",
                "n_plaquettes_total",
                "n_plaquettes_kept",
                "C_loc",
                "C_loc_full",
            ]
        )
        for r in results:
            full = proxy_chern(r.rho_p)  # ignoring mask, healthy reference
            w.writerow(
                [
                    r.label,
                    r.L,
                    r.m,
                    r.rho_p.size,
                    int(r.mask.sum()),
                    f"{r.C_loc:.6f}",
                    f"{full:.6f}",
                ]
            )
            _save_rho_csv(r.label, r.rho_p)
            _plot_rho(r.label, r.rho_p, r.mask, r.C_loc)

    _plot_bar(results)

    print("=== Local Plaquette-Flux Proxy Chern ===")
    print(f"L = {L}, m = {m} (topological phase, |C| = 1; sign is gauge-dependent)")
    for r in results:
        print(
            f"  {r.label:<14s}  C_loc = {r.C_loc:+.6f}  "
            f"({int(r.mask.sum())}/{r.rho_p.size} plaquettes kept)"
        )
    print(f"summary -> {summary.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
