"""Hero animation for the Local Plaquette-Flux Proxy Chern README.

Render (preview):
    manim -pqm animations/proxy_chern_hero.py ProxyChernHero

Render high quality MP4:
    manim -qh animations/proxy_chern_hero.py ProxyChernHero

Render a GIF for the README (smaller, loops on GitHub):
    manim -qm --format=gif animations/proxy_chern_hero.py ProxyChernHero

Concept
-------
1. Title + equation reveal.
2. A 24x24 BZ-mesh of plaquettes flashes in, each tinted by its actual rho_p
   value (computed via the repo's plaquette_chern.py with eigh-based gauge),
   in a diagonal wave.
3. A running C_loc counter on the right ticks down plaquette-by-plaquette
   and lands exactly on -1.000 (FHS quantisation).
4. "DAMAGE" event: a central horizontal strip turns black; the masked
   plaquettes' rho_p values lift off and dissolve. The counter drops live
   to ~-0.452.
5. "HEAL" event: damage fades, counter returns to -1.000.
6. Closing frame: equation, result table.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Make the repo's simulation module importable.
_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "simulations"))

from manim import (  # noqa: E402
    BLACK,
    DOWN,
    LEFT,
    ORIGIN,
    RIGHT,
    UP,
    WHITE,
    Create,
    DecimalNumber,
    FadeIn,
    FadeOut,
    Indicate,
    LaggedStart,
    ManimColor,
    MathTex,
    Rectangle,
    Scene,
    Square,
    Tex,
    Text,
    VGroup,
    Write,
    color_gradient,
    config,
    interpolate_color,
    linear,
    smooth,
)

from plaquette_chern import (  # noqa: E402
    build_chern_field,
    central_strip_mask,
    proxy_chern,
)

# ----------------------------------------------------------------------------
# Render configuration: 16:9, dark background suitable for GitHub README.
# ----------------------------------------------------------------------------
config.background_color = "#0b1020"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

# Diverging palette (RdBu_r-ish) hand-picked for visibility on dark BG.
_PALETTE = [
    ManimColor("#3b6bff"),   # strong negative
    ManimColor("#7a9bff"),
    ManimColor("#c0d0ff"),
    ManimColor("#1a1f3a"),   # ~zero, near background
    ManimColor("#ffd0a8"),
    ManimColor("#ff8a5b"),
    ManimColor("#ff3b3b"),   # strong positive
]


def _color_for(value: float, vmax: float) -> ManimColor:
    """Map a real value in [-vmax, vmax] to a diverging palette colour."""
    if vmax <= 0:
        t = 0.5
    else:
        t = 0.5 * (1.0 + max(-1.0, min(1.0, value / vmax)))
    n = len(_PALETTE) - 1
    idx = t * n
    lo = int(np.floor(idx))
    hi = min(n, lo + 1)
    frac = idx - lo
    return interpolate_color(_PALETTE[lo], _PALETTE[hi], frac)


def _wave_order(L: int) -> list[tuple[int, int]]:
    """Diagonal-wave ordering of plaquette indices (looks great in a sweep)."""
    pairs = [(i, j) for i in range(L) for j in range(L)]
    pairs.sort(key=lambda ij: (ij[0] + ij[1], ij[0]))
    return pairs


# ----------------------------------------------------------------------------
# Scene
# ----------------------------------------------------------------------------

class ProxyChernHero(Scene):
    def construct(self) -> None:
        # ---- 1. Pre-compute physics ---------------------------------------
        L = 24
        m = -1.5
        rho_p, _Ux, _Uy = build_chern_field(L=L, m=m)
        c_full = proxy_chern(rho_p)
        mask = central_strip_mask(L, half_width=3)
        c_damaged = proxy_chern(rho_p, mask=mask)
        vmax = float(np.max(np.abs(rho_p))) or 1.0e-6

        # Sign convention may give -1; show |C| in the storytelling but
        # display the signed live counter so the math is honest.
        order = _wave_order(L)

        # ---- 2. Title + equation -----------------------------------------
        title = Text(
            "Local Plaquette-Flux Proxy Chern",
            font="Inter",
            weight="BOLD",
            color=WHITE,
        ).scale(0.7).to_edge(UP, buff=0.4)

        equation = MathTex(
            r"\rho_p \;=\; \tfrac{1}{2\pi}\,\arg\!\Big("
            r"g_{p,0}\,g_{p,1}\,\overline{g_{p,2}}\,\overline{g_{p,3}}\Big)"
            r"\,,\quad C_{\mathrm{loc}} \;=\; \sum_p \rho_p",
            color=WHITE,
        ).scale(0.75).next_to(title, DOWN, buff=0.35)

        self.play(Write(title, run_time=1.0))
        self.play(Write(equation, run_time=1.6))
        self.wait(0.4)

        # ---- 3. Build lattice --------------------------------------------
        lattice_size = 5.0   # in scene units
        cell = lattice_size / L
        lattice_origin_x = -3.6
        lattice_origin_y = -lattice_size / 2 - 0.6

        squares: list[list[Square]] = [[None] * L for _ in range(L)]  # type: ignore
        for i in range(L):
            for j in range(L):
                sq = Square(side_length=cell)
                sq.set_stroke(width=0)
                sq.set_fill(_color_for(0.0, vmax), opacity=0.0)
                sq.move_to(
                    [
                        lattice_origin_x + (i + 0.5) * cell,
                        lattice_origin_y + (j + 0.5) * cell,
                        0.0,
                    ]
                )
                squares[i][j] = sq

        lattice_group = VGroup(*(squares[i][j] for i in range(L) for j in range(L)))

        # Frame around the lattice
        frame = Rectangle(
            width=lattice_size + 0.18,
            height=lattice_size + 0.18,
            stroke_color=WHITE,
            stroke_opacity=0.55,
            stroke_width=2.0,
        ).move_to(
            [lattice_origin_x + lattice_size / 2,
             lattice_origin_y + lattice_size / 2,
             0.0]
        )
        bz_label = Tex(r"Brillouin zone, $L = 24$", color=WHITE).scale(0.55)
        bz_label.next_to(frame, DOWN, buff=0.15)

        # Side panel: live counter + scenario label + ideal reference
        panel_x = 4.0
        panel_top_y = -0.4

        scenario_label = Text("HEALTHY", font="Inter", weight="BOLD", color="#9be7a3").scale(0.55)
        scenario_label.move_to([panel_x, panel_top_y, 0.0])

        c_label = MathTex(r"C_{\mathrm{loc}} = ", color=WHITE).scale(0.95)
        c_value = DecimalNumber(
            0.0,
            num_decimal_places=4,
            include_sign=True,
            color=WHITE,
        ).scale(1.4)
        c_row = VGroup(c_label, c_value).arrange(RIGHT, buff=0.18)
        c_row.move_to([panel_x, panel_top_y - 1.05, 0.0])

        # Reference markers
        ref_text = Tex(r"topological: $|C| = 1$", color="#aaa").scale(0.5)
        ref_text.move_to([panel_x, panel_top_y - 1.95, 0.0])

        kept_label = Tex(r"plaquettes kept: $576/576$", color="#cfd6e8").scale(0.5)
        kept_label.move_to([panel_x, panel_top_y - 2.4, 0.0])

        # Reveal lattice frame + side panel
        self.play(
            Create(frame, run_time=0.8),
            FadeIn(bz_label, run_time=0.6),
            FadeIn(scenario_label, run_time=0.6),
            Write(c_label, run_time=0.8),
            FadeIn(c_value, run_time=0.6),
            FadeIn(ref_text, run_time=0.6),
            FadeIn(kept_label, run_time=0.6),
        )

        # ---- 4. Wave fill: precompute target colors and partial sums -----
        partial = 0.0
        targets: list[tuple[Square, ManimColor, float]] = []
        for (i, j) in order:
            v = float(rho_p[i, j])
            partial += v
            targets.append((squares[i][j], _color_for(v, vmax), partial))

        # Animate in chunks so we can drive the counter alongside the fills.
        chunk_size = max(1, len(targets) // 60)  # ~60 frames of fills
        cumulative = 0.0
        for k in range(0, len(targets), chunk_size):
            chunk = targets[k : k + chunk_size]
            anims = []
            for sq, col, _ in chunk:
                anims.append(sq.animate.set_fill(col, opacity=1.0))
            cumulative = chunk[-1][2]
            self.play(
                LaggedStart(*anims, lag_ratio=0.05, run_time=0.18),
                c_value.animate.set_value(cumulative).set_color(WHITE),
                run_time=0.18,
            )

        # Snap to exact final value (handles float drift in display).
        self.play(
            c_value.animate.set_value(c_full),
            run_time=0.4,
        )
        self.play(Indicate(c_value, scale_factor=1.18, color="#9be7a3", run_time=0.8))
        self.wait(0.3)

        # ---- 5. DAMAGE event ---------------------------------------------
        damage_label = Text("CENTRAL STRIP DAMAGE",
                            font="Inter", weight="BOLD", color="#ff8a5b").scale(0.55)
        damage_label.move_to(scenario_label.get_center())

        damaged_squares = [squares[i][j] for i in range(L) for j in range(L) if not mask[i, j]]
        n_damaged = len(damaged_squares)
        n_kept = L * L - n_damaged
        kept_label_new = Tex(rf"plaquettes kept: ${n_kept}/{L*L}$", color="#cfd6e8").scale(0.5)
        kept_label_new.move_to(kept_label.get_center())

        # Build a "shockwave" overlay: a thin rectangle that sweeps across the
        # damaged strip vertically, leaving black behind it.
        cy_idx = L // 2
        strip_y0 = lattice_origin_y + (cy_idx - 3) * cell
        strip_y1 = lattice_origin_y + (cy_idx + 3) * cell
        strip_overlay = Rectangle(
            width=lattice_size,
            height=strip_y1 - strip_y0,
            fill_color=BLACK,
            fill_opacity=0.0,
            stroke_color="#ff8a5b",
            stroke_opacity=0.0,
            stroke_width=1.5,
        ).move_to(
            [lattice_origin_x + lattice_size / 2,
             (strip_y0 + strip_y1) / 2,
             0.0]
        )
        self.add(strip_overlay)

        damage_anims = [sq.animate.set_fill(BLACK, opacity=0.95) for sq in damaged_squares]
        self.play(
            FadeOut(scenario_label, run_time=0.3),
            FadeIn(damage_label, run_time=0.3),
            strip_overlay.animate.set_fill(BLACK, opacity=0.55).set_stroke(opacity=0.9),
            LaggedStart(*damage_anims, lag_ratio=0.01, run_time=1.2),
            c_value.animate.set_value(c_damaged).set_color("#ff8a5b"),
            kept_label.animate.become(kept_label_new),
            run_time=1.2,
        )
        self.play(Indicate(c_value, scale_factor=1.1, color="#ff8a5b", run_time=0.7))
        self.wait(0.5)

        # ---- 6. HEAL event ------------------------------------------------
        heal_label = Text("HEALED", font="Inter", weight="BOLD", color="#9be7a3").scale(0.55)
        heal_label.move_to(scenario_label.get_center())
        kept_label_full = Tex(rf"plaquettes kept: ${L*L}/{L*L}$",
                              color="#cfd6e8").scale(0.5).move_to(kept_label.get_center())

        # Restore the original fills for the damaged plaquettes.
        restore_anims = []
        for (i, j) in order:
            if not mask[i, j]:
                v = float(rho_p[i, j])
                restore_anims.append(
                    squares[i][j].animate.set_fill(_color_for(v, vmax), opacity=1.0)
                )
        self.play(
            FadeOut(damage_label, run_time=0.3),
            FadeIn(heal_label, run_time=0.3),
            strip_overlay.animate.set_fill(BLACK, opacity=0.0).set_stroke(opacity=0.0),
            LaggedStart(*restore_anims, lag_ratio=0.005, run_time=1.4),
            c_value.animate.set_value(c_full).set_color("#9be7a3"),
            kept_label.animate.become(kept_label_full),
            run_time=1.4,
        )
        self.play(Indicate(c_value, scale_factor=1.15, color="#9be7a3", run_time=0.7))
        self.wait(0.4)

        # ---- 7. Closing badge --------------------------------------------
        badge = VGroup(
            Tex(r"FHS-quantised: $\sum_p \rho_p \in \mathbb{Z}$",
                color=WHITE).scale(0.6),
            Tex(r"healthy $|C_{\mathrm{loc}}| = 1.000000$",
                color="#9be7a3").scale(0.55),
            Tex(r"damage drops it to $|C_{\mathrm{loc}}| \approx 0.452$",
                color="#ff8a5b").scale(0.55),
        ).arrange(DOWN, buff=0.16, aligned_edge=LEFT)
        badge.move_to([panel_x, panel_top_y - 3.25, 0.0])

        self.play(FadeIn(badge, run_time=1.0))
        self.wait(1.6)
