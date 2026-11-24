"""Analyse projected branching geometry and simple electrochemical scalings."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

import charging_models as cm
from geometry import (
    box_counting_dimension,
    branching_multipliers,
    generation_surface_volume,
    grow,
    ideal_similarity_dimension,
    lateral_surface_area,
    material_volume,
    material_weighted_root_path,
    surface_efficient_radius_window,
    total_length,
)

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "reports" / "figures"


def render(net, path: Path) -> None:
    segs = [[(s.x0, s.y0), (s.x1, s.y1)] for s in net.segments]
    max_r = max(s.radius for s in net.segments)
    widths = [0.45 + 7.0 * s.radius / max_r for s in net.segments]
    depths = np.array([s.depth for s in net.segments])
    colors = plt.cm.YlOrRd(0.25 + 0.7 * depths / depths.max())

    fig, ax = plt.subplots(figsize=(8, 9))
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#0d1117")
    ax.add_collection(LineCollection(segs, linewidths=widths, colors=colors,
                                     capstyle="round"))
    ax.autoscale()
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Projected root-network electrode architecture",
                 color="w", fontsize=13, pad=12)
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=fig.get_facecolor())
    plt.close(fig)


def depth_study(depths, path: Path, **grow_kwargs):
    length_ratio, path_proxy, area_to_volume = [], [], []
    for d in depths:
        net = grow(depth=d, **grow_kwargs)
        trunk = next(s for s in net.segments if s.depth == 1)
        length_ratio.append(total_length(net) / trunk.length)
        path_proxy.append(material_weighted_root_path(net) / trunk.length)
        area_to_volume.append(lateral_surface_area(net) / material_volume(net))

    fig, ax1 = plt.subplots(figsize=(7.5, 4.8))
    ax1.plot(depths, length_ratio, "o-", lw=1.8,
             label="centreline length / trunk length")
    ax1.set_yscale("log")
    ax1.set(xlabel="Branching generations",
            ylabel="Centreline-length ratio (log)")
    ax2 = ax1.twinx()
    ax2.plot(depths, path_proxy, "s--", lw=1.8,
             label="material-weighted root path / trunk length")
    ax2.set_ylabel("Material-weighted root path ratio")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [ln.get_label() for ln in lines], loc="upper left")
    ax1.grid(alpha=0.3)
    ax1.set_title("Finite-network growth: more centreline, longer connected paths")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return length_ratio, path_proxy, area_to_volume


def scaling_window_figure(path: Path, branching_factor=3, length_shrink=0.65,
                          default_radius_shrink=0.65):
    gammas = np.linspace(0.30, 0.85, 250)
    surface = branching_factor * length_shrink * gammas
    volume = branching_factor * length_shrink * gammas**2
    low, high = surface_efficient_radius_window(branching_factor, length_shrink)
    s0, v0 = branching_multipliers(branching_factor, length_shrink,
                                    default_radius_shrink)

    fig, ax = plt.subplots(figsize=(7.3, 4.8))
    ax.plot(gammas, surface, label="lateral-surface multiplier / generation")
    ax.plot(gammas, volume, label="material-volume multiplier / generation")
    ax.axhline(1.0, ls="--", lw=1.0, alpha=0.65)
    ax.axvspan(low, high, alpha=0.16,
               label="surface grows while added material falls")
    ax.scatter([default_radius_shrink, default_radius_shrink], [s0, v0], zorder=3)
    ax.set(xlabel="Radius reduction factor γ",
           ylabel="Ideal per-generation multiplier",
           title="Cylindrical branching scalings, b=3 and length shrink=0.65")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return low, high, s0, v0


def box_slope_stability(path: Path, depths=range(4, 9), seeds=range(10),
                        **grow_kwargs):
    rows = []
    for d in depths:
        vals = []
        for seed in seeds:
            net = grow(depth=d, seed=seed, **grow_kwargs)
            slope, _, _ = box_counting_dimension(net.points)
            vals.append(slope)
            rows.append({"depth": int(d), "seed": int(seed), "slope": float(slope)})
        vals = np.asarray(vals)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for d in depths:
        vals = np.array([r["slope"] for r in rows if r["depth"] == d])
        x = np.full_like(vals, float(d))
        ax.scatter(x, vals, alpha=0.45, s=20)
        ax.plot([d - 0.16, d + 0.16], [np.median(vals)] * 2, lw=2.3)
    ax.set(xlabel="Finite branching depth",
           ylabel="Projected finite-scale box-counting slope",
           title="The projected box slope is not a depth-invariant material constant")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return rows


def charging_figures(figdir: Path):
    t = np.logspace(-1, 3, 200)
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for D in (2.0, 2.25, 2.5, 2.75, 3.0):
        ax.loglog(t, cm.degennes_charge(t, D), lw=1.8,
                  label=f"D = {D:.2f} (slope {(3-D)/2:.2f})")
    ax.set(xlabel="Time (arbitrary units)", ylabel="Relative charge",
           title="De-Gennes rough-surface scaling: Q ∝ t^((3-D)/2)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(figdir / "charge_scaling.png", dpi=150)
    plt.close(fig)

    gens = np.arange(0, 8)
    fig2, ax2 = plt.subplots(figsize=(7, 4.8))
    for gamma in (0.55, 0.65, 0.75):
        ax2.semilogy(gens, cm.charging_time_ratio(gens, gamma), "o-",
                     lw=1.8, label=f"radius reduction γ={gamma}")
    ax2.set(xlabel="Branching generation i",
            ylabel="Branch-local charge-time ratio tᵢ/t₀ (log)",
            title="Reaction-zone scaling: finer branch radii shorten local charge time")
    ax2.legend()
    ax2.grid(alpha=0.3, which="both")
    fig2.tight_layout()
    fig2.savefig(figdir / "charging_time.png", dpi=150)
    plt.close(fig2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", type=int, default=7)
    ap.add_argument("--children", type=int, default=3)
    ap.add_argument("--spread", type=float, default=55.0)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--length-shrink", type=float, default=0.65)
    ap.add_argument("--radius-shrink", type=float, default=0.65)
    args = ap.parse_args()

    FIGURES.mkdir(parents=True, exist_ok=True)
    kwargs = dict(n_children=args.children, spread_deg=args.spread,
                  length_shrink=args.length_shrink,
                  radius_shrink=args.radius_shrink)

    net = grow(depth=args.depth, seed=args.seed, **kwargs)
    render(net, FIGURES / "electrode_render.png")

    dim, eps, counts = box_counting_dimension(net.points)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(np.log(1 / eps), np.log(counts), "o-", lw=1.6)
    ax.set(xlabel="log(1/box size)", ylabel="log(occupied boxes)",
           title=f"Finite projected box-counting slope = {dim:.2f}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "box_counting.png", dpi=150)
    plt.close(fig)

    depths = list(range(3, args.depth + 1))
    length_ratio, path_proxy, area_to_volume = depth_study(
        depths, FIGURES / "depth_tradeoff.png", **kwargs)

    low, high, surf_mult, vol_mult = scaling_window_figure(
        FIGURES / "branching_scaling_window.png",
        branching_factor=args.children,
        length_shrink=args.length_shrink,
        default_radius_shrink=args.radius_shrink)

    stability_rows = box_slope_stability(
        FIGURES / "box_slope_stability.png",
        n_children=args.children,
        spread_deg=args.spread,
        length_shrink=args.length_shrink,
        radius_shrink=args.radius_shrink)

    charging_figures(FIGURES)

    sim_dim = ideal_similarity_dimension(args.children, args.length_shrink)
    depth4 = [r["slope"] for r in stability_rows if r["depth"] == 4]
    depth8 = [r["slope"] for r in stability_rows if r["depth"] == 8]
    generation5_ratio = float(cm.charging_time_ratio(5, args.radius_shrink))
    summary = {
        "scope": "2-D projected geometry plus separate illustrative electrochemical scalings",
        "generations": args.depth,
        "children_per_node": args.children,
        "n_segments": len(net.segments),
        "length_shrink": args.length_shrink,
        "radius_shrink": args.radius_shrink,
        "projected_box_slope_default": round(dim, 3),
        "projected_box_slope_depth4_median_10seeds": round(float(np.median(depth4)), 3),
        "projected_box_slope_depth8_median_10seeds": round(float(np.median(depth8)), 3),
        "ideal_similarity_dimension_from_branch_count_and_length_shrink": round(sim_dim, 3),
        "exceeds_2d_embedding_dimension": bool(sim_dim > 2.0),
        "surface_multiplier_per_generation": round(surf_mult, 4),
        "material_volume_multiplier_per_generation": round(vol_mult, 4),
        "surface_efficient_radius_window": [round(low, 4), round(high, 4)],
        "generation5_local_charge_time_ratio": round(generation5_ratio, 4),
        "centreline_length_ratio_at_max_depth": round(length_ratio[-1], 3),
        "material_weighted_root_path_ratio_at_max_depth": round(path_proxy[-1], 3),
        "surface_to_volume_proxy_at_max_depth": round(area_to_volume[-1], 3),
    }
    out = ROOT / "reports" / "geometry_summary.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    (ROOT / "reports" / "box_slope_stability.json").write_text(
        json.dumps(stability_rows, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"Figures -> {FIGURES}")


if __name__ == "__main__":
    main()
