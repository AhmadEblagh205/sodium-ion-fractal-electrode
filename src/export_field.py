"""Rasterise a projected branching network into a conservative material field.

The field is a generic 2-D material-per-cell representation. It is *not* a
ready-made PyBaMM input: mapping a 2-D architecture into a particular battery
PDE model is a separate modelling step.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from geometry import grow, material_volume

ROOT = Path(__file__).resolve().parents[1]


def _analytic_centroid(net) -> np.ndarray:
    weights = np.array([np.pi * s.radius**2 * s.length for s in net.segments])
    mids = np.array([[(s.x0 + s.x1) / 2, (s.y0 + s.y1) / 2]
                     for s in net.segments])
    return np.average(mids, axis=0, weights=weights)


def density_field(net, resolution: int = 200, oversample: float = 1.5):
    """Conservative bilinear raster of cylindrical material-volume proxy.

    Each segment's volume ``pi*r^2*L`` is split among sampled centreline points
    and deposited bilinearly into neighbouring grid cells. The field sum is
    therefore conserved by construction; centroid error is also reported to
    test spatial fidelity rather than conservation alone.
    """
    if resolution < 16:
        raise ValueError("resolution must be >= 16")
    if oversample <= 0:
        raise ValueError("oversample must be positive")
    pts = np.array([v for s in net.segments for v in ((s.x0, s.y0),
                                                      (s.x1, s.y1))], dtype=float)
    span = np.ptp(pts, axis=0)
    pad = max(float(span.max()) * 0.02, 1e-6)
    lo = pts.min(axis=0) - pad
    hi = pts.max(axis=0) + pad

    field = np.zeros((resolution, resolution), dtype=float)
    dx = (hi[0] - lo[0]) / resolution
    dy = (hi[1] - lo[1]) / resolution
    ds_target = min(dx, dy) / oversample

    for s in net.segments:
        n = max(2, int(np.ceil(s.length / ds_target)))
        # Midpoint sampling avoids duplicating branch endpoints.
        u = (np.arange(n) + 0.5) / n
        xs = s.x0 + u * (s.x1 - s.x0)
        ys = s.y0 + u * (s.y1 - s.y0)
        mass = np.pi * s.radius**2 * s.length / n

        gx = (xs - lo[0]) / dx - 0.5
        gy = (ys - lo[1]) / dy - 0.5
        i0 = np.floor(gx).astype(int)
        j0 = np.floor(gy).astype(int)
        fx = gx - i0
        fy = gy - j0

        for di, wx in ((0, 1 - fx), (1, fx)):
            for dj, wy in ((0, 1 - fy), (1, fy)):
                ii = np.clip(i0 + di, 0, resolution - 1)
                jj = np.clip(j0 + dj, 0, resolution - 1)
                np.add.at(field, (jj, ii), mass * wx * wy)

    x_centres = lo[0] + (np.arange(resolution) + 0.5) * dx
    y_centres = lo[1] + (np.arange(resolution) + 0.5) * dy
    return field, (lo, hi), (x_centres, y_centres)


def field_diagnostics(net, field, axes) -> dict:
    x, y = axes
    yy, xx = np.meshgrid(y, x, indexing="ij")
    total = float(field.sum())
    raster_centroid = np.array([(field * xx).sum() / total,
                                (field * yy).sum() / total])
    analytic = _analytic_centroid(net)
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    cell_diag = float(np.hypot(dx, dy))
    return {
        "analytic_material_volume": float(material_volume(net)),
        "field_integral": total,
        "relative_conservation_error": float(abs(total - material_volume(net)) /
                                             material_volume(net)),
        "analytic_centroid": analytic.tolist(),
        "raster_centroid": raster_centroid.tolist(),
        "centroid_error": float(np.linalg.norm(raster_centroid - analytic)),
        "centroid_error_cell_diagonals": float(np.linalg.norm(raster_centroid - analytic) /
                                                cell_diag),
        "occupied_cell_fraction": float((field > 0).mean()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", type=int, default=7)
    ap.add_argument("--children", type=int, default=3)
    ap.add_argument("--resolution", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    net = grow(depth=args.depth, n_children=args.children, seed=args.seed)
    field, (lo, hi), axes = density_field(net, args.resolution)
    diag = field_diagnostics(net, field, axes)

    out_dir = ROOT / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(exist_ok=True)
    np.savez_compressed(out_dir / "density_field.npz", field=field,
                        x=axes[0], y=axes[1],
                        x_min=lo[0], y_min=lo[1], x_max=hi[0], y_max=hi[1])
    meta = {"resolution": args.resolution,
            "segments": len(net.segments),
            **{k: (round(v, 8) if isinstance(v, float) else v)
               for k, v in diag.items()}}
    (out_dir / "density_field.json").write_text(json.dumps(meta, indent=2))

    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    im = ax.imshow(field, origin="lower", cmap="cividis",
                   extent=[lo[0], hi[0], lo[1], hi[1]])
    fig.colorbar(im, ax=ax, label="cylindrical material-volume proxy / cell")
    ax.set(title="Conservative 2-D material field")
    fig.tight_layout()
    fig.savefig(out_dir / "figures" / "density_field.png", dpi=150)
    plt.close(fig)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
