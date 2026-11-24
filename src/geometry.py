"""Root-network geometry and auditable scaling metrics.

The network is a *2-D projection* of a branching architecture. Cylindrical
surface/volume calculations are deliberately kept separate from the projected
box-counting slope so a 3-D electrochemical fractal dimension is not confused
with a 2-D drawing.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import math

import numpy as np


@dataclass(frozen=True)
class Segment:
    x0: float
    y0: float
    x1: float
    y1: float
    radius: float
    depth: int
    path_start: float

    @property
    def length(self) -> float:
        return math.hypot(self.x1 - self.x0, self.y1 - self.y0)

    @property
    def path_midpoint(self) -> float:
        return self.path_start + 0.5 * self.length


@dataclass
class RootNetwork:
    segments: list[Segment] = field(default_factory=list)

    def sampled_points(self, samples_per_unit: float = 40.0) -> np.ndarray:
        if not self.segments:
            raise ValueError("network has no segments")
        if samples_per_unit <= 0:
            raise ValueError("samples_per_unit must be positive")
        pts = []
        for s in self.segments:
            n = max(2, int(math.ceil(s.length * samples_per_unit)))
            pts.append(np.column_stack([
                np.linspace(s.x0, s.x1, n),
                np.linspace(s.y0, s.y1, n),
            ]))
        return np.vstack(pts)

    @property
    def points(self) -> np.ndarray:
        """Compatibility alias for the default projected point cloud."""
        return self.sampled_points()


def _validate_growth(depth: int, n_children: int, length: float,
                     length_shrink: float, radius0: float,
                     radius_shrink: float, jitter: float) -> None:
    if depth < 1:
        raise ValueError("depth must be >= 1")
    if n_children < 1:
        raise ValueError("n_children must be >= 1")
    if length <= 0 or radius0 <= 0:
        raise ValueError("length and radius0 must be positive")
    if not (0 < length_shrink < 1):
        raise ValueError("length_shrink must lie in (0, 1)")
    if not (0 < radius_shrink < 1):
        raise ValueError("radius_shrink must lie in (0, 1)")
    if jitter < 0:
        raise ValueError("jitter must be non-negative")


def grow(depth: int = 7, n_children: int = 3, length: float = 1.0,
         length_shrink: float = 0.65, radius0: float = 0.06,
         radius_shrink: float = 0.65, spread_deg: float = 55.0,
         jitter: float = 0.12, seed: int = 3) -> RootNetwork:
    """Grow a 2-D projected branching network.

    ``length_shrink`` controls centreline length; ``radius_shrink`` controls
    the cylindrical branch radius used only for proxy surface/volume metrics.
    Keeping them separate avoids the old ambiguity where a plotting linewidth
    was also treated as a physical electrode width.
    """
    _validate_growth(depth, n_children, length, length_shrink, radius0,
                     radius_shrink, jitter)
    rng = np.random.default_rng(seed)
    net = RootNetwork()

    def branch(x: float, y: float, angle: float, seg_len: float,
               radius: float, d: int, path_start: float) -> None:
        x1 = x + seg_len * math.sin(angle)
        y1 = y - seg_len * math.cos(angle)
        net.segments.append(Segment(x, y, x1, y1, radius, d, path_start))
        if d >= depth:
            return
        for i in range(n_children):
            frac = (i - (n_children - 1) / 2) / max(n_children - 1, 1)
            child_angle = (angle + math.radians(spread_deg) * 2 * frac
                           + rng.normal(0, jitter))
            child_len = seg_len * length_shrink * (1 + rng.normal(0, jitter / 2))
            if child_len <= 0:
                # Extremely large jitter could otherwise create an invalid
                # branch. The default jitter is far from this regime.
                child_len = np.finfo(float).eps
            branch(x1, y1, child_angle, child_len,
                   radius * radius_shrink, d + 1,
                   path_start + seg_len)

    branch(0.0, 0.0, 0.0, length, radius0, 1, 0.0)
    return net


def total_length(net: RootNetwork) -> float:
    return float(sum(s.length for s in net.segments))


def lateral_surface_area(net: RootNetwork) -> float:
    """Cylindrical side-area proxy, ignoring end caps."""
    return float(sum(2.0 * math.pi * s.radius * s.length
                     for s in net.segments))


def material_volume(net: RootNetwork) -> float:
    """Cylindrical material-volume proxy."""
    return float(sum(math.pi * s.radius**2 * s.length
                     for s in net.segments))


def material_weighted_root_path(net: RootNetwork) -> float:
    """Mean network distance from the root to branch material.

    Each segment is weighted by its cylindrical material-volume proxy and is
    represented at its midpoint along the connected branch path.
    """
    weights = np.array([math.pi * s.radius**2 * s.length for s in net.segments])
    paths = np.array([s.path_midpoint for s in net.segments])
    if weights.sum() <= 0:
        raise ValueError("network material volume is zero")
    return float(np.average(paths, weights=weights))


def length_by_depth(net: RootNetwork) -> dict[int, float]:
    out: dict[int, float] = {}
    for s in net.segments:
        out[s.depth] = out.get(s.depth, 0.0) + s.length
    return dict(sorted(out.items()))


def generation_surface_volume(net: RootNetwork) -> dict[int, dict[str, float]]:
    out: dict[int, dict[str, float]] = {}
    for s in net.segments:
        row = out.setdefault(s.depth, {"length": 0.0, "surface": 0.0,
                                       "volume": 0.0})
        row["length"] += s.length
        row["surface"] += 2 * math.pi * s.radius * s.length
        row["volume"] += math.pi * s.radius**2 * s.length
    return dict(sorted(out.items()))


def branching_multipliers(branching_factor: int, length_shrink: float,
                          radius_shrink: float) -> tuple[float, float]:
    """Ideal per-generation multipliers for cylindrical branches.

    For deterministic self-similar branching (ignoring overlap and end caps):
      lateral surface ~ b * lambda * gamma
      material volume ~ b * lambda * gamma^2
    """
    if branching_factor < 1:
        raise ValueError("branching_factor must be >= 1")
    if not (0 < length_shrink < 1 and 0 < radius_shrink < 1):
        raise ValueError("shrink factors must lie in (0, 1)")
    return (float(branching_factor * length_shrink * radius_shrink),
            float(branching_factor * length_shrink * radius_shrink**2))


def surface_efficient_radius_window(branching_factor: int,
                                    length_shrink: float) -> tuple[float, float]:
    """Radius-shrink interval where surface grows but material added falls.

    Conditions are b*lambda*gamma > 1 and b*lambda*gamma^2 < 1.  The
    interval exists only when b*lambda > 1.
    """
    if branching_factor < 1 or not (0 < length_shrink < 1):
        raise ValueError("invalid branching factor or length shrink")
    b_lambda = branching_factor * length_shrink
    if b_lambda <= 1:
        raise ValueError("no such interval when branching_factor*length_shrink <= 1")
    return float(1.0 / b_lambda), float(1.0 / math.sqrt(b_lambda))


def ideal_similarity_dimension(branching_factor: int,
                               length_shrink: float) -> float:
    """Ideal similarity dimension log(b)/log(1/lambda).

    This is an IFS-style geometric scaling quantity, not an electrochemical
    De-Gennes surface dimension. If it exceeds the embedding dimension, a
    non-overlapping self-similar construction cannot realise it there.
    """
    if branching_factor < 1:
        raise ValueError("branching_factor must be >= 1")
    if not (0 < length_shrink < 1):
        raise ValueError("length_shrink must lie in (0, 1)")
    return float(math.log(branching_factor) / math.log(1.0 / length_shrink))


def rescale_radii(net: RootNetwork, factor: float) -> RootNetwork:
    if factor <= 0:
        raise ValueError("factor must be positive")
    return RootNetwork([replace(s, radius=s.radius * factor)
                        for s in net.segments])


def box_counting_dimension(points: np.ndarray, n_scales: int = 10
                           ) -> tuple[float, np.ndarray, np.ndarray]:
    """Finite-scale projected box-counting slope for a 2-D point cloud.

    This is intentionally called a *finite-scale slope* in the public report:
    a finite stochastic branching tree need not possess a stable asymptotic
    fractal dimension over the sampled scales.
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2 or len(pts) < 2:
        raise ValueError("points must have shape (n, 2) with n >= 2")
    if n_scales < 4:
        raise ValueError("n_scales must be >= 4")
    if not np.isfinite(pts).all():
        raise ValueError("points must be finite")
    pts = pts - pts.min(axis=0)
    span = float(np.ptp(pts, axis=0).max())
    if span <= 0:
        raise ValueError("point cloud must span a non-zero distance")
    pts = pts / span
    epsilons = np.logspace(-0.5, -2.2, n_scales)
    counts = []
    for eps in epsilons:
        boxes = np.unique(np.floor(pts / eps).astype(int), axis=0)
        counts.append(len(boxes))
    counts = np.asarray(counts, dtype=float)
    slope, _ = np.polyfit(np.log(1 / epsilons), np.log(counts), 1)
    return float(slope), epsilons, counts
