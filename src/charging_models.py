"""Small electrochemical scaling relations used in the design-study rebuild.

These are illustrative scaling models, not a validated sodium-ion cell model.
The De-Gennes fractal dimension belongs to a rough *surface in 3-D* and is
kept conceptually separate from the 2-D projected branching geometry.
"""
from __future__ import annotations

import numpy as np


def degennes_charge(t: np.ndarray, D: float) -> np.ndarray:
    """Relative charge Q(t) ~ t^((3-D)/2), normalised to Q(1)=1.

    The conventional surface-roughness interpretation used here has 2 <= D <= 3.
    """
    t = np.asarray(t, dtype=float)
    if np.any(t <= 0) or not np.isfinite(t).all():
        raise ValueError("t must contain finite positive values")
    if not (2.0 <= D <= 3.0):
        raise ValueError("De-Gennes surface dimension D must lie in [2, 3]")
    return t ** ((3.0 - D) / 2.0)


def charging_time(generation, r0: float = 10e-6, gamma: float = 0.65,
                  eps: float = 0.4, q_plus: float = 3.0e9,
                  j: float = 20.0) -> np.ndarray:
    """Illustrative branch-axis charging time under the reaction-zone model.

    t_i = r0 * gamma^i * (1-eps) * q_plus / (2*j)
    """
    generation = np.asarray(generation, dtype=float)
    if np.any(generation < 0) or not np.isfinite(generation).all():
        raise ValueError("generation must be finite and non-negative")
    if r0 <= 0 or q_plus <= 0 or j <= 0:
        raise ValueError("r0, q_plus and j must be positive")
    if not (0 < gamma < 1):
        raise ValueError("gamma must lie in (0, 1)")
    if not (0 <= eps < 1):
        raise ValueError("eps must lie in [0, 1)")
    r_k = r0 * gamma ** generation
    return r_k * (1 - eps) * q_plus / (2 * j)


def charging_time_ratio(generation, gamma: float = 0.65) -> np.ndarray:
    """t_i/t_0 = gamma^i for the branch-local reaction-zone model."""
    generation = np.asarray(generation, dtype=float)
    if np.any(generation < 0):
        raise ValueError("generation must be non-negative")
    if not (0 < gamma < 1):
        raise ValueError("gamma must lie in (0, 1)")
    return gamma ** generation
