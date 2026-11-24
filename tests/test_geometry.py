import math

import numpy as np
import pytest

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


def test_straight_line_has_projected_box_slope_near_one():
    pts = np.column_stack([np.linspace(0, 1, 4000), np.zeros(4000)])
    dim, _, _ = box_counting_dimension(pts)
    assert 0.9 < dim < 1.1


def test_growth_segment_count_and_positive_geometry():
    net = grow(depth=4, n_children=2, seed=0)
    assert len(net.segments) == 2**4 - 1
    assert total_length(net) > 0
    assert material_volume(net) > 0
    assert lateral_surface_area(net) > 0
    assert all(s.radius > 0 and s.length > 0 for s in net.segments)


def test_similarity_dimension_is_not_the_projected_box_slope():
    d = ideal_similarity_dimension(3, 0.65)
    assert d == pytest.approx(2.55, abs=0.01)
    assert d > 2.0  # impossible as a non-overlapping planar similarity dimension


def test_default_surface_volume_multipliers_define_efficiency_window():
    surface, volume = branching_multipliers(3, 0.65, 0.65)
    low, high = surface_efficient_radius_window(3, 0.65)
    assert surface == pytest.approx(1.2675)
    assert volume == pytest.approx(0.823875)
    assert low < 0.65 < high
    assert surface > 1 and volume < 1


def test_deterministic_generation_ratios_match_analytic_multipliers():
    net = grow(depth=5, n_children=3, length_shrink=0.65,
               radius_shrink=0.65, jitter=0.0, seed=0)
    by = generation_surface_volume(net)
    surface_ratio = by[2]["surface"] / by[1]["surface"]
    volume_ratio = by[2]["volume"] / by[1]["volume"]
    expected_surface, expected_volume = branching_multipliers(3, 0.65, 0.65)
    assert surface_ratio == pytest.approx(expected_surface)
    assert volume_ratio == pytest.approx(expected_volume)


def test_material_weighted_root_path_is_within_network_extent():
    net = grow(depth=5, seed=1)
    mean_path = material_weighted_root_path(net)
    max_path = max(s.path_start + s.length for s in net.segments)
    assert 0 < mean_path < max_path


def test_projected_box_slope_increases_with_added_deterministic_depth():
    shallow = grow(depth=4, n_children=3, jitter=0.0, seed=0)
    deep = grow(depth=8, n_children=3, jitter=0.0, seed=0)
    d_shallow, _, _ = box_counting_dimension(shallow.points)
    d_deep, _, _ = box_counting_dimension(deep.points)
    assert d_deep > d_shallow + 0.2


def test_growth_validation():
    with pytest.raises(ValueError):
        grow(depth=0)
    with pytest.raises(ValueError):
        grow(length_shrink=1.2)
    with pytest.raises(ValueError):
        grow(radius_shrink=0)
