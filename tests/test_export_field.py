import numpy as np
import pytest

from export_field import density_field, field_diagnostics
from geometry import grow, material_volume


def test_density_field_conserves_material_to_roundoff():
    net = grow(depth=5, n_children=2, seed=1)
    field, _, axes = density_field(net, resolution=120)
    diag = field_diagnostics(net, field, axes)
    assert field.sum() == pytest.approx(material_volume(net), rel=1e-12, abs=1e-14)
    assert diag["relative_conservation_error"] < 1e-12


def test_raster_centroid_tracks_analytic_centroid_within_one_cell_diagonal():
    net = grow(depth=5, n_children=2, seed=2)
    field, _, axes = density_field(net, resolution=160)
    diag = field_diagnostics(net, field, axes)
    assert diag["centroid_error_cell_diagonals"] < 1.0


def test_material_integral_is_resolution_invariant():
    net = grow(depth=4, n_children=3, seed=4)
    vals = []
    for res in (80, 160):
        field, _, _ = density_field(net, resolution=res)
        vals.append(field.sum())
    assert vals[0] == pytest.approx(vals[1], rel=1e-12)


def test_density_field_validation():
    net = grow(depth=3, seed=0)
    with pytest.raises(ValueError):
        density_field(net, resolution=8)
