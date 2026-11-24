import numpy as np
import pytest

from charging_models import charging_time, charging_time_ratio, degennes_charge


def test_degennes_plane_recovers_sqrt_time_and_d3_is_flat():
    t = np.array([1.0, 4.0, 100.0])
    assert np.allclose(degennes_charge(t, 2.0), np.sqrt(t))
    assert np.allclose(degennes_charge(t, 3.0), 1.0)


def test_degennes_dimension_domain_is_explicit():
    with pytest.raises(ValueError):
        degennes_charge(np.array([1.0, 2.0]), 1.5)
    with pytest.raises(ValueError):
        degennes_charge(np.array([1.0, 2.0]), 3.2)


def test_finer_generations_charge_faster_geometrically():
    gens = np.arange(6)
    ratio = charging_time_ratio(gens, 0.65)
    assert np.allclose(ratio, 0.65**gens)
    assert ratio[5] == pytest.approx(0.1160290625)
    assert np.all(np.diff(charging_time(gens, gamma=0.65)) < 0)


def test_charging_parameter_validation():
    with pytest.raises(ValueError):
        charging_time([0, 1], gamma=1.0)
    with pytest.raises(ValueError):
        charging_time([-1], gamma=0.65)
    with pytest.raises(ValueError):
        charging_time([0], eps=1.0)
