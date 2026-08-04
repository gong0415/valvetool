import math
from dataclasses import replace

from solenoid_model import magnetics
from solenoid_model.baseline import BASELINE_PARAMS


def test_reluctance_matches_air_gap_dominated_analytical():
    params = replace(BASELINE_PARAMS, mu_r_core=1e9)  # make core reluctance negligible
    gap = params.g0
    expected = gap / (magnetics.MU_0 * params.A_gap)
    assert math.isclose(magnetics.reluctance(gap, params), expected, rel_tol=1e-6)


def test_inductance_matches_air_gap_dominated_analytical():
    params = replace(BASELINE_PARAMS, mu_r_core=1e9)
    gap = params.g0
    expected = params.N_turns ** 2 * magnetics.MU_0 * params.A_gap / gap
    assert math.isclose(magnetics.inductance(gap, params), expected, rel_tol=1e-6)


def test_reluctance_includes_core_contribution():
    params = BASELINE_PARAMS
    gap = params.g0
    air_gap_only = gap / (magnetics.MU_0 * params.A_gap)
    assert magnetics.reluctance(gap, params) > air_gap_only


def test_dinductance_dgap_matches_numerical_derivative():
    params = BASELINE_PARAMS
    gap = params.g0
    h = gap * 1e-6
    numerical = (magnetics.inductance(gap + h, params) - magnetics.inductance(gap - h, params)) / (2 * h)
    analytical = magnetics.dinductance_dgap(gap, params)
    assert math.isclose(analytical, numerical, rel_tol=1e-6)


def test_dinductance_dgap_is_negative():
    params = BASELINE_PARAMS
    assert magnetics.dinductance_dgap(params.g0, params) < 0


def test_magnetic_force_matches_maxwell_stress_tensor():
    params = replace(BASELINE_PARAMS, mu_r_core=1e9)
    gap = params.g0
    i = 0.3
    coenergy_force = magnetics.magnetic_force_closing(gap, i, params)
    B = magnetics.flux_density(gap, i, params)
    maxwell_force = B ** 2 * params.A_gap / (2 * magnetics.MU_0)
    assert math.isclose(coenergy_force, maxwell_force, rel_tol=1e-6)


def test_magnetic_force_is_positive():
    params = BASELINE_PARAMS
    force = magnetics.magnetic_force_closing(params.g0, 0.3, params)
    assert force > 0


def test_saturation_check_flags_high_current():
    params = BASELINE_PARAMS
    assert magnetics.saturation_check(params.g0, 50.0, params) is True


def test_saturation_check_false_at_low_current():
    params = BASELINE_PARAMS
    assert magnetics.saturation_check(params.g0, 1e-6, params) is False
