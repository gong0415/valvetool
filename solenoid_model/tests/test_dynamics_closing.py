import math

import pytest

from solenoid_model import dynamics, flyback
from solenoid_model.baseline import BASELINE_PARAMS

RC_CIRCUIT = flyback.RCFlyback(R_snub=14655.9, C_snub=1e-7)  # zeta ~= 1.5


def test_diode_closing_matches_verified_baseline_timing():
    result = dynamics.simulate_closing(BASELINE_PARAMS, flyback.DiodeFlyback())
    assert math.isclose(result.t_release, 0.061402777555826515, rel_tol=1e-6)
    assert math.isclose(result.t_close, 0.06274842379361782, rel_tol=1e-6)


def test_zener_closing_matches_verified_baseline_timing():
    result = dynamics.simulate_closing(BASELINE_PARAMS, flyback.ZenerFlyback())
    assert math.isclose(result.t_release, 0.01696499443027462, rel_tol=1e-6)
    assert math.isclose(result.t_close, 0.017456210526240436, rel_tol=1e-6)


def test_rc_closing_matches_verified_baseline_timing():
    result = dynamics.simulate_closing(BASELINE_PARAMS, RC_CIRCUIT)
    assert math.isclose(result.t_release, 0.0002936551130574534, rel_tol=1e-5)
    assert math.isclose(result.t_close, 0.00047571971038468336, rel_tol=1e-5)


def test_zener_closes_faster_than_diode():
    t_diode = dynamics.simulate_closing(BASELINE_PARAMS, flyback.DiodeFlyback()).t_close
    t_zener = dynamics.simulate_closing(BASELINE_PARAMS, flyback.ZenerFlyback()).t_close
    assert t_zener < t_diode


def test_stroke_position_monotonically_decreases_to_zero():
    result = dynamics.simulate_closing(BASELINE_PARAMS, flyback.ZenerFlyback())
    x_t = result.sol_stroke.y[1]
    assert all(a >= b for a, b in zip(x_t, x_t[1:]))
    assert math.isclose(x_t[-1], 0.0, abs_tol=1e-9)


def test_diode_current_rises_during_stroke_due_to_flux_conservation():
    # The flyback clamp only controls the held-phase electrical decay. Once
    # released, gap increases from gap_min to g0 as the armature returns to
    # x=0, and inductance collapses ~5.8x over the stroke (L(gap_min)=2.4127 H,
    # L(g0)=0.4160 H). Coil flux linkage L*i is approximately conserved on the
    # stroke's fast timescale, so current actually rises during the release
    # stroke rather than continuing to decay toward zero -- it does not reach
    # 0 by the time the valve closes.
    result = dynamics.simulate_closing(BASELINE_PARAMS, flyback.DiodeFlyback())
    i_release = result.sol_stroke.y[0][0]
    i_close = result.sol_stroke.y[0][-1]
    assert i_close > i_release
    assert math.isclose(i_close, 0.19429423446295552, rel_tol=1e-5)


def test_simulate_closing_rejects_never_releasing_circuit():
    with pytest.raises(ValueError):
        dynamics.simulate_closing(BASELINE_PARAMS, flyback.DiodeFlyback(), t_max_hold=1e-6)
