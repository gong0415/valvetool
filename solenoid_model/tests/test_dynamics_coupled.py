import numpy as np

from solenoid_model import dynamics, magnetics
from solenoid_model.baseline import BASELINE_PARAMS


def test_armature_stays_pinned_before_pull_in_threshold():
    params = BASELINE_PARAMS
    sol = dynamics.simulate_opening(params, t_max=0.0005)  # 0.5 ms, well before the ~4.26 ms pull-in point
    assert np.all(sol.y[1] == 0.0)
    assert np.all(sol.y[2] == 0.0)


def test_current_matches_rl_charging_while_pinned():
    params = BASELINE_PARAMS
    sol = dynamics.simulate_opening(params, t_max=0.0005)
    L0 = magnetics.inductance(params.g0, params)
    analytical = (params.V_bus / params.R_coil_20C) * (
        1 - np.exp(-sol.t * params.R_coil_20C / L0)
    )
    assert np.allclose(sol.y[0], analytical, rtol=1e-4, atol=1e-6)


def test_armature_reaches_full_stroke_within_10ms():
    params = BASELINE_PARAMS
    sol = dynamics.simulate_opening(params, t_max=0.05)
    assert len(sol.t_events[0]) == 1, "armature never reached full stroke within t_max"
    t_open = sol.t_events[0][0]
    assert 0.0 < t_open < 0.010  # spec target: t_open < 10 ms


def test_current_signature_shows_back_emf_dip_when_armature_moves():
    params = BASELINE_PARAMS
    sol = dynamics.simulate_opening(params, t_max=0.05)
    assert len(sol.t_events[0]) == 1
    i_t = sol.y[0]
    diffs = np.diff(i_t)
    assert np.any(diffs < 0), "current trace never dips; expected back-EMF signature once armature accelerates"
