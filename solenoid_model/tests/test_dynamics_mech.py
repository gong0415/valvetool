import math
from dataclasses import replace

import numpy as np
from scipy.integrate import solve_ivp

from solenoid_model import dynamics
from solenoid_model.baseline import BASELINE_PARAMS


def test_spring_force_is_preload_at_rest_and_increases_linearly():
    params = BASELINE_PARAMS
    assert dynamics.spring_force(0.0, params) == params.F_preload
    assert math.isclose(
        dynamics.spring_force(1e-4, params),
        params.F_preload + params.k_spring * 1e-4,
    )


def test_mechanical_dv_dt_matches_simple_harmonic_motion():
    # Isolate pure spring-mass: no preload, no pressure load, no damping, no magnetic force.
    params = replace(BASELINE_PARAMS, F_preload=0.0, delta_P=0.0, damping_coeff=0.0)
    x0 = 5e-5
    omega = math.sqrt(params.k_spring / params.m_arm)

    def rhs(t, y):
        x, v = y
        return [v, dynamics.mechanical_dv_dt(x, v, 0.0, params)]

    t_eval = np.linspace(0.0, 0.002, 200)
    sol = solve_ivp(rhs, (0.0, 0.002), [x0, 0.0], t_eval=t_eval, rtol=1e-10, atol=1e-14)

    analytical = x0 * np.cos(omega * t_eval)
    assert np.allclose(sol.y[0], analytical, rtol=1e-4, atol=1e-9)
