import numpy as np
from scipy.integrate import solve_ivp

from solenoid_model import dynamics, magnetics
from solenoid_model.baseline import BASELINE_PARAMS


def test_electrical_di_dt_matches_rl_step_response_at_fixed_gap():
    params = BASELINE_PARAMS
    gap = params.g0  # fixed gap: no motion, so back-EMF term is zero (v=0)

    def rhs(t, y):
        return [dynamics.electrical_di_dt(y[0], gap, 0.0, params)]

    t_eval = np.linspace(0.0, 0.02, 200)
    sol = solve_ivp(rhs, (0.0, 0.02), [0.0], t_eval=t_eval, rtol=1e-10, atol=1e-12)

    L = magnetics.inductance(gap, params)
    analytical = (params.V_bus / params.R_coil_20C) * (1 - np.exp(-t_eval * params.R_coil_20C / L))
    assert np.allclose(sol.y[0], analytical, rtol=1e-4, atol=1e-6)


def test_electrical_di_dt_v_override_replaces_v_bus():
    params = BASELINE_PARAMS
    gap = params.g0 - params.x_stroke  # fixed gap: no motion, back-EMF term zero
    L = magnetics.inductance(gap, params)
    R = params.R_coil_20C
    I_ss = params.V_bus / R
    V = -30.0  # arbitrary override, unrelated to params.V_bus

    def rhs(t, y):
        return [dynamics.electrical_di_dt(y[0], gap, 0.0, params, V=V)]

    t_eval = np.linspace(0.0, 0.02, 200)
    sol = solve_ivp(rhs, (0.0, 0.02), [I_ss], t_eval=t_eval, rtol=1e-10, atol=1e-12)

    i_ss_new = V / R
    analytical = i_ss_new + (I_ss - i_ss_new) * np.exp(-t_eval * R / L)
    assert np.allclose(sol.y[0], analytical, rtol=1e-4, atol=1e-6)
