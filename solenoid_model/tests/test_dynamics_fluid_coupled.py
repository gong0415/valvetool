import math

import pytest

from solenoid_model import dynamics, fluid
from solenoid_model.baseline import BASELINE_PARAMS

N2_VACUUM = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)


def test_coupled_gas_baseline_t_open_reference():
    sol = dynamics.simulate_opening(BASELINE_PARAMS, medium=fluid.N2, cond=N2_VACUUM)
    assert math.isclose(sol.t_events[0][0], 5.0312899433304515e-3, rel_tol=1e-6)


def test_coupled_t_open_not_faster_than_dry():
    dry = dynamics.simulate_opening(BASELINE_PARAMS).t_events[0][0]
    coupled = dynamics.simulate_opening(BASELINE_PARAMS, medium=fluid.N2, cond=N2_VACUUM).t_events[0][0]
    assert coupled >= dry


def test_medium_without_conditions_raises():
    with pytest.raises(ValueError):
        dynamics.simulate_opening(BASELINE_PARAMS, medium=fluid.N2)


def test_coupled_liquid_medium_opens_and_is_not_faster():
    cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
    dry = dynamics.simulate_opening(BASELINE_PARAMS).t_events[0][0]
    sol = dynamics.simulate_opening(BASELINE_PARAMS, medium=fluid.WATER_20C, cond=cond)
    assert len(sol.t_events[0]) == 1
    assert sol.t_events[0][0] >= dry
