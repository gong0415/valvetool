"""P2 thermal coupling: T_coil=None must be bit-identical to the dry run
(characterization tests are the umbrella); a given T_coil swaps in
thermal.R_coil(T_coil, params) as a constant for the whole simulation."""
import math

import pytest

from solenoid_model import dynamics, fluid
from solenoid_model.baseline import BASELINE_PARAMS


def test_t_coil_at_reference_temp_is_bit_identical_to_dry_run():
    sol_dry = dynamics.simulate_opening(BASELINE_PARAMS)
    sol_ref = dynamics.simulate_opening(BASELINE_PARAMS, T_coil=293.15)
    # R(293.15) == R_coil_20C 精確 → 整條軌跡逐位相同
    assert sol_ref.t_events[0][0] == sol_dry.t_events[0][0]
    assert sol_ref.y[0][-1] == sol_dry.y[0][-1]


def test_hot_coil_70c_slows_opening_to_anchor():
    # planning 原型定案:R=95.72 Ω → t_open = 5.683742 ms
    sol = dynamics.simulate_opening(BASELINE_PARAMS, T_coil=343.15)
    assert math.isclose(sol.t_events[0][0], 5.683742e-3, rel_tol=1e-6)


def test_ln2_coil_77k_speeds_opening_to_anchor():
    # planning 原型定案:R=10.4 Ω → t_open = 3.572129 ms
    sol = dynamics.simulate_opening(BASELINE_PARAMS, T_coil=77.0)
    assert math.isclose(sol.t_events[0][0], 3.572129e-3, rel_tol=1e-6)


def test_thermal_and_fluid_coupling_are_orthogonal():
    # 熱只動電氣路徑的 R,流體只動機械路徑的力 → 可同時開
    cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
    sol = dynamics.simulate_opening(
        BASELINE_PARAMS, medium=fluid.N2, cond=cond, T_coil=343.15)
    assert sol.t_events[0].size == 1
    # 熱態耦合比常溫耦合更慢——跟常溫耦合案例（同一 medium/cond,無 T_coil）比較,
    # 不用硬編常數,才不會在 baseline 流道幾何變動時悄悄過時
    room_temp = dynamics.simulate_opening(BASELINE_PARAMS, medium=fluid.N2, cond=cond)
    assert sol.t_events[0][0] > room_temp.t_events[0][0]


def test_out_of_domain_temperature_propagates_value_error():
    with pytest.raises(ValueError):
        dynamics.simulate_opening(BASELINE_PARAMS, T_coil=76.0)
