import math

from solenoid_model import winding
from solenoid_model.baseline import BASELINE_PARAMS


def test_baseline_params_have_expected_values():
    p = BASELINE_PARAMS
    assert p.V_bus == 28.0
    assert p.R_coil_20C == 80.0
    assert p.N_turns == 2000.0
    assert p.A_gap == 3.0e-5
    assert p.g0 == 3.5e-4
    assert p.x_stroke == 3.0e-4
    assert p.m_arm == 0.001
    assert p.k_spring == 20000.0
    assert p.F_preload == 10.0
    assert p.delta_P == 2.4e6
    assert p.A_seat == 5.0e-6
    assert p.damping_coeff == 0.1


def test_valve_params_is_a_plain_dataclass_so_replace_works():
    from dataclasses import replace
    p2 = replace(BASELINE_PARAMS, N_turns=3000.0)
    assert p2.N_turns == 3000.0
    assert BASELINE_PARAMS.N_turns == 2000.0  # original untouched


def test_baseline_winding_defaults_are_calibrated_to_r_coil():
    # the winding defaults must reproduce the independent R_coil_20C at N=2000
    p = BASELINE_PARAMS
    R = winding.coil_resistance(p.N_turns, p.A_winding, p.k_fill, p.l_turn_mean)
    assert math.isclose(R, p.R_coil_20C, rel_tol=1e-9)
    # exact default values
    assert p.k_fill == 0.5
    assert math.isclose(p.A_winding, 5.2190904529497215e-05, rel_tol=1e-12)
    assert math.isclose(p.l_turn_mean, 0.031066014600891194, rel_tol=1e-12)


def test_baseline_sizing_defaults():
    p = BASELINE_PARAMS
    assert p.k_pack == 3.0
    assert math.isclose(p.J_max, 20e6, rel_tol=1e-12)
