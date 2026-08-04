import math

import pytest

from solenoid_model import flyback


def test_diode_clamp_voltage_is_negative_v_fwd_while_conducting():
    circuit = flyback.DiodeFlyback(V_fwd=0.7)
    assert flyback.clamp_voltage(circuit, i=0.1) == -0.7


def test_zener_clamp_voltage_is_negative_v_clamp_while_conducting():
    circuit = flyback.ZenerFlyback(V_clamp=30.0)
    assert flyback.clamp_voltage(circuit, i=0.1) == -30.0


def test_clamp_voltage_is_zero_once_current_reaches_zero():
    circuit = flyback.DiodeFlyback(V_fwd=0.7)
    assert flyback.clamp_voltage(circuit, i=0.0) == 0.0
    assert flyback.clamp_voltage(circuit, i=-0.01) == 0.0


def test_clamp_voltage_rejects_rc_circuit():
    circuit = flyback.RCFlyback(R_snub=100.0, C_snub=1e-7)
    with pytest.raises(TypeError):
        flyback.clamp_voltage(circuit, i=0.1)


def test_rc_terminal_voltage_formula():
    circuit = flyback.RCFlyback(R_snub=100.0, C_snub=1e-7)
    # V = -(i*R_snub + V_c) = -(0.05*100.0 + 2.0) = -7.0
    assert math.isclose(
        flyback.rc_terminal_voltage(circuit, i=0.05, V_c=2.0), -7.0
    )


def test_default_diode_and_zener_voltages():
    assert flyback.DiodeFlyback().V_fwd == 0.7
    assert flyback.ZenerFlyback().V_clamp == 30.0
