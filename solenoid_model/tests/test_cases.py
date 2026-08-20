"""Design-case regression tests: the recorded numbers must stay re-derivable.

These lock the *sizing claims* in cases.py (flow target, pull-in margin, hold
power), so an edit to the physics or to a case constant that breaks the design
intent fails here rather than silently shipping a valve that cannot open.
"""
import math

import pytest

from solenoid_model import fluid
from solenoid_model.cases import CASES, N2_25BAR
from solenoid_model.dynamics import simulate_opening
from solenoid_model.magnetics import flux_density, magnetic_force_closing
from solenoid_model.thermal import equilibrium_temp
from solenoid_model.winding import coil_resistance


def test_cases_registry_keyed_by_name():
    for name, case in CASES.items():
        assert case.name == name


def test_n2_case_meets_mass_flow_target():
    """Bore was reverse-solved from 1.28 g/s; full-open flow must still hit it."""
    p = N2_25BAR.params
    mdot = fluid.mdot_gas(p.x_stroke, p, N2_25BAR.medium, N2_25BAR.cond)
    assert mdot == pytest.approx(1.28e-3, rel=0.05)


def test_n2_case_is_choked():
    r = N2_25BAR.cond.P_down / N2_25BAR.cond.P_up
    assert r < fluid.critical_pressure_ratio(N2_25BAR.medium)


def test_n2_stroke_is_bore_limited():
    """x_stroke = D/4 is where curtain area meets bore area — more travel is wasted."""
    p = N2_25BAR.params
    assert p.x_stroke == pytest.approx(p.D_seat_bore / 4, rel=1e-6)
    a_full = fluid.effective_area(p.x_stroke, p)
    assert a_full == pytest.approx(math.pi * p.D_seat_bore ** 2 / 4, rel=1e-6)


def test_n2_resistance_consistent_with_winding_window():
    """The app warns if these diverge >5%; keep them exactly co-derived."""
    p = N2_25BAR.params
    assert coil_resistance(p.N_turns, p.A_winding, p.k_fill,
                           p.l_turn_mean) == pytest.approx(p.R_coil_20C, rel=1e-9)


def test_n2_pulls_in_with_margin_when_hot():
    """Worst case is the hot coil: higher R, less current, less force."""
    p = N2_25BAR.params
    T_eq = equilibrium_temp(p.V_bus / p.R_coil_20C, N2_25BAR.T_coil, p)
    from solenoid_model.thermal import R_coil
    i_hot = p.V_bus / R_coil(T_eq, p)
    resisting = p.F_preload + N2_25BAR.cond.P_up * p.A_seat
    assert magnetic_force_closing(p.g0, i_hot, p) > 2.0 * resisting


def test_n2_hold_power_under_one_watt():
    p = N2_25BAR.params
    assert p.V_bus ** 2 / p.R_coil_20C < 1.0


def test_n2_not_saturated_at_closed_gap():
    """Force is scarcest at the closed gap; that is where saturation would hurt."""
    p = N2_25BAR.params
    assert flux_density(p.g0, p.V_bus / p.R_coil_20C, p) < p.B_sat


def test_n2_opens_within_10ms():
    p = N2_25BAR.params
    res = simulate_opening(p, medium=N2_25BAR.medium, cond=N2_25BAR.cond,
                           T_coil=N2_25BAR.T_coil)
    x = res.y[1]
    reached = [t for t, xv in zip(res.t, x) if xv >= p.x_stroke * 0.999]
    assert reached, "case must reach full stroke"
    assert reached[0] < 10e-3
