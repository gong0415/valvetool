"""Drive-point solver tests: the reverse-solved currents must round-trip.

drive.py answers "what current does this gap need?" by inverting the same
linear-magnetics force law dynamics integrates. These tests pin the inversion
against the forward model (feed the solved current back in, get the force
asked for) rather than against remembered constants, so a change to the force
law surfaces here instead of silently shifting a design point.
"""
import math

import pytest

from solenoid_model import drive, limits
from solenoid_model.baseline import BASELINE_PARAMS
from solenoid_model.cases import N2_25BAR_PH, N2_25BAR_PH_PARAMS
from solenoid_model.dynamics import spring_force
from solenoid_model.magnetics import magnetic_force_closing


def test_current_for_force_margin_round_trips_through_force_law():
    """The solved current must produce exactly margin * F_resist of force."""
    p = N2_25BAR_PH_PARAMS
    gap = p.g0 - p.x_stroke
    F_resist = 3.5
    for margin in (1.0, 2.0, 5.0):
        i = drive.current_for_force_margin(gap, margin, p, F_resist)
        assert magnetic_force_closing(gap, i, p) == pytest.approx(margin * F_resist)


def test_current_scales_as_sqrt_margin():
    """F ~ i^2, so doubling the force margin costs sqrt(2) in current."""
    p = N2_25BAR_PH_PARAMS
    gap = p.g0 - p.x_stroke
    i1 = drive.current_for_force_margin(gap, 1.0, p, 3.5)
    i2 = drive.current_for_force_margin(gap, 2.0, p, 3.5)
    assert i2 / i1 == pytest.approx(math.sqrt(2.0))


def test_hold_current_balances_spring_plus_pressure_at_open_gap():
    """At margin=1.0 the magnetic force must exactly equal the resisting load
    the armature sees while held open (spring at full stroke + static pressure)."""
    p = N2_25BAR_PH_PARAMS
    i = drive.hold_current(p, margin=1.0)
    gap_open = p.g0 - p.x_stroke
    F_resist = spring_force(p.x_stroke, p) + p.delta_P * p.A_seat
    assert magnetic_force_closing(gap_open, i, p) == pytest.approx(F_resist)


def test_pull_in_current_matches_limits_motion_threshold():
    """pull_in_current(margin=1.0) is the same quantity limits.py already
    exposed as motion_threshold_current; the refactor must not move it."""
    for p in (BASELINE_PARAMS, N2_25BAR_PH_PARAMS):
        assert drive.pull_in_current(p, margin=1.0) == pytest.approx(
            limits.motion_threshold_current(p))


def test_pull_in_costs_more_current_than_hold():
    """The peak-and-hold premise: the g0 gap is expensive, the closed gap is
    cheap. This headroom ratio is what constant-voltage drive wastes as heat."""
    p = N2_25BAR_PH_PARAMS
    ratio = drive.pull_in_current(p, margin=1.0) / drive.hold_current(p, margin=1.0)
    assert ratio == pytest.approx(3.14, abs=0.01)


def test_ph_case_hold_current_is_derived_not_hardcoded():
    """The case's recorded i_hold must be reproducible from the force balance
    at its stated 2.0x margin."""
    assert N2_25BAR_PH.i_hold == pytest.approx(
        drive.hold_current(N2_25BAR_PH_PARAMS, margin=2.0))
    assert N2_25BAR_PH.i_hold == pytest.approx(0.011440, abs=1e-6)


def test_hold_duty_round_trips_to_current():
    """D = i*R/V_bus is just KVL; multiplying back by V_bus/R returns i."""
    p = N2_25BAR_PH_PARAMS
    i = drive.hold_current(p, margin=2.0)
    D = drive.hold_duty(i, p)
    assert D * p.V_bus / p.R_coil_20C == pytest.approx(i)
    assert 0.0 < D < 1.0


def test_hold_duty_honours_explicit_resistance():
    """A hot-coil R must raise the duty needed for the same current."""
    p = N2_25BAR_PH_PARAMS
    i = drive.hold_current(p, margin=2.0)
    D_cold = drive.hold_duty(i, p)
    D_hot = drive.hold_duty(i, p, R=p.R_coil_20C * 1.2)
    assert D_hot == pytest.approx(D_cold * 1.2)


def test_rejects_nonpositive_margin():
    p = N2_25BAR_PH_PARAMS
    with pytest.raises(ValueError):
        drive.hold_current(p, margin=0.0)
    with pytest.raises(ValueError):
        drive.pull_in_current(p, margin=-1.0)


def test_rejects_duty_above_unity():
    """A demanded current the bus cannot source is a design error, not a
    silently-clipped duty."""
    p = N2_25BAR_PH_PARAMS
    i_unreachable = 10.0 * p.V_bus / p.R_coil_20C
    with pytest.raises(ValueError):
        drive.hold_duty(i_unreachable, p)
