import math

import pytest

from solenoid_model import drive, magnetization
from solenoid_model.cases import N2_25BAR_PH_PARAMS

PH = N2_25BAR_PH_PARAMS


def test_mag_none_matches_existing_hold_current():
    """相容性：mag=None 必須與現行解析反解相同。"""
    mag_free = drive.hold_current(PH, margin=2.0)
    assert mag_free == pytest.approx(0.01143954024087197, rel=1e-12)


def test_linear_hold_current_overstates_the_margin():
    """spec 3A.3 發現 2：線性算出的 11.44 mA 在飽和下只有 1.60x。"""
    mag = magnetization.from_params(PH)
    i_linear = drive.hold_current(PH, margin=2.0)
    true_margin = drive.force_margin_at_current(PH, i_linear, mag)
    assert true_margin == pytest.approx(1.60, abs=0.02)


def test_saturating_hold_current_restores_true_margin():
    """真 2.0x 需 13.07 mA。"""
    mag = magnetization.from_params(PH)
    i = drive.hold_current(PH, margin=2.0, mag=mag)
    assert i == pytest.approx(0.01307, abs=5e-5)
    assert drive.force_margin_at_current(PH, i, mag) == pytest.approx(
        2.0, rel=1e-6)


def test_saturating_hold_current_exceeds_linear():
    """飽和使吸力變小，所以需要更多電流。"""
    mag = magnetization.from_params(PH)
    assert drive.hold_current(PH, margin=2.0, mag=mag) > drive.hold_current(
        PH, margin=2.0)


def test_hold_power_penalty_is_reported_correctly():
    """0.037 W -> 0.048 W（+31%）。"""
    mag = magnetization.from_params(PH)
    P_lin = drive.hold_current(PH, 2.0) ** 2 * PH.R_coil_20C
    P_sat = drive.hold_current(PH, 2.0, mag=mag) ** 2 * PH.R_coil_20C
    assert P_lin == pytest.approx(0.0368, abs=5e-4)
    assert P_sat == pytest.approx(0.0481, abs=5e-4)


def test_force_margin_round_trips_under_saturation():
    mag = magnetization.from_params(PH)
    for margin in (1.5, 2.0, 3.0):
        i = drive.hold_current(PH, margin, mag=mag)
        assert drive.force_margin_at_current(PH, i, mag) == pytest.approx(
            margin, rel=1e-6)


def test_hold_current_rejects_nonpositive_margin():
    mag = magnetization.from_params(PH)
    with pytest.raises(ValueError, match="margin"):
        drive.hold_current(PH, margin=0.0, mag=mag)
