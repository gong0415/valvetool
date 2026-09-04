import math

import pytest

from solenoid_model import magnetics, magnetization
from solenoid_model.baseline import BASELINE_PARAMS
from solenoid_model.cases import N2_25BAR_PH_PARAMS

PH = N2_25BAR_PH_PARAMS
I_PEAK = PH.V_bus / PH.R_coil_20C


def test_mag_none_is_bit_identical_to_linear():
    """相容性鐵律：mag=None 必須與現行線性結果逐位元相同。"""
    for params in (BASELINE_PARAMS, PH):
        for gap in (params.g0, params.g0 - params.x_stroke):
            for i in (0.01, 0.1, 0.35):
                linear_B = (params.N_turns * i
                            / magnetics.reluctance(gap, params) / params.A_gap)
                assert magnetics.flux_density(gap, i, params) == linear_B
                assert magnetics.flux_density(gap, i, params, None) == linear_B


def test_solver_satisfies_the_mmf_balance():
    """解回代 H(B)*l_core + B*gap/mu_0 必須等於 N*i。"""
    mag = magnetization.from_params(PH)
    gap = PH.g0
    B = magnetics.solve_flux_density(gap, I_PEAK, PH, mag)
    mmf = mag.H_of_B(B) * PH.l_core + B * gap / magnetics.MU_0
    assert math.isclose(mmf, PH.N_turns * I_PEAK, rel_tol=1e-6)


def test_saturating_B_stays_below_B_sat():
    """線性模型在閉合氣隙報 8.16 T；飽和模型必須低於 B_sat。"""
    mag = magnetization.from_params(PH)
    gap = PH.g0 - PH.x_stroke
    assert magnetics.flux_density(gap, I_PEAK, PH) > 8.0        # 線性：非物理
    assert magnetics.flux_density(gap, I_PEAK, PH, mag) < PH.B_sat


def test_saturating_force_never_exceeds_linear_force():
    """飽和只會減少吸力，不會增加。"""
    mag = magnetization.from_params(PH)
    for k in range(6):
        gap = PH.g0 - PH.x_stroke * k / 5
        F_lin = magnetics.magnetic_force_closing(gap, I_PEAK, PH)
        F_sat = magnetics.magnetic_force_closing(gap, I_PEAK, PH, mag)
        assert F_sat <= F_lin


def test_force_at_rest_gap_matches_verified_value():
    """對照 spec 3A.3 已驗證數值（線性 44.58 N -> 飽和 26.47 N）。"""
    mag = magnetization.from_params(PH)
    F = magnetics.magnetic_force_closing(PH.g0, I_PEAK, PH, mag)
    assert F == pytest.approx(26.47, abs=0.05)


def test_force_at_closed_gap_matches_verified_value():
    """spec 3A.3：閉合氣隙 530.30 N -> 32.70 N。"""
    mag = magnetization.from_params(PH)
    F = magnetics.magnetic_force_closing(PH.g0 - PH.x_stroke, I_PEAK, PH, mag)
    assert F == pytest.approx(32.70, abs=0.05)


def test_saturating_force_is_nearly_flat_across_stroke():
    """spec 3A.3 發現 1：真實吸力沿行程近乎持平（線性模型報 12 倍成長）。"""
    mag = magnetization.from_params(PH)
    F_rest = magnetics.magnetic_force_closing(PH.g0, I_PEAK, PH, mag)
    F_closed = magnetics.magnetic_force_closing(
        PH.g0 - PH.x_stroke, I_PEAK, PH, mag)
    assert F_closed / F_rest < 1.5          # 飽和：+24%
    F_lin_rest = magnetics.magnetic_force_closing(PH.g0, I_PEAK, PH)
    F_lin_closed = magnetics.magnetic_force_closing(
        PH.g0 - PH.x_stroke, I_PEAK, PH)
    assert F_lin_closed / F_lin_rest > 10.0  # 線性：12 倍，對照組


def test_low_current_converges_to_linear_model():
    """電流小到不飽和時，兩模型必須收斂。"""
    mag = magnetization.from_params(PH)
    i = 1e-4
    B_lin = magnetics.flux_density(PH.g0, i, PH)
    B_sat = magnetics.flux_density(PH.g0, i, PH, mag)
    assert math.isclose(B_lin, B_sat, rel_tol=1e-3)


def test_zero_current_gives_zero_flux():
    mag = magnetization.from_params(PH)
    assert magnetics.solve_flux_density(PH.g0, 0.0, PH, mag) == 0.0


def test_saturation_check_honours_mag():
    mag = magnetization.from_params(PH)
    gap = PH.g0 - PH.x_stroke
    assert magnetics.saturation_check(gap, I_PEAK, PH) is True
    assert magnetics.saturation_check(gap, I_PEAK, PH, mag) is False
