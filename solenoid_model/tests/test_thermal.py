import math
from dataclasses import replace

import pytest

from solenoid_model import thermal
from solenoid_model.baseline import BASELINE_PARAMS


# --- R_coil:規格域線性 α ---

def test_r_coil_at_reference_temp_returns_r20_exactly():
    # 1 + α·0.0 乘法在 IEEE 754 下精確 → GUI 預設 293.15 K 不改變任何結果的根據
    assert thermal.R_coil(293.15, BASELINE_PARAMS) == BASELINE_PARAMS.R_coil_20C


def test_r_coil_hot_70c_matches_linear_alpha():
    # 80 × (1 + 0.00393 × 50) = 95.72
    assert math.isclose(thermal.R_coil(343.15, BASELINE_PARAMS), 95.72, rel_tol=1e-12)


# --- R_coil:低溫查表段 ---

def test_r_coil_ln2_77k_table_value():
    # 80 × 0.13 = 10.4(工程表值,RRR 影響 ±10% 量級,文件標註)
    assert math.isclose(thermal.R_coil(77.0, BASELINE_PARAMS), 10.4, rel_tol=1e-12)


def test_r_coil_mid_table_interpolates():
    # 100 K 表點:80 × 0.21 = 16.8
    assert math.isclose(thermal.R_coil(100.0, BASELINE_PARAMS), 16.8, rel_tol=1e-12)


def test_r_coil_seam_continuous_at_minus_40c():
    # 表錨點 = 線性模型在 233.15 K 的值 → 兩支連續
    r_linear = thermal.R_coil(233.15, BASELINE_PARAMS)
    r_table_side = thermal.R_coil(233.15 - 1e-6, BASELINE_PARAMS)
    assert math.isclose(r_linear, 61.136, rel_tol=1e-9)
    assert math.isclose(r_table_side, r_linear, rel_tol=1e-6)


def test_r_coil_below_table_floor_raises():
    with pytest.raises(ValueError):
        thermal.R_coil(76.0, BASELINE_PARAMS)


# --- equilibrium_temp ---

def test_equilibrium_temp_zero_current_returns_ambient():
    assert thermal.equilibrium_temp(0.0, 293.15, BASELINE_PARAMS) == 293.15


def test_equilibrium_temp_baseline_anchor():
    # planning 階段原型定案:I=0.35 A、T_amb=293.15 K → 350.5829 K(ΔT≈57.4 K)
    T_eq = thermal.equilibrium_temp(0.35, 293.15, BASELINE_PARAMS)
    assert math.isclose(T_eq, 350.5829, abs_tol=1e-3)


def test_equilibrium_temp_satisfies_energy_balance():
    p = BASELINE_PARAMS
    T_amb, I = 293.15, 0.35
    T_eq = thermal.equilibrium_temp(I, T_amb, p)
    losses = (p.G_th_cond * (T_eq - T_amb)
              + p.emissivity * thermal.STEFAN_BOLTZMANN * p.A_rad
              * (T_eq**4 - T_amb**4))
    assert math.isclose(I**2 * thermal.R_coil(T_eq, p), losses, rel_tol=1e-6)


def test_equilibrium_temp_thermal_runaway_raises():
    # 定電流 + 幾乎無散熱路徑 → T_max 內無平衡
    p = replace(BASELINE_PARAMS, G_th_cond=0.0, A_rad=1e-9)
    with pytest.raises(ValueError, match="runaway"):
        thermal.equilibrium_temp(0.35, 293.15, p)


# --- curie_margin ---

def test_curie_margin_room_temp():
    assert math.isclose(thermal.curie_margin(293.15), 750.0, rel_tol=1e-12)


# --- CoilMaterial 表驗證 + R_coil NaN 守衛(P2) ---

def test_coil_material_rejects_non_ascending_table():
    with pytest.raises(ValueError, match="ascending"):
        thermal.CoilMaterial(alpha=0.004, T_ref=293.15,
                             r_table=((100.0, 0.21), (77.0, 0.13), (233.15, 0.76)))


def test_coil_material_rejects_single_row_table():
    with pytest.raises(ValueError, match="at least"):
        thermal.CoilMaterial(alpha=0.004, T_ref=293.15, r_table=((77.0, 0.13),))


def test_coil_material_rejects_non_finite_or_non_positive_values():
    # NaN 表值、inf α、非正 R/R₀ 各自都要擋下
    with pytest.raises(ValueError, match="finite"):
        thermal.CoilMaterial(alpha=0.004, T_ref=293.15,
                             r_table=((77.0, float("nan")), (233.15, 0.76)))
    with pytest.raises(ValueError, match="finite"):
        thermal.CoilMaterial(alpha=float("inf"), T_ref=293.15,
                             r_table=((77.0, 0.13), (233.15, 0.76)))
    with pytest.raises(ValueError, match="> 0"):
        thermal.CoilMaterial(alpha=0.004, T_ref=293.15,
                             r_table=((77.0, -0.1), (233.15, 0.76)))


def test_r_coil_nan_temperature_raises_value_error():
    # 現行為:NaN 穿過所有區間比較,落到 AssertionError——錯誤類型誤導
    with pytest.raises(ValueError, match="NaN"):
        thermal.R_coil(float("nan"), BASELINE_PARAMS)
