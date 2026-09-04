import math
from dataclasses import replace

import pytest

from solenoid_model import layout, sizing
from solenoid_model.cases import N2_25BAR_PH_PARAMS

PH = N2_25BAR_PH_PARAMS


def test_layout_does_not_import_physics_layers():
    """分層不變式：layout 是純代數層。"""
    import solenoid_model.layout as m
    src = open(m.__file__, encoding="utf-8").read()
    for banned in ("import dynamics", "import magnetics", "import fluid",
                   "from solenoid_model import dynamics",
                   "from solenoid_model import magnetics",
                   "from solenoid_model import fluid"):
        assert banned not in src


def test_yoke_area_at_B_sat_equals_gap_area():
    """磁通連續 + B_sat 鉗位：同材質下 A_yoke >= A_gap。"""
    assert layout.yoke_area(PH) == pytest.approx(PH.A_gap, rel=1e-12)


def test_yoke_area_scales_with_real_pole_flux():
    """真實 B 較低時所需截面等比例縮小（spec 3: 1.824 T -> 16.97 mm2）。"""
    A = layout.yoke_area(PH, B_pole=1.824)
    assert A == pytest.approx(16.97e-6, rel=1e-3)
    assert A < layout.yoke_area(PH)


def test_yoke_area_rejects_B_pole_above_B_sat():
    with pytest.raises(ValueError, match="B_sat"):
        layout.yoke_area(PH, B_pole=3.0)


def test_core_radius_and_pole_diameter_match_A_gap():
    assert layout.pole_diameter(PH) == pytest.approx(5.05e-3, abs=5e-6)
    assert layout.core_radius(PH) == pytest.approx(
        layout.pole_diameter(PH) / 2, rel=1e-12)


def test_shell_wall_magnetic_gives_annulus_of_yoke_area():
    R_i = 9.0e-3
    t = layout.shell_wall_magnetic(PH, R_i)
    assert t == pytest.approx(0.35e-3, abs=1e-5)
    A_annulus = math.pi * ((R_i + t) ** 2 - R_i ** 2)
    assert A_annulus == pytest.approx(layout.yoke_area(PH), rel=1e-9)


def test_hoop_wall_thickness_matches_thin_wall_formula():
    t = layout.hoop_wall_thickness(P=25e5, R_inner=9.0e-3,
                                   S_yield=205e6, SF=3.0)
    assert t == pytest.approx(329.3e-6, abs=1e-6)


def test_hoop_wall_scales_with_safety_factor():
    a = layout.hoop_wall_thickness(25e5, 9.0e-3, 205e6, 3.0)
    b = layout.hoop_wall_thickness(25e5, 9.0e-3, 205e6, 6.0)
    assert b == pytest.approx(2.0 * a, rel=1e-12)


def test_shell_wall_adopts_the_largest_of_three():
    w = layout.shell_wall_thickness(PH, R_inner=9.0e-3, P=25e5,
                                    S_yield=205e6, SF=3.0,
                                    t_manufacturing=1.0e-3)
    assert w.adopted == pytest.approx(1.0e-3, rel=1e-12)
    assert w.adopted == max(w.magnetic, w.structural, w.manufacturing)
    assert "manufacturing" in w.reason


def test_shell_wall_reason_names_the_binding_constraint():
    w = layout.shell_wall_thickness(PH, R_inner=9.0e-3, P=250e5,
                                    S_yield=205e6, SF=3.0,
                                    t_manufacturing=0.1e-3)
    assert w.adopted == pytest.approx(w.structural, rel=1e-12)
    assert "structural" in w.reason


def test_end_plate_thickness_carries_yoke_flux_at_core_radius():
    t = layout.end_plate_thickness(PH)
    assert t == pytest.approx(1.26e-3, abs=1e-5)
    A_radial = 2 * math.pi * layout.core_radius(PH) * t
    assert A_radial == pytest.approx(layout.yoke_area(PH), rel=1e-9)


def test_spring_geometry_round_trips_to_k_spring():
    """由 (d, D, n) 正算回 k 必須等於 params.k_spring。"""
    g = layout.spring_geometry(PH, d_wire=0.35e-3, D_coil=1.8e-3, G=79e9)
    k_back = 79e9 * (0.35e-3) ** 4 / (8 * (1.8e-3) ** 3 * g.n_active)
    assert k_back == pytest.approx(PH.k_spring, rel=1e-9)


def test_spring_geometry_matches_verified_design_point():
    """spec 6 已驗證：n=6.35, C=5.14, tau=389 MPa, L_free=3.72 mm。"""
    g = layout.spring_geometry(PH, d_wire=0.35e-3, D_coil=1.8e-3, G=79e9)
    assert g.n_active == pytest.approx(6.35, abs=0.02)
    assert g.index == pytest.approx(5.14, abs=0.02)
    assert g.tau_max == pytest.approx(389e6, rel=0.02)
    assert g.L_free == pytest.approx(3.72e-3, abs=5e-5)


def test_spring_wahl_factor_approaches_one_for_large_index():
    """C -> infinity 時 Wahl 修正 -> 1。"""
    assert layout.wahl_factor(1e6) == pytest.approx(1.0, abs=1e-5)
    assert layout.wahl_factor(5.0) > 1.2


def test_spring_solid_length_is_below_free_length():
    g = layout.spring_geometry(PH, d_wire=0.35e-3, D_coil=1.8e-3, G=79e9)
    assert g.L_solid < g.L_free


def test_armature_thickness_derives_from_declared_mass():
    """銜鐵厚度由 m_arm 反解，非 EMPIRICAL。"""
    t = layout.armature_thickness(PH, rho=7870.0,
                                  D=layout.pole_diameter(PH))
    assert t == pytest.approx(5.08e-3, abs=2e-5)


def test_armature_mass_consistency_flags_mismatch():
    """4.0 mm 只給 0.631 g，與宣告的 0.8 g 不符。"""
    m_dec, m_geo, rel = layout.armature_mass_consistency(
        PH, t_arm=4.0e-3, rho=7870.0, D=layout.pole_diameter(PH))
    assert m_dec == pytest.approx(0.8e-3, rel=1e-12)
    assert m_geo == pytest.approx(0.631e-3, abs=5e-6)
    assert rel > 0.2


def test_armature_mass_consistency_passes_at_derived_thickness():
    t = layout.armature_thickness(PH, rho=7870.0, D=layout.pole_diameter(PH))
    _, _, rel = layout.armature_mass_consistency(
        PH, t_arm=t, rho=7870.0, D=layout.pole_diameter(PH))
    assert rel < 1e-9


def test_axial_stack_total_equals_sum_of_segments():
    stack = layout.axial_stack(PH)
    assert stack["total"] == pytest.approx(
        sum(v for k, v in stack["segments"]), rel=1e-12)


def test_axial_stack_matches_verified_total():
    """spec 7：合計 18.68 mm。"""
    assert layout.axial_stack(PH)["total"] == pytest.approx(18.679e-3, abs=2e-5)


def test_axial_stack_includes_stroke_and_working_gap_from_params():
    stack = dict(layout.axial_stack(PH)["segments"])
    assert stack["行程 x_stroke"] == pytest.approx(PH.x_stroke, rel=1e-12)
    assert stack["工作氣隙 g0"] == pytest.approx(PH.g0, rel=1e-12)


def test_envelope_outer_diameter_uses_adopted_wall():
    env = layout.envelope(PH, coil_OD=18.0e-3, t_manufacturing=1.0e-3)
    assert env["OD"] == pytest.approx(20.0e-3, abs=1e-5)


def test_lengths_scale_linearly_under_self_similar_scaling():
    """自相似縮放下長度 ~ s、截面 ~ s^2。"""
    s = 2.0
    scaled = sizing.scale_params(PH, s)
    assert layout.core_radius(scaled) == pytest.approx(
        s * layout.core_radius(PH), rel=1e-9)
    assert layout.yoke_area(scaled) == pytest.approx(
        s ** 2 * layout.yoke_area(PH), rel=1e-9)


def test_every_empirical_entry_is_fully_declared():
    assert layout.LAYOUT_EMPIRICAL
    for name, entry in layout.LAYOUT_EMPIRICAL.items():
        for field in ("value", "unit", "reason", "source"):
            assert field in entry, f"{name} missing {field}"
        assert entry["source"] == "EMPIRICAL"
        assert entry["reason"].strip()
