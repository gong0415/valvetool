import math

from solenoid_model import fluid
from solenoid_model.baseline import BASELINE_PARAMS


def test_valve_params_has_seat_bore_and_cd_with_baseline_values():
    assert BASELINE_PARAMS.D_seat_bore == 0.14e-3
    assert BASELINE_PARAMS.C_d == 0.8


def test_valve_params_new_fields_have_defaults():
    from solenoid_model.params import ValveParams
    import dataclasses
    fields = {f.name: f for f in dataclasses.fields(ValveParams)}
    assert fields["D_seat_bore"].default == 0.52e-3
    assert fields["C_d"].default == 0.8


def test_gas_presets_match_documented_constants():
    assert math.isclose(fluid.N2.gamma, 1.4)
    assert math.isclose(fluid.N2.R_specific, 296.8)
    assert math.isclose(fluid.XE.gamma, 1.667)
    assert math.isclose(fluid.XE.R_specific, 63.33)
    assert math.isclose(fluid.HE.R_specific, 2077.1)


def test_liquid_presets_match_documented_constants():
    assert math.isclose(fluid.WATER_20C.rho, 998.2)
    assert math.isclose(fluid.WATER_20C.c_sound, 1482.0)
    assert math.isclose(fluid.LN2_77K.rho, 806.1)
    assert math.isclose(fluid.LN2_77K.P_vap, 101325.0)


def test_effective_area_zero_at_zero_lift():
    assert fluid.effective_area(0.0, BASELINE_PARAMS) == 0.0


def test_effective_area_curtain_bore_crossover_at_quarter_diameter():
    params = BASELINE_PARAMS
    x_cross = params.D_seat_bore / 4
    curtain = math.pi * params.D_seat_bore * x_cross
    bore = math.pi * params.D_seat_bore ** 2 / 4
    assert math.isclose(curtain, bore, rel_tol=1e-12)
    assert math.isclose(fluid.effective_area(x_cross, params), bore, rel_tol=1e-12)
    # beyond crossover: saturates at bore area
    assert math.isclose(fluid.effective_area(params.x_stroke, params), bore, rel_tol=1e-12)
    # below crossover: curtain-limited (relative to D_seat_bore so this stays
    # below the crossover regardless of which bore size baseline uses)
    x1 = params.D_seat_bore / 8
    assert math.isclose(fluid.effective_area(x1, params), math.pi * params.D_seat_bore * x1, rel_tol=1e-12)


def test_critical_pressure_ratio_n2():
    assert math.isclose(fluid.critical_pressure_ratio(fluid.N2), 0.5282817877171742, rel_tol=1e-12)


def test_mdot_gas_choked_matches_hand_computed_reference():
    cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
    mdot = fluid.mdot_gas(BASELINE_PARAMS.x_stroke, BASELINE_PARAMS, fluid.N2, cond)
    assert math.isclose(mdot, 6.862805689539614e-5, rel_tol=1e-9)


def test_mdot_gas_continuous_at_critical_ratio():
    params = BASELINE_PARAMS
    r_crit = fluid.critical_pressure_ratio(fluid.N2)
    P_up = 2.4e6
    just_below = fluid.mdot_gas(params.x_stroke, params, fluid.N2,
                                fluid.FlowConditions(P_up, P_up * (r_crit - 1e-9), 293.0))
    just_above = fluid.mdot_gas(params.x_stroke, params, fluid.N2,
                                fluid.FlowConditions(P_up, P_up * (r_crit + 1e-9), 293.0))
    assert math.isclose(just_below, just_above, rel_tol=1e-6)


def test_mdot_gas_subsonic_reference_and_limits():
    params = BASELINE_PARAMS
    P_up = 2.4e6
    mdot_09 = fluid.mdot_gas(params.x_stroke, params, fluid.N2,
                             fluid.FlowConditions(P_up, 0.9 * P_up, 293.0))
    assert math.isclose(mdot_09, 4.2353689603347495e-5, rel_tol=1e-9)
    # r -> 1: flow vanishes
    mdot_tiny = fluid.mdot_gas(params.x_stroke, params, fluid.N2,
                               fluid.FlowConditions(P_up, P_up * 0.999999, 293.0))
    assert mdot_tiny < 1e-5
    # monotonic in P_up (choked regime: mdot proportional to P_up)
    mdot_2x = fluid.mdot_gas(params.x_stroke, params, fluid.N2,
                             fluid.FlowConditions(2 * P_up, 0.0, 293.0))
    choked_ref = fluid.mdot_gas(params.x_stroke, params, fluid.N2,
                                fluid.FlowConditions(P_up, 0.0, 293.0))
    assert math.isclose(mdot_2x, 2 * choked_ref, rel_tol=1e-12)


def test_mdot_liquid_matches_hand_computed_reference():
    cond = fluid.FlowConditions(P_up=2.0e5, P_down=1.0e5, T0=293.0)
    mdot = fluid.mdot_liquid(1.0e-4, BASELINE_PARAMS, fluid.WATER_20C, cond)
    assert math.isclose(mdot, 1.74004195666987e-4, rel_tol=1e-9)


def test_water_hammer_joukowsky_reference():
    dp = fluid.water_hammer_dp(fluid.LN2_77K, delta_v=2.0)
    assert math.isclose(dp, 806.1 * 850.0 * 2.0, rel_tol=1e-12)


def test_flashing_risk_ln2_at_boiling_true_subcooled_false():
    at_boiling = fluid.FlowConditions(P_up=5.0e5, P_down=101325.0, T0=77.0)
    assert fluid.flashing_risk(fluid.LN2_77K, at_boiling) is True
    subcooled = fluid.FlowConditions(P_up=5.0e5, P_down=5.0e5, T0=77.0)
    assert fluid.flashing_risk(fluid.LN2_77K, subcooled) is False


def test_full_open_cv_meets_baseline_target():
    # Baseline flow target (spec §6.4 baseline case flow item): a
    # small-thrust-class attitude-control thruster valve (Cv~0.0007,
    # commercial ~5 lbf/22N redundant-seat range).
    cv = fluid.cv_from_geometry(BASELINE_PARAMS.x_stroke, BASELINE_PARAMS)
    assert math.isclose(cv, 0.0007, rel_tol=0.04)   # 0.0007248 vs target 0.0007
    kv = fluid.kv_from_geometry(BASELINE_PARAMS.x_stroke, BASELINE_PARAMS)
    assert math.isclose(cv, 1.156 * kv, rel_tol=1e-12)


def test_kv_definition_water_reference_consistency():
    # Kv must equal the water volumetric flow (m^3/h) at its defining condition (1 bar, water).
    params = BASELINE_PARAMS
    x = params.x_stroke
    rho_ref = 1000.0
    mdot_ref = params.C_d * fluid.effective_area(x, params) * math.sqrt(2 * rho_ref * 1e5)
    q_m3_per_h = mdot_ref / rho_ref * 3600
    assert math.isclose(fluid.kv_from_geometry(x, params), q_m3_per_h, rel_tol=1e-12)


def test_flow_force_liquid_reference():
    cond = fluid.FlowConditions(P_up=2.0e5, P_down=1.0e5, T0=293.0)
    f = fluid.flow_force_liquid(1.0e-4, BASELINE_PARAMS, cond)
    assert math.isclose(f, 2.4630086404143978e-3, rel_tol=1e-9)


def test_schematic_script_produces_png():
    from solenoid_model.report import generate_schematic
    out = generate_schematic.generate()
    assert out.exists()
    assert out.stat().st_size > 10_000  # a real rendered figure, not an empty file


def test_guards_negative_lift_clamps_to_zero_area():
    assert fluid.effective_area(-1e-5, BASELINE_PARAMS) == 0.0


def test_guards_gas_zero_or_reversed_pressure_gives_zero_flow():
    p = BASELINE_PARAMS
    assert fluid.mdot_gas(p.x_stroke, p, fluid.N2, fluid.FlowConditions(0.0, 0.0, 293.0)) == 0.0
    assert fluid.mdot_gas(p.x_stroke, p, fluid.N2, fluid.FlowConditions(-1e5, 0.0, 293.0)) == 0.0
    # 逆壓差（P_down > P_up）：不模擬逆流，回傳 0
    assert fluid.mdot_gas(p.x_stroke, p, fluid.N2, fluid.FlowConditions(1e5, 2e5, 293.0)) == 0.0
    # 負的下游壓力視為 0（真空下限）
    choked_ref = fluid.mdot_gas(p.x_stroke, p, fluid.N2, fluid.FlowConditions(2.4e6, 0.0, 293.0))
    same = fluid.mdot_gas(p.x_stroke, p, fluid.N2, fluid.FlowConditions(2.4e6, -1e4, 293.0))
    assert math.isclose(same, choked_ref, rel_tol=1e-12)


def test_guards_liquid_reversed_dp_gives_zero_flow_and_force():
    p = BASELINE_PARAMS
    rev = fluid.FlowConditions(1e5, 2e5, 293.0)
    assert fluid.mdot_liquid(p.x_stroke, p, fluid.WATER_20C, rev) == 0.0
    assert fluid.flow_force_liquid(p.x_stroke, p, rev) == 0.0
    # 負的下游壓力視為 0（真空下限），與氣體路徑鉗位語意一致
    vac = fluid.FlowConditions(2.0e5, 0.0, 293.0)
    neg = fluid.FlowConditions(2.0e5, -1e4, 293.0)
    assert math.isclose(fluid.mdot_liquid(1.0e-4, p, fluid.WATER_20C, neg),
                        fluid.mdot_liquid(1.0e-4, p, fluid.WATER_20C, vac), rel_tol=1e-12)
    assert math.isclose(fluid.flow_force_liquid(1.0e-4, p, neg),
                        fluid.flow_force_liquid(1.0e-4, p, vac), rel_tol=1e-12)


def test_guards_do_not_change_legal_domain_values():
    p = BASELINE_PARAMS
    cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
    assert math.isclose(fluid.mdot_gas(p.x_stroke, p, fluid.N2, cond),
                        6.862805689539614e-5, rel_tol=1e-9)
    liq = fluid.FlowConditions(P_up=2.0e5, P_down=1.0e5, T0=293.0)
    assert math.isclose(fluid.mdot_liquid(1.0e-4, p, fluid.WATER_20C, liq),
                        1.74004195666987e-4, rel_tol=1e-9)
    assert math.isclose(fluid.cv_from_geometry(p.x_stroke, p), 0.0007, rel_tol=0.04)


def test_flow_force_gas_choked_vacuum_reference():
    cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
    f = fluid.flow_force_gas(BASELINE_PARAMS.x_stroke, BASELINE_PARAMS, fluid.N2, cond)
    assert math.isclose(f, 0.04137697092552667, rel_tol=1e-9)


def test_flow_force_gas_subsonic_reference_and_vanishing_limit():
    p = BASELINE_PARAMS
    f09 = fluid.flow_force_gas(p.x_stroke, p, fluid.N2,
                               fluid.FlowConditions(2.4e6, 0.9 * 2.4e6, 293.0))
    assert math.isclose(f09, 0.005690504306347273, rel_tol=1e-9)
    tiny = fluid.flow_force_gas(p.x_stroke, p, fluid.N2,
                                fluid.FlowConditions(2.4e6, 2.4e6 * 0.999999, 293.0))
    assert tiny < 1e-4


def test_flow_force_gas_continuous_at_critical_ratio():
    p = BASELINE_PARAMS
    r_crit = fluid.critical_pressure_ratio(fluid.N2)
    below = fluid.flow_force_gas(p.x_stroke, p, fluid.N2,
                                 fluid.FlowConditions(2.4e6, 2.4e6 * (r_crit - 1e-9), 293.0))
    above = fluid.flow_force_gas(p.x_stroke, p, fluid.N2,
                                 fluid.FlowConditions(2.4e6, 2.4e6 * (r_crit + 1e-9), 293.0))
    assert math.isclose(below, above, rel_tol=1e-6)


def test_flow_force_gas_zero_at_zero_lift_and_reversed_pressure():
    p = BASELINE_PARAMS
    cond = fluid.FlowConditions(2.4e6, 0.0, 293.0)
    assert fluid.flow_force_gas(0.0, p, fluid.N2, cond) == 0.0
    assert fluid.flow_force_gas(p.x_stroke, p, fluid.N2,
                                fluid.FlowConditions(1e5, 2e5, 293.0)) == 0.0


def test_p1_coupled_report_produces_comparison():
    from pathlib import Path
    from solenoid_model.report import generate_p1_coupled_report
    out = generate_p1_coupled_report.generate()
    assert out["t_open_coupled"] >= out["t_open_dry"]
    assert 0.0 < out["t_open_dry"] < 0.010
    png = Path(generate_p1_coupled_report.__file__).parent / "dry_vs_coupled_x_t.png"
    assert png.exists() and png.stat().st_size > 10_000
