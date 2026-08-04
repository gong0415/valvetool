import math
from pathlib import Path

from solenoid_model.report import generate_p0_report, generate_p2_l6_envelope, generate_p3_flyback_report, generate_p4_l2_pareto_report, generate_p4_l1_report, generate_p5_l3_report, generate_p6_report


def test_generate_report_produces_plots_and_reasonable_t_open():
    t_open = generate_p0_report.generate()
    assert 0.0 < t_open < 0.010

    report_dir = Path(generate_p0_report.__file__).parent
    assert (report_dir / "current_i_t.png").exists()
    assert (report_dir / "position_x_t.png").exists()
    assert (report_dir / "force_F_t.png").exists()


def test_generate_l6_envelope_report_bare_selfheat_and_ln2(tmp_path):
    # 粗網格/粗容差控制測試時長(~30 秒);單調性不受影響
    # out_path 走 tmp:測試不覆寫版本庫的全解析度 PNG
    out_png = tmp_path / "l6_envelope_v_t.png"
    out = generate_p2_l6_envelope.generate(n_points=3, tol=0.5, t_max=0.02,
                                           out_path=out_png)
    assert out_png.exists()

    v_bare = out["v_pi_bare"]
    assert all(a < b for a, b in zip(v_bare, v_bare[1:]))  # L6:V_pull-in 隨 T 上升
    # 含自熱曲線必在裸曲線上方(T_hot > T_amb 且 V_pi 對 T 單調)
    assert all(s >= b for s, b in zip(out["v_pi_selfheat"], v_bare))
    assert all(t_hot > t_amb for t_hot, t_amb in zip(out["T_selfheat"], out["T_grid"]))
    # LN₂ 77 K:域外低電阻 → 拉入電壓遠低;穩態電流大 → 飽和旗標必觸發
    assert out["v_pi_ln2"] < v_bare[0]
    assert out["ln2_saturated"] is True
    # 響應特性 4 表:三溫 × 兩電壓全部開得了
    assert len(out["resp4"]) == 6
    for _T, _V, t_open in out["resp4"]:
        assert t_open is not None


def test_generate_p3_flyback_report_zener_faster_than_diode_and_in_spec_range(tmp_path):
    out_png = tmp_path / "p3_flyback_current_t.png"
    out = generate_p3_flyback_report.generate(out_path=out_png)
    assert out_png.exists()

    results = out["results"]
    assert results["zener"].t_close < results["diode"].t_close
    assert 3.0 <= out["diode_to_zener_ratio"] <= 5.0


def test_generate_p4_l2_pareto_report_baseline_saturated_and_ceiling(tmp_path):
    out_png = tmp_path / "p4_l2_pareto.png"
    out = generate_p4_l2_pareto_report.generate(out_path=out_png)
    assert out_png.exists()

    assert math.isclose(out["F_sat"], 43.09120084213066, rel_tol=1e-9)
    # baseline hold sits in the saturated flat region
    assert out["baseline_saturated"] is True
    # ceiling starts binding at a hold power far below the baseline's 9.8 W
    assert 0.0 < out["P_cross"] < 1.0
    # front is index-aligned and every F_hold respects the ceiling
    front = out["front"]
    assert len(front["N"]) == len(front["F_hold"])
    assert all(f <= out["F_sat"] + 1e-9 for f in front["F_hold"])


def test_generate_p4_l1_report_bound_valid_and_naive_is_not(tmp_path):
    out_png = tmp_path / "p4_l1_sweeps.png"
    out = generate_p4_l1_report.generate(out_path=out_png)
    assert out_png.exists()

    # the bound holds everywhere and is tight
    lo, hi = out["bound_ratio_range"]
    assert 0.70 <= lo <= hi <= 1.00

    # the naive estimator lands on BOTH sides of the truth -> not a bound
    n_lo, n_hi = out["naive_ratio_range"]
    assert n_lo < 1.0 < n_hi

    # the two turns curves disagree: window-consistent N=4000 cannot pull in
    # while the fixed-R curve claims it opens
    assert out["turns_window"]["pulls_in"][-1] is False
    assert out["turns_fixed_R"]["pulls_in"][-1] is True


def test_generate_p5_l3_report_limits_and_pi_2_trend(tmp_path):
    out_png = tmp_path / "p5_l3_region.png"
    out = generate_p5_l3_report.generate(out_path=out_png)
    assert out_png.exists()

    assert math.isclose(out["s_star"], 1.0139509758476493, rel_tol=1e-9)
    assert math.isclose(out["dP_max"], 8618240.168426132, rel_tol=1e-9)
    assert math.isclose(out["groups"]["pi_2_time_ratio"], 23.25455910965206,
                        rel_tol=1e-9)

    # pi_2 falls monotonically with scale -- smaller valves are less
    # dead-time dominated
    pi2 = [v for _s, v in out["pi_2_by_scale"]]
    assert all(a > b for a, b in zip(pi2, pi2[1:]))

    # the feasible region is a rectangle: every row is a prefix of Trues
    for row in out["region"]["feasible"]:
        assert row == sorted(row, reverse=True)


def test_generate_p6_report_land_limits_and_the_zener_impact_bill(tmp_path):
    out = generate_p6_report.generate(
        l4_path=tmp_path / "p6_l4_leak.png",
        l5_path=tmp_path / "p6_l5_life.png",
        bounce_path=tmp_path / "p6_bounce_x_t.png")
    for key in ("l4", "l5", "bounce"):
        assert out[key].exists()
    s = out["summary"]
    # soft seat seals with a manufacturable land, metal-to-metal does not
    assert 60.0 < s["pctfe_land_limit_um"] < 70.0
    assert s["metal_land_limit_um"] < 2.0
    # P3's zener speedup is paid for in impact energy
    assert s["zener_energy_ratio"] > 10.0
    assert s["v_zener"] > s["v_diode"]
    assert s["n_bounce"] >= 2


def test_generate_winding_window_schematic_produces_png():
    from solenoid_model.report import generate_schematic
    out = generate_schematic.generate_winding_window()
    assert out.name == "winding_window_schematic.png"
    assert out.exists() and out.stat().st_size > 0


def test_generate_seal_land_schematic_produces_png():
    from solenoid_model.report import generate_schematic
    out = generate_schematic.generate_seal_land()
    assert out.name == "seal_land_schematic.png"
    assert out.exists() and out.stat().st_size > 0
