"""L6 pull-in-voltage search + L2 power-vs-holding-force Pareto (P4).
L6 anchors are deterministic bisection endpoints: with tol=0.5, t_max=0.02
the sequence from [0, 28] lands on exact binary fractions. L2 tests check
the B_sat force ceiling and the winding-slope-vs-ceiling Pareto front.

Also covers L1 response-time analysis (P4): the analytic
building blocks (electrical_time_constant, motion_threshold_current,
response_time_bound, naive_response_time_estimate), the l1_sweep helper,
and the module-consistency guards that check the B_sat force ceiling and
the resisting-force thresholds never bind ahead of where the ODE actually
operates."""
import math
from dataclasses import replace

import pytest

from solenoid_model import dynamics, limits, sizing, winding
from solenoid_model.baseline import BASELINE_PARAMS


def test_pullin_voltage_room_temp_anchor():
    v = limits.find_pullin_voltage(BASELINE_PARAMS, tol=0.5, t_max=0.02)
    assert math.isclose(v, 16.625, abs_tol=1e-9)


def test_pullin_voltage_rises_with_temperature():
    # planning 定案:13.125 (−40°C) < 16.625 (20°C) < 19.25 (+70°C)
    v_cold = limits.find_pullin_voltage(BASELINE_PARAMS, T_coil=233.15, tol=0.5, t_max=0.02)
    v_room = limits.find_pullin_voltage(BASELINE_PARAMS, tol=0.5, t_max=0.02)
    v_hot = limits.find_pullin_voltage(BASELINE_PARAMS, T_coil=343.15, tol=0.5, t_max=0.02)
    assert math.isclose(v_cold, 13.125, abs_tol=1e-9)
    assert math.isclose(v_hot, 19.25, abs_tol=1e-9)
    assert v_cold < v_room < v_hot


def test_pullin_voltage_returns_successful_endpoint():
    # 回傳值本身開得了;差一個 tol 的下方開不了(判準語意)
    v = limits.find_pullin_voltage(BASELINE_PARAMS, tol=0.5, t_max=0.02)
    from solenoid_model import dynamics
    sol_at = dynamics.simulate_opening(replace(BASELINE_PARAMS, V_bus=v), t_max=0.02)
    sol_below = dynamics.simulate_opening(replace(BASELINE_PARAMS, V_bus=v - 0.5), t_max=0.02)
    assert sol_at.t_events[0].size == 1
    assert sol_below.t_events[0].size == 0


def test_pullin_voltage_raises_when_valve_cannot_pull_in():
    p = replace(BASELINE_PARAMS, F_preload=10000.0)
    with pytest.raises(ValueError, match="no pull-in"):
        limits.find_pullin_voltage(p, t_max=0.005)


def test_pullin_voltage_rejects_nonpositive_v_bus():
    # 現行為:hi = V_bus 的倍增迴圈在 0 卡死成無窮迴圈、負值行為未定義
    for v in (0.0, -28.0):
        with pytest.raises(ValueError, match="V_bus"):
            limits.find_pullin_voltage(replace(BASELINE_PARAMS, V_bus=v),
                                       t_max=0.005)


def test_saturation_force_ceiling_value():
    F_sat = limits.saturation_force_ceiling(BASELINE_PARAMS)
    assert math.isclose(F_sat, 43.09120084213066, rel_tol=1e-9)


def test_holding_force_is_capped_at_saturation_not_fictitious_linear():
    # baseline hold current 0.35 A: linear magnetics predicts ~2364 N, but
    # the real force is capped by B_sat at ~43 N (coil is deep in saturation)
    p = BASELINE_PARAMS
    i = p.V_bus / p.R_coil_20C  # 0.35 A
    F = limits.holding_force(p, i)
    assert math.isclose(F, 43.09120084213066, rel_tol=1e-9)
    assert F < 100.0  # emphatically not the fictitious 2364 N


def test_holding_force_unsaturated_returns_linear():
    # a tiny current stays below the ceiling -> returns the linear force
    p = BASELINE_PARAMS
    i_small = 1e-3
    F = limits.holding_force(p, i_small)
    F_sat = limits.saturation_force_ceiling(p)
    assert F < F_sat
    # matches 0.5*i^2*|dL/dx| at the pulled-in gap
    from solenoid_model import magnetics
    dLdx = -magnetics.dinductance_dgap(p.g0 - p.x_stroke, p)
    assert math.isclose(F, 0.5 * i_small**2 * abs(dLdx), rel_tol=1e-12)


def test_l2_pareto_front_shape_and_saturation_region():
    p = BASELINE_PARAMS
    N_values = [2000.0, 8000.0, 20000.0]
    front = limits.l2_pareto_front(p, N_values)
    assert front["N"] == N_values
    F_sat = limits.saturation_force_ceiling(p)
    # each F_hold is min(F_linear, F_sat)
    for F_hold, F_lin in zip(front["F_hold"], front["F_linear"]):
        assert math.isclose(F_hold, min(F_lin, F_sat), rel_tol=1e-9)
    # N=2000 baseline sits in the saturated (capped) region
    assert math.isclose(front["F_hold"][0], F_sat, rel_tol=1e-9)
    # N=20000 (fine wire, low power) drops below the ceiling
    assert front["F_hold"][2] < F_sat
    # power P = i^2 R with i = V/R and R from the winding model
    R0 = winding.coil_resistance(2000.0, p.A_winding, p.k_fill, p.l_turn_mean)
    assert math.isclose(front["R"][0], R0, rel_tol=1e-12)
    assert math.isclose(front["P_hold"][0], (p.V_bus / R0) ** 2 * R0, rel_tol=1e-9)


def test_electrical_time_constant_at_rest_and_pulled_in_gaps():
    p = BASELINE_PARAMS
    assert math.isclose(limits.electrical_time_constant(p),
                        0.005199877495596899, rel_tol=1e-9)
    assert math.isclose(limits.electrical_time_constant(p, gap=p.g0 - p.x_stroke),
                        0.030159289474462003, rel_tol=1e-9)


def test_electrical_time_constant_is_invariant_with_turns_at_fixed_window():
    # L ~ N^2 and (winding-derived) R ~ N^2, so tau_e = L/R cancels N entirely:
    # turns is NOT a lever on the electrical time constant at a fixed window.
    p = BASELINE_PARAMS
    for N in (1000.0, 2000.0, 4000.0, 8000.0):
        p_n = replace(p, N_turns=N)
        R = winding.coil_resistance(N, p.A_winding, p.k_fill, p.l_turn_mean)
        assert math.isclose(limits.electrical_time_constant(p_n, R=R),
                            0.005199877495596899, rel_tol=1e-9)


def test_motion_threshold_current_balances_spring_and_pressure():
    p = BASELINE_PARAMS
    i_th = limits.motion_threshold_current(p)
    assert math.isclose(i_th, 0.19581177043586317, rel_tol=1e-9)
    # at i_th the magnetic force at the rest gap exactly equals the resistance
    from solenoid_model import magnetics
    F_mag = 0.5 * i_th**2 * abs(magnetics.dinductance_dgap(p.g0, p))
    F_resist = dynamics.spring_force(0.0, p) + p.delta_P * p.A_seat
    assert math.isclose(F_mag, F_resist, rel_tol=1e-9)
    assert math.isclose(F_resist, 22.0, rel_tol=1e-12)


def test_response_time_bound_components_and_total():
    p = BASELINE_PARAMS
    b = limits.response_time_bound(p)
    assert math.isclose(b["t_i"], 0.0042626465242997305, rel_tol=1e-9)
    assert math.isclose(b["t_x"], 0.00011799981508206165, rel_tol=1e-9)
    assert math.isclose(b["t_bound"], 0.004380646339381792, rel_tol=1e-9)
    assert math.isclose(b["t_bound"], b["t_i"] + b["t_x"], rel_tol=1e-12)


def test_response_time_bound_is_below_simulated_t_open():
    p = BASELINE_PARAMS
    b = limits.response_time_bound(p)
    t_open = dynamics.simulate_opening(p).t_events[0][0]
    assert b["t_bound"] < t_open
    assert math.isclose(b["t_bound"] / t_open, 0.8714345527671004, rel_tol=1e-6)


def test_response_time_bound_returns_none_when_valve_cannot_pull_in():
    # i_ss below the motion threshold -> never opens, no bound to report
    p = replace(BASELINE_PARAMS, V_bus=10.0)
    assert limits.response_time_bound(p) is None


def test_naive_estimate_value_and_that_it_is_not_a_bound():
    # The spec's tau_e + sqrt(2mx/F) decomposition is an order-of-magnitude
    # estimator, NOT a bound: it lands on both sides of the simulated t_open.
    p = BASELINE_PARAMS
    assert math.isclose(limits.naive_response_time_estimate(p),
                        0.005311346986978533, rel_tol=1e-9)

    p_low = replace(p, V_bus=18.0)
    p_high = replace(p, V_bus=50.0)
    n_low = limits.naive_response_time_estimate(p_low)
    n_high = limits.naive_response_time_estimate(p_high)
    t_low = dynamics.simulate_opening(p_low).t_events[0][0]
    t_high = dynamics.simulate_opening(p_high).t_events[0][0]
    assert n_low < t_low    # under-predicts at low voltage
    assert n_high > t_high  # over-predicts at high voltage


def _l1_sweep_cases():
    """Three reduced sweeps (3 points each) covering voltage, turns
    (winding-consistent R) and stroke. The report script runs these at full
    resolution; the test keeps 9 simulations to stay fast."""
    p = BASELINE_PARAMS
    cases = []
    for V in (22.0, 28.0, 36.0):
        cases.append((f"V={V:g}", replace(p, V_bus=V), p.R_coil_20C))
    for N in (1000.0, 2000.0, 3000.0):
        R = winding.coil_resistance(N, p.A_winding, p.k_fill, p.l_turn_mean)
        cases.append((f"N={N:g}", replace(p, N_turns=N), R))
    for x in (2.0e-4, 3.0e-4, 4.0e-4):
        cases.append((f"x={x*1e3:g}mm",
                      replace(p, x_stroke=x, g0=x + 0.5e-4), p.R_coil_20C))
    return cases


def test_l1_sweep_returns_index_aligned_lists():
    cases = _l1_sweep_cases()
    out = limits.l1_sweep(cases)
    n = len(cases)
    for key in ("label", "t_open", "t_bound", "t_naive", "pulls_in"):
        assert len(out[key]) == n
    assert out["label"] == [c[0] for c in cases]
    assert all(out["pulls_in"])  # every reduced case opens


def test_bound_holds_below_t_open_across_every_sweep():
    # The central correctness claim of the L1 analysis: the bound is a real
    # lower bound on the simulated opening time in every swept case, and is
    # tight enough to be useful (0.70-1.00 of the actual).
    out = limits.l1_sweep(_l1_sweep_cases())
    for label, t_open, t_bound in zip(out["label"], out["t_open"], out["t_bound"]):
        assert t_bound < t_open, f"bound not below t_open at {label}"
        assert 0.70 <= t_bound / t_open <= 1.00, f"bound too loose at {label}"


def test_l1_sweep_marks_non_opening_cases():
    p = BASELINE_PARAMS
    # window-consistent N=4000 -> R=320 ohm -> steady current far too low
    R = winding.coil_resistance(4000.0, p.A_winding, p.k_fill, p.l_turn_mean)
    out = limits.l1_sweep([("N=4000", replace(p, N_turns=4000.0), R)])
    assert out["pulls_in"] == [False]
    assert out["t_open"] == [None]
    assert out["t_bound"] == [None]


def test_saturation_cap_does_not_bind_at_the_decisive_thresholds():
    # Guard: limits.holding_force caps force at B_sat, but magnetics/dynamics
    # (which drive the ODE) do not. That division of labour is only sound
    # while the forces that actually decide the transients stay below the
    # ceiling -- the force needed to start moving, and the force at which the
    # armature releases. If someone lowers B_sat or raises preload/delta_P
    # past this point, the ODE starts using fictitious uncapped forces and
    # this test must fail loudly rather than let it pass silently.
    p = BASELINE_PARAMS
    F_sat = limits.saturation_force_ceiling(p)
    F_start = dynamics.spring_force(0.0, p) + p.delta_P * p.A_seat
    F_release = dynamics.spring_force(p.x_stroke, p) + p.delta_P * p.A_seat
    assert math.isclose(F_start, 22.0, rel_tol=1e-12)
    assert math.isclose(F_release, 28.0, rel_tol=1e-12)
    assert math.isclose(F_sat, 43.09120084213066, rel_tol=1e-9)
    assert F_start < F_sat
    assert F_release < F_sat

    # The quantity that actually gates the bound's validity is the peak
    # magnetic force the (uncapped) ODE develops during the stroke -- not the
    # resisting-force thresholds above. Assert the ceiling clears that too,
    # otherwise response_time_bound's t_x leg (which uses the capped force)
    # stops being the optimistic term the bound argument requires.
    from solenoid_model import magnetics
    sol = dynamics.simulate_opening(p)
    F_peak = max(magnetics.magnetic_force_closing(p.g0 - x, i, p)
                 for i, x in zip(sol.y[0], sol.y[1]))
    assert math.isclose(F_peak, 28.878, rel_tol=1e-3)
    assert F_peak < F_sat


def test_dimensionless_groups_baseline_values():
    g = limits.dimensionless_groups(BASELINE_PARAMS)
    assert math.isclose(g["pi_1_stroke_gap"], 0.8571428571428571, rel_tol=1e-9)
    assert math.isclose(g["pi_2_time_ratio"], 23.25455910965206, rel_tol=1e-9)
    assert math.isclose(g["pi_3_magnetic_number"], 4.39299995317669, rel_tol=1e-9)
    assert math.isclose(g["pi_4_pressure_load"], 0.7500000000000001, rel_tol=1e-9)
    assert math.isclose(g["pi_5_saturation_margin"], 3.590933403510888, rel_tol=1e-9)


def test_electrical_and_mechanical_time_constants_scale_as_d2_and_d():
    # tau_e ~ D^2 and tau_m ~ D under self-similar scaling, so their ratio ~ D.
    p = BASELINE_PARAMS
    N = p.N_turns
    for s in (1.0, 0.5, 0.2):
        ps = sizing.scale_params(p, s)
        R = winding.coil_resistance(N, ps.A_winding, ps.k_fill, ps.l_turn_mean)
        tau_e = limits.electrical_time_constant(ps, R=R)
        tau_m = 1.0 / math.sqrt(ps.k_spring / ps.m_arm)
        assert math.isclose(tau_e / s**2, 0.005199877495596899, rel_tol=1e-9)
        assert math.isclose(tau_m / s, 0.00022360679774997895, rel_tol=1e-9)


def test_pi_2_is_the_only_scale_dependent_group_and_falls_with_size():
    # pi_2 ~ D: a smaller valve shifts from dead-time dominated toward
    # mechanically dominated. This predicts the L1 analytic bound (whose
    # tightness rests on dead-time domination) should loosen at small scale.
    p = BASELINE_PARAMS
    N = p.N_turns
    ratios = []
    for s in (1.0, 0.5, 0.2):
        ps = sizing.scale_params(p, s)
        # R MUST come from the scaled winding window, not the untouched
        # params.R_coil_20C default: scale_params deliberately leaves
        # R_coil_20C unscaled (see its docstring), so a scaled params object
        # is winding-inconsistent by construction. pi_3's invariance below
        # holds only under this winding-consistent R -- swap in the default
        # R_coil_20C and pi_3 blows up as D**-2 instead of staying put (see
        # limits.dimensionless_groups's docstring for the contrast).
        R = winding.coil_resistance(N, ps.A_winding, ps.k_fill, ps.l_turn_mean)
        g = limits.dimensionless_groups(ps, R=R)
        ratios.append(g["pi_2_time_ratio"])
        # the other four groups stay put
        assert math.isclose(g["pi_1_stroke_gap"], 0.8571428571428571, rel_tol=1e-9)
        assert math.isclose(g["pi_3_magnetic_number"], 4.39299995317669, rel_tol=1e-9)
        assert math.isclose(g["pi_4_pressure_load"], 0.7500000000000001, rel_tol=1e-9)
        assert math.isclose(g["pi_5_saturation_margin"], 3.590933403510888, rel_tol=1e-9)
    assert math.isclose(ratios[0], 23.25455910965206, rel_tol=1e-9)
    assert math.isclose(ratios[1], 11.62727955482603, rel_tol=1e-9)
    assert math.isclose(ratios[2], 4.650911821930412, rel_tol=1e-9)
    assert ratios[0] > ratios[1] > ratios[2]
    for s, r in zip((1.0, 0.5, 0.2), ratios):
        assert math.isclose(r / s, 23.25455910965206, rel_tol=1e-9)


def test_l4_leak_curve_brackets_the_percolation_threshold():
    from solenoid_model import fluid, limits, materials
    from solenoid_model.baseline import BASELINE_PARAMS as P
    cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
    seats = [materials.PCTFE_ON_440C, materials.METAL_17_4PH_ON_440C]
    w_lands = [20e-6, 50e-6, 100e-6, 200e-6]
    out = limits.l4_leak_curve(P, seats, w_lands, cond)
    assert set(out) == {s.name for s in seats}
    for rows in out.values():
        assert len(rows) == len(w_lands)
        for r in rows:
            assert (r["leak"] == 0.0) == r["sealed"]
        # Among unsealed rows the leak FALLS as the land widens: the gap
        # closure saturates at Rq_c while the flow path keeps lengthening.
        # That is the design doc's finding-2 artefact, pinned here so the
        # known limitation stays visible rather than being mistaken for
        # physics.
        unsealed = [r["leak"] for r in rows if not r["sealed"]]
        assert unsealed == sorted(unsealed, reverse=True)
    pctfe = out[materials.PCTFE_ON_440C.name]
    assert pctfe[0]["sealed"] and pctfe[1]["sealed"]      # 20, 50 um seal
    assert not pctfe[3]["sealed"]                          # 200 um does not
    metal = out[materials.METAL_17_4PH_ON_440C.name]
    assert not any(r["sealed"] for r in metal)


def test_l5_life_curve_is_monotonic_and_flags_shakedown():
    from solenoid_model import limits, materials
    from solenoid_model.baseline import BASELINE_PARAMS as P
    rows = limits.l5_life_curve(
        P, materials.PCTFE_ON_440C, [0.05, 0.1, 0.4, 1.5])
    assert len(rows) == 4
    lives = [r["N_cycle"] for r in rows]
    assert lives == sorted(lives, reverse=True)   # faster impact -> shorter life
    assert not rows[-1]["shakedown"]              # 1.5 m/s is plastic


def test_l5_life_curve_reports_infinite_life_under_shakedown():
    import math

    from solenoid_model import limits, materials
    from solenoid_model.baseline import BASELINE_PARAMS as P
    rows = limits.l5_life_curve(P, materials.METAL_17_4PH_ON_440C, [1e-3])
    assert rows[0]["shakedown"]
    assert math.isinf(rows[0]["N_cycle"])


def test_bounce_leak_spike_reports_off_seat_leakage():
    from solenoid_model import fluid, flyback, limits, materials
    from solenoid_model.baseline import BASELINE_PARAMS as P
    cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
    out = limits.bounce_leak_spike(
        P, materials.PCTFE_ON_440C, flyback.ZenerFlyback(), cond, fluid.N2)
    assert out["n_bounce"] >= 2
    assert out["t_open_total"] > 0.0
    assert out["mass_leaked"] > 0.0
