import math
from dataclasses import replace

from solenoid_model import sizing, thermal, winding
from solenoid_model.baseline import BASELINE_PARAMS


def test_active_mass_is_core_plus_coil():
    p = BASELINE_PARAMS
    m = sizing.active_mass(p)
    m_core = p.A_gap * p.l_core * sizing.RHO_FE
    m_coil = p.A_winding * p.l_turn_mean * sizing.RHO_CU
    assert math.isclose(m_core, 0.011805, rel_tol=1e-9)
    assert math.isclose(m_coil, 0.014527416083237826, rel_tol=1e-9)
    assert math.isclose(m, m_core + m_coil, rel_tol=1e-12)
    assert math.isclose(m, 0.026332416083237824, rel_tol=1e-9)


def test_valve_mass_applies_packaging_factor_and_lands_in_spec_range():
    p = BASELINE_PARAMS
    m = sizing.valve_mass(p, k_pack=3.0)
    assert math.isclose(m, 0.07899724824971346, rel_tol=1e-9)
    # spec §4 lists m_valve 30-500 g; the calibrated baseline must land inside
    assert 0.030 <= m <= 0.500


def test_pressure_ceiling_matches_maxwell_saturation_limit():
    p = BASELINE_PARAMS
    dP = sizing.pressure_ceiling(p)
    assert math.isclose(dP, 8618240.168426132, rel_tol=1e-9)
    # equivalently: the ceiling force divided by the sealing area
    expected = p.B_sat**2 * p.A_gap / (2.0 * sizing.MU_0 * p.A_seat)
    assert math.isclose(dP, expected, rel_tol=1e-12)


def test_min_feasible_scale_falls_with_allowed_current_density():
    p = BASELINE_PARAMS
    assert math.isclose(sizing.min_feasible_scale(p, 10e6),
                        2.0279019516952985, rel_tol=1e-9)
    assert math.isclose(sizing.min_feasible_scale(p, 20e6),
                        1.0139509758476493, rel_tol=1e-9)
    assert math.isclose(sizing.min_feasible_scale(p, 30e6),
                        0.6759673172317662, rel_tol=1e-9)


def test_min_feasible_scale_is_where_available_mmf_equals_required():
    # at s = s*, MMF available from the window exactly meets what the gap needs
    p = BASELINE_PARAMS
    J_max = 20e6
    s_star = sizing.min_feasible_scale(p, J_max)
    NI_need = p.B_sat * (p.g0 * s_star) / sizing.MU_0
    NI_avail = J_max * (p.A_winding * s_star**2) * p.k_fill
    assert math.isclose(NI_avail, NI_need, rel_tol=1e-9)


def test_scale_params_applies_the_documented_exponents():
    p = BASELINE_PARAMS
    s = 0.5
    ps = sizing.scale_params(p, s)
    assert math.isclose(ps.g0, p.g0 * s, rel_tol=1e-12)
    assert math.isclose(ps.x_stroke, p.x_stroke * s, rel_tol=1e-12)
    assert math.isclose(ps.l_core, p.l_core * s, rel_tol=1e-12)
    assert math.isclose(ps.l_turn_mean, p.l_turn_mean * s, rel_tol=1e-12)
    assert math.isclose(ps.D_seat_bore, p.D_seat_bore * s, rel_tol=1e-12)
    assert math.isclose(ps.A_gap, p.A_gap * s**2, rel_tol=1e-12)
    assert math.isclose(ps.A_seat, p.A_seat * s**2, rel_tol=1e-12)
    assert math.isclose(ps.A_winding, p.A_winding * s**2, rel_tol=1e-12)
    assert math.isclose(ps.A_rad, p.A_rad * s**2, rel_tol=1e-12)
    assert math.isclose(ps.m_arm, p.m_arm * s**3, rel_tol=1e-12)
    assert math.isclose(ps.k_spring, p.k_spring * s, rel_tol=1e-12)
    assert math.isclose(ps.F_preload, p.F_preload * s**2, rel_tol=1e-12)
    assert math.isclose(ps.G_th_cond, p.G_th_cond * s, rel_tol=1e-12)
    # untouched
    assert math.isclose(ps.B_sat, p.B_sat, rel_tol=1e-12)
    assert math.isclose(ps.V_bus, p.V_bus, rel_tol=1e-12)
    assert math.isclose(ps.k_fill, p.k_fill, rel_tol=1e-12)


def test_pressure_ceiling_is_scale_invariant():
    # A_gap and A_seat both go as D^2 and cancel: shrinking the valve does
    # not lower its pressure rating.
    p = BASELINE_PARAMS
    base = sizing.pressure_ceiling(p)
    for s in (1.0, 0.5, 0.25):
        assert math.isclose(sizing.pressure_ceiling(sizing.scale_params(p, s)),
                            base, rel_tol=1e-12)
    assert math.isclose(base, 8618240.168426132, rel_tol=1e-9)


def test_force_margin_is_scale_invariant():
    # The spec names seat dP force and magnetic saturation as L3 mechanisms,
    # but both go as D^2, so their ratio is the same at every scale and they
    # produce NO miniaturisation limit on their own.
    p = BASELINE_PARAMS
    for s in (1.0, 0.5, 0.25, 0.1):
        ps = sizing.scale_params(p, s)
        F_avail = ps.B_sat**2 * ps.A_gap / (2.0 * sizing.MU_0)
        F_req = ps.delta_P * ps.A_seat
        assert math.isclose(F_avail / F_req, 3.590933403510888, rel_tol=1e-9)


def test_active_mass_scales_as_the_cube():
    p = BASELINE_PARAMS
    for s in (1.0, 0.5, 0.25):
        m = sizing.active_mass(sizing.scale_params(p, s))
        assert math.isclose(m / s**3, 0.026332416083237824, rel_tol=1e-9)


def test_l3_feasible_region_shape_and_alignment():
    p = BASELINE_PARAMS
    scales = [1.5, 1.0, 0.5]
    pressures = [1.0e6, 5.0e6]
    out = sizing.l3_feasible_region(p, J_max=20e6, k_pack=3.0,
                                    scales=scales, pressures=pressures)
    assert out["scales"] == scales
    assert out["pressures"] == pressures
    assert len(out["feasible"]) == len(scales)
    assert all(len(row) == len(pressures) for row in out["feasible"])
    assert len(out["mass"]) == len(scales)
    # mass depends only on scale, and grows with it
    assert out["mass"][0] > out["mass"][1] > out["mass"][2]


def test_l3_feasible_region_matches_the_closed_form_boundaries():
    # WIRING check, not independent validation: confirms l3_feasible_region
    # threads min_feasible_scale and pressure_ceiling into the grid correctly
    # (including boundary cells). By construction, this test cannot catch an
    # error inside those two closed forms — only a mismatch between what they
    # return and what appears in the grid. For genuine independent validation
    # against the underlying physics, see test_l3_feasible_region_agrees_with_a_direct_mmf_and_force_check.
    p = BASELINE_PARAMS
    J_max = 20e6
    s_star = sizing.min_feasible_scale(p, J_max)
    dP_max = sizing.pressure_ceiling(p)
    scales = [1.5, 1.1, 1.0139509758476493, 0.9, 0.5]
    pressures = [1.0e6, 8.0e6, 8618240.168426132, 9.0e6]
    out = sizing.l3_feasible_region(p, J_max=J_max, k_pack=3.0,
                                    scales=scales, pressures=pressures)
    for i, s in enumerate(scales):
        for j, dP in enumerate(pressures):
            assert out["feasible"][i][j] is bool(s >= s_star and dP <= dP_max)


def test_l3_feasible_region_agrees_with_a_direct_mmf_and_force_check():
    # Genuinely independent validation. Recompute feasibility from scale_params
    # and the underlying MMF (NI_avail >= NI_need) and Maxwell-force
    # (dP*A_seat <= F_avail) physics, touching neither min_feasible_scale nor
    # pressure_ceiling. This path is independent arithmetic; agreement proves
    # the closed-form limits capture the physics correctly.
    p = BASELINE_PARAMS
    J_max = 20e6
    scales = [1.5, 1.0, 0.5]
    pressures = [1.0e6, 9.0e6]
    out = sizing.l3_feasible_region(p, J_max=J_max, k_pack=3.0,
                                    scales=scales, pressures=pressures)
    for i, s in enumerate(scales):
        ps = sizing.scale_params(p, s)
        NI_need = ps.B_sat * ps.g0 / sizing.MU_0
        NI_avail = J_max * ps.A_winding * ps.k_fill
        F_avail = ps.B_sat**2 * ps.A_gap / (2.0 * sizing.MU_0)
        for j, dP in enumerate(pressures):
            direct = bool(NI_avail >= NI_need and dP * ps.A_seat <= F_avail)
            assert out["feasible"][i][j] is direct


# --- thermal-limited winding-window MMF: pins the load-bearing prose in
# docs/spec/PARAMS.md and the report's "force margin is scale-invariant"
# claim, which otherwise covers only 2 of the 3 spec-named L3 mechanisms
# (saturation, seat dP) and leaves the winding-window one resting on
# unshipped prose. ---

_T_AMB = 293.15   # K, matches the P2 thermal-model test convention
_T_CEIL = 373.15  # K = 100 degC, the design-doc prototype's coil ceiling


def _thermal_limited_ni_ratio(p, s):
    """Back-solve the thermal-limited steady current at scale s from the P2
    conduction+radiation model, and return (NI_avail/NI_need, current
    density [A/m**2]).

    scale_params deliberately does not rescale R_coil_20C (see Finding 1 /
    sizing.scale_params's docstring), so a scaled params object is
    winding-inconsistent by construction. R here is therefore derived from
    the SCALED winding window via winding.coil_resistance -- not read from
    the params default -- before it is handed to the thermal model.

    NI_need uses the same gap-only convention as min_feasible_scale (no
    iron term; see Finding 5), so the ratio is directly comparable to s*.
    """
    ps = sizing.scale_params(p, s)
    R20_scaled = winding.coil_resistance(ps.N_turns, ps.A_winding,
                                         ps.k_fill, ps.l_turn_mean)
    ps_r = replace(ps, R_coil_20C=R20_scaled)
    losses = (ps_r.G_th_cond * (_T_CEIL - _T_AMB)
              + ps_r.emissivity * thermal.STEFAN_BOLTZMANN * ps_r.A_rad
              * (_T_CEIL ** 4 - _T_AMB ** 4))
    I_thermal = math.sqrt(losses / thermal.R_coil(_T_CEIL, ps_r))
    NI_avail = ps.N_turns * I_thermal
    NI_need = p.B_sat * ps.g0 / sizing.MU_0
    a_wire = ps.A_winding * ps.k_fill / ps.N_turns
    J = I_thermal / a_wire
    return NI_avail / NI_need, J


def test_thermal_limited_mmf_ratio_stays_above_one_across_a_decade_of_scale():
    # Reconstructs the PARAMS.md claim: back-solving the thermal-limited
    # current at each scale and forming NI_avail/NI_need never drops below
    # 1 across a 10x scale range (s=1.0 down to s=0.1) -- i.e. the
    # winding-window mechanism alone gives no L3 limit either, matching the
    # scale-invariant force-margin result for saturation and seat dP.
    p = BASELINE_PARAMS
    ratio_s1, _ = _thermal_limited_ni_ratio(p, 1.0)
    ratio_s01, _ = _thermal_limited_ni_ratio(p, 0.1)
    assert math.isclose(ratio_s1, 1.5114027906905578, rel_tol=1e-9)
    assert math.isclose(ratio_s01, 1.4780021269186137, rel_tol=1e-9)
    for s in (1.0, 0.8, 0.6, 0.4, 0.2, 0.1):
        ratio, _ = _thermal_limited_ni_ratio(p, s)
        assert 1.4780021269186137 - 1e-6 <= ratio <= 1.5114027906905578 + 1e-6
        assert ratio > 1.0


def test_thermal_limited_current_density_at_s_0p1_is_an_extrapolation_artifact():
    # The "never below 1" band above is an extrapolation artifact, not
    # evidence the winding-window mechanism is truly limitless: the
    # back-solved current density at s=0.1 is ~299.72 A/mm**2, wildly
    # outside any real continuous copper rating (~5-20 A/mm**2 -- the
    # baseline itself already runs pulsed duty at 26.82 A/mm**2, see
    # test_min_feasible_scale_falls_with_allowed_current_density's J_max=20
    # case). The P2 thermal model has no notion of insulation breakdown or
    # conductor fusing -- the mechanisms that would actually fail up here --
    # so its equilibrium prediction in this region is not credible physics.
    # That extrapolation gap is precisely why the L3 winding-window limit
    # this repo ships comes from the EMPIRICAL J_max ceiling, not from this
    # thermal back-solve.
    p = BASELINE_PARAMS
    _, J = _thermal_limited_ni_ratio(p, 0.1)
    assert math.isclose(J, 299724339.7788058, rel_tol=1e-9)
    J_A_per_mm2 = J / 1e6
    assert math.isclose(J_A_per_mm2, 299.7243397788058, rel_tol=1e-9)
    assert J_A_per_mm2 > 20.0 * 10  # order of magnitude beyond a generous ceiling
