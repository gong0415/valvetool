import math

from solenoid_model import impact, sealing
from solenoid_model.baseline import BASELINE_PARAMS as P

PCTFE_ON_440C = sealing.SeatPair(
    name="PCTFE/440C", H=1.0e8, E_star=1.696042e9,
    Rq_c=4.123106e-7, e_restitution=0.4)
METAL_17_4_ON_440C = sealing.SeatPair(
    name="17-4PH/440C", H=4.0e9, E_star=1.073642e11,
    Rq_c=1.414214e-7, e_restitution=0.6)

M = 1e-3        # baseline armature mass, kg
R = 1.0e-3      # baseline poppet tip radius, m
V_DIODE = 0.3934    # measured seat-impact speed, diode flyback, m/s
V_ZENER = 1.5338    # measured seat-impact speed, zener flyback, m/s
V_STOP = 0.8543     # measured stop-impact speed on opening, m/s


def test_hertz_stiffness_matches_closed_form():
    k = impact.hertz_stiffness(PCTFE_ON_440C, R)
    assert math.isclose(
        k, (4.0 / 3.0) * PCTFE_ON_440C.E_star * math.sqrt(R), rel_tol=1e-12)


def test_indentation_max_satisfies_the_energy_balance():
    # 1/2 m v^2 must equal the Hertz strain energy (2/5) k d^2.5
    d = impact.indentation_max(M, V_DIODE, PCTFE_ON_440C, R)
    k = impact.hertz_stiffness(PCTFE_ON_440C, R)
    assert math.isclose(0.5 * M * V_DIODE**2, 0.4 * k * d**2.5, rel_tol=1e-12)


def test_contact_force_max_is_stiffness_times_indentation_to_1p5():
    d = impact.indentation_max(M, V_ZENER, METAL_17_4_ON_440C, R)
    k = impact.hertz_stiffness(METAL_17_4_ON_440C, R)
    F = impact.contact_force_max(M, V_ZENER, METAL_17_4_ON_440C, R)
    assert math.isclose(F, k * d**1.5, rel_tol=1e-12)


def test_contact_pressure_matches_prototyped_values():
    p = impact.contact_pressure_max(M, V_DIODE, PCTFE_ON_440C, R)
    assert math.isclose(p, 165.866e6, rel_tol=1e-4)
    p = impact.contact_pressure_max(M, V_ZENER, PCTFE_ON_440C, R)
    assert math.isclose(p, 285.845e6, rel_tol=1e-4)
    p = impact.contact_pressure_max(M, V_DIODE, METAL_17_4_ON_440C, R)
    assert math.isclose(p, 4580.319e6, rel_tol=1e-4)


def test_contact_pressure_scales_as_velocity_to_the_two_fifths():
    a = impact.contact_pressure_max(M, 1.0, PCTFE_ON_440C, R)
    b = impact.contact_pressure_max(M, 32.0, PCTFE_ON_440C, R)
    assert math.isclose(b / a, 32.0**0.4, rel_tol=1e-9)


def test_contact_duration_matches_prototyped_values():
    assert math.isclose(
        impact.contact_duration(M, V_DIODE, PCTFE_ON_440C, R), 157.458e-6,
        rel_tol=1e-4)
    assert math.isclose(
        impact.contact_duration(M, V_ZENER, METAL_17_4_ON_440C, R), 22.825e-6,
        rel_tol=1e-4)


def test_every_baseline_impact_is_plastic_on_both_seat_types():
    # design doc finding 4: p_max/H > 1 in all six combinations, so each
    # cycle indents the seat -- this is the L5 wear mechanism itself
    for seat in (PCTFE_ON_440C, METAL_17_4_ON_440C):
        for v in (V_DIODE, V_ZENER, V_STOP):
            assert not impact.is_shakedown(M, v, seat, R)


def test_shakedown_is_reached_at_sufficiently_low_speed():
    assert impact.is_shakedown(M, 1e-3, METAL_17_4_ON_440C, R)


def test_bounce_sequence_is_geometric():
    speeds = impact.bounce_sequence(1.0, 0.5, v_min=0.05, n_max=50)
    assert speeds[0] == 1.0
    for a, b in zip(speeds, speeds[1:]):
        assert math.isclose(b / a, 0.5, rel_tol=1e-12)


def test_bounce_sequence_counts_match_prototyped_values():
    # measured at the opening stop speed with v_min = 1 cm/s
    assert len(impact.bounce_sequence(V_STOP, 0.2)) == 3
    assert len(impact.bounce_sequence(V_STOP, 0.5)) == 7
    assert len(impact.bounce_sequence(V_STOP, 0.8)) == 20


def test_bounce_sequence_stops_when_rebound_falls_below_v_min():
    # e*v0 <= v_min -> the first impact is the only one
    assert impact.bounce_sequence(V_STOP, 1e-6) == [V_STOP]


def test_bounce_sequence_honours_the_zeno_guard():
    assert len(impact.bounce_sequence(1.0, 0.999999, n_max=12)) == 12


def test_impact_energy_total_approaches_the_initial_kinetic_energy():
    E0 = 0.5 * M * V_STOP**2
    for e in (0.2, 0.5, 0.8):
        E = impact.impact_energy_total(M, V_STOP, e)
        assert E <= E0
        assert math.isclose(E, E0, rel_tol=2e-3)


def test_contact_time_stays_well_below_the_stroke_time():
    # Guard for the instantaneous-restitution assumption (design doc
    # finding 4). Measured stroke phases: 0.6181 ms opening motion,
    # 1.3456 ms diode closing stroke, 0.4912 ms zener closing stroke.
    # Measured ratios span 4.10x (PCTFE + zener) to 44.91x (metal + diode);
    # the threshold is 3x, so drifting parameters trip this before the
    # idealisation silently stops being valid.
    cases = ((V_STOP, 0.6181e-3), (V_DIODE, 1.3456e-3), (V_ZENER, 0.4912e-3))
    for seat in (PCTFE_ON_440C, METAL_17_4_ON_440C):
        for v, t_stroke in cases:
            t_c = impact.contact_duration(M, v, seat, R)
            assert t_stroke / t_c >= 3.0, (seat.name, v, t_stroke / t_c)
