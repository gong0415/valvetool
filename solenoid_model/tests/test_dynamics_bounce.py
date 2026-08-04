import math

from solenoid_model import dynamics, impact, materials
from solenoid_model.baseline import BASELINE_PARAMS as P

SEAT = materials.PCTFE_ON_440C


def test_seat_none_is_bit_identical_to_the_unbounced_solution():
    a = dynamics.simulate_opening(P)
    b = dynamics.simulate_opening(P, seat=None)
    assert a.t_events[0][0] == b.t_events[0][0]
    assert list(a.t) == list(b.t)
    for ra, rb in zip(a.y, b.y):
        assert list(ra) == list(rb)


def test_bounce_mode_returns_a_bounce_result_with_the_same_first_contact():
    plain = dynamics.simulate_opening(P)
    r = dynamics.simulate_opening(P, seat=SEAT)
    assert isinstance(r, dynamics.BounceResult)
    assert math.isclose(r.t_first, plain.t_events[0][0], rel_tol=1e-12)


def test_bouncing_delays_settling_beyond_first_contact():
    r = dynamics.simulate_opening(P, seat=SEAT)
    assert r.n_bounce >= 2
    assert r.t_settle > r.t_first
    assert len(r.v_impacts) == r.n_bounce
    assert len(r.segments) == r.n_bounce


def test_impact_speeds_converge_to_the_restitution_coefficient_from_above():
    # The ODE bounce train does NOT decay by exactly e, and that is physics,
    # not a modelling error. The armature reaches the travel stop while the
    # coil current is still only ~11% of steady state (0.0387 A of 0.35 A),
    # so the magnetic force there (~28.9 N) barely exceeds the spring plus
    # pressure load (~28 N). With that little margin the first rebound flies
    # far and long (~264 us). At fixed current the magnetic force derives
    # from the co-energy potential and does ZERO net work over a closed
    # excursion, so the sole energy source available during that flight is
    # the coil current's continued rise; position-dependence does not add
    # energy on its own, it only amplifies the effect by widening the range
    # of x over which the still-rising force is sampled. That is what
    # returns the armature FASTER than e*v. Later flights shrink
    # geometrically (264->77->28->11->4 us), narrowing that range until the
    # effect dies out, so the ratio converges to e from above.
    # Measured: PCTFE (e=0.4) 0.605, 0.430, 0.409, 0.403, 0.401;
    #           metal (e=0.6) 0.898, 0.665, 0.628, 0.614, ..., 0.601.
    for seat in (materials.PCTFE_ON_440C, materials.METAL_17_4PH_ON_440C):
        r = dynamics.simulate_opening(P, seat=seat)
        e = seat.e_restitution
        ratios = [b / a for a, b in zip(r.v_impacts, r.v_impacts[1:])]
        assert len(ratios) >= 3, seat.name
        assert all(b < a for a, b in zip(r.v_impacts, r.v_impacts[1:])), seat.name
        assert all(x >= e for x in ratios), (seat.name, ratios)
        assert ratios[0] > 1.2 * e, (seat.name, ratios[0])
        assert math.isclose(ratios[-1], e, rel_tol=0.03), (seat.name, ratios[-1])


def test_dead_seat_settles_immediately():
    from dataclasses import replace
    dead = replace(SEAT, e_restitution=1e-6)
    plain = dynamics.simulate_opening(P)
    r = dynamics.simulate_opening(P, seat=dead)
    assert r.n_bounce == 1
    assert r.t_settle == r.t_first
    assert math.isclose(r.t_settle, plain.t_events[0][0], rel_tol=1e-12)


def test_zeno_guard_caps_the_impact_count():
    from dataclasses import replace
    springy = replace(SEAT, e_restitution=0.999)
    r = dynamics.simulate_opening(P, seat=springy, n_max=4)
    assert r.n_bounce <= 4


def test_ode_bounce_count_agrees_with_the_analytic_train():
    # The ODE train runs systematically LONGER than the geometric-decay
    # prediction, never shorter: impact speeds converge to e from ABOVE (see
    # test_impact_speeds_converge_to_the_restitution_coefficient_from_above),
    # so each real impact is a bit faster than e**k * v0 would predict, which
    # delays the train crossing v_min by one extra impact. Measured: PCTFE
    # 6 vs analytic 5, metal 10 vs analytic 9 -- both exactly +1, so the old
    # <=1 bound passed with zero margin on both seats. Widened to <=2 so the
    # next parameter tweak has room before this becomes a tripwire.
    for seat in (materials.PCTFE_ON_440C, materials.METAL_17_4PH_ON_440C):
        r = dynamics.simulate_opening(P, seat=seat)
        analytic = impact.bounce_sequence(r.v_impacts[0], seat.e_restitution)
        assert abs(r.n_bounce - len(analytic)) <= 2, (
            seat.name, r.n_bounce, len(analytic))


from solenoid_model import flyback


def test_closing_without_seat_reports_no_bounce():
    r = dynamics.simulate_closing(P, flyback.DiodeFlyback())
    assert r.bounce is None


def test_closing_t_close_is_unchanged_by_enabling_bounce():
    plain = dynamics.simulate_closing(P, flyback.DiodeFlyback())
    bounced = dynamics.simulate_closing(P, flyback.DiodeFlyback(), seat=SEAT)
    assert math.isclose(bounced.t_close, plain.t_close, rel_tol=1e-12)
    assert math.isclose(bounced.t_release, plain.t_release, rel_tol=1e-12)


def test_closing_bounce_settles_after_first_contact():
    r = dynamics.simulate_closing(P, flyback.DiodeFlyback(), seat=SEAT)
    assert r.bounce.n_bounce >= 2
    assert r.bounce.t_settle > r.bounce.t_first


def test_zener_slams_the_seat_harder_than_the_diode():
    # the P6 headline: P3's zener speedup is paid for in impact energy
    d = dynamics.simulate_closing(P, flyback.DiodeFlyback(), seat=SEAT)
    z = dynamics.simulate_closing(P, flyback.ZenerFlyback(), seat=SEAT)
    assert z.bounce.v_impacts[0] > 3.0 * d.bounce.v_impacts[0]
    assert z.t_close < 0.5 * d.t_close
