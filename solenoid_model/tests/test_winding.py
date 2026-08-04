import math

from solenoid_model import winding

# calibrated winding geometry (see docs/spec/PARAMS.md, P4);
# chosen so coil_resistance at N=2000 equals the baseline R_coil_20C=80 ohm
A_WINDING = 5.2190904529497215e-05
K_FILL = 0.5
L_TURN_MEAN = 0.031066014600891194


def test_coil_resistance_calibrated_to_baseline():
    R = winding.coil_resistance(2000.0, A_WINDING, K_FILL, L_TURN_MEAN)
    assert math.isclose(R, 80.0, rel_tol=1e-12)


def test_coil_resistance_scales_with_n_squared_over_window():
    # R = rho*N^2*l_turn/(A_w*k_fill): fixed window -> N^2/R invariant
    ratios = []
    for N in (1000.0, 2000.0, 4000.0, 8000.0):
        R = winding.coil_resistance(N, A_WINDING, K_FILL, L_TURN_MEAN)
        ratios.append(N**2 / R)
    assert all(math.isclose(r, 50000.0, rel_tol=1e-9) for r in ratios)


def test_wire_diameter_matches_closed_form():
    N = 2000.0
    d = winding.wire_diameter(N, A_WINDING, K_FILL)
    a_wire = A_WINDING * K_FILL / N
    expected = math.sqrt(4.0 * a_wire / math.pi)
    assert math.isclose(d, expected, rel_tol=1e-12)
    assert math.isclose(d, 0.0001288909650852741, rel_tol=1e-9)


def test_turns_for_resistance_inverts_coil_resistance():
    N = winding.turns_for_resistance(80.0, A_WINDING, K_FILL, L_TURN_MEAN)
    assert math.isclose(N, 2000.0, rel_tol=1e-9)
    # round-trip: feeding N back yields the same R
    R = winding.coil_resistance(N, A_WINDING, K_FILL, L_TURN_MEAN)
    assert math.isclose(R, 80.0, rel_tol=1e-9)


def test_coil_resistance_consistency_is_exact_at_baseline():
    from solenoid_model.baseline import BASELINE_PARAMS
    R_declared, R_geometric, rel_error = winding.coil_resistance_consistency(
        BASELINE_PARAMS)
    assert math.isclose(R_declared, 80.0, rel_tol=1e-12)
    assert math.isclose(R_geometric, 80.0, rel_tol=1e-12)
    assert math.isclose(rel_error, 0.0, abs_tol=1e-12)


def test_coil_resistance_consistency_detects_stale_declared_resistance():
    # Changing N_turns alone leaves params.R_coil_20C stale: the geometry says
    # 320 ohm at N=4000 while the declared value is still 80 ohm. dynamics
    # would silently simulate a coil that cannot exist.
    from dataclasses import replace
    from solenoid_model.baseline import BASELINE_PARAMS
    p = replace(BASELINE_PARAMS, N_turns=4000.0)
    R_declared, R_geometric, rel_error = winding.coil_resistance_consistency(p)
    assert math.isclose(R_declared, 80.0, rel_tol=1e-12)
    assert math.isclose(R_geometric, 320.0, rel_tol=1e-9)
    assert math.isclose(rel_error, 3.0, rel_tol=1e-9)
