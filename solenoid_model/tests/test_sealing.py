import math

import pytest

from solenoid_model import sealing
from solenoid_model.baseline import BASELINE_PARAMS as P

PCTFE_ON_440C = sealing.SeatPair(
    name="PCTFE/440C", H=1.0e8, E_star=1.696042e9,
    Rq_c=4.123106e-7, e_restitution=0.4)
METAL_17_4_ON_440C = sealing.SeatPair(
    name="17-4PH/440C", H=4.0e9, E_star=1.073642e11,
    Rq_c=1.414214e-7, e_restitution=0.6)


def test_seal_diameter_derives_from_seat_area():
    D = sealing.seal_diameter(P)
    assert math.isclose(D, math.sqrt(4.0 * P.A_seat / math.pi), rel_tol=1e-12)
    assert math.isclose(D, 2.523132522e-3, rel_tol=1e-9)


def test_seal_force_is_preload_plus_pressure_load():
    F = sealing.seal_force(P)
    assert math.isclose(F, 10.0 + 12.0, rel_tol=1e-12)
    assert math.isclose(F, 22.0, rel_tol=1e-12)


def test_conditions_supersede_delta_p_when_given():
    class Cond:
        P_up, P_down = 1.2e6, 0.2e6
    F = sealing.seal_force(P, cond=Cond())
    assert math.isclose(F, 10.0 + 1.0e6 * P.A_seat, rel_tol=1e-12)


def test_contact_stress_at_baseline_land_width():
    p = sealing.contact_stress(P)
    assert math.isclose(p, 55.5089e6, rel_tol=1e-5)


def test_contact_area_fraction_is_capped_at_unity():
    from dataclasses import replace
    tiny_land = replace(P, w_land=1e-7)  # absurdly narrow -> p >> H
    assert sealing.contact_area_fraction(tiny_land, PCTFE_ON_440C) == 1.0


def test_percolation_separates_the_two_seat_types_at_baseline():
    assert math.isclose(
        sealing.contact_area_fraction(P, PCTFE_ON_440C), 0.555089, rel_tol=1e-5)
    assert math.isclose(
        sealing.contact_area_fraction(P, METAL_17_4_ON_440C), 0.013877, rel_tol=1e-4)
    assert sealing.is_percolated(P, PCTFE_ON_440C)
    assert not sealing.is_percolated(P, METAL_17_4_ON_440C)


def test_residual_gap_is_zero_once_percolated():
    assert sealing.residual_gap(P, PCTFE_ON_440C) == 0.0


def test_residual_gap_below_threshold_matches_closed_form():
    phi = sealing.contact_area_fraction(P, METAL_17_4_ON_440C)
    h = sealing.residual_gap(P, METAL_17_4_ON_440C)
    assert math.isclose(h, METAL_17_4_ON_440C.Rq_c * (1.0 - phi), rel_tol=1e-12)


def test_land_width_for_seal_matches_prototyped_values():
    assert math.isclose(
        sealing.land_width_for_seal(P, PCTFE_ON_440C), 66.0820e-6, rel_tol=1e-5)
    assert math.isclose(
        sealing.land_width_for_seal(P, METAL_17_4_ON_440C), 1.6521e-6, rel_tol=1e-4)


def test_land_width_for_seal_inverts_is_percolated():
    # at exactly the returned land width the pair is (marginally) percolated;
    # 1% wider and it is not
    from dataclasses import replace
    for seat in (PCTFE_ON_440C, METAL_17_4_ON_440C):
        w = sealing.land_width_for_seal(P, seat)
        assert sealing.is_percolated(replace(P, w_land=w * 0.99), seat)
        assert not sealing.is_percolated(replace(P, w_land=w * 1.01), seat)


class _Cond:
    def __init__(self, P_up, P_down):
        self.P_up, self.P_down = P_up, P_down


BASELINE_COND = _Cond(2.4e6, 0.0)


def test_mean_free_path_scales_inversely_with_pressure():
    assert math.isclose(sealing.mean_free_path(101325.0), 1.94e-7, rel_tol=1e-12)
    assert math.isclose(sealing.mean_free_path(2 * 101325.0), 0.97e-7, rel_tol=1e-12)


def test_knudsen_number_is_mfp_over_gap():
    assert math.isclose(
        sealing.knudsen_number(1e-7, 101325.0), 1.94, rel_tol=1e-12)


def test_viscous_branch_scales_as_gap_cubed():
    a = sealing.leak_viscous(1e-7, 1e-3, 1e-4, 2.4e6, 0.0)
    b = sealing.leak_viscous(2e-7, 1e-3, 1e-4, 2.4e6, 0.0)
    assert math.isclose(b / a, 8.0, rel_tol=1e-12)


def test_molecular_branch_scales_as_gap_squared():
    a = sealing.leak_molecular(1e-7, 1e-3, 1e-4, 2.4e6, 0.0)
    b = sealing.leak_molecular(2e-7, 1e-3, 1e-4, 2.4e6, 0.0)
    assert math.isclose(b / a, 4.0, rel_tol=1e-12)


def test_leak_rate_is_zero_when_percolated():
    assert sealing.leak_rate(P, PCTFE_ON_440C, BASELINE_COND) == 0.0


def test_leak_rate_below_threshold_is_orders_above_the_spec_line():
    # design doc finding 1: sealing is near-binary -- an unconformed seat
    # misses the 1e-4 scc/s internal-leak spec by several orders
    q = sealing.leak_rate(P, METAL_17_4_ON_440C, BASELINE_COND)
    assert q > 1e-3
    assert q / 1e-4 > 10.0


def test_leak_rate_lies_between_the_two_pure_branches():
    h = sealing.residual_gap(P, METAL_17_4_ON_440C, BASELINE_COND)
    b = math.pi * sealing.seal_diameter(P)
    q_v = sealing.leak_viscous(h, b, P.w_land, 2.4e6, 0.0)
    q_m = sealing.leak_molecular(h, b, P.w_land, 2.4e6, 0.0)
    q = sealing.leak_rate(P, METAL_17_4_ON_440C, BASELINE_COND)
    assert min(q_v, q_m) <= q <= max(q_v, q_m)
