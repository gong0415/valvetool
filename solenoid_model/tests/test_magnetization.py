import math

import pytest

from solenoid_model import magnetization
from solenoid_model.cases import N2_25BAR_PH_PARAMS


def test_analytic_low_H_reduces_to_linear():
    """H->0 時 Frohlich-Kennelly 必須退化為 B = mu_r*mu_0*H。"""
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    H = 0.1
    linear = 4000.0 * magnetization.MU_0 * H
    assert math.isclose(bh.B_of_H(H), linear, rel_tol=1e-3)


def test_analytic_high_H_approaches_B_sat_without_exceeding():
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    assert bh.B_of_H(1e6) == pytest.approx(2.1491, abs=1e-3)
    for H in (1e3, 1e6, 1e9):
        assert bh.B_of_H(H) < 2.15


def test_analytic_is_monotonic():
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    prev = -1.0
    for k in range(0, 9):
        B = bh.B_of_H(10.0 ** k)
        assert B > prev
        prev = B


def test_analytic_round_trip_H_of_B():
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    for H in (1.0, 1e2, 1e4):
        assert math.isclose(bh.H_of_B(bh.B_of_H(H)), H, rel_tol=1e-9)


def test_H_of_B_at_or_above_B_sat_is_infinite():
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    assert math.isinf(bh.H_of_B(2.15))
    assert math.isinf(bh.H_of_B(3.0))


def test_from_params_uses_existing_fields_only():
    bh = magnetization.from_params(N2_25BAR_PH_PARAMS)
    assert bh.mu_r == N2_25BAR_PH_PARAMS.mu_r_core
    assert bh.B_sat == N2_25BAR_PH_PARAMS.B_sat


def test_tabulated_matches_analytic_on_its_own_points():
    """同一組點上，查表與解析式必須一致（介面可互換）。"""
    analytic = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    Hs = (1.0, 10.0, 100.0, 1000.0, 10000.0)
    table = magnetization.TabulatedBH(
        tuple((H, analytic.B_of_H(H)) for H in Hs))
    for H in Hs:
        assert math.isclose(table.B_of_H(H), analytic.B_of_H(H), rel_tol=1e-12)


def test_tabulated_interpolates_between_points():
    table = magnetization.TabulatedBH(((0.0, 0.0), (100.0, 1.0)))
    assert math.isclose(table.B_of_H(50.0), 0.5, rel_tol=1e-12)


def test_tabulated_clamps_above_last_point():
    """超出表格上界時鉗位在最後一點，不外插（外插會給出非物理的 B）。"""
    table = magnetization.TabulatedBH(((0.0, 0.0), (100.0, 1.0)))
    assert math.isclose(table.B_of_H(1e6), 1.0, rel_tol=1e-12)


def test_tabulated_rejects_non_ascending_points():
    with pytest.raises(ValueError, match="ascending"):
        magnetization.TabulatedBH(((0.0, 0.0), (100.0, 1.0), (50.0, 0.5)))


def test_tabulated_rejects_too_few_points():
    with pytest.raises(ValueError, match="at least 2"):
        magnetization.TabulatedBH(((0.0, 0.0),))


def test_analytic_rejects_nonpositive_B_sat():
    with pytest.raises(ValueError, match="B_sat"):
        magnetization.AnalyticBH(mu_r=4000.0, B_sat=0.0)
