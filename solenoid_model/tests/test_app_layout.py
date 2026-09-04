import pytest

from solenoid_model import app, layout
from solenoid_model.cases import N2_25BAR_PH_PARAMS


def test_render_layout_tab_exists():
    assert hasattr(app, "render_layout_tab")


def test_compute_layout_returns_full_chain():
    result = app.compute_layout(N2_25BAR_PH_PARAMS)
    assert result["envelope"]["OD"] == pytest.approx(20.0e-3, abs=1e-5)
    assert result["envelope"]["L"] == pytest.approx(18.679e-3, abs=2e-5)
    assert result["wall"].adopted == pytest.approx(1.0e-3, rel=1e-12)
    assert result["spring"].n_active == pytest.approx(6.35, abs=0.02)


def test_compute_layout_reports_armature_consistency():
    result = app.compute_layout(N2_25BAR_PH_PARAMS)
    assert result["armature_rel_error"] < 1e-9


def test_compute_layout_surfaces_empirical_entries():
    result = app.compute_layout(N2_25BAR_PH_PARAMS)
    assert result["empirical"] is layout.LAYOUT_EMPIRICAL
