import math

import pytest

from solenoid_model import materials, sealing


def test_every_material_declares_a_source_tag():
    assert materials.MATERIALS
    for name, mat in materials.MATERIALS.items():
        assert mat.source in ("ASTM_E595", "LITERATURE"), name


def test_astm_e595_thresholds_are_the_spec_values():
    assert materials.TML_LIMIT == 1.0
    assert materials.CVCM_LIMIT == 0.1


def test_outgassing_screen_passes_pctfe_and_fails_vespel_sp1():
    assert materials.outgassing_ok(materials.PCTFE)
    # Vespel SP-1 absorbs moisture; its TML exceeds the 1% screen
    assert not materials.outgassing_ok(materials.VESPEL_SP1)


def test_seat_pair_takes_the_softer_hardness_and_composes_the_modulus():
    pair = materials.seat_pair(materials.PCTFE, materials.SS_440C, 0.4)
    assert isinstance(pair, sealing.SeatPair)
    assert math.isclose(pair.H, materials.PCTFE.H, rel_tol=1e-12)
    expected = 1.0 / ((1 - materials.PCTFE.nu**2) / materials.PCTFE.E
                      + (1 - materials.SS_440C.nu**2) / materials.SS_440C.E)
    assert math.isclose(pair.E_star, expected, rel_tol=1e-12)
    assert math.isclose(pair.E_star, 1.696042e9, rel_tol=1e-5)
    assert math.isclose(pair.Rq_c, 4.123106e-7, rel_tol=1e-5)


def test_preset_pairs_match_the_values_used_across_p6():
    assert math.isclose(materials.PCTFE_ON_440C.E_star, 1.696042e9, rel_tol=1e-5)
    assert math.isclose(materials.PCTFE_ON_440C.H, 1.0e8, rel_tol=1e-12)
    assert math.isclose(
        materials.METAL_17_4PH_ON_440C.E_star, 1.073642e11, rel_tol=1e-5)
    assert math.isclose(materials.METAL_17_4PH_ON_440C.H, 4.0e9, rel_tol=1e-12)


def test_cold_weld_risk_is_high_for_like_metals_under_high_stress():
    assert materials.cold_weld_risk(
        materials.SS_440C, materials.SS_440C, 3.0e9) == "high"


def test_cold_weld_risk_is_low_whenever_a_polymer_is_involved():
    assert materials.cold_weld_risk(
        materials.PCTFE, materials.SS_440C, 3.0e9) == "low"


def test_cold_weld_risk_is_moderate_for_unlike_metals():
    assert materials.cold_weld_risk(
        materials.SS_17_4PH, materials.SS_440C, 3.0e9) == "moderate"


def test_radiation_margin_is_limit_over_dose():
    m = materials.radiation_margin(materials.PTFE, 1.0e3)
    assert math.isclose(m, materials.PTFE.dose_limit / 1.0e3, rel_tol=1e-12)
    assert m > 1.0


def test_ptfe_fails_radiation_where_vespel_passes():
    dose = 1.0e6
    assert materials.radiation_margin(materials.PTFE, dose) < 1.0
    assert materials.radiation_margin(materials.VESPEL_SP1, dose) > 1.0


def test_propellant_compatibility_returns_a_grade():
    assert materials.propellant_compatibility(materials.PCTFE, "hydrazine") == "A"


def test_unknown_propellant_raises():
    with pytest.raises(KeyError):
        materials.propellant_compatibility(materials.PCTFE, "unobtainium")
