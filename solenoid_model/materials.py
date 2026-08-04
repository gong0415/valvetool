"""Materials and space-environment domain (spec §2.6).

Mostly a lookup table, not a dynamics model -- and deliberately a small one,
covering only what the baseline valve is built from.

SOURCE DISCIPLINE. Every entry carries a `source` tag with exactly two
possible values:

  ASTM_E595   the TML < 1% / CVCM < 0.1% screening thresholds themselves.
              These are specification limits, i.e. hard numbers.
  LITERATURE  everything else -- per-material TML/CVCM figures, radiation
              dose limits, hardness, compatibility grades. These are
              textbook/handbook-level values recorded to give the model
              something to run on. They are NOT measurements, and before any
              flight use they must be re-verified against the NASA
              outgassing database and the relevant material certificates.

Imports sealing (for SeatPair); sealing must never import this module.
"""
import math
from dataclasses import dataclass

from solenoid_model import sealing

# ASTM E595 screening thresholds -- specification values, not estimates.
TML_LIMIT = 1.0    # total mass loss, %
CVCM_LIMIT = 0.1   # collected volatile condensable material, %

# EMPIRICAL: contact stress above this fraction of the softer member's
# hardness is treated as intimate enough for like-metal cold welding.
COLD_WELD_STRESS_FRACTION = 0.5


@dataclass(frozen=True)
class MaterialProperties:
    """Immutable material record. See the module docstring on `source`."""
    name: str
    E: float                    # Young's modulus, Pa
    nu: float                   # Poisson ratio, -
    H: float                    # indentation hardness, Pa
    Rq_typical: float           # typical finished roughness, m
    TML: float                  # total mass loss, %
    CVCM: float                 # condensable volatiles, %
    dose_limit: float           # radiation dose limit, rad(Si)
    compatibility: dict         # propellant -> grade
    source: str
    is_polymer: bool


# All per-material numbers below are LITERATURE grade -- see module docstring.
PCTFE = MaterialProperties(
    "PCTFE", E=1.5e9, nu=0.35, H=1.0e8, Rq_typical=0.40e-6,
    TML=0.02, CVCM=0.00, dose_limit=1.0e6,
    compatibility={"hydrazine": "A", "MMH": "A", "NTO": "B", "H2O2": "B",
                   "GN2": "A", "GHe": "A"},
    source="LITERATURE", is_polymer=True)

VESPEL_SP1 = MaterialProperties(
    "Vespel SP-1", E=3.1e9, nu=0.41, H=2.5e8, Rq_typical=0.50e-6,
    TML=1.09, CVCM=0.00, dose_limit=1.0e9,
    compatibility={"hydrazine": "B", "MMH": "B", "NTO": "C", "H2O2": "C",
                   "GN2": "A", "GHe": "A"},
    source="LITERATURE", is_polymer=True)

PTFE = MaterialProperties(
    "PTFE", E=0.5e9, nu=0.46, H=3.0e7, Rq_typical=0.60e-6,
    TML=0.03, CVCM=0.00, dose_limit=1.0e4,
    compatibility={"hydrazine": "A", "MMH": "A", "NTO": "A", "H2O2": "A",
                   "GN2": "A", "GHe": "A"},
    source="LITERATURE", is_polymer=True)

SS_17_4PH = MaterialProperties(
    "17-4PH H900", E=197e9, nu=0.27, H=4.0e9, Rq_typical=0.10e-6,
    TML=0.0, CVCM=0.0, dose_limit=1.0e12,
    compatibility={"hydrazine": "A", "MMH": "A", "NTO": "B", "H2O2": "B",
                   "GN2": "A", "GHe": "A"},
    source="LITERATURE", is_polymer=False)

SS_440C = MaterialProperties(
    "440C", E=200e9, nu=0.28, H=6.0e9, Rq_typical=0.10e-6,
    TML=0.0, CVCM=0.0, dose_limit=1.0e12,
    compatibility={"hydrazine": "B", "MMH": "B", "NTO": "C", "H2O2": "C",
                   "GN2": "A", "GHe": "A"},
    source="LITERATURE", is_polymer=False)

SS_430F = MaterialProperties(
    "430F", E=200e9, nu=0.28, H=2.0e9, Rq_typical=0.20e-6,
    TML=0.0, CVCM=0.0, dose_limit=1.0e12,
    compatibility={"hydrazine": "B", "MMH": "B", "NTO": "C", "H2O2": "C",
                   "GN2": "A", "GHe": "A"},
    source="LITERATURE", is_polymer=False)

PURE_IRON = MaterialProperties(
    "ARMCO pure iron", E=200e9, nu=0.29, H=1.0e9, Rq_typical=0.20e-6,
    TML=0.0, CVCM=0.0, dose_limit=1.0e12,
    compatibility={"hydrazine": "C", "MMH": "C", "NTO": "D", "H2O2": "D",
                   "GN2": "A", "GHe": "A"},
    source="LITERATURE", is_polymer=False)

COPPER = MaterialProperties(
    "Cu (magnet wire)", E=117e9, nu=0.34, H=0.87e9, Rq_typical=0.40e-6,
    TML=0.0, CVCM=0.0, dose_limit=1.0e12,
    compatibility={"hydrazine": "D", "MMH": "D", "NTO": "D", "H2O2": "D",
                   "GN2": "A", "GHe": "A"},
    source="LITERATURE", is_polymer=False)

MATERIALS = {m.name: m for m in (
    PCTFE, VESPEL_SP1, PTFE, SS_17_4PH, SS_440C, SS_430F, PURE_IRON, COPPER)}


def seat_pair(mat_a, mat_b, e_restitution):
    """Compose a sealing.SeatPair from two materials.

    Hardness is the softer member's (it is what yields); the reduced modulus
    follows 1/E* = (1-nu_a**2)/E_a + (1-nu_b**2)/E_b; roughness combines in
    quadrature. e_restitution is an EMPIRICAL input, not derived here."""
    E_star = 1.0 / ((1 - mat_a.nu**2) / mat_a.E + (1 - mat_b.nu**2) / mat_b.E)
    Rq_c = math.sqrt(mat_a.Rq_typical**2 + mat_b.Rq_typical**2)
    return sealing.SeatPair(
        name=f"{mat_a.name}/{mat_b.name}",
        H=min(mat_a.H, mat_b.H), E_star=E_star, Rq_c=Rq_c,
        e_restitution=e_restitution)


# EMPIRICAL restitution coefficients: polymer seats damp far more than steel.
PCTFE_ON_440C = seat_pair(PCTFE, SS_440C, 0.4)
METAL_17_4PH_ON_440C = seat_pair(SS_17_4PH, SS_440C, 0.6)


def outgassing_ok(mat):
    """ASTM E595 screen: TML < 1% and CVCM < 0.1%."""
    return mat.TML < TML_LIMIT and mat.CVCM < CVCM_LIMIT


def cold_weld_risk(mat_a, mat_b, p_contact):
    """Criteria-based cold-welding risk: "low" | "moderate" | "high".

    Not a derived model. Polymers do not cold weld, so any polymer member
    gives "low". Like metals in vacuum pressed to an appreciable fraction of
    their hardness are the classic risk case; unlike metals, which carry
    dissimilar oxide films, are treated as intermediate."""
    if mat_a.is_polymer or mat_b.is_polymer:
        return "low"
    H_soft = min(mat_a.H, mat_b.H)
    if mat_a.name == mat_b.name and p_contact >= COLD_WELD_STRESS_FRACTION * H_soft:
        return "high"
    return "moderate"


def radiation_margin(mat, dose):
    """Dose margin = dose_limit/dose. Greater than 1 passes."""
    return mat.dose_limit / dose


def propellant_compatibility(mat, propellant):
    """Compatibility grade: "A" compatible, "B" limited, "C" not recommended,
    "D" incompatible. Raises KeyError for a propellant not in the table --
    silently returning a default would be the dangerous behaviour here."""
    return mat.compatibility[propellant]
