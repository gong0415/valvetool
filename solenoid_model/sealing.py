"""Sealing domain and the L4 leak-rate limit (spec §2.4, §3 L4).

Never imports dynamics/magnetics/fluid -- pure algebraic helpers, the same
invariant winding.py and sizing.py hold. leak_rate covers the CLOSED state
only (the contact-slit problem); once the poppet lifts off the seat the flow
is an orifice problem already covered by fluid.py, and the two are combined
at the analysis layer (limits.py), never here.
"""
import math
from dataclasses import dataclass

# EMPIRICAL: 2-D continuum percolation threshold for the real-contact area
# fraction at which the leak path is pinched off. Theory gives 0.42;
# reported values span 0.40-0.50.
PHI_PERCOLATION = 0.42

# GHe properties at 293 K, used for the internal-leak spec (scc/s GHe)
MU_HE = 1.99e-5           # dynamic viscosity, Pa*s
M_HE = 4.0026e-3          # molar mass, kg/mol
R_UNIVERSAL = 8.314       # J/(mol*K)
T_REF = 293.0             # K
LAMBDA_HE_1ATM = 1.94e-7  # mean free path at 101325 Pa, 293 K, m
P_STD = 101325.0          # Pa, standard pressure defining scc


@dataclass(frozen=True)
class SeatPair:
    """A seat material pairing. Deliberately separate from ValveParams for
    the same reason fluid.GasProperties is: one valve can be built with
    different seat materials (ARCHITECTURE.md layer-1 principle)."""
    name: str
    H: float              # indentation hardness of the softer member, Pa
    E_star: float         # reduced modulus, Pa
    Rq_c: float           # combined rms roughness sqrt(Rq1**2+Rq2**2), m
    e_restitution: float  # coefficient of restitution, - (EMPIRICAL)


def _delta_p(params, cond):
    """Seat pressure differential [Pa]. Same dual-mode rule the rest of the
    project uses: explicit conditions supersede params.delta_P, so the two
    never double-count."""
    if cond is None:
        return params.delta_P
    return cond.P_up - cond.P_down


def seal_diameter(params):
    """Seal ring equivalent diameter [m], derived from the pressure-acting
    area as A_seat = pi/4 * D_seal**2. Derived rather than an independent
    input, so no hidden degree of freedom is introduced (spec §2.0)."""
    return math.sqrt(4.0 * params.A_seat / math.pi)


def seal_force(params, cond=None):
    """Force pressing the poppet onto the seat [N] in the closed state:
    spring preload plus the pressure load, which also acts to close."""
    return params.F_preload + _delta_p(params, cond) * params.A_seat


def nominal_contact_area(params):
    """Nominal seal land area [m^2]: an annulus of width w_land at D_seal."""
    return math.pi * seal_diameter(params) * params.w_land


def contact_stress(params, cond=None):
    """Nominal contact stress on the seal land [Pa]."""
    return seal_force(params, cond) / nominal_contact_area(params)


def contact_area_fraction(params, seat, cond=None):
    """Real/nominal contact area fraction, Bowden-Tabor fully-plastic
    closure A_real = F/H. Capped at 1 (full conformity)."""
    return min(contact_stress(params, cond) / seat.H, 1.0)


def is_percolated(params, seat, cond=None):
    """True when the contact patches percolate and pinch off the leak path."""
    return contact_area_fraction(params, seat, cond) >= PHI_PERCOLATION


def residual_gap(params, seat, cond=None):
    """Mean residual interfacial gap [m], zero once percolated.

    EMPIRICAL engineering closure h = Rq_c*(1-phi). Known weakness: for
    p << H this degenerates to h -> Rq_c, i.e. the gap barely responds to
    load, so the magnitude of the leaking branch is NOT trustworthy -- only
    the threshold location is. See the design doc's finding 2."""
    if is_percolated(params, seat, cond):
        return 0.0
    return seat.Rq_c * (1.0 - contact_area_fraction(params, seat, cond))


def land_width_for_seal(params, seat, cond=None):
    """Largest land width [m] that still reaches the percolation threshold.
    The practical L4 design answer: how narrow the land must be made."""
    p_needed = PHI_PERCOLATION * seat.H
    return seal_force(params, cond) / (p_needed * math.pi * seal_diameter(params))


def mean_free_path(p_mean):
    """GHe mean free path [m] at p_mean, 293 K (ideal-gas 1/p scaling)."""
    return LAMBDA_HE_1ATM * P_STD / p_mean


def knudsen_number(gap, p_mean):
    """Kn = lambda/gap. <0.01 viscous, >10 free-molecular, between is the
    transitional regime this module blends across."""
    return mean_free_path(p_mean) / gap


def _mean_thermal_speed():
    return math.sqrt(8.0 * R_UNIVERSAL * T_REF / (math.pi * M_HE))


def leak_viscous(gap, b, w, P_up, P_down):
    """Compressible Poiseuille slit flow [scc/s]: width b, height gap,
    flow-path length w. Scales as gap**3."""
    Q = b * gap**3 * (P_up**2 - P_down**2) / (24.0 * MU_HE * w * P_STD)
    return Q * 1e6


def leak_molecular(gap, b, w, P_up, P_down):
    """Free-molecular (Knudsen) slit conductance [scc/s]. Scales as gap**2."""
    C = (2.0 / 3.0) * _mean_thermal_speed() * b * gap**2 / w
    return C * (P_up - P_down) / P_STD * 1e6


def leak_rate(params, seat, cond):
    """Internal helium leak rate [scc/s GHe] in the CLOSED state.

    Zero once the contact percolates. Otherwise the residual gap is treated
    as an annular slit of circumference pi*D_seal and flow length w_land,
    with the two flow regimes blended by f = 1/(1+Kn) -- an EMPIRICAL
    transitional-regime treatment, not a first-principles result.

    Only the closed state: a lifted poppet is an orifice problem belonging
    to fluid.py, combined with this at the analysis layer.

    NOTE: per the design doc's finding 2, the magnitude returned on the
    leaking branch is not trustworthy (the gap closure saturates at Rq_c);
    the trustworthy output is whether the seat percolates at all, and
    land_width_for_seal.
    """
    gap = residual_gap(params, seat, cond)
    if gap <= 0.0:
        return 0.0
    b = math.pi * seal_diameter(params)
    w = params.w_land
    Kn = knudsen_number(gap, 0.5 * (cond.P_up + cond.P_down))
    f = 1.0 / (1.0 + Kn)
    return (f * leak_viscous(gap, b, w, cond.P_up, cond.P_down)
            + (1.0 - f) * leak_molecular(gap, b, w, cond.P_up, cond.P_down))
