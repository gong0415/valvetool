"""Coil thermal library (P2): resistance-temperature model, steady-state
self-heating equilibrium, and Curie-margin check.

Model scope (spec §2.5):
- Isothermal per opening event: a single ~5 ms opening is far shorter than
  the coil thermal time constant (seconds), so T is an input parameter and
  R is constant during one simulation.
- R(T): linear alpha in the spec range T >= 233.15 K (-40 degC); below,
  an engineering copper-resistivity table (piecewise linear, seam anchored
  to the linear branch at -40 degC for continuity) down to 77 K (LN2).
  Below 77 K -> ValueError (outside table domain).
- Self-heating: steady-state balance I^2 R(T) = conduction + radiation
  (vacuum: no convection, spec §2.5), constant-current worst case.
"""
import math
from dataclasses import dataclass

from scipy.optimize import brentq

STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m^2 K^4)


@dataclass(frozen=True)
class CoilMaterial:
    alpha: float    # temperature coefficient of resistance, 1/K (linear branch)
    T_ref: float    # reference temperature of params.R_coil_20C, K
    r_table: tuple  # ((T [K], R/R0), ...) ascending; low-temperature branch

    def __post_init__(self):
        if not (math.isfinite(self.alpha) and math.isfinite(self.T_ref)):
            raise ValueError(
                f"alpha={self.alpha}, T_ref={self.T_ref}: both must be finite")
        if len(self.r_table) < 2:
            raise ValueError(
                f"r_table has {len(self.r_table)} entry; needs at least 2 "
                "(T, R/R0) points to interpolate")
        for T, r in self.r_table:
            if not (math.isfinite(T) and math.isfinite(r)):
                raise ValueError(f"r_table entry ({T}, {r}) must be finite")
            if r <= 0.0:
                raise ValueError(f"r_table ratio {r} at T={T} K must be > 0")
        temps = [T for T, _ in self.r_table]
        if any(t1 <= t0 for t0, t1 in zip(temps, temps[1:])):
            raise ValueError(
                f"r_table temperatures {temps} must be strictly ascending")


_ALPHA_CU = 0.00393  # 1/K, copper (spec §2.5)
_T_REF = 293.15      # K = 20 degC, definition point of R_coil_20C
_T_SEAM = 233.15     # K = -40 degC, linear branch ends / table begins

# Engineering copper R/R0 values, ~2 significant figures. The 77 K point
# depends on wire purity (RRR) at the +-10% level. The seam entry equals
# the linear branch at -40 degC so the two branches join continuously.
COPPER = CoilMaterial(
    alpha=_ALPHA_CU,
    T_ref=_T_REF,
    r_table=(
        (77.0, 0.13),
        (100.0, 0.21),
        (150.0, 0.41),
        (200.0, 0.62),
        (_T_SEAM, 1.0 + _ALPHA_CU * (_T_SEAM - _T_REF)),
    ),
)


def R_coil(T, params, material=COPPER):
    """Coil resistance at temperature T [K], ohm.

    Linear-alpha branch above the -40 degC seam; piecewise-linear table
    below, floored at 77 K (ValueError outside). Accuracy of the linear
    branch degrades above ~473 K; it is left unbounded only so that
    equilibrium_temp can probe for thermal runaway. At T == material.T_ref
    the multiplier is exactly 1.0, returning params.R_coil_20C bit-exactly.
    """
    if math.isnan(T):
        raise ValueError("T is NaN; temperature must be a real number in K")
    table = material.r_table
    if T < table[0][0]:
        raise ValueError(
            f"T={T} K is below the resistivity table floor "
            f"{table[0][0]} K (model domain)")
    if T >= table[-1][0]:
        return params.R_coil_20C * (1.0 + material.alpha * (T - material.T_ref))
    for (t0, r0), (t1, r1) in zip(table, table[1:]):
        if t0 <= T <= t1:
            return params.R_coil_20C * (r0 + (r1 - r0) * (T - t0) / (t1 - t0))
    raise AssertionError("unreachable: ascending table and T within range")


def equilibrium_temp(I_hold, T_amb, params, material=COPPER, T_max=1000.0):
    """Steady-state coil temperature [K] under constant hold current.

    Solves I^2 R(T) = G_th_cond (T - T_amb)
                      + emissivity * sigma * A_rad * (T^4 - T_amb^4)
    (vacuum: conduction + radiation only). Constant current is the
    conservative worst case — heating grows with T, so runaway is possible;
    constant-voltage drive (P = V^2/R) falls with T and self-stabilizes.
    Raises ValueError when no equilibrium exists below T_max.
    """
    if I_hold == 0.0:
        return T_amb

    def residual(T):
        losses = (params.G_th_cond * (T - T_amb)
                  + params.emissivity * STEFAN_BOLTZMANN * params.A_rad
                  * (T**4 - T_amb**4))
        return I_hold**2 * R_coil(T, params, material) - losses

    if residual(T_max) > 0.0:
        raise ValueError(
            f"thermal runaway: no equilibrium below T_max={T_max} K for "
            f"I_hold={I_hold} A (check G_th_cond/emissivity/A_rad)")
    return brentq(residual, T_amb, T_max, xtol=1e-6)


def curie_margin(T, T_curie=1043.15):
    """Margin to the Curie point [K] (soft iron ~770 degC = 1043.15 K).

    Check-only helper: B_sat degradation with temperature is NOT modeled
    (operating-range margin is huge; documented honesty note).
    """
    return T_curie - T
