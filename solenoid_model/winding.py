"""Winding-window coil model (P4): couples turns, wire gauge, and winding
window into coil resistance. Closes ontology primitive 1 (spec §2.0) — the
wire-gauge / winding-window degrees of freedom the independent N_turns and
R_coil_20C inputs left implicit.

Never imports dynamics/magnetics — these are pure algebraic helpers. The
coil resistance derived here is an analysis layer (used by limits.py L2);
dynamics.py keeps using the independent R_coil_20C input.
"""
import math

RHO_CU_20C = 1.68e-8  # ohm*m, annealed copper at 20 degC


def coil_resistance(N, A_winding, k_fill, l_turn_mean, rho_cu=RHO_CU_20C):
    """Coil resistance [ohm] from winding geometry.

    Each turn gets copper cross-section a_wire = A_winding*k_fill/N; total
    wire length is N*l_turn_mean; so
        R = rho_cu * (N*l_turn_mean) / a_wire
          = rho_cu * N**2 * l_turn_mean / (A_winding * k_fill).
    For a fixed window this makes N**2/R constant — the invariant behind
    the L2 power-vs-force slope.
    """
    return rho_cu * N**2 * l_turn_mean / (A_winding * k_fill)


def wire_diameter(N, A_winding, k_fill):
    """Round-wire diameter [m] for N turns filling k_fill of the window."""
    a_wire = A_winding * k_fill / N
    return math.sqrt(4.0 * a_wire / math.pi)


def turns_for_resistance(R, A_winding, k_fill, l_turn_mean, rho_cu=RHO_CU_20C):
    """Inverse of coil_resistance: turns N giving resistance R [ohm]."""
    return math.sqrt(R * A_winding * k_fill / (rho_cu * l_turn_mean))


def coil_resistance_consistency(params):
    """Compare the declared coil resistance against the winding geometry.

    Returns (R_declared, R_geometric, rel_error). params.R_coil_20C is an
    independent input that dynamics uses directly, while this module derives
    resistance from the winding window; they are calibrated to agree at the
    baseline but nothing enforces that. Changing N_turns alone, for instance,
    leaves R_coil_20C stale -- 80 ohm declared against 320 ohm of geometry at
    N=4000 -- and the simulation would silently model a coil that cannot be
    wound. Reports rather than raises: ODE solvers probe wide parameter
    ranges, so callers decide what divergence is acceptable."""
    R_declared = params.R_coil_20C
    R_geometric = coil_resistance(params.N_turns, params.A_winding,
                                  params.k_fill, params.l_turn_mean)
    return R_declared, R_geometric, abs(R_geometric - R_declared) / R_declared
