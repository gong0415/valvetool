"""Magnetic circuit: reluctance, inductance, coenergy force, saturation.

Two modes, following the dual-mode rule in ARCHITECTURE.md. `mag=None`
(the default) is the original linear circuit with a constant
mu_r_core -- bit-identical to every result this module produced before
P8, which the characterization tests freeze. Passing a
magnetization.AnalyticBH or TabulatedBH switches to a nonlinear solve
where the core saturates.

Why this matters: the linear model reports B = 8.16 T at the closed gap
of the N2_25BAR_PH case, which is not a physical flux density. It
overstates force where force is scarcest, and it makes the force look
like it explodes as the gap closes when the real curve is nearly flat.

reluctance/inductance/dinductance_dgap keep their original signatures
and stay linear-only on purpose: with a saturating core there is no
single reluctance -- it depends on the operating point -- so a `mag`
argument there would return a number that quietly means something else.
Nonlinear callers go through solve_flux_density instead.
"""
import math

from solenoid_model import magnetization

MU_0 = 4 * math.pi * 1e-7  # vacuum permeability, H/m

# Bisection bracket is [0, B_sat) and the residual is monotonic in B, so
# a fixed iteration count converges to machine precision without needing
# a tolerance argument: 200 halvings of a ~2 T bracket is far below
# double precision.
_BISECT_ITERS = 200


def reluctance(gap, params):
    r_gap = gap / (MU_0 * params.A_gap)
    r_core = params.l_core / (params.mu_r_core * MU_0 * params.A_gap)
    return r_gap + r_core


def inductance(gap, params):
    return params.N_turns ** 2 / reluctance(gap, params)


def dinductance_dgap(gap, params):
    r_total = reluctance(gap, params)
    dr_dgap = 1.0 / (MU_0 * params.A_gap)
    return -params.N_turns ** 2 * dr_dgap / r_total ** 2


def solve_flux_density(gap, i, params, mag):
    """Gap flux density [T] with a saturating core.

    Solves the MMF balance around the loop,

        H(B)*l_core + B*gap/MU_0 = N*i

    for B by bisection. The left side is strictly increasing in B (H(B)
    is monotonic and the gap term is linear), so the root is unique and
    bisection cannot miss it. The bracket is [0, B_sat): H(B_sat) is
    infinite, so the upper end is approached, never reached.

    Sign convention: current magnitude only. Force goes as B**2, so the
    sign of i does not change the force, and a negative i would only
    break the bracket.
    """
    NI = abs(params.N_turns * i)
    if NI == 0.0:
        return 0.0
    lo, hi = 0.0, mag.B_sat * (1.0 - 1e-12)
    for _ in range(_BISECT_ITERS):
        mid = 0.5 * (lo + hi)
        if mag.H_of_B(mid) * params.l_core + mid * gap / MU_0 < NI:
            lo = mid
        else:
            hi = mid
    return lo


def flux_density(gap, i, params, mag=None):
    """Gap flux density [T]. mag=None -> linear circuit (unchanged)."""
    if mag is None:
        flux = params.N_turns * i / reluctance(gap, params)
        return flux / params.A_gap
    return solve_flux_density(gap, i, params, mag)


def magnetic_force_closing(gap, i, params, mag=None):
    """Force [N] pulling the armature toward the pole, always >= 0.

    Linear mode keeps the coenergy form 0.5*i^2*|dL/dgap| untouched.
    Saturating mode uses the Maxwell stress form B^2*A/(2*mu_0), which
    is the same physics: with a linear circuit the two agree exactly,
    but only the B-form stays valid once the core saturates, because
    dL/dgap is no longer a constant of the operating point.
    """
    if mag is None:
        return -0.5 * i ** 2 * dinductance_dgap(gap, params)
    B = solve_flux_density(gap, i, params, mag)
    return B ** 2 * params.A_gap / (2.0 * MU_0)


def saturation_check(gap, i, params, mag=None):
    """True when the core is driven past B_sat.

    With `mag` supplied this is always False by construction -- the
    solver cannot return B >= B_sat -- so it reports whether the LINEAR
    model would have exceeded B_sat, i.e. whether using `mag` changes
    the answer materially at this operating point.
    """
    return flux_density(gap, i, params, mag) > params.B_sat
