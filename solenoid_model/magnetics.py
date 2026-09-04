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
    bisection cannot miss it. The bracket is [0, B_sat]: both curve
    types resolve H(B_sat) to math.inf exactly at B_sat (not just as B
    approaches it), so hi=B_sat is a valid, never-crossed upper bound
    for either AnalyticBH (H(B) diverges continuously as B -> B_sat) or
    TabulatedBH (H(B) is finite everywhere below B_sat, since _interp
    clamps, and only becomes inf via the B >= B_sat check in H_of_B
    itself). Anchoring hi exactly at B_sat, rather than pulling it back
    by an epsilon, is what makes the bracket sound for a clamped table:
    a table can supply only as much MMF as H_of_B(B_sat-eps) gives, and
    that can be far short of NI, so the residual would never cross NI
    inside a bracket that stops short of B_sat.

    After the loop we back-substitute and check the MMF balance holds.
    If NI exceeds what the curve can supply at all (a demand beyond
    what even B=B_sat delivers -- only possible for a clamped table,
    since an analytic law's H diverges first), bisection still returns
    something (lo walks up to ~B_sat), but that value does not actually
    balance the loop, so callers would silently get a physically wrong
    B. Raise instead of returning a number that looks converged but
    isn't; this also guards any future curve type with the same
    clamp-instead-of-diverge shape.

    Sign convention: current magnitude only. Force goes as B**2, so the
    sign of i does not change the force, and a negative i would only
    break the bracket.
    """
    NI = abs(params.N_turns * i)
    if NI == 0.0:
        return 0.0
    lo, hi = 0.0, mag.B_sat
    for _ in range(_BISECT_ITERS):
        mid = 0.5 * (lo + hi)
        if mag.H_of_B(mid) * params.l_core + mid * gap / MU_0 < NI:
            lo = mid
        else:
            hi = mid
    mmf_achieved = mag.H_of_B(lo) * params.l_core + lo * gap / MU_0
    if not math.isclose(mmf_achieved, NI, rel_tol=1e-6, abs_tol=1e-6):
        raise ValueError(
            f"solve_flux_density: demanded MMF {NI:.3f} A-t exceeds what "
            f"this B-H curve can supply (achieved {mmf_achieved:.3f} A-t "
            f"at B={lo:.6g} T, curve B_sat={mag.B_sat:.6g} T); the curve's "
            "table or law does not reach far enough to balance this "
            "current at this gap")
    return lo


def flux_density(gap, i, params, mag=None):
    """Gap flux density [T]. mag=None -> linear circuit (unchanged).

    Sign asymmetry: with mag=None, a negative i returns a signed B of
    matching sign. With mag supplied, solve_flux_density takes abs(N*i),
    so the result is always the non-negative magnitude regardless of the
    sign of i. Harmless for this codebase since force uses B**2, but a
    trap for a future caller that reads B directly and expects a sign.
    """
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

    With `mag` supplied this is always False by construction: the
    result comes from the NONLINEAR path (flux_density delegates to
    solve_flux_density), whose bisection bracket cannot return B >= B_sat.
    So passing `mag` here always reports "not saturated" regardless of
    operating point; the useful comparison is calling this twice, once
    with mag=None (the linear model, which can exceed B_sat) and once
    with mag supplied, to see whether the nonlinear model changes the
    answer at this operating point.
    """
    return flux_density(gap, i, params, mag) > params.B_sat
