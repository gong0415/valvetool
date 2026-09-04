"""Drive-point solvers: what coil current does a given gap require?

Every function here inverts the same force law that dynamics.py integrates
forward, in one of two modes selected by the optional `mag` argument.
`mag=None` (the default) is the linear-magnetics closed form,

    F_mag(gap, i) = 0.5 * i**2 * |dL/dgap|   ->   i = sqrt(2*F / |dL/dgap|)

so a design point solved here and a linear point simulated there cannot
disagree. Passing a B-H curve (see magnetization.py) switches to the
saturating circuit instead: B is nonlinear in i there, so there is no closed
form, and the current is found by bracketing and solving numerically against
magnetics.magnetic_force_closing. There are no empirical constants beyond
that curve: the force law, the spring law, and KVL are the whole content.
Force margin is a *design choice* and is therefore always an explicit
argument, never a default hiding a judgement call.

Why this module exists separately from limits.py: limits.py answers "how good
can this valve possibly be" (bounds, Pareto fronts). These functions answer
"what current does this operating point need" -- a sizing question. The two
share the inversion, so motion_threshold_current now delegates to
pull_in_current rather than repeating the algebra.

The B_sat caveat: the linear form (`mag=None`) is optimistic above B_sat --
the real core produces less force than the formula promises, so a current
solved without `mag` is a lower bound on what the core actually needs, not
an exact answer. This bites hardest at HOLD currents, not peak ones: the
held-open gap is the CLOSED magnetic gap, which is exactly where the linear
model most overstates force, so a hold current sized without `mag` buys less
margin than it claims (see hold_current's docstring for the concrete case).
A solved PEAK current should still be checked with
magnetics.saturation_check; a solved HOLD current should be solved with
`mag` supplied in the first place rather than checked after the fact.
"""
import math

from scipy.optimize import brentq

from solenoid_model import dynamics, magnetics


def current_for_force_margin(gap, margin, params, F_resist, mag=None):
    """Coil current [A] whose force at `gap` is `margin` * `F_resist`.

    Linear mode inverts the force law directly. Saturating mode cannot:
    B is nonlinear in i, so there is no closed form -- it brackets and
    solves numerically instead. Force is monotonic in current in both
    modes, so the root is unique.
    """
    if margin <= 0.0:
        raise ValueError(f"margin must be positive, got {margin}")
    if F_resist <= 0.0:
        raise ValueError(f"F_resist must be positive, got {F_resist}")
    F_target = margin * F_resist
    if mag is None:
        dL_dx = abs(magnetics.dinductance_dgap(gap, params))
        return math.sqrt(2.0 * F_target / dL_dx)

    def residual(i):
        return magnetics.magnetic_force_closing(gap, i, params, mag) - F_target

    # The linear solution is an underestimate of what a saturating core
    # needs (saturation only removes force), so it is a safe lower
    # bracket. Double upward until the force target is cleared; the
    # B_sat ceiling means an unreachable target must terminate rather
    # than loop, hence the explicit cap.
    lo = math.sqrt(2.0 * F_target
                   / abs(magnetics.dinductance_dgap(gap, params)))
    hi = lo
    for _ in range(60):
        if residual(hi) > 0.0:
            return brentq(residual, lo, hi, xtol=1e-12, rtol=1e-12)
        hi *= 2.0
    raise ValueError(
        f"force {F_target:.3f} N at gap {gap:.2e} m is unreachable: the core "
        f"saturates at B_sat={mag.B_sat} T before the coil can produce it")


def hold_current(params, margin, mag=None):
    """Current [A] needed to hold the armature open at `margin` margin.

    The held-open position is x = x_stroke, so the gap is g0 - x_stroke
    (the CLOSED magnetic gap -- the valve is open when the magnetic
    circuit is shut) and the load is the spring at full compression plus
    static pressure.

    Saturation bites hardest exactly here. The closed gap is where the
    linear model most overstates force, so a hold current solved without
    `mag` buys less margin than it claims: on the N2_25BAR_PH case the
    linear 11.44 mA delivers 1.60x, not the 2.0x it was sized for.
    """
    gap_open = params.g0 - params.x_stroke
    F_resist = (dynamics.spring_force(params.x_stroke, params)
                + params.delta_P * params.A_seat)
    return current_for_force_margin(gap_open, margin, params, F_resist, mag)


def pull_in_current(params, margin, mag=None):
    """Current [A] needed to start the armature moving from rest."""
    F_resist = (dynamics.spring_force(0.0, params)
                + params.delta_P * params.A_seat)
    return current_for_force_margin(params.g0, margin, params, F_resist, mag)


def force_margin_at_current(params, i, mag=None):
    """Force margin (F_mag / F_resist) at the held-open position."""
    gap_open = params.g0 - params.x_stroke
    F_resist = (dynamics.spring_force(params.x_stroke, params)
                + params.delta_P * params.A_seat)
    return magnetics.magnetic_force_closing(
        gap_open, i, params, mag) / F_resist


def max_ripple_fraction(margin):
    """Largest fractional current sag a hold designed at `margin` survives.

    Force goes as i**2, so a coil sized for `margin` reaches bare balance when
    the current falls to 1/sqrt(margin) of nominal:

        r_max = 1 - 1/sqrt(margin)

    Note what is NOT in that expression: the valve. Gap, turns, seat area,
    spring rate all cancel with the force ratio, so this is a property of the
    chosen design margin alone and holds for any solenoid. A 2.0x force margin
    tolerates 29.3% sag and no more -- which is a good deal less headroom than
    "2x" sounds like, because the margin is quadratic in the quantity that
    actually ripples.

    This is the number to check a driver's ripple spec against. It does not
    predict ripple: that takes the switching model this module omits.
    """
    if margin < 1.0:
        raise ValueError(
            f"margin={margin} is already below drop-out; a ripple budget "
            f"is undefined")
    return 1.0 - 1.0 / math.sqrt(margin)


def hold_duty(i_hold, params, R=None):
    """Duty ratio needed to average `i_hold` [A] from the bus, as a fraction.

    Straight from KVL on the cycle average: the coil's electrical time
    constant is far longer than any sane chopper period, so the winding
    integrates the switched waveform and sees only its mean, D * V_bus. Then
    D = i*R/V_bus.

    This is the AVERAGE, and says nothing about ripple: the instantaneous
    current swings either side of i_hold, and it is the LOW excursion that
    decides whether the armature stays held. Sizing that band needs a
    switching model, which this module deliberately does not have.
    """
    if R is None:
        R = params.R_coil_20C
    D = i_hold * R / params.V_bus
    if D > 1.0:
        raise ValueError(
            f"i_hold={i_hold} A needs duty {D:.3f} > 1: the {params.V_bus} V "
            f"bus cannot source it through {R:.1f} ohm")
    return D
