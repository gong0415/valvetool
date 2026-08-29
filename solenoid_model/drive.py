"""Drive-point solvers: what coil current does a given gap require?

Every function here inverts the same linear-magnetics force law that
dynamics.py integrates forward,

    F_mag(gap, i) = 0.5 * i**2 * |dL/dgap|   ->   i = sqrt(2*F / |dL/dgap|)

so a design point solved here and a point simulated there cannot disagree.
There are no empirical constants in this module: the force law, the spring
law, and KVL are the whole content. Force margin is a *design choice* and is
therefore always an explicit argument, never a default hiding a judgement
call.

Why this module exists separately from limits.py: limits.py answers "how good
can this valve possibly be" (bounds, Pareto fronts). These functions answer
"what current does this operating point need" -- a sizing question. The two
share the inversion, so motion_threshold_current now delegates to
pull_in_current rather than repeating the algebra.

The B_sat caveat: this inversion is linear-magnetics only. Above B_sat the
real core produces less force than the formula promises, so a solved current
is optimistic there. Hold currents sit at the far end of that concern (the
closed gap is cheap, the currents are milliamps), but a solved PEAK current
should be checked with magnetics.saturation_check before it is trusted.
"""
import math

from solenoid_model import dynamics, magnetics


def current_for_force_margin(gap, margin, params, F_resist):
    """Coil current [A] whose magnetic force at `gap` is `margin` * `F_resist`.

    The direct inversion of magnetics.magnetic_force_closing. `margin` is the
    force safety factor (1.0 = bare balance, nothing to spare); `F_resist` is
    the load being overcome, in newtons.
    """
    if margin <= 0.0:
        raise ValueError(f"margin must be positive, got {margin}")
    if F_resist <= 0.0:
        raise ValueError(f"F_resist must be positive, got {F_resist}")
    dL_dx = abs(magnetics.dinductance_dgap(gap, params))
    return math.sqrt(2.0 * margin * F_resist / dL_dx)


def hold_current(params, margin):
    """Current [A] needed to hold the armature open at `margin` force margin.

    The held-open position is x = x_stroke, so the gap is g0 - x_stroke (the
    CLOSED magnetic gap -- the valve is open when the magnetic circuit is
    shut) and the load is the spring at full compression plus static pressure.

    This is the cheap half of peak-and-hold: the same force that costs
    pull_in_current at the rest gap costs a fraction of it here, because
    |dL/dgap| grows sharply as the gap closes.
    """
    gap_open = params.g0 - params.x_stroke
    F_resist = (dynamics.spring_force(params.x_stroke, params)
                + params.delta_P * params.A_seat)
    return current_for_force_margin(gap_open, margin, params, F_resist)


def pull_in_current(params, margin):
    """Current [A] needed to start the armature moving from rest, at `margin`.

    Evaluated at the rest gap g0 against spring preload + static pressure.
    At margin=1.0 this is exactly limits.motion_threshold_current -- the
    current below which the valve cannot open however long you wait.
    """
    F_resist = (dynamics.spring_force(0.0, params)
                + params.delta_P * params.A_seat)
    return current_for_force_margin(params.g0, margin, params, F_resist)


def force_margin_at_current(params, i):
    """Force margin (F_mag / F_resist) at the held-open position for current i.

    The forward reading of hold_current: pass the current back in and get the
    margin out. Below 1.0 the armature drops out.
    """
    gap_open = params.g0 - params.x_stroke
    F_resist = (dynamics.spring_force(params.x_stroke, params)
                + params.delta_P * params.A_seat)
    return magnetics.magnetic_force_closing(gap_open, i, params) / F_resist


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
