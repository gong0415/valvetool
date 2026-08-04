"""Seat-impact contact mechanics and bounce (spec §2.2 seat bounce, §3 L5).

Never imports dynamics -- pure algebraic helpers. The dependency runs the
other way: dynamics.py imports this module for its bounce event handling,
and limits.py imports it for the L5 life curve.

All formulae are the standard Hertz elastic-impact closed forms for a sphere
of radius R_tip striking a flat, with F = k*delta**1.5.
"""
import math


def hertz_stiffness(seat, R_tip):
    """Hertz contact stiffness k [N/m**1.5]: k = (4/3)*E*·sqrt(R_tip)."""
    return (4.0 / 3.0) * seat.E_star * math.sqrt(R_tip)


def indentation_max(m, v, seat, R_tip):
    """Peak elastic indentation [m] from the energy balance
    (1/2)*m*v**2 = (2/5)*k*delta**2.5  ->  delta = (5*m*v**2/(4*k))**0.4."""
    k = hertz_stiffness(seat, R_tip)
    return (5.0 * m * v**2 / (4.0 * k)) ** 0.4


def contact_force_max(m, v, seat, R_tip):
    """Peak contact force [N]: F = k*delta_max**1.5."""
    k = hertz_stiffness(seat, R_tip)
    return k * indentation_max(m, v, seat, R_tip) ** 1.5


def contact_pressure_max(m, v, seat, R_tip):
    """Peak Hertzian contact pressure [Pa]:
    p_max = (6*F_max*E*^2/(pi**3*R_tip**2))**(1/3), so p_max ~ v**0.4."""
    F = contact_force_max(m, v, seat, R_tip)
    return (6.0 * F * seat.E_star**2 / (math.pi**3 * R_tip**2)) ** (1.0 / 3.0)


def contact_duration(m, v, seat, R_tip):
    """Hertz impact contact time [s]: t_c = 2.87*(m**2/(k**2*v))**0.2.

    Used as the modelling self-consistency check behind the instantaneous
    coefficient-of-restitution assumption: t_c must stay well below the
    stroke time for that idealisation to hold."""
    k = hertz_stiffness(seat, R_tip)
    return 2.87 * (m**2 / (k**2 * v)) ** 0.2


def is_shakedown(m, v, seat, R_tip):
    """True when p_max < H -- H being seat.H, the FULLY-PLASTIC indentation
    hardness, not the elastic limit. Classical Hertz contact puts first
    subsurface yield at p0 ~= 1.60*Y while H ~= 3*Y (Tabor's rule), so first
    yield actually starts around p0 ~= 0.53*H: every impact with p_max
    between roughly 0.53*H and H IS already accumulating plastic
    indentation, this function just does not say so. Using the fully-plastic
    threshold rather than the first-yield pressure as the branch test is an
    EMPIRICAL simplification, and it is the unconservative one -- callers
    (l5_life_curve) report N_cycle = +inf for the whole "True" band,
    including the part that is already yielding. See MODEL_NOTES.md's
    known-limitations item 9 for the measured extent of that band on the
    baseline valve."""
    return contact_pressure_max(m, v, seat, R_tip) < seat.H


# EMPIRICAL cut-offs for the bounce train: below V_MIN_DEFAULT the rebound
# carries negligible energy, and N_MAX_DEFAULT is the Zeno guard.
V_MIN_DEFAULT = 0.01   # m/s
N_MAX_DEFAULT = 50


def bounce_sequence(v0, e, v_min=V_MIN_DEFAULT, n_max=N_MAX_DEFAULT):
    """Impact speeds [m/s] of the bounce train: v0, e*v0, e**2*v0, ...

    Terminates when the NEXT rebound would fall to or below v_min, so an
    e -> 0 seat yields exactly one impact. dynamics.py's bounce loop uses
    this same rule, which is what lets the two be cross-checked.
    """
    speeds = [abs(v0)]
    while len(speeds) < n_max:
        v_next = e * speeds[-1]
        if v_next <= v_min:
            break
        speeds.append(v_next)
    return speeds


def impact_energy_total(m, v0, e, v_min=V_MIN_DEFAULT, n_max=N_MAX_DEFAULT):
    """Total kinetic energy [J] dissipated across the bounce train. Tends to
    the initial kinetic energy, the residue being the truncated tail."""
    return sum(0.5 * m * v**2 * (1.0 - e**2)
               for v in bounce_sequence(v0, e, v_min, n_max))
