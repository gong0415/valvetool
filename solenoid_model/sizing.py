"""Valve sizing and the L3 miniaturisation limit (spec §3 L3).

Never imports dynamics -- these are pure algebraic helpers, the same
invariant winding.py holds. m_arm remains an independent ValveParams input
consumed by the ODE; nothing here overrides it.

Prototyping showed the three mechanisms the spec names for L3 (winding
window, magnetic saturation, seat dP force) produce no miniaturisation limit
on their own: the force balance is exactly scale-invariant (both the required
seat force and the available magnetic force go as D**2), and the MMF balance
stays flat once the current is set by the thermal model. The real limit comes
from an engineering ceiling on current density -- see J_max in ValveParams
and docs/spec/PARAMS.md.
"""
import math
from dataclasses import replace

MU_0 = 4.0 * math.pi * 1e-7  # vacuum permeability, H/m
RHO_FE = 7870.0              # kg/m^3, soft iron (magnetic circuit)
RHO_CU = 8960.0              # kg/m^3, copper (winding)


def active_mass(params):
    """Magnetically and electrically active mass [kg]: the iron of the
    magnetic circuit plus the copper of the winding. Excludes housing,
    seat, spring and fasteners -- see valve_mass for the packaged figure."""
    m_core = params.A_gap * params.l_core * RHO_FE
    m_coil = params.A_winding * params.l_turn_mean * RHO_CU
    return m_core + m_coil


def valve_mass(params, k_pack):
    """Packaged valve mass [kg] = active_mass * k_pack.

    k_pack is an EMPIRICAL packaging/housing factor, not a derived quantity;
    it is calibrated so the baseline lands inside the 30-500 g band spec §4
    quotes for m_valve."""
    return active_mass(params) * k_pack


def pressure_ceiling(params):
    """Largest seat pressure differential [Pa] the magnetic circuit can hold
    shut: dP_max = B_sat**2 * A_gap / (2*mu_0*A_seat).

    Scale-invariant under self-similar scaling -- A_gap and A_seat both go as
    D**2 and cancel -- so shrinking the valve does not lower its pressure
    rating."""
    return params.B_sat**2 * params.A_gap / (2.0 * MU_0 * params.A_seat)


def min_feasible_scale(params, J_max):
    """Smallest self-similar scale factor s at which the winding can still
    drive the working gap to saturation, given a current-density ceiling
    J_max [A/m**2]:

        s* = B_sat * g0 / (mu_0 * J_max * A_winding * k_fill)

    The MMF the gap needs goes as D (B_sat*gap/mu_0) while the MMF the window
    can supply goes as D**2 (J_max*A_winding*k_fill), so the ratio goes as D
    and there is a hard floor. Below s* no winding fits enough ampere-turns.

    This is a gap-only magnetic circuit: NI_need omits the iron path
    (magnetics.reluctance's l_core/(mu_r_core*mu_0*A_gap) term), which this
    module does not import. At baseline that term is 12.5 um against a
    350 um gap, i.e. about 3.6% of it, so NI_need here is understated by the
    same ~3.6% and s* is understated too: including the iron path would move
    s* from 1.0140 to 1.0502 at J_max=20 A/mm**2. Both terms scale as D, so
    no scaling law above changes -- only the numeric prefactor would. Not
    corrected here by design: re-pinning s* would ripple through every
    literal derived from it (tests, docs, the regenerated figure) for a
    correction smaller than the EMPIRICAL J_max ceiling's own uncertainty."""
    NI_need_unit = params.B_sat * params.g0 / MU_0
    NI_avail_unit = J_max * params.A_winding * params.k_fill
    return NI_need_unit / NI_avail_unit


def scale_params(params, s):
    """Self-similar geometric scaling by factor s, for scaling-law analysis.

    Lengths go as s, areas as s**2, armature mass as s**3. A geometrically
    similar helical spring has k = G*d**4/(8*D**3*n), so k_spring goes as s
    and the preload k*x0 goes as s**2. Conduction conductance goes as s
    (cross-section over path length) and radiating area as s**2.

    Electrical drive and material properties (V_bus, N_turns, B_sat,
    mu_r_core, delta_P, k_fill, ...) are deliberately left untouched: this
    scales the hardware, not the operating point.

    R_coil_20C is hardware too, but it is deliberately NOT rescaled here: the
    winding window (A_winding, l_turn_mean) does scale, so a params object
    returned by this function is winding-inconsistent by construction --
    R_coil_20C still describes the *unscaled* winding. Any caller doing a
    tau_e or pi_2 (dimensionless_groups) analysis on a scaled params object
    must supply R explicitly, derived from the scaled window via
    winding.coil_resistance(N_turns, ps.A_winding, ps.k_fill,
    ps.l_turn_mean) -- passing the stale R_coil_20C default silently answers
    a different, physically inconsistent question."""
    return replace(
        params,
        g0=params.g0 * s,
        x_stroke=params.x_stroke * s,
        l_core=params.l_core * s,
        l_turn_mean=params.l_turn_mean * s,
        D_seat_bore=params.D_seat_bore * s,
        A_gap=params.A_gap * s**2,
        A_seat=params.A_seat * s**2,
        A_winding=params.A_winding * s**2,
        A_rad=params.A_rad * s**2,
        m_arm=params.m_arm * s**3,
        k_spring=params.k_spring * s,
        F_preload=params.F_preload * s**2,
        G_th_cond=params.G_th_cond * s,
    )


def l3_feasible_region(params, J_max, k_pack, scales, pressures):
    """L3 feasible region over a (scale, pressure) grid.

    A cell is feasible when the winding can still drive the gap
    (scale >= min_feasible_scale) and the magnetic circuit can hold the seat
    shut (pressure <= pressure_ceiling). The two limits are orthogonal --
    one constrains size, the other pressure, and neither depends on the
    other -- so the feasible region is a rectangle.

    That rectangle is a consequence of the size criterion, not an
    independent physical result: min_feasible_scale requires the winding to
    drive the gap all the way to B_sat regardless of the pressure the valve
    actually has to hold, which is exactly what makes the size axis
    ΔP-independent. A demand-matched criterion (MMF sufficient only for
    F_preload + delta_P*A_seat at the rest gap, not full saturation) would
    couple the two axes -- s* becomes ΔP-dependent (0.598 @ 1 MPa, 0.846 @
    4 MPa, 1.092 @ 8 MPa) and the boundary becomes a curve, not a rectangle;
    at the baseline delta_P it comes out to s* ~ 0.72, about 1.4x smaller
    than the to-saturation s* = 1.014 this module ships. The to-saturation
    criterion is the deliberate choice (it is the worst case, independent of
    the operating pressure); the rectangle is downstream of that choice.

    Returns index-aligned structures: feasible[i][j] for scales[i] x
    pressures[j], and mass[i] for scales[i] (mass depends only on scale)."""
    s_star = min_feasible_scale(params, J_max)
    dP_max = pressure_ceiling(params)
    feasible = [[bool(s >= s_star and dP <= dP_max) for dP in pressures]
                for s in scales]
    mass = [valve_mass(scale_params(params, s), k_pack) for s in scales]
    return {"scales": list(scales), "pressures": list(pressures),
            "feasible": feasible, "mass": mass}
