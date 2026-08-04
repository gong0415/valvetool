"""L1-L6 physical-limit scans (spec §3). L6 pull-in voltage arrived with P2,
L1/L2 with P4, L3's groups with P5, and L4/L5 with P6 -- the set is now
complete."""
import math
from dataclasses import replace

from solenoid_model import (dynamics, fluid, impact, magnetics, sealing,
                            winding)

MU_0 = 4 * math.pi * 1e-7  # vacuum permeability, H/m


def find_pullin_voltage(params, T_coil=None, medium=None, cond=None,
                        t_max=0.05, tol=0.05):
    """Minimum bus voltage [V] at which the valve completes its stroke
    within t_max [s].

    Dynamic criterion: success means the stroke event of a full coupled
    simulation fires (back-EMF included). Bisection over V; the upper
    bracket starts at params.V_bus and doubles until the valve opens
    (ValueError beyond 4x V_bus). Returns the successful upper endpoint,
    within tol [V] of the threshold. Cost: ~log2(V_bus/tol) simulations;
    sub-threshold probes run the full t_max, so pass a reduced t_max for
    interactive or test use.
    """
    if params.V_bus <= 0.0:
        raise ValueError(
            f"V_bus={params.V_bus} V must be > 0: it seeds the doubling "
            "bracket of the bisection search")

    def opens(V):
        p = replace(params, V_bus=V)
        sol = dynamics.simulate_opening(p, t_max=t_max, medium=medium,
                                        cond=cond, T_coil=T_coil)
        return sol.t_events[0].size > 0

    hi = params.V_bus
    while not opens(hi):
        hi *= 2.0
        if hi > 4.0 * params.V_bus:
            raise ValueError(
                f"no pull-in below {4.0 * params.V_bus} V within "
                f"t_max={t_max} s")
    lo = 0.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if opens(mid):
            hi = mid
        else:
            lo = mid
    return hi


def saturation_force_ceiling(params):
    """Maxwell holding-force ceiling [N] at magnetic saturation:
    F_sat = B_sat**2 * A_gap / (2 * mu_0). Real holding force cannot exceed
    this no matter how much coil power is applied."""
    return params.B_sat**2 * params.A_gap / (2.0 * MU_0)


def _linear_holding_force(params, i):
    """Unsaturated (linear-magnetics) holding force 0.5*i**2*|dL/dx| at the
    pulled-in gap (gap = g0 - x_stroke). Fictitious above B_sat — callers cap
    it with saturation_force_ceiling."""
    gap = params.g0 - params.x_stroke
    dL_dx = magnetics.dinductance_dgap(gap, params)
    return 0.5 * i**2 * abs(dL_dx)


def holding_force(params, i):
    """Holding force [N] at the pulled-in gap (gap = g0 - x_stroke) for coil
    current i, capped at the B_sat ceiling.

    Below saturation the linear-magnetics force 0.5*i**2*|dL/dx| applies;
    above it that formula is fictitious (linear magnetics ignores B_sat), so
    the returned value is min(linear, F_sat). For the baseline the coil runs
    deep in saturation at hold (linear ~2364 N vs F_sat ~43 N)."""
    return min(_linear_holding_force(params, i), saturation_force_ceiling(params))


def l2_pareto_front(params, N_values):
    """L2 power-vs-holding-force front over a turns sweep.

    For each N: derive R from the winding window (winding.coil_resistance),
    take the steady hold current i = V_bus / R, hold power P = i**2 * R, the
    fictitious linear force, and the real (B_sat-capped) holding force. The
    winding window fixes F/P slope (N**2/R invariant); the B_sat ceiling caps
    the front. Returns index-aligned lists under keys N/R/P_hold/F_hold/
    F_linear."""
    R_list, P_list, F_list, F_lin_list = [], [], [], []
    for N in N_values:
        R = winding.coil_resistance(N, params.A_winding, params.k_fill,
                                    params.l_turn_mean)
        i = params.V_bus / R
        p_n = replace(params, N_turns=N)
        F_linear = _linear_holding_force(p_n, i)
        R_list.append(R)
        P_list.append(i**2 * R)
        F_lin_list.append(F_linear)
        F_list.append(min(F_linear, saturation_force_ceiling(p_n)))
    return {"N": list(N_values), "R": R_list, "P_hold": P_list,
            "F_hold": F_list, "F_linear": F_lin_list}


def _resisting_force_at_rest(params):
    """Spring preload + static pressure opposing the armature at x=0 [N]."""
    return dynamics.spring_force(0.0, params) + params.delta_P * params.A_seat


def electrical_time_constant(params, R=None, gap=None):
    """Coil electrical time constant tau_e = L(gap)/R [s] (spec §2.1).

    gap defaults to the rest gap g0; R defaults to params.R_coil_20C. At a
    fixed winding window tau_e is invariant with turns: L ~ N**2 (L = N**2/R_m)
    and the winding-derived R ~ N**2, so N cancels. Turns is therefore not a
    lever on the electrical time constant -- the sister result to L2's
    N**2/R invariant."""
    if R is None:
        R = params.R_coil_20C
    if gap is None:
        gap = params.g0
    return magnetics.inductance(gap, params) / R


def motion_threshold_current(params):
    """Coil current [A] at which the magnetic force at the rest gap first
    balances spring preload + static pressure. Below this the armature cannot
    begin to move, however long the current is applied."""
    dL_dx = abs(magnetics.dinductance_dgap(params.g0, params))
    return math.sqrt(2.0 * _resisting_force_at_rest(params) / dL_dx)


def response_time_bound(params, R=None):
    """L1 lower bound on t_open [s], as {"t_i", "t_x", "t_bound"}.

    Built from two terms whose sum is a genuine bound (verified against
    simulation across the voltage/turns/stroke sweeps, tightness ~0.80-0.91),
    but the two terms are not on equal footing:

      t_i = -tau_e * ln(1 - i_th/i_ss)   RL rise to the motion-threshold
                                          current, ignoring back-EMF (which
                                          only ever slows the rise)
      t_x = sqrt(2*m*x_stroke/F_max)     constant-acceleration transit at the
                                          largest force available anywhere in
                                          the stroke (B_sat-capped)

    t_i is EXACT within this model, not merely optimistic: coupled_rhs pins
    the armature at x=0 while dv_dt <= 0, so the pre-motion interval is a
    pure RL rise at L(g0)/R that ends precisely when i = i_th -- there is no
    approximation here for back-EMF to slacken. t_x is the only optimistic
    leg (constant acceleration at the largest force available anywhere in
    the stroke), and it is a loose one: measured against simulation it
    recovers only about 6-26% of the real transit time. Essentially all of
    the bound's slack therefore lives in t_x. Consequently, when the bound
    comes out tight (0.796-0.913 on the shipped sweeps), that tightness is
    mainly because this valve is dead-time dominated (t_i dominates t_open),
    not because the constant-acceleration transit model is an accurate
    description of the real stroke -- it is not evidence for t_x's accuracy.

    Precondition: the bound holds only while F_sat (the saturation force
    ceiling used to cap F_max) exceeds the peak magnetic force the
    (uncapped) ODE actually develops during the stroke; otherwise t_x would
    be computed from a force below what the ODE really sees and would stop
    being optimistic. See the guard test
    test_saturation_cap_does_not_bind_at_the_decisive_thresholds.

    Returns None when the steady current cannot reach the motion threshold
    (the valve never opens, so there is no response time to bound)."""
    if R is None:
        R = params.R_coil_20C
    i_ss = params.V_bus / R
    i_th = motion_threshold_current(params)
    if i_ss <= i_th:
        return None
    tau_e = electrical_time_constant(params, R=R)
    t_i = -tau_e * math.log(1.0 - i_th / i_ss)
    F_max = holding_force(params, i_ss)   # already min(linear, F_sat)
    t_x = math.sqrt(2.0 * params.m_arm * params.x_stroke / F_max)
    return {"t_i": t_i, "t_x": t_x, "t_bound": t_i + t_x}


def naive_response_time_estimate(params, R=None):
    """The spec §3 L1 decomposition tau_e + sqrt(2*m*x/F_net) [s].

    WARNING: this is an order-of-magnitude estimator, NOT a bound in either
    direction. Measured against simulation it ranges from 0.43x (under) to
    3.11x (over) the true t_open, because the armature starts moving long
    before the current settles -- charging the full tau_e as dead time
    double-counts. Use response_time_bound for a claim that actually holds.

    Returns None when the net force at the rest gap is non-positive (the
    valve never opens)."""
    if R is None:
        R = params.R_coil_20C
    i_ss = params.V_bus / R
    F_mag = 0.5 * i_ss**2 * abs(magnetics.dinductance_dgap(params.g0, params))
    F_net = F_mag - _resisting_force_at_rest(params)
    if F_net <= 0.0:
        return None
    tau_e = electrical_time_constant(params, R=R)
    return tau_e + math.sqrt(2.0 * params.m_arm * params.x_stroke / F_net)


def l1_sweep(cases):
    """L1 response-time sweep over explicit parameter variants.

    cases is an iterable of (label, params, R): the caller supplies both the
    parameter variant and the resistance that physically belongs with it. That
    is deliberate -- for a turns sweep R must come from the winding window
    (winding.coil_resistance), and making the choice explicit at the call site
    keeps the physically-inconsistent fixed-R variant clearly labelled rather
    than hidden inside this function.

    Returns index-aligned lists under label/t_open/t_bound/t_naive/pulls_in;
    t_open, t_bound and t_naive are None for cases that never open."""
    labels, t_opens, t_bounds, t_naives, pulls_in = [], [], [], [], []
    for label, params, R in cases:
        sol = dynamics.simulate_opening(replace(params, R_coil_20C=R))
        # bool(): sol.t_events[0].size is a numpy int, so the comparison
        # yields numpy.bool_, which fails callers' `is True` / `is False`
        opened = bool(sol.t_events[0].size > 0)
        bound = response_time_bound(params, R=R)
        labels.append(label)
        t_opens.append(sol.t_events[0][0] if opened else None)
        t_bounds.append(bound["t_bound"] if bound is not None else None)
        t_naives.append(naive_response_time_estimate(params, R=R))
        pulls_in.append(opened)
    return {"label": labels, "t_open": t_opens, "t_bound": t_bounds,
            "t_naive": t_naives, "pulls_in": pulls_in}


def dimensionless_groups(params, R=None):
    """Buckingham-pi groups governing the valve (spec §3).

    pi_1 stroke/gap          x_stroke/g0
    pi_2 time-constant ratio tau_e/tau_m       (electrical vs mechanical)
    pi_3 magnetic number     F_mag/F_spring    (at rest gap, steady current)
    pi_4 pressure load       dP*A_seat/F_spring
    pi_5 saturation margin   F_sat/(dP*A_seat)

    The scaling-law paragraph below holds only when R is supplied from the
    scaled winding window, e.g. R=winding.coil_resistance(N_turns,
    ps.A_winding, ps.k_fill, ps.l_turn_mean) for ps = sizing.scale_params(p,
    s). With the default R=params.R_coil_20C -- which sizing.scale_params
    deliberately does not rescale -- the natural call
    dimensionless_groups(scale_params(p, s)) gives the opposite law: tau_e
    (and hence pi_2) stays constant with s, while pi_3 goes as D**-2 instead
    of being invariant. The default is correct for an UNSCALED params
    object; it is a trap for a scaled one.

    Under self-similar scaling with a winding-consistent R, tau_e goes as
    D**2 and tau_m as D, so pi_2 is the only scale-dependent group: it goes
    as D, meaning a smaller valve shifts from dead-time dominated toward
    mechanically dominated. The other four are scale-invariant."""
    if R is None:
        R = params.R_coil_20C
    tau_e = electrical_time_constant(params, R=R)
    tau_m = 1.0 / math.sqrt(params.k_spring / params.m_arm)
    i_ss = params.V_bus / R
    F_mag = 0.5 * i_ss**2 * abs(magnetics.dinductance_dgap(params.g0, params))
    F_spring = params.F_preload + params.k_spring * params.x_stroke
    F_pressure = params.delta_P * params.A_seat
    return {
        "pi_1_stroke_gap": params.x_stroke / params.g0,
        "pi_2_time_ratio": tau_e / tau_m,
        "pi_3_magnetic_number": F_mag / F_spring,
        "pi_4_pressure_load": F_pressure / F_spring,
        "pi_5_saturation_margin": saturation_force_ceiling(params) / F_pressure,
    }


# EMPIRICAL life-model coefficients. None of these are derived; they set the
# scale of the L5 estimate and should be treated as order-of-magnitude only.
ARCHARD_K = 1.0e-4          # dimensionless wear coefficient, polymer-on-steel
FAILURE_DEPTH_FRACTION = 0.1  # indentation depth reaching this fraction of
                              # w_land is taken as seal-geometry failure


def l4_leak_curve(params, seats, w_lands, cond):
    """L4: contact stress vs helium leak rate over a range of land widths.

    Returns {seat.name: [{w_land, stress, phi, sealed, gap, leak}, ...]}.
    The trustworthy content is the `sealed` flag and the land width at which
    it flips; see sealing.leak_rate on why the magnitudes are not.
    """
    out = {}
    for seat in seats:
        rows = []
        for w in w_lands:
            p = replace(params, w_land=w)
            rows.append({
                "w_land": w,
                "stress": sealing.contact_stress(p, cond),
                "phi": sealing.contact_area_fraction(p, seat, cond),
                "sealed": sealing.is_percolated(p, seat, cond),
                "gap": sealing.residual_gap(p, seat, cond),
                "leak": sealing.leak_rate(p, seat, cond),
            })
        out[seat.name] = rows
    return out


def l5_life_curve(params, seat, velocities, R_tip=None):
    """L5: impact speed vs cycle-life estimate.

    Shakedown first: if the impact stays elastic (p_max < H) no plastic
    indentation accumulates and the wear route reports infinite life -- the
    honest statement being that life is then fatigue-limited by a mechanism
    this model does not carry. Otherwise an Archard-type wear volume gives a
    per-cycle indentation depth, and life is the cycle count at which that
    depth reaches FAILURE_DEPTH_FRACTION of the land width.

    Every coefficient here is EMPIRICAL; treat the result as an
    order-of-magnitude estimate, never as a qualification number.
    """
    R_tip = params.R_tip if R_tip is None else R_tip
    rows = []
    for v in velocities:
        shakedown = impact.is_shakedown(params.m_arm, v, seat, R_tip)
        p_max = impact.contact_pressure_max(params.m_arm, v, seat, R_tip)
        E_imp = 0.5 * params.m_arm * v**2
        if shakedown:
            N = float("inf")
            depth = 0.0
        else:
            # Archard-type wear: worn volume = K * (normal load * slide
            # distance)/H, with the impact's kinetic energy standing in for
            # the load-distance product. The full impact energy is used, not
            # just the plastically dissipated fraction: for a coefficient of
            # restitution e the plastic share is a constant factor (1 - e**2)
            # that ARCHARD_K, itself EMPIRICAL and uncalibrated, absorbs
            # entirely. The distinction would matter only if K were ever
            # calibrated against measurements.
            volume = ARCHARD_K * E_imp / seat.H
            area = math.pi * sealing.seal_diameter(params) * params.w_land
            depth = volume / area
            N = FAILURE_DEPTH_FRACTION * params.w_land / depth
        rows.append({
            "v_impact": v, "p_max": p_max, "shakedown": shakedown,
            "depth_per_cycle": depth, "N_cycle": N,
        })
    return rows


def bounce_leak_spike(params, seat, circuit, cond, medium):
    """Leakage during the closing bounce train (the L4-L5 coupling).

    The only cross-domain composite in P6, and it belongs here rather than in
    sealing.py: while the poppet is off the seat the flow is an orifice
    problem owned by fluid.py, and only the seated intervals are a contact
    problem owned by sealing.py.

    Returns a flat dict whose t_first/t_settle are release-relative (the
    closing BounceResult's own clock, t=0 at stop release), NOT
    de-energization-relative like ClosingResult.t_close/.t_release -- see
    BounceResult's docstring in dynamics.py.
    """
    result = dynamics.simulate_closing(params, circuit, seat=seat)
    bounce = result.bounce
    t_open_total = 0.0
    mass_leaked = 0.0
    for seg in bounce.segments[1:]:
        t_seg = seg.t
        x_seg = seg.y[1]
        for k in range(1, len(t_seg)):
            dt = t_seg[k] - t_seg[k - 1]
            x_mid = 0.5 * (x_seg[k] + x_seg[k - 1])
            if x_mid <= 0.0:
                continue
            t_open_total += dt
            mass_leaked += fluid.mdot_gas(x_mid, params, medium, cond) * dt
    return {
        "n_bounce": bounce.n_bounce,
        "v_impacts": bounce.v_impacts,
        "t_first": bounce.t_first,
        "t_settle": bounce.t_settle,
        "t_open_total": t_open_total,
        "mass_leaked": mass_leaked,
        "sealed_leak": sealing.leak_rate(params, seat, cond),
    }
