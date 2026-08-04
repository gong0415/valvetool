from dataclasses import dataclass

from scipy.integrate import solve_ivp

from solenoid_model import flyback as flyback_lib
from solenoid_model import fluid as fluid_lib
from solenoid_model import impact as impact_lib
from solenoid_model import magnetics
from solenoid_model import thermal


@dataclass
class ClosingResult:
    t_release: float   # time from de-energization to stop-release, s
    t_close: float      # time from de-energization to valve fully closed, s
    sol_hold: object    # OdeResult for the held-open phase (t=0 at de-energization)
    sol_stroke: object  # OdeResult for the released stroke (t=0 at the release instant)
    bounce: object = None  # BounceResult when seat bounce is enabled, else None
                           # t_close keeps meaning FIRST seat contact either way


@dataclass
class BounceResult:
    """Result of a simulation run with seat bounce enabled.

    Clock base is NOT the same on the two calling paths. On the opening path
    (simulate_opening) there is only one clock for the whole call: t=0 is
    energization, and t_first/t_settle/segments are all relative to it. On
    the closing path (simulate_closing) this object wraps only Phase B (the
    released stroke), whose own integration restarts at t=0 = the release
    instant -- so a closing BounceResult's t_first/t_settle are
    RELEASE-relative, whereas ClosingResult.t_close and .t_release are
    DE-ENERGIZATION-relative. The two differ by exactly t_release: a caller
    comparing e.g. bounce.t_settle against t_close directly, without adding
    t_release first, will be off by t_release.

    t_settle is the time of the LAST recorded impact, not necessarily the
    time bouncing actually decayed below v_min -- those coincide when the
    train ends normally (the next rebound would be <= v_min), but not when
    the n_max guard cuts the train off first, in which case t_settle is
    simply wherever the loop stopped, possibly still bouncing above v_min."""
    t_first: float     # first contact with the travel limit, s (clock base
                       # is call-dependent -- see docstring)
    t_settle: float    # time of the LAST recorded impact, s, same clock
                       # base as t_first; NOT necessarily when |v| actually
                       # fell below v_min -- see docstring
    n_bounce: int      # number of impacts, >= 1
    v_impacts: list    # impact speed of each impact, m/s
    segments: list     # OdeResult per flight segment, in order


def electrical_di_dt(i, gap, v, params, R=None, V=None):
    if R is None:
        R = params.R_coil_20C
    if V is None:
        V = params.V_bus
    L = magnetics.inductance(gap, params)
    dL_dx = -magnetics.dinductance_dgap(gap, params)  # chain rule: gap = g0 - x, d(gap)/dx = -1
    back_emf_term = i * dL_dx * v
    return (V - i * R - back_emf_term) / L


def spring_force(x, params):
    return params.F_preload + params.k_spring * x


def net_mechanical_force(x, v, F_mag, F_pressure, F_flow, params):
    return F_mag - spring_force(x, params) - F_pressure - F_flow - params.damping_coeff * v


def mechanical_dv_dt(x, v, F_mag, params):
    F_pressure = params.delta_P * params.A_seat
    return net_mechanical_force(x, v, F_mag, F_pressure, 0.0, params) / params.m_arm


def _fluid_forces(x, params, medium, cond):
    """Coupled-mode forces: static pressure term from conditions (supersedes
    params.delta_P) plus the medium-appropriate flow force."""
    F_pressure = (cond.P_up - cond.P_down) * params.A_seat
    if isinstance(medium, fluid_lib.GasProperties):
        F_flow = fluid_lib.flow_force_gas(x, params, medium, cond)
    else:
        F_flow = fluid_lib.flow_force_liquid(x, params, cond)
    return F_pressure, F_flow


def coupled_rhs(t, state, params, medium=None, cond=None, R_coil_eff=None):
    i, x, v = state
    gap = params.g0 - x
    di_dt = electrical_di_dt(i, gap, v, params, R=R_coil_eff)
    F_mag = magnetics.magnetic_force_closing(gap, i, params)
    if medium is None:
        dv_dt = mechanical_dv_dt(x, v, F_mag, params)
    else:
        F_pressure, F_flow = _fluid_forces(x, params, medium, cond)
        dv_dt = net_mechanical_force(x, v, F_mag, F_pressure, F_flow, params) / params.m_arm

    if x <= 0.0 and dv_dt <= 0.0:
        return [di_dt, 0.0, 0.0]  # held against the rest-position mechanical stop

    return [di_dt, v, dv_dt]


# Restart offset for the bounce-train contact-plane nudge below. Needs to be
# small enough to be a numerical formality rather than a physically
# significant displacement, and far below the shallowest real flight
# excursion this loop ever integrates -- the shallowest measured is ~12 nm
# (PCTFE opening bounce, last flight before settling), four orders of
# magnitude above this value. Coincidentally equal to the ODE solver's atol;
# that is not why this value was chosen, and the two are free to diverge.
_CONTACT_EPS = 1e-12


def _bounce_train(first_sol, rhs, event, t_end, x_contact, seat, v_min, n_max,
                  solve_kwargs):
    """Shared bounce loop for both opening and closing.

    Reflects velocity at each contact (v -> -e*v) and re-integrates until the
    next rebound would fall to or below v_min -- the same termination rule
    impact.bounce_sequence uses, so the two stay comparable. The restart
    position is nudged just off the contact plane (by _CONTACT_EPS) because
    solve_ivp would otherwise see the event function sitting exactly at zero
    at t0.
    """
    t_first = first_sol.t_events[0][0]
    segments = [first_sol]
    v_impacts = [abs(first_sol.y_events[0][0][2])]
    t_now = t_first
    state = list(first_sol.y_events[0][0])
    while len(v_impacts) < n_max:
        v_rebound = seat.e_restitution * v_impacts[-1]
        if v_rebound <= v_min:
            break
        direction = -1.0 if x_contact > 0.0 else 1.0
        state[1] = x_contact + direction * _CONTACT_EPS
        state[2] = direction * v_rebound
        seg = solve_ivp(rhs, (t_now, t_end), state, events=event,
                        **solve_kwargs)
        segments.append(seg)
        if not seg.t_events[0].size:
            break
        t_now = seg.t_events[0][0]
        state = list(seg.y_events[0][0])
        v_impacts.append(abs(state[2]))
    return BounceResult(t_first=t_first, t_settle=t_now,
                        n_bounce=len(v_impacts), v_impacts=v_impacts,
                        segments=segments)


def simulate_opening(params, t_max=0.05, medium=None, cond=None, T_coil=None,
                     seat=None, v_min=None, n_max=None):
    """Simulate the opening transient.

    seat=None (default) keeps the historical behaviour exactly: the armature
    is hard-clamped at the travel stop and the raw OdeResult is returned.
    Passing a sealing.SeatPair switches on bounce and returns a BounceResult
    instead -- the fourth orthogonal mode switch alongside medium/cond and
    T_coil.
    """
    if medium is not None and cond is None:
        raise ValueError("coupled mode requires explicit FlowConditions (cond=...)")
    # isothermal per event: R evaluated once, constant for the whole simulation
    R_coil_eff = None if T_coil is None else thermal.R_coil(T_coil, params)

    def stroke_event(t, state):
        return state[1] - params.x_stroke

    stroke_event.terminal = True
    stroke_event.direction = 1

    rhs = lambda t, y: coupled_rhs(t, y, params, medium=medium, cond=cond,
                                   R_coil_eff=R_coil_eff)
    solve_kwargs = dict(max_step=t_max / 5000, rtol=1e-9, atol=1e-12,
                        dense_output=True)
    sol = solve_ivp(rhs, (0.0, t_max), [0.0, 0.0, 0.0], events=stroke_event,
                    **solve_kwargs)
    if seat is None:
        return sol
    if not sol.t_events[0].size:
        raise ValueError(
            "valve does not reach the travel stop; bounce is undefined")
    v_min = impact_lib.V_MIN_DEFAULT if v_min is None else v_min
    n_max = impact_lib.N_MAX_DEFAULT if n_max is None else n_max
    return _bounce_train(sol, rhs, stroke_event, t_max, params.x_stroke,
                         seat, v_min, n_max, solve_kwargs)


def _release_dv(i, params):
    """Net mechanical acceleration at the held-open position (x=x_stroke,
    v=0) for the given coil current. Negative once the magnetic force drops
    below spring+pressure resistance — that is the stop-release condition."""
    gap = params.g0 - params.x_stroke
    F_mag = magnetics.magnetic_force_closing(gap, i, params)
    return mechanical_dv_dt(params.x_stroke, 0.0, F_mag, params)


def _closing_phase_a_rhs_clamped(t, y, params, circuit):
    i = y[0]
    gap = params.g0 - params.x_stroke
    V = flyback_lib.clamp_voltage(circuit, i)
    di_dt = electrical_di_dt(i, gap, 0.0, params, V=V) if i > 0.0 else 0.0
    return [di_dt]


def _closing_phase_a_rhs_rc(t, y, params, circuit):
    i, V_c = y
    gap = params.g0 - params.x_stroke
    V = flyback_lib.rc_terminal_voltage(circuit, i, V_c)
    di_dt = electrical_di_dt(i, gap, 0.0, params, V=V)
    dVc_dt = i / circuit.C_snub
    return [di_dt, dVc_dt]


def _closing_stroke_rhs_clamped(t, state, params, circuit):
    i, x, v = state
    gap = params.g0 - x
    V = flyback_lib.clamp_voltage(circuit, i)
    di_dt = electrical_di_dt(i, gap, v, params, V=V) if i > 0.0 else 0.0
    F_mag = magnetics.magnetic_force_closing(gap, i, params)
    dv_dt = mechanical_dv_dt(x, v, F_mag, params)
    if x <= 0.0 and dv_dt <= 0.0:
        return [di_dt, 0.0, 0.0]  # held against the fully-closed seat
    return [di_dt, v, dv_dt]


def _closing_stroke_rhs_rc(t, state, params, circuit):
    i, x, v, V_c = state
    gap = params.g0 - x
    V = flyback_lib.rc_terminal_voltage(circuit, i, V_c)
    di_dt = electrical_di_dt(i, gap, v, params, V=V)
    dVc_dt = i / circuit.C_snub
    F_mag = magnetics.magnetic_force_closing(gap, i, params)
    dv_dt = mechanical_dv_dt(x, v, F_mag, params)
    if x <= 0.0 and dv_dt <= 0.0:
        return [di_dt, 0.0, 0.0, dVc_dt]  # held against the fully-closed seat
    return [di_dt, v, dv_dt, dVc_dt]


def simulate_closing(params, circuit, t_max_hold=0.3, t_max_stroke=0.05,
                     seat=None, v_min=None, n_max=None):
    """Simulate the closing (de-energized) transient: from full-open steady
    holding current, through the flyback path's electrical decay, to the
    armature returning to x=0 (valve fully closed). Dry mode only — no
    fluid/thermal coupling.

    Two-phase integration: Phase A holds x=x_stroke, v=0 fixed (the armature
    is pinned against the fully-open mechanical stop while the holding
    current — ~2.4 kN at the baseline case — vastly exceeds spring+pressure
    resistance, ~28 N) and integrates only the electrical decay until the
    magnetic force drops enough to release the stop. Phase B then integrates
    the full state from the release point to x=0. A single combined ODE
    with inline clamps was tried first and chatters at the release boundary
    (adaptive step size collapses); this two-phase split avoids that.

    seat=None keeps the historical hard-clamp behaviour. Passing a
    sealing.SeatPair additionally integrates the bounce train after first
    seat contact and attaches it as .bounce; t_close and t_release keep
    their original meanings so P3's results are unaffected.
    """
    I_ss = params.V_bus / params.R_coil_20C
    is_rc = isinstance(circuit, flyback_lib.RCFlyback)

    def release_event(t, y):
        return _release_dv(y[0], params)

    release_event.terminal = True
    release_event.direction = -1

    if is_rc:
        sol_hold = solve_ivp(
            lambda t, y: _closing_phase_a_rhs_rc(t, y, params, circuit),
            (0.0, t_max_hold), [I_ss, 0.0], events=release_event,
            max_step=t_max_hold / 5000, rtol=1e-9, atol=1e-12, dense_output=True,
        )
    else:
        sol_hold = solve_ivp(
            lambda t, y: _closing_phase_a_rhs_clamped(t, y, params, circuit),
            (0.0, t_max_hold), [I_ss], events=release_event,
            max_step=t_max_hold / 5000, rtol=1e-9, atol=1e-12, dense_output=True,
        )

    if not sol_hold.t_events[0].size:
        raise ValueError(
            f"circuit {circuit!r} does not release the armature within "
            f"t_max_hold={t_max_hold}s"
        )

    t_release = sol_hold.t_events[0][0]
    y_release = sol_hold.y_events[0][0]

    def closed_event(t, state):
        return state[1]

    closed_event.terminal = True
    closed_event.direction = -1

    if is_rc:
        i_r, Vc_r = y_release
        y0_stroke = [i_r, params.x_stroke, 0.0, Vc_r]
        rhs_stroke = lambda t, y: _closing_stroke_rhs_rc(t, y, params, circuit)
    else:
        i_r = y_release[0]
        y0_stroke = [i_r, params.x_stroke, 0.0]
        rhs_stroke = lambda t, y: _closing_stroke_rhs_clamped(t, y, params, circuit)

    sol_stroke = solve_ivp(
        rhs_stroke, (0.0, t_max_stroke), y0_stroke, events=closed_event,
        max_step=t_max_stroke / 5000, rtol=1e-9, atol=1e-12, dense_output=True,
    )

    if not sol_stroke.t_events[0].size:
        raise ValueError(
            f"circuit {circuit!r} released at t={t_release}s but did not "
            f"reach x=0 within t_max_stroke={t_max_stroke}s"
        )

    t_close = t_release + sol_stroke.t_events[0][0]
    result = ClosingResult(t_release=t_release, t_close=t_close,
                           sol_hold=sol_hold, sol_stroke=sol_stroke)
    if seat is None:
        return result

    v_min = impact_lib.V_MIN_DEFAULT if v_min is None else v_min
    n_max = impact_lib.N_MAX_DEFAULT if n_max is None else n_max
    solve_kwargs = dict(max_step=t_max_stroke / 5000, rtol=1e-9, atol=1e-12,
                        dense_output=True)
    result.bounce = _bounce_train(sol_stroke, rhs_stroke, closed_event,
                                  t_max_stroke, 0.0, seat, v_min, n_max,
                                  solve_kwargs)
    return result
