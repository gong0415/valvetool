"""Standalone fluid library (P1): valve flow physics as pure functions.

Never imports dynamics/magnetics — coupling into the armature ODE lives in
dynamics.py.
"""
import math
from dataclasses import dataclass


@dataclass
class GasProperties:
    gamma: float        # specific heat ratio
    R_specific: float   # specific gas constant, J/(kg K)


@dataclass
class LiquidProperties:
    rho: float          # density, kg/m^3 (at stated temperature)
    c_sound: float      # speed of sound, m/s (water hammer)
    P_vap: float        # vapor pressure, Pa (flashing flag)


@dataclass
class FlowConditions:
    P_up: float         # upstream stagnation pressure, Pa
    P_down: float       # downstream pressure, Pa
    T0: float           # upstream stagnation temperature, K (gas branch)


# Preset property constants (engineering tables, ~3 sig figs; sources in docs/spec/PARAMS.md)
N2 = GasProperties(gamma=1.4, R_specific=296.8)
XE = GasProperties(gamma=1.667, R_specific=63.33)
HE = GasProperties(gamma=1.667, R_specific=2077.1)

WATER_20C = LiquidProperties(rho=998.2, c_sound=1482.0, P_vap=2339.0)
LN2_77K = LiquidProperties(rho=806.1, c_sound=850.0, P_vap=101325.0)  # at normal boiling point

_KV_RHO_WATER = 1000.0   # kg/m^3, Kv definition condition
_KV_DP_REF = 1.0e5       # Pa (1 bar), Kv definition condition
_CV_PER_KV = 1.156       # standard Cv/Kv conversion constant


def effective_area(x, params):
    """Flat-poppet curtain-area model; negative lift clamps to zero (closed).

    Clamping (not raising) keeps ODE solvers safe when they probe x<0.
    """
    x_eff = max(x, 0.0)
    a_curtain = math.pi * params.D_seat_bore * x_eff
    a_bore = math.pi * params.D_seat_bore ** 2 / 4
    return min(a_curtain, a_bore)


def critical_pressure_ratio(gas):
    g = gas.gamma
    return (2 / (g + 1)) ** (g / (g - 1))


def mdot_gas(x, params, gas, cond):
    """Ideal-gas mass flow through the seat: choked below r_crit, isentropic subsonic above.

    Guards (documented physical clamps, no backflow modeling): P_up<=0 or
    reversed/zero pressure difference returns 0; negative P_down treated as vacuum.
    """
    if cond.P_up <= 0.0:
        return 0.0
    a_eff = effective_area(x, params)
    g, R, T0, P_up = gas.gamma, gas.R_specific, cond.T0, cond.P_up
    r = max(cond.P_down, 0.0) / P_up
    if r >= 1.0:
        return 0.0
    if r <= critical_pressure_ratio(gas):
        flux = math.sqrt(g / (R * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
    else:
        flux = math.sqrt(2 * g / (R * T0 * (g - 1)) * (r ** (2 / g) - r ** ((g + 1) / g)))
    return params.C_d * a_eff * P_up * flux


def mdot_liquid(x, params, liquid, cond):
    """Incompressible Bernoulli orifice flow; reversed dP returns 0 (no backflow modeling)."""
    dp = cond.P_up - max(cond.P_down, 0.0)
    if dp <= 0.0:
        return 0.0
    return params.C_d * effective_area(x, params) * math.sqrt(2 * liquid.rho * dp)


def water_hammer_dp(liquid, delta_v):
    """Joukowsky surge pressure for an instantaneous velocity change."""
    return liquid.rho * liquid.c_sound * delta_v


def flashing_risk(liquid, cond):
    """Conservative flag: downstream pressure at/below vapor pressure means flashing.

    Vena-contracta pressure dips below P_down, so absence of this flag does not
    guarantee no cavitation (refined criteria are empirical; see PARAMS.md).
    """
    return cond.P_down <= liquid.P_vap


def kv_from_geometry(x, params):
    """Kv (m^3/h of water at 1 bar) equivalent of the current opening."""
    velocity = math.sqrt(2 * _KV_DP_REF / _KV_RHO_WATER)
    return params.C_d * effective_area(x, params) * velocity * 3600


def cv_from_geometry(x, params):
    return _CV_PER_KV * kv_from_geometry(x, params)


def flow_force_liquid(x, params, cond):
    """Momentum-flux reaction on the poppet, magnitude: mdot*v_jet = 2*C_d*A_eff*dP.

    Direction convention: tends to CLOSE the valve at small lift. Reversed dP
    returns 0 (no backflow modeling). Omits jet-angle cos(theta) correction
    (sharp seat ~69 deg, empirical) — overestimates the closing force.
    """
    dp = cond.P_up - max(cond.P_down, 0.0)
    if dp <= 0.0:
        return 0.0
    return 2 * params.C_d * effective_area(x, params) * dp


def flow_force_gas(x, params, gas, cond):
    """Gas flow force on the poppet: jet momentum + throat pressure imbalance.

    Choked (r <= r_crit): F = mdot*v_throat + (P_throat - P_down)*A_eff, with sonic
    throat conditions T_t = T0*2/(gamma+1), v_t = sqrt(gamma*R*T_t), P_t = P_up*r_crit.
    Subsonic: exit expands to P_down, so F = mdot*v_exit only. The two branches are
    analytically equal at r_crit. Direction convention: tends to CLOSE the valve.
    Omits jet-angle cos(theta) correction (empirical) — overestimates closing force.
    Same guards as mdot_gas (no backflow modeling).
    """
    if cond.P_up <= 0.0:
        return 0.0
    g, R, T0 = gas.gamma, gas.R_specific, cond.T0
    r = max(cond.P_down, 0.0) / cond.P_up
    if r >= 1.0:
        return 0.0
    a_eff = effective_area(x, params)
    m = mdot_gas(x, params, gas, cond)
    if r <= critical_pressure_ratio(gas):
        t_throat = T0 * 2 / (g + 1)
        v_throat = math.sqrt(g * R * t_throat)
        p_throat = cond.P_up * critical_pressure_ratio(gas)
        return m * v_throat + (p_throat - cond.P_down) * a_eff
    v_exit = math.sqrt(2 * g * R * T0 / (g - 1) * (1 - r ** ((g - 1) / g)))
    return m * v_exit
