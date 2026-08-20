"""Design cases: named (ValveParams + operating conditions) bundles.

`baseline.BASELINE_PARAMS` stays the reference validation case. Cases added
here are *design results* — each one records the requirement it was sized
against, so the numbers can be re-derived rather than trusted blindly.

Operating conditions travel with the case because a design point is only
meaningful at the duty it was sized for: N2_25BAR's bore is reverse-solved
from a choked-flow mass-flow target at P_up=25 bar, and reading it at a
different upstream pressure silently invalidates that bore.
"""
import math
from dataclasses import dataclass

from solenoid_model import fluid
from solenoid_model.baseline import BASELINE_PARAMS
from solenoid_model.params import ValveParams
from solenoid_model.winding import coil_resistance


@dataclass
class DesignCase:
    """A ValveParams plus the duty point it was sized against."""
    name: str
    params: ValveParams
    medium: object              # fluid.GasProperties | fluid.LiquidProperties
    cond: fluid.FlowConditions
    T_coil: float               # K, isothermal coil assumption
    notes: str = ""


# --- N2 25 bar case ---------------------------------------------------------
# Requirement (user-supplied): N2, P_up=25 bar, P_down~0 (vacuum), T0=25 C,
# mdot=1.28 g/s, 28 V bus, normally-closed direct-acting, minimize hold power.
#
# Sizing chain (each step verified against the model, not hand algebra):
#  1. r = P_down/P_up = 0 << r_crit(N2)=0.528  -> always choked. Mass flow then
#     depends only on A_eff and P_up, so the bore alone sets flow:
#     A_req = mdot / (C_d * P_up * flux), flux = sqrt(g/(R*T0))*(2/(g+1))^((g+1)/(2(g-1)))
#     -> A_req = 0.278 mm^2 -> D_SEAT_BORE = 0.60 mm (gives 1.302 g/s = 101.7% of target).
#  2. Curtain area pi*D*x saturates at the bore area when x = D/4, so
#     x_stroke = 0.15 mm; further travel buys no flow.
#  3. Static seat load is only 25 bar * A_seat = 0.71 N, which is what makes a
#     direct-acting armature viable here at all.
#  4. For a fixed winding window, coil_resistance gives R ~ N^2, so hold power
#     P = V^2/R = rho*l_turn*(N*i)^2/(A_winding*k_fill) is independent of N at
#     fixed MMF. Turns were then swept against the coupled ODE for the knee of
#     the power-vs-speed curve: N=6000 -> 0.76 W at t_open=5.35 ms with 2.57x
#     pull-in force margin (2.40x with the coil at its 28.8 C self-heat
#     equilibrium). N=10000 (0.27 W) fails to pull in at all.
_N2_OD, _N2_ID, _N2_H = 13.0e-3, 6.0e-3, 10.0e-3   # coil outer/inner dia, height, m
_N2_K_FILL = 0.5
_N2_L_TURN = 2 * math.pi * ((_N2_OD + _N2_ID) / 4)
_N2_A_WINDING = ((_N2_OD - _N2_ID) / 2) * _N2_H
_N2_TURNS = 6000.0

N2_25BAR_PARAMS = ValveParams(
    # Electrical — R derived from the winding window so the two stay consistent
    # (the app flags a >5% mismatch between R_coil_20C and the window-implied R).
    V_bus=28.0,
    R_coil_20C=coil_resistance(_N2_TURNS, _N2_A_WINDING, _N2_K_FILL, _N2_L_TURN),
    N_turns=_N2_TURNS,

    # Magnetic circuit — 5.05 mm pole dia; B=0.97 T at the closed gap keeps
    # margin to B_sat where force is scarcest.
    A_gap=20.0e-6,
    g0=0.20e-3,
    x_stroke=0.15e-3,
    l_core=0.045,
    mu_r_core=4000.0,
    B_sat=2.15,          # soft iron

    # Mechanical — preload sized ~3x the 0.71 N static load for seat sealing.
    m_arm=0.8e-3,
    k_spring=4000.0,
    F_preload=2.2,

    # Simplified loads — delta_P/A_seat only bite in dry-run mode; with fluid
    # coupling on, the static load is taken from (P_up-P_down)*A_seat instead.
    delta_P=2.5e6,
    A_seat=0.2827e-6,    # = bore area, poppet seals on the bore
    damping_coeff=0.05,  # EMPIRICAL placeholder, as in baseline

    # Flow path — the sized result
    D_seat_bore=0.60e-3,
    C_d=0.8,             # EMPIRICAL

    # Winding window — co-consistent with R_coil_20C above
    l_turn_mean=_N2_L_TURN,
    k_fill=_N2_K_FILL,
    A_winding=_N2_A_WINDING,
)

N2_25BAR = DesignCase(
    name="N₂ 25bar 1.28g/s",
    params=N2_25BAR_PARAMS,
    medium=fluid.N2,
    cond=fluid.FlowConditions(P_up=25.0e5, P_down=0.0, T0=298.15),
    T_coil=298.15,
    notes=("常閉直動式，以最小保持功耗為目標。0.76 W / t_open 5.35 ms / "
           "吸合裕度 2.57×。孔徑 0.60 mm 由 1.28 g/s 壅塞流反推；"
           "接管 1.88 mm 為孔徑面積的 21 倍，不構成節流。"),
)

BASELINE = DesignCase(
    name="基準案例 (baseline)",
    params=BASELINE_PARAMS,
    medium=fluid.N2,
    cond=fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0),
    T_coil=293.15,
    notes="原始驗證案例，非特定商用產品；見 baseline.py",
)

CASES = {c.name: c for c in (BASELINE, N2_25BAR)}
