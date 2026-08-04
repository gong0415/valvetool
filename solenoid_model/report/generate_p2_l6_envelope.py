"""P2 L6 report: voltage-temperature operating envelope (spec §3 L6),
response-characteristic-4 table (t_open / V_pull-in drift at -40/+70 degC),
and the LN2 77 K out-of-envelope case.

Full-resolution run (defaults) takes a few minutes: each envelope point is
~10 dynamic simulations and sub-threshold probes run the full t_max. The
test suite calls generate() with a coarse grid instead.
"""
import sys
from dataclasses import replace
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from solenoid_model import limits, magnetics, thermal
from solenoid_model.baseline import BASELINE_PARAMS
from solenoid_model.dynamics import simulate_opening

V_BUS_MIN, V_BUS_MAX = 22.0, 36.0  # spec §3 L6: typical satellite bus band


def generate(n_points=12, tol=0.1, t_max=0.03, out_path=None):
    p = BASELINE_PARAMS
    T_grid = list(np.linspace(233.15, 343.15, n_points))

    # bare envelope: coil at ambient temperature
    v_pi_bare = [limits.find_pullin_voltage(p, T_coil=T, t_max=t_max, tol=tol)
                 for T in T_grid]

    # self-heated envelope: continuous hold at V_bus, constant-current
    # conservative assumption I = V_bus / R(T_amb) (same convention as GUI)
    T_selfheat, v_pi_selfheat = [], []
    for T_amb in T_grid:
        I_hs = p.V_bus / thermal.R_coil(T_amb, p)
        T_hot = thermal.equilibrium_temp(I_hs, T_amb, p)
        T_selfheat.append(T_hot)
        v_pi_selfheat.append(
            limits.find_pullin_voltage(p, T_coil=T_hot, t_max=t_max, tol=tol))

    # LN2 77 K: outside the -40..+70 degC envelope, annotated separately
    R77 = thermal.R_coil(77.0, p)
    v_pi_ln2 = limits.find_pullin_voltage(p, T_coil=77.0, t_max=t_max, tol=tol)
    I_ss_ln2 = p.V_bus / R77
    ln2_saturated = magnetics.saturation_check(p.g0 - p.x_stroke, I_ss_ln2, p)

    # response characteristic 4: t_open drift over temperature and voltage
    resp4 = []
    for T in (233.15, 293.15, 343.15):
        for V in (22.0, 28.0):
            sol = simulate_opening(replace(p, V_bus=V), T_coil=T)
            t_open = sol.t_events[0][0] if sol.t_events[0].size else None
            resp4.append((T, V, t_open))

    _plot(T_grid, v_pi_bare, v_pi_selfheat, v_pi_ln2, ln2_saturated, out_path)
    _print_summary(T_grid, v_pi_bare, T_selfheat, v_pi_selfheat,
                   R77, v_pi_ln2, I_ss_ln2, ln2_saturated, resp4)

    return {
        "T_grid": T_grid, "v_pi_bare": v_pi_bare,
        "T_selfheat": T_selfheat, "v_pi_selfheat": v_pi_selfheat,
        "v_pi_ln2": v_pi_ln2, "ln2_saturated": ln2_saturated,
        "resp4": resp4,
    }


def _plot(T_grid, v_pi_bare, v_pi_selfheat, v_pi_ln2, ln2_saturated,
          out_path=None):
    fig, ax = plt.subplots()
    ax.axhspan(V_BUS_MIN, V_BUS_MAX, alpha=0.15, color="tab:green",
               label=f"bus band {V_BUS_MIN:.0f}-{V_BUS_MAX:.0f} V")
    ax.plot(T_grid, v_pi_bare, marker="o", label="V_pull-in (coil at ambient)")
    ax.plot(T_grid, v_pi_selfheat, marker="s",
            label="V_pull-in (continuous hold, self-heated)")
    suffix = ", B_sat!" if ln2_saturated else ""
    ax.annotate(f"LN2 77 K: {v_pi_ln2:.1f} V (out of envelope{suffix})",
                xy=(T_grid[0], v_pi_bare[0]),
                xytext=(0.02, 0.03), textcoords="axes fraction", fontsize=8)
    ax.set_xlabel("coil temperature (K)")
    ax.set_ylabel("pull-in voltage (V)")
    ax.set_title("L6: voltage-temperature operating envelope")
    ax.legend(fontsize=8)
    if out_path is None:
        out_path = Path(__file__).parent / "l6_envelope_v_t.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _print_summary(T_grid, v_pi_bare, T_selfheat, v_pi_selfheat,
                   R77, v_pi_ln2, I_ss_ln2, ln2_saturated, resp4):
    print("== L6 envelope ==")
    for T, vb, Th, vs in zip(T_grid, v_pi_bare, T_selfheat, v_pi_selfheat):
        print(f"  T_amb={T:7.2f} K  V_pi={vb:6.2f} V | "
              f"self-heated T={Th:7.2f} K  V_pi={vs:6.2f} V")
    print(f"  margin at hot soak: {V_BUS_MIN - v_pi_selfheat[-1]:+.2f} V "
          f"vs bus min {V_BUS_MIN} V")
    print("== LN2 77 K (outside envelope) ==")
    print(f"  R(77 K)={R77:.2f} ohm  V_pi={v_pi_ln2:.2f} V  "
          f"I_ss={I_ss_ln2:.2f} A")
    if ln2_saturated:
        print("  WARNING: B > B_sat at steady state — linear magnetics "
              "invalid here; result is a flag, not a prediction")
    print("== response characteristic 4: t_open (ms) ==")
    for T, V, t_open in resp4:
        t_ms = "no pull-in" if t_open is None else f"{t_open * 1e3:.3f}"
        print(f"  T={T:7.2f} K  V={V:4.1f} V  t_open={t_ms}")


if __name__ == "__main__":
    generate()
