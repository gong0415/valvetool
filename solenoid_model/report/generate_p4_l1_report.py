"""P4 report: L1 response-time limit (spec §3 L1). Sweeps t_open over bus
voltage, turns, and stroke, overlaying the analytic lower bound and the
spec's own tau_e + sqrt(2mx/F) decomposition.

Two headline findings: the analytic bound holds in every case at ~0.80-0.91
tightness, while the spec's decomposition is not a bound in either direction;
and the turns sweep runs in opposite directions depending on whether coil
resistance follows the winding window or is held fixed -- the fixed-R curve
is the hidden degree of freedom the winding-window model closed, drawn here
for contrast.
"""
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from solenoid_model import limits, winding
from solenoid_model.baseline import BASELINE_PARAMS

V_VALUES = [18.0, 22.0, 28.0, 36.0, 50.0]
N_VALUES = [1000.0, 1500.0, 2000.0, 3000.0, 4000.0]
X_VALUES = [1.0e-4, 2.0e-4, 3.0e-4, 4.0e-4, 5.0e-4]
RESIDUAL_GAP = 0.5e-4   # held constant so the stroke sweep moves one DOF


def generate(out_path=None):
    p = BASELINE_PARAMS

    voltage = limits.l1_sweep(
        [(f"{V:g}", replace(p, V_bus=V), p.R_coil_20C) for V in V_VALUES])
    turns_window = limits.l1_sweep(
        [(f"{N:g}", replace(p, N_turns=N),
          winding.coil_resistance(N, p.A_winding, p.k_fill, p.l_turn_mean))
         for N in N_VALUES])
    turns_fixed_R = limits.l1_sweep(
        [(f"{N:g}", replace(p, N_turns=N), p.R_coil_20C) for N in N_VALUES])
    stroke = limits.l1_sweep(
        [(f"{x*1e3:g}", replace(p, x_stroke=x, g0=x + RESIDUAL_GAP),
          p.R_coil_20C) for x in X_VALUES])

    sweeps = {"voltage": voltage, "turns_window": turns_window,
              "turns_fixed_R": turns_fixed_R, "stroke": stroke}
    bound_ratio_range = _ratio_range(sweeps, "t_bound")
    naive_ratio_range = _ratio_range(sweeps, "t_naive")

    _plot(voltage, turns_window, turns_fixed_R, stroke, out_path)
    _print_summary(p, sweeps, bound_ratio_range, naive_ratio_range)

    out = dict(sweeps)
    out["bound_ratio_range"] = bound_ratio_range
    out["naive_ratio_range"] = naive_ratio_range
    return out


def _ratio_range(sweeps, key):
    ratios = [v / t for s in sweeps.values()
              for v, t in zip(s[key], s["t_open"])
              if v is not None and t is not None]
    return (min(ratios), max(ratios))


def _panel(ax, sweep, xs, title, xlabel, label=None, style="o-"):
    xs_ok = [x for x, ok in zip(xs, sweep["pulls_in"]) if ok]
    ax.plot(xs_ok, [t * 1e3 for t in sweep["t_open"] if t is not None],
            style, label=label or "simulated t_open")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)


def _plot(voltage, turns_window, turns_fixed_R, stroke, out_path=None):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    ax = axes[0]
    xs = [float(s) for s in voltage["label"]]
    _panel(ax, voltage, xs, "t_open vs bus voltage", "V_bus (V)")
    ax.plot([x for x, ok in zip(xs, voltage["pulls_in"]) if ok],
            [t * 1e3 for t in voltage["t_bound"] if t is not None], "s--",
            label="analytic lower bound")
    ax.plot([x for x, ok in zip(xs, voltage["pulls_in"]) if ok],
            [t * 1e3 for t in voltage["t_naive"] if t is not None], "^:",
            label="spec naive estimate (not a bound)")
    ax.set_ylabel("t_open (ms)")
    ax.legend(fontsize=7)

    ax = axes[1]
    xs = [float(s) for s in turns_window["label"]]
    _panel(ax, turns_window, xs, "t_open vs turns", "N_turns",
           label="winding-consistent R (physical)")
    _panel(ax, turns_fixed_R, xs, "t_open vs turns", "N_turns",
           label="fixed R (inconsistent)", style="x--")
    ax.legend(fontsize=7)

    ax = axes[2]
    xs = [float(s) for s in stroke["label"]]
    _panel(ax, stroke, xs, "t_open vs stroke", "x_stroke (mm)")
    ax.plot([x for x, ok in zip(xs, stroke["pulls_in"]) if ok],
            [t * 1e3 for t in stroke["t_bound"] if t is not None], "s--",
            label="analytic lower bound")
    ax.legend(fontsize=7)

    fig.suptitle("L1: response-time limit sweeps", fontsize=11)
    fig.tight_layout()
    if out_path is None:
        out_path = Path(__file__).parent / "p4_l1_sweeps.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _print_summary(p, sweeps, bound_ratio_range, naive_ratio_range):
    print("== L1: response-time limit ==")
    print(f"  tau_e at rest gap      = {limits.electrical_time_constant(p)*1e3:.4f} ms")
    print(f"  tau_e at pulled-in gap = "
          f"{limits.electrical_time_constant(p, gap=p.g0-p.x_stroke)*1e3:.4f} ms")
    print(f"  motion-threshold current = {limits.motion_threshold_current(p):.5f} A")
    print(f"  analytic bound / t_open : {bound_ratio_range[0]:.3f} .. "
          f"{bound_ratio_range[1]:.3f}  (all < 1.0 -> a real bound)")
    print(f"  spec naive  / t_open    : {naive_ratio_range[0]:.3f} .. "
          f"{naive_ratio_range[1]:.3f}  (straddles 1.0 -> NOT a bound)")
    tw, tf = sweeps["turns_window"], sweeps["turns_fixed_R"]
    print("  turns sweep, window-consistent vs fixed R:")
    for label, a, b in zip(tw["label"], tw["t_open"], tf["t_open"]):
        fa = "no pull-in" if a is None else f"{a*1e3:8.4f} ms"
        fb = "no pull-in" if b is None else f"{b*1e3:8.4f} ms"
        # The fixed-R variant is exactly winding.coil_resistance_consistency's
        # target case: N_turns changes but R_coil_20C stays at the baseline
        # declared value, so it drifts from the winding-window geometric R as
        # N grows. Surface that divergence alongside the fixed-R timing so
        # the "inconsistent" label above is backed by a number.
        _, _, rel_error = winding.coil_resistance_consistency(
            replace(p, N_turns=float(label)))
        print(f"    N={label:>5}: window {fa}   fixed-R {fb}"
              f"   (fixed-R declared vs geometric: rel_err={rel_error:.3f})")


if __name__ == "__main__":
    generate()
