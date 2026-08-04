"""P4 report: L2 power-vs-holding-force Pareto (spec §3 L2). The winding
window fixes the F/P slope; the B_sat force ceiling caps the front. Headline
finding: the baseline coil runs deep in saturation at hold, so extra hold
power buys zero extra force (the spike-and-hold argument, quantified).
"""
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from solenoid_model import limits
from solenoid_model.baseline import BASELINE_PARAMS

# turns sweep spanning the saturated baseline (N=2000) up through the
# fine-wire, low-power region where the ceiling stops binding (~14815 turns)
N_VALUES = list(np.linspace(2000.0, 24000.0, 45))


def generate(out_path=None):
    p = BASELINE_PARAMS
    front = limits.l2_pareto_front(p, N_VALUES)
    F_sat = limits.saturation_force_ceiling(p)

    # baseline (N=2000) hold power and saturation state
    P_baseline = front["P_hold"][0]
    baseline_saturated = bool(front["F_hold"][0] >= F_sat - 1e-9)

    # crossover power: smallest hold power in the sweep that is still capped
    # (the ceiling binds at and below this power / above this N)
    capped_powers = [P for P, F in zip(front["P_hold"], front["F_hold"])
                     if F >= F_sat - 1e-9]
    P_cross = min(capped_powers) if capped_powers else float("nan")

    _plot(front, F_sat, P_baseline, out_path)
    _print_summary(front, F_sat, P_baseline, P_cross, baseline_saturated)

    return {"front": front, "F_sat": F_sat, "P_cross": P_cross,
            "baseline_saturated": baseline_saturated}


def _plot(front, F_sat, P_baseline, out_path=None):
    fig, ax = plt.subplots()
    ax.plot(front["P_hold"], front["F_linear"], "--", color="tab:gray",
            label="linear magnetics (fictitious above B_sat)")
    ax.plot(front["P_hold"], front["F_hold"], marker="o", ms=3,
            color="tab:blue", label="real force (B_sat-capped front)")
    ax.axhline(F_sat, color="tab:red", lw=1,
               label=f"B_sat ceiling {F_sat:.1f} N")
    ax.axvline(P_baseline, color="tab:green", lw=1, ls=":",
               label=f"baseline hold {P_baseline:.1f} W")
    ax.set_xlabel("hold power P = i^2 R (W)")
    ax.set_ylabel("holding force (N)")
    ax.set_title("L2: power vs holding force (winding slope + B_sat ceiling)")
    ax.legend(fontsize=8)
    if out_path is None:
        out_path = Path(__file__).parent / "p4_l2_pareto.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _print_summary(front, F_sat, P_baseline, P_cross, baseline_saturated):
    print("== L2: power vs holding-force Pareto ==")
    print(f"  B_sat force ceiling F_sat = {F_sat:.3f} N")
    print(f"  baseline (N={front['N'][0]:.0f}) hold: "
          f"P={P_baseline:.3f} W, F_hold={front['F_hold'][0]:.3f} N, "
          f"F_linear(fictitious)={front['F_linear'][0]:.1f} N")
    print(f"  baseline in saturated region: {baseline_saturated}")
    print(f"  ceiling binds at/above hold power ~{P_cross:.4f} W "
          f"-> hold power above this buys zero extra force")


if __name__ == "__main__":
    generate()
