"""P3 report: closing transient / flyback circuit comparison (spec §5
response characteristic 2). Compares diode / zener / RC-snubber freewheel
paths and checks the spec's claimed "zener clamp accelerates the closing
transient 3-5x vs diode" figure against the baseline parameter set.
"""
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from solenoid_model import dynamics, flyback
from solenoid_model.baseline import BASELINE_PARAMS

# near-critical damping (zeta ~= 1.5) at the baseline's full-open
# inductance; see docs/spec/PARAMS.md for the damping-ratio formula
RC_CIRCUIT = flyback.RCFlyback(R_snub=14655.9, C_snub=1e-7)


def generate(out_path=None):
    p = BASELINE_PARAMS
    circuits = [
        ("diode", flyback.DiodeFlyback()),
        ("zener", flyback.ZenerFlyback()),
        ("RC snubber", RC_CIRCUIT),
    ]

    results = {name: dynamics.simulate_closing(p, circuit) for name, circuit in circuits}
    ratio = results["diode"].t_close / results["zener"].t_close

    _plot(results, out_path)
    _print_summary(results, ratio)

    return {"results": results, "diode_to_zener_ratio": ratio}


def _plot(results, out_path=None):
    fig, ax = plt.subplots()
    for name, result in results.items():
        t = np.concatenate([result.sol_hold.t, result.t_release + result.sol_stroke.t])
        i = np.concatenate([result.sol_hold.y[0], result.sol_stroke.y[0]])
        ax.plot(t * 1e3, i, label=f"{name} (t_close={result.t_close * 1e3:.2f} ms)")
    ax.set_xlabel("time since de-energization (ms)")
    ax.set_ylabel("coil current (A)")
    ax.set_title("P3: closing-transient current decay by flyback circuit")
    ax.legend(fontsize=8)
    if out_path is None:
        out_path = Path(__file__).parent / "p3_flyback_current_t.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _print_summary(results, ratio):
    print("== P3: closing transient / flyback comparison ==")
    for name, result in results.items():
        print(f"  {name:12s} t_release={result.t_release * 1e3:8.4f} ms  "
              f"t_close={result.t_close * 1e3:8.4f} ms")
    print(f"  diode/zener t_close ratio = {ratio:.3f} "
          f"(spec claims zener clamp accelerates 3-5x)")
    print(f"  within spec's claimed 3-5x range: {3.0 <= ratio <= 5.0}")


if __name__ == "__main__":
    generate()
