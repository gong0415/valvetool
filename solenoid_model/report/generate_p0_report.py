from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from solenoid_model import dynamics, magnetics
from solenoid_model.baseline import BASELINE_PARAMS

REPORT_DIR = Path(__file__).parent


def generate():
    params = BASELINE_PARAMS
    sol = dynamics.simulate_opening(params, t_max=0.05)
    if len(sol.t_events[0]) == 0:
        raise RuntimeError("armature did not reach full stroke within t_max; check baseline parameters")
    t_open = sol.t_events[0][0]

    t_ms = sol.t * 1e3
    i_t, x_t, _v_t = sol.y

    fig, ax = plt.subplots()
    ax.plot(t_ms, i_t)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("coil current (A)")
    ax.set_title(f"i(t), t_open={t_open * 1e3:.2f} ms")
    fig.savefig(REPORT_DIR / "current_i_t.png")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(t_ms, x_t * 1e3)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("armature displacement (mm)")
    ax.set_title("x(t)")
    fig.savefig(REPORT_DIR / "position_x_t.png")
    plt.close(fig)

    F_mag_t = [
        magnetics.magnetic_force_closing(params.g0 - x, i, params)
        for i, x in zip(i_t, x_t)
    ]
    fig, ax = plt.subplots()
    ax.plot(t_ms, F_mag_t)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("magnetic force (N)")
    ax.set_title("F_mag(t)")
    fig.savefig(REPORT_DIR / "force_F_t.png")
    plt.close(fig)

    return t_open


if __name__ == "__main__":
    result = generate()
    print(f"t_open = {result * 1e3:.3f} ms")
