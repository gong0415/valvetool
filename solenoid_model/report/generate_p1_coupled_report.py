"""P1 baseline comparison: dry run vs N2-coupled opening transient."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from solenoid_model import dynamics, fluid
from solenoid_model.baseline import BASELINE_PARAMS

REPORT_DIR = Path(__file__).parent
N2_VACUUM = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)


def generate():
    params = BASELINE_PARAMS
    sol_dry = dynamics.simulate_opening(params)
    sol_cpl = dynamics.simulate_opening(params, medium=fluid.N2, cond=N2_VACUUM)
    t_dry = sol_dry.t_events[0][0]
    t_cpl = sol_cpl.t_events[0][0]

    fig, ax = plt.subplots()
    ax.plot(sol_dry.t * 1e3, sol_dry.y[1] * 1e3, label=f"dry: t_open={t_dry*1e3:.3f} ms")
    ax.plot(sol_cpl.t * 1e3, sol_cpl.y[1] * 1e3, "--",
            label=f"N2 coupled: t_open={t_cpl*1e3:.3f} ms")
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("armature displacement (mm)")
    ax.set_title("x(t): dry vs N2-coupled (2.4 MPa -> vacuum)")
    ax.legend()
    fig.savefig(REPORT_DIR / "dry_vs_coupled_x_t.png", dpi=150)
    plt.close(fig)

    return {
        "t_open_dry": t_dry,
        "t_open_coupled": t_cpl,
        "cv_full_open": fluid.cv_from_geometry(params.x_stroke, params),
        "mdot_full_open": fluid.mdot_gas(params.x_stroke, params, fluid.N2, N2_VACUUM),
        "flow_force_full_open": fluid.flow_force_gas(params.x_stroke, params, fluid.N2, N2_VACUUM),
    }


if __name__ == "__main__":
    r = generate()
    print(f"dry     t_open = {r['t_open_dry']*1e3:.3f} ms")
    print(f"coupled t_open = {r['t_open_coupled']*1e3:.3f} ms  (+{(r['t_open_coupled']-r['t_open_dry'])*1e6:.1f} µs)")
    print(f"Cv={r['cv_full_open']:.4f}  mdot={r['mdot_full_open']*1e3:.3f} g/s  F_flow={r['flow_force_full_open']:.3f} N")
