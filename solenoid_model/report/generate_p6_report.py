"""P6 report: the L4 leak limit, the L5 impact/life limit, and the seat-bounce
leak spike (spec §2.2 bounce, §2.4, §2.6, §3 L4/L5).

Two headline findings. First, sealing at this force level is near-binary: a
PCTFE soft seat percolates with a comfortable ~66 um land while metal-to-metal
would need ~1.65 um -- effectively a knife edge -- and below the threshold the
leak misses the 1e-4 scc/s spec by 2-4 orders (4-6 above 1e-6). Second, P3's zener
flyback speedup carries a bill this report finally computes: closing ~3.6x
faster means striking the seat at ~3.9x the speed and ~15x the energy.
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

from solenoid_model import (dynamics, fluid, flyback, limits, materials,
                            sealing)
from solenoid_model.baseline import BASELINE_PARAMS

_HERE = Path(__file__).resolve().parent

SEATS = [materials.PCTFE_ON_440C, materials.METAL_17_4PH_ON_440C]
W_LANDS = [float(w) for w in np.logspace(-6.3, -3.3, 40)]    # 0.5 .. 500 um
VELOCITIES = [float(v) for v in np.logspace(-1.6, 0.4, 40)]  # 0.025 .. 2.5 m/s
COND = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
SPEC_LEAK_LOOSE = 1e-4     # scc/s GHe, spec §4.3 internal-leak band
SPEC_LEAK_TIGHT = 1e-6
SPEC_CYCLES = (1e5, 1e7)   # spec §4.3 N_cycle band


def _seat_impact_speed(params, circuit):
    """Speed at which the poppet reaches the seat, no bounce."""
    r = dynamics.simulate_closing(params, circuit)
    return abs(r.sol_stroke.y_events[0][0][2])


def _plot_l4(params, path):
    curves = limits.l4_leak_curve(params, SEATS, W_LANDS, COND)
    land_limits = {s.name: sealing.land_width_for_seal(params, s, COND)
                   for s in SEATS}
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for seat in SEATS:
        rows = [r for r in curves[seat.name] if not r["sealed"]]
        ax.plot([r["w_land"] * 1e6 for r in rows],
                [r["leak"] for r in rows], marker=".",
                label=f"{seat.name}: leaking branch")
        ax.axvline(land_limits[seat.name] * 1e6, linestyle="--", linewidth=1,
                   label=f"{seat.name}: seals below "
                         f"{land_limits[seat.name] * 1e6:.2f} um")
    ax.axhline(SPEC_LEAK_LOOSE, color="k", linestyle=":", linewidth=1)
    ax.axhline(SPEC_LEAK_TIGHT, color="k", linestyle=":", linewidth=1)
    ax.text(0.02, SPEC_LEAK_LOOSE, " spec 1e-4", va="bottom", fontsize=7)
    ax.text(0.02, SPEC_LEAK_TIGHT, " spec 1e-6", va="bottom", fontsize=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("seal land width (um)")
    ax.set_ylabel("internal leak (scc/s GHe)")
    ax.set_title("L4: sealing is near-binary -- the threshold is the result, "
                 "not the magnitude")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return land_limits


def _plot_l5(params, path, v_ops):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    shakedown_edges = {}
    for seat in SEATS:
        all_rows = limits.l5_life_curve(params, seat, VELOCITIES)
        rows = [r for r in all_rows if np.isfinite(r["N_cycle"])]
        line, = ax.plot([r["v_impact"] for r in rows],
                        [r["N_cycle"] for r in rows], marker=".",
                        label=f"{seat.name} (plastic branch)")
        v_finite = [r["v_impact"] for r in rows]
        if v_finite:
            shakedown_edges[seat.name] = (min(v_finite), line.get_color())
    ax.axhspan(SPEC_CYCLES[0], SPEC_CYCLES[1], alpha=0.12, color="tab:green",
               label="spec 1e5-1e7 cycles")
    for label, v in v_ops.items():
        ax.axvline(v, linestyle="--", linewidth=1,
                   label=f"{label} close: {v:.3f} m/s")
    ax.set_xscale("log")
    ax.set_yscale("log")

    # I4 fix: l5_life_curve reports N_cycle = +inf under shakedown (p_max <
    # H, the FULLY-PLASTIC hardness -- NOT first yield, see
    # impact.is_shakedown's docstring and MODEL_NOTES.md #9) and those rows
    # are excluded from the curves above. Left unmarked, the resulting gap
    # just reads as missing data and a reader naturally interpolates across
    # it -- exactly wrong, since the model is actually asserting unlimited
    # life there. Shade it per seat instead so that assertion is explicit.
    x_lo, x_hi = ax.get_xlim()
    for seat_name, (edge, color) in shakedown_edges.items():
        ax.axvspan(x_lo, edge, color=color, alpha=0.10,
                   label=f"{seat_name}: shakedown -- claims UNLIMITED life")
    ax.set_xlim(x_lo, x_hi)
    ax.text(0.02, 0.03,
           "shading = p_max < H (fully-plastic threshold, NOT first yield --\n"
           "impacts down to ~0.53*H are already yielding but untracked here)",
           transform=ax.transAxes, fontsize=6, va="bottom", ha="left",
           style="italic")

    ax.set_xlabel("seat impact speed (m/s)")
    ax.set_ylabel("estimated cycles to seal-geometry failure")
    ax.set_title("L5: impact speed vs cycle life (EMPIRICAL wear "
                 "coefficients -- order of magnitude only)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _peak_bounce_lift(params, circuit, seat):
    """Max poppet position [m] reached during any closing-bounce flight."""
    res = dynamics.simulate_closing(params, circuit, seat=seat)
    return max(seg.y[1].max() for seg in res.bounce.segments[1:])


def _plot_bounce(params, path, seat):
    res = dynamics.simulate_closing(params, flyback.ZenerFlyback(), seat=seat)
    b = res.bounce
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for k, seg in enumerate(b.segments):
        ax.plot(seg.t * 1e3, seg.y[1] * 1e6, color="tab:blue",
                label="poppet position" if k == 0 else None)
    ax.axvline(b.t_first * 1e3, linestyle="--", linewidth=1, color="tab:red",
               label=f"first contact {b.t_first * 1e3:.3f} ms")
    ax.axvline(b.t_settle * 1e3, linestyle=":", linewidth=1, color="tab:green",
               label=f"settled {b.t_settle * 1e3:.3f} ms "
                     f"({b.n_bounce} impacts)")
    ax.set_xlabel("time since stop release (ms)")
    ax.set_ylabel("poppet position (um)")
    ax.set_title(f"Seat bounce on closing, zener flyback, {seat.name}")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return b


def generate(l4_path=None, l5_path=None, bounce_path=None):
    p = BASELINE_PARAMS
    l4_path = _HERE / "p6_l4_leak.png" if l4_path is None else Path(l4_path)
    l5_path = _HERE / "p6_l5_life.png" if l5_path is None else Path(l5_path)
    bounce_path = (_HERE / "p6_bounce_x_t.png" if bounce_path is None
                   else Path(bounce_path))

    land_limits = _plot_l4(p, l4_path)
    v_ops = {"diode": _seat_impact_speed(p, flyback.DiodeFlyback()),
             "zener": _seat_impact_speed(p, flyback.ZenerFlyback())}
    _plot_l5(p, l5_path, v_ops)
    soft = materials.PCTFE_ON_440C
    _plot_bounce(p, bounce_path, soft)
    spike = limits.bounce_leak_spike(p, soft, flyback.ZenerFlyback(), COND,
                                     fluid.N2)
    # Third arm of the "zener speedup bill" alongside t_close x0.28 and
    # impact energy x15.2 (see design doc finding 3): the harder zener
    # bounce also leaks more mass, even though it spends LESS total time
    # off the seat, because it lifts the poppet further and a wider orifice
    # dominates a shorter duration.
    spike_diode = limits.bounce_leak_spike(p, soft, flyback.DiodeFlyback(),
                                           COND, fluid.N2)
    diode_lift = _peak_bounce_lift(p, flyback.DiodeFlyback(), soft)
    zener_lift = _peak_bounce_lift(p, flyback.ZenerFlyback(), soft)

    # design doc finding 6: the opening bounce ratios converge to e from
    # above, because the armature reaches the stop with almost no magnetic
    # force margin over the spring+pressure load
    open_bounce = dynamics.simulate_opening(p, seat=soft)
    open_ratios = [b / a for a, b in
                   zip(open_bounce.v_impacts, open_bounce.v_impacts[1:])]

    summary = {
        "open_ratios": open_ratios,
        "open_n_bounce": open_bounce.n_bounce,
        "open_t_first_ms": open_bounce.t_first * 1e3,
        "open_t_settle_ms": open_bounce.t_settle * 1e3,
        "pctfe_land_limit_um": land_limits[soft.name] * 1e6,
        "metal_land_limit_um":
            land_limits[materials.METAL_17_4PH_ON_440C.name] * 1e6,
        "v_diode": v_ops["diode"],
        "v_zener": v_ops["zener"],
        "zener_speed_ratio": v_ops["zener"] / v_ops["diode"],
        "zener_energy_ratio": (v_ops["zener"] / v_ops["diode"]) ** 2,
        "n_bounce": spike["n_bounce"],
        "t_open_total": spike["t_open_total"],
        "mass_leaked": spike["mass_leaked"],
        "diode_n_bounce": spike_diode["n_bounce"],
        "diode_t_open_total": spike_diode["t_open_total"],
        "diode_mass_leaked": spike_diode["mass_leaked"],
        "bounce_mass_ratio": spike["mass_leaked"] / spike_diode["mass_leaked"],
        "bounce_time_ratio": spike["t_open_total"] / spike_diode["t_open_total"],
        "diode_bounce_lift_um": diode_lift * 1e6,
        "zener_bounce_lift_um": zener_lift * 1e6,
    }
    return {"l4": l4_path, "l5": l5_path, "bounce": bounce_path,
            "summary": summary, "land_limits": land_limits,
            "spike": spike, "spike_diode": spike_diode}


if __name__ == "__main__":
    out = generate()
    s = out["summary"]
    print("== L4: land width needed to seal (percolation threshold) ==")
    print(f"  PCTFE/440C   : <= {s['pctfe_land_limit_um']:8.3f} um "
          f"(manufacturable)")
    print(f"  17-4PH/440C  : <= {s['metal_land_limit_um']:8.3f} um "
          f"(effectively a knife edge)")
    print("  below the threshold the leak misses the 1e-4 scc/s spec by "
          "2-4 orders (4-6 above 1e-6): sealing is near-binary, so only the "
          "threshold is a trustworthy output of this model (see design doc "
          "finding 2)")
    print("== L5: the bill for P3's zener speedup ==")
    print(f"  diode close  : {s['v_diode']:.4f} m/s")
    print(f"  zener close  : {s['v_zener']:.4f} m/s "
          f"({s['zener_speed_ratio']:.2f}x speed, "
          f"{s['zener_energy_ratio']:.2f}x impact energy)")
    print("== Seat bounce on closing, PCTFE seat -- zener speedup bill, "
          "arm 3 (arms 1-2: t_close x0.28, impact energy x15.2) ==")
    print(f"  diode : {s['diode_n_bounce']} impacts, off-seat "
          f"{s['diode_t_open_total'] * 1e3:.4f} ms, leaked "
          f"{s['diode_mass_leaked'] * 1e9:.4f} ug, lifts poppet "
          f"{s['diode_bounce_lift_um']:.2f} um (N2 at MEOP)")
    print(f"  zener : {s['n_bounce']} impacts, off-seat "
          f"{s['t_open_total'] * 1e3:.4f} ms, leaked "
          f"{s['mass_leaked'] * 1e9:.4f} ug, lifts poppet "
          f"{s['zener_bounce_lift_um']:.2f} um (N2 at MEOP)")
    print(f"  zener/diode: {s['bounce_mass_ratio']:.2f}x the leaked mass in "
          f"{s['bounce_time_ratio']:.3f}x the off-seat time -- the harder "
          "zener bounce lifts the poppet further, opening a wider orifice "
          "that outweighs the shorter total off-seat time")
    print("== Bounce on opening (design doc finding 6) ==")
    print(f"  impacts      : {s['open_n_bounce']}, first contact "
          f"{s['open_t_first_ms']:.4f} ms, settled "
          f"{s['open_t_settle_ms']:.4f} ms")
    print("  speed ratios : "
          + ", ".join(f"{x:.3f}" for x in s["open_ratios"]))
    print("  these converge to e from ABOVE, not to e exactly: the armature "
          "reaches the stop at ~11% of steady-state current, where magnetic "
          "force barely exceeds the spring+pressure load, so the first "
          "rebound flies far and long. The energy it gains comes from the "
          "current still rising during that flight -- at fixed current the "
          "magnetic force is conservative and would do no net work over a "
          "closed excursion; the field's position-dependence only amplifies "
          "the effect by widening the range of x over which the rising "
          "current is sampled.")
