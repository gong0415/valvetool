"""P5 report: L3 miniaturisation limit and Buckingham-pi scaling (spec §3).

Two headline findings. First, the mechanisms the spec names for L3 -- winding
window, magnetic saturation, seat dP force -- produce no miniaturisation limit
on their own: the force margin is exactly scale-invariant. The limit that does
exist comes from an EMPIRICAL current-density ceiling, giving two orthogonal
closed-form bounds and hence a rectangular feasible region. Second, pi_2 =
tau_e/tau_m is the only scale-dependent group and goes as D, so smaller valves
shift from dead-time dominated toward mechanically dominated.
"""
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import math
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from solenoid_model import limits, sizing, thermal, winding
from solenoid_model.baseline import BASELINE_PARAMS

SCALES = [round(float(s), 4) for s in np.linspace(1.6, 0.4, 13)]
PRESSURES = [float(x) for x in np.linspace(1.0e6, 12.0e6, 12)]
PI2_SCALES = [1.0, 0.8, 0.6, 0.4, 0.2]

_T_AMB = 293.15   # K, P2 thermal-model convention
_T_CEIL = 373.15  # K = 100 degC, design-doc prototype coil ceiling


def _thermal_limited_mmf_ratio(p, s):
    """NI_avail/NI_need at scale s, with NI_avail from the P2 thermal
    model's back-solved steady current (not J_max) -- see
    tests/test_sizing.py's _thermal_limited_ni_ratio for the derivation and
    the honesty note on why this is not the shipped limit. R is derived
    from the SCALED winding window (winding.coil_resistance), because
    scale_params deliberately leaves R_coil_20C at its unscaled value."""
    ps = sizing.scale_params(p, s)
    R20_scaled = winding.coil_resistance(ps.N_turns, ps.A_winding,
                                         ps.k_fill, ps.l_turn_mean)
    ps_r = replace(ps, R_coil_20C=R20_scaled)
    losses = (ps_r.G_th_cond * (_T_CEIL - _T_AMB)
              + ps_r.emissivity * thermal.STEFAN_BOLTZMANN * ps_r.A_rad
              * (_T_CEIL ** 4 - _T_AMB ** 4))
    I_thermal = math.sqrt(losses / thermal.R_coil(_T_CEIL, ps_r))
    NI_avail = ps.N_turns * I_thermal
    NI_need = p.B_sat * ps.g0 / sizing.MU_0
    return NI_avail / NI_need


def generate(out_path=None):
    p = BASELINE_PARAMS
    region = sizing.l3_feasible_region(p, J_max=p.J_max, k_pack=p.k_pack,
                                       scales=SCALES, pressures=PRESSURES)
    s_star = sizing.min_feasible_scale(p, p.J_max)
    dP_max = sizing.pressure_ceiling(p)
    groups = limits.dimensionless_groups(p)

    pi_2_by_scale = []
    for s in PI2_SCALES:
        ps = sizing.scale_params(p, s)
        R = winding.coil_resistance(p.N_turns, ps.A_winding, ps.k_fill,
                                    ps.l_turn_mean)
        pi_2_by_scale.append(
            (s, limits.dimensionless_groups(ps, R=R)["pi_2_time_ratio"]))

    mmf_ratio_baseline = _thermal_limited_mmf_ratio(p, 1.0)

    _plot(p, region, s_star, dP_max, pi_2_by_scale, out_path)
    _print_summary(p, region, s_star, dP_max, groups, pi_2_by_scale,
                   mmf_ratio_baseline)

    return {"region": region, "s_star": s_star, "dP_max": dP_max,
            "groups": groups, "pi_2_by_scale": pi_2_by_scale,
            "mmf_ratio_baseline": mmf_ratio_baseline}


def _plot(p, region, s_star, dP_max, pi_2_by_scale, out_path=None):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    ax = axes[0]
    grid = np.array([[1 if f else 0 for f in row] for row in region["feasible"]])
    # SCALES descends, so grid row 0 is the LARGEST scale -> origin="upper"
    # puts it at the top, matching extent's top=scales[0].
    ax.imshow(grid, aspect="auto", origin="upper", cmap="Greens", vmin=0, vmax=1,
              extent=[region["pressures"][0] / 1e6, region["pressures"][-1] / 1e6,
                      region["scales"][-1], region["scales"][0]])
    ax.axvline(dP_max / 1e6, color="tab:red", lw=1.5,
               label=f"dP_max = {dP_max/1e6:.1f} MPa")
    ax.axhline(s_star, color="tab:blue", lw=1.5,
               label=f"s* = {s_star:.3f} (J_max EMPIRICAL)")
    # Baseline marker: s=1.0 at the baseline seat dP. s=1.0 < s* here, so the
    # star lands just outside (below) the green region -- the shipped valve
    # sits outside its own feasible region under the EMPIRICAL J_max ceiling.
    ax.plot(p.delta_P / 1e6, 1.0, marker="*", markersize=14, color="black",
            markeredgecolor="white", markeredgewidth=0.6, linestyle="none",
            zorder=5, label="baseline valve (s=1)")
    ax.set_xlabel("seat pressure differential (MPa)")
    ax.set_ylabel("scale factor s")
    ax.set_title("L3 feasible region (green)", fontsize=9)
    ax.legend(fontsize=6.5, loc="upper right")

    ax = axes[1]
    ax.loglog(region["scales"], [m * 1e3 for m in region["mass"]], "o-",
              label="mass (k_pack EMPIRICAL)")
    ax.axvline(s_star, color="tab:blue", lw=1.5, ls="--", label=f"s* = {s_star:.3f}")
    ax.set_xlabel("scale factor s")
    ax.set_ylabel("valve mass (g)")
    ax.set_title("mass vs scale (slope 3 => mass ~ D^3)", fontsize=9)
    ax.legend(fontsize=7)

    ax = axes[2]
    ax.plot([s for s, _ in pi_2_by_scale], [v for _, v in pi_2_by_scale], "s-")
    ax.axhline(1.0, color="tab:gray", lw=1, ls=":",
               label="pi_2 = 1 (equal time constants)")
    ax.set_xlabel("scale factor s")
    ax.set_ylabel("pi_2 = tau_e / tau_m")
    ax.set_title("pi_2 ~ D: small valves turn mechanical", fontsize=9)
    ax.legend(fontsize=7)

    fig.suptitle("L3 miniaturisation limit and pi scaling", fontsize=11)
    fig.tight_layout()
    if out_path is None:
        out_path = Path(__file__).parent / "p5_l3_region.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _print_summary(p, region, s_star, dP_max, groups, pi_2_by_scale,
                   mmf_ratio_baseline):
    a_wire = p.A_winding * p.k_fill / p.N_turns
    J_baseline = (p.V_bus / p.R_coil_20C) / a_wire
    print("== L3 miniaturisation limit ==")
    print(f"  minimum feasible scale s* = {s_star:.4f} "
          f"(at J_max = {p.J_max/1e6:.0f} A/mm^2, EMPIRICAL)")
    print(f"  pressure ceiling dP_max   = {dP_max/1e6:.3f} MPa "
          f"({dP_max/1e5:.2f} bar), scale-invariant")
    print(f"  baseline valve mass       = {sizing.valve_mass(p, p.k_pack)*1e3:.2f} g "
          f"(k_pack = {p.k_pack}, EMPIRICAL)")
    print("  note: the force margin (saturation vs. seat dP) is "
          "scale-invariant, and the thermal-limited winding-window MMF "
          "ratio below never drops below 1 either -- so none of the "
          "spec's three named mechanisms alone give a limit; s* comes "
          "from the EMPIRICAL J_max ceiling instead")
    print(f"  winding-window MMF ratio  = {mmf_ratio_baseline:.4f} "
          f"(NI_avail/NI_need at baseline, thermal-limited current per P2 "
          f"model, T_ceiling={_T_CEIL:.2f} K -- see test_sizing.py: this "
          f"ratio's 'never below 1' band is itself an extrapolation "
          f"artifact, see the current-density line below)")
    if s_star > 1.0:
        print(f"  ** baseline s=1.0 is BELOW s*={s_star:.4f}: the shipped "
              f"baseline valve sits OUTSIDE its own feasible region under "
              f"the EMPIRICAL J_max={p.J_max/1e6:.0f} A/mm^2 continuous "
              f"ceiling -- the baseline coil actually runs at "
              f"{J_baseline/1e6:.2f} A/mm^2, above that ceiling, consistent "
              f"with it being a pulsed-duty valve, not a continuous one")
    print("== Buckingham pi (baseline) ==")
    for k, v in groups.items():
        print(f"  {k:24s} = {v:10.4f}")
    print("== pi_2 vs scale ==")
    for s, v in pi_2_by_scale:
        print(f"  s={s:4.2f}: pi_2 = {v:8.4f}")


if __name__ == "__main__":
    generate()
