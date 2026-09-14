import sys
from dataclasses import dataclass, replace
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from solenoid_model import (drive, dynamics, flyback, fluid, impact, limits,
                            magnetics, materials, sealing, sizing, thermal,
                            winding)
from solenoid_model import cases
from solenoid_model import layout
from solenoid_model.report import generate_layout
from solenoid_model.baseline import BASELINE_PARAMS
from solenoid_model.params import ValveParams


@dataclass
class SimulationResult:
    t_open: float
    sol: object
    saturated: bool


def compute_result(params, t_max=0.05, medium=None, cond=None, T_coil=None):
    sol = dynamics.simulate_opening(params, t_max=t_max, medium=medium,
                                    cond=cond, T_coil=T_coil)
    if len(sol.t_events[0]) == 0:
        return None
    t_open = sol.t_events[0][0]
    saturated = any(
        magnetics.saturation_check(params.g0 - x, i, params)
        for i, x in zip(sol.y[0], sol.y[1])
    )
    return SimulationResult(t_open=t_open, sol=sol, saturated=saturated)


def compute_opening_bounce(params, seat, t_max=0.05, medium=None, cond=None,
                           T_coil=None):
    return dynamics.simulate_opening(params, t_max=t_max, medium=medium,
                                     cond=cond, T_coil=T_coil, seat=seat)


@dataclass
class ClosingSummary:
    t_release: float
    t_close: float
    v_impact: float
    result: object


def compute_closing(params, circuit, seat=None):
    r = dynamics.simulate_closing(params, circuit, seat=seat)
    return ClosingSummary(t_release=r.t_release, t_close=r.t_close,
                          v_impact=abs(r.sol_stroke.y[2][-1]), result=r)


def compute_bounce_spike(params, seat, circuit, cond, medium):
    return limits.bounce_leak_spike(params, seat, circuit, cond, medium)


def parse_float_list(text):
    return [float(tok) for tok in text.split(",") if tok.strip()]


def l1_sweep_cases(kind, values, params):
    if kind == "voltage":
        return [(f"{v:g} V", replace(params, V_bus=v), params.R_coil_20C)
                for v in values]
    if kind == "turns_window":
        return [(f"N={n:g}", replace(params, N_turns=n),
                 winding.coil_resistance(n, params.A_winding, params.k_fill,
                                         params.l_turn_mean))
                for n in values]
    raise ValueError(f"unknown sweep kind: {kind}")


def render_limits_tab(params):
    st.subheader("L1 響應時間下限")
    tau_e = limits.electrical_time_constant(params)
    i_th = limits.motion_threshold_current(params)
    bound = limits.response_time_bound(params)
    naive = limits.naive_response_time_estimate(params)
    c1, c2, c3 = st.columns(3)
    c1.metric("τ_e 電氣時間常數", f"{tau_e * 1e3:.3f} ms")
    c2.metric("運動門檻電流 i_th", f"{i_th:.4f} A")
    if bound is None:
        st.error("穩態電流達不到運動門檻——此參數下閥不會開，L1 無界可算")
    else:
        c3.metric("L1 響應時間下限", f"{bound['t_bound'] * 1e3:.3f} ms",
                  help="t_i（RL 升流死區，模型內精確）+ t_x（定加速度渡越，樂觀）")
        naive_txt = (f"{naive * 1e3:.3f} ms" if naive is not None
                    else "N/A（此工況閥不會開）")
        st.caption(f"t_i = {bound['t_i'] * 1e3:.3f} ms、"
                   f"t_x = {bound['t_x'] * 1e3:.3f} ms。"
                   f"規格分解式 naive 估計 = "
                   f"{naive_txt}——⚠️ 它**不是上下界**"
                   "（實測落在真值 0.43–3.11 倍之間），有效宣稱請用左側下限")

    sweep_kind = st.selectbox("掃描維度", ["電壓 (V)", "匝數（窗口一致 R）"])
    default_vals = ("18, 22, 28, 36, 50" if sweep_kind == "電壓 (V)"
                    else "1000, 1500, 2000, 3000, 4000")
    vals_text = st.text_input("掃描點（逗號分隔）", value=default_vals)
    if st.button("執行 L1 掃描（每點跑一次 ODE）", key="btn_l1_sweep"):
        kind = "voltage" if sweep_kind == "電壓 (V)" else "turns_window"
        sw = limits.l1_sweep(l1_sweep_cases(kind, parse_float_list(vals_text),
                                            params))
        st.session_state["l1_sweep"] = (params, sw)
    if "l1_sweep" in st.session_state:
        stored_params, sw = st.session_state["l1_sweep"]
        if stored_params != params:
            st.caption("⚠️ 顯示的結果對應舊參數——請重按掃描")
        fmt = lambda t: f"{t * 1e3:.3f}" if t is not None else "—"
        st.dataframe({
            "case": sw["label"],
            "t_open (ms)": [fmt(t) for t in sw["t_open"]],
            "L1 下限 (ms)": [fmt(t) for t in sw["t_bound"]],
            "naive (ms)": [fmt(t) for t in sw["t_naive"]],
            "吸得上": ["✅" if p else "❌" for p in sw["pulls_in"]],
        })

    st.subheader("L2 保持力–功率 Pareto（繞線窗口）")
    c1, c2, c3 = st.columns(3)
    N_lo = c1.number_input("N 起", value=2000.0, min_value=1.0)
    N_hi = c2.number_input("N 迄", value=24000.0, min_value=1.0)
    N_pts = c3.number_input("點數", value=45, min_value=2, step=1)
    front = limits.l2_pareto_front(params, list(np.linspace(N_lo, N_hi,
                                                            int(N_pts))))
    F_sat = limits.saturation_force_ceiling(params)
    st.metric("B_sat 保持力天花板", f"{F_sat:.2f} N",
              help="F_sat = B_sat²·A_gap/(2μ₀)：不論線圈功率多大，保持力到此為止")
    fig, ax = plt.subplots()
    ax.plot(front["P_hold"], front["F_linear"], "--", color="gray",
            label="linear (fictitious above B_sat)")
    ax.plot(front["P_hold"], front["F_hold"], color="tab:blue",
            label="F_hold (B_sat-capped)")
    ax.axhline(F_sat, color="tab:red", lw=1, label="F_sat ceiling")
    ax.set_xlabel("hold power (W)")
    ax.set_ylabel("holding force (N)")
    ax.set_xscale("log")
    ax.set_title("L2 Pareto front over turns sweep")
    ax.legend()
    st.pyplot(fig)
    plt.close(fig)
    st.caption("窗口固定時 N²/R 不變 → 前緣斜率固定，天花板由 B_sat 提供；"
               "基準閥保持工況深度飽和（spike-and-hold 論點，P4）")

    st.subheader("L3 小型化可行域 + Buckingham π")
    groups = limits.dimensionless_groups(params)
    st.table({"π 群": ["π₁ 行程/氣隙", "π₂ τ_e/τ_m", "π₃ 磁力數",
                       "π₄ 壓力負載", "π₅ 飽和裕度"],
              "值": [f"{v:.4g}" for v in groups.values()]})
    st.caption("π₂ 是唯一隨尺度變的群（∝D）：小閥從死區主導轉為機械主導（P5）。"
               "本表對未縮放參數以預設 R 計算——縮放律曲線見 P5 報告腳本")
    c1, c2 = st.columns(2)
    k_pack = c1.number_input("包裝質量係數 k_pack", value=params.k_pack,
                             help="整閥/活性質量比")
    J_max = c2.number_input("電流密度上限 J_max (A/m²)", value=params.J_max,
                            format="%.2e", help="連續工作的繞組電流密度天花板")
    st.caption("⚠️ k_pack 與 J_max 皆為 EMPIRICAL 工程係數，非第一性推導；"
               "ODE 模擬不使用它們")
    c1, c2, c3 = st.columns(3)
    c1.metric("最小可行縮放 s*", f"{sizing.min_feasible_scale(params, J_max):.3f}",
              help="窗口 MMF ∝D² vs 氣隙需求 ∝D → 縮小有硬地板")
    c2.metric("壓力天花板", f"{sizing.pressure_ceiling(params) / 1e6:.2f} MPa",
              help="磁路壓得住的最大閥座壓差，尺度不變")
    c3.metric("整閥質量估計", f"{sizing.valve_mass(params, k_pack) * 1e3:.1f} g")
    scales = [round(float(s), 4) for s in np.linspace(1.6, 0.4, 13)]
    pressures = [float(x) for x in np.linspace(1.0e6, 12.0e6, 12)]
    region = sizing.l3_feasible_region(params, J_max=J_max, k_pack=k_pack,
                                       scales=scales, pressures=pressures)
    grid = np.array([[1 if f else 0 for f in row]
                     for row in region["feasible"]])
    fig, ax = plt.subplots()
    ax.imshow(grid, origin="upper", aspect="auto", cmap="RdYlGn",
              extent=[pressures[0] / 1e6, pressures[-1] / 1e6,
                      scales[-1], scales[0]])
    ax.set_xlabel("seat pressure differential (MPa)")
    ax.set_ylabel("scale factor s")
    ax.set_title("L3 feasible region (green = feasible)")
    st.pyplot(fig)
    plt.close(fig)


def render_sealing_tab(params, seat, seat_mats, medium, cond):
    st.subheader("L4 密封門檻與氦漏率")
    w_max = sealing.land_width_for_seal(params, seat, cond)
    sealed = sealing.is_percolated(params, seat, cond)
    leak = sealing.leak_rate(params, seat, cond)
    c1, c2, c3 = st.columns(3)
    c1.metric("密封 land 寬度上限", f"{w_max * 1e6:.2f} µm",
              help="此閥座材料對在當前密封力下，實接觸面積分數達滲流門檻 φ_c=0.42 的最大 land 寬度")
    c2.metric("接觸應力", f"{sealing.contact_stress(params, cond) / 1e6:.1f} MPa")
    c3.metric("目前 w_land 判定",
              "密封" if sealed else "不密封",
              delta=None if sealed else f"漏率 ≈ {leak:.2e} scc/s",
              delta_color="inverse")
    st.warning("⚠️ 本模型可信的是**密封/不密封的門檻位置**，不是漏率絕對值"
               "（門檻以下實測比規格線高 2–4 個數量級；P6 發現 2）")

    w_lands = [float(w) for w in np.logspace(-6.3, -3.3, 40)]  # 0.5–500 µm
    curves = limits.l4_leak_curve(params, [materials.PCTFE_ON_440C,
                                           materials.METAL_17_4PH_ON_440C],
                                  w_lands, cond)
    fig, ax = plt.subplots()
    for name, rows in curves.items():
        xs = [r["w_land"] * 1e6 for r in rows if not r["sealed"]]
        ys = [r["leak"] for r in rows if not r["sealed"]]
        ax.plot(xs, ys, marker=".", label=f"{name} (unsealed)")
    ax.axhline(1e-4, color="tab:red", lw=1, ls="--", label="spec 1e-4 scc/s")
    ax.axhline(1e-6, color="tab:red", lw=1, ls=":", label="spec 1e-6 scc/s")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("land width (µm)")
    ax.set_ylabel("He leak rate (scc/s)")
    ax.set_title("L4: land width vs leak, two built-in seats")
    ax.legend(fontsize=8)
    st.pyplot(fig)
    plt.close(fig)
    st.caption("曲線只畫未密封點（密封點漏率為 0，對數軸不可畫）；"
               "兩條內建座對照，現值指標用側欄選的閥座")

    st.subheader("材料與環境準則")
    dose = st.number_input("任務輻射劑量 (rad(Si))", value=1.0e5, format="%.2e")
    prop = st.selectbox("推進劑", ["GN2", "GHe", "hydrazine", "MMH", "NTO",
                                    "H2O2"])
    rows = []
    for mat in seat_mats:
        rows.append({
            "材料": mat.name,
            "出氣 (TML/CVCM)": ("✅" if materials.outgassing_ok(mat) else "❌")
                               + f" {mat.TML}/{mat.CVCM} %",
            "輻射裕度": f"{materials.radiation_margin(mat, dose):.1f}×",
            f"相容性:{prop}": materials.propellant_compatibility(mat, prop),
        })
    st.table(rows)
    risk = materials.cold_weld_risk(seat_mats[0], seat_mats[1],
                                    sealing.contact_stress(params, cond))
    st.caption(f"冷焊風險（真空金屬對金屬、當前接觸應力）：**{risk}**。"
               "材料數字皆 LITERATURE 等級（ASTM E595 門檻除外），非量測值")

    st.subheader("L5 撞擊壽命")
    closing = st.session_state.get("closing")
    v_default = closing[1].v_impact if closing else 0.39
    v_probe = st.number_input("評估撞擊速度 (m/s)", value=float(v_default),
                              min_value=0.001,
                              help="預設取頁籤 2 最近一次關閉模擬的觸座速度；無則 0.39（二極體基準）")
    c1, c2, c3 = st.columns(3)
    c1.metric("Hertz 峰值接觸應力",
              f"{impact.contact_pressure_max(params.m_arm, v_probe, seat, params.R_tip) / 1e6:.1f} MPa")
    c2.metric("接觸時間",
              f"{impact.contact_duration(params.m_arm, v_probe, seat, params.R_tip) * 1e6:.1f} µs")
    c3.metric("Shakedown（純彈性）",
              "是" if impact.is_shakedown(params.m_arm, v_probe, seat,
                                          params.R_tip) else "否")
    velocities = [float(v) for v in np.logspace(-1.6, 0.4, 40)]
    rows = limits.l5_life_curve(params, seat, velocities)
    xs = [r["v_impact"] for r in rows if not r["shakedown"]]
    ys = [r["N_cycle"] for r in rows if not r["shakedown"]]
    fig, ax = plt.subplots()
    ax.plot(xs, ys, marker=".")
    if any(r["shakedown"] for r in rows):
        v_sd = max(r["v_impact"] for r in rows if r["shakedown"])
        ax.axvspan(velocities[0], v_sd, color="tab:green", alpha=0.15,
                   label="shakedown (elastic; fatigue not modelled)")
        ax.legend(fontsize=8)
    if closing:
        ax.axvline(closing[1].v_impact, color="tab:red", lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("impact speed (m/s)")
    ax.set_ylabel("cycle life estimate")
    ax.set_title("L5: impact speed vs life")
    st.pyplot(fig)
    plt.close(fig)
    st.warning("⚠️ 壽命模型每個係數都是 EMPIRICAL（Archard K、失效深度分數），"
               "結果只能當**數量級**參考，不得作為鑑定數字；"
               "shakedown 區內壽命由本模型未涵蓋的疲勞機制決定")

    st.subheader("關閉回彈漏量尖峰（L4×L5 耦合）")
    if not isinstance(medium, fluid.GasProperties):
        st.info("回彈尖峰的離座流量走氣體壅塞流（mdot_gas）——"
                "請在側欄「流體耦合」改選氣體介質")
    else:
        spike_circ = st.selectbox("續流電路（回彈用）", ["稽納", "二極體"],
                                  help="預設稽納：關閉最快、回彈最兇的工況（P6 報告基準）")
        if st.button("計算回彈漏量尖峰（跑關閉 ODE）", key="btn_spike"):
            circuit = (flyback.ZenerFlyback() if spike_circ == "稽納"
                       else flyback.DiodeFlyback())
            st.session_state["spike"] = (params, compute_bounce_spike(
                params, seat, circuit, cond, medium))
        if "spike" in st.session_state:
            stored_params, spike = st.session_state["spike"]
            if stored_params != params:
                st.caption("⚠️ 顯示的結果對應舊參數——請重按計算")
            c1, c2, c3, c4 = st.columns(4)
            # 標籤刻意不同於頁籤 1 的「回彈次數」：測試以 label 為 key 收集
            # 全域 metric，同名會互相蓋掉
            c1.metric("關閉回彈次數", f"{spike['n_bounce']}")
            c2.metric("回彈離座總時間", f"{spike['t_open_total'] * 1e6:.0f} µs")
            c3.metric("回彈漏出質量", f"{spike['mass_leaked'] * 1e6:.3f} mg")
            c4.metric("靜態密封漏率", f"{spike['sealed_leak']:.2e} scc/s")
            st.caption("回彈期間 poppet 離座，漏量走孔口壅塞流——"
                       "一次關閉的尖峰漏出質量 vs 密封後的穩態漏率對照")


def compute_layout(params, coil_OD=18.0e-3):
    """Pure function: the whole physical dimension chain for one params."""
    env = layout.envelope(params, coil_OD=coil_OD)
    spring = layout.spring_geometry(params, **layout.SPRING_DESIGN)
    t_arm = layout.armature_thickness(params)
    _, _, rel = layout.armature_mass_consistency(params, t_arm)
    return {
        "envelope": env,
        "wall": env["wall"],
        "stack": env["stack"],
        "spring": spring,
        "armature_thickness": t_arm,
        "armature_rel_error": rel,
        "yoke_area": layout.yoke_area(params),
        "pole_diameter": layout.pole_diameter(params),
        "end_plate": layout.end_plate_thickness(params),
        "empirical": layout.LAYOUT_EMPIRICAL,
    }


def render_layout_tab(params):
    r = compute_layout(params)
    st.subheader("實體佈局尺寸鏈")
    c1, c2, c3 = st.columns(3)
    c1.metric("外徑", f"{r['envelope']['OD']*1e3:.1f} mm")
    c2.metric("總長", f"{r['envelope']['L']*1e3:.2f} mm")
    c3.metric("極面直徑", f"{r['pole_diameter']*1e3:.2f} mm")

    st.markdown("**軸向尺寸鏈**")
    st.table({"段": [n for n, _ in r["stack"]["segments"]],
              "厚度 (mm)": [f"{t*1e3:.2f}"
                            for _, t in r["stack"]["segments"]]})

    st.markdown("**外殼壁厚：三個獨立下限**")
    w = r["wall"]
    st.table({"下限": ["磁性", "結構", "製造"],
              "值 (mm)": [f"{w.magnetic*1e3:.2f}",
                          f"{w.structural*1e3:.2f}",
                          f"{w.manufacturing*1e3:.2f}"]})
    st.caption(f"採用 {w.adopted*1e3:.2f} mm — {w.reason}")

    st.markdown("**彈簧**")
    s = r["spring"]
    L_installed = s.L_free - params.F_preload / params.k_spring
    st.write(f"線徑 {s.d_wire*1e3:.2f} mm / 中徑 {s.D_coil*1e3:.2f} mm / "
             f"外徑 {(s.D_coil + s.d_wire)*1e3:.2f} mm / "
             f"有效圈數 {s.n_active:.2f} / 彈簧指數 {s.index:.2f} / "
             f"最大剪應力 {s.tau_max/1e6:.0f} MPa")
    st.caption(
        f"自由長 {s.L_free*1e3:.2f} mm（未裝入）→ 閥關安裝長 "
        f"{L_installed*1e3:.2f} mm（已被 F_preload 壓縮 "
        f"{params.F_preload/params.k_spring*1e3:.2f} mm）→ 併圈長 "
        f"{s.L_solid*1e3:.2f} mm。尺寸鏈的「彈簧腔」用自由長，是組裝時"
        "須容得下的腔體尺寸，非裝配後長度。")

    if r["armature_rel_error"] > 0.05:
        st.warning(
            f"銜鐵厚度與宣告的 m_arm 不一致（偏差 "
            f"{r['armature_rel_error']*100:.1f}%）")

    st.info(
        "⚠️ 濕式銜鐵的隔離套厚度落在徑向磁路上，但 magnetics.reluctance "
        "未建模此項 → 實際吸力低於模型報告的裕度。")

    st.markdown("**EMPIRICAL 項目**（非第一性推導）")
    st.table({
        "項目": list(r["empirical"].keys()),
        "值": [f"{e['value']*1e3:.2f} mm" for e in r["empirical"].values()],
        "理由": [e["reason"] for e in r["empirical"].values()],
    })

    for path, cap in ((generate_layout.SECTION_OUT, "按比例剖面圖"),
                      (generate_layout.ACTUATION_OUT,
                       "作動前後對照（閥關 ↔ 閥開）"),
                      (generate_layout.ARCH_OUT, "架構圖")):
        if path.exists():
            st.image(str(path), caption=cap)


_MEDIUM_ORDER = [fluid.N2, fluid.XE, fluid.HE, fluid.WATER_20C, fluid.LN2_77K]


def _medium_index(medium):
    """Index of a preset medium in the selectbox; custom media fall back to N2."""
    for i, preset in enumerate(_MEDIUM_ORDER):
        if medium is preset:
            return i
    return 0


def _key_suffix(case_name):
    """Stable per-case widget-key suffix (Streamlit keys must be identifiers)."""
    return "".join(ch if ch.isalnum() else "_" for ch in case_name)


def _render_hold_phase(params, case):
    """Static force-balance summary for the hold phase of a peak-and-hold case."""
    gap_open = params.g0 - params.x_stroke
    i_peak = params.V_bus / params.R_coil_20C
    # Same force balance drive.hold_current inverts, evaluated forward here so
    # the panel and the sizing solver cannot disagree about the margin.
    resisting = (dynamics.spring_force(params.x_stroke, params)
                 + params.delta_P * params.A_seat)
    f_hold = magnetics.magnetic_force_closing(gap_open, case.i_hold, params)
    p_hold = case.i_hold ** 2 * params.R_coil_20C

    margin = f_hold / resisting
    r_max = drive.max_ripple_fraction(margin) if margin >= 1.0 else 0.0

    with st.expander("峰值-保持驅動：保持相（靜態力平衡，非 ODE）", expanded=True):
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("保持電流", f"{case.i_hold * 1e3:.2f} mA",
                  delta=f"峰值 {i_peak * 1e3:.1f} mA")
        c2.metric("保持功耗", f"{p_hold * 1e3:.1f} mW",
                  delta=f"峰值 {params.V_bus * i_peak:.2f} W")
        c3.metric("保持裕度", f"{margin:.2f}×")
        c4.metric("等效 duty",
                  f"{drive.hold_duty(case.i_hold, params) * 100:.1f} %")
        c5.metric("最大容許漣波", f"{r_max * 100:.1f} %",
                  help="電流下擺超過此比例，磁力就低於阻力而掉閥。"
                       "力∝電流²，所以 2.0× 力裕度只容許 29.3% 電流下擺。")
        if f_hold < resisting:
            st.error(f"保持電流不足：磁力 {f_hold:.2f} N < 阻力 {resisting:.2f} N，"
                     "閥會在保持相掉落")
        st.caption(
            f"保持相在氣隙 {gap_open * 1e3:.2f} mm（已吸合）計算："
            f"磁力 {f_hold:.2f} N vs 阻力 {resisting:.2f} N。"
            f"上列電流與 duty 都是**週期平均值**；真正決定掉不掉閥的是漣波"
            f"下擺的瞬時最低電流：下擺幅度須小於 {r_max * 100:.1f}%"
            "（本模型不含切換模型，無法預測漣波大小，請查驅動 IC 規格）。"
            "下方四個頁籤全部是**峰值相**（ODE 以定電壓 V_bus 驅動），"
            "不代表保持相的電流與功耗。")



def main():
    st.title("電磁閥 參數輸入介面")

    # Case selector: picking a case reseeds every sidebar default, including the
    # duty point (P_up/T0/medium). Streamlit keeps a widget's old value when its
    # default changes, so each widget's key is suffixed with the case name —
    # switching cases makes them fresh widgets that adopt the new defaults.
    case_name = st.sidebar.selectbox(
        "設計案例", list(cases.CASES),
        help="切換案例會一併帶入該案例的閥件參數與工況（P_up／T₀／流體）。"
             "選定後仍可在下方各欄位手動微調")
    case = cases.CASES[case_name]
    b = case.params
    k = _key_suffix(case_name)
    if case.notes:
        st.sidebar.caption(case.notes)

    _SPEC_DIR = _REPO_ROOT / "docs" / "spec"
    with st.sidebar.expander("參數示意圖", expanded=False):
        st.image(str(_SPEC_DIR / "valve_schematic.png"),
                 caption="整閥剖面（詳見 docs/spec/PARAMS.md）")

    with st.sidebar.expander("電氣", expanded=True):
        V_bus = st.number_input("供電電壓 (V)", value=b.V_bus, key=f"V_bus_{k}", help="電磁閥線圈的驅動電壓，典型衛星匯流排電壓 22–36V")
        R_coil_20C = st.number_input("線圈電阻 (Ω，20°C)", value=b.R_coil_20C, key=f"R_coil_{k}", help="常溫下的線圈電阻，影響穩態電流大小與電氣時間常數 τ_e=L/R")
        N_turns = st.number_input("線圈匝數", value=b.N_turns, key=f"N_turns_{k}", help="線圈繞線圈數，越多電感與磁力越大，但電氣時間常數也隨之變長")

    with st.sidebar.expander("磁路", expanded=True):
        A_gap = st.number_input("氣隙截面積 (m²)", value=b.A_gap, format="%.2e", key=f"A_gap_{k}", help="銜鐵與極面之間氣隙的有效截面積，影響磁阻與吸力大小")
        g0 = st.number_input("靜止氣隙 (m)", value=b.g0, format="%.2e", key=f"g0_{k}", help="未通電（銜鐵在靜止位置）時的氣隙距離")
        x_stroke = st.number_input("行程 (m)", value=b.x_stroke, format="%.2e", key=f"x_stroke_{k}", help="銜鐵從靜止位置到完全吸合所需移動的距離")
        l_core = st.number_input("鐵芯磁路長度 (m)", value=b.l_core, key=f"l_core_{k}", help="磁通在鐵芯內部行進的路徑長度，用於計算鐵芯磁阻")
        mu_r_core = st.number_input("鐵芯相對導磁率", value=b.mu_r_core, key=f"mu_r_{k}", help="數值越大，鐵芯磁阻相對氣隙磁阻越可忽略（本模型假設固定值，不做飽和非線性）")
        B_sat = st.number_input("飽和磁通密度 (T)", value=b.B_sat, key=f"B_sat_{k}", help="鐵芯材料的磁飽和上限，模擬中磁通密度若超過此值只會顯示警告，不會改變計算結果")

    with st.sidebar.expander("機械", expanded=True):
        m_arm = st.number_input("銜鐵質量 (kg)", value=b.m_arm, format="%.4f", key=f"m_arm_{k}", help="運動部件（銜鐵）的質量，影響加速度與響應速度")
        k_spring = st.number_input("彈簧剛性 (N/m)", value=b.k_spring, key=f"k_spring_{k}", help="彈簧的線性勁度係數，數值越大彈簧越硬")
        F_preload = st.number_input("彈簧預載力 (N)", value=b.F_preload, key=f"F_preload_{k}", help="銜鐵在靜止位置（x=0）時彈簧已經施加的作用力")

    with st.sidebar.expander("簡化負載", expanded=True):
        delta_P = st.number_input("壓力差 (Pa)", value=b.delta_P, format="%.2e", key=f"delta_P_{k}", help="跨閥座的靜態壓力負載，視為固定值，非完整流體域模擬")
        A_seat = st.number_input("閥座面積 (m²)", value=b.A_seat, format="%.2e", key=f"A_seat_{k}", help="用來把壓力差換算成靜態作用力（壓力力 = delta_P × A_seat）")
        damping_coeff = st.number_input("阻尼係數 (N·s/m)", value=b.damping_coeff, key=f"damping_{k}", help="簡化的黏滯阻尼項，非第一性推導，僅為模型簡化假設")
        st.caption("⚠️ 阻尼為簡化項，非第一性推導")

    with st.sidebar.expander("流道幾何", expanded=True):
        D_seat_bore = st.number_input("閥座流道孔徑 (m)", value=b.D_seat_bore, format="%.2e", key=f"D_bore_{k}", help="決定流量的孔徑；A_eff(x)=min(π·D·x, π·D²/4)。基準值 0.52 mm 由「全開 Cv≈0.01」反推定案。與 A_seat（密封環受壓面積）是不同物理量")
        C_d = st.number_input("流量係數 C_d", value=b.C_d, key=f"C_d_{k}", help="⚠️ 經驗係數（銳緣孔口典型 0.6–0.9），非第一性推導，精度約 ±15%。不要與 Cv（閥容量係數）混淆")

    with st.sidebar.expander("流體耦合", expanded=True):
        fluid_enabled = st.checkbox("啟用流體耦合", value=True, help="啟用後：流體力進入銜鐵運動方程式，且開啟動態／極限掃描／密封壽命三頁籤的靜態壓力負載一律改由下方工況（P_up−P_down）推導，取代上方「壓力差」欄位。關閉＝乾跑模式（P0 行為），此時上方「壓力差」欄位才是實際生效值")
        medium_name = st.selectbox("流體", ["N₂", "Xe", "He", "水 (20°C)", "LN₂ (77K)", "自訂氣體", "自訂液體"],
                                   index=_medium_index(case.medium), key=f"medium_{k}", help="內建物性常數見 docs/spec/PARAMS.md；自訂時輸入物性")
        if medium_name == "自訂氣體":
            gamma_in = st.number_input("比熱比 γ", value=1.4)
            R_in = st.number_input("比氣體常數 R (J/kg·K)", value=296.8)
            medium = fluid.GasProperties(gamma=gamma_in, R_specific=R_in)
        elif medium_name == "自訂液體":
            rho_in = st.number_input("密度 ρ (kg/m³)", value=998.2)
            c_in = st.number_input("音速 c (m/s)", value=1482.0)
            pvap_in = st.number_input("蒸氣壓 P_vap (Pa)", value=2339.0)
            medium = fluid.LiquidProperties(rho=rho_in, c_sound=c_in, P_vap=pvap_in)
        else:
            medium = {"N₂": fluid.N2, "Xe": fluid.XE, "He": fluid.HE,
                      "水 (20°C)": fluid.WATER_20C, "LN₂ (77K)": fluid.LN2_77K}[medium_name]
        P_up = st.number_input("上游壓力 P_up (Pa)", value=case.cond.P_up, format="%.2e", key=f"P_up_{k}", help="上游滯止壓力（絕對壓）。預設＝基準 MEOP 2.4 MPa")
        P_down = st.number_input("下游壓力 P_down (Pa)", value=case.cond.P_down, format="%.2e", key=f"P_down_{k}", help="下游壓力（絕對壓）。預設 0＝排真空")
        T0 = st.number_input("上游溫度 T₀ (K)", value=case.cond.T0, key=f"T0_{k}", help="氣體分支使用的滯止溫度；液體分支不使用")
        cond = fluid.FlowConditions(P_up=P_up, P_down=P_down, T0=T0)

    with st.sidebar.expander("熱域", expanded=True):
        T_coil = st.number_input(
            "線圈溫度 T_coil (K)", value=case.T_coil, min_value=77.0, max_value=500.0,
            key=f"T_coil_{k}",
            help="等溫假設：單次開啟（~5 ms）遠短於線圈熱時間常數（秒級），"
                 "溫度視為輸入參數。對照：233.15 K = −40°C、293.15 K = 20°C、"
                 "343.15 K = +70°C、77 K = LN₂。77 K 以下超出銅電阻表域。"
                 "預設 293.15 K 時熱態電阻恰等於 20°C 電阻，結果與未耦合完全相同")

    with st.sidebar.expander("繞線窗口", expanded=False):
        st.image(str(_SPEC_DIR / "winding_window_schematic.png"))
        A_winding = st.number_input("繞線窗口面積 A_winding (m²)",
                                    value=b.A_winding, format="%.4e", key=f"A_win_{k}",
                                    help="鐵芯窗口可繞線的截面積，決定 R=ρ·N²/(A_w·k_fill)·l_turn")
        k_fill = st.number_input("銅填充率 k_fill", value=b.k_fill, key=f"k_fill_{k}",
                                 help="銅截面佔窗口面積比例，圓線+絕緣典型 0.4–0.6")
        l_turn_mean = st.number_input("平均匝長 l_turn_mean (m)",
                                      value=b.l_turn_mean, format="%.4e", key=f"l_turn_{k}",
                                      help="一匝的平均周長")
        R_derived = winding.coil_resistance(N_turns, A_winding, k_fill,
                                            l_turn_mean)
        R_nominal = R_coil_20C
        rel_err = (R_derived - R_nominal) / R_nominal if R_nominal else 0.0
        if abs(rel_err) > 0.05:
            st.warning(f"⚠️ 繞線窗口推導電阻 {R_derived:.1f} Ω 與輸入的 "
                       f"R_coil_20C={R_nominal:.1f} Ω 不一致（差 {rel_err:+.0%}）")
        else:
            st.caption(f"窗口推導 R = {R_derived:.2f} Ω（與 R_coil_20C 一致）")

    with st.sidebar.expander("閥座/密封", expanded=False):
        st.image(str(_SPEC_DIR / "seal_land_schematic.png"))
        seat_label = st.selectbox("閥座材料對",
                                  ["PCTFE/440C", "17-4PH/440C", "自訂"],
                                  help="決定接觸力學（E*/H/Rq_c）與恢復係數 e；內建組合見 materials.py")
        if seat_label == "自訂":
            mat_names = list(materials.MATERIALS)
            mat_a = st.selectbox("閥座材料 A", mat_names, index=mat_names.index("PCTFE"))
            mat_b = st.selectbox("閥座材料 B", mat_names, index=mat_names.index("440C"))
            e_rest = st.number_input("恢復係數 e", value=0.4, min_value=0.0,
                                     max_value=1.0,
                                     help="撞擊恢復係數（0=完全塑性、1=完全彈性），EMPIRICAL")
            seat_mats = (materials.MATERIALS[mat_a], materials.MATERIALS[mat_b])
            seat = materials.seat_pair(*seat_mats, e_rest)
        else:
            seat, seat_mats = {
                "PCTFE/440C": (materials.PCTFE_ON_440C,
                               (materials.PCTFE, materials.SS_440C)),
                "17-4PH/440C": (materials.METAL_17_4PH_ON_440C,
                                (materials.SS_17_4PH, materials.SS_440C)),
            }[seat_label]
        w_land = st.number_input("密封 land 寬度 w_land (m)", value=b.w_land,
                                 format="%.2e",
                                 help="⚠️ EMPIRICAL 幾何參數；密封門檻與壽命深度基準都以它為尺")
        R_tip = st.number_input("poppet 端面曲率半徑 R_tip (m)", value=b.R_tip,
                                format="%.2e", help="⚠️ EMPIRICAL；Hertz 撞擊接觸的等效曲率半徑")

    params = ValveParams(
        V_bus=V_bus, R_coil_20C=R_coil_20C, N_turns=N_turns,
        A_gap=A_gap, g0=g0, x_stroke=x_stroke, l_core=l_core,
        mu_r_core=mu_r_core, B_sat=B_sat,
        m_arm=m_arm, k_spring=k_spring, F_preload=F_preload,
        delta_P=delta_P, A_seat=A_seat, damping_coeff=damping_coeff,
        D_seat_bore=D_seat_bore, C_d=C_d,
        A_winding=A_winding, k_fill=k_fill, l_turn_mean=l_turn_mean,
        w_land=w_land, R_tip=R_tip,
    )

    # Peak-and-hold cases: the ODE drives a constant V_bus, so every tab below
    # shows the PEAK phase. Report the hold phase separately rather than let the
    # peak numbers read as the whole duty cycle.
    if case.i_hold is not None and abs(params.V_bus - b.V_bus) < 1e-9:
        _render_hold_phase(params, case)

    tab_open, tab_close, tab_limits, tab_seal, tab_layout = st.tabs(
        ["開啟動態", "關閉/續流", "極限掃描 L1–L3", "密封/壽命 L4–L5",
         "實體佈局"])
    with tab_open:
        try:
            render_opening_tab(params, medium, cond, T_coil, fluid_enabled, seat)
        except ValueError as e:
            st.error(f"開啟動態計算失敗：{e}")
    with tab_close:
        try:
            render_closing_tab(params, seat)
        except ValueError as e:
            st.error(f"關閉/續流計算失敗：{e}")
    with tab_limits:
        try:
            # Dual-mode rule (ARCHITECTURE.md): in fluid-coupled mode the
            # operating pressure comes from `cond`, not `params.delta_P`.
            # `limits.response_time_bound`/`dimensionless_groups` take no
            # `cond` argument, so without this they'd silently keep using
            # stale `delta_P` while Tab 1's t_open already reflects `cond`
            # -- producing an L1 "lower bound" that isn't actually a bound
            # on the displayed t_open whenever P_up-P_down != delta_P.
            limits_params = (replace(params, delta_P=cond.P_up - cond.P_down)
                             if fluid_enabled else params)
            render_limits_tab(limits_params)
        except ValueError as e:
            st.error(f"極限掃描計算失敗：{e}")
    with tab_seal:
        try:
            render_sealing_tab(params, seat, seat_mats, medium, cond)
        except ValueError as e:
            st.error(f"密封/壽命計算失敗：{e}")
    with tab_layout:
        try:
            render_layout_tab(params)
        except ValueError as e:
            st.error(f"佈局計算失敗：{e}")


def render_opening_tab(params, medium, cond, T_coil, fluid_enabled, seat):
    result_dry = compute_result(params, T_coil=T_coil)
    if fluid_enabled:
        result = compute_result(params, medium=medium, cond=cond, T_coil=T_coil)
    else:
        result = result_dry

    if result is None:
        st.error("此參數組合下電磁力不足以克服彈簧/壓力負載，閥件不會吸合。")
        return

    if fluid_enabled and result_dry is not None:
        delta_ms = (result.t_open - result_dry.t_open) * 1e3
        st.metric("t_open", f"{result.t_open * 1e3:.3f} ms",
                  delta=f"{delta_ms:+.3f} ms vs 乾跑", delta_color="inverse")
    else:
        st.metric("t_open", f"{result.t_open * 1e3:.3f} ms")

    if fluid_enabled:
        col1, col2 = st.columns(2)
        col1.metric("全開 Cv", f"{fluid.cv_from_geometry(params.x_stroke, params):.4f}")
        if isinstance(medium, fluid.GasProperties):
            mdot_full = fluid.mdot_gas(params.x_stroke, params, medium, cond)
        else:
            mdot_full = fluid.mdot_liquid(params.x_stroke, params, medium, cond)
        col2.metric("全開質量流率", f"{mdot_full * 1e3:.3f} g/s")
        if isinstance(medium, fluid.LiquidProperties) and fluid.flashing_risk(medium, cond):
            st.warning("⚠️ 閃蒸風險：下游壓力 ≤ 蒸氣壓，流體會閃蒸（兩相流）——本模型液體公式在此工況失效，動態結果不可信。")

    R_hot = thermal.R_coil(T_coil, params)
    col_r, col_v = st.columns(2)
    col_r.metric("熱態線圈電阻", f"{R_hot:.2f} Ω",
                 delta=f"{R_hot - params.R_coil_20C:+.2f} Ω vs 20°C")
    try:
        I_ss = params.V_bus / R_hot
        T_eq = thermal.equilibrium_temp(I_ss, T_coil, params)
        st.caption(f"持續通電保持（定電流 {I_ss:.3f} A 保守假設，真空僅傳導+輻射散熱）："
                   f"平衡線圈溫度 ≈ {T_eq:.1f} K（溫升 +{T_eq - T_coil:.1f} K）")
    except ValueError:
        st.error("熱失控：此參數下持續通電無穩態平衡（檢查 G_th_cond / emissivity / A_rad）")
    if st.button("計算拉入電壓（數秒）", key="btn_pullin"):
        try:
            with st.spinner("二分搜尋拉入電壓中…"):
                v_pi = limits.find_pullin_voltage(params, T_coil=T_coil,
                                                  tol=0.1, t_max=0.03)
        except ValueError as e:
            st.error(f"拉入電壓搜尋失敗：{e}")
        else:
            col_v.metric("V_pull-in", f"{v_pi:.2f} V")
            if v_pi > params.V_bus:
                st.warning("⚠️ 拉入電壓高於供電電壓——超出工作包絡，此溫度下閥不會吸合")
            else:
                st.caption(f"拉入裕度 {params.V_bus - v_pi:.2f} V（對當前供電電壓）")

    if result.saturated:
        st.warning("模擬中磁通密度曾超過飽和值 B_sat（模型未做非線性 B-H 曲線，此結果僅供警示）")

    bounce_on = st.checkbox("啟用閥座回彈", value=False,
                            help="傳入閥座材料對（側欄），銜鐵到達止擋後以恢復係數 e 回彈積分；"
                                 "關閉＝P0–P5 硬鉗制行為。頭條 t_open 為首次到達止擋時間，不受影響")
    if bounce_on:
        br = compute_opening_bounce(params, seat,
                                    medium=medium if fluid_enabled else None,
                                    cond=cond if fluid_enabled else None,
                                    T_coil=T_coil)
        c1, c2, c3 = st.columns(3)
        c1.metric("回彈次數", f"{br.n_bounce}")
        c2.metric("首次觸擊速度", f"{br.v_impacts[0]:.3f} m/s")
        c3.metric("末次觸擊速度", f"{br.v_impacts[-1]:.3f} m/s")
        st.caption(f"t_settle = {br.t_settle * 1e3:.3f} ms（通電起算）。"
                   "首飛期間電流仍在爬升，撞速比從上方收斂到 e（P6 發現 6）")

    sol = result.sol
    t_ms = sol.t * 1e3
    i_t, x_t, _v_t = sol.y

    fig_i, ax_i = plt.subplots()
    ax_i.plot(t_ms, i_t)
    ax_i.set_xlabel("time (ms)")
    ax_i.set_ylabel("coil current (A)")
    ax_i.set_title("i(t)")
    st.pyplot(fig_i)
    plt.close(fig_i)

    fig_x, ax_x = plt.subplots()
    ax_x.plot(t_ms, x_t * 1e3)
    ax_x.set_xlabel("time (ms)")
    ax_x.set_ylabel("armature displacement (mm)")
    ax_x.set_title("x(t)")
    st.pyplot(fig_x)
    plt.close(fig_x)

    F_mag_t = [
        magnetics.magnetic_force_closing(params.g0 - x, i, params)
        for i, x in zip(i_t, x_t)
    ]
    fig_f, ax_f = plt.subplots()
    ax_f.plot(t_ms, F_mag_t)
    ax_f.set_xlabel("time (ms)")
    ax_f.set_ylabel("magnetic force (N)")
    ax_f.set_title("F_mag(t)")
    st.pyplot(fig_f)
    plt.close(fig_f)


def render_closing_tab(params, seat):
    circ_name = st.selectbox("續流電路", ["二極體", "稽納", "RC 緩衝"],
                             help="斷電後線圈電流的洩放路徑，決定關閉速度與電壓應力")
    if circ_name == "二極體":
        V_fwd = st.number_input("順向壓降 V_fwd (V)", value=0.7,
                                help="⚠️ 簡化為定值（實際隨電流呈對數變化），非第一性")
        circuit = flyback.DiodeFlyback(V_fwd=V_fwd)
    elif circ_name == "稽納":
        V_clamp = st.number_input("鉗位電壓 V_clamp (V)", value=30.0,
                                  help="稽納鉗位電壓越高，電流洩放越快、關閉越快，"
                                       "但觸座速度也越高（P6 發現 3）")
        circuit = flyback.ZenerFlyback(V_clamp=V_clamp)
    else:
        st.caption("R_snub/C_snub 刻意無預設：欠阻尼組合不是合理的 snubber 設計，"
                   "阻尼比公式見 docs/spec/PARAMS.md")
        R_snub = st.number_input("緩衝電阻 R_snub (Ω)", value=0.0, min_value=0.0)
        C_snub = st.number_input("緩衝電容 C_snub (F)", value=0.0,
                                 min_value=0.0, format="%.2e")
        if R_snub <= 0.0 or C_snub <= 0.0:
            st.warning("請輸入正的 R_snub 與 C_snub 才能模擬 RC 緩衝電路")
            return
        circuit = flyback.RCFlyback(R_snub=R_snub, C_snub=C_snub)

    compare = st.checkbox("與二極體對照",
                          help="加跑一次二極體續流，顯示關閉時間比與觸座速度比")
    bounce_on = st.checkbox("閥座回彈", key="cb_close_bounce",
                            help="關閉觸座後以恢復係數 e 回彈積分（側欄閥座材料對）")
    if st.button("模擬關閉", key="btn_closing"):
        summary = compute_closing(params, circuit,
                                  seat=seat if bounce_on else None)
        ref = compute_closing(params, flyback.DiodeFlyback()) if compare else None
        st.session_state["closing"] = (params, summary)
        st.session_state["closing_ref"] = ref
    if "closing" not in st.session_state:
        st.info("按「模擬關閉」執行（斷電 → 止擋保持期 → 釋放衝程期）")
        return
    stored_params, summary = st.session_state["closing"]
    if stored_params != params:
        st.caption("⚠️ 顯示的結果對應舊參數——請重按「模擬關閉」")
    ref = st.session_state.get("closing_ref")
    c1, c2, c3 = st.columns(3)
    c1.metric("t_release", f"{summary.t_release * 1e3:.3f} ms")
    c2.metric("t_close", f"{summary.t_close * 1e3:.3f} ms")
    c3.metric("觸座速度", f"{summary.v_impact:.4f} m/s")
    if ref is not None:
        st.caption(f"對照二極體：關閉時間 {ref.t_close / summary.t_close:.2f} 倍快、"
                   f"觸座速度 {summary.v_impact / ref.v_impact:.2f} 倍高"
                   "（稽納加速的代價，P6 發現 3）")

    r = summary.result
    fig_i, ax_i = plt.subplots()
    ax_i.plot(r.sol_hold.t * 1e3, r.sol_hold.y[0], label="hold")
    ax_i.plot((r.t_release + r.sol_stroke.t) * 1e3, r.sol_stroke.y[0],
              label="stroke")
    ax_i.set_xlabel("time since de-energization (ms)")
    ax_i.set_ylabel("coil current (A)")
    ax_i.set_title("i(t) decay")
    ax_i.legend()
    st.pyplot(fig_i)
    plt.close(fig_i)

    fig_x, ax_x = plt.subplots()
    ax_x.plot((r.t_release + r.sol_stroke.t) * 1e3, r.sol_stroke.y[1] * 1e3,
              label="stroke")
    if r.bounce is not None:
        for seg in r.bounce.segments[1:]:
            ax_x.plot((r.t_release + seg.t) * 1e3, seg.y[1] * 1e3,
                      color="tab:red")
        st.caption(f"回彈 {r.bounce.n_bounce} 次，首次觸座速度 "
                   f"{r.bounce.v_impacts[0]:.3f} m/s")
    ax_x.set_xlabel("time since de-energization (ms)")
    ax_x.set_ylabel("poppet position (mm)")
    ax_x.set_title("x(t) closing stroke")
    st.pyplot(fig_x)
    plt.close(fig_x)


if __name__ == "__main__":
    main()
