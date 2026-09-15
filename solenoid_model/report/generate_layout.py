"""按比例剖面圖、作動前後對照圖與架構圖（P8）。

剖面圖用真實 mm 比例（set_aspect("equal")），與 docs/spec/valve_schematic.png
的符號示意圖不同——後者是參數對照用、不按比例。

作動對照圖的上排為真實比例，下排放大圖只把「導杆＋銜鐵中段」摺疊
（不按比例）以突顯 0.05-0.20 mm 的間隙；間隙本身仍按 1:1 且兩狀態共用
刻度，所以圖上的長度比等於真實比例。

拓樸依 valve_schematic.png：閥座 → 閥芯 → 導杆 → 銜鐵 → 工作氣隙 →
固定極。閥芯／導杆／銜鐵是同一動件組，彈簧套在導杆上、夾在固定彈簧座
與閥芯頂面之間，故動件上移 x 即壓縮 x（見 state_geometry）。
"""
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

from solenoid_model import fluid as fluid_lib
from solenoid_model import layout
from solenoid_model.cases import N2_25BAR_PH_PARAMS

_DOCS = Path(__file__).resolve().parents[2] / "docs" / "spec"
SECTION_OUT = _DOCS / "layout_section.png"
ARCH_OUT = _DOCS / "layout_architecture.png"
ACTUATION_OUT = _DOCS / "layout_actuation.png"

COL_MAG = "#4477aa"     # 磁性件
COL_NONMAG = "#bbbbbb"  # 非磁性件
COL_COIL = "#f4a261"    # 線圈
COL_SEAL = "#cc3311"    # 密封/流體
COL_MECH = "#228833"    # 機械

# 比照 generate_schematic.py：以 rc_context 區域套用，不改全域 rcParams
_CJK_FONT_RC = {
    "font.sans-serif": ["Arial Unicode MS", "PingFang TC", "Heiti TC"],
    "axes.unicode_minus": False,
}

# 外殼內徑：封裝外形輸入（envelope 的 coil_OD 參數沿用此名），不是繞線
# 窗口本身。線圈斷面一律由 params.A_winding / params.l_turn_mean 推導，
# 見 coil_box()，故此處不再宣告 COIL_ID / COIL_H 常數——先前版本的
# 18/6/12 mm 與 params 數值一致純屬巧合（A_winding 恰為 72 mm^2），
# 一旦 params 變動圖面便會無聲說謊。
COIL_OD = 18.0e-3


def _mm(x):
    return x * 1e3


@dataclass(frozen=True)
class CoilBox:
    """線圈繞線斷面（半剖，單位 m）。r_in/r_out 為半徑，y_lo/y_hi 為軸向。"""
    r_in: float
    r_out: float
    y_lo: float
    y_hi: float

    @property
    def width(self):
        return self.r_out - self.r_in

    @property
    def height(self):
        return self.y_hi - self.y_lo

    @property
    def area(self):
        """繞線窗口面積 = 全寬 x 高（全寬 = 2 x 半剖寬度）。"""
        return 2.0 * self.width * self.height


def coil_box(params, y_win_lo, y_win_hi, t_sleeve, R_shell):
    """由 params 推導線圈斷面，而非以常數硬寫。

    內緣貼在隔離套外側（r_core + t_sleeve），中間不留不明空隙。軸向填滿
    可用窗口（閥座頂面到彈簧腔底面）——線圈不能與座體或彈簧腔爭同一空間。
    寬度再由 A_winding = 2*w*h 反解，故繞線面積與 ValveParams 恆等：
    params 一變，圖面跟著變，不會像常數版（18/6/12 mm）在 A_winding
    改動後無聲說謊。

    先前版本改以 r_mean 定寬、置中於鐵芯跨距，結果線圈底緣侵入座體
    0.787 mm——同一空間被兩個零件佔用。以窗口為準即無此衝突。
    """
    r_in = core_radius_with_sleeve(params, t_sleeve)
    h = y_win_hi - y_win_lo
    if h <= 0:
        raise ValueError(f"線圈軸向窗口 {h * 1e3:.3f} mm 非正：軸向鏈不相容")
    w = params.A_winding / (2.0 * h)
    if r_in + w > R_shell:
        raise ValueError(
            f"繞線窗口 A_winding={params.A_winding * 1e6:.1f} mm² 在高 "
            f"{h * 1e3:.2f} mm 下需寬 {w * 1e3:.2f} mm，外緣 "
            f"{(r_in + w) * 1e3:.2f} mm 超出外殼內徑 {R_shell * 1e3:.2f} mm")
    return CoilBox(r_in=r_in, r_out=r_in + w, y_lo=y_win_lo, y_hi=y_win_hi)


def core_radius_with_sleeve(params, t_sleeve):
    """線圈內緣半徑：鐵芯外加非磁性隔離套，中間不留不明空隙。"""
    return layout.core_radius(params) + t_sleeve


def spring_installed_length(params, spring):
    """彈簧「安裝長」：閥關（銜鐵坐封）時的實際長度 [m]。

    spring.L_free 是自由長——未裝入時的長度。裝配後彈簧永遠至少被預載
    壓縮 F_preload/k_spring，閥關時即為此狀態；閥全開再多壓 x_stroke。
    axial_stack 的「彈簧腔」用 L_free 是腔體尺寸（組裝時須容得下），
    但剖面圖畫的是已裝配狀態，兩者不同，故此處另算不動 layout.py。
    """
    return spring.L_free - params.F_preload / params.k_spring


def axial_parts(params, t_seat, t_pole, t_plate, spring):
    """軸向零件鏈 [(name, y_lo, y_hi)]，座面朝上，含兩片端板。

    順序與厚度對齊 layout.axial_stack 的 segments（端板 x2），因此圖面
    總長恆等於 envelope()['L']；先前版本漏畫座端端板，標註 18.679 mm
    卻只畫到 17.417 mm。

    段序依 valve_schematic.png：閥座在下，往上是行程、彈簧腔、銜鐵、
    工作氣隙、固定極（線圈在其外側）。彈簧在銜鐵「下方」是力學要求，
    見 layout.axial_stack 的說明。
    """
    parts, y = [], 0.0
    t_pop = layout.LAYOUT_EMPIRICAL["t_poppet"]["value"]
    for name, t in (("端板（座端）", t_plate),
                    ("閥座座體", t_seat),
                    ("行程", params.x_stroke),
                    ("閥芯", t_pop),
                    ("彈簧腔", spring.L_free),
                    ("銜鐵", layout.armature_thickness(params)),
                    ("工作氣隙", params.g0),
                    ("固定極", t_pole),
                    ("端板（極端）", t_plate)):
        parts.append((name, y, y + t))
        y += t
    return parts


def section_geometry(params=N2_25BAR_PH_PARAMS):
    """剖面圖的全部幾何 [m]，與繪圖共用同一組數字。

    抽成獨立函式讓測試能直接斷言幾何（總長是否等於 envelope L、線圈
    面積是否等於 A_winding、零件是否互相穿插），而不是只驗 PNG 檔案
    大小——舊測試僅斷言 st_size > 1000，六個幾何錯誤因此全數無聲通過。
    """
    env = layout.envelope(params, coil_OD=COIL_OD)
    spring = layout.spring_geometry(params, **layout.SPRING_DESIGN)
    t_plate = layout.end_plate_thickness(params)
    t_pole = layout.LAYOUT_EMPIRICAL["t_fixed_pole"]["value"]
    t_seat = layout.LAYOUT_EMPIRICAL["t_seat_body"]["value"]
    t_sleeve = layout.LAYOUT_EMPIRICAL["t_sleeve"]["value"]
    parts = axial_parts(params, t_seat, t_pole, t_plate, spring)
    span = dict((n, (lo, hi)) for n, lo, hi in parts)
    # 線圈軸向窗口：繞線必須包住磁路的可動段，即銜鐵底面到固定極頂面。
    # 線圈不在閥座側（那裡是彈簧與閥芯，非磁路），也不侵入端板。
    coil = coil_box(params, span["銜鐵"][0], span["固定極"][1],
                    t_sleeve, COIL_OD / 2.0)
    return {
        "env": env, "spring": spring, "parts": parts, "span": span,
        "coil": coil, "t_sleeve": t_sleeve, "t_plate": t_plate,
        "R_shell": COIL_OD / 2.0, "wall": env["wall"].adopted,
        "r_core": layout.core_radius(params),
        "L_total": parts[-1][2],
        "spring_installed": spring_installed_length(params, spring),
        "r_spring_out": (spring.D_coil + spring.d_wire) / 2.0,
    }


def state_geometry(params=N2_25BAR_PH_PARAMS, x=0.0):
    """作動狀態幾何 [m]：銜鐵位移 x 時各界面的位置。

    遵循 dynamics 的唯一約定 gap = params.g0 - x（見 dynamics.coupled_rhs
    與 simulate_closing 的 "armature returning to x=0"）：

      x = 0          閥關：閥芯坐封，氣隙最大 = g0
      x = x_stroke   閥開：動件組被吸向固定極，氣隙最小 = g0 - x_stroke

    開度與工作氣隙是同一動件組的兩端，此消彼長，和恆為
    x_stroke + (g0 - x_stroke) = g0。早期剖面圖把兩者同時畫成最大，
    是物理上不可能的狀態（動件不能同時離座又離極最遠）。

    動件組＝閥芯＋導杆＋銜鐵，整體平移 x：密封由下方閥芯完成，磁吸由
    上方銜鐵承受。彈簧套在導杆上，上端頂固定彈簧座（不動）、下端頂閥芯
    頂面（隨 x 上移），因此壓縮量恰為 x，F = F_preload + k*x 與
    dynamics.spring_force 恆等。
    """
    if not 0.0 <= x <= params.x_stroke + 1e-15:
        raise ValueError(
            f"x={x * 1e3:.3f} mm 超出行程 0..{params.x_stroke * 1e3:.3f} mm")
    G = section_geometry(params)
    span = G["span"]
    y_seat_top = span["閥座座體"][1]
    t_arm = layout.armature_thickness(params)
    t_pop = layout.LAYOUT_EMPIRICAL["t_poppet"]["value"]
    # 動件組（閥芯＋導杆＋銜鐵）整體平移 x：閥芯底面離座 x，銜鐵頂面
    # 逼近固定極，故 gap = g0 - x。固定極與彈簧座都不動。
    # 閥關(x=0)時閥芯坐在座面上；行程是它上方的淨空，隨 x 被吃掉。
    # 用 span["閥芯"][0] 當基準會讓 x=0 的閥芯浮在座上 0.15 mm——
    # 尺寸鏈把「行程」排在閥芯之下（組裝淨空），但實體位置以座面為準。
    y_pop_lo = y_seat_top + x
    y_pop_hi = y_pop_lo + t_pop
    y_arm_hi = span["固定極"][0] - (params.g0 - x)
    y_arm_lo = y_arm_hi - t_arm
    # 彈簧：上端頂在固定彈簧座（不動），下端頂在閥芯頂面（隨 x 上移），
    # 故彈簧長度 = 安裝長 - x，壓縮量隨行程增加 → F = F_preload + k*x
    L_spring = G["spring_installed"] - x
    y_spring_lo = y_pop_hi
    y_spring_hi = y_spring_lo + L_spring
    # 固定彈簧座：貼在彈簧腔頂（= 銜鐵底面）之下，不動。錨定在腔頂而非
    # 彈簧頂端，否則它會隨行程移動（彈簧的上支承必須是固定面），也不會
    # 像先前那樣長進銜鐵與線圈裡。
    t_seat_spr = layout.LAYOUT_EMPIRICAL["t_spring_seat"]["value"]
    y_seat_spr_hi = span["彈簧腔"][1]
    y_seat_spr_lo = y_seat_spr_hi - t_seat_spr
    return {
        "base": G, "x": x, "gap": params.g0 - x,
        "seat_gap": x,
        "y_seat_top": y_seat_top,
        "y_pop_lo": y_pop_lo, "y_pop_hi": y_pop_hi,
        "y_arm_lo": y_arm_lo, "y_arm_hi": y_arm_hi,
        "y_pole_lo": span["固定極"][0],
        "y_spring_lo": y_spring_lo, "y_spring_hi": y_spring_hi,
        "L_spring": L_spring,
        "y_seat_spr_lo": y_seat_spr_lo, "y_seat_spr_hi": y_seat_spr_hi,
        "t_spring_seat": t_seat_spr,
        "F_spring": params.F_preload + params.k_spring * x,
        "energised": x > 0.0,
    }


def generate_section(params=N2_25BAR_PH_PARAMS):
    """按比例軸對稱半剖圖，單位 mm（閥關狀態）。"""
    G = section_geometry(params)
    S = state_geometry(params, x=0.0)
    env, spring, span, coil = G["env"], G["spring"], G["span"], G["coil"]
    wall, R_shell, r_core = G["wall"], G["R_shell"], G["r_core"]
    L_total = G["L_total"]
    r_stem = layout.LAYOUT_EMPIRICAL["d_stem"]["value"] / 2.0
    r_pop = layout.LAYOUT_EMPIRICAL["d_poppet"]["value"] / 2.0
    r_spr = G["r_spring_out"]

    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, ax = plt.subplots(figsize=(8.4, 10.5))
        x_lead = _mm(R_shell + wall) + 2.4
        R_out = R_shell + wall

        def yl(n):
            return span[n][0]

        def yh(n):
            return span[n][1]

        # --- 不動件：端板 x2、外殼
        for name in ("端板（座端）", "端板（極端）"):
            ax.add_patch(Rectangle((0, _mm(yl(name))), _mm(R_out),
                                   _mm(yh(name) - yl(name)),
                                   fc=COL_MAG, ec="k"))
        ax.add_patch(Rectangle((_mm(R_shell), 0), _mm(wall), _mm(L_total),
                               fc=COL_MAG, ec="k"))
        ax.text(_mm(R_out) + 0.6, _mm(L_total) / 2,
                f"外殼 {_mm(wall):.2f} mm", fontsize=8, color=COL_MAG,
                rotation=90, va="center")
        # --- 閥座座體 + 流道孔
        ax.add_patch(Rectangle((0, _mm(yl("閥座座體"))), _mm(R_out),
                               _mm(yh("閥座座體") - yl("閥座座體")),
                               fc=COL_NONMAG, ec="k"))
        ax.text(_mm(R_shell) * 0.72, _mm(sum(span["閥座座體"]) / 2),
                "閥座座體", ha="center", va="center", fontsize=8)
        r_bore = _mm(params.D_seat_bore / 2)
        ax.add_patch(Rectangle((0, 0), r_bore, _mm(yh("閥座座體")),
                               fc="w", ec="none", zorder=2))
        ax.plot([0, r_bore, r_bore], [0, 0, _mm(yh("閥座座體"))],
                color=COL_SEAL, lw=1.5, zorder=3)
        y_bore_label = _mm(yh("閥座座體") / 2)
        ax.plot([r_bore, x_lead], [y_bore_label, y_bore_label],
                color=COL_SEAL, lw=0.9)
        ax.text(x_lead, y_bore_label,
                f" 流道孔徑 D_seat_bore = {_mm(params.D_seat_bore):.2f} mm",
                fontsize=7.5, color=COL_SEAL, ha="left", va="center")
        # --- 固定彈簧座：彈簧的上支承（不動），套在導杆外
        ax.add_patch(Rectangle((_mm(r_stem), _mm(S["y_seat_spr_lo"])),
                               _mm(R_out - r_stem), _mm(S["t_spring_seat"]),
                               fc=COL_NONMAG, ec="k"))
        ax.plot([_mm(R_out), x_lead - 0.3],
                [_mm(S["y_seat_spr_lo"]) + 0.2] * 2, color="k", lw=0.7)
        ax.text(x_lead - 0.2, _mm(S["y_seat_spr_lo"]) - 0.9,
                " 固定彈簧座（不動）：彈簧的上支承",
                fontsize=7, color="k", ha="left", va="center")
        # --- 動件組：閥芯 + 導杆 + 銜鐵（閥關狀態 x=0）
        ax.add_patch(Rectangle((0, _mm(S["y_pop_lo"])), _mm(r_pop),
                               _mm(S["y_pop_hi"] - S["y_pop_lo"]),
                               fc=COL_MAG, ec="k", lw=1.4))
        ax.text(_mm(r_pop) + 0.35, _mm((S["y_pop_lo"] + S["y_pop_hi"]) / 2),
                f"閥芯 ⌀{_mm(2 * r_pop):.2f}", fontsize=7, ha="left",
                va="center")
        ax.add_patch(Rectangle((0, _mm(S["y_pop_hi"])), _mm(r_stem),
                               _mm(S["y_arm_lo"] - S["y_pop_hi"]),
                               fc=COL_MAG, ec="k", lw=1.4))
        ax.add_patch(Rectangle((0, _mm(S["y_arm_lo"])), _mm(r_core),
                               _mm(S["y_arm_hi"] - S["y_arm_lo"]),
                               fc=COL_MAG, ec="k", lw=1.4))
        ax.text(_mm(r_core) / 2, _mm((S["y_arm_lo"] + S["y_arm_hi"]) / 2),
                "銜鐵", ha="center", va="center", fontsize=8, color="w")
        # --- 彈簧：套在導杆上，夾在固定彈簧座與閥芯頂面之間
        ax.add_patch(Rectangle((_mm(r_stem), _mm(S["y_spring_lo"])),
                               _mm(r_spr - r_stem), _mm(S["L_spring"]),
                               fc="none", ec=COL_MECH, lw=1.6, ls="--"))
        ax.text(_mm(r_spr) + 0.3, _mm((S["y_spring_lo"] + S["y_spring_hi"]) / 2),
                f"彈簧 d{_mm(spring.d_wire):.2f}/D{_mm(spring.D_coil):.2f}\n"
                f"安裝長 {_mm(S['L_spring']):.2f}（自由長 "
                f"{_mm(spring.L_free):.2f}）\n"
                f"預載 {S['F_spring']:.2f} N ↓ 壓閥芯坐封",
                fontsize=6.5, color=COL_MECH, ha="left", va="center")
        # --- 固定極（磁路上端，不動）
        ax.add_patch(Rectangle((0, _mm(yl("固定極"))), _mm(r_core),
                               _mm(yh("固定極") - yl("固定極")),
                               fc=COL_MAG, ec="k"))
        ax.text(_mm(r_core) / 2, _mm(sum(span["固定極"]) / 2), "固定極",
                ha="center", va="center", fontsize=8, color="w")
        # --- 線圈（包住磁路）與隔離套
        ax.add_patch(Rectangle((_mm(coil.r_in), _mm(coil.y_lo)),
                               _mm(coil.width), _mm(coil.height),
                               fc=COL_COIL, ec="k"))
        ax.text(_mm((coil.r_in + coil.r_out) / 2),
                _mm((coil.y_lo + coil.y_hi) / 2),
                f"線圈\n{params.N_turns:.0f} 匝\n{coil.area * 1e6:.1f} mm²",
                ha="center", va="center", fontsize=8)
        ax.add_patch(Rectangle((_mm(coil.r_in - G["t_sleeve"]),
                                _mm(coil.y_lo)), _mm(G["t_sleeve"]),
                               _mm(coil.height),
                               fc="none", ec=COL_SEAL, lw=1.5, hatch="//"))
        # --- 兩個次毫米間隙的引線
        for name, color, fmt, y_end in (
                ("行程", COL_MECH, "行程 x_stroke = {:.2f} mm（閥芯離座）",
                 _mm(S["y_pop_lo"]) - 1.1),
                ("工作氣隙", COL_MAG, "工作氣隙 g0 = {:.2f} mm",
                 _mm(coil.y_hi) + 0.8)):
            lo, hi = span[name]
            y_mid = _mm((lo + hi) / 2)
            r_from = _mm(r_pop) if name == "行程" else _mm(r_core)
            r_stub = r_from + 0.35
            ax.plot([r_from, r_stub, r_stub, x_lead],
                    [y_mid, y_mid, y_end, y_end],
                    color=color, lw=1.0, solid_capstyle="butt")
            ax.text(x_lead, y_end, " " + fmt.format(_mm(hi - lo)),
                    fontsize=8, color=color, ha="left", va="center")
        # --- 中心線、總長
        ax.axvline(0, color="k", lw=0.8, ls="-.")
        ax.text(0.1, _mm(L_total) + 0.6, "軸心", fontsize=7)
        ax.add_patch(FancyArrowPatch((-1.2, 0), (-1.2, _mm(L_total)),
                                     arrowstyle="<->", mutation_scale=12,
                                     color=COL_MECH))
        ax.text(-1.6, _mm(L_total) / 2, f"總長 {_mm(env['L']):.2f} mm",
                rotation=90, ha="center", va="center", fontsize=9,
                color=COL_MECH)
        ax.text(x_lead, _mm((coil.y_lo + coil.y_hi) / 2) + 2.0,
                f" 線圈斷面由 A_winding = {params.A_winding * 1e6:.1f} mm²\n"
                f" 與磁路窗口高 {_mm(coil.height):.2f} mm\n"
                f" 反解寬 {_mm(coil.width):.2f} mm",
                fontsize=6.5, color=COL_COIL, ha="left", va="top")
        ax.set_xlim(-3, _mm(R_out) + 13)
        ax.set_ylim(-1.5, _mm(L_total) + 2)
        ax.set_aspect("equal")
        ax.set_xlabel("半徑 (mm)")
        ax.set_ylabel("軸向 (mm)")
        ax.set_title(
            f"電磁閥實體剖面（半剖，按比例，閥關 x=0）\n"
            f"外徑 {_mm(env['OD']):.1f} mm × 長 {_mm(env['L']):.1f} mm",
            fontsize=11)
        SECTION_OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(SECTION_OUT, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return SECTION_OUT


def flow_state(params, x, gas=None, cond=None):
    """流體狀態 [SI]：有效流通面積、質量流率與瓶頸所在。

    effective_area 是「簾幕」模型 min(pi*D*x, pi*D^2/4)：氣體沿孔上行，
    在閥芯端面與座面之間的環狀縫隙轉為徑向流出，故節流面是這道簾幕而
    非孔本身——直到行程大到簾幕面積追上孔面積為止。臨界點 x = D/4，
    而本案 x_stroke 恰等於 D/4（cases.py sizing chain 第 2 步：再多行程
    也買不到流量），這是圖上值得畫出來的機制。
    """
    if gas is None:
        gas = fluid_lib.N2
    if cond is None:
        from solenoid_model.cases import N2_25BAR
        cond = N2_25BAR.cond
    a_curtain = math.pi * params.D_seat_bore * max(x, 0.0)
    a_bore = math.pi * params.D_seat_bore ** 2 / 4.0
    A_eff = fluid_lib.effective_area(x, params)
    r = cond.P_down / cond.P_up if cond.P_up > 0 else 0.0
    return {
        "A_eff": A_eff, "a_curtain": a_curtain, "a_bore": a_bore,
        "mdot": fluid_lib.mdot_gas(x, params, gas, cond),
        "choked": r < fluid_lib.critical_pressure_ratio(gas),
        "limiter": "簾幕" if a_curtain < a_bore else "孔口",
        "cond": cond, "open": x > 0.0,
    }


def _draw_flow(ax, params, x, G, S, scale=1.0):
    """在半剖圖上畫流體路徑：孔內上行 → 簾幕徑向轉折 → 徑向流出。

    閥關時不畫流線（沒有流動），改以壓力箭頭表示 25 bar 壓在已坐封的
    閥芯上——那是坐封狀態下流體唯一的作用，也是失效關閉的一部分。
    """
    F = flow_state(params, x)
    r_bore = _mm(params.D_seat_bore / 2)
    r_pop = _mm(layout.LAYOUT_EMPIRICAL["d_poppet"]["value"] / 2)
    y_seat = _mm(S["y_seat_top"])
    y_in = _mm(0.0)

    if not F["open"]:
        # 坐封：壓力箭頭在孔內頂到閥芯底面，無流線。孔半徑僅 0.30 mm，
        # 故箭頭置於孔中線、標籤一律拉到孔外。
        ax.annotate("", xy=(r_bore * 0.5, y_seat - 0.10),
                    xytext=(r_bore * 0.5, y_in + 0.35),
                    arrowprops=dict(arrowstyle="-|>", color=COL_SEAL,
                                    lw=1.4 * scale, shrinkA=0, shrinkB=0))
        ax.plot([r_bore * 0.5, r_bore + 0.9], [y_seat * 0.45] * 2,
                color=COL_SEAL, lw=0.7)
        ax.text(r_bore + 1.0, y_seat * 0.45,
                f"N₂ {F['cond'].P_up / 1e5:.0f} bar 壓在已坐封的閥芯上\n"
                f"{params.delta_P * params.A_seat:.2f} N，與彈簧同向助封",
                fontsize=6.5, color=COL_SEAL, ha="left", va="center")
        return F

    # 開啟：孔內上行（只畫在孔內，到閥芯底面為止）
    ax.annotate("", xy=(r_bore * 0.5, _mm(S["y_pop_lo"]) - 0.05),
                xytext=(r_bore * 0.5, y_in + 0.35),
                arrowprops=dict(arrowstyle="-|>", color=COL_SEAL,
                                lw=1.6 * scale, shrinkA=0, shrinkB=0))
    # 簾幕轉折 + 徑向流出。以折線明畫 90 度轉折，不用 annotate 的
    # "angle" connectionstyle——垂直分量為零時它會除以零（本模組先前
    # 已在流道孔標註處踩過同一個雷）。
    y_curtain = _mm(S["y_seat_top"]) + _mm(S["seat_gap"]) / 2
    R_out_mm = _mm(G["R_shell"] + G["wall"])
    ax.annotate("", xy=(r_pop + 1.3, y_curtain),
                xytext=(r_bore * 0.5, y_curtain),
                arrowprops=dict(arrowstyle="-|>", color=COL_SEAL,
                                lw=1.6 * scale, shrinkA=0, shrinkB=0))
    # 標籤拉到外殼右側淨空，並以折線引出：貼在閥芯右邊會壓到彈簧標籤
    ax.plot([r_pop + 1.3, R_out_mm + 0.5, R_out_mm + 0.5],
            [y_curtain, y_curtain, y_curtain - 1.5],
            color=COL_SEAL, lw=0.7)
    ax.text(R_out_mm + 0.65, y_curtain - 1.6,
            f"簾幕 π·D·x = {F['A_eff'] * 1e6:.3f} mm²\n"
            f"{'壅塞流' if F['choked'] else '次音速'} "
            f"{F['mdot'] * 1e3:.2f} g/s\n徑向流出",
            fontsize=6.5, color=COL_SEAL, ha="left", va="top")
    return F


def _draw_state(ax, params, x, G, label, sub):
    """在 ax 上畫單一作動狀態的半剖圖（真實比例）。"""
    S = state_geometry(params, x=x)
    span, coil = G["span"], G["coil"]
    wall, R_shell, r_core = G["wall"], G["R_shell"], G["r_core"]
    R_out = R_shell + wall
    L = G["L_total"]
    r_stem = layout.LAYOUT_EMPIRICAL["d_stem"]["value"] / 2.0
    r_pop = layout.LAYOUT_EMPIRICAL["d_poppet"]["value"] / 2.0
    r_spr = G["r_spring_out"]

    # --- 不動件
    for name in ("端板（座端）", "端板（極端）"):
        lo, hi = span[name]
        ax.add_patch(Rectangle((0, _mm(lo)), _mm(R_out), _mm(hi - lo),
                               fc=COL_MAG, ec="k", lw=0.8))
    ax.add_patch(Rectangle((_mm(R_shell), 0), _mm(wall), _mm(L),
                           fc=COL_MAG, ec="k", lw=0.8))
    lo, hi = span["閥座座體"]
    ax.add_patch(Rectangle((0, _mm(lo)), _mm(R_out), _mm(hi - lo),
                           fc=COL_NONMAG, ec="k", lw=0.8))
    r_bore = _mm(params.D_seat_bore / 2)
    ax.add_patch(Rectangle((0, 0), r_bore, _mm(hi), fc="w", ec="none",
                           zorder=2))
    ax.plot([0, r_bore, r_bore], [0, 0, _mm(hi)], color=COL_SEAL, lw=1.2,
            zorder=3)
    # 固定彈簧座（不動）
    ax.add_patch(Rectangle((_mm(r_stem), _mm(S["y_seat_spr_lo"])),
                           _mm(R_out - r_stem), _mm(S["t_spring_seat"]),
                           fc=COL_NONMAG, ec="k", lw=0.8))
    # 固定極（不動）
    lo, hi = span["固定極"]
    ax.add_patch(Rectangle((0, _mm(lo)), _mm(r_core), _mm(hi - lo),
                           fc=COL_MAG, ec="k", lw=0.8))
    ax.text(_mm(r_core) / 2, _mm((lo + hi) / 2), "固定極", ha="center",
            va="center", fontsize=7, color="w")
    # 線圈與隔離套（不動）
    ax.add_patch(Rectangle((_mm(coil.r_in), _mm(coil.y_lo)),
                           _mm(coil.width), _mm(coil.height),
                           fc=COL_COIL, ec="k", lw=0.8))
    ax.text(_mm((coil.r_in + coil.r_out) / 2),
            _mm((coil.y_lo + coil.y_hi) / 2), "線圈", ha="center",
            va="center", fontsize=7)
    ax.add_patch(Rectangle((_mm(coil.r_in - G["t_sleeve"]), _mm(coil.y_lo)),
                           _mm(G["t_sleeve"]), _mm(coil.height),
                           fc="none", ec=COL_SEAL, lw=1.0, hatch="//"))
    # --- 動件組：閥芯 + 導杆 + 銜鐵（整體隨 x 平移）
    ax.add_patch(Rectangle((0, _mm(S["y_pop_lo"])), _mm(r_pop),
                           _mm(S["y_pop_hi"] - S["y_pop_lo"]),
                           fc=COL_MAG, ec="k", lw=1.4))
    ax.add_patch(Rectangle((0, _mm(S["y_pop_hi"])), _mm(r_stem),
                           _mm(S["y_arm_lo"] - S["y_pop_hi"]),
                           fc=COL_MAG, ec="k", lw=1.4))
    ax.add_patch(Rectangle((0, _mm(S["y_arm_lo"])), _mm(r_core),
                           _mm(S["y_arm_hi"] - S["y_arm_lo"]),
                           fc=COL_MAG, ec="k", lw=1.4))
    ax.text(_mm(r_core) / 2, _mm((S["y_arm_lo"] + S["y_arm_hi"]) / 2),
            "銜鐵", ha="center", va="center", fontsize=7, color="w")
    # --- 彈簧：下端隨閥芯上移，上端固定 → 長度 = 安裝長 - x
    ax.add_patch(Rectangle((_mm(r_stem), _mm(S["y_spring_lo"])),
                           _mm(r_spr - r_stem), _mm(S["L_spring"]),
                           fc="none", ec=COL_MECH, lw=1.6, ls="--"))
    ax.text(_mm(r_spr) + 0.35, _mm((S["y_spring_lo"] + S["y_spring_hi"]) / 2),
            f"彈簧 {_mm(S['L_spring']):.2f}\n{S['F_spring']:.2f} N ↓",
            fontsize=6.5, color=COL_MECH, ha="left", va="center")
    # --- 流體路徑（閥關畫壓力，閥開畫流線）
    _draw_flow(ax, params, x, G, S)
    # --- 座面基準線與運動箭頭
    ax.plot([0, _mm(R_out)], [_mm(S["y_seat_top"])] * 2, color=COL_SEAL,
            lw=0.7, ls=":", zorder=5)
    if S["energised"]:
        x_arrow = _mm(R_out) + 1.4
        y_mid = _mm((S["y_pop_lo"] + S["y_arm_hi"]) / 2)
        ax.annotate("", xy=(x_arrow, y_mid + 1.6),
                    xytext=(x_arrow, y_mid - 1.6),
                    arrowprops=dict(arrowstyle="-|>", color=COL_MECH, lw=2.2))
        ax.text(x_arrow + 0.35, y_mid,
                f"動件組\n上移 {_mm(params.x_stroke):.2f}", fontsize=7,
                color=COL_MECH, ha="left", va="center")
    ax.set_xlim(-1.0, _mm(R_out) + 4.2)
    ax.set_ylim(-1.0, _mm(L) + 1.0)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"{label}\n{sub}", fontsize=9.5)
    return S


def _draw_gap_detail(ax, params, x, G, color):
    """兩個間隙的放大圖：真實比例下 0.05-0.20 mm 看不見，故另開放大視窗。

    銜鐵與導杆的中段不按比例摺疊（已標示），但兩個間隙本身按 1:1 且
    兩狀態共用刻度，所以圖上的長度比等於真實比例。
    """
    S = state_geometry(params, x=x)
    r_core = G["r_core"]
    r_show = _mm(r_core)
    r_bore = _mm(params.D_seat_bore / 2)
    r_pop = _mm(layout.LAYOUT_EMPIRICAL["d_poppet"]["value"] / 2)
    r_stem = _mm(layout.LAYOUT_EMPIRICAL["d_stem"]["value"] / 2)

    # 摺疊後的繪圖座標：間隙按真實 mm，動件中段壓成固定高度
    body = 0.85
    stub = 0.20
    y_seat = 0.0
    y_pop_lo = y_seat + _mm(S["seat_gap"])
    y_pop_hi = y_pop_lo + stub
    y_arm_lo = y_pop_hi + body
    y_arm_hi = y_arm_lo + stub
    y_pole_lo = y_arm_hi + _mm(S["gap"])
    y_lo, y_hi = y_seat - 0.42, y_pole_lo + 0.38

    # 閥座
    ax.add_patch(Rectangle((0, y_lo), r_show * 1.35, y_seat - y_lo,
                           fc=COL_NONMAG, ec="k", lw=0.8))
    ax.add_patch(Rectangle((0, y_lo), r_bore, y_seat - y_lo, fc="w",
                           ec="none", zorder=2))
    ax.plot([0, r_bore, r_bore], [y_lo, y_lo, y_seat], color=COL_SEAL,
            lw=1.2, zorder=3)
    ax.text(r_show * 0.85, y_lo + (y_seat - y_lo) / 2, "閥座", fontsize=6.5,
            ha="center", va="center")
    # 閥芯（下端面按比例）
    ax.add_patch(Rectangle((0, y_pop_lo), r_pop, stub, fc=COL_MAG, ec="k",
                           lw=1.4))
    ax.text(r_pop / 2, y_pop_lo + stub / 2, "閥芯", fontsize=6,
            ha="center", va="center", color="w")
    # 動件中段（摺疊，含導杆與銜鐵本體）
    ax.add_patch(Rectangle((0, y_pop_hi), r_stem, body, fc=COL_MAG,
                           ec="none", alpha=0.35))
    for yb in (y_pop_hi, y_arm_lo):
        ax.plot([0, r_show], [yb] * 2, color="k", lw=0.7, ls=(0, (4, 2.5)),
                zorder=4)
    ax.text(r_show * 0.52, y_pop_hi + body / 2, "導杆＋銜鐵\n（中段不按比例）",
            fontsize=5.8, ha="center", va="center",
            bbox=dict(fc="w", ec="none", pad=0.6))
    # 銜鐵極面（按比例）
    ax.add_patch(Rectangle((0, y_arm_lo), r_show, stub, fc=COL_MAG, ec="k",
                           lw=1.4))
    # 固定極
    ax.add_patch(Rectangle((0, y_pole_lo), r_show, y_hi - y_pole_lo,
                           fc=COL_MAG, ec="k", lw=0.8))
    ax.text(r_show * 0.5, y_pole_lo + (y_hi - y_pole_lo) / 2, "固定極",
            fontsize=6.5, color="w", ha="center", va="center")

    x_dim = r_show * 1.14
    ax.annotate("", xy=(x_dim, y_arm_hi), xytext=(x_dim, y_pole_lo),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.3))
    ax.text(x_dim + 0.05, (y_arm_hi + y_pole_lo) / 2,
            f" 工作氣隙 {_mm(S['gap']):.2f}", fontsize=8.5, color=color,
            ha="left", va="center")
    if S["seat_gap"] > 0:
        ax.annotate("", xy=(x_dim, y_seat), xytext=(x_dim, y_pop_lo),
                    arrowprops=dict(arrowstyle="<->", color=COL_MECH, lw=1.3))
        ax.text(x_dim + 0.05, (y_seat + y_pop_lo) / 2,
                f" 開度 {_mm(S['seat_gap']):.2f}", fontsize=8.5,
                color=COL_MECH, ha="left", va="center")
        # 孔內上行 → 簾幕徑向轉出：放大圖才看得清這個 90 度轉折
        F = flow_state(params, x)
        ax.annotate("", xy=(r_bore * 0.45, y_pop_lo - 0.015),
                    xytext=(r_bore * 0.45, y_lo + 0.05),
                    arrowprops=dict(arrowstyle="-|>", color=COL_SEAL, lw=1.6))
        y_curtain = y_seat + (y_pop_lo - y_seat) / 2
        ax.plot([r_bore * 0.45, r_bore * 0.45], [y_seat, y_curtain],
                color=COL_SEAL, lw=1.6, solid_capstyle="butt")
        ax.annotate("", xy=(r_show * 1.05, y_curtain),
                    xytext=(r_bore * 0.45, y_curtain),
                    arrowprops=dict(arrowstyle="-|>", color=COL_SEAL,
                                    lw=1.6, shrinkA=0, shrinkB=0))
        ax.text(r_bore + 0.06, y_lo + 0.08, "入口 ↑", fontsize=6,
                color=COL_SEAL, ha="left", va="bottom")
    else:
        ax.plot([0, r_pop], [y_seat] * 2, color=COL_SEAL, lw=3.0,
                solid_capstyle="butt", zorder=5)
        ax.text(x_dim + 0.05, y_seat + stub * 0.8,
                " 坐封：開度 0（無流動）", fontsize=8, color=COL_SEAL,
                ha="left", va="center")
    ax.set_xlim(-0.05, r_show * 2.6)
    ax.set_ylim(y_lo - 0.05, y_hi + 0.05)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def generate_actuation(params=N2_25BAR_PH_PARAMS):
    """作動前後對照圖：同一動件組的兩個位置，真實比例 + 間隙放大。

    上排真實比例看封裝與力路，下排放大看兩個間隙——0.15 mm 行程在
    18.7 mm 的零件上只有 0.8% 高度，不放大則兩個狀態看起來完全一樣。
    """
    G = section_geometry(params)
    F_open = flow_state(params, params.x_stroke)
    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.6),
                                 gridspec_kw={"height_ratios": [2.4, 1.0]})
        _draw_state(axes[0][0], params, 0.0, G, "作動前（未通電）",
                    "閥關：彈簧預載壓閥芯坐封")
        _draw_state(axes[0][1], params, params.x_stroke, G,
                    "作動後（通電）", "閥開：磁吸力克服彈簧＋壓差")
        _draw_gap_detail(axes[1][0], params, 0.0, G, COL_MAG)
        _draw_gap_detail(axes[1][1], params, params.x_stroke, G, COL_MAG)
        for a in axes[1]:
            a.set_title("間隙放大（同一刻度）", fontsize=8)

        fig.suptitle("電磁閥作動前後對照（半剖，單位 mm）", fontsize=12,
                     y=0.975)
        fig.text(0.5, 0.048,
                 f"動件組＝閥芯＋導杆＋銜鐵，整體平移：開度 + 工作氣隙 "
                 f"恆等於 g0 = {_mm(params.g0):.2f} mm（dynamics 的 "
                 f"gap = g0 − x）\n"
                 f"彈簧套在導杆上，上端頂固定彈簧座、下端頂閥芯頂面："
                 f"動件上移 x 即壓縮 x，故 F = F_preload + k·x = "
                 f"{params.F_preload:.2f} → "
                 f"{params.F_preload + params.k_spring * params.x_stroke:.2f} N，"
                 f"方向恆向下（坐封方向）\n"
                 f"失效關閉：斷電後彈簧與壓差 "
                 f"{params.delta_P * params.A_seat:.2f} N 同向，"
                 f"將動件組推回坐封\n"
                 f"流體：氣體沿孔上行，在閥芯與座面之間的簾幕轉為徑向流出——"
                 f"節流面是簾幕 π·D·x，不是孔本身\n"
                 f"x_stroke = D/4 = {_mm(params.x_stroke):.2f} mm 時簾幕恰"
                 f"追上孔面積 {F_open['a_bore'] * 1e6:.3f} mm²，"
                 f"再多行程也買不到流量（壅塞 "
                 f"{F_open['mdot'] * 1e3:.2f} g/s）",
                 ha="center", fontsize=8.5)
        fig.subplots_adjust(top=0.90, bottom=0.135, hspace=0.20)
        ACTUATION_OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(ACTUATION_OUT, dpi=150)
        plt.close(fig)
    return ACTUATION_OUT


def generate_architecture(params=N2_25BAR_PH_PARAMS, V_bus=28.0,
                          cond=None, mdot_req=1.28e-3):
    """四個功能域方塊與跨域耦合點。

    方塊上的數字由 params 推導，不寫死字面值。流量同時標出「需求」與
    「模型實算」：1.28 g/s 是使用者給定的需求（cases.py 的孔徑就是由它
    反推，案例名稱即 "N₂ 25bar 1.28g/s"），1.302 g/s 是 fluid.mdot_gas
    在該案例宣告條件下的實現值（101.7% of target）。只標其中一個會讓
    讀者分不清哪個是規格、哪個是結果。

    cond 預設取 N2_25BAR 案例的宣告條件（P_up=25 bar, P_down=0 真空,
    T0=298.15 K）——不可自行假設 26→1 bar / 293 K，那會給 1.37 g/s 並
    與 cases.py 的 1.302 g/s 對不上。V_bus 為參數，因匯流排電壓屬於
    circuit 而非 ValveParams。
    """
    if cond is None:
        from solenoid_model.cases import N2_25BAR
        cond = N2_25BAR.cond
    mdot = fluid_lib.mdot_gas(params.x_stroke, params, fluid_lib.N2, cond)
    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, ax = plt.subplots(figsize=(10, 7))

        def box(x, y, w, h, label, color):
            ax.add_patch(Rectangle((x, y), w, h, fc=color, ec="k",
                                   alpha=0.25, lw=1.5))
            # CJK 後援字型（Arial Unicode MS 等）沒有 bold 字重，
            # weight="bold" 會靜默退回 regular（產生 findfont 警告卻無視覺
            # 效果），故改用放大字級來做標題強調，而非依賴 bold
            ax.text(x + w / 2, y + h - 0.3, label, ha="center",
                    va="top", fontsize=11.5)

        def node(x, y, text, w=2.0):
            ax.add_patch(Rectangle((x, y), w, 0.7, fc="w", ec="k"))
            ax.text(x + w / 2, y + 0.35, text, ha="center", va="center",
                    fontsize=8)

        def arrow(a, b, color="k", style="-|>"):
            ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style,
                                         mutation_scale=14, color=color,
                                         lw=1.6))

        box(0.3, 7.2, 9.4, 1.6, "電氣域", COL_COIL)
        node(0.8, 7.5, f"{V_bus:.0f} V 母線")
        node(3.4, 7.5, "峰值-保持驅動")
        node(6.6, 7.5, f"線圈 {params.N_turns:.0f} 匝")
        arrow((2.8, 7.85), (3.4, 7.85)); arrow((5.4, 7.85), (6.6, 7.85))

        box(0.3, 4.6, 9.4, 2.2, "磁域", COL_MAG)
        node(0.8, 5.6, "MMF = N·i"); node(3.4, 5.6, "鐵芯 B-H\n（飽和）")
        node(6.6, 5.6, f"工作氣隙 g0 = {_mm(params.g0):.2f} mm")
        node(3.4, 4.75, "磁軛回路\nA≥A_gap")
        arrow((1.8, 5.6), (1.8, 5.05)); arrow((2.8, 5.95), (3.4, 5.95))
        arrow((5.4, 5.95), (6.6, 5.95))
        arrow((7.6, 5.6), (7.6, 5.1)); arrow((7.6, 5.1), (5.4, 5.1))

        # 機械域／流體域：加高方塊（1.9→2.9）讓標題與第一排節點之間留出
        # 清楚間距（標題錨點 y=4.35，第一排節點上緣 y=3.85，間隙 0.5），
        # 且節點改用較窄寬度（1.85）並在兩節點間留 0.3 間隙，避免節點
        # 方塊互相貼邊；第二排節點（閥座止擋／壅塞流）與第一排之間也留
        # 0.25 間隙，不與其上緣重疊
        box(0.3, 1.7, 4.6, 2.9, "機械域", COL_MECH)
        node(0.55, 3.15, f"銜鐵 {params.m_arm * 1e3:.1f} g", w=1.85)
        node(2.7, 3.15, f"彈簧 {params.k_spring:.0f} N/m", w=1.85)
        node(1.6, 2.2, f"閥芯坐封 ⌀{_mm(params.D_seat_bore):.2f}", w=1.85)

        box(5.1, 1.7, 4.6, 2.9, "流體域", COL_SEAL)
        node(5.35, 3.15, f"N₂ {params.delta_P / 1e5:.0f} bar", w=1.85)
        node(7.5, 3.15, f"孔徑 {_mm(params.D_seat_bore):.2f} mm", w=1.85)
        # 需求 vs 實算：兩者都標，讀者才分得清規格與結果
        node(6.15, 2.2,
             f"壅塞流 {mdot * 1e3:.2f} g/s\n"
             f"（需求 {mdot_req * 1e3:.2f}，"
             f"{mdot / mdot_req * 100:.1f}%）", w=2.3)

        # 跨域耦合點
        arrow((7.6, 5.6), (7.6, 4.3), color=COL_MAG)
        ax.text(7.75, 4.9, "耦合①\n氣隙：磁↔機械", fontsize=8, color=COL_MAG)
        # 耦合②：改在兩個方塊的交界正下方（低於兩者的下緣 y=1.7）水平
        # 連接，不再貫穿機械域／流體域方塊內部；箭頭、標籤、下方註腳三者
        # 之間各留 >=0.3 的垂直間距，避免彼此貼在一起
        arrow((2.6, 1.4), (5.4, 1.4), color=COL_SEAL, style="<->")
        ax.text(4.0, 1.52, "耦合②　閥座：機械↔流體", ha="center",
                fontsize=8, color=COL_SEAL)

        ax.text(5.0, 0.05,
                "耦合①：氣隙同時決定磁阻與機械位置（dynamics.coupled_rhs）\n"
                "耦合②：閥座開度決定流量，噴流反作用力回饋進力平衡\n"
                f"流量：需求 {mdot_req * 1e3:.2f} g/s（孔徑由此反推），"
                f"fluid.mdot_gas 在全開 x={_mm(params.x_stroke):.2f} mm、"
                f"{cond.P_up / 1e5:.0f}→{cond.P_down / 1e5:.0f} bar、"
                f"{cond.T0:.1f} K 下實現 {mdot * 1e3:.3f} g/s",
                ha="center", fontsize=8.5)
        ax.set_xlim(0, 10); ax.set_ylim(-1.1, 9.2)
        ax.axis("off")
        ax.set_title("電磁閥架構圖：四個物理域與跨域耦合點", fontsize=12)
        ARCH_OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(ARCH_OUT, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return ARCH_OUT


if __name__ == "__main__":
    print(generate_section())
    print(generate_actuation())
    print(generate_architecture())
