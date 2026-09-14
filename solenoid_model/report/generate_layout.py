"""按比例剖面圖、作動前後對照圖與架構圖（P8）。

剖面圖用真實 mm 比例（set_aspect("equal")），與 docs/spec/valve_schematic.png
的符號示意圖不同——後者是參數對照用、不按比例。

作動對照圖的上排為真實比例，下排放大圖只把銜鐵「中段」摺疊（不按比例）
以突顯 0.05-0.20 mm 的間隙；間隙本身仍按 1:1 且兩狀態共用刻度，所以
圖上的長度比等於真實比例。
"""
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

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
    """
    parts, y = [], 0.0
    for name, t in (("端板（座端）", t_plate),
                    ("閥座座體", t_seat),
                    ("行程", params.x_stroke),
                    ("銜鐵", layout.armature_thickness(params)),
                    ("工作氣隙", params.g0),
                    ("固定極", t_pole),
                    ("彈簧腔", spring.L_free),
                    ("端板（彈簧端）", t_plate)):
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
    # 線圈軸向窗口：閥座頂面到彈簧腔底面，兩端都是實體零件的界面，
    # 線圈不得侵入（以鐵芯跨距置中會壓進座體 0.787 mm）。
    coil = coil_box(params, span["閥座座體"][1], span["彈簧腔"][0],
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

      x = 0          閥關：銜鐵坐在閥座上，氣隙最大 = g0
      x = x_stroke   閥開：銜鐵被吸向固定極，氣隙最小 = g0 - x_stroke

    行程間隙與工作氣隙是同一個剛體的兩端，此消彼長，和恆為
    x_stroke + (g0 - x_stroke) = g0。static 剖面圖把兩者同時畫成最大，
    是物理上不可能的狀態（銜鐵不能同時離座又離極最遠）。
    """
    if not 0.0 <= x <= params.x_stroke + 1e-15:
        raise ValueError(
            f"x={x * 1e3:.3f} mm 超出行程 0..{params.x_stroke * 1e3:.3f} mm")
    G = section_geometry(params)
    y_seat_top = G["span"]["閥座座體"][1]
    t_arm = layout.armature_thickness(params)
    # 銜鐵下緣 = 座面 + x（x=0 時貼座）；固定極下緣固定不動
    y_arm_lo = y_seat_top + x
    y_arm_hi = y_arm_lo + t_arm
    y_pole_lo = y_seat_top + params.x_stroke + t_arm + (params.g0
                                                        - params.x_stroke)
    gap = params.g0 - x
    return {
        "base": G, "x": x, "gap": gap,
        "seat_gap": x,
        "y_seat_top": y_seat_top,
        "y_arm_lo": y_arm_lo, "y_arm_hi": y_arm_hi,
        "y_pole_lo": y_pole_lo,
        "energised": x > 0.0,
    }


def generate_section(params=N2_25BAR_PH_PARAMS):
    """按比例軸對稱半剖圖，單位 mm。"""
    G = section_geometry(params)
    env, spring, span, coil = G["env"], G["spring"], G["span"], G["coil"]
    wall, R_shell, r_core = G["wall"], G["R_shell"], G["r_core"]
    L_total = G["L_total"]

    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, ax = plt.subplots(figsize=(8, 10))
        # 淨空標註用的引線終點：避開外殼標籤（位於 R_shell+wall+0.6），
        # 拉遠到 R_shell+wall+2.4 之外的空白區
        x_lead = _mm(R_shell + wall) + 2.4
        R_out = R_shell + wall

        def yl(name):
            return span[name][0]

        def yh(name):
            return span[name][1]

        # --- 端板 x2：axial_stack 明寫「端板 x2」，座端與彈簧端各一片。
        for name in ("端板（座端）", "端板（彈簧端）"):
            ax.add_patch(Rectangle((0, _mm(yl(name))), _mm(R_out),
                                   _mm(yh(name) - yl(name)),
                                   fc=COL_MAG, ec="k"))
        # --- 閥座座體（先畫實體，流道孔再疊上「挖空」）
        ax.add_patch(Rectangle((0, _mm(yl("閥座座體"))), _mm(R_out),
                               _mm(yh("閥座座體") - yl("閥座座體")),
                               fc=COL_NONMAG, ec="k"))
        ax.text(_mm(R_shell) * 0.7, _mm(sum(span["閥座座體"]) / 2),
                "閥座座體", ha="center", va="center", fontsize=8)
        # 流道孔：白色填色挖空座體，僅畫孔口與孔壁，讀作貫穿孔而非色塊。
        # 孔自座端端板底面一路貫穿到閥座頂面（流體入口）。
        r_bore = _mm(params.D_seat_bore / 2)
        y_bore_lo, y_bore_hi = 0.0, yh("閥座座體")
        ax.add_patch(Rectangle((0, _mm(y_bore_lo)), r_bore,
                               _mm(y_bore_hi - y_bore_lo),
                               fc="w", ec="none", zorder=2))
        ax.plot([0, r_bore, r_bore],
                [_mm(y_bore_lo), _mm(y_bore_lo), _mm(y_bore_hi)],
                color=COL_SEAL, lw=1.5, zorder=3)
        y_bore_label = _mm((y_bore_lo + y_bore_hi) / 2)
        ax.plot([r_bore, x_lead], [y_bore_label, y_bore_label],
                color=COL_SEAL, lw=0.9)
        ax.text(x_lead, y_bore_label,
                f" 流道孔徑 D_seat_bore = {_mm(params.D_seat_bore):.2f} mm",
                fontsize=7.5, color=COL_SEAL, ha="left", va="center")
        # --- 銜鐵、固定極
        for name, label in (("銜鐵", "銜鐵"), ("固定極", "固定極")):
            ax.add_patch(Rectangle((0, _mm(yl(name))), _mm(r_core),
                                   _mm(yh(name) - yl(name)),
                                   fc=COL_MAG, ec="k"))
            ax.text(_mm(r_core) / 2, _mm(sum(span[name]) / 2), label,
                    ha="center", va="center", fontsize=8, color="w")
        # --- 行程與工作氣隙：次毫米間隙，引線繞過線圈拉到右側淨空帶。
        # 垂直段走在鐵芯外緣與線圈內緣之間的隔離套帶，不穿過線圈色塊。
        y_coil_top, y_coil_bot = _mm(coil.y_hi), _mm(coil.y_lo)
        r_stub = _mm(r_core) + 0.12
        # 兩條引線都折到線圈上方的淨空並上下錯開，不再往線圈下方走——
        # 那裡是座體與其標籤，舊版的行程引線會橫穿「閥座座體」文字。
        for name, color, size, fmt in (
                ("行程", COL_MECH, 8, "行程 x_stroke = {:.2f} mm"),
                ("工作氣隙", COL_MAG, 9, "工作氣隙 g0 = {:.2f} mm")):
            lo, hi = span[name]
            y_mid = _mm((lo + hi) / 2)
            y_end = (y_coil_top + 0.6 if name == "工作氣隙"
                     else y_coil_top + 2.2)
            ax.plot([_mm(r_core), r_stub, r_stub, x_lead],
                    [y_mid, y_mid, y_end, y_end],
                    color=color, lw=1.0, solid_capstyle="butt")
            ax.text(x_lead, y_end, " " + fmt.format(_mm(hi - lo)),
                    fontsize=size, color=color, ha="left", va="center")
        # --- 線圈：斷面由 params.A_winding / l_turn_mean 推導（見 coil_box），
        # 內緣貼隔離套外側，不留不明空隙。高度大於鐵芯跨距時對稱溢出，
        # 屬真實幾何，另加註說明而非隱藏。
        ax.add_patch(Rectangle((_mm(coil.r_in), _mm(coil.y_lo)),
                               _mm(coil.width), _mm(coil.height),
                               fc=COL_COIL, ec="k"))
        ax.text(_mm((coil.r_in + coil.r_out) / 2),
                _mm((coil.y_lo + coil.y_hi) / 2),
                f"線圈\n{params.N_turns:.0f} 匝\n"
                f"{coil.area * 1e6:.1f} mm²",
                ha="center", va="center", fontsize=8)
        # --- 隔離套：沿線圈內緣，貼在鐵芯外側
        ax.add_patch(Rectangle((_mm(coil.r_in - G["t_sleeve"]),
                                _mm(coil.y_lo)), _mm(G["t_sleeve"]),
                               _mm(coil.height),
                               fc="none", ec=COL_SEAL, lw=1.5, hatch="//"))
        # --- 彈簧：真實外徑 D_coil + d_wire（舊版誤畫到 r_core，寬 2.35 倍），
        # 長度用安裝長（閥關時已被預載壓縮），非自由長。
        r_spr = G["r_spring_out"]
        L_inst = G["spring_installed"]
        y_spr = yl("彈簧腔")
        ax.add_patch(Rectangle((0, _mm(y_spr)), _mm(r_spr), _mm(L_inst),
                               fc="none", ec=COL_MECH, lw=1.5, ls="--"))
        # 標籤置於彈簧正上方的腔內淨空（安裝長之上到腔頂），純垂直讓位，
        # 不往右伸：往右會與線圈左上角及右側的行程引線打結。
        ax.text(_mm(r_spr) + 0.3, _mm(yh("彈簧腔")) - 0.2,
                f"彈簧 d{_mm(spring.d_wire):.2f} / D{_mm(spring.D_coil):.2f}\n"
                f"安裝長 {_mm(L_inst):.2f} mm"
                f"（自由長 {_mm(spring.L_free):.2f}）",
                ha="left", va="top", fontsize=6.5, color=COL_MECH)
        # --- 外殼
        ax.add_patch(Rectangle((_mm(R_shell), 0), _mm(wall), _mm(L_total),
                               fc=COL_MAG, ec="k"))
        ax.text(_mm(R_out) + 0.6, _mm(L_total) / 2,
                f"外殼 {_mm(wall):.2f} mm", fontsize=8, color=COL_MAG,
                rotation=90, va="center")
        # --- 中心線
        ax.axvline(0, color="k", lw=0.8, ls="-.")
        ax.text(0.1, _mm(L_total) + 0.6, "軸心", fontsize=7)
        # --- 總尺寸標註：畫到 L_total，且 L_total 恆等於 envelope()['L']
        ax.add_patch(FancyArrowPatch((-1.2, 0), (-1.2, _mm(L_total)),
                                     arrowstyle="<->", mutation_scale=12,
                                     color=COL_MECH))
        ax.text(-1.6, _mm(L_total) / 2, f"總長 {_mm(env['L']):.2f} mm",
                rotation=90, ha="center", va="center", fontsize=9,
                color=COL_MECH)
        # 線圈斷面來源註記：讓讀者知道這個框不是畫好看的，是 A_winding 反解。
        # 放在線圈中段右側的大片淨空，避開下方 D_seat_bore 標籤所在的 y。
        ax.text(x_lead, _mm((coil.y_lo + coil.y_hi) / 2) - 2.0,
                f" 線圈斷面由 A_winding = {params.A_winding * 1e6:.1f} mm²\n"
                f" 與軸向窗口高 {_mm(coil.height):.2f} mm\n"
                f" 反解寬 {_mm(coil.width):.2f} mm",
                fontsize=6.5, color=COL_COIL, ha="left", va="top")
        ax.set_xlim(-3, _mm(R_out) + 11)
        ax.set_ylim(-1.5, _mm(L_total) + 2)
        ax.set_aspect("equal")
        ax.set_xlabel("半徑 (mm)")
        ax.set_ylabel("軸向 (mm)")
        ax.set_title(
            f"電磁閥實體剖面（半剖，按比例）\n"
            f"外徑 {_mm(env['OD']):.1f} mm × 長 {_mm(env['L']):.1f} mm",
            fontsize=11)
        SECTION_OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(SECTION_OUT, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return SECTION_OUT


def _draw_state(ax, params, x, G, label, sub):
    """在 ax 上畫單一作動狀態的半剖圖（真實比例）。"""
    S = state_geometry(params, x=x)
    span, coil = G["span"], G["coil"]
    wall, R_shell, r_core = G["wall"], G["R_shell"], G["r_core"]
    R_out = R_shell + wall
    L = G["L_total"]
    spring = G["spring"]
    # 彈簧安裝長隨行程變化：閥開時再被多壓 x（腔體不變，彈簧更短）
    L_spr = G["spring_installed"] - x

    # 不動件：端板 x2、外殼、閥座座體、固定極
    for name in ("端板（座端）", "端板（彈簧端）"):
        lo, hi = span[name]
        ax.add_patch(Rectangle((0, _mm(lo)), _mm(R_out), _mm(hi - lo),
                               fc=COL_MAG, ec="k", lw=0.8))
    ax.add_patch(Rectangle((_mm(R_shell), 0), _mm(wall), _mm(L),
                           fc=COL_MAG, ec="k", lw=0.8))
    lo, hi = span["閥座座體"]
    ax.add_patch(Rectangle((0, _mm(lo)), _mm(R_out), _mm(hi - lo),
                           fc=COL_NONMAG, ec="k", lw=0.8))
    # 流道孔（挖空）：閥關時被銜鐵封住，閥開時連通
    r_bore = _mm(params.D_seat_bore / 2)
    ax.add_patch(Rectangle((0, 0), r_bore, _mm(hi),
                           fc="w", ec="none", zorder=2))
    ax.plot([0, r_bore, r_bore], [0, 0, _mm(hi)],
            color=COL_SEAL, lw=1.2, zorder=3)
    # 固定極（不動）
    ax.add_patch(Rectangle((0, _mm(S["y_pole_lo"])), _mm(r_core),
                           _mm(span["固定極"][1] - S["y_pole_lo"]),
                           fc=COL_MAG, ec="k", lw=0.8))
    ax.text(_mm(r_core) / 2, _mm((S["y_pole_lo"] + span["固定極"][1]) / 2),
            "固定極", ha="center", va="center", fontsize=7, color="w")
    # 線圈與隔離套（不動）
    ax.add_patch(Rectangle((_mm(coil.r_in), _mm(coil.y_lo)),
                           _mm(coil.width), _mm(coil.height),
                           fc=COL_COIL, ec="k", lw=0.8))
    ax.text(_mm((coil.r_in + coil.r_out) / 2),
            _mm((coil.y_lo + coil.y_hi) / 2), "線圈",
            ha="center", va="center", fontsize=7)
    ax.add_patch(Rectangle((_mm(coil.r_in - G["t_sleeve"]), _mm(coil.y_lo)),
                           _mm(G["t_sleeve"]), _mm(coil.height),
                           fc="none", ec=COL_SEAL, lw=1.0, hatch="//"))
    # 動件：銜鐵（依 x 位移）
    ax.add_patch(Rectangle((0, _mm(S["y_arm_lo"])), _mm(r_core),
                           _mm(S["y_arm_hi"] - S["y_arm_lo"]),
                           fc=COL_MAG, ec="k", lw=1.4))
    ax.text(_mm(r_core) / 2, _mm((S["y_arm_lo"] + S["y_arm_hi"]) / 2),
            "銜鐵", ha="center", va="center", fontsize=7, color="w")
    # 彈簧（安裝長隨 x 縮短：閥開時被多壓 x）
    r_spr = G["r_spring_out"]
    y_spr = span["彈簧腔"][0]
    ax.add_patch(Rectangle((0, _mm(y_spr)), _mm(r_spr), _mm(L_spr),
                           fc="none", ec=COL_MECH, lw=1.6, ls="--"))
    ax.text(_mm(r_spr) + 0.3, _mm(y_spr + L_spr / 2),
            f"彈簧\n{_mm(L_spr):.2f}", fontsize=6.5, color=COL_MECH,
            ha="left", va="center")
    # 座面位置與銜鐵運動方向：兩圖唯一的差別就在這裡
    ax.plot([0, _mm(R_out)], [_mm(S["y_seat_top"])] * 2,
            color=COL_SEAL, lw=0.7, ls=":", zorder=5)
    if S["energised"]:
        ax.annotate("", xy=(_mm(r_core) * 0.55, _mm(S["y_arm_lo"]) - 0.15),
                    xytext=(_mm(r_core) * 0.55,
                            _mm(S["y_arm_lo"]) - 1.5),
                    arrowprops=dict(arrowstyle="-|>", color=COL_MECH,
                                    lw=2.0))
        ax.text(_mm(r_core) * 0.55 + 0.25, _mm(S["y_arm_lo"]) - 0.9,
                f"上移 {_mm(params.x_stroke):.2f}", fontsize=7,
                color=COL_MECH, ha="left", va="center")
    ax.set_xlim(-1.0, _mm(R_out) + 3.4)
    ax.set_ylim(-1.0, _mm(L) + 1.0)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"{label}\n{sub}", fontsize=9.5)
    return S


def _draw_gap_detail(ax, params, x, G, color):
    """氣隙區域放大圖：真實比例下 0.15 mm 看不見，故另開放大視窗。

    視窗只取座面與固定極下緣之間，銜鐵本體僅露出上下各一小截（以剖斷
    線示意其延續），否則 5.08 mm 厚的銜鐵會把 0.05-0.20 mm 的間隙壓成
    看不見的細線——放大圖的重點是間隙，不是銜鐵。
    """
    S = state_geometry(params, x=x)
    r_core = G["r_core"]
    r_show = _mm(r_core)
    r_bore = _mm(params.D_seat_bore / 2)

    # 座標為「摺疊後」的繪圖座標（mm）：銜鐵中段 5.08 mm 不按比例，壓成
    # 固定的 0.9 繪圖單位並以剖斷線標示。間隙（0.05-0.20 mm）與銜鐵端面
    # 則按 1:1 真實比例，兩狀態共用同一刻度，所以 0.20 與 0.05 的長度比
    # 在圖上仍是 4:1，可直接目視比較。
    stub = 0.20             # 銜鐵端面露出高度（繪圖單位 = mm 真實）
    body = 0.9              # 中段摺疊後的固定高度（不按比例）
    y_seat = 0.0
    y_arm_lo = y_seat + _mm(S["seat_gap"])
    y_arm_hi_drawn = y_arm_lo + stub + body + stub
    y_pole_lo = y_arm_hi_drawn + _mm(S["gap"])
    y_lo, y_hi = y_seat - 0.45, y_pole_lo + 0.40

    # 閥座座體
    ax.add_patch(Rectangle((0, y_lo), r_show * 1.30, y_seat - y_lo,
                           fc=COL_NONMAG, ec="k", lw=0.8))
    ax.add_patch(Rectangle((0, y_lo), r_bore, y_seat - y_lo,
                           fc="w", ec="none", zorder=2))
    ax.plot([0, r_bore, r_bore], [y_lo, y_lo, y_seat],
            color=COL_SEAL, lw=1.2, zorder=3)
    ax.text(r_show * 0.72, y_lo + (y_seat - y_lo) / 2, "閥座",
            fontsize=6.5, ha="center", va="center")
    # 銜鐵：上下端面按真實比例，中段摺疊
    ax.add_patch(Rectangle((0, y_arm_lo), r_show, stub,
                           fc=COL_MAG, ec="k", lw=1.4))
    ax.add_patch(Rectangle((0, y_arm_hi_drawn - stub), r_show, stub,
                           fc=COL_MAG, ec="k", lw=1.4))
    ax.add_patch(Rectangle((0, y_arm_lo + stub), r_show, body,
                           fc=COL_MAG, ec="none", alpha=0.35))
    for yb in (y_arm_lo + stub, y_arm_hi_drawn - stub):
        ax.plot([0, r_show], [yb] * 2, color="k", lw=0.7,
                ls=(0, (4, 2.5)), zorder=4)
    ax.text(r_show * 0.5, y_arm_lo + stub + body / 2,
            f"銜鐵 {_mm(layout.armature_thickness(params)):.2f}\n（中段不按比例）",
            fontsize=6, ha="center", va="center")
    # 固定極
    ax.add_patch(Rectangle((0, y_pole_lo), r_show, y_hi - y_pole_lo,
                           fc=COL_MAG, ec="k", lw=0.8))
    ax.text(r_show * 0.5, y_pole_lo + (y_hi - y_pole_lo) / 2, "固定極",
            fontsize=6.5, color="w", ha="center", va="center")

    x_dim = r_show * 1.10
    # 工作氣隙（真實比例）
    ax.annotate("", xy=(x_dim, y_arm_hi_drawn), xytext=(x_dim, y_pole_lo),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.3))
    ax.text(x_dim + 0.05, (y_arm_hi_drawn + y_pole_lo) / 2,
            f" 工作氣隙 {_mm(S['gap']):.2f}", fontsize=8.5, color=color,
            ha="left", va="center")
    # 座開度
    if S["seat_gap"] > 0:
        ax.annotate("", xy=(x_dim, y_seat), xytext=(x_dim, y_arm_lo),
                    arrowprops=dict(arrowstyle="<->", color=COL_MECH, lw=1.3))
        ax.text(x_dim + 0.05, (y_seat + y_arm_lo) / 2,
                f" 開度 {_mm(S['seat_gap']):.2f}", fontsize=8.5,
                color=COL_MECH, ha="left", va="center")
        ax.annotate("", xy=(r_bore * 0.45, y_arm_lo - 0.02),
                    xytext=(r_bore * 0.45, y_lo + 0.05),
                    arrowprops=dict(arrowstyle="-|>", color=COL_SEAL, lw=1.6))
        ax.text(r_bore + 0.06, (y_lo + y_seat) / 2 + 0.08, "流動",
                fontsize=7, color=COL_SEAL, ha="left", va="center")
    else:
        ax.plot([0, r_show], [y_seat] * 2, color=COL_SEAL, lw=3.0,
                solid_capstyle="butt", zorder=5)
        # 標籤上移半個 stub，避免與座體內的「閥座」字樣同高相撞
        ax.text(x_dim + 0.05, y_seat + stub * 0.7,
                " 坐封：開度 0（無流動）", fontsize=8, color=COL_SEAL,
                ha="left", va="center")
    ax.set_xlim(-0.05, r_show * 2.45)
    ax.set_ylim(y_lo - 0.05, y_hi + 0.05)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def generate_actuation(params=N2_25BAR_PH_PARAMS):
    """作動前後對照圖：同一剛體的兩個位置，真實比例 + 氣隙放大。

    上排真實比例看封裝，下排放大看間隙——0.15 mm 行程在 18.7 mm 的
    零件上只有 0.8% 高度，不放大則兩個狀態看起來完全一樣。
    """
    G = section_geometry(params)
    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, axes = plt.subplots(2, 2, figsize=(9.5, 9.0),
                                 gridspec_kw={"height_ratios": [2.4, 1.0]})
        _draw_state(axes[0][0], params, 0.0, G, "作動前（未通電）",
                    "閥關：彈簧預載壓銜鐵坐封")
        _draw_state(axes[0][1], params, params.x_stroke, G,
                    "作動後（通電）", "閥開：磁吸力克服彈簧+壓差")
        _draw_gap_detail(axes[1][0], params, 0.0, G, COL_MAG)
        _draw_gap_detail(axes[1][1], params, params.x_stroke, G, COL_MAG)
        axes[1][0].set_title("氣隙區放大", fontsize=8)
        axes[1][1].set_title("氣隙區放大", fontsize=8)

        fig.suptitle("電磁閥作動前後對照（半剖，單位 mm）", fontsize=12,
                     y=0.975)
        fig.text(0.5, 0.055,
                 f"銜鐵是單一剛體：行程開度 + 工作氣隙 恆等於 g0 = "
                 f"{_mm(params.g0):.2f} mm（dynamics 的 gap = g0 − x）\n"
                 f"作動前 x=0：開度 0 ／ 氣隙 {_mm(params.g0):.2f} mm　→　"
                 f"作動後 x={_mm(params.x_stroke):.2f} mm：開度 "
                 f"{_mm(params.x_stroke):.2f} mm ／ 氣隙 "
                 f"{_mm(params.g0 - params.x_stroke):.2f} mm\n"
                 f"失效關閉：斷電後彈簧預載 {params.F_preload:.2f} N 與壓差 "
                 f"{params.delta_P * params.A_seat:.2f} N 同向，將銜鐵推回坐封",
                 ha="center", fontsize=8.5)
        fig.subplots_adjust(top=0.90, bottom=0.135, hspace=0.18)
        ACTUATION_OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(ACTUATION_OUT, dpi=150)
        plt.close(fig)
    return ACTUATION_OUT


def generate_architecture():
    """四個功能域方塊與跨域耦合點。"""
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
        node(0.8, 7.5, "28 V 母線"); node(3.4, 7.5, "峰值-保持驅動")
        node(6.6, 7.5, "線圈 4000 匝")
        arrow((2.8, 7.85), (3.4, 7.85)); arrow((5.4, 7.85), (6.6, 7.85))

        box(0.3, 4.6, 9.4, 2.2, "磁域", COL_MAG)
        node(0.8, 5.6, "MMF = N·i"); node(3.4, 5.6, "鐵芯 B-H\n（飽和）")
        node(6.6, 5.6, "工作氣隙 g0")
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
        node(0.55, 3.15, "銜鐵 0.8 g", w=1.85)
        node(2.7, 3.15, "彈簧 4000 N/m", w=1.85)
        node(1.6, 2.2, "閥座止擋", w=1.85)

        box(5.1, 1.7, 4.6, 2.9, "流體域", COL_SEAL)
        node(5.35, 3.15, "N₂ 25 bar", w=1.85)
        node(7.5, 3.15, "孔徑 0.60 mm", w=1.85)
        node(6.4, 2.2, "1.28 g/s 壅塞流", w=1.85)

        # 跨域耦合點
        arrow((7.6, 5.6), (7.6, 4.3), color=COL_MAG)
        ax.text(7.75, 4.9, "耦合①\n氣隙：磁↔機械", fontsize=8, color=COL_MAG)
        # 耦合②：改在兩個方塊的交界正下方（低於兩者的下緣 y=1.7）水平
        # 連接，不再貫穿機械域／流體域方塊內部；箭頭、標籤、下方註腳三者
        # 之間各留 >=0.3 的垂直間距，避免彼此貼在一起
        arrow((2.6, 1.4), (5.4, 1.4), color=COL_SEAL, style="<->")
        ax.text(4.0, 1.1, "耦合②　閥座：機械↔流體", ha="center",
                fontsize=8, color=COL_SEAL)

        ax.text(5.0, 0.35,
                "耦合①：氣隙同時決定磁阻與機械位置（dynamics.coupled_rhs）\n"
                "耦合②：閥座開度決定流量，噴流反作用力回饋進力平衡",
                ha="center", fontsize=8.5)
        ax.set_xlim(0, 10); ax.set_ylim(-0.5, 9.2)
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
