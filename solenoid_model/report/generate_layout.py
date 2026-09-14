"""按比例剖面圖與架構圖（P8）。

剖面圖用真實 mm 比例（set_aspect("equal")），與 docs/spec/valve_schematic.png
的符號示意圖不同——後者是參數對照用、不按比例。
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
    print(generate_architecture())
