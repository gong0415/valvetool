"""按比例剖面圖與架構圖（P8）。

剖面圖用真實 mm 比例（set_aspect("equal")），與 docs/spec/valve_schematic.png
的符號示意圖不同——後者是參數對照用、不按比例。
"""
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

COIL_OD, COIL_ID, COIL_H = 18.0e-3, 6.0e-3, 12.0e-3


def _mm(x):
    return x * 1e3


def generate_section(params=N2_25BAR_PH_PARAMS):
    """按比例軸對稱半剖圖，單位 mm。"""
    env = layout.envelope(params, coil_OD=COIL_OD)
    wall = env["wall"].adopted
    r_core = layout.core_radius(params)
    t_arm = layout.armature_thickness(params)
    t_plate = layout.end_plate_thickness(params)
    t_pole = layout.LAYOUT_EMPIRICAL["t_fixed_pole"]["value"]
    t_seat = layout.LAYOUT_EMPIRICAL["t_seat_body"]["value"]
    t_sleeve = layout.LAYOUT_EMPIRICAL["t_sleeve"]["value"]
    R_shell = COIL_OD / 2

    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, ax = plt.subplots(figsize=(8, 10))
        # 淨空標註用的引線終點：避開外殼標籤（位於 R_shell+wall+0.6），
        # 拉遠到 R_shell+wall+2.4 之外的空白區
        x_lead = _mm(R_shell + wall) + 2.4
        y = 0.0
        # 閥座座體（先畫實體，流道孔再疊在上面「挖空」，避免孔看起來像脫離的碎片）
        ax.add_patch(Rectangle((0, _mm(y)), _mm(R_shell + wall),
                               _mm(t_seat), fc=COL_NONMAG, ec="k"))
        ax.text(_mm(R_shell) * 0.7, _mm(y + t_seat / 2), "閥座座體",
                ha="center", va="center", fontsize=8)
        # 流道孔：白色填色挖空座體，僅在孔口（頂邊）與孔壁（右邊）畫線，
        # 底邊與軸心邊不畫線，讀作貫穿座體的孔而非獨立色塊
        r_bore = _mm(params.D_seat_bore / 2)
        ax.add_patch(Rectangle((0, _mm(y)), r_bore, _mm(t_seat),
                               fc="w", ec="none", zorder=2))
        ax.plot([0, r_bore, r_bore], [_mm(y), _mm(y), _mm(y + t_seat)],
                color=COL_SEAL, lw=1.5, zorder=3)
        # 標籤文字固定在座體中段高度（y=t_seat/2），落在 x_lead 右側淨空
        # 帶內、x 軸刻度標籤上方，不落到座標軸下緣以外與刻度文字重疊。
        # 直接以純水平線連接（xy/xytext 同一 y），改用 ax.plot 而非
        # annotate 的 "angle" connectionstyle，避免其在零垂直分量下
        # 除以零產生 RuntimeWarning
        y_bore_label = _mm(y + t_seat / 2)
        ax.plot([r_bore, x_lead], [y_bore_label, y_bore_label],
                color=COL_SEAL, lw=0.9)
        ax.text(x_lead, y_bore_label,
                f" 流道孔徑 D_seat_bore = {_mm(params.D_seat_bore):.2f} mm",
                fontsize=7.5, color=COL_SEAL, ha="left", va="center")
        y += t_seat
        y_stroke_lo, y_stroke_hi = y, y + params.x_stroke
        y += params.x_stroke      # 行程（閥開時的間隙）
        # 銜鐵
        ax.add_patch(Rectangle((0, _mm(y)), _mm(r_core), _mm(t_arm),
                               fc=COL_MAG, ec="k"))
        ax.text(_mm(r_core) / 2, _mm(y + t_arm / 2), "銜鐵",
                ha="center", va="center", fontsize=8, color="w")
        y += t_arm
        y_gap_lo, y_gap_hi = y, y + params.g0
        y += params.g0            # 工作氣隙
        # 固定極
        ax.add_patch(Rectangle((0, _mm(y)), _mm(r_core), _mm(t_pole),
                               fc=COL_MAG, ec="k"))
        ax.text(_mm(r_core) / 2, _mm(y + t_pole / 2), "固定極",
                ha="center", va="center", fontsize=8, color="w")
        # 行程與工作氣隙：次毫米級間隙，用引線標到右側淨空帶。兩個間隙的
        # y 範圍都落在線圈徑向涵蓋之內（線圈 r=3.0-9.0mm），若引線直接
        # 水平拉到 x_lead 會整條穿過線圈色塊與「線圈 4000 匝」文字，看起
        # 來像是線圈的一部分。改成在間隙本體的淨空半徑帶（r_core 與線圈
        # 內緣 COIL_ID/2 之間，此帶內僅有薄薄的隔離套，仍留有可見空隙）
        # 先水平拉一小段，再垂直折到線圈外緣正上／正下方的淨空 y，最後
        # 水平拉到標籤：三段都不與線圈本體重疊。
        # 垂直段刻意離開 r_core 一點（不緊貼銜鐵/固定極的黑色描邊），
        # 避免細線被同一 x 上的元件邊框線蓋掉而不可見。
        # （y_coil 提前算出：線圈實際繪製於本函式稍後，但幾何在此已確定）
        y_coil = y - COIL_H / 2
        y_coil_top = _mm(y_coil + COIL_H)
        y_coil_bot = _mm(y_coil)
        y_stroke_mid = _mm((y_stroke_lo + y_stroke_hi) / 2)
        y_gap_mid = _mm((y_gap_lo + y_gap_hi) / 2)
        y_stroke_lead = y_coil_bot - 0.6
        y_gap_lead = y_coil_top + 0.6
        r_stub = _mm(r_core) + 0.35  # 離開元件描邊，落在隔離套內緣之前的淨空帶
        ax.plot([_mm(r_core), r_stub, r_stub, x_lead],
                [y_stroke_mid, y_stroke_mid, y_stroke_lead, y_stroke_lead],
                color=COL_MECH, lw=1.0, solid_capstyle="butt")
        ax.text(x_lead, y_stroke_lead,
                f" 行程 x_stroke = {_mm(y_stroke_hi - y_stroke_lo):.2f} mm",
                fontsize=8, color=COL_MECH, ha="left", va="center")
        ax.plot([_mm(r_core), r_stub, r_stub, x_lead],
                [y_gap_mid, y_gap_mid, y_gap_lead, y_gap_lead],
                color=COL_MAG, lw=1.0, solid_capstyle="butt")
        ax.text(x_lead, y_gap_lead,
                f" 工作氣隙 g0 = {_mm(y_gap_hi - y_gap_lo):.2f} mm",
                fontsize=9, color=COL_MAG, ha="left", va="center")
        # 線圈（在固定極與銜鐵徑向外側；y_coil 已於上方引線區塊算出）
        ax.add_patch(Rectangle((_mm(COIL_ID / 2), _mm(y_coil)),
                               _mm((COIL_OD - COIL_ID) / 2), _mm(COIL_H),
                               fc=COL_COIL, ec="k"))
        ax.text(_mm((COIL_OD + COIL_ID) / 4), _mm(y_coil + COIL_H / 2),
                f"線圈\n{params.N_turns:.0f} 匝", ha="center", va="center",
                fontsize=8)
        # 隔離套
        ax.add_patch(Rectangle((_mm(COIL_ID / 2 - t_sleeve), _mm(y_coil)),
                               _mm(t_sleeve), _mm(COIL_H),
                               fc="none", ec=COL_SEAL, lw=1.5, hatch="//"))
        y += t_pole
        # 彈簧腔 + 端板
        spring = layout.spring_geometry(params, **layout.SPRING_DESIGN)
        ax.add_patch(Rectangle((0, _mm(y)), _mm(r_core), _mm(spring.L_free),
                               fc="none", ec=COL_MECH, lw=1.5, ls="--"))
        ax.text(_mm(r_core) / 2, _mm(y + spring.L_free / 2),
                f"彈簧\nd{_mm(spring.d_wire):.2f}", ha="center", va="center",
                fontsize=7, color=COL_MECH)
        y += spring.L_free
        ax.add_patch(Rectangle((0, _mm(y)), _mm(R_shell + wall),
                               _mm(t_plate), fc=COL_MAG, ec="k"))
        y += t_plate
        # 外殼
        ax.add_patch(Rectangle((_mm(R_shell), 0), _mm(wall), _mm(y),
                               fc=COL_MAG, ec="k"))
        ax.text(_mm(R_shell + wall) + 0.6, _mm(y) / 2,
                f"外殼 {_mm(wall):.2f} mm", fontsize=8, color=COL_MAG,
                rotation=90, va="center")
        # 中心線
        ax.axvline(0, color="k", lw=0.8, ls="-.")
        ax.text(0.1, _mm(y) + 0.6, "軸心", fontsize=7)
        # 總尺寸標註
        ax.add_patch(FancyArrowPatch((-1.2, 0), (-1.2, _mm(y)),
                                     arrowstyle="<->", mutation_scale=12,
                                     color=COL_MECH))
        ax.text(-1.6, _mm(y) / 2, f"總長 {_mm(env['L']):.2f} mm",
                rotation=90, ha="center", va="center", fontsize=9,
                color=COL_MECH)
        ax.set_xlim(-3, _mm(R_shell + wall) + 10)
        ax.set_ylim(-1.5, _mm(y) + 2)
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
