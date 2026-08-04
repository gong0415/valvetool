"""Generate the annotated valve cross-section schematic for docs/spec/PARAMS.md."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

OUT_PATH = Path(__file__).resolve().parents[2] / "docs" / "spec" / "valve_schematic.png"

COL_MAG = "#4477aa"    # magnetic quantities
COL_MECH = "#228833"   # mechanical quantities
COL_FLUID = "#cc3311"  # fluid quantities

# macOS matplotlib default font lacks Traditional Chinese glyphs (renders as
# tofu boxes + "Glyph missing from current font" warnings). Force a CJK-capable
# fallback chain. Scoped via rc_context inside generate() (not set globally at
# import time) so importing this module never changes matplotlib's rcParams
# for unrelated callers (e.g. generate_p0_report.py) sharing the process.
_CJK_FONT_RC = {
    "font.sans-serif": ["Arial Unicode MS", "PingFang TC", "Heiti TC"],
    "axes.unicode_minus": False,
}


def _dim_arrow(ax, xy_a, xy_b, text, color, text_dx=0.0, text_dy=0.0):
    ax.add_patch(FancyArrowPatch(xy_a, xy_b, arrowstyle="<->", mutation_scale=12, color=color, lw=1.4))
    mid = ((xy_a[0] + xy_b[0]) / 2 + text_dx, (xy_a[1] + xy_b[1]) / 2 + text_dy)
    ax.text(*mid, text, color=color, fontsize=9, ha="left", va="center")


def generate():
    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, ax = plt.subplots(figsize=(7.5, 9))

        # Coil windows (left/right of core)
        ax.add_patch(Rectangle((1.0, 6.0), 1.4, 3.2, fc="#f4a261", ec="k", label="coil"))
        ax.add_patch(Rectangle((5.6, 6.0), 1.4, 3.2, fc="#f4a261", ec="k"))
        ax.text(1.7, 7.6, "線圈\nN_turns\nR_coil_20C", ha="center", va="center", fontsize=8.5)
        ax.text(6.3, 7.6, "線圈\nN_turns\nR_coil_20C", ha="center", va="center", fontsize=8.5)

        # Core yoke (top bridge + side legs)
        ax.add_patch(Rectangle((1.0, 9.2), 6.0, 0.6, fc="#b0b0b0", ec="k"))
        ax.add_patch(Rectangle((0.4, 6.0), 0.6, 3.8, fc="#b0b0b0", ec="k"))
        ax.add_patch(Rectangle((7.0, 6.0), 0.6, 3.8, fc="#b0b0b0", ec="k"))
        ax.text(4.0, 9.5, "鐵芯磁軛（頂部軛鐵）", ha="center", va="center", fontsize=9, color=COL_MAG)
        ax.text(4.0, 8.6, "鐵芯磁軛\nl_core, mu_r_core, B_sat", ha="center", fontsize=9, color=COL_MAG)

        # Fixed pole face
        ax.add_patch(Rectangle((2.4, 6.0), 3.2, 0.5, fc="#8d99ae", ec="k"))
        ax.text(4.0, 6.25, "極面（面積 A_gap）", ha="center", fontsize=9, color=COL_MAG)

        # Armature below pole, separated by air gap
        ax.add_patch(Rectangle((2.4, 4.6), 3.2, 0.9, fc="#6c757d", ec="k"))
        ax.text(4.0, 5.05, "銜鐵 m_arm", ha="center", fontsize=10, color="w")

        # Air gap dimension: g0 (rest) with x measured downward-to-up
        _dim_arrow(ax, (6.0, 5.5), (6.0, 6.0), " gap = g0 − x", COL_MAG, text_dx=0.1)
        ax.text(6.1, 5.2, "x：位移 0→x_stroke\n（吸合方向：向上）", fontsize=8, color=COL_MECH)

        # Spring (zigzag) beside armature
        zx = [2.0] * 7
        zy = [4.6, 4.35, 4.1, 3.85, 3.6, 3.35, 3.1]
        for i in range(len(zy) - 1):
            ax.plot([zx[i] - 0.15 + 0.3 * (i % 2), zx[i] - 0.15 + 0.3 * ((i + 1) % 2)],
                    [zy[i], zy[i + 1]], color=COL_MECH, lw=1.6)
        ax.text(0.6, 3.8, "彈簧\nk_spring\nF_preload", fontsize=9, color=COL_MECH)

        # Poppet stem + poppet
        ax.add_patch(Rectangle((3.8, 3.0), 0.4, 1.6, fc="#6c757d", ec="k"))
        ax.add_patch(Rectangle((3.3, 2.5), 1.4, 0.5, fc="#6c757d", ec="k"))
        ax.text(4.0, 2.72, "poppet", ha="center", fontsize=9, color="w")

        # Seat with bore
        ax.add_patch(Rectangle((2.2, 1.6), 1.3, 0.9, fc="#adb5bd", ec="k"))
        ax.add_patch(Rectangle((4.5, 1.6), 1.3, 0.9, fc="#adb5bd", ec="k"))
        ax.text(4.0, 1.2, "閥座", ha="center", fontsize=9)

        # Sealing land annotation (A_seat) and bore (D_seat_bore)
        _dim_arrow(ax, (3.5, 2.5), (4.5, 2.5), "", COL_FLUID)
        ax.text(4.62, 2.62, "密封環受壓面積 A_seat（力平衡用）", fontsize=8, color=COL_FLUID)
        _dim_arrow(ax, (3.5, 1.75), (4.5, 1.75), "", COL_FLUID)
        ax.text(4.62, 1.75, "流道孔徑 D_seat_bore（流量用）\nA_eff(x)=min(π·D·x, π·D²/4)", fontsize=8, color=COL_FLUID)

        # Pressure difference arrows
        ax.annotate("", xy=(4.0, 0.55), xytext=(4.0, 1.35),
                    arrowprops=dict(arrowstyle="-|>", color=COL_FLUID, lw=2))
        ax.text(4.2, 0.8, "ΔP = delta_P（P_up 上游 → P_down 下游）", fontsize=9, color=COL_FLUID)

        ax.set_xlim(0, 8.4)
        ax.set_ylim(0, 10)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title("電磁閥剖面示意（幾何量對照 ValveParams 欄位）", fontsize=12)

        OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return OUT_PATH


WINDING_OUT = OUT_PATH.parent / "winding_window_schematic.png"
SEAL_OUT = OUT_PATH.parent / "seal_land_schematic.png"


def generate_winding_window():
    """繞線窗口剖面：A_winding / k_fill / l_turn_mean 定義圖（P4）。"""
    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, ax = plt.subplots(figsize=(7, 5))
        # 鐵芯中柱 + 窗口外框
        ax.add_patch(Rectangle((0.5, 1.0), 1.2, 6.0, fc="#b0b0b0", ec="k"))
        ax.text(1.1, 4.0, "鐵芯\n中柱", ha="center", fontsize=9)
        ax.add_patch(Rectangle((2.0, 1.5), 3.6, 5.0, fc="none", ec=COL_MAG, lw=2))
        ax.text(3.8, 6.9, "繞線窗口面積 A_winding", ha="center",
                fontsize=10, color=COL_MAG)
        # 線材圓形示意（部分填充 → k_fill）
        for row in range(5):
            for col in range(6):
                ax.add_patch(plt.Circle((2.4 + col * 0.6, 2.0 + row * 0.9),
                                        0.26, fc="#f4a261", ec="k", lw=0.6))
        ax.text(3.8, 0.9, "銅填充率 k_fill = 銅截面 / A_winding"
                          "（圓形線材 + 絕緣層 → 必小於 1）",
                ha="center", fontsize=9)
        # 平均匝長（俯視小圖）
        circ = plt.Circle((7.6, 4.0), 1.3, fc="none", ec=COL_MECH, lw=2,
                          linestyle="--")
        ax.add_patch(circ)
        ax.text(7.6, 4.0, "俯視", ha="center", fontsize=8)
        ax.text(7.6, 2.2, "平均匝長 l_turn_mean\n（一匝繞行的平均周長）",
                ha="center", fontsize=9, color=COL_MECH)
        ax.set_xlim(0, 9.4)
        ax.set_ylim(0, 7.6)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title("繞線窗口示意（對照 ValveParams：A_winding, k_fill, l_turn_mean）",
                     fontsize=11)
        fig.savefig(WINDING_OUT, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return WINDING_OUT


def generate_seal_land():
    """密封 land 放大剖面：w_land / R_tip / D_seal / 殘留間隙與漏流路徑（P6）。"""
    with matplotlib.rc_context(_CJK_FONT_RC):
        fig, ax = plt.subplots(figsize=(7.5, 5))
        # 閥座（左右兩塊，中間是流道孔）
        ax.add_patch(Rectangle((0.5, 1.0), 3.0, 1.5, fc="#adb5bd", ec="k"))
        ax.add_patch(Rectangle((6.5, 1.0), 3.0, 1.5, fc="#adb5bd", ec="k"))
        ax.text(5.0, 1.3, "流道 D_seat_bore", ha="center", fontsize=9,
                color=COL_FLUID)
        # poppet 圓弧端（半徑 R_tip）
        arc = plt.Circle((5.0, 6.2), 3.4, fc="#6c757d", ec="k")
        ax.add_patch(arc)
        ax.set_clip_on(True)
        ax.text(5.0, 5.4, "poppet 端面\n曲率半徑 R_tip", ha="center",
                fontsize=10, color="w")
        # 密封 land（接觸帶）；label 移到左側全開放區（x<1.6，poppet 圓弧
        # 在此完全不下探，最安全），並拉高與下方漏流文字保持大間距
        _dim_arrow(ax, (2.5, 2.75), (3.5, 2.75), " 密封 land 寬度 w_land",
                   COL_FLUID, text_dx=-4.5, text_dy=2.0)
        # 殘留間隙 + 漏流路徑；箭頭留在原位（貼近 land 實際位置），文字
        # 移到左側、與上方 w_land 標籤保持約 1.7 個資料單位的垂直間距
        ax.annotate("", xy=(3.6, 2.6), xytext=(0.9, 2.6),
                    arrowprops=dict(arrowstyle="<-", color=COL_FLUID, lw=1.6,
                                    linestyle=":"))
        ax.text(-1.0, 3.0, "殘留間隙 h = Rq_c·(1−φ)\n漏流沿 land 寬度方向流過",
                fontsize=9, color=COL_FLUID)
        # D_seal 標註
        _dim_arrow(ax, (3.0, 0.55), (7.0, 0.55),
                   " 密封環等效直徑 D_seal = √(4·A_seat/π)", COL_MECH,
                   text_dx=-2.6, text_dy=-0.35)
        ax.set_xlim(0, 10)
        ax.set_ylim(-0.4, 6.4)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title("密封 land 放大示意（對照 ValveParams：w_land, R_tip, A_seat）",
                     fontsize=11)
        fig.savefig(SEAL_OUT, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return SEAL_OUT


if __name__ == "__main__":
    print(generate())
    print(generate_winding_window())
    print(generate_seal_land())
