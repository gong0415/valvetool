# P8：B-H 飽和曲線 + 實體佈局推導層 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 導入單值 B-H 飽和曲線修正被高估的電磁吸力，並在其上建立從 `ValveParams` 推導實體尺寸鏈的 `layout.py`，產出按比例剖面圖／架構圖與 GUI 第五頁籤。

**Architecture:** B-H 曲線作為**第五個正交模式旋鈕**加進 `magnetics.py`：新增 `mag=None` 引數，`None` 時走原線性碼路徑、與現有結果 bit-level 相同（沿用 `medium`／`T_coil`／`seat` 的既有雙模式原則）。`layout.py` 是純函式推導層，不 import `dynamics`／`magnetics`／`fluid`（沿用 `sizing.py`、`winding.py` 的不變式），把磁通連續、環向應力、彈簧封閉解三組公式組成尺寸鏈。

**Tech Stack:** Python 3.12、pytest、matplotlib（Agg）、Streamlit、scipy（`brentq`，`thermal.py` 已用）

**Spec:** `docs/superpowers/specs/2026-09-04-physical-layout-design.md`

## Global Constraints

- **相容性鐵律**：`mag=None` 時 `magnetics` 所有函式與現行結果 bit-level 相同。`solenoid_model/tests/test_dynamics_characterization.py`（乾跑基準凍結）全程綠燈，**不得修改該檔**。
- **分層不變式**：`layout.py` 不得 import `dynamics`／`magnetics`／`fluid`。`magnetics.py` 維持不 import `dynamics`。
- **本體論（規格 §2.0）**：不得引入無法對應的隱藏自由度。Fröhlich-Kennelly 的 μ_i、B_sat 一律由現有 `ValveParams` 欄位推得。
- **EMPIRICAL 標註**：非第一性推導的數值一律標 `EMPIRICAL`，比照現有 `k_pack`／`C_d`／`damping_coeff` 的註解風格。
- **中文字型**：所有 matplotlib 產圖沿用 `generate_schematic.py` 的 `_CJK_FONT_RC`，以 `matplotlib.rc_context` 區域套用，不得全域改 rcParams。
- **執行環境**：測試一律 `.venv/bin/python -m pytest`。
- **基準案例**：`cases.N2_25BAR_PH`（`A_gap=20 mm²`、`g0=0.20 mm`、`x_stroke=0.15 mm`、`N=4000`、`B_sat=2.15 T`、`mu_r_core=4000`、`l_core=0.045 m`、`k_spring=4000 N/m`、`F_preload=2.2 N`、`m_arm=0.8 g`）。

## File Structure

| 檔案 | 責任 |
|---|---|
| `solenoid_model/magnetization.py`（新增） | B-H 曲線資料型別：`AnalyticBH`／`TabulatedBH`，只管 `B_of_H`／`H_of_B`。獨立於 `magnetics.py` 以免磁路求解與曲線來源糾纏。 |
| `solenoid_model/magnetics.py`（修改） | 新增 `solve_flux_density()`；`reluctance`／`flux_density`／`magnetic_force_closing`／`saturation_check` 加 `mag=None`。 |
| `solenoid_model/layout.py`（新增） | 實體尺寸鏈推導：磁路截面、外殼壁厚、彈簧幾何、軸向堆疊、外型包絡、一致性檢查。 |
| `solenoid_model/report/generate_layout.py`（新增） | 產生 `layout_section.png`（按比例剖面）與 `layout_architecture.png`（架構圖）。 |
| `solenoid_model/app.py`（修改） | 新增 `render_layout_tab()`，`st.tabs` 由四頁擴為五頁。 |
| `solenoid_model/tests/test_magnetization.py`（新增） | 曲線本身：極限、單調、往返一致、兩來源一致。 |
| `solenoid_model/tests/test_magnetics_bh.py`（新增） | 磁路求解與相容性鐵律。 |
| `solenoid_model/tests/test_layout.py`（新增） | 尺寸鏈推導與不變式。 |

---

### Task 1: B-H 曲線資料型別

**Files:**
- Create: `solenoid_model/magnetization.py`
- Test: `solenoid_model/tests/test_magnetization.py`

**Interfaces:**
- Consumes: 無（本任務為葉節點，只用標準庫）
- Produces:
  - `AnalyticBH(mu_r: float, B_sat: float)` — frozen dataclass，方法 `B_of_H(H: float) -> float`、`H_of_B(B: float) -> float`
  - `TabulatedBH(points: tuple)` — frozen dataclass，`points` 為 `((H, B), ...)` 嚴格遞增，同樣兩個方法
  - `from_params(params) -> AnalyticBH` — 由 `params.mu_r_core`／`params.B_sat` 建構
  - 模組常數 `MU_0 = 4*math.pi*1e-7`

- [ ] **Step 1: 寫失敗測試**

建立 `solenoid_model/tests/test_magnetization.py`：

```python
import math

import pytest

from solenoid_model import magnetization
from solenoid_model.cases import N2_25BAR_PH_PARAMS


def test_analytic_low_H_reduces_to_linear():
    """H->0 時 Frohlich-Kennelly 必須退化為 B = mu_r*mu_0*H。"""
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    H = 0.1
    linear = 4000.0 * magnetization.MU_0 * H
    assert math.isclose(bh.B_of_H(H), linear, rel_tol=1e-3)


def test_analytic_high_H_approaches_B_sat_without_exceeding():
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    assert bh.B_of_H(1e6) == pytest.approx(2.1491, abs=1e-3)
    for H in (1e3, 1e6, 1e9):
        assert bh.B_of_H(H) < 2.15


def test_analytic_is_monotonic():
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    prev = -1.0
    for k in range(0, 9):
        B = bh.B_of_H(10.0 ** k)
        assert B > prev
        prev = B


def test_analytic_round_trip_H_of_B():
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    for H in (1.0, 1e2, 1e4):
        assert math.isclose(bh.H_of_B(bh.B_of_H(H)), H, rel_tol=1e-9)


def test_H_of_B_at_or_above_B_sat_is_infinite():
    bh = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    assert math.isinf(bh.H_of_B(2.15))
    assert math.isinf(bh.H_of_B(3.0))


def test_from_params_uses_existing_fields_only():
    bh = magnetization.from_params(N2_25BAR_PH_PARAMS)
    assert bh.mu_r == N2_25BAR_PH_PARAMS.mu_r_core
    assert bh.B_sat == N2_25BAR_PH_PARAMS.B_sat


def test_tabulated_matches_analytic_on_its_own_points():
    """同一組點上，查表與解析式必須一致（介面可互換）。"""
    analytic = magnetization.AnalyticBH(mu_r=4000.0, B_sat=2.15)
    Hs = (1.0, 10.0, 100.0, 1000.0, 10000.0)
    table = magnetization.TabulatedBH(
        tuple((H, analytic.B_of_H(H)) for H in Hs))
    for H in Hs:
        assert math.isclose(table.B_of_H(H), analytic.B_of_H(H), rel_tol=1e-12)


def test_tabulated_interpolates_between_points():
    table = magnetization.TabulatedBH(((0.0, 0.0), (100.0, 1.0)))
    assert math.isclose(table.B_of_H(50.0), 0.5, rel_tol=1e-12)


def test_tabulated_clamps_above_last_point():
    """超出表格上界時鉗位在最後一點，不外插（外插會給出非物理的 B）。"""
    table = magnetization.TabulatedBH(((0.0, 0.0), (100.0, 1.0)))
    assert math.isclose(table.B_of_H(1e6), 1.0, rel_tol=1e-12)


def test_tabulated_rejects_non_ascending_points():
    with pytest.raises(ValueError, match="ascending"):
        magnetization.TabulatedBH(((0.0, 0.0), (100.0, 1.0), (50.0, 0.5)))


def test_tabulated_rejects_too_few_points():
    with pytest.raises(ValueError, match="at least 2"):
        magnetization.TabulatedBH(((0.0, 0.0),))


def test_analytic_rejects_nonpositive_B_sat():
    with pytest.raises(ValueError, match="B_sat"):
        magnetization.AnalyticBH(mu_r=4000.0, B_sat=0.0)
```

- [ ] **Step 2: 跑測試確認失敗**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_magnetization.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'solenoid_model.magnetization'`

- [ ] **Step 3: 實作**

建立 `solenoid_model/magnetization.py`：

```python
"""Single-valued B-H magnetization curves (P8).

Pure data + interpolation: these types answer only "what B does this H
give, and vice versa". Solving a magnetic circuit with one of them lives
in magnetics.py, so the curve source and the circuit solver stay
independent -- a tabulated measurement and the analytic law are
interchangeable at the same interface.

Hysteresis is deliberately absent. These curves are single-valued, so
they model saturation but not remanence: B depends on H alone, with no
magnetization history. For this valve that is the right trade -- the
remanent holding force in a gap-dominated circuit is ~0.02 N against a
2.2 N spring preload, two orders of magnitude short of a non-reseating
risk, while a hysteresis model would turn the ODE from memoryless into
history-dependent. See the design spec section 3A.5.
"""
import math
from dataclasses import dataclass

MU_0 = 4.0 * math.pi * 1e-7  # vacuum permeability, H/m


@dataclass(frozen=True)
class AnalyticBH:
    """Frohlich-Kennelly saturation law.

        B(H) = mu_i*H / (1 + mu_i*H/B_sat),   mu_i = mu_r * MU_0

    Introduces NO new degrees of freedom: both constants come from
    existing ValveParams fields (mu_r_core, B_sat), which is what keeps
    this compatible with the spec's ontology rule (section 2.0). It
    reduces to the linear model as H->0 and asymptotes to B_sat as
    H->infinity, so it is a strict refinement of the constant-mu_r model
    rather than a different one.
    """
    mu_r: float
    B_sat: float

    def __post_init__(self):
        if not (math.isfinite(self.mu_r) and self.mu_r > 0.0):
            raise ValueError(f"mu_r={self.mu_r} must be finite and positive")
        if not (math.isfinite(self.B_sat) and self.B_sat > 0.0):
            raise ValueError(f"B_sat={self.B_sat} must be finite and positive")

    @property
    def mu_i(self):
        """Initial permeability [H/m]."""
        return self.mu_r * MU_0

    def B_of_H(self, H):
        """Flux density [T] at field strength H [A/m]."""
        if H < 0.0:
            raise ValueError(f"H={H} must be non-negative")
        return self.mu_i * H / (1.0 + self.mu_i * H / self.B_sat)

    def H_of_B(self, B):
        """Field strength [A/m] for flux density B [T].

        Infinite at and above B_sat: no finite field reaches saturation
        under this law. Callers bisecting on B must keep the bracket
        strictly below B_sat.
        """
        if B < 0.0:
            raise ValueError(f"B={B} must be non-negative")
        if B >= self.B_sat:
            return math.inf
        return B / (self.mu_i * (1.0 - B / self.B_sat))


@dataclass(frozen=True)
class TabulatedBH:
    """Measured (H, B) points with piecewise-linear interpolation.

    Same interface as AnalyticBH, so a solver cannot tell them apart.
    Follows thermal.py's low-temperature table: strictly ascending
    points, validated at construction, and CLAMPED rather than
    extrapolated outside the table -- extrapolating a B-H curve past its
    last measured point produces unbounded, non-physical B.
    """
    points: tuple

    def __post_init__(self):
        if len(self.points) < 2:
            raise ValueError(
                f"points has {len(self.points)} entry; needs at least 2 "
                "(H, B) points to interpolate")
        for H, B in self.points:
            if not (math.isfinite(H) and math.isfinite(B)):
                raise ValueError(f"point ({H}, {B}) must be finite")
            if H < 0.0 or B < 0.0:
                raise ValueError(f"point ({H}, {B}) must be non-negative")
        Hs = [H for H, _ in self.points]
        Bs = [B for _, B in self.points]
        if any(h1 <= h0 for h0, h1 in zip(Hs, Hs[1:])):
            raise ValueError(f"H values {Hs} must be strictly ascending")
        if any(b1 <= b0 for b0, b1 in zip(Bs, Bs[1:])):
            raise ValueError(f"B values {Bs} must be strictly ascending")

    @property
    def B_sat(self):
        """Largest representable B [T] -- the table's last point."""
        return self.points[-1][1]

    def B_of_H(self, H):
        if H < 0.0:
            raise ValueError(f"H={H} must be non-negative")
        return _interp(H, self.points)

    def H_of_B(self, B):
        if B < 0.0:
            raise ValueError(f"B={B} must be non-negative")
        if B >= self.B_sat:
            return math.inf
        return _interp(B, tuple((b, h) for h, b in self.points))


def _interp(x, points):
    """Piecewise-linear interpolation, clamped at both ends."""
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return points[-1][1]


def from_params(params):
    """AnalyticBH built from a ValveParams' existing magnetic fields."""
    return AnalyticBH(mu_r=params.mu_r_core, B_sat=params.B_sat)
```

- [ ] **Step 4: 跑測試確認通過**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_magnetization.py -v`
Expected: 全數 PASS

- [ ] **Step 5: Commit**

```bash
git add solenoid_model/magnetization.py solenoid_model/tests/test_magnetization.py
git commit -m "feat(magnetics): add single-valued B-H curve types

Frohlich-Kennelly analytic law and a tabulated variant behind one
interface. Zero new degrees of freedom: mu_i and B_sat both come from
existing ValveParams fields.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: 非線性磁路求解 + 相容性鐵律

**Files:**
- Modify: `solenoid_model/magnetics.py`
- Test: `solenoid_model/tests/test_magnetics_bh.py`

**Interfaces:**
- Consumes: Task 1 的 `magnetization.AnalyticBH`／`from_params`／`MU_0`
- Produces:
  - `magnetics.solve_flux_density(gap, i, params, mag) -> float` — 對 B 二分求解 `H(B)*l_core + B*gap/MU_0 = N*i`
  - `magnetics.flux_density(gap, i, params, mag=None) -> float`
  - `magnetics.magnetic_force_closing(gap, i, params, mag=None) -> float`
  - `magnetics.saturation_check(gap, i, params, mag=None) -> bool`
  - `reluctance`／`inductance`／`dinductance_dgap` 簽章**不變**（線性磁路概念，非線性下無單一值）

- [ ] **Step 1: 寫失敗測試**

建立 `solenoid_model/tests/test_magnetics_bh.py`：

```python
import math

import pytest

from solenoid_model import magnetics, magnetization
from solenoid_model.baseline import BASELINE_PARAMS
from solenoid_model.cases import N2_25BAR_PH_PARAMS

PH = N2_25BAR_PH_PARAMS
I_PEAK = PH.V_bus / PH.R_coil_20C


def test_mag_none_is_bit_identical_to_linear():
    """相容性鐵律：mag=None 必須與現行線性結果逐位元相同。"""
    for params in (BASELINE_PARAMS, PH):
        for gap in (params.g0, params.g0 - params.x_stroke):
            for i in (0.01, 0.1, 0.35):
                linear_B = (params.N_turns * i
                            / magnetics.reluctance(gap, params) / params.A_gap)
                assert magnetics.flux_density(gap, i, params) == linear_B
                assert magnetics.flux_density(gap, i, params, None) == linear_B


def test_solver_satisfies_the_mmf_balance():
    """解回代 H(B)*l_core + B*gap/mu_0 必須等於 N*i。"""
    mag = magnetization.from_params(PH)
    gap = PH.g0
    B = magnetics.solve_flux_density(gap, I_PEAK, PH, mag)
    mmf = mag.H_of_B(B) * PH.l_core + B * gap / magnetics.MU_0
    assert math.isclose(mmf, PH.N_turns * I_PEAK, rel_tol=1e-6)


def test_saturating_B_stays_below_B_sat():
    """線性模型在閉合氣隙報 8.16 T；飽和模型必須低於 B_sat。"""
    mag = magnetization.from_params(PH)
    gap = PH.g0 - PH.x_stroke
    assert magnetics.flux_density(gap, I_PEAK, PH) > 8.0        # 線性：非物理
    assert magnetics.flux_density(gap, I_PEAK, PH, mag) < PH.B_sat


def test_saturating_force_never_exceeds_linear_force():
    """飽和只會減少吸力，不會增加。"""
    mag = magnetization.from_params(PH)
    for k in range(6):
        gap = PH.g0 - PH.x_stroke * k / 5
        F_lin = magnetics.magnetic_force_closing(gap, I_PEAK, PH)
        F_sat = magnetics.magnetic_force_closing(gap, I_PEAK, PH, mag)
        assert F_sat <= F_lin


def test_force_at_rest_gap_matches_verified_value():
    """對照 spec 3A.3 已驗證數值（線性 44.58 N -> 飽和 26.47 N）。"""
    mag = magnetization.from_params(PH)
    F = magnetics.magnetic_force_closing(PH.g0, I_PEAK, PH, mag)
    assert F == pytest.approx(26.47, abs=0.05)


def test_force_at_closed_gap_matches_verified_value():
    """spec 3A.3：閉合氣隙 530.30 N -> 32.70 N。"""
    mag = magnetization.from_params(PH)
    F = magnetics.magnetic_force_closing(PH.g0 - PH.x_stroke, I_PEAK, PH, mag)
    assert F == pytest.approx(32.70, abs=0.05)


def test_saturating_force_is_nearly_flat_across_stroke():
    """spec 3A.3 發現 1：真實吸力沿行程近乎持平（線性模型報 12 倍成長）。"""
    mag = magnetization.from_params(PH)
    F_rest = magnetics.magnetic_force_closing(PH.g0, I_PEAK, PH, mag)
    F_closed = magnetics.magnetic_force_closing(
        PH.g0 - PH.x_stroke, I_PEAK, PH, mag)
    assert F_closed / F_rest < 1.5          # 飽和：+24%
    F_lin_rest = magnetics.magnetic_force_closing(PH.g0, I_PEAK, PH)
    F_lin_closed = magnetics.magnetic_force_closing(
        PH.g0 - PH.x_stroke, I_PEAK, PH)
    assert F_lin_closed / F_lin_rest > 10.0  # 線性：12 倍，對照組


def test_low_current_converges_to_linear_model():
    """電流小到不飽和時，兩模型必須收斂。"""
    mag = magnetization.from_params(PH)
    i = 1e-4
    B_lin = magnetics.flux_density(PH.g0, i, PH)
    B_sat = magnetics.flux_density(PH.g0, i, PH, mag)
    assert math.isclose(B_lin, B_sat, rel_tol=1e-3)


def test_zero_current_gives_zero_flux():
    mag = magnetization.from_params(PH)
    assert magnetics.solve_flux_density(PH.g0, 0.0, PH, mag) == 0.0


def test_saturation_check_honours_mag():
    mag = magnetization.from_params(PH)
    gap = PH.g0 - PH.x_stroke
    assert magnetics.saturation_check(gap, I_PEAK, PH) is True
    assert magnetics.saturation_check(gap, I_PEAK, PH, mag) is False
```

- [ ] **Step 2: 跑測試確認失敗**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_magnetics_bh.py -v`
Expected: FAIL — `AttributeError: module 'solenoid_model.magnetics' has no attribute 'solve_flux_density'`

- [ ] **Step 3: 實作**

修改 `solenoid_model/magnetics.py`。在檔案頂端 `import math` 之後加入模組 docstring 與 import：

```python
"""Magnetic circuit: reluctance, inductance, coenergy force, saturation.

Two modes, following the dual-mode rule in ARCHITECTURE.md. `mag=None`
(the default) is the original linear circuit with a constant
mu_r_core -- bit-identical to every result this module produced before
P8, which the characterization tests freeze. Passing a
magnetization.AnalyticBH or TabulatedBH switches to a nonlinear solve
where the core saturates.

Why this matters: the linear model reports B = 8.16 T at the closed gap
of the N2_25BAR_PH case, which is not a physical flux density. It
overstates force where force is scarcest, and it makes the force look
like it explodes as the gap closes when the real curve is nearly flat.

reluctance/inductance/dinductance_dgap keep their original signatures
and stay linear-only on purpose: with a saturating core there is no
single reluctance -- it depends on the operating point -- so a `mag`
argument there would return a number that quietly means something else.
Nonlinear callers go through solve_flux_density instead.
"""
import math

from solenoid_model import magnetization

MU_0 = 4 * math.pi * 1e-7  # vacuum permeability, H/m

# Bisection bracket is [0, B_sat) and the residual is monotonic in B, so
# a fixed iteration count converges to machine precision without needing
# a tolerance argument: 200 halvings of a ~2 T bracket is far below
# double precision.
_BISECT_ITERS = 200
```

然後**取代** `flux_density`、`magnetic_force_closing`、`saturation_check` 三個函式，並新增 `solve_flux_density`：

```python
def solve_flux_density(gap, i, params, mag):
    """Gap flux density [T] with a saturating core.

    Solves the MMF balance around the loop,

        H(B)*l_core + B*gap/MU_0 = N*i

    for B by bisection. The left side is strictly increasing in B (H(B)
    is monotonic and the gap term is linear), so the root is unique and
    bisection cannot miss it. The bracket is [0, B_sat): H(B_sat) is
    infinite, so the upper end is approached, never reached.

    Sign convention: current magnitude only. Force goes as B**2, so the
    sign of i does not change the force, and a negative i would only
    break the bracket.
    """
    NI = abs(params.N_turns * i)
    if NI == 0.0:
        return 0.0
    lo, hi = 0.0, mag.B_sat * (1.0 - 1e-12)
    for _ in range(_BISECT_ITERS):
        mid = 0.5 * (lo + hi)
        if mag.H_of_B(mid) * params.l_core + mid * gap / MU_0 < NI:
            lo = mid
        else:
            hi = mid
    return lo


def flux_density(gap, i, params, mag=None):
    """Gap flux density [T]. mag=None -> linear circuit (unchanged)."""
    if mag is None:
        flux = params.N_turns * i / reluctance(gap, params)
        return flux / params.A_gap
    return solve_flux_density(gap, i, params, mag)


def magnetic_force_closing(gap, i, params, mag=None):
    """Force [N] pulling the armature toward the pole, always >= 0.

    Linear mode keeps the coenergy form 0.5*i^2*|dL/dgap| untouched.
    Saturating mode uses the Maxwell stress form B^2*A/(2*mu_0), which
    is the same physics: with a linear circuit the two agree exactly,
    but only the B-form stays valid once the core saturates, because
    dL/dgap is no longer a constant of the operating point.
    """
    if mag is None:
        return -0.5 * i ** 2 * dinductance_dgap(gap, params)
    B = solve_flux_density(gap, i, params, mag)
    return B ** 2 * params.A_gap / (2.0 * MU_0)


def saturation_check(gap, i, params, mag=None):
    """True when the core is driven past B_sat.

    With `mag` supplied this is always False by construction -- the
    solver cannot return B >= B_sat -- so it reports whether the LINEAR
    model would have exceeded B_sat, i.e. whether using `mag` changes
    the answer materially at this operating point.
    """
    return flux_density(gap, i, params, mag) > params.B_sat
```

- [ ] **Step 4: 跑測試確認通過（含相容性鐵律）**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_magnetics_bh.py -v`
Expected: 全數 PASS

再跑全套確認沒破壞既有行為，**特別是特徵測試**：

Run: `.venv/bin/python -m pytest solenoid_model/tests/ -q`
Expected: 全數 PASS，`test_dynamics_characterization.py` 綠燈

- [ ] **Step 5: Commit**

```bash
git add solenoid_model/magnetics.py solenoid_model/tests/test_magnetics_bh.py
git commit -m "feat(magnetics): solve the circuit with a saturating core

Adds solve_flux_density (bisection on the MMF balance) and a mag=None
argument on flux_density/magnetic_force_closing/saturation_check. The
default path is bit-identical to the linear model, so the dry-run
characterization freeze is untouched.

Corrects an unphysical result: the linear model reports 8.16 T at the
closed gap and a 12x force rise across the stroke; the real figures are
2.03 T and +24%.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: 修正保持電流裕度（drive.py）

**Files:**
- Modify: `solenoid_model/drive.py`
- Test: `solenoid_model/tests/test_drive_bh.py`

**Interfaces:**
- Consumes: Task 2 的 `magnetics.magnetic_force_closing(gap, i, params, mag)`
- Produces:
  - `drive.hold_current(params, margin, mag=None) -> float`
  - `drive.force_margin_at_current(params, i, mag=None) -> float`
  - `drive.current_for_force_margin(gap, margin, params, F_resist, mag=None) -> float`

**背景**：`cases.py` 宣稱 `i_hold=11.44 mA` 給 2.0× 裕度，飽和模型下實際只有 1.60×，需 13.07 mA 才是真 2.0×（spec §3A.3 發現 2）。線性模型可解析反解，飽和模型不行（B 對 i 非線性），故改用二分。

- [ ] **Step 1: 寫失敗測試**

建立 `solenoid_model/tests/test_drive_bh.py`：

```python
import math

import pytest

from solenoid_model import drive, magnetization
from solenoid_model.cases import N2_25BAR_PH_PARAMS

PH = N2_25BAR_PH_PARAMS


def test_mag_none_matches_existing_hold_current():
    """相容性：mag=None 必須與現行解析反解相同。"""
    mag_free = drive.hold_current(PH, margin=2.0)
    assert mag_free == pytest.approx(0.01143954024087197, rel=1e-12)


def test_linear_hold_current_overstates_the_margin():
    """spec 3A.3 發現 2：線性算出的 11.44 mA 在飽和下只有 1.60x。"""
    mag = magnetization.from_params(PH)
    i_linear = drive.hold_current(PH, margin=2.0)
    true_margin = drive.force_margin_at_current(PH, i_linear, mag)
    assert true_margin == pytest.approx(1.60, abs=0.02)


def test_saturating_hold_current_restores_true_margin():
    """真 2.0x 需 13.07 mA。"""
    mag = magnetization.from_params(PH)
    i = drive.hold_current(PH, margin=2.0, mag=mag)
    assert i == pytest.approx(0.01307, abs=5e-5)
    assert drive.force_margin_at_current(PH, i, mag) == pytest.approx(
        2.0, rel=1e-6)


def test_saturating_hold_current_exceeds_linear():
    """飽和使吸力變小，所以需要更多電流。"""
    mag = magnetization.from_params(PH)
    assert drive.hold_current(PH, margin=2.0, mag=mag) > drive.hold_current(
        PH, margin=2.0)


def test_hold_power_penalty_is_reported_correctly():
    """0.037 W -> 0.048 W（+31%）。"""
    mag = magnetization.from_params(PH)
    P_lin = drive.hold_current(PH, 2.0) ** 2 * PH.R_coil_20C
    P_sat = drive.hold_current(PH, 2.0, mag=mag) ** 2 * PH.R_coil_20C
    assert P_lin == pytest.approx(0.0368, abs=5e-4)
    assert P_sat == pytest.approx(0.0481, abs=5e-4)


def test_force_margin_round_trips_under_saturation():
    mag = magnetization.from_params(PH)
    for margin in (1.5, 2.0, 3.0):
        i = drive.hold_current(PH, margin, mag=mag)
        assert drive.force_margin_at_current(PH, i, mag) == pytest.approx(
            margin, rel=1e-6)


def test_hold_current_rejects_nonpositive_margin():
    mag = magnetization.from_params(PH)
    with pytest.raises(ValueError, match="margin"):
        drive.hold_current(PH, margin=0.0, mag=mag)
```

- [ ] **Step 2: 跑測試確認失敗**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_drive_bh.py -v`
Expected: FAIL — `TypeError: hold_current() got an unexpected keyword argument 'mag'`

- [ ] **Step 3: 實作**

修改 `solenoid_model/drive.py`。頂端加 import：

```python
from solenoid_model import dynamics, magnetics
```
改為
```python
from scipy.optimize import brentq

from solenoid_model import dynamics, magnetics
```

**取代** `current_for_force_margin`、`hold_current`、`force_margin_at_current` 三個函式：

```python
def current_for_force_margin(gap, margin, params, F_resist, mag=None):
    """Coil current [A] whose force at `gap` is `margin` * `F_resist`.

    Linear mode inverts the force law directly. Saturating mode cannot:
    B is nonlinear in i, so there is no closed form -- it brackets and
    solves numerically instead. Force is monotonic in current in both
    modes, so the root is unique.
    """
    if margin <= 0.0:
        raise ValueError(f"margin must be positive, got {margin}")
    if F_resist <= 0.0:
        raise ValueError(f"F_resist must be positive, got {F_resist}")
    F_target = margin * F_resist
    if mag is None:
        dL_dx = abs(magnetics.dinductance_dgap(gap, params))
        return math.sqrt(2.0 * F_target / dL_dx)

    def residual(i):
        return magnetics.magnetic_force_closing(gap, i, params, mag) - F_target

    # The linear solution is an underestimate of what a saturating core
    # needs (saturation only removes force), so it is a safe lower
    # bracket. Double upward until the force target is cleared; the
    # B_sat ceiling means an unreachable target must terminate rather
    # than loop, hence the explicit cap.
    lo = math.sqrt(2.0 * F_target
                   / abs(magnetics.dinductance_dgap(gap, params)))
    hi = lo
    for _ in range(60):
        if residual(hi) > 0.0:
            return brentq(residual, lo, hi, xtol=1e-12, rtol=1e-12)
        hi *= 2.0
    raise ValueError(
        f"force {F_target:.3f} N at gap {gap:.2e} m is unreachable: the core "
        f"saturates at B_sat={mag.B_sat} T before the coil can produce it")


def hold_current(params, margin, mag=None):
    """Current [A] needed to hold the armature open at `margin` margin.

    The held-open position is x = x_stroke, so the gap is g0 - x_stroke
    (the CLOSED magnetic gap -- the valve is open when the magnetic
    circuit is shut) and the load is the spring at full compression plus
    static pressure.

    Saturation bites hardest exactly here. The closed gap is where the
    linear model most overstates force, so a hold current solved without
    `mag` buys less margin than it claims: on the N2_25BAR_PH case the
    linear 11.44 mA delivers 1.60x, not the 2.0x it was sized for.
    """
    gap_open = params.g0 - params.x_stroke
    F_resist = (dynamics.spring_force(params.x_stroke, params)
                + params.delta_P * params.A_seat)
    return current_for_force_margin(gap_open, margin, params, F_resist, mag)


def force_margin_at_current(params, i, mag=None):
    """Force margin (F_mag / F_resist) at the held-open position."""
    gap_open = params.g0 - params.x_stroke
    F_resist = (dynamics.spring_force(params.x_stroke, params)
                + params.delta_P * params.A_seat)
    return magnetics.magnetic_force_closing(
        gap_open, i, params, mag) / F_resist
```

同時把 `pull_in_current` 加上 `mag` 透傳：

```python
def pull_in_current(params, margin, mag=None):
    """Current [A] needed to start the armature moving from rest."""
    F_resist = (dynamics.spring_force(0.0, params)
                + params.delta_P * params.A_seat)
    return current_for_force_margin(params.g0, margin, params, F_resist, mag)
```

- [ ] **Step 4: 跑測試確認通過**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_drive_bh.py solenoid_model/tests/ -q`
Expected: 全數 PASS

- [ ] **Step 5: Commit**

```bash
git add solenoid_model/drive.py solenoid_model/tests/test_drive_bh.py
git commit -m "feat(drive): solve hold current against a saturating core

hold_current/pull_in_current/force_margin_at_current take mag=None.
Saturating mode brackets and solves numerically because B is nonlinear
in i.

The correction matters: the linear 11.44 mA hold current on
N2_25BAR_PH buys 1.60x margin, not the 2.0x it was sized for. A true
2.0x needs 13.07 mA, raising hold dissipation 0.037 W -> 0.048 W.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: 佈局推導層 layout.py

**Files:**
- Create: `solenoid_model/layout.py`
- Test: `solenoid_model/tests/test_layout.py`

**Interfaces:**
- Consumes: 無（純函式層，**不得** import `magnetics`／`dynamics`／`fluid`）
- Produces:
  - `yoke_area(params, B_pole=None) -> float`
  - `core_radius(params) -> float`、`pole_diameter(params) -> float`
  - `shell_wall_magnetic(params, R_inner) -> float`
  - `hoop_wall_thickness(P, R_inner, S_yield, SF) -> float`
  - `shell_wall_thickness(params, R_inner, P, S_yield, SF, t_manufacturing) -> WallThickness`（dataclass：`magnetic`／`structural`／`manufacturing`／`adopted`／`reason`）
  - `end_plate_thickness(params) -> float`
  - `spring_geometry(params, d_wire, D_coil, G) -> SpringGeometry`（dataclass：`n_active`／`index`／`tau_max`／`L_solid`／`L_free`）
  - `armature_thickness(params, rho, D) -> float`
  - `armature_mass_consistency(params, t_arm, rho, D) -> tuple`
  - `axial_stack(params, **overrides) -> dict`
  - `envelope(params, ...) -> dict`
  - `LAYOUT_EMPIRICAL: dict` — 每項 `{"value", "unit", "reason", "source"}`

- [ ] **Step 1: 寫失敗測試**

建立 `solenoid_model/tests/test_layout.py`：

```python
import math
from dataclasses import replace

import pytest

from solenoid_model import layout, sizing
from solenoid_model.cases import N2_25BAR_PH_PARAMS

PH = N2_25BAR_PH_PARAMS


def test_layout_does_not_import_physics_layers():
    """分層不變式：layout 是純代數層。"""
    import solenoid_model.layout as m
    src = open(m.__file__, encoding="utf-8").read()
    for banned in ("import dynamics", "import magnetics", "import fluid",
                   "from solenoid_model import dynamics",
                   "from solenoid_model import magnetics",
                   "from solenoid_model import fluid"):
        assert banned not in src


def test_yoke_area_at_B_sat_equals_gap_area():
    """磁通連續 + B_sat 鉗位：同材質下 A_yoke >= A_gap。"""
    assert layout.yoke_area(PH) == pytest.approx(PH.A_gap, rel=1e-12)


def test_yoke_area_scales_with_real_pole_flux():
    """真實 B 較低時所需截面等比例縮小（spec 3: 1.824 T -> 16.97 mm2）。"""
    A = layout.yoke_area(PH, B_pole=1.824)
    assert A == pytest.approx(16.97e-6, rel=1e-3)
    assert A < layout.yoke_area(PH)


def test_yoke_area_rejects_B_pole_above_B_sat():
    with pytest.raises(ValueError, match="B_sat"):
        layout.yoke_area(PH, B_pole=3.0)


def test_core_radius_and_pole_diameter_match_A_gap():
    assert layout.pole_diameter(PH) == pytest.approx(5.05e-3, abs=5e-6)
    assert layout.core_radius(PH) == pytest.approx(
        layout.pole_diameter(PH) / 2, rel=1e-12)


def test_shell_wall_magnetic_gives_annulus_of_yoke_area():
    R_i = 9.0e-3
    t = layout.shell_wall_magnetic(PH, R_i)
    assert t == pytest.approx(0.35e-3, abs=1e-5)
    A_annulus = math.pi * ((R_i + t) ** 2 - R_i ** 2)
    assert A_annulus == pytest.approx(layout.yoke_area(PH), rel=1e-9)


def test_hoop_wall_thickness_matches_thin_wall_formula():
    t = layout.hoop_wall_thickness(P=25e5, R_inner=9.0e-3,
                                   S_yield=205e6, SF=3.0)
    assert t == pytest.approx(329.3e-6, abs=1e-6)


def test_hoop_wall_scales_with_safety_factor():
    a = layout.hoop_wall_thickness(25e5, 9.0e-3, 205e6, 3.0)
    b = layout.hoop_wall_thickness(25e5, 9.0e-3, 205e6, 6.0)
    assert b == pytest.approx(2.0 * a, rel=1e-12)


def test_shell_wall_adopts_the_largest_of_three():
    w = layout.shell_wall_thickness(PH, R_inner=9.0e-3, P=25e5,
                                    S_yield=205e6, SF=3.0,
                                    t_manufacturing=1.0e-3)
    assert w.adopted == pytest.approx(1.0e-3, rel=1e-12)
    assert w.adopted == max(w.magnetic, w.structural, w.manufacturing)
    assert "manufacturing" in w.reason


def test_shell_wall_reason_names_the_binding_constraint():
    w = layout.shell_wall_thickness(PH, R_inner=9.0e-3, P=250e5,
                                    S_yield=205e6, SF=3.0,
                                    t_manufacturing=0.1e-3)
    assert w.adopted == pytest.approx(w.structural, rel=1e-12)
    assert "structural" in w.reason


def test_end_plate_thickness_carries_yoke_flux_at_core_radius():
    t = layout.end_plate_thickness(PH)
    assert t == pytest.approx(1.26e-3, abs=1e-5)
    A_radial = 2 * math.pi * layout.core_radius(PH) * t
    assert A_radial == pytest.approx(layout.yoke_area(PH), rel=1e-9)


def test_spring_geometry_round_trips_to_k_spring():
    """由 (d, D, n) 正算回 k 必須等於 params.k_spring。"""
    g = layout.spring_geometry(PH, d_wire=0.35e-3, D_coil=1.8e-3, G=79e9)
    k_back = 79e9 * (0.35e-3) ** 4 / (8 * (1.8e-3) ** 3 * g.n_active)
    assert k_back == pytest.approx(PH.k_spring, rel=1e-9)


def test_spring_geometry_matches_verified_design_point():
    """spec 6 已驗證：n=6.35, C=5.14, tau=389 MPa, L_free=3.72 mm。"""
    g = layout.spring_geometry(PH, d_wire=0.35e-3, D_coil=1.8e-3, G=79e9)
    assert g.n_active == pytest.approx(6.35, abs=0.02)
    assert g.index == pytest.approx(5.14, abs=0.02)
    assert g.tau_max == pytest.approx(389e6, rel=0.02)
    assert g.L_free == pytest.approx(3.72e-3, abs=5e-5)


def test_spring_wahl_factor_approaches_one_for_large_index():
    """C -> infinity 時 Wahl 修正 -> 1。"""
    assert layout.wahl_factor(1e6) == pytest.approx(1.0, abs=1e-5)
    assert layout.wahl_factor(5.0) > 1.2


def test_spring_solid_length_is_below_free_length():
    g = layout.spring_geometry(PH, d_wire=0.35e-3, D_coil=1.8e-3, G=79e9)
    assert g.L_solid < g.L_free


def test_armature_thickness_derives_from_declared_mass():
    """銜鐵厚度由 m_arm 反解，非 EMPIRICAL。"""
    t = layout.armature_thickness(PH, rho=7870.0,
                                  D=layout.pole_diameter(PH))
    assert t == pytest.approx(5.08e-3, abs=2e-5)


def test_armature_mass_consistency_flags_mismatch():
    """4.0 mm 只給 0.631 g，與宣告的 0.8 g 不符。"""
    m_dec, m_geo, rel = layout.armature_mass_consistency(
        PH, t_arm=4.0e-3, rho=7870.0, D=layout.pole_diameter(PH))
    assert m_dec == pytest.approx(0.8e-3, rel=1e-12)
    assert m_geo == pytest.approx(0.631e-3, abs=5e-6)
    assert rel > 0.2


def test_armature_mass_consistency_passes_at_derived_thickness():
    t = layout.armature_thickness(PH, rho=7870.0, D=layout.pole_diameter(PH))
    _, _, rel = layout.armature_mass_consistency(
        PH, t_arm=t, rho=7870.0, D=layout.pole_diameter(PH))
    assert rel < 1e-9


def test_axial_stack_total_equals_sum_of_segments():
    stack = layout.axial_stack(PH)
    assert stack["total"] == pytest.approx(
        sum(v for k, v in stack["segments"]), rel=1e-12)


def test_axial_stack_matches_verified_total():
    """spec 7：合計 18.68 mm。"""
    assert layout.axial_stack(PH)["total"] == pytest.approx(18.679e-3, abs=2e-5)


def test_axial_stack_includes_stroke_and_working_gap_from_params():
    stack = dict(layout.axial_stack(PH)["segments"])
    assert stack["行程 x_stroke"] == pytest.approx(PH.x_stroke, rel=1e-12)
    assert stack["工作氣隙 g0"] == pytest.approx(PH.g0, rel=1e-12)


def test_envelope_outer_diameter_uses_adopted_wall():
    env = layout.envelope(PH, coil_OD=18.0e-3, t_manufacturing=1.0e-3)
    assert env["OD"] == pytest.approx(20.0e-3, abs=1e-5)


def test_lengths_scale_linearly_under_self_similar_scaling():
    """自相似縮放下長度 ~ s、截面 ~ s^2。"""
    s = 2.0
    scaled = sizing.scale_params(PH, s)
    assert layout.core_radius(scaled) == pytest.approx(
        s * layout.core_radius(PH), rel=1e-9)
    assert layout.yoke_area(scaled) == pytest.approx(
        s ** 2 * layout.yoke_area(PH), rel=1e-9)


def test_every_empirical_entry_is_fully_declared():
    assert layout.LAYOUT_EMPIRICAL
    for name, entry in layout.LAYOUT_EMPIRICAL.items():
        for field in ("value", "unit", "reason", "source"):
            assert field in entry, f"{name} missing {field}"
        assert entry["source"] == "EMPIRICAL"
        assert entry["reason"].strip()
```

- [ ] **Step 2: 跑測試確認失敗**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_layout.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'solenoid_model.layout'`

- [ ] **Step 3: 實作**

建立 `solenoid_model/layout.py`：

```python
"""Physical layout: from ValveParams to a manufacturable dimension chain (P8).

Never imports dynamics/magnetics/fluid -- the same invariant sizing.py and
winding.py hold. Everything here is closed-form algebra over an existing
ValveParams; nothing feeds back into the ODE.

The split this module is careful about is DERIVED vs EMPIRICAL. Yoke
sections, wall thicknesses, spring geometry and armature thickness all
follow from physics already in ValveParams. Seat-body and fixed-pole
thicknesses do not -- they are packaging judgements, declared in
LAYOUT_EMPIRICAL with a stated reason rather than buried as literals, so
that reading a dimension tells you how much to trust it.
"""
import math
from dataclasses import dataclass

MU_0 = 4.0 * math.pi * 1e-7  # vacuum permeability, H/m
RHO_FE = 7870.0              # kg/m^3, soft iron

# EMPIRICAL layout inputs: packaging judgements, not derived quantities.
# Declared here rather than defaulted inside functions so that changing a
# judgement is a visible edit, following how params.py labels k_pack/C_d.
LAYOUT_EMPIRICAL = {
    "t_seat_body": {
        "value": 3.0e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "容納 0.60 mm 孔、密封 land 與閥座壓入配合的最小座體厚度；"
                  "無一致性檢查可攔截，須由 CAD 階段回頭確認",
    },
    "t_fixed_pole": {
        "value": 4.0e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "固定極需容納磁通轉向與繞線骨架端面支撐；厚度不足會在極面"
                  "根部先飽和，但本模型未解 3D 磁通分佈，故取封裝經驗值",
    },
    "t_manufacturing": {
        "value": 1.0e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "外殼可加工與剛性下限。磁性下限 0.35 mm 與結構下限 0.33-0.44 mm "
                  "皆遠低於此，故壁厚實際由製造決定而非物理",
    },
    "t_sleeve": {
        "value": 0.25e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "濕式銜鐵的非磁性隔離套（316L）。25 bar 下結構僅需 0.11 mm，"
                  "餘為加工餘裕。⚠️ 此厚度落在徑向磁路上但 magnetics.reluctance "
                  "未建模，故實際吸力低於模型報告值",
    },
    "clearance_spring": {
        "value": 0.1e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "彈簧全壓縮時與腔壁的餘隙，避免併圈干涉",
    },
}


# --- Magnetic circuit sections -------------------------------------------

def core_radius(params):
    """Core / pole radius [m] from the pole face area."""
    return math.sqrt(params.A_gap / math.pi)


def pole_diameter(params):
    """Pole face diameter [m]."""
    return 2.0 * core_radius(params)


def yoke_area(params, B_pole=None):
    """Minimum yoke cross-section [m^2] that will not saturate before the pole.

    Flux continuity around the circuit: the yoke carries the same flux the
    pole face does, so

        B_pole * A_gap = B_yoke * A_yoke,  B_yoke <= B_sat
        => A_yoke >= A_gap * B_pole / B_sat

    `B_pole=None` uses the worst case B_pole = B_sat, giving
    A_yoke >= A_gap. That is deliberately the shipped default: it is
    independent of the operating point, so it stays valid when the bus
    voltage or turns count changes later. Passing the actual solved pole
    flux (magnetics.solve_flux_density) gives the tighter,
    operating-point-specific bound -- 16.97 mm^2 against the conservative
    20.0 mm^2 on the N2_25BAR_PH case.

    Assumes yoke and pole share a material. For a mixed-material circuit,
    scale by the ratio of the two B_sat values.
    """
    if B_pole is None:
        return params.A_gap
    if B_pole <= 0.0:
        raise ValueError(f"B_pole={B_pole} must be positive")
    if B_pole > params.B_sat:
        raise ValueError(
            f"B_pole={B_pole} T exceeds B_sat={params.B_sat} T: the pole "
            f"cannot carry that flux density, so sizing a yoke to it is "
            f"meaningless")
    return params.A_gap * B_pole / params.B_sat


def shell_wall_magnetic(params, R_inner, B_pole=None):
    """Shell wall thickness [m] whose annulus carries the yoke flux.

    Solves pi*((R_i+t)^2 - R_i^2) = A_yoke for t.
    """
    A = yoke_area(params, B_pole)
    return math.sqrt(A / math.pi + R_inner ** 2) - R_inner


def end_plate_thickness(params, B_pole=None):
    """End plate thickness [m] for radially outward flux.

    A washer carrying flux radially has cross-section A = 2*pi*r*t, which
    is SMALLEST at the smallest radius. The tightest point is therefore
    the core outer surface, and sizing there covers the whole plate.
    """
    A = yoke_area(params, B_pole)
    return A / (2.0 * math.pi * core_radius(params))


# --- Pressure boundary ----------------------------------------------------

def hoop_wall_thickness(P, R_inner, S_yield, SF):
    """Thin-wall hoop-stress thickness [m]: t = P*R/(S_yield/SF)."""
    if SF <= 0.0:
        raise ValueError(f"SF={SF} must be positive")
    if S_yield <= 0.0:
        raise ValueError(f"S_yield={S_yield} must be positive")
    return P * R_inner / (S_yield / SF)


@dataclass(frozen=True)
class WallThickness:
    """Three independent lower bounds on a wall, and which one binds."""
    magnetic: float
    structural: float
    manufacturing: float
    adopted: float
    reason: str


def shell_wall_thickness(params, R_inner, P, S_yield, SF,
                         t_manufacturing=None, B_pole=None):
    """Shell wall [m] as the largest of three independent lower bounds.

    Reporting all three rather than just the winner is the point: on this
    valve the magnetic (0.35 mm) and structural (0.33-0.44 mm) bounds are
    the same order and BOTH are below what can be machined, so the wall is
    set by manufacturing, not by physics. A single returned number would
    hide that.
    """
    if t_manufacturing is None:
        t_manufacturing = LAYOUT_EMPIRICAL["t_manufacturing"]["value"]
    t_mag = shell_wall_magnetic(params, R_inner, B_pole)
    t_str = hoop_wall_thickness(P, R_inner, S_yield, SF)
    candidates = {"magnetic": t_mag, "structural": t_str,
                  "manufacturing": t_manufacturing}
    name = max(candidates, key=candidates.get)
    return WallThickness(
        magnetic=t_mag, structural=t_str, manufacturing=t_manufacturing,
        adopted=candidates[name],
        reason=f"{name} bound governs "
               f"(magnetic {t_mag*1e3:.2f} mm, structural {t_str*1e3:.2f} mm, "
               f"manufacturing {t_manufacturing*1e3:.2f} mm)")


# --- Spring ---------------------------------------------------------------

def wahl_factor(C):
    """Wahl stress-correction factor for spring index C = D/d.

    Corrects for curvature and direct shear, both of which raise the real
    stress above the straight-torsion estimate. Tends to 1 as C grows.
    """
    if C <= 1.0:
        raise ValueError(f"spring index C={C} must exceed 1")
    return (4.0 * C - 1.0) / (4.0 * C - 4.0) + 0.615 / C


@dataclass(frozen=True)
class SpringGeometry:
    """Helical compression spring realising params.k_spring."""
    d_wire: float     # m
    D_coil: float     # m
    n_active: float   # active coils, -
    index: float      # C = D/d, -
    tau_max: float    # Pa, at full stroke, Wahl-corrected
    L_solid: float    # m
    L_free: float     # m


def spring_geometry(params, d_wire, D_coil, G, n_dead=2.0, clearance=None):
    """Solve the coil count that realises params.k_spring, and its stresses.

    k = G*d^4/(8*D^3*n)  =>  n = G*d^4/(8*D^3*k)

    Free length is solid height plus the working deflection (preload plus
    stroke) plus a clearance, so the spring never reaches solid at full
    stroke.
    """
    if clearance is None:
        clearance = LAYOUT_EMPIRICAL["clearance_spring"]["value"]
    if d_wire <= 0.0 or D_coil <= 0.0:
        raise ValueError(f"d_wire={d_wire}, D_coil={D_coil} must be positive")
    C = D_coil / d_wire
    n = G * d_wire ** 4 / (8.0 * D_coil ** 3 * params.k_spring)
    if n <= 0.0:
        raise ValueError(f"k_spring={params.k_spring} gives n={n} coils")
    F_max = params.F_preload + params.k_spring * params.x_stroke
    tau = wahl_factor(C) * 8.0 * F_max * D_coil / (math.pi * d_wire ** 3)
    L_solid = (n + n_dead) * d_wire
    deflection = params.F_preload / params.k_spring + params.x_stroke
    return SpringGeometry(
        d_wire=d_wire, D_coil=D_coil, n_active=n, index=C, tau_max=tau,
        L_solid=L_solid, L_free=L_solid + deflection + clearance)


# --- Armature -------------------------------------------------------------

def armature_thickness(params, rho=RHO_FE, D=None):
    """Armature disc thickness [m] implied by params.m_arm.

    DERIVED, not a packaging guess: t = m_arm / (rho * A). Defaults to a
    disc at the pole diameter.
    """
    if D is None:
        D = pole_diameter(params)
    A = math.pi * D ** 2 / 4.0
    return params.m_arm / (rho * A)


def armature_mass_consistency(params, t_arm, rho=RHO_FE, D=None):
    """Compare declared m_arm against the mass a given thickness implies.

    Reports rather than raises, following
    winding.coil_resistance_consistency: an armature may legitimately stop
    being a plain disc (a guide stem, a lightening bore), which decouples
    thickness from mass. The caller decides what divergence is acceptable.

    Returns (m_declared, m_geometric, rel_error).
    """
    if D is None:
        D = pole_diameter(params)
    A = math.pi * D ** 2 / 4.0
    m_geom = rho * A * t_arm
    return params.m_arm, m_geom, abs(m_geom - params.m_arm) / params.m_arm


# --- Assembly -------------------------------------------------------------

def axial_stack(params, t_seat_body=None, t_fixed_pole=None,
                d_wire=0.35e-3, D_coil=1.8e-3, G=79e9, rho=RHO_FE):
    """Ordered axial dimension chain [m], seat face upward.

    Returns {"segments": [(name, thickness), ...], "total": float}. Ordered
    so the list reads as the physical stack, not as a dict of parts.
    """
    if t_seat_body is None:
        t_seat_body = LAYOUT_EMPIRICAL["t_seat_body"]["value"]
    if t_fixed_pole is None:
        t_fixed_pole = LAYOUT_EMPIRICAL["t_fixed_pole"]["value"]
    spring = spring_geometry(params, d_wire, D_coil, G)
    segments = [
        ("閥座座體", t_seat_body),
        ("行程 x_stroke", params.x_stroke),
        ("銜鐵", armature_thickness(params, rho)),
        ("工作氣隙 g0", params.g0),
        ("固定極", t_fixed_pole),
        ("彈簧腔", spring.L_free),
        ("端板 x2", 2.0 * end_plate_thickness(params)),
    ]
    return {"segments": segments, "total": sum(t for _, t in segments)}


def envelope(params, coil_OD, R_inner=None, P=25e5, S_yield=205e6, SF=3.0,
             t_manufacturing=None, **stack_kwargs):
    """Overall package envelope [m]: outer diameter and length."""
    if R_inner is None:
        R_inner = coil_OD / 2.0
    wall = shell_wall_thickness(params, R_inner, P, S_yield, SF,
                                t_manufacturing)
    stack = axial_stack(params, **stack_kwargs)
    return {"OD": 2.0 * (R_inner + wall.adopted),
            "L": stack["total"],
            "wall": wall,
            "stack": stack}
```

- [ ] **Step 4: 跑測試確認通過**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_layout.py -v`
Expected: 全數 PASS

- [ ] **Step 5: Commit**

```bash
git add solenoid_model/layout.py solenoid_model/tests/test_layout.py
git commit -m "feat(layout): derive the physical dimension chain

Yoke sections from flux continuity, wall thickness as the largest of
magnetic/structural/manufacturing bounds, spring geometry solved from
k_spring, and armature thickness derived from m_arm with a consistency
check. Packaging judgements live in LAYOUT_EMPIRICAL with stated
reasons rather than as bare literals.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: 剖面圖與架構圖

**Files:**
- Create: `solenoid_model/report/generate_layout.py`
- Test: `solenoid_model/tests/test_layout_report.py`

**Interfaces:**
- Consumes: Task 4 的 `layout.envelope`／`axial_stack`／`spring_geometry`／`pole_diameter`／`core_radius`
- Produces:
  - `generate_section() -> Path`（輸出 `docs/spec/layout_section.png`）
  - `generate_architecture() -> Path`（輸出 `docs/spec/layout_architecture.png`）
  - 模組常數 `SECTION_OUT`、`ARCH_OUT`

- [ ] **Step 1: 寫失敗測試**

建立 `solenoid_model/tests/test_layout_report.py`：

```python
from solenoid_model.report import generate_layout


def test_generate_section_writes_a_png():
    path = generate_layout.generate_section()
    assert path.exists()
    assert path.suffix == ".png"
    assert path.stat().st_size > 1000


def test_generate_architecture_writes_a_png():
    path = generate_layout.generate_architecture()
    assert path.exists()
    assert path.suffix == ".png"
    assert path.stat().st_size > 1000


def test_importing_does_not_mutate_global_rcparams():
    """中文字型必須以 rc_context 區域套用，不得污染全域。"""
    import matplotlib
    before = list(matplotlib.rcParams["font.sans-serif"])
    import importlib
    importlib.reload(generate_layout)
    assert list(matplotlib.rcParams["font.sans-serif"]) == before
```

- [ ] **Step 2: 跑測試確認失敗**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_layout_report.py -v`
Expected: FAIL — `ImportError: cannot import name 'generate_layout'`

- [ ] **Step 3: 實作**

建立 `solenoid_model/report/generate_layout.py`：

```python
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
        y = 0.0
        # 閥座座體
        ax.add_patch(Rectangle((0, _mm(y)), _mm(R_shell + wall),
                               _mm(t_seat), fc=COL_NONMAG, ec="k"))
        ax.text(_mm(R_shell) * 0.55, _mm(y + t_seat / 2), "閥座座體",
                ha="center", va="center", fontsize=8)
        # 流道孔
        ax.add_patch(Rectangle((0, _mm(y)), _mm(params.D_seat_bore / 2),
                               _mm(t_seat), fc="w", ec=COL_SEAL, lw=1.5))
        y += t_seat
        y += params.x_stroke      # 行程（閥開時的間隙）
        # 銜鐵
        ax.add_patch(Rectangle((0, _mm(y)), _mm(r_core), _mm(t_arm),
                               fc=COL_MAG, ec="k"))
        ax.text(_mm(r_core) / 2, _mm(y + t_arm / 2), "銜鐵",
                ha="center", va="center", fontsize=8, color="w")
        y += t_arm
        y += params.g0            # 工作氣隙
        # 固定極
        ax.add_patch(Rectangle((0, _mm(y)), _mm(r_core), _mm(t_pole),
                               fc=COL_MAG, ec="k"))
        ax.text(_mm(r_core) / 2, _mm(y + t_pole / 2), "固定極",
                ha="center", va="center", fontsize=8, color="w")
        # 線圈（在固定極與銜鐵徑向外側）
        y_coil = y - COIL_H / 2
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
        spring = layout.spring_geometry(params, 0.35e-3, 1.8e-3, 79e9)
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
        ax.set_xlim(-3, _mm(R_shell + wall) + 4)
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
            ax.text(x + w / 2, y + h - 0.35, label, ha="center",
                    va="top", fontsize=10, weight="bold")

        def node(x, y, text):
            ax.add_patch(Rectangle((x, y), 2.0, 0.7, fc="w", ec="k"))
            ax.text(x + 1.0, y + 0.35, text, ha="center", va="center",
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

        box(0.3, 2.4, 4.6, 1.9, "機械域", COL_MECH)
        node(0.6, 3.3, "銜鐵 0.8 g"); node(2.7, 3.3, "彈簧 4000 N/m")
        node(1.6, 2.55, "閥座止擋")
        arrow((2.6, 3.65), (2.7, 3.65))

        box(5.1, 2.4, 4.6, 1.9, "流體域", COL_SEAL)
        node(5.4, 3.3, "N₂ 25 bar"); node(7.5, 3.3, "孔徑 0.60 mm")
        node(6.4, 2.55, "1.28 g/s 壅塞流")
        arrow((7.4, 3.65), (7.5, 3.65))

        # 跨域耦合點
        arrow((7.6, 5.6), (7.6, 4.3), color=COL_MAG)
        ax.text(7.75, 4.9, "耦合①\n氣隙：磁↔機械", fontsize=8, color=COL_MAG)
        arrow((2.6, 2.55), (5.4, 2.55), color=COL_SEAL, style="<->")
        ax.text(3.0, 2.2, "耦合②　閥座：機械↔流體", fontsize=8, color=COL_SEAL)

        ax.text(5.0, 0.9,
                "耦合①：氣隙同時決定磁阻與機械位置（dynamics.coupled_rhs）\n"
                "耦合②：閥座開度決定流量，噴流反作用力回饋進力平衡",
                ha="center", fontsize=8.5)
        ax.set_xlim(0, 10); ax.set_ylim(0.4, 9.2)
        ax.axis("off")
        ax.set_title("電磁閥架構圖：四個物理域與跨域耦合點", fontsize=12)
        ARCH_OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(ARCH_OUT, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return ARCH_OUT


if __name__ == "__main__":
    print(generate_section())
    print(generate_architecture())
```

- [ ] **Step 4: 跑測試並目視檢查產出的圖**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_layout_report.py -v`
Expected: 全數 PASS

Run: `.venv/bin/python -m solenoid_model.report.generate_layout`
Expected: 印出兩個路徑

**目視檢查**（用 Read 工具開啟 `docs/spec/layout_section.png` 與
`docs/spec/layout_architecture.png`）：確認無中文豆腐字、無元件重疊、
比例看起來合理（外徑約為總長的一半）。若有重疊，調整座標後重跑。

- [ ] **Step 5: Commit**

```bash
git add solenoid_model/report/generate_layout.py \
        solenoid_model/tests/test_layout_report.py \
        docs/spec/layout_section.png docs/spec/layout_architecture.png
git commit -m "feat(layout): add to-scale section and architecture figures

Unlike valve_schematic.png (a not-to-scale parameter key), the section
is drawn at true mm scale from layout.py's dimension chain.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: GUI 第五頁籤

**Files:**
- Modify: `solenoid_model/app.py`（`main()` 的 `st.tabs`，約 `app.py:512-513`）
- Test: `solenoid_model/tests/test_app_layout.py`

**Interfaces:**
- Consumes: Task 4 的 `layout` 全部、Task 5 的 `SECTION_OUT`／`ARCH_OUT`
- Produces: `app.render_layout_tab(params)`

- [ ] **Step 1: 寫失敗測試**

建立 `solenoid_model/tests/test_app_layout.py`：

```python
import pytest

from solenoid_model import app, layout
from solenoid_model.cases import N2_25BAR_PH_PARAMS


def test_render_layout_tab_exists():
    assert hasattr(app, "render_layout_tab")


def test_compute_layout_returns_full_chain():
    result = app.compute_layout(N2_25BAR_PH_PARAMS)
    assert result["envelope"]["OD"] == pytest.approx(20.0e-3, abs=1e-5)
    assert result["envelope"]["L"] == pytest.approx(18.679e-3, abs=2e-5)
    assert result["wall"].adopted == pytest.approx(1.0e-3, rel=1e-12)
    assert result["spring"].n_active == pytest.approx(6.35, abs=0.02)


def test_compute_layout_reports_armature_consistency():
    result = app.compute_layout(N2_25BAR_PH_PARAMS)
    assert result["armature_rel_error"] < 1e-9


def test_compute_layout_surfaces_empirical_entries():
    result = app.compute_layout(N2_25BAR_PH_PARAMS)
    assert result["empirical"] is layout.LAYOUT_EMPIRICAL
```

- [ ] **Step 2: 跑測試確認失敗**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_app_layout.py -v`
Expected: FAIL — `AttributeError: module 'solenoid_model.app' has no attribute 'render_layout_tab'`

- [ ] **Step 3: 實作**

在 `solenoid_model/app.py` 的 import 區加入：

```python
from solenoid_model import layout
from solenoid_model.report import generate_layout
```

新增純函式與頁籤（放在 `render_sealing_tab` 之後）：

```python
def compute_layout(params, coil_OD=18.0e-3):
    """Pure function: the whole physical dimension chain for one params."""
    env = layout.envelope(params, coil_OD=coil_OD)
    spring = layout.spring_geometry(params, 0.35e-3, 1.8e-3, 79e9)
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
    st.write(f"線徑 {s.d_wire*1e3:.2f} mm / 中徑 {s.D_coil*1e3:.2f} mm / "
             f"有效圈數 {s.n_active:.2f} / 彈簧指數 {s.index:.2f} / "
             f"最大剪應力 {s.tau_max/1e6:.0f} MPa / "
             f"自由長 {s.L_free*1e3:.2f} mm")

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
                      (generate_layout.ARCH_OUT, "架構圖")):
        if path.exists():
            st.image(str(path), caption=cap)
```

**修改 `app.py:512-513`**，把四頁改為五頁：

```python
    tab_open, tab_close, tab_limits, tab_seal = st.tabs(
        ["開啟動態", "關閉/續流", "極限掃描 L1–L3", "密封/壽命 L4–L5"])
```
改為
```python
    tab_open, tab_close, tab_limits, tab_seal, tab_layout = st.tabs(
        ["開啟動態", "關閉/續流", "極限掃描 L1–L3", "密封/壽命 L4–L5",
         "實體佈局"])
```

並在 `with tab_seal:` 區塊之後加入：

```python
    with tab_layout:
        try:
            render_layout_tab(params)
        except ValueError as e:
            st.error(f"佈局計算失敗：{e}")
```

- [ ] **Step 4: 跑測試確認通過**

Run: `.venv/bin/python -m pytest solenoid_model/tests/test_app_layout.py -v`
Expected: 全數 PASS

Run: `.venv/bin/python -m pytest solenoid_model/tests/ -q`
Expected: 全數 PASS（含既有 `test_app.py` 的 AppTest UI 測試）

- [ ] **Step 5: Commit**

```bash
git add solenoid_model/app.py solenoid_model/tests/test_app_layout.py
git commit -m "feat(gui): add the physical layout tab

compute_layout as a pure function plus render_layout_tab, following the
existing compute_*/render_*_tab split. Surfaces all three wall-thickness
bounds and the EMPIRICAL table rather than just the adopted numbers.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: 文件同步

**Files:**
- Modify: `ARCHITECTURE.md`、`ROADMAP.md`、`solenoid_model/report/MODEL_NOTES.md`、`docs/spec/PARAMS.md`、`USAGE.md`

**Interfaces:**
- Consumes: Task 1–6 全部
- Produces: 無程式介面

**背景**：`ARCHITECTURE.md` 開頭寫明「新增模組時同步更新此檔案」，`ROADMAP.md`
自稱工作項目編號的唯一定義處。

- [ ] **Step 1: 更新 ARCHITECTURE.md**

在資料夾樹加入三個新模組（`magnetization.py`、`layout.py`、
`report/generate_layout.py`）與三個新測試檔。

在「分層架構與雙模式設計」段落，把既有的「四個正交旋鈕」敘述改為**五個**，
新增一段：

```markdown
磁化模型沿用同一原則：`mag=None`（預設）維持既有線性磁路，對現有呼叫者
bit-level 相同；傳入 `magnetization.AnalyticBH` 或 `TabulatedBH` 才切換到
飽和磁路（`magnetics.solve_flux_density` 對 B 二分求解 MMF 平衡）。連同
`medium`/`cond`、`T_coil` 與 `seat`，`mag` 是**第五個正交模式旋鈕**。

`layout.py` 屬第 2 層（物理函式層）的純代數分支，與 `sizing.py`／`winding.py`
同級：不 import `dynamics`／`magnetics`／`fluid`，只從 `ValveParams` 推導
實體尺寸鏈。
```

- [ ] **Step 2: 更新 ROADMAP.md**

在模組總表末尾加入 P8 列：

```markdown
| P8 | B-H 飽和曲線 + 實體佈局 | `feat(magnetics)`／`feat(layout)` | `magnetization.py`（Fröhlich-Kennelly 解析式 + 查表，零新增自由度）；`magnetics.solve_flux_density`（對 B 二分解 MMF 平衡）+ `mag=None` 第五旋鈕（線性路徑 bit-level 不變）；`drive.py` 全面支援 `mag`；`layout.py`（磁通連續定磁軛截面、壁厚取磁性/結構/製造三下限之最大、彈簧封閉解、銜鐵厚度由 m_arm 反解 + 一致性檢查、LAYOUT_EMPIRICAL 集中宣告）；`generate_layout.py` 兩張圖（按比例剖面 + 架構）；GUI 第五頁籤。**發現 1**：線性模型「吸力隨氣隙關閉暴增」是假象——真實吸力沿行程僅 +24%（26.5→32.7 N），線性報 12 倍（44.6→530 N），影響 P6 撞擊速度與 L5 壽命估算；**發現 2**：`N2_25BAR_PH` 宣稱的 2.0× 保持裕度實為 **1.60×**，真 2.0× 需 13.07 mA（非 11.44 mA），保持功耗 0.037→0.048 W；**發現 3**：外殼壁厚的磁性下限 0.35 mm 與結構下限 0.33–0.44 mm 皆低於可加工厚度，壁厚實由製造決定；**刻意不做**：磁滯迴線（殘留吸力 ~0.02 N vs 彈簧預載 2.2 N，差兩個數量級） |
```

- [ ] **Step 3: 更新 MODEL_NOTES.md**

在已知簡化項清單加入：

```markdown
## P8 已知限制

1. **隔離套磁阻未計入**：`magnetics.reluctance` 只有軸向氣隙 + 鐵芯路徑，
   濕式銜鐵的徑向非磁性隔離套（0.25 mm）等效氣隙未建模 → 實際吸力低於
   模型報告值。
2. **B_sat 未含溫度效應**：B-H 曲線沿用常數 `B_sat`，高溫下實際飽和磁通
   較低（`thermal.curie_margin` 已標註此項未建模）→ 高溫吸力被高估。
3. **單值曲線不含損耗**：無磁滯損耗、無渦流。快速暫態下渦流會延遲磁通
   建立（P4 已列「仍缺渦流延遲修正」），本層未改善。
4. **既有報告數字仍為線性模式**：`cases.py` 註記與 P0–P7 各報告的裕度
   數字皆在 `mag=None` 下產生，未以飽和模型重算。引用時須註明模式。
5. **軸向尺寸鏈部分為 EMPIRICAL**：閥座座體與固定極厚度無一致性檢查可
   攔截（銜鐵厚度有）。
6. **未含公差與配合**：所有尺寸為標稱值。
```

- [ ] **Step 4: 更新 PARAMS.md 與 USAGE.md**

`PARAMS.md`：新增一段說明 `mag` 旋鈕與 `LAYOUT_EMPIRICAL` 各項的物理意義，
並指向 `docs/spec/layout_section.png`。

`USAGE.md`：新增「實體佈局」頁籤的操作說明與 `generate_layout.py` 的執行方式：

```bash
.venv/bin/python -m solenoid_model.report.generate_layout
```

- [ ] **Step 5: 跑全套測試並 Commit**

Run: `.venv/bin/python -m pytest solenoid_model/tests/ -q`
Expected: 全數 PASS

```bash
git add ARCHITECTURE.md ROADMAP.md solenoid_model/report/MODEL_NOTES.md \
        docs/spec/PARAMS.md USAGE.md
git commit -m "docs: record P8 in the roadmap, architecture and model notes

Registers the fifth orthogonal mode knob (mag) and layout.py's layer
position, and lists the P8 limitations -- chief among them that the
existing reports' margin figures were all produced in linear mode and
have not been recomputed.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Self-Review

**1. Spec 覆蓋**

| Spec 章節 | 對應任務 |
|---|---|
| §3A.1 線性模型問題 | Task 2（測試對照 8.16 T） |
| §3A.2 Fröhlich-Kennelly | Task 1 |
| §3A.3 設計衝擊（三發現） | Task 2（發現 1、3）、Task 3（發現 2） |
| §3A.4 第五旋鈕介面 | Task 1 + Task 2 |
| §3A.5 磁滯不做 | Task 1 docstring + Task 7 ROADMAP |
| §3 磁軛截面 | Task 4 `yoke_area` |
| §4 壁厚取大者 | Task 4 `shell_wall_thickness` |
| §5 承壓邊界／隔離套 | Task 4 `LAYOUT_EMPIRICAL["t_sleeve"]`、Task 6 警示 |
| §6 彈簧 | Task 4 `spring_geometry` |
| §7 軸向尺寸鏈 | Task 4 `axial_stack` |
| §8 模組設計 | Task 1、2、4 |
| §9 兩張圖 | Task 5 |
| §10 GUI | Task 6 |
| §11 測試 | Task 1–6 各自的測試檔 |
| §12 已知限制 | Task 7 MODEL_NOTES.md |

無缺口。

**2. Placeholder 掃描**：無 TBD／TODO／「類似 Task N」；每個程式步驟都有完整可貼上的程式碼。

**3. 型別一致性**：`mag` 引數名在 `magnetization`／`magnetics`／`drive` 三處一致；`layout.WallThickness` 的欄位名（`magnetic`／`structural`／`manufacturing`／`adopted`／`reason`）在 Task 4 定義、Task 6 使用一致；`SpringGeometry` 的 `n_active`／`index`／`tau_max`／`L_solid`／`L_free` 同此；`generate_layout.SECTION_OUT`／`ARCH_OUT` 在 Task 5 定義、Task 6 引用一致；`axial_stack` 回傳 `{"segments": [(name, t)], "total"}` 在 Task 4 與 Task 6 用法一致。
