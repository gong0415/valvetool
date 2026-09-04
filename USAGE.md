# USAGE — 使用說明書

> 電磁閥開發工具的使用指南。新增功能後同步更新本文件（維護方式同 `ARCHITECTURE.md`）。
> 涵蓋：P0 電磁-機械模擬與網頁介面、P1 流體函式庫與耦合、P2 熱域耦合、
> P3 關閉瞬態/flyback、P4 繞線窗口+L1/L2 掃描、P5 L3 小型化+Buckingham π、
> P6 密封域/L4+閥座回彈/L5+材料域、P7 GUI 四頁籤整合、
> P8 B-H 飽和曲線+實體佈局（GUI 第五頁籤）。

---

## 方式一：網頁介面（最簡單，不用寫程式）

```bash
cd /Users/eric/Desktop/solenoid-valve
streamlit run solenoid_model/app.py
```

瀏覽器開啟 http://localhost:8501 。

介面分成**側欄（參數輸入）**與**主畫面五頁籤（結果/掃描）**：改側欄任何數字，五個頁籤在同一次 rerun 內全部重算並重新渲染（Streamlit 的執行模型：每次互動重跑整支 script）。

### 側欄：10 個折疊區塊（由上到下）

1. **參數示意圖**（預設收合）——閥件剖面圖 `valve_schematic.png`，對照 `docs/spec/PARAMS.md`
2. **電氣**（預設展開）——`V_bus`／`R_coil_20C`／`N_turns`
3. **磁路**（預設展開）——`A_gap`／`g0`／`x_stroke`／`l_core`／`mu_r_core`／`B_sat`
4. **機械**（預設展開）——`m_arm`／`k_spring`／`F_preload`
5. **簡化負載**（預設展開）——`delta_P`／`A_seat`／`damping_coeff`（⚠️ 阻尼非第一性推導）
6. **流道幾何**（預設展開）——`D_seat_bore`／`C_d`（⚠️ `C_d` 為經驗係數，不要與 Cv 混淆）
7. **流體耦合**（預設展開）——「啟用流體耦合」勾選（預設開）+ 流體選擇（N₂/Xe/He/水/LN₂/自訂氣體/自訂液體）+ `P_up`／`P_down`／`T0`
8. **熱域**（預設展開）——線圈溫度 `T_coil`（K，預設 293.15 = 20°C；等溫假設，單次開啟遠短於熱時間常數）
9. **繞線窗口**（預設收合）——`winding_window_schematic.png` 示意圖 + `A_winding`／`k_fill`／`l_turn_mean`；即時顯示由窗口反推的電阻 `R_derived`，與上面 `R_coil_20C` 差超過 5% 會出現黃色不一致警告
10. **閥座/密封**（預設收合）——`seal_land_schematic.png` 示意圖 + 閥座材料對選單（`PCTFE/440C`／`17-4PH/440C`／自訂材料+恢復係數 `e`）+ `w_land`／`R_tip`（皆 ⚠️ EMPIRICAL）

側欄湊出的 `ValveParams` 會同時餵給下面五個頁籤；各頁籤若吸不動/算不出（如彈簧預載過大），只有該頁籤顯示紅色錯誤訊息，其餘頁籤仍正常渲染（各頁籤各自包一層 `try/except ValueError`）。

### 主畫面：五個頁籤

**頁籤 1「開啟動態」**——原始開啟模擬，每次側欄改動即時重算（單次 ODE 積分到吸合為止，成本低不需按鈕）：
- **t_open** 開啟時間；流體耦合開啟時旁邊顯示與乾跑（不耦合）的差值
- 流體耦合開啟時：全開 Cv、全開質量流率；選液體且下游壓力 ≤ 蒸氣壓時顯示黃色**閃蒸警告**
- 熱態線圈電阻（相對 20°C 差值）+ 持續通電保持的穩態自熱平衡溫度（真空傳導+輻射，定電流保守假設；無穩態解則顯示熱失控錯誤）
- 按鈕「計算拉入電壓（數秒）」——二分搜尋 `V_pull-in`（`limits.find_pullin_voltage`，內部要跑多次 ODE 才收斂，故按鈕觸發，不隨側欄變動自動重跑）；算出後顯示對供電電壓的裕度，若高於供電電壓則警告超出工作包絡
- 磁通密度超過 `B_sat` 時的黃色飽和警告
- 勾選框「啟用閥座回彈」（預設關）——開啟後銜鐵到止擋改以恢復係數彈回（側欄閥座材料對），顯示回彈次數、首次/末次觸擊速度、`t_settle`；關閉＝既有硬鉗制行為，頭條 t_open 不受影響
- 電流 i(t)／位移 x(t)／磁力 F_mag(t) 三張圖

**頁籤 2「關閉/續流」**——斷電後的關閉瞬態：
- 續流電路選擇：二極體（`V_fwd`）／稽納（`V_clamp`）／RC 緩衝（`R_snub`／`C_snub`，刻意無預設值，輸入非正值會擋下並提示）
- 勾選框「與二極體對照」——加跑一次二極體續流基準，顯示關閉時間比與觸座速度比
- 勾選框「閥座回彈」——關閉觸座後以恢復係數回彈積分
- 按鈕「模擬關閉」——關閉模擬（止擋保持期+釋放衝程期兩段積分）按鈕觸發，結果存進 `st.session_state["closing"]`；未按過時只顯示提示訊息，其餘頁籤（如頁籤 4 的 L5）會沿用上一次按下的結果作為預設撞擊速度
- 按下後：`t_release`／`t_close`／觸座速度三個 metric；勾了對照時顯示與二極體的倍率比較
- 電流衰減圖（保持期+衝程期）、位移圖（回彈時紅色標出回彈段）

**頁籤 3「極限掃描 L1–L3」**——三組極限，皆為閉式解或快速掃描（除 L1 掃描外都不需按鈕）：
- **L1**：`τ_e`／`i_th`／L1 響應時間下限即時算出（若穩態電流吸不動閥則顯示無界錯誤）；規格分解式 naive 估計旁附警語（**不是**上下界，實測落在真值 0.43–3.11 倍之間）；掃描維度選單（電壓／匝數-窗口一致 R）+ 掃描點輸入 + 按鈕「執行 L1 掃描（每點跑一次 ODE）」——每個掃描點都要積分一次 ODE，成本隨點數線性增加，故按鈕觸發，結果存 `st.session_state["l1_sweep"]`，側欄參數變了但沒重按會提示「顯示的結果對應舊參數」
- **L2**：`N_lo`／`N_hi`／點數三個輸入即時重算保持力-功率 Pareto 前緣（純解析公式，非 ODE，不需按鈕）；`B_sat` 保持力天花板 metric；Pareto 前緣圖
- **L3**：Buckingham π 五群表；`k_pack`／`J_max` 輸入（⚠️ 皆 EMPIRICAL 工程係數，ODE 不使用）；最小可行縮放 `s*`／壓力天花板／整閥質量估計三個 metric；縮放-壓力可行域熱力圖（皆為解析尺寸公式掃描，不需按鈕）

**頁籤 4「密封/壽命 L4–L5」**——密封與撞擊壽命，同樣以閉式解為主：
- **L4**：密封 land 寬度上限、接觸應力、目前 `w_land` 密封判定（不密封時附漏率）三個 metric；警語：模型可信的是門檻位置，不是漏率絕對值；land 寬度 vs 漏率曲線（兩內建閥座對照 + 規格線）
- 材料/環境準則：輻射劑量、推進劑選單；各閥座材料的出氣（TML/CVCM）、輻射裕度、推進劑相容性表；冷焊風險評語
- **L5**：評估撞擊速度輸入（預設取頁籤 2 最近一次關閉模擬的觸座速度，無則 0.39 m/s）；Hertz 峰值接觸應力／接觸時間／Shakedown 三個 metric；撞擊速度 vs 壽命曲線（shakedown 區間標色，疊上頁籤 2 的關閉工況點）；警語：壽命係數皆 EMPIRICAL，只能當量級參考
- 「關閉回彈漏量尖峰」（L4×L5 耦合）：僅氣體介質可算（液體會顯示提示改選側欄氣體介質）；續流電路選單（稽納/二極體，預設稽納）+ 按鈕「計算回彈漏量尖峰（跑關閉 ODE）」——這段要重跑一次帶回彈的關閉 ODE，故按鈕觸發，結果存 `st.session_state["spike"]`；顯示關閉回彈次數、回彈離座總時間、回彈漏出質量、靜態密封漏率四個 metric

**頁籤 5「實體佈局」**（P8）——把 `ValveParams` 展開成一條可製造的軸向尺寸鏈，皆為閉式解，不需按鈕：
- 三個 metric：外徑、總長、極面直徑
- 「軸向尺寸鏈」表：閥座座體／行程／銜鐵／工作氣隙／固定極／彈簧腔／端板 x2，依序疊加
- 「外殼壁厚：三個獨立下限」表：磁性/結構/製造三個下限並列，並顯示採用值與判斷理由（哪一項下限最大）
- 彈簧一行摘要：線徑/中徑/有效圈數/彈簧指數/最大剪應力/自由長（彈簧設計點取自規格 §6 的 `layout.SPRING_DESIGN`，DERIVED）
- 銜鐵厚度與宣告的 `m_arm` 偏差超過 5% 時顯示黃色不一致警告
- 固定提示：濕式銜鐵隔離套厚度落在徑向磁路上但 `magnetics.reluctance` 未建模，實際吸力低於模型報告的裕度
- 「EMPIRICAL 項目」表：`layout.LAYOUT_EMPIRICAL` 五項（`t_seat_body`／`t_fixed_pole`／`t_manufacturing`／`t_sleeve`／`clearance_spring`）連同數值與理由
- 按比例剖面圖（`layout_section.png`）與四域架構圖（`layout_architecture.png`）：兩張圖若已產生則直接顯示；未產生時頁籤留白，見下方指令自行產生

### 為什麼有些東西要按鈕才跑

Streamlit 的執行模型是**每次互動整支 script 全部重跑**——包括所有五個頁籤的內容，不是只重跑被改到的那個頁籤。純解析／查表計算（L2 Pareto、L3 可行域、L4 漏率曲線、L5 壽命曲線、Hertz 接觸量、τ_e/i_th 等）成本可忽略，每次重跑都算沒問題，直接顯示。但凡是要積分 ODE 的計算——開啟模擬本身例外（單次積分到吸合，成本仍低，故也是即時重算）——只要是**多次**積分或**額外**一次較貴的積分（關閉模擬的兩段式積分、L1 掃描逐點跑 ODE、拉入電壓二分搜尋要反覆試電壓、回彈漏量尖峰要重跑一次帶回彈的關閉模擬），就會拖慢每一次側欄互動，所以一律做成按鈕，並把結果快取進 `st.session_state`，按下之後即使側欄後續改了別的參數也不會自動重算——因此每個按鈕結果旁都有「顯示的結果對應舊參數，請重按」的提示邏輯，讀者看到這行提示代表要重新按按鈕才能拿到跟當前側欄一致的數字。

### EMPIRICAL 警語怎麼讀

畫面上出現「⚠️ EMPIRICAL」或「⚠️ 非第一性推導」字樣時，代表該數字/係數不是從電磁-機械-流體第一性物理推導出來的，而是工程經驗值或簡化假設（例如 `damping_coeff`、`C_d`、`k_pack`、`J_max`、`w_land`、`R_tip`、恢復係數 `e`、Archard 磨耗係數）。這類警語不代表功能有 bug，而是誠實揭露模型的可信範圍：**結論可信的通常是門檻位置或量級（例如密封 land 的密封/不密封分界、壽命的數量級），不是絕對數值**——引用這些結果時應照著警語旁的措辭，不要把 EMPIRICAL 數字當成量測值或鑑定數字使用。

## 方式二：跑基準案例報告（一行指令）

```bash
python3 -m solenoid_model.report.generate_p0_report
```

印出基準閥件的 `t_open = 5.027 ms`，並在 `solenoid_model/report/` 更新三張 PNG（電流、位移、磁力曲線）。

乾跑 vs 流體耦合對照報告（P1）：

```bash
python3 -m solenoid_model.report.generate_p1_coupled_report
```

印出乾跑 5.027 ms / N₂ 耦合 5.031 ms（+4.4 µs）、全開 Cv/質量流率/流體力，並更新 `dry_vs_coupled_x_t.png` 位移對照圖。

L6 電壓-溫度工作邊界報告（P2，全解析度需**數分鐘**）：

```bash
python3 -m solenoid_model.report.generate_p2_l6_envelope
```

印出包絡表（−40°C 到 +70°C 掃描，含自熱前後兩條 `V_pull-in` 曲線）、LN₂ 77 K 域外案例（含 B_sat 飽和警告）、響應特性 4 表（`t_open` 隨溫度/電壓漂移），並更新 `l6_envelope_v_t.png` 包絡圖。基準案例：+70°C 熱浸 + 持續自熱下 `V_pull-in ≈ 21.5 V`，距 22 V 匯流排下限 ~0.5 V 裕度。

基準閥件規格：28 V、MEOP 2.4 MPa、D_seat_bore 0.14 mm（全開 Cv = 0.0007，小推力等級姿控推力器閥）、t_open < 10 ms。完整參數值見 `solenoid_model/baseline.py`。

實體佈局兩張圖（按比例剖面 + 四域架構，P8，用 `cases.py` 的 `N2_25BAR_PH_PARAMS`）：

```bash
.venv/bin/python -m solenoid_model.report.generate_layout
```

更新 `docs/spec/layout_section.png`（按真實 mm 比例的軸對稱半剖圖）與 `docs/spec/layout_architecture.png`（電氣/磁/機械/流體四域方塊圖）；GUI 頁籤 5「實體佈局」若找不到這兩個檔案會留白，需先跑一次此指令。

## 方式三：寫 Python（功能最完整）

### A. 電磁-機械模擬（算開啟時間與響應曲線）

```python
from dataclasses import replace
from solenoid_model import dynamics
from solenoid_model.baseline import BASELINE_PARAMS

# 以基準閥件為底，改你要的參數（欄位定義見 docs/spec/PARAMS.md）
custom = replace(BASELINE_PARAMS, N_turns=3000.0, k_spring=15000.0)

sol = dynamics.simulate_opening(custom)
if len(sol.t_events[0]):
    print(f"t_open = {sol.t_events[0][0]*1e3:.3f} ms")
else:
    print("此參數組合下閥件不會吸合")

# sol.t = 時間軸；sol.y = [電流 i, 位移 x, 速度 v] 三條軌跡，可自行畫圖
```

### B. 流體計算（P1）

```python
from solenoid_model import fluid
from solenoid_model.baseline import BASELINE_PARAMS as P

# 步驟 1：選流體（內建預設常數）
#   氣體：fluid.N2、fluid.XE、fluid.HE
#   液體：fluid.WATER_20C、fluid.LN2_77K
# 步驟 2：給工況（SI 單位：Pa、K）
cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)

# 步驟 3：呼叫計算函式（x = 閥開度，0 到 x_stroke）
fluid.effective_area(P.x_stroke, P)                    # 等效流通面積 m²
fluid.cv_from_geometry(P.x_stroke, P)                  # 該開度的 Cv（全開=0.0007）
fluid.kv_from_geometry(P.x_stroke, P)                  # Kv（公制流量係數）
fluid.mdot_gas(P.x_stroke, P, fluid.N2, cond)          # 氣體質量流率 kg/s（自動判斷壅塞/亞音速）
fluid.critical_pressure_ratio(fluid.N2)                # 壅塞臨界壓力比（N₂=0.528）
fluid.mdot_liquid(P.x_stroke, P, fluid.WATER_20C, cond)  # 液體質量流率 kg/s
fluid.water_hammer_dp(fluid.LN2_77K, delta_v=2.0)      # 急關水錘壓升 Pa（Joukowsky）
fluid.flashing_risk(fluid.LN2_77K, cond)               # 閃蒸風險旗標 True/False
fluid.flow_force_liquid(P.x_stroke, P, cond)           # 液體流體力 N（方向：傾向關閉）
fluid.flow_force_gas(P.x_stroke, P, fluid.N2, cond)    # 氣體流體力 N（壅塞/亞音速自動分支）
```

**流體耦合模擬（P1）**——把流體力接進銜鐵運動方程式：

```python
sol = dynamics.simulate_opening(params, medium=fluid.N2,
                                cond=fluid.FlowConditions(2.4e6, 0.0, 293.0))
```

- `medium` 可為任何 `GasProperties` / `LiquidProperties`；傳了 `medium` 就必須明確傳 `cond`（否則 `ValueError`）
- **耦合模式下 `params.delta_P` 被工況取代**：靜態壓力力 = `(P_up−P_down)·A_seat`，避免雙重帳本
- 不傳 `medium` 即為乾跑模式（向後相容）

**自訂流體**不用改程式碼，直接建物性物件：

```python
hydrazine = fluid.LiquidProperties(rho=1004.0, c_sound=2100.0, P_vap=1900.0)  # 自行查證物性來源
fluid.mdot_liquid(P.x_stroke, P, hydrazine, cond)
```

**LN₂ 等低溫液體注意**：務必檢查 `flashing_risk`——下游壓力貼近蒸氣壓會閃蒸（兩相流），本模型的液體公式即失效。旗標未觸發也不保證完全無空蝕（縮流頸部壓力更低），精細判別尚未實作。

### C. 自訂閥件的完整流程建議

1. 從 `replace(BASELINE_PARAMS, ...)` 改參數（參數意義查 `docs/spec/PARAMS.md`）
2. 跑 `dynamics.simulate_opening` 確認 t_open 與吸合行為
3. 跑 `fluid.cv_from_geometry` / `mdot_gas` 確認流量能力
4. 低溫/液體應用加查 `flashing_risk` 與 `water_hammer_dp`

### D. 熱域耦合（P2）

```python
from solenoid_model import dynamics, limits, thermal
from solenoid_model.baseline import BASELINE_PARAMS

# 熱態開啟（+70°C 熱浸）
sol = dynamics.simulate_opening(BASELINE_PARAMS, T_coil=343.15)

# LN₂ 77 K（注意：穩態電流大，B_sat 飽和警告會觸發）
sol_ln2 = dynamics.simulate_opening(BASELINE_PARAMS, T_coil=77.0)

# 熱態拉入電壓（L6）
v_pi = limits.find_pullin_voltage(BASELINE_PARAMS, T_coil=343.15)

# 持續保持的穩態自熱平衡（真空：傳導+輻射）
T_eq = thermal.equilibrium_temp(I_hold=0.35, T_amb=293.15, params=BASELINE_PARAMS)
```

- `T_coil=None`（預設）= 乾跑（線圈電阻視為 `R_coil_20C` 定值，行為與 P0/P1 完全相同）；傳入具體 K 值即啟用熱耦合，`R_coil` 依溫度改變 → 影響電流 → 影響拉力/開啟時間
- `T_coil` 與流體耦合（`medium`/`cond`）互相正交，可同時使用：`simulate_opening(params, medium=fluid.N2, cond=cond, T_coil=343.15)`
- 報告腳本 `python -m solenoid_model.report.generate_p2_l6_envelope` 全解析度需**數分鐘**（見上方方式二）

### E. 密封與壽命（P6）

```python
from solenoid_model import dynamics, fluid, flyback, limits, materials, sealing
from solenoid_model.baseline import BASELINE_PARAMS as P

# L4：目前的 land 寬度密封嗎？多窄才密封？（滲流門檻 phi_c=0.42）
seat = materials.PCTFE_ON_440C   # 軟閥座；金屬對照見 materials.METAL_17_4PH_ON_440C
sealing.is_percolated(P, seat)                # True/False：目前 w_land 是否密封
sealing.land_width_for_seal(P, seat)          # 密封門檻的 land 寬度上限 [m]（PCTFE ≈ 66 um）

# L4 掃描：land 寬度 x 多閥座對照的洩漏率曲線（cond 型別同流體域 B）
cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
curves = limits.l4_leak_curve(
    P, [materials.PCTFE_ON_440C, materials.METAL_17_4PH_ON_440C],
    w_lands=[10e-6, 50e-6, 100e-6], cond=cond)

# L5：關閉時的閥座撞擊 + 回彈（seat 是第四個正交模式旋鈕，見 ARCHITECTURE.md）
closed = dynamics.simulate_closing(P, flyback.ZenerFlyback(), seat=seat)
closed.t_close           # 首次觸座時間（P3 既有定義，不受 seat 影響）
closed.bounce.n_bounce   # 回彈次數
closed.bounce.t_settle   # 真正密封（回彈衰減完畢）的時間

# L5 壽命估計：撞擊速度 vs 循環壽命（EMPIRICAL 磨耗係數，僅供量級參考）
life = limits.l5_life_curve(P, seat, velocities=[0.4, 0.8, 1.5])
```

- `sealing.leak_rate`/`l4_leak_curve` 只涵蓋**閉合狀態**的洩漏（接觸-狹縫問題）；離座後的流量已由 `fluid.py` 涵蓋（孔口問題），兩者在 `limits.bounce_leak_spike` 於分析層組合
- **洩漏分支的絕對值不可信，可信的是密封/不密封的門檻位置**：工程閉合式 `h=Rq_c·(1−φ)` 在 `p≪H` 時退化為 `h→Rq_c`，漏率對載荷不敏感，詳見 `MODEL_NOTES.md`
- `seat=None`（預設）維持既有硬鉗制止擋行為，對既有呼叫者與 P0–P5 bit-level 相同；傳入 `sealing.SeatPair` 才啟用回彈積分，回彈期間 `t_close`/`t_release` 語意不變，新增 `.bounce.n_bounce`/`.bounce.t_settle`
- 報告腳本 `python -m solenoid_model.report.generate_p6_report` 印出 L4 門檻（PCTFE ≈66 um、17-4PH/440C ≈1.65 um）、L5 帳單（稽納撞擊動能約 15 倍於二極體）、開/關兩側回彈摘要，並更新 `p6_l4_leak.png`／`p6_l5_life.png`／`p6_bounce_x_t.png`

---

## 查資料的地方

| 想查什麼 | 看哪裡 |
|---|---|
| 每個參數的意義、單位、常見誤解警告 | `docs/spec/PARAMS.md` + 對照圖 `docs/spec/valve_schematic.png` |
| 三個「面積」的差別（A_gap / A_seat / D_seat_bore） | `docs/spec/PARAMS.md` 的三面積辨析章節 |
| 工作項目編號（P0–P8）與模組總表 | `ROADMAP.md` |
| 資料夾結構 | `ARCHITECTURE.md` |
| 規格書（任務需求來源） | `docs/spec/solenoid_valve_first_principles.md` |
| 已知模型簡化項（誠實原則清單） | `solenoid_model/report/MODEL_NOTES.md` |

## 驗證安裝／改動後一切正常

```bash
.venv/bin/python -m pytest solenoid_model/tests/ -q   # 應顯示 316 passed
```

⚠️ 直接跑裸 `pytest` 可能解析到系統直譯器（缺 scipy/streamlit），測試模組會直接 collection error。務必用上面 `.venv/bin/python -m pytest` 的完整寫法，確保吃到專案 venv。

## 已知限制（依誠實原則揭露）

- `C_d`（流量係數）為經驗係數（預設 0.8，典型 0.6–0.9，精度約 ±15%），非第一性推導
- `damping_coeff`（阻尼）為工程佔位值，非第一性推導
- 磁飽和：`mag=None`（預設）只做警告旗標，不做非線性 B-H 曲線；傳入 `magnetization.AnalyticBH`/`TabulatedBH`（P8）才真的求解飽和磁路，但單值曲線仍不含磁滯（刻意排除，見 `docs/spec/PARAMS.md`「磁化模型參數」）；渦流延遲未建模
- 流體力已耦合進運動方程式（P1）；未含射流角 cosθ 修正（高估關閉力）；不模擬逆流與兩相流動力學
- 熱域（P2）：`G_th_cond`/`emissivity`/`A_rad` 為 EMPIRICAL 量級估計，非量測值；`equilibrium_temp` 僅定電流保守假設（定電壓自穩定情境未建模）；B_sat 隨溫度的退化（居里點附近磁飽和下降）未建模（P8 的 B-H 曲線同樣未含溫度效應）；暫態 T(t)（熱累積、脈衝模式）明確排除於 P2
- 實體佈局（P8）：濕式銜鐵隔離套（0.25 mm）未計入 `magnetics.reluctance`（實際吸力低於模型報告值）；閥座座體/固定極厚度為無一致性檢查的 EMPIRICAL 值；所有尺寸為標稱值，未含公差與配合；`cases.py` 與 P0–P7 各報告的裕度數字皆在 `mag=None` 線性模式下產生，未以飽和模型重算——完整清單見 `solenoid_model/report/MODEL_NOTES.md`「P8 已知限制」
