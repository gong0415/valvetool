# 實體電磁閥架構佈局設計（P8）

> 設計文件。基準案例：`N₂ 25bar 峰值保持`（`cases.N2_25BAR_PH`）。
> 本文所有數值皆由模型或本文所列封閉解算出並驗證，未經驗證者一律標 EMPIRICAL。

## 1. 目標與缺口

現有 P0–P7 建立的是**磁路等效參數**（`A_gap`、`g0`、`l_core`…），不是可製造的實體。
`docs/spec/valve_schematic.png` 是不按比例的參數對照圖，非實體佈局。

本工作項目補上從 `ValveParams` 到實體尺寸鏈的推導層：外殼與磁軛截面、極面幾何、
彈簧實體尺寸、軸向尺寸鏈、外型包絡，並產出按比例剖面圖與架構圖。

**非目標**：3D CAD 模型、公差配合表、製程與加工工藝、材料採購規格。本層只到
「主要尺寸與其物理依據」為止。

## 2. 設計基準（來自 cases.py，已驗證）

| 項目 | 數值 | 來源 |
|---|---|---|
| 極面面積 A_gap | 20.0 mm² | `N2_25BAR_PH_PARAMS` |
| 極面直徑 D_pole | 5.05 mm | √(4·A_gap/π) |
| 靜止氣隙 g0 | 0.20 mm | 同上 |
| 行程 x_stroke | 0.15 mm | 曲率面積飽和 x=D/4 |
| 閥座孔徑 | 0.60 mm | 1.28 g/s 壅塞流反推 |
| 線圈 OD/ID/高 | 18 / 6 / 12 mm | `_PH_OD/_PH_ID/_PH_H` |
| 匝數／線徑 | 4000 匝／107 µm | `winding.wire_diameter` |
| 線圈電阻 | 281.5 Ω | `winding.coil_resistance` |
| 峰值／保持電流 | 99.5 / 11.4 mA | `drive.hold_current` |
| 彈簧剛度／預載 | 4000 N/m / 2.2 N | 同上 |
| 活動件質量 | 0.8 g | `m_arm` |
| B_sat（軟鐵） | 2.15 T | 同上 |

## 3. 磁路截面：為何不用線性模型的 B

線性磁路模型在峰值 99.5 mA、閉合氣隙下算出 **B = 8.16 T**，物理上不可能
（`cases.py` 註解已記錄此限制）。直接拿它定磁軛截面會得到 75.9 mm² 的荒謬結果。

**採用的規則**：磁通連續 + B 被 B_sat 鉗住。

極面實際 B 最高只能到 B_sat。磁軛與極面同材質時，磁軛只要不比極面先飽和即可：

    Φ = B_pole · A_gap = B_yoke · A_yoke ,  B_pole ≤ B_sat
    ⇒ A_yoke ≥ A_gap · (B_pole/B_sat) ,  最嚴苛在 B_pole = B_sat
    ⇒ A_yoke ≥ A_gap = 20.0 mm²

這是保守的封閉解，且與 `sizing.pressure_ceiling` 同樣的「B_sat 鉗位」思路一致。
異材質（例如外殼用不鏽鋼）時 B_sat 不同，公式保留 B_sat 比值項。

**推導出的尺寸**（已驗證）：

| 尺寸 | 值 | 依據 |
|---|---|---|
| 鐵芯半徑 | 2.52 mm | √(A_gap/π) |
| 磁軛截面 A_yoke | ≥ 20.0 mm² | 磁通連續 |
| 外殼磁性壁厚 | 0.35 mm | 環形回路 π(Ro²−Ri²)=A_yoke，Ri=9.0 mm |
| 端板厚度 | 1.26 mm | 徑向磁通 A=2πr·t，最嚴苛在 r=鐵芯 OD |
| 外型外徑（磁性下限時） | 18.7 mm | 線圈 OD 18 + 2×0.35；**採用壁厚後為 20.0 mm，見 §4** |

端板厚度取最小半徑處（鐵芯外緣）是因為徑向磁通的可用截面隨 r 增加，
r 最小處最容易飽和。

## 4. 外殼壁厚：磁性下限 vs 結構下限

外殼壁厚有兩個獨立下限，**取大者**：

**磁性下限** 0.35 mm（上節）。

**結構下限**：薄壁環向應力 t = P·Ri/(Sy/SF)，P=25 bar、Ri=9.0 mm：

| 材料 | Sy | SF=3 | SF=4 |
|---|---|---|---|
| 304 SS | 205 MPa | 0.33 mm | 0.44 mm |
| 17-4PH H900 | 1000 MPa | 0.068 mm | 0.090 mm |

結構計算假設外殼即承壓邊界。實際上磁性外殼未必承壓（見 §5 承壓邊界討論）。

**結論**：兩者同量級（0.33–0.44 vs 0.35 mm），皆遠低於可製造下限。
實務壁厚由**加工與剛性**而非物理決定 → 取 **1.0 mm（EMPIRICAL）**，
標註方式比照現有 `k_pack`／`C_d`。模組同時回報三個值（磁性、結構、採用）
與採用理由，不隱藏判斷。

採用 1.0 mm 後，**外型外徑 = 18 + 2×1.0 = 20.0 mm**（§3 的 18.7 mm 是磁性
下限情形，非最終值）。壁厚加厚只增裕度、不損磁路：A_yoke 隨之大於下限。

## 5. 承壓邊界（架構決策）

常閉直動式、P_up=25 bar、P_down=真空。承壓邊界必須把 25 bar 與外界隔開，
且線圈在真空側或大氣側會影響散熱路徑（`G_th_cond`／`A_rad` 已在 `thermal.py`）。

**採用：濕式銜鐵（wet armature）**——流體充滿銜鐵腔，線圈在承壓邊界外，
以非磁性隔離套（isolation sleeve）分隔。理由：

- 動密封只在靜止的閥座處，銜鐵不需通過動密封 → 消除滑動密封的洩漏與摩擦
- 線圈可維修／更換，不接觸 N₂
- 代價：隔離套厚度直接加進磁路氣隙

**隔離套對磁路的影響是本設計最關鍵的耦合項**：套壁厚 t_sleeve 落在徑向磁路上，
但**不落在軸向工作氣隙 g0 上**（銜鐵與極面之間仍是 g0）。徑向路徑的等效氣隙
會提高總磁阻、降低吸力。此項在 `magnetics.reluctance` 中未建模
（現有模型只有軸向 gap + 鐵芯路徑）→ **列為已知限制**，於 MODEL_NOTES.md 記錄，
並在模組回報中明示「隔離套磁阻未計入吸力裕度」。

隔離套厚度取 **0.25 mm（EMPIRICAL）**：非磁性 316L，25 bar 下結構足夠
（同 §4 公式，Ri=3.0 mm 時結構下限僅 0.11 mm @ SF=4），其餘為加工餘裕。

## 6. 彈簧實體尺寸（封閉解）

由 `k_spring=4000 N/m`、`F_preload=2.2 N` 反解幾何。
螺旋壓縮彈簧 k = G·d⁴/(8·D³·n)，G=79 GPa（302 不鏽鋼）。

工作條件：預載變形 0.550 mm，全行程變形 0.700 mm，F_max=2.80 N。

含 Wahl 修正的候選解（已驗證）：

| d (mm) | D (mm) | n | C=D/d | τ_max (MPa) | 實高 (mm) | 自由長 (mm) |
|---|---|---|---|---|---|---|
| 0.30 | 1.50 | 5.92 | 5.00 | 519 | 2.38 | 3.18 |
| **0.35** | **1.80** | **6.35** | **5.14** | **389** | **2.92** | **3.72** |
| 0.40 | 2.00 | 7.90 | 5.00 | 292 | 3.96 | 4.76 |

**採用 d=0.35 / D=1.80 / n=6.35**：τ=389 MPa 對 302 不鏽鋼許用
700–900 MPa 有 1.8–2.3× 裕度；彈簧指數 C=5.14 落在易製造的 4–12 範圍；
自由長 3.72 mm 符合軸向預算（§7）。

註：k=4000 N/m 在此力量級偏硬，導致細線無解（d=0.15 mm 只得 n=0.37 匝，
非實體）。若後續放寬 k，彈簧可再細。

## 7. 軸向尺寸鏈

自閥座面往上累加（EMPIRICAL 項已標示）：

| 段 | 厚度 (mm) | 依據 |
|---|---|---|
| 閥座座體 | 3.0 | EMPIRICAL，容納 0.60 mm 孔與密封 land |
| 閥座↔銜鐵行程 | 0.15 | x_stroke |
| 銜鐵 | 5.08 | **推導值**：由 m_arm=0.8 g 與軟鐵密度反算 |
| 工作氣隙 g0 | 0.20 | 模型 |
| 固定極 | 4.0 | EMPIRICAL |
| 彈簧腔 | 3.72 | 彈簧自由長 |
| 端板 ×2 | 2×1.26 | §3 磁通連續 |
| **合計** | **18.67** | 上列之和；未含端蓋與接管 |

銜鐵厚度**不是 EMPIRICAL 值，而是由 `m_arm` 反解**：以極面直徑 5.05 mm 的
圓柱、軟鐵密度 7870 kg/m³ 計，t = m_arm/(ρ·A_pole) = **5.08 mm**。
（起草時誤填 4.0 mm，該值只給 0.631 g，與宣告的 0.8 g 不符——正是下述
一致性檢查要攔截的情形。）

模組提供 `armature_mass_consistency()` 反向檢查，比照
`winding.coil_resistance_consistency` 的「回報而不拋錯」做法：銜鐵若因結構
需要改為非等徑（例如加導向段），厚度與質量會脫鉤，此時由檢查函式回報偏差，
由設計者決定是否回頭修正 `m_arm`。

## 8. 模組設計：`solenoid_model/layout.py`

純函式，**不 import `dynamics`／`magnetics`／`fluid`**（沿用 `sizing.py`、
`winding.py` 的既有不變式），輸入 `ValveParams` 輸出尺寸鏈。

```
layout.py
├─ 磁路推導（第一性）
│   yoke_area(params)                 -> A_yoke [m²]
│   shell_wall_magnetic(params, R_i)  -> t [m]
│   end_plate_thickness(params)       -> t [m]
│   core_radius(params)               -> r [m]
├─ 結構推導（第一性）
│   hoop_wall_thickness(P, R_i, Sy, SF) -> t [m]
│   shell_wall_thickness(params, ...) -> {magnetic, structural, adopted, reason}
├─ 彈簧推導（第一性）
│   spring_geometry(params, d, D, G)  -> {n, C, tau_wahl, L_solid, L_free}
│   spring_stress_ok(...)             -> bool
├─ 一致性檢查（比照 winding.coil_resistance_consistency）
│   armature_mass_consistency(params, t_arm, rho) -> (m_declared, m_geom, rel_err)
└─ 組裝
    axial_stack(params, **empirical)  -> 有序尺寸鏈 + 總長
    envelope(params, ...)             -> {OD, L, 各段}
```

**EMPIRICAL 參數處置**：不寫死在函式內，改為帶預設值的具名引數，
並由一個 `LAYOUT_EMPIRICAL` dict 集中宣告（值 + 單位 + 理由 + 來源標籤），
使「哪些是推導、哪些是判斷」在程式碼層面可見。

## 9. 圖：`report/generate_layout.py`

沿用 `generate_schematic.py` 的 `_CJK_FONT_RC` 中文字型處理與配色常數。

1. **按比例剖面圖** `docs/spec/layout_section.png`
   真實 mm 比例（`set_aspect('equal')`）、軸對稱半剖、尺寸鏈標註、
   材料分區著色（磁性件／非磁性件／密封件／線圈）。

2. **架構圖** `docs/spec/layout_architecture.png`
   四個功能域方塊與介面：電氣（驅動→線圈）、磁（線圈→鐵芯→氣隙→銜鐵→
   磁軛→回路）、流體（入口→閥座→出口）、機械（彈簧／止擋／銜鐵）。
   標出跨域耦合點（氣隙＝磁↔機械、閥座＝機械↔流體）。

## 10. GUI：第五頁籤

`app.py` 新增 `render_layout_tab()`，比照既有 `render_*_tab` 模式：
尺寸鏈表格（含推導／EMPIRICAL 標記）、兩張圖、一致性檢查結果、
以及 §5 的隔離套限制警示。皆為閉式解 → 即時重算，不需按鈕觸發。

## 11. 測試：`tests/test_layout.py`

- 磁通連續不變式：`yoke_area(p) >= p.A_gap`；異材質時比值正確
- 壁厚取大者邏輯：三種情境（磁性大／結構大／EMPIRICAL 大）各一案
- 彈簧封閉解：由 (d, D, n) 正算回 k，須與 `params.k_spring` 相符（往返一致）
- Wahl 修正在 C→∞ 時趨近 1
- 自相似縮放：`sizing.scale_params(p, s)` 下長度尺寸按 s、截面按 s² 縮放
- 銜鐵質量一致性檢查能抓出不一致
- EMPIRICAL 宣告完整性：每項都有單位、理由、來源標籤
- 產圖腳本可執行且輸出檔案存在（比照 `test_report_script.py`）

## 12. 已知限制（同步寫入 MODEL_NOTES.md）

1. **隔離套磁阻未計入**：`magnetics.reluctance` 只有軸向氣隙 + 鐵芯路徑，
   濕式銜鐵的徑向隔離套等效氣隙未建模 → 實際吸力低於模型報告值。
2. **峰值相 B 超過 B_sat**：沿用 `cases.py` 既有記錄；本層以 B_sat 鉗位
   規避，但不修正吸力預測本身。
3. **軸向尺寸鏈部分為 EMPIRICAL**：閥座座體與固定極厚度由封裝經驗定，非推導值，
   且無一致性檢查可攔截。銜鐵厚度由 `m_arm` 反解且有檢查（§7）。
4. **未含公差與配合**：所有尺寸為標稱值。
