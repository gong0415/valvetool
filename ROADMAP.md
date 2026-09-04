# ROADMAP — 第一版 v1 工作項目與模組總表

> 本檔案是工作項目編號的**唯一定義處**。對話、文件一律用此處的 ID 稱呼工作項目。

第一版 v1 由 P0–P7 八個工作項目構成,規格 §3 的 L1–L6 物理極限至此全數到位。P8 為第一版之後的加做項目:B-H 飽和曲線 + 實體佈局。

## 命名規則

- **`P<n>`**:工作項目編號（P0–P8）
- **文件檔名**:設計/計畫文件檔名包含 ID
- **Commit 訊息**:依規格書 §6.0.5 用物理域/極限標籤（`feat(fluid):`、`feat(L3):`）,ID ↔ 標籤對照見下表

## 模組總表

| ID | 名稱 | commit 標籤 | 主要產物與關鍵結論 |
|---|---|---|---|
| P0 | 電磁+機械耦合核心 + Streamlit 參數介面 | `feat(core)`／`feat(gui)` | `magnetics.py`（磁路:reluctance／inductance／coenergy 力／飽和標記,只吃物理氣隙 gap）;`dynamics.py`（電路–機械耦合 ODE:`electrical_di_dt`／`spring_force`／`net_mechanical_force`／`mechanical_dv_dt`／`coupled_rhs`／`simulate_opening`）;`app.py`（參數輸入介面）;報告 `generate_p0_report.py`（基準案例 `t_open` = 5.027 ms） |
| P1 | 流體域:獨立函式庫 + 動力學耦合 | `feat(fluid)` | `fluid.py`（`effective_area`／`mdot_gas`／`mdot_liquid`／`water_hammer_dp`／`flashing_risk`／`Cv`/`Kv`／`flow_force_gas`／`flow_force_liquid`;輸入域鉗位確保 ODE 求解器安全;不匯入 dynamics/magnetics）;雙模式 `simulate_opening`（乾跑 bit-level 由特徵測試凍結,耦合模式 `delta_P` 被工況取代,避免雙重帳本）;GUI 流體面板;對照報告 `generate_p1_coupled_report.py`（N₂ 基準 5.027→5.031 ms,+4.4 µs）;**已知限制**:射流角 cosθ 未修正 |
| P2 | 熱域 | `feat(thermal)` | `thermal.py`（R(T) 線性 α + 77 K 低溫查表、穩態自熱平衡、居里餘裕檢查）;雙模式 `simulate_opening(T_coil=None)` 乾跑 bit-level 不變;`limits.py` 起檔（L6 動態拉入電壓二分搜尋 `find_pullin_voltage`）;L6 包絡報告 `generate_p2_l6_envelope.py`（+70°C 熱浸+自熱 V_pull-in = 21.49 V,距 22 V 母線下限僅 0.51 V;LN₂ 77 K 附 B_sat 警示）;解鎖 L6 與響應特性 4 |
| P3 | 關閉瞬態 + 續流電路 | `feat(dynamics)` | `flyback.py`（二極體／稽納／RC 三種電路）;`dynamics.simulate_closing()`（兩階段積分:止擋保持期 + 釋放衝程期）;`electrical_di_dt` 泛化支援外部 `V`,乾跑開啟相 bit-level 不變;比較報告 `generate_p3_flyback_report.py`（基準參數組二極體→稽納關閉時間比 ≈3.59 倍,落在規格宣稱的 3–5 倍區間）;涵蓋響應特性 2 |
| P4 | 繞線窗口 + L2 Pareto + L1 響應時間下限 | `feat(L2)`／`feat(L1)` | `winding.py`（`R = ρ·N²/(A_w·k_fill)·l_turn`,補上本體論原始項 1）;`ValveParams` 新增 `l_turn_mean`／`k_fill`／`A_winding`（校準至 R ≈ 80 Ω,不動 dynamics）;`limits.l2_pareto_front`（窗口斜率 + B_sat 天花板;基準深度飽和 → spike-and-hold 論點）;`limits.py` L1 組（`electrical_time_constant`／`motion_threshold_current`／`response_time_bound`／`naive_response_time_estimate`／`l1_sweep`）;**τ_e 對匝數不變**（L2 的 N²/R 不變之直接推論:τ_e = L/R = (N²/R)/ℛ）;真下界在 18 個可吸合案例 100% 成立、鬆緊度 0.796–0.913;**如實標示規格的 τ_e + √(2mx/F) 分解式不是界**（實測 0.43–3.11 倍,兩側皆會跑）;報告 `generate_p4_l2_pareto_report.py`／`generate_p4_l1_report.py`;**仍缺**:渦流延遲修正 |
| P5 | L3 小型化 + Buckingham π | `feat(L3)` | `sizing.py`（質量模型 + `scale_params` 自相似縮放 + 兩個正交封閉解 → 矩形可行域）;`ValveParams` 新增 `k_pack`／`J_max`（皆 EMPIRICAL,ODE 不使用）;`limits.dimensionless_groups`（五個 π 群 + 縮放律）;**發現規格點名的三個 L3 機制尺度不變、單獨給不出極限**,極限改由電流密度上限 `J_max` 提供;π₂ ∝ D 預測小閥由死時間主導轉為機械主導（基準 23.25,s=0.2 時降至 4.65）;報告 `generate_p5_l3_report.py` |
| P6 | 密封域/L4、閥座回彈/L5、材料域 | `feat(L4)`／`feat(L5)` | `sealing.py`（接觸應力→實接觸面積分數→滲流判別→殘留間隙→Knudsen 分流→氦漏率;`leak_rate`／`land_width_for_seal`）;`impact.py`（Hertz 彈性撞擊閉式解 `indentation_max`／`contact_pressure_max`／`contact_duration`／`is_shakedown` + 恢復係數回彈序列 `bounce_sequence`,被 `dynamics`／`limits` 共用）;`materials.py`（PCTFE/Vespel/PTFE/17-4PH/440C/430F/純鐵/銅材料表,`source` 區分 `ASTM_E595`／`LITERATURE`;`seat_pair()` 組出 `SeatPair`）;`seat=None` 雙模式（第四個正交旋鈕,`None` 時與既有呼叫者 bit-level 相同）;`limits.l4_leak_curve`／`l5_life_curve`／`bounce_leak_spike`;`ValveParams` 新增 `w_land`／`R_tip`（皆 EMPIRICAL）;**發現 1**:密封在此力量級近乎二元——PCTFE/440C 軟閥座 land ≤66.082 µm 即密封,17-4PH/440C 金屬對金屬需 ≤1.652 µm（近乎刀刃線接觸）,門檻以下漏率高於規格 1e-4 scc/s 達 2–4 個數量級;**發現 3**:稽納加速的代價——關閉時間降至 0.28 倍的同時,觸座速度升至 3.90 倍、撞擊動能升至 15.20 倍;**發現 6**:開啟回彈的撞擊速度比不等於恢復係數 e,而是從上方收斂（PCTFE e=0.4 實測 0.605→0.430→0.409→0.403→0.401）,因銜鐵到止擋時磁力餘裕僅約 0.9 N（~28.9 N vs ~28 N 阻力）,首次長飛行內電流仍在爬升;報告 `generate_p6_report.py` + 三張圖;**已知限制**:洩漏絕對值不可信（`h→Rq_c` 退化）,可信的是密封/不密封的門檻位置 |
| P7 | GUI 整合 P3–P6 | `feat(gui)` | `app.py` 拆為 `compute_*` 純函式（`compute_result`／`compute_opening_bounce`／`compute_closing`／`compute_bounce_spike`）+ `render_*_tab`（`render_opening_tab`／`render_closing_tab`／`render_limits_tab`／`render_sealing_tab`）+ `main()`;主畫面 `st.tabs` 四頁（開啟動態｜關閉/續流｜極限掃描 L1–L3｜密封/壽命 L4–L5）;側欄 10 個 expander,含繞線窗口推導電阻與 `R_coil_20C` 的一致性檢查;`generate_schematic.py` 新增 `generate_winding_window()`／`generate_seal_land()`,連同既有 `valve_schematic.png` 共三張示意圖;ODE 密集運算（關閉模擬、L1 掃描、拉入電壓二分搜尋、回彈漏量尖峰）改為按鈕觸發並存進 `st.session_state`（Streamlit 每次互動會重跑整支 script）,閉式解運算（L2 Pareto、L3 可行域、L4 漏率、L5 壽命、Hertz 接觸量）維持即時重算 |
| P8 | B-H 飽和曲線 + 實體佈局 | `feat(magnetics)`／`feat(layout)` | `magnetization.py`（Fröhlich-Kennelly 解析式 + 查表，零新增自由度）；`magnetics.solve_flux_density`（對 B 二分解 MMF 平衡）+ `mag=None` 第五旋鈕（線性路徑 bit-level 不變）；`drive.py` 全面支援 `mag`；`layout.py`（磁通連續定磁軛截面、壁厚取磁性/結構/製造三下限之最大、彈簧封閉解、銜鐵厚度由 m_arm 反解 + 一致性檢查、LAYOUT_EMPIRICAL 集中宣告）；`generate_layout.py` 兩張圖（按比例剖面 + 架構）；GUI 第五頁籤。**發現 1**：線性模型「吸力隨氣隙關閉暴增」是假象——真實吸力沿行程僅 +24%（26.5→32.7 N），線性報 12 倍（44.6→530 N），影響 P6 撞擊速度與 L5 壽命估算；**發現 2**：`N2_25BAR_PH` 宣稱的 2.0× 保持裕度實為 **1.60×**，真 2.0× 需 13.07 mA（非 11.44 mA），保持功耗 0.037→0.048 W；**發現 3**：外殼壁厚的磁性下限 0.35 mm 與結構下限 0.33–0.44 mm 皆低於可加工厚度，壁厚實由製造決定；**刻意不做**：磁滯迴線（殘留吸力 ~0.02 N vs 彈簧預載 2.2 N，差兩個數量級） |

## 待決事項

| 項目 | 選項 |
|---|---|
| `damping_coeff` 處置 | (a) 保留並持續標註為簡化項 (b) 以擠膜阻尼（squeeze-film damping）第一性模型取代 |

## 相關文件

- 規格書:`docs/spec/solenoid_valve_first_principles.md`
- 參數手冊:`docs/spec/PARAMS.md`
- 資料夾結構:`ARCHITECTURE.md`
- 使用說明:`USAGE.md`
- 已知模型簡化項:`solenoid_model/report/MODEL_NOTES.md`
