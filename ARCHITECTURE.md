# 專案資料架構

> 本檔案記錄電磁閥開發工具的資料夾結構。新增模組時同步更新此檔案。

```
solenoid-valve/
├── USAGE.md                  # 使用說明書（網頁介面/報告腳本/Python API，新功能同步更新）
├── ROADMAP.md                # 工作項目編號（P0–P7）與模組總表——唯一定義處
│
├── docs/spec/                 # 任務規格書、參數手冊與示意圖
│   ├── README.md              # 本資料夾內容指引
│   ├── solenoid_valve_first_principles.md        # 任務規格書（含 §2.0 本體論、§6.0 增量指引）
│   ├── PARAMS.md              # ValveParams/fluid.py 參數定義手冊（P1）；對照 valve_schematic.png
│   ├── valve_schematic.png    # 閥件剖面示意圖（由 generate_schematic.py 的 generate() 產生）
│   ├── winding_window_schematic.png  # 繞線窗口示意圖（A_winding/k_fill/l_turn_mean，GUI 側欄顯示）
│   └── seal_land_schematic.png       # 密封 land 放大示意圖（w_land/R_tip/A_seat，GUI 側欄顯示）
│
├── solenoid_model/                              # 核心程式套件
│   ├── __init__.py
│   ├── params.py           # ValveParams dataclass — 所有可調參數的定義
│   ├── baseline.py         # BASELINE_PARAMS — 預先驗證過的基準案例數值
│   ├── magnetics.py        # 磁路：reluctance/inductance/coenergy力/飽和標記（只吃物理氣隙 gap）
│   ├── dynamics.py         # 電路–機械耦合 ODE：electrical_di_dt、spring_force、
│   │                       #   net_mechanical_force、mechanical_dv_dt、coupled_rhs、
│   │                       #   simulate_opening（雙模式：medium=None 乾跑／medium+cond 流體耦合）
│   ├── fluid.py            # 流體域函式庫（P1）：effective_area、mdot_gas/mdot_liquid、
│   │                       #   water_hammer_dp、flashing_risk、Cv/Kv、flow_force_liquid、
│   │                       #   flow_force_gas；輸入域鉗位（ODE 求解器安全）；不匯入 dynamics/magnetics
│   ├── thermal.py          # 熱域函式庫（P2）：R_coil（線性 α + 低溫查表，77–473 K）、
│   │                       #   equilibrium_temp（定電流保守假設，真空傳導+輻射穩態平衡）、
│   │                       #   curie_margin（檢查用，B_sat 溫度退化未建模）
│   ├── sealing.py          # 密封域（P6）：seal_force/contact_stress/contact_area_fraction、
│   │                       #   is_percolated/land_width_for_seal（滲流門檻）、Knudsen 分流氦
│   │                       #   漏率 leak_rate；純函式，不 import dynamics/magnetics/fluid
│   ├── impact.py           # 閥座撞擊接觸力學（P6）：Hertz 彈性撞擊閉式解（indentation_max／
│   │                       #   contact_pressure_max／contact_duration／is_shakedown）+ 恢復
│   │                       #   係數回彈序列 bounce_sequence；被 dynamics.py 與 limits.py 共用
│   ├── materials.py        # 材料與環境域（P6）：材料屬性表（source 標註 ASTM_E595／
│   │                       #   LITERATURE）+ outgassing_ok／cold_weld_risk／radiation_margin／
│   │                       #   propellant_compatibility 查表準則；seat_pair() 組出 SeatPair
│   ├── limits.py           # L1–L6 物理極限（P2 起 L6；P6 新增 L4/L5：l4_leak_curve／
│   │                       #   l5_life_curve／bounce_leak_spike）：find_pullin_voltage
│   │                       #   （二分搜尋拉入電壓，`T_coil` 可選）；L1/L2/L3 掃描屬 P4/P5，
│   │                       #   完整進度見 ROADMAP.md 進度總表
│   ├── app.py              # Streamlit 四頁籤介面：compute_*（純函式：compute_result／
│   │                       #   compute_opening_bounce／compute_closing／compute_bounce_spike）+
│   │                       #   render_*_tab（render_opening_tab／render_closing_tab／
│   │                       #   render_limits_tab／render_sealing_tab）+ main()（側欄組裝
│   │                       #   ValveParams、st.tabs 分派四頁：開啟動態｜關閉/續流｜
│   │                       #   極限掃描 L1–L3｜密封/壽命 L4–L5，P7）
│   │                       #   執行方式：streamlit run solenoid_model/app.py
│   │
│   ├── tests/               # pytest 測試（同時也是分層驗證）
│   │   ├── test_params.py
│   │   ├── test_magnetics.py
│   │   ├── test_dynamics_rl.py
│   │   ├── test_dynamics_mech.py
│   │   ├── test_dynamics_coupled.py
│   │   ├── test_dynamics_characterization.py  # 乾跑基準凍結（規格 §6.0.3，全程綠燈鐵律）
│   │   ├── test_dynamics_fluid_coupled.py     # 耦合模式（氣/液流體力進 ODE、medium 無 cond 拋錯）
│   │   ├── test_fluid.py    # fluid.py 全套（Cv/A_eff、壅塞流、液體流+水錘、鉗位防護、
│   │   │                    #   flow_force_gas、schematic/對照報告腳本）
│   │   ├── test_report_script.py
│   │   ├── test_app.py      # compute_result 純函式測試 + AppTest UI 測試（含流體面板）
│   │   ├── test_thermal.py             # R_coil（線性+查表）、equilibrium_temp（含熱失控 ValueError）、curie_margin
│   │   ├── test_dynamics_thermal_coupled.py  # T_coil=None 乾跑 bit-identical、熱耦合改變 i(t)/t_open
│   │   └── test_limits.py              # find_pullin_voltage（二分搜尋、逾 4×V_bus 拋錯）
│   │
│   └── report/               # 基準案例報告腳本與輸出圖
│       ├── generate_p0_report.py       # 基準案例模擬 + 三張瞬態圖
│       ├── generate_schematic.py       # 產生 docs/spec/ 三張示意圖：valve_schematic.png
│       │                               #   （PARAMS.md 對照圖）+ generate_winding_window()／
│       │                               #   generate_seal_land()（供 app.py 側欄顯示）
│       ├── generate_p1_coupled_report.py   # 乾跑 vs N₂ 耦合對照報告（P1）
│       ├── generate_p2_l6_envelope.py  # L6 電壓-溫度包絡 + 響應特性 4 + LN₂ 77 K 域外案例（P2）
│       ├── generate_p3_flyback_report.py   # 二極體/稽納/RC 三種續流電路關閉對照（P3）
│       ├── generate_p4_l2_pareto_report.py # L2 功率-保持力 Pareto 前緣（P4）
│       ├── generate_p4_l1_report.py    # L1 響應時間下界三組掃描（電壓/匝數/行程）（P4）
│       ├── generate_p5_l3_report.py    # L3 小型化可行域 + Buckingham π 群（P5）
│       ├── generate_p6_report.py       # L4 洩漏門檻 + L5 撞擊/壽命 + 關閉回彈漏量尖峰（P6）
│       ├── MODEL_NOTES.md              # 驗證結果、已知簡化項清單
│       ├── current_i_t.png
│       ├── position_x_t.png
│       ├── force_F_t.png
│       ├── dry_vs_coupled_x_t.png      # 乾跑 vs 耦合位移對照圖
│       ├── l6_envelope_v_t.png         # L6 電壓-溫度工作邊界圖（P2）
│       ├── p3_flyback_current_t.png    # 三種續流電路的電流衰減對照（P3）
│       ├── p4_l2_pareto.png            # L2 功率-保持力 Pareto（P4）
│       ├── p4_l1_sweeps.png            # L1 三組掃描 + 下界疊圖（P4）
│       ├── p5_l3_region.png            # L3 矩形可行域（P5）
│       ├── p6_l4_leak.png              # land 寬度 vs 漏率，兩閥座對照 + 規格線（P6）
│       ├── p6_l5_life.png              # 撞擊速度 vs 估計壽命 + 二極體/稽納工作點（P6）
│       └── p6_bounce_x_t.png           # 稽納續流下的關閉回彈位移軌跡，PCTFE 閥座（P6）
│
├── pyproject.toml            # pytest 設定（讓 solenoid_model 可被匯入）+ 相依宣告（numpy/scipy/matplotlib/streamlit）
├── .gitignore
└── .claude/settings.json     # 權限白名單設定
```

## 資料流向

`params.py`（定義）→ `baseline.py`（具體數值）→ `magnetics.py`（磁路物理）→ `dynamics.py`（耦合 ODE，呼叫 magnetics）→ `report/`（跑基準案例、產圖、寫報告）；`tests/` 貫穿每一層，各自對照解析解驗證。

## 分層架構與雙模式設計

```
第 1 層｜參數層    ValveParams（閥件硬體）｜fluid 物性 dataclass（流體身分+工況）
                  ↑ 刻意分離：同一顆閥可在不同流體/工況下使用
第 2 層｜物理函式層 magnetics.py（磁路）｜fluid.py（流量/流體力）｜dynamics 元件函式
第 3 層｜模擬層    coupled_rhs + simulate_opening（把第 2 層組成 ODE）
第 4 層｜使用者層  app.py（GUI）｜report/ 腳本
```

**雙模式規則**（乾跑 vs 流體耦合）：

| | 乾跑模式 | 耦合模式 |
|---|---|---|
| 第 3 層呼叫 | `simulate_opening(params)`（medium=None，向後相容） | `simulate_opening(params, medium=..., cond=...)`，明確傳入流體與工況 |
| 靜態壓力力 | `delta_P·A_seat`（簡化模式） | 由工況推導 `(P_up−P_down)·A_seat`；**`delta_P` 在此模式被忽略**（避免雙重帳本） |
| 流體力 F_flow | 無 | 疊加噴流動量反作用（氣/液分流） |

**「預設開啟」層級原則**：模式預設由**第 4 層（使用者層）決定**——GUI 與報告腳本預設啟用耦合（N₂、P_up=MEOP、P_down=真空、T₀=293K 攤在畫面可改，可一鍵切乾跑對照）；第 3 層以下無隱藏假設、無魔法預設（規格 §2.0 本體論：不得引入無法對應的隱藏自由度）。

熱域沿用同一雙模式原則：庫層（`thermal.py`、`dynamics.simulate_opening` 的 `T_coil` 參數）無隱藏溫度假設，`T_coil` 的預設值（293.15 K，即 20°C 定義點）落在第 4 層（GUI/報告腳本），庫層呼叫者必須明確傳入。

密封／回彈沿用同一原則：`seat=None`（預設）維持既有硬鉗制止擋行為，對現有呼叫者 bit-level 相同；傳入 `sealing.SeatPair` 才切換到回彈積分（`dynamics.simulate_opening`/`simulate_closing` 的 `seat` 參數，事件點 `v → −e·v` 後續積分至 `|v| < v_min` 或達 `n_max`）。連同 `medium`/`cond` 與 `T_coil`，`seat` 是**第四個正交模式旋鈕**：四者互相獨立、可任意組合，庫層本身不預設 `None` 以外的任何模式，「預設開啟」與否的選擇權一律留給第 4 層。

## 模組總表

各工作項目的產物與關鍵結論一律見 **`ROADMAP.md` 模組總表**（唯一來源，本檔不重複維護，避免兩表不同步）。

結構性備註：`validation/` 獨立資料夾為刻意省略——分層驗證已併入 `tests/`。
