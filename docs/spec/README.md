# 任務規格書

電磁閥開發工具的第一性原理分析規格,以及參數定義手冊與示意圖。

## 內容

| 檔案 | 說明 |
|---|---|
| `solenoid_valve_first_principles.md` | 任務規格書:任務目標、第一性原理分解（§2,含 §2.0 模型本體論）、L1–L6 物理極限分析要求（§3）、關鍵參數清單（§4）、響應特性分析（§5）、執行與交付要求（§6,含 §6.0 增量更新指引）、參考錨點（§7） |
| `PARAMS.md` | `ValveParams`／`fluid.py`／`winding.py`／密封參數的定義手冊,對照下列示意圖 |
| `valve_schematic.png` | 整閥剖面示意圖 |
| `winding_window_schematic.png` | 繞線窗口示意圖（`A_winding`／`k_fill`／`l_turn_mean`） |
| `seal_land_schematic.png` | 密封 land 放大示意圖（`w_land`／`R_tip`／`A_seat`） |
| `layout_section.png` | 實體佈局剖面圖（P8,按真實 mm 比例） |
| `layout_actuation.png` | 作動前後對照圖（閥關 x=0 ↔ 閥開 x=x_stroke,含氣隙放大） |
| `layout_architecture.png` | 四物理域架構圖與跨域耦合點 |

上列前三張示意圖由 `solenoid_model/report/generate_schematic.py` 產生,
後三張 P8 佈局圖由 `solenoid_model/report/generate_layout.py` 產生,皆不手工編輯。

佈局圖的幾何全部由 `ValveParams`／`layout.py` 推導（線圈斷面由 `A_winding`
反解、軸向鏈等於 `layout.axial_stack`、作動狀態遵循 `gap = g0 - x`）,
不得在繪圖端硬寫尺寸常數,並由測試斷言幾何而非僅斷言 PNG 產出。

## 修改規格書

依規格書 §6.0 增量更新指引執行:先盤點現有實作產出對照表 → 確認優先順序 → 在特徵測試保護下逐項實作 → 本體論回溯檢查 → 逐項 commit。
