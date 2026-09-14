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
| `layout_actuation.png` | 作動前後對照圖（閥關 x=0 ↔ 閥開 x=x_stroke,含間隙放大與流體路徑） |
| `layout_architecture.png` | 四物理域架構圖與跨域耦合點（方塊數字由 params 推導,流量為 `fluid.mdot_gas` 實算） |

上列前三張示意圖由 `solenoid_model/report/generate_schematic.py` 產生,
後三張 P8 佈局圖由 `solenoid_model/report/generate_layout.py` 產生,皆不手工編輯。

佈局圖的幾何全部由 `ValveParams`／`layout.py` 推導（線圈斷面由 `A_winding`
反解、軸向鏈等於 `layout.axial_stack`、作動狀態遵循 `gap = g0 - x`）,
不得在繪圖端硬寫尺寸常數,並由測試斷言幾何而非僅斷言 PNG 產出。

### 動件組拓樸（與 `valve_schematic.png` 一致）

由下而上:閥座 → **閥芯 poppet** → 導杆 → **銜鐵** → 工作氣隙 → 固定極
（線圈包在銜鐵與固定極外側）。閥芯、導杆、銜鐵是同一個動件組,整體平移
`x`;密封由閥芯完成,磁吸由銜鐵承受。

**彈簧套在導杆上,位於閥芯上方、銜鐵下方**:上端頂在固定彈簧座（不動）,
下端頂在閥芯頂面（隨 `x` 上移）,故動件上移 `x` 即壓縮 `x`。這不是任意的
封裝選擇,而是 `dynamics` 的力學要求——`spring_force = F_preload + k*x` 在
`net_mechanical_force` 中取負,即彈簧力向下且隨行程增大,壓縮彈簧只能被動件
頂向上方的固定面。若把彈簧放在固定極之後（早期版本如此）,它會夾在兩個
不動面之間:既不能改變長度,也隔著 4 mm 實心極碰不到銜鐵。

`d_stem`／`d_poppet`／`t_poppet` 是 `LAYOUT_EMPIRICAL` 的封裝判斷（動力學
只需要 `m_arm` 與 `A_seat`,從未給定這三個尺寸）,但少了它們就畫不出彈簧的
力路。

### 流體路徑（作動圖）

氣體沿 0.60 mm 孔上行,在閥芯端面與座面之間的**簾幕**（環狀縫隙）轉為徑向
流出。節流面是簾幕 `π·D·x`,不是孔本身——`fluid.effective_area` 取
`min(π·D·x, π·D²/4)`。兩者在 `x = D/4` 相等,而本案 `x_stroke` 恰為
0.15 mm = D/4,所以全開時簾幕剛好追上孔面積 0.283 mm²,**再多行程也買不到
流量**（`cases.py` sizing chain 第 2 步）。下游真空使壓比遠低於
`r_crit(N₂)=0.528`,恆為壅塞流 1.302 g/s。

閥關時圖上不畫流線（流量為零）,改標 25 bar 壓在已坐封閥芯上的 0.71 N——
該力與彈簧同向,是失效關閉的一部分。

## 修改規格書

依規格書 §6.0 增量更新指引執行:先盤點現有實作產出對照表 → 確認優先順序 → 在特徵測試保護下逐項實作 → 本體論回溯檢查 → 逐項 commit。
