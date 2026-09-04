from dataclasses import replace

from streamlit.testing.v1 import AppTest

from solenoid_model import app
from solenoid_model.baseline import BASELINE_PARAMS


def test_compute_result_succeeds_for_baseline_params():
    result = app.compute_result(BASELINE_PARAMS)
    assert result is not None
    assert 0.0049 < result.t_open < 0.0052  # known ~5.027 ms
    assert result.saturated is False


def test_compute_result_returns_none_when_armature_never_pulls_in():
    params = replace(BASELINE_PARAMS, F_preload=10000.0)
    result = app.compute_result(params, t_max=0.01)
    assert result is None


def test_app_ui_default_coupled_shows_slower_t_open():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    assert not at.exception
    assert len(at.sidebar.number_input) >= 20  # 15 既有 + 流道幾何 2 + 工況 3
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["t_open"] == "5.031 ms"
    assert len(at.error) == 0


def test_app_ui_dry_toggle_restores_baseline():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    coupling = next(cb for cb in at.sidebar.checkbox if cb.label == "啟用流體耦合")
    coupling.set_value(False)
    at.run()
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["t_open"] == "5.027 ms"


def test_app_ui_liquid_selection_shows_flashing_warning():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    medium_select = next(sb for sb in at.sidebar.selectbox if sb.label == "流體")
    medium_select.set_value("LN₂ (77K)")
    at.run()
    assert not at.exception
    assert any("閃蒸" in w.value for w in at.warning)


def test_app_ui_shows_error_for_extreme_preload():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    preload_input = next(
        ni for ni in at.sidebar.number_input if ni.label == "彈簧預載力 (N)"
    )
    preload_input.set_value(10000.0)
    at.run()
    assert not at.exception
    # 開啟動態頁籤本身無 metric；τ_e/i_th 屬於 L1 頁籤，獨立於 t_open 計算，
    # 兩頁籤同一次 run() 都會渲染，故仍會出現（Task 6 起）。
    metrics = {m.label: m.value for m in at.metric}
    assert "t_open" not in metrics
    errors = [e.value for e in at.error]
    assert any("不會吸合" in e for e in errors)          # 開啟動態頁籤
    assert any("L1 無界可算" in e for e in errors)       # 極限掃描頁籤：此參數下 L1 判斷同樣不會開


def test_app_ui_thermal_default_shows_r20_and_keeps_baseline():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["熱態線圈電阻"] == "80.00 Ω"   # T=293.15 → 恰為 R_coil_20C
    assert metrics["t_open"] == "5.031 ms"        # 預設行為不變
    assert len(at.error) == 0


def test_app_ui_hot_coil_slows_opening_and_updates_r():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    t_default = float({m.label: m.value for m in at.metric}["t_open"].split()[0])
    T_input = next(ni for ni in at.sidebar.number_input
                   if ni.label == "線圈溫度 T_coil (K)")
    T_input.set_value(343.15)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["熱態線圈電阻"] == "95.72 Ω"
    assert float(metrics["t_open"].split()[0]) > t_default


def test_app_ui_pullin_button_shows_voltage_metric():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    button = next(b for b in at.button if b.key == "btn_pullin")
    button.set_value(True)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert "V_pull-in" in metrics
    # 常溫基準 V_pull-in ≈ 15.7 V,遠低於 28 V → 不該出現包絡警告
    assert float(metrics["V_pull-in"].split()[0]) < 28.0


def test_app_ui_pullin_button_search_failure_shows_error(monkeypatch):
    # ValueError 路徑無法用真實參數觸發(按鈕僅在主模擬成功時渲染,
    # 而 t_open 對電壓單調)——monkeypatch 注入,驗證 try/except 存在
    from solenoid_model import limits

    def raise_no_pullin(*args, **kwargs):
        raise ValueError("no pull-in below 112.0 V within t_max=0.03 s")

    monkeypatch.setattr(limits, "find_pullin_voltage", raise_no_pullin)
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    button = next(b for b in at.button if b.key == "btn_pullin")
    button.set_value(True)
    at.run()
    assert not at.exception
    assert any("拉入電壓搜尋失敗" in e.value for e in at.error)


def test_app_ui_sidebar_has_winding_and_seat_expanders():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    assert not at.exception
    labels = [ni.label for ni in at.sidebar.number_input]
    for expected in ("繞線窗口面積 A_winding (m²)", "銅填充率 k_fill",
                     "平均匝長 l_turn_mean (m)", "密封 land 寬度 w_land (m)",
                     "poppet 端面曲率半徑 R_tip (m)"):
        assert expected in labels
    seat_select = next(sb for sb in at.sidebar.selectbox
                       if sb.label == "閥座材料對")
    assert seat_select.value == "PCTFE/440C"
    # 基準值下窗口推導 R 與 R_coil_20C 一致 → 不得出現不一致警告
    assert not any("繞線窗口" in w.value for w in at.warning)
    # 頭條行為不變
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["t_open"] == "5.031 ms"


def test_app_ui_custom_seat_pair_shows_material_inputs():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    seat_select = next(sb for sb in at.sidebar.selectbox
                       if sb.label == "閥座材料對")
    seat_select.set_value("自訂")
    at.run()
    assert not at.exception
    sb_labels = [sb.label for sb in at.sidebar.selectbox]
    assert "閥座材料 A" in sb_labels and "閥座材料 B" in sb_labels
    assert any(ni.label == "恢復係數 e" for ni in at.sidebar.number_input)


def test_app_ui_has_five_tabs_and_baseline_behavior_unchanged():
    # P8 Task 6 added a fifth tab ("實體佈局"); the label list below pins
    # the new correct baseline rather than a length/subset check, so
    # accidental tab drift is still caught.
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    assert not at.exception
    assert [t.label for t in at.tabs] == ["開啟動態", "關閉/續流",
                                          "極限掃描 L1–L3", "密封/壽命 L4–L5",
                                          "實體佈局"]
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["t_open"] == "5.031 ms"


def test_compute_opening_bounce_returns_decaying_impact_train():
    from solenoid_model import materials
    br = app.compute_opening_bounce(BASELINE_PARAMS, materials.PCTFE_ON_440C)
    assert br.n_bounce >= 2
    # P6 發現 6：連續撞擊速度比（v[n+1]/v[n]）從上方收斂到 e=0.4
    # （量測序列 0.605 → 0.430 → 0.409 → 0.403 → 0.401，與任務簡報一致）。
    # 簡報 Step 1 原將此比值序列的首項 0.605 誤植為 v_impacts[0] 本身；
    # 實際首次觸擊速度（與硬鉗制模擬在 t_open 當下的速度完全吻合）為
    # ~0.854 m/s，見 task-4-report.md 的量測證據。
    assert br.v_impacts[0] > br.v_impacts[-1]
    assert 0.80 < br.v_impacts[0] < 0.90
    ratios = [b / a for a, b in zip(br.v_impacts, br.v_impacts[1:])]
    assert ratios[0] > ratios[-1] > 0.4  # 比值從上方收斂到 e=0.4


def test_compute_closing_diode_vs_zener_matches_p6_report():
    from solenoid_model import flyback
    diode = app.compute_closing(BASELINE_PARAMS, flyback.DiodeFlyback())
    zener = app.compute_closing(BASELINE_PARAMS, flyback.ZenerFlyback())
    assert 0.062 < diode.t_close < 0.064      # 62.748 ms
    assert 0.017 < zener.t_close < 0.018      # 17.456 ms
    assert 3.0 < diode.t_close / zener.t_close < 5.0   # 規格宣稱 3–5 倍
    assert 0.38 < diode.v_impact < 0.41       # 0.3934 m/s
    assert 1.50 < zener.v_impact < 1.57       # 1.5338 m/s


def test_app_ui_closing_button_runs_and_compares_to_diode():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    circ = next(sb for sb in at.selectbox if sb.label == "續流電路")
    circ.set_value("稽納")
    cmp_cb = next(cb for cb in at.checkbox if cb.label == "與二極體對照")
    cmp_cb.set_value(True)
    next(b for b in at.button if b.key == "btn_closing").set_value(True)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["t_close"] == "17.456 ms"
    assert "觸座速度" in metrics
    assert any("倍" in c.value for c in at.caption)  # 對照比值已顯示


def test_app_ui_opening_bounce_toggle_adds_metrics():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=30)
    at.run()
    assert "回彈次數" not in {m.label for m in at.metric}  # 預設關
    toggle = next(cb for cb in at.checkbox if cb.label == "啟用閥座回彈")
    toggle.set_value(True)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert "回彈次數" in metrics
    assert "首次觸擊速度" in metrics
    assert metrics["t_open"] == "5.031 ms"  # 頭條不因回彈旋鈕改變
    assert metrics["熱態線圈電阻"] == "80.00 Ω"


def test_parse_float_list_and_l1_cases():
    assert app.parse_float_list("18, 22,28") == [18.0, 22.0, 28.0]
    import pytest as _pytest
    with _pytest.raises(ValueError):
        app.parse_float_list("18, abc")
    cases = app.l1_sweep_cases("turns_window", [1000.0], BASELINE_PARAMS)
    label, p, R = cases[0]
    assert p.N_turns == 1000.0
    from solenoid_model import winding
    assert R == winding.coil_resistance(1000.0, p.A_winding, p.k_fill,
                                        p.l_turn_mean)
    cases_v = app.l1_sweep_cases("voltage", [18.0], BASELINE_PARAMS)
    assert cases_v[0][1].V_bus == 18.0
    assert cases_v[0][2] == BASELINE_PARAMS.R_coil_20C


def test_app_ui_l1_instant_metrics_and_sweep_button():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=60)
    at.run()
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["τ_e 電氣時間常數"] == "5.200 ms"
    assert metrics["L1 響應時間下限"] == "4.381 ms"
    assert any("不是上下界" in c.value for c in at.caption)
    next(b for b in at.button if b.key == "btn_l1_sweep").set_value(True)
    at.run()
    assert not at.exception
    assert len(at.dataframe) >= 1


def test_app_ui_l1_bound_tracks_cond_pressure_not_stale_delta_p():
    # Finding 1 (P7 review): Tab 3's L1 bound must describe the SAME
    # operating point as Tab 1's t_open when fluid coupling is on. Set
    # P_up far below the default delta_P (2.4 MPa baseline == default
    # P_up so they agree only by coincidence today) and confirm the L1
    # bound is still a genuine lower bound on the coupled t_open.
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=60)
    at.run()
    P_up_input = next(ni for ni in at.sidebar.number_input
                      if ni.label == "上游壓力 P_up (Pa)")
    P_up_input.set_value(0.2e6)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    t_open = float(metrics["t_open"].split()[0])
    bound = float(metrics["L1 響應時間下限"].split()[0])
    assert bound < t_open, (
        f"L1 bound {bound} ms is not a lower bound on t_open {t_open} ms "
        "-- Tab 3 is using stale delta_P instead of the cond pressure"
    )


def test_app_ui_l2_l3_sections_render_with_baseline_values():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=60)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["B_sat 保持力天花板"] == "43.09 N"
    assert metrics["最小可行縮放 s*"] == "1.014"
    assert metrics["壓力天花板"] == "8.62 MPa"
    assert metrics["整閥質量估計"] == "79.0 g"
    # π 群表格存在（st.table → at.table）
    assert len(at.table) >= 1
    # EMPIRICAL 警語存在
    assert any("EMPIRICAL" in c.value for c in at.caption)


def test_app_ui_l4_metrics_and_materials_criteria():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=60)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["密封 land 寬度上限"] == "66.08 µm"   # PCTFE/440C 門檻
    assert metrics["接觸應力"] == "55.5 MPa"
    # 注意:不可寫 `"密封" in ...`——「不密封」也含這兩字(已實測 True)
    assert metrics["目前 w_land 判定"] == "密封"
    # 發現 2 的誠實警語
    assert any("門檻位置" in w.value for w in at.warning)
    # 材料準則:內建 PCTFE/440C 出氣皆過 → 有 ✅;推進劑相容性表存在
    assert any("出氣" in c.value or "outgassing" in c.value.lower()
               for c in at.caption) or len(at.table) >= 2


def test_app_ui_l4_metal_seat_shows_not_sealed():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=60)
    at.run()
    seat_select = next(sb for sb in at.sidebar.selectbox
                       if sb.label == "閥座材料對")
    seat_select.set_value("17-4PH/440C")
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["密封 land 寬度上限"] == "1.65 µm"
    assert metrics["目前 w_land 判定"] == "不密封"


def test_compute_bounce_spike_zener_n2_matches_prototype():
    from solenoid_model import flyback, fluid, materials
    cond = fluid.FlowConditions(P_up=2.4e6, P_down=0.0, T0=293.0)
    spike = app.compute_bounce_spike(BASELINE_PARAMS, materials.PCTFE_ON_440C,
                                     flyback.ZenerFlyback(), cond, fluid.N2)
    assert spike["n_bounce"] >= 2
    assert 2e-4 < spike["t_open_total"] < 5e-4       # 實測 3.270e-4 s
    assert 5e-9 < spike["mass_leaked"] < 2e-8        # 實測 1.076e-8 kg
    assert spike["sealed_leak"] == 0.0               # PCTFE 基準 w_land 密封


def test_app_ui_l5_curve_and_spike_button():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=60)
    at.run()
    metrics = {m.label: m.value for m in at.metric}
    assert "Hertz 峰值接觸應力" in metrics
    assert any("數量級" in w.value or "EMPIRICAL" in w.value
               for w in at.warning)
    next(b for b in at.button if b.key == "btn_spike").set_value(True)
    at.run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert "關閉回彈次數" in metrics
    assert "回彈離座總時間" in metrics
    assert "回彈漏出質量" in metrics


def test_app_ui_spike_button_disabled_for_liquid_medium():
    at = AppTest.from_file("solenoid_model/app.py", default_timeout=60)
    at.run()
    medium_select = next(sb for sb in at.sidebar.selectbox if sb.label == "流體")
    medium_select.set_value("水 (20°C)")
    at.run()
    assert not at.exception
    spike_buttons = [b for b in at.button if b.key == "btn_spike"]
    assert spike_buttons == [] or spike_buttons[0].disabled
