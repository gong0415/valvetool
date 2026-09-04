# Model Notes — Baseline Validation and Known Simplifications

## Baseline case

| Parameter | Value | Source |
|---|---|---|
| V_bus | 28 V | spec baseline |
| R_coil_20C | 80 Ω | spec 4.1 typical range (20-200 Ω) |
| N_turns | 2000 | spec 4.2 typical range (500-5000) |
| A_gap | 30 mm² | chosen for a small thruster-valve-sized pole face |
| g0 | 0.35 mm | spec 4.2 stroke range + small residual gap |
| x_stroke | 0.30 mm | spec 4.2 typical range (0.05-0.5 mm) |
| m_arm | 1 g | spec 4.2 typical range (0.5-10 g) |
| k_spring | 20 N/mm | chosen to give a preload-dominated static balance |
| F_preload | 10 N | spec 4.2 typical range (5-50 N) |
| delta_P | 2.4 MPa | spec baseline MEOP |
| A_seat | 5 mm² | consistent with a small-thrust-class (Cv≈0.0007) attitude-control thruster valve |
| damping_coeff | 0.1 N·s/m | placeholder, not first-principles (see Known Simplifications) |

## Results

- Static pull-in threshold current at rest (gap=g0): 0.196 A (analytical force balance against F_preload + ΔP·A_seat = 22 N)
- Full coupled simulation: t_open = 5.03 ms — **meets the spec target of < 10 ms**
- Current signature: rises per RL charging until ~4.26 ms (matches the 0.196 A threshold almost exactly), then dips from a peak of ~0.198 A down to ~0.039 A as the armature completes its stroke — reproduces the back-EMF dip required by response characteristic 1 (current-signature diagnostics)
- No saturation flagged at the operating currents seen in this run (peak current ~0.2 A vs. the 0.274 A threshold for B=B_sat at gap=g0)

## Validation summary (5 layered checks, all passing)

1. Coenergy force vs. Maxwell stress tensor: exact match (< 1e-6 relative)
2. Decoupled RL circuit vs. analytical step response: matches within 1e-4 relative tolerance
3. Decoupled spring-mass vs. analytical SHM: matches within 1e-4 relative tolerance
4. Static equilibrium / rest-position clamp: armature stays pinned at x=0 exactly while current charges per the RL curve, releases at the analytically predicted threshold
5. Full coupled run: t_open within target, back-EMF dip present

## Known simplifications (not first-principles, per the design spec's honesty principle)

- `damping_coeff` is a placeholder viscous term, not derived from friction/seal physics
- Eddy current delay is excluded entirely (deferred to the L1 limit-analysis phase)
- Core reluctance uses a fixed `mu_r_core`, not a saturating B-H curve
- Magnetic saturation is only flagged (`saturation_check`), not modeled nonlinearly
- The rest-position mechanical stop is a simple clamp (freezes x, v at the boundary) by default (`seat=None`); passing a `sealing.SeatPair` switches on a real collision/bounce model (P6) — see below for what that model still leaves simplified
- Pressure force is a static load (`ΔP·A_seat`); dynamic flow force `F_flow(x, ṁ)` is deferred to the fluid-domain phase

## Sealing, impact, and materials — known limitations

`sealing.py` (L4 leak model), `impact.py` (Hertz seat impact + bounce), and `materials.py` add the sealing, seat-impact-life, and materials/environment domains (spec §2.4, §2.6, §3 L4/L5) and the seat-bounce mode of `simulate_opening`/`simulate_closing` (spec §2.2). Per the same honesty principle as above, this section records nine known limitations:

1. **The leak-branch magnitude is not trustworthy.** The engineering closure `h = Rq_c*(1-phi)` degenerates to `h -> Rq_c` as `p << H`, so the residual gap — and hence the leak rate — barely responds to load once you are away from the percolation threshold. The trustworthy output of `sealing.py` is limited to the sealed/not-sealed threshold location and the land width required to reach it (`sealing.land_width_for_seal`), never the leak-rate magnitude on the leaking branch. Confirmed with the user to ship with this caveat rather than block on a better contact model; revisit in a future round.
2. `PHI_PERCOLATION = 0.42` is the 2-D continuum percolation theory value; reported experimental values span 0.40-0.50, so it is EMPIRICAL.
3. `contact_area_fraction` (`A_real = F/H`, Bowden-Tabor) assumes fully plastic contact; the elastic-plastic transition is not modelled.
4. The transitional flow regime between viscous and free-molecular conduction is blended with `f = 1/(1+Kn)` — an engineering treatment, not a first-principles derivation.
5. Not modelled: surface waviness, particle embedment, and PCTFE viscoelastic cold-flow (real cold-flow would tend to improve sealing over time, so this model is conservative in that respect, not optimistic).
6. The L5 life-model coefficients (`ARCHARD_K`, `FAILURE_DEPTH_FRACTION`) are all EMPIRICAL; treat `l5_life_curve`'s cycle counts as order-of-magnitude only, never as a qualification number.
7. Every `materials.py` entry beyond the ASTM E595 TML/CVCM thresholds themselves (which are specification limits) is `LITERATURE` grade — textbook/handbook values, not measurements. Re-verify against the NASA outgassing database and the actual material certificates before any flight use.
8. The coefficient of restitution `e` (`SeatPair.e_restitution`) is an EMPIRICAL constant; its dependence on impact speed is not modelled.
9. **The shakedown/plastic branch test (`impact.is_shakedown`) uses the fully-plastic hardness `H`, not the first-yield pressure, and that is the unconservative direction.** Classical Hertz contact theory puts first subsurface yield at `p0 ~= 1.60*Y` while `H ~= 3*Y` (Tabor's rule), i.e. first yield actually starts around `p0 ~= 0.53*H`; every impact with `p_max` between roughly `0.53*H` and `H` is already accumulating plastic indentation, but because the test is `p_max < H` (not `p_max < 0.53*H`), `l5_life_curve` reports `N_cycle = +inf` for all of it -- this model does not track that band at all. Concretely, the shipped model claims infinite life for metal-seat impacts up to 0.280 m/s (71.3% of the diode's 0.393 m/s closing speed) and for PCTFE up to 0.111 m/s (28.2% of it) -- both verified via `scipy.optimize.brentq` on `impact.contact_pressure_max(m, v, seat, R) - seat.H`, independently confirmed by direct bisection on `impact.is_shakedown` and by a from-scratch Hertz-formula reimplementation, all three agreeing to 10 significant figures. This is the same branch criterion item 6's `ARCHARD_K`/`FAILURE_DEPTH_FRACTION` only ever apply downstream of: item 6 covers whether those coefficients are trustworthy once the plastic branch is entered, this item covers whether the model enters it at all.

**Consequence of design-doc finding 6** (the ODE bounce train's impact-speed ratios converge to `e` from above, not decay by exactly `e` — see `test_dynamics_bounce.py`): the armature reaches the travel stop with almost no magnetic force margin (~28.9 N against a ~28 N spring-plus-pressure load). That near-balance means bounce excursions are a significant fraction of the stroke — the first opening-side rebound alone can fly for tens of microns: 29.41 um (9.8% of the 300 um `x_stroke`) for the PCTFE seat, 54.53 um (18.2%) for the metal seat — and any response-time estimate that assumes the valve settles immediately on first contact will understate `t_settle`.

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
