from solenoid_model.params import ValveParams

# Baseline validation case. Values chosen within the typical ranges
# of the source spec (section 4.1-4.3), not a specific commercial product.
# Numerically pre-verified: static pull-in threshold current ~0.196 A at
# gap=g0 is reachable under steady-state I=V/R=0.35 A without exceeding
# B_sat at that current; full coupled simulation reaches t_open ~5.0 ms
# with a visible back-EMF current dip (see Task 6/7 notes).
# Fluid path: D_seat_bore=0.14 mm with C_d=0.8 gives full-open Cv=0.0007
# (bore-limited above x=0.035 mm) — re-targeted from the original
# Cv~0.01 spec §6.4 case to a small-thrust-class attitude-control thruster
# valve, in the flow-coefficient range of commercial ~5 lbf (22 N)
# redundant-seat monopropellant thruster valves (Moog datasheet reference).
BASELINE_PARAMS = ValveParams(
    V_bus=28.0,
    R_coil_20C=80.0,
    N_turns=2000.0,
    A_gap=3.0e-5,
    g0=3.5e-4,
    x_stroke=3.0e-4,
    l_core=0.05,
    mu_r_core=4000.0,
    B_sat=1.9,
    m_arm=0.001,
    k_spring=20000.0,
    F_preload=10.0,
    delta_P=2.4e6,
    A_seat=5.0e-6,
    damping_coeff=0.1,
    D_seat_bore=0.14e-3,
    C_d=0.8,
)
