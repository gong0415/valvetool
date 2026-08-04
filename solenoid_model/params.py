from dataclasses import dataclass


@dataclass
class ValveParams:
    # Electrical
    V_bus: float          # supply voltage, V
    R_coil_20C: float     # coil resistance at 20°C, ohm
    N_turns: float        # number of coil turns

    # Magnetic circuit
    A_gap: float          # pole face / air gap cross-sectional area, m^2
    g0: float             # total air gap at rest (de-energized), m
    x_stroke: float       # armature travel from rest to fully pulled-in, m
    l_core: float         # core magnetic path length, m
    mu_r_core: float      # core relative permeability
    B_sat: float          # core saturation flux density, T

    # Mechanical
    m_arm: float          # armature moving mass, kg
    k_spring: float       # spring stiffness, N/m
    F_preload: float      # spring preload force at rest (x=0), N

    # Simplified loads (placeholders, not first-principles)
    delta_P: float        # static pressure differential across seat, Pa
    A_seat: float         # seat area for static pressure force, m^2
    damping_coeff: float  # viscous damping coefficient, N*s/m

    # Fluid path geometry (P1)
    D_seat_bore: float = 0.52e-3   # seat flow bore diameter, m — sets flow capacity;
                                   # distinct from A_seat (static-pressure sealing area)
    C_d: float = 0.8               # discharge coefficient — EMPIRICAL (typ. 0.6-0.9),
                                   # not first-principles; see docs/spec/PARAMS.md

    # Thermal path (P2) — EMPIRICAL order-of-magnitude placeholders, not
    # measured values; see docs/spec/PARAMS.md
    G_th_cond: float = 0.2      # coil-to-mount conduction conductance, W/K
    emissivity: float = 0.8     # outer surface emissivity (potted/oxidized), -
    A_rad: float = 1.5e-3       # radiating outer surface area, m^2

    # Winding window (P4) — EMPIRICAL geometry, co-calibrated so
    # winding.coil_resistance(N_turns, ...) == R_coil_20C at N=2000; see
    # docs/spec/PARAMS.md. Analysis-layer only (L2); dynamics uses R_coil_20C.
    l_turn_mean: float = 0.031066014600891194    # mean turn length, m
    k_fill: float = 0.5                          # copper fill factor, -
    A_winding: float = 5.2190904529497215e-05    # winding window area, m^2

    # Sizing / L3 (P5) — EMPIRICAL, analysis-layer only; the ODE uses none of
    # these. See docs/spec/PARAMS.md.
    k_pack: float = 3.0      # packaging/housing mass factor, -
    J_max: float = 20e6      # continuous current-density ceiling, A/m^2

    # Sealing / impact geometry (P6) — EMPIRICAL, analysis-layer only; the ODE
    # uses neither unless a SeatPair is passed to simulate_opening/closing.
    # See docs/spec/PARAMS.md.
    w_land: float = 50e-6    # seal land width, m
    R_tip: float = 1.0e-3    # poppet tip radius of curvature, m
