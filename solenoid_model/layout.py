"""Physical layout: from ValveParams to a manufacturable dimension chain (P8).

Never imports dynamics/magnetics/fluid -- the same invariant sizing.py and
winding.py hold. Everything here is closed-form algebra over an existing
ValveParams; nothing feeds back into the ODE.

The split this module is careful about is DERIVED vs EMPIRICAL. Yoke
sections, wall thicknesses, spring geometry and armature thickness all
follow from physics already in ValveParams. Seat-body and fixed-pole
thicknesses do not -- they are packaging judgements, declared in
LAYOUT_EMPIRICAL with a stated reason rather than buried as literals, so
that reading a dimension tells you how much to trust it.
"""
import math
from dataclasses import dataclass

RHO_FE = 7870.0              # kg/m^3, soft iron

# EMPIRICAL layout inputs: packaging judgements, not derived quantities.
# Declared here rather than defaulted inside functions so that changing a
# judgement is a visible edit, following how params.py labels k_pack/C_d.
LAYOUT_EMPIRICAL = {
    "t_seat_body": {
        "value": 3.0e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "容納 0.60 mm 孔、密封 land 與閥座壓入配合的最小座體厚度；"
                  "無一致性檢查可攔截，須由 CAD 階段回頭確認",
    },
    "t_fixed_pole": {
        "value": 4.0e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "固定極需容納磁通轉向與繞線骨架端面支撐；厚度不足會在極面"
                  "根部先飽和，但本模型未解 3D 磁通分佈，故取封裝經驗值",
    },
    "t_manufacturing": {
        "value": 1.0e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "外殼可加工與剛性下限。磁性下限 0.35 mm 與結構下限 0.33-0.44 mm "
                  "皆遠低於此，故壁厚實際由製造決定而非物理",
    },
    "t_sleeve": {
        "value": 0.25e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "濕式銜鐵的非磁性隔離套（316L）。25 bar 下結構僅需 0.11 mm，"
                  "餘為加工餘裕。⚠️ 此厚度落在徑向磁路上但 magnetics.reluctance "
                  "未建模，故實際吸力低於模型報告值",
    },
    "clearance_spring": {
        "value": 0.1e-3, "unit": "m", "source": "EMPIRICAL",
        "reason": "彈簧全壓縮時與腔壁的餘隙，避免併圈干涉",
    },
}

# The spring design point solved in spec section 6: d=0.35 mm on a 1.8 mm
# coil diameter in 302 stainless (G=79 GPa) realises k_spring=4000 N/m at
# 6.35 active coils, with tau=389 MPa against a 700-900 MPa allowable and
# a spring index of 5.14 (4-12 is the easily-wound range). DERIVED from
# k_spring/F_preload, not an EMPIRICAL packaging guess -- which is why it
# is not in LAYOUT_EMPIRICAL.
SPRING_DESIGN = {"d_wire": 0.35e-3, "D_coil": 1.8e-3, "G": 79e9}


# --- Magnetic circuit sections -------------------------------------------

def core_radius(params):
    """Core / pole radius [m] from the pole face area."""
    return math.sqrt(params.A_gap / math.pi)


def pole_diameter(params):
    """Pole face diameter [m]."""
    return 2.0 * core_radius(params)


def yoke_area(params, B_pole=None):
    """Minimum yoke cross-section [m^2] that will not saturate before the pole.

    Flux continuity around the circuit: the yoke carries the same flux the
    pole face does, so

        B_pole * A_gap = B_yoke * A_yoke,  B_yoke <= B_sat
        => A_yoke >= A_gap * B_pole / B_sat

    `B_pole=None` uses the worst case B_pole = B_sat, giving
    A_yoke >= A_gap. That is deliberately the shipped default: it is
    independent of the operating point, so it stays valid when the bus
    voltage or turns count changes later. Passing the actual solved pole
    flux (magnetics.solve_flux_density) gives the tighter,
    operating-point-specific bound -- 16.97 mm^2 against the conservative
    20.0 mm^2 on the N2_25BAR_PH case.

    Assumes yoke and pole share a material. For a mixed-material circuit,
    scale by the ratio of the two B_sat values.
    """
    if B_pole is None:
        return params.A_gap
    if B_pole <= 0.0:
        raise ValueError(f"B_pole={B_pole} must be positive")
    if B_pole > params.B_sat:
        raise ValueError(
            f"B_pole={B_pole} T exceeds B_sat={params.B_sat} T: the pole "
            f"cannot carry that flux density, so sizing a yoke to it is "
            f"meaningless")
    return params.A_gap * B_pole / params.B_sat


def shell_wall_magnetic(params, R_inner, B_pole=None):
    """Shell wall thickness [m] whose annulus carries the yoke flux.

    Solves pi*((R_i+t)^2 - R_i^2) = A_yoke for t.
    """
    A = yoke_area(params, B_pole)
    return math.sqrt(A / math.pi + R_inner ** 2) - R_inner


def end_plate_thickness(params, B_pole=None):
    """End plate thickness [m] for radially outward flux.

    A washer carrying flux radially has cross-section A = 2*pi*r*t, which
    is SMALLEST at the smallest radius. The tightest point is therefore
    the core outer surface, and sizing there covers the whole plate.
    """
    A = yoke_area(params, B_pole)
    return A / (2.0 * math.pi * core_radius(params))


# --- Pressure boundary ----------------------------------------------------

def hoop_wall_thickness(P, R_inner, S_yield, SF):
    """Thin-wall hoop-stress thickness [m]: t = P*R/(S_yield/SF)."""
    if SF <= 0.0:
        raise ValueError(f"SF={SF} must be positive")
    if S_yield <= 0.0:
        raise ValueError(f"S_yield={S_yield} must be positive")
    return P * R_inner / (S_yield / SF)


@dataclass(frozen=True)
class WallThickness:
    """Three independent lower bounds on a wall, and which one binds."""
    magnetic: float
    structural: float
    manufacturing: float
    adopted: float
    reason: str


def shell_wall_thickness(params, R_inner, P, S_yield, SF,
                         t_manufacturing=None, B_pole=None):
    """Shell wall [m] as the largest of three independent lower bounds.

    Reporting all three rather than just the winner is the point: on this
    valve the magnetic (0.35 mm) and structural (0.33-0.44 mm) bounds are
    the same order and BOTH are below what can be machined, so the wall is
    set by manufacturing, not by physics. A single returned number would
    hide that.
    """
    if t_manufacturing is None:
        t_manufacturing = LAYOUT_EMPIRICAL["t_manufacturing"]["value"]
    t_mag = shell_wall_magnetic(params, R_inner, B_pole)
    t_str = hoop_wall_thickness(P, R_inner, S_yield, SF)
    candidates = {"magnetic": t_mag, "structural": t_str,
                  "manufacturing": t_manufacturing}
    name = max(candidates, key=candidates.get)
    return WallThickness(
        magnetic=t_mag, structural=t_str, manufacturing=t_manufacturing,
        adopted=candidates[name],
        reason=f"{name} bound governs "
               f"(magnetic {t_mag*1e3:.2f} mm, structural {t_str*1e3:.2f} mm, "
               f"manufacturing {t_manufacturing*1e3:.2f} mm)")


# --- Spring ---------------------------------------------------------------

def wahl_factor(C):
    """Wahl stress-correction factor for spring index C = D/d.

    Corrects for curvature and direct shear, both of which raise the real
    stress above the straight-torsion estimate. Tends to 1 as C grows.
    """
    if C <= 1.0:
        raise ValueError(f"spring index C={C} must exceed 1")
    return (4.0 * C - 1.0) / (4.0 * C - 4.0) + 0.615 / C


@dataclass(frozen=True)
class SpringGeometry:
    """Helical compression spring realising params.k_spring."""
    d_wire: float     # m
    D_coil: float     # m
    n_active: float   # active coils, -
    index: float      # C = D/d, -
    tau_max: float    # Pa, at full stroke, Wahl-corrected
    L_solid: float    # m
    L_free: float     # m


def spring_geometry(params, d_wire, D_coil, G, n_dead=2.0, clearance=None):
    """Solve the coil count that realises params.k_spring, and its stresses.

    k = G*d^4/(8*D^3*n)  =>  n = G*d^4/(8*D^3*k)

    Free length is solid height plus the working deflection (preload plus
    stroke) plus a clearance, so the spring never reaches solid at full
    stroke.
    """
    if clearance is None:
        clearance = LAYOUT_EMPIRICAL["clearance_spring"]["value"]
    if d_wire <= 0.0 or D_coil <= 0.0:
        raise ValueError(f"d_wire={d_wire}, D_coil={D_coil} must be positive")
    C = D_coil / d_wire
    n = G * d_wire ** 4 / (8.0 * D_coil ** 3 * params.k_spring)
    if n <= 0.0:
        raise ValueError(f"k_spring={params.k_spring} gives n={n} coils")
    F_max = params.F_preload + params.k_spring * params.x_stroke
    tau = wahl_factor(C) * 8.0 * F_max * D_coil / (math.pi * d_wire ** 3)
    L_solid = (n + n_dead) * d_wire
    deflection = params.F_preload / params.k_spring + params.x_stroke
    return SpringGeometry(
        d_wire=d_wire, D_coil=D_coil, n_active=n, index=C, tau_max=tau,
        L_solid=L_solid, L_free=L_solid + deflection + clearance)


# --- Armature -------------------------------------------------------------

def armature_thickness(params, rho=RHO_FE, D=None):
    """Armature disc thickness [m] implied by params.m_arm.

    DERIVED, not a packaging guess: t = m_arm / (rho * A). Defaults to a
    disc at the pole diameter.
    """
    if D is None:
        D = pole_diameter(params)
    A = math.pi * D ** 2 / 4.0
    return params.m_arm / (rho * A)


def armature_mass_consistency(params, t_arm, rho=RHO_FE, D=None):
    """Compare declared m_arm against the mass a given thickness implies.

    Reports rather than raises, following
    winding.coil_resistance_consistency: an armature may legitimately stop
    being a plain disc (a guide stem, a lightening bore), which decouples
    thickness from mass. The caller decides what divergence is acceptable.

    Returns (m_declared, m_geometric, rel_error).
    """
    if D is None:
        D = pole_diameter(params)
    A = math.pi * D ** 2 / 4.0
    m_geom = rho * A * t_arm
    return params.m_arm, m_geom, abs(m_geom - params.m_arm) / params.m_arm


# --- Assembly -------------------------------------------------------------

def axial_stack(params, t_seat_body=None, t_fixed_pole=None,
                d_wire=None, D_coil=None, G=None, rho=RHO_FE):
    """Ordered axial dimension chain [m], seat face upward.

    Returns {"segments": [(name, thickness), ...], "total": float}. Ordered
    so the list reads as the physical stack, not as a dict of parts.

    Spring parameters default to SPRING_DESIGN, the single verified design
    point (spec section 6) that the figure generator and GUI also draw
    from -- so the dimension chain here and the section drawing it feeds
    can never silently disagree about which spring they mean.
    """
    if t_seat_body is None:
        t_seat_body = LAYOUT_EMPIRICAL["t_seat_body"]["value"]
    if t_fixed_pole is None:
        t_fixed_pole = LAYOUT_EMPIRICAL["t_fixed_pole"]["value"]
    if d_wire is None:
        d_wire = SPRING_DESIGN["d_wire"]
    if D_coil is None:
        D_coil = SPRING_DESIGN["D_coil"]
    if G is None:
        G = SPRING_DESIGN["G"]
    spring = spring_geometry(params, d_wire, D_coil, G)
    segments = [
        ("閥座座體", t_seat_body),
        ("行程 x_stroke", params.x_stroke),
        ("銜鐵", armature_thickness(params, rho)),
        ("工作氣隙 g0", params.g0),
        ("固定極", t_fixed_pole),
        ("彈簧腔", spring.L_free),
        ("端板 x2", 2.0 * end_plate_thickness(params)),
    ]
    return {"segments": segments, "total": sum(t for _, t in segments)}


def envelope(params, coil_OD, R_inner=None, P=25e5, S_yield=205e6, SF=3.0,
             t_manufacturing=None, **stack_kwargs):
    """Overall package envelope [m]: outer diameter and length."""
    if R_inner is None:
        R_inner = coil_OD / 2.0
    wall = shell_wall_thickness(params, R_inner, P, S_yield, SF,
                                t_manufacturing)
    stack = axial_stack(params, **stack_kwargs)
    return {"OD": 2.0 * (R_inner + wall.adopted),
            "L": stack["total"],
            "wall": wall,
            "stack": stack}
