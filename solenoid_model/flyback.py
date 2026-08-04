"""Flyback circuit models (P3): coil terminal voltage during de-energization.

Never imports dynamics/magnetics — these are pure algebraic voltage helpers;
the ODE coupling (electrical decay + armature motion) lives in dynamics.py.
"""
from dataclasses import dataclass


@dataclass
class DiodeFlyback:
    V_fwd: float = 0.7      # forward voltage drop, V — simplified constant
                             # (real diodes have a log(I)-dependent drop;
                             # not first-principles, same treatment as C_d)


@dataclass
class ZenerFlyback:
    V_clamp: float = 30.0   # zener clamp voltage, V


@dataclass
class RCFlyback:
    R_snub: float            # snubber resistance, ohm — no default: an
                              # underdamped (too-small R_snub) combination
                              # is not a physically sane snubber design; see
                              # docs/spec/PARAMS.md for the damping-ratio
                              # formula to pick a sensible value
    C_snub: float             # snubber capacitance, F


def clamp_voltage(circuit, i):
    """Coil terminal voltage during diode/zener freewheel conduction.

    Returns the (negative) voltage clamp while current is still flowing
    (i > 0). Returns 0.0 once current has reached zero or below — the diode
    stops conducting (no reverse current), so the coil is open-circuit.
    """
    if i <= 0.0:
        return 0.0
    if isinstance(circuit, DiodeFlyback):
        return -circuit.V_fwd
    if isinstance(circuit, ZenerFlyback):
        return -circuit.V_clamp
    raise TypeError(
        f"clamp_voltage expects DiodeFlyback or ZenerFlyback, "
        f"got {type(circuit).__name__}"
    )


def rc_terminal_voltage(circuit, i, V_c):
    """Coil terminal voltage for the RC snubber path.

    KVL around the coil-snubber loop: L di/dt = -i*(R_coil+R_snub) - V_c,
    so the terminal voltage contributed by the snubber branch is
    -(i*R_snub + V_c).
    """
    return -(i * circuit.R_snub + V_c)
