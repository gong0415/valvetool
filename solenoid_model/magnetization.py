"""Single-valued B-H magnetization curves (P8).

Pure data + interpolation: these types answer only "what B does this H
give, and vice versa". Solving a magnetic circuit with one of them lives
in magnetics.py, so the curve source and the circuit solver stay
independent -- a tabulated measurement and the analytic law are
interchangeable at the same interface.

Hysteresis is deliberately absent. These curves are single-valued, so
they model saturation but not remanence: B depends on H alone, with no
magnetization history. For this valve that is the right trade -- the
remanent holding force in a gap-dominated circuit is ~0.02 N against a
2.2 N spring preload, two orders of magnitude short of a non-reseating
risk, while a hysteresis model would turn the ODE from memoryless into
history-dependent. See the design spec section 3A.5.
"""
import math
from dataclasses import dataclass

MU_0 = 4.0 * math.pi * 1e-7  # vacuum permeability, H/m


@dataclass(frozen=True)
class AnalyticBH:
    """Frohlich-Kennelly saturation law.

        B(H) = mu_i*H / (1 + mu_i*H/B_sat),   mu_i = mu_r * MU_0

    Introduces NO new degrees of freedom: both constants come from
    existing ValveParams fields (mu_r_core, B_sat), which is what keeps
    this compatible with the spec's ontology rule (section 2.0). It
    reduces to the linear model as H->0 and asymptotes to B_sat as
    H->infinity, so it is a strict refinement of the constant-mu_r model
    rather than a different one.
    """
    mu_r: float
    B_sat: float

    def __post_init__(self):
        if not (math.isfinite(self.mu_r) and self.mu_r > 0.0):
            raise ValueError(f"mu_r={self.mu_r} must be finite and positive")
        if not (math.isfinite(self.B_sat) and self.B_sat > 0.0):
            raise ValueError(f"B_sat={self.B_sat} must be finite and positive")

    @property
    def mu_i(self):
        """Initial permeability [H/m]."""
        return self.mu_r * MU_0

    def B_of_H(self, H):
        """Flux density [T] at field strength H [A/m]."""
        if H < 0.0:
            raise ValueError(f"H={H} must be non-negative")
        return self.mu_i * H / (1.0 + self.mu_i * H / self.B_sat)

    def H_of_B(self, B):
        """Field strength [A/m] for flux density B [T].

        Infinite at and above B_sat: no finite field reaches saturation
        under this law. Callers bisecting on B must keep the bracket
        strictly below B_sat.
        """
        if B < 0.0:
            raise ValueError(f"B={B} must be non-negative")
        if B >= self.B_sat:
            return math.inf
        return B / (self.mu_i * (1.0 - B / self.B_sat))


@dataclass(frozen=True)
class TabulatedBH:
    """Measured (H, B) points with piecewise-linear interpolation.

    Same interface as AnalyticBH, so a solver cannot tell them apart.
    Follows thermal.py's low-temperature table: strictly ascending
    points, validated at construction, and CLAMPED rather than
    extrapolated outside the table -- extrapolating a B-H curve past its
    last measured point produces unbounded, non-physical B.
    """
    points: tuple

    def __post_init__(self):
        if len(self.points) < 2:
            raise ValueError(
                f"points has {len(self.points)} entry; needs at least 2 "
                "(H, B) points to interpolate")
        for H, B in self.points:
            if not (math.isfinite(H) and math.isfinite(B)):
                raise ValueError(f"point ({H}, {B}) must be finite")
            if H < 0.0 or B < 0.0:
                raise ValueError(f"point ({H}, {B}) must be non-negative")
        Hs = [H for H, _ in self.points]
        Bs = [B for _, B in self.points]
        if any(h1 <= h0 for h0, h1 in zip(Hs, Hs[1:])):
            raise ValueError(f"H values {Hs} must be strictly ascending")
        if any(b1 <= b0 for b0, b1 in zip(Bs, Bs[1:])):
            raise ValueError(f"B values {Bs} must be strictly ascending")

    @property
    def B_sat(self):
        """Largest representable B [T] -- the table's last point."""
        return self.points[-1][1]

    def B_of_H(self, H):
        if H < 0.0:
            raise ValueError(f"H={H} must be non-negative")
        return _interp(H, self.points)

    def H_of_B(self, B):
        if B < 0.0:
            raise ValueError(f"B={B} must be non-negative")
        if B >= self.B_sat:
            return math.inf
        return _interp(B, tuple((b, h) for h, b in self.points))


def _interp(x, points):
    """Piecewise-linear interpolation, clamped at both ends."""
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return points[-1][1]


def from_params(params):
    """AnalyticBH built from a ValveParams' existing magnetic fields."""
    return AnalyticBH(mu_r=params.mu_r_core, B_sat=params.B_sat)
