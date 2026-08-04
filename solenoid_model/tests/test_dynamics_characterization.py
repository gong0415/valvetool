"""Characterization tests (spec §6.0.3): freeze the dry-run baseline
before any dynamics.py modification. These must stay green for the entire
the fluid-coupling effort — a failure here means the refactor changed dry
behavior."""
import math

import numpy as np

from solenoid_model import dynamics
from solenoid_model.baseline import BASELINE_PARAMS


def test_dry_baseline_t_open_frozen():
    sol = dynamics.simulate_opening(BASELINE_PARAMS)
    assert math.isclose(sol.t_events[0][0], 5.026936705082159e-3, rel_tol=1e-9)


def test_dry_baseline_current_signature_frozen():
    sol = dynamics.simulate_opening(BASELINE_PARAMS)
    i_t = sol.y[0]
    assert math.isclose(i_t.max(), 0.198312, rel_tol=1e-4)
    assert math.isclose(i_t[-1], 0.038680, rel_tol=1e-4)
