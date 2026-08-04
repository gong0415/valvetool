import math

MU_0 = 4 * math.pi * 1e-7  # vacuum permeability, H/m


def reluctance(gap, params):
    r_gap = gap / (MU_0 * params.A_gap)
    r_core = params.l_core / (params.mu_r_core * MU_0 * params.A_gap)
    return r_gap + r_core


def inductance(gap, params):
    return params.N_turns ** 2 / reluctance(gap, params)


def dinductance_dgap(gap, params):
    r_total = reluctance(gap, params)
    dr_dgap = 1.0 / (MU_0 * params.A_gap)
    return -params.N_turns ** 2 * dr_dgap / r_total ** 2


def flux_density(gap, i, params):
    flux = params.N_turns * i / reluctance(gap, params)
    return flux / params.A_gap


def magnetic_force_closing(gap, i, params):
    return -0.5 * i ** 2 * dinductance_dgap(gap, params)


def saturation_check(gap, i, params):
    return flux_density(gap, i, params) > params.B_sat
