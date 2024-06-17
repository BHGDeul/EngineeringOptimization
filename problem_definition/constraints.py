from constants.constants import *
import numpy as np
"""Constraints"""
def constraint_traction(x):
    # TODO: rewrite in negative-null form
    # Assume gamma_in and gamma_out are the nominal values
    gamma_in, gamma_out = x[1], x[2]
    a_elev_in = x[3]
    a_elev_out = x[4]
    # v_w_n = x[5]
    A_proj = x[0]

    T_out_elevation = 0.5 * rho * v_w_n ** 2 * A_proj * (
            np.cos(a_elev_out) - gamma_out) ** 2 * F_out

    T_in_elevation = 0.5 * rho * v_w_n ** 2 * A_proj * (
            1 + 2 * gamma_in * np.cos(a_elev_in) + gamma_in**2) * F_in

    return T_out_elevation/40E3 - 1, T_in_elevation/40E3 - 1

def constraint_reel_speed(x):
    v_w_n = x[5]
    return x[1] / max_reel_speed / v_w_n - 1