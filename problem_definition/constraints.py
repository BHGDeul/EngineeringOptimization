from constants import data
import numpy as np
"""Constraints"""
def constraint_traction(x):
    # TODO: rewrite in negative-null form
    # Assume gamma_in and gamma_out are the nominal values
    gamma_in, gamma_out = x[0], x[1]
    a_elev_out, a_elev_in = x[2], x[3]
    v_w_n = x[5]

    T_out_elevation = 0.5 * data['rho'] * v_w_n ** 2 * x[0] * (
            np.cos(a_elev_out) - gamma_out) ** 2 * data['F_out']

    T_in_elevation = 0.5 * data['rho'] * v_w_n ** 2 * x[0] * (
            1 + 2 * gamma_in * np.cos(a_elev_in) + gamma_in**2) * data['F_in']

    return T_out_elevation, T_in_elevation

def constraint_reel_speed(x):
    gamma_in, gamma_out = x[1], x[2]
    v_w_n = x[5]
    return gamma_in, gamma_out / data['max_reel_speed'] / v_w_n - 1