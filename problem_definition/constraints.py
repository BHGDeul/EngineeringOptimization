from constants.constants import *
import numpy as np
from problem_definition.objective import power_output, aero_module

"""Constraints"""
def constraint_traction(x):
    gamma_in, gamma_out = x[1], x[2]
    a_elev_in = x[3]
    a_elev_out = x[4]
    v_w_n = x[5]
    A_proj = x[0]

    # F_out, F_in = aero_module(x)

    T_out_elevation = 0.5 * rho * v_w_n ** 2 * A_proj * (
            np.cos(a_elev_out) - gamma_out) ** 2 * F_out

    T_in_elevation = 0.5 * rho * v_w_n ** 2 * A_proj * (
            1 + 2 * gamma_in * np.cos(a_elev_in) + gamma_in**2) * F_in

    return T_out_elevation/max_traction - 1, T_in_elevation/max_traction - 1

### For linear constraints
A = np.zeros((6, 6))
# A[5, 5] = 1 / max_reel_speed / v_w_n - 1
A[4, 4] = 1 / (np.pi / 2) - 1

lb = np.zeros((6, ))
lb[5] = -1
lb[4] = -1
ub = np.zeros((6, ))

def constraint_reel_speed(x):
    return x[1] / max_reel_speed - 1

def constraint_power(x):
    return - power_output(x) / np.sqrt(power_output(x)**2) + 1