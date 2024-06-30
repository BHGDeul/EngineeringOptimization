from constants.constants import rho, eff_out, eff_in, baseline_power, baseline_area, F_out, F_in
import numpy as np
from aero_module.aero_main_function import main_aero_function


def aero_module(x):
    """
    Module to make use of the VSM aerodynamic analysis module to determine F_out and F_in
    for the current design.

    THIS FUNCTION IS NOT IMPLEMENTED IN THE OPTIMIZATION TO SAVE RUN TIME
    :param x: vector of design variables
    :return: F_out, F_in: the force factors when reeling out/in the kite.
    """
    A_proj = x[0]

    AoA_range_out = np.array([18])
    AoA_range_in = np.array([x[6]])

    Points = 1000
    Kite_segments = 6
    N_split = 5
    V_wind = 10

    F_out, F_in = main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_out, AoA_range_in, V_wind, False)[2:4]
    return F_out, F_in

def power_output(x):
    """
    Calculate the power output of the design
    :param x: vector of design variables
    :return: power output in Watt
    """
    gamma_in, gamma_out = x[1], x[2]
    a_elev_in = x[3]
    a_elev_out = x[4]
    A_proj = x[0]
    v_w_n = x[5]

    # F_out, F_in = aero_module(x)

    P_w = 0.5 * v_w_n ** 3 * rho  # Wind Power


    power = P_w * A_proj * (
            eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
            (F_in / eff_in) * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) * (
                    (gamma_out * gamma_in) / (gamma_out + gamma_in))

    return power

def objective(x):
    '''
    Objective function
    :param x: vector of design variables
    :return: objective function value
    '''
    P_ref = baseline_power
    A_ref = baseline_area
    return (P_ref / (1 + power_output(x))) + x[0] / A_ref