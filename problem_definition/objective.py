from constants.constants import rho, eff_out, eff_in, baseline_power, baseline_area, F_out, F_in
import numpy as np
from aero_module.aero_main_function import main_aero_function


def aero_module(x):
    A_proj = x[0]

    Points = 1000
    Kite_segments = 12
    N_split = 5
    V_wind = 10

    AoA_range_out = np.arange(8, x[6], 1)
    AoA_range_in = np.arange(-1, x[7], 1)

    F_out, F_in = main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_out, AoA_range_in, V_wind, Print=False)[2:4]
    plt.close()
    return F_out, F_in

def power_output(x):
    gamma_in, gamma_out = x[1], x[2]
    a_elev_in = x[3]
    a_elev_out = x[4]
    A_proj = x[0]

    v_w_n = x[5]

    # F_out, F_in = aero_module(x)


    # print("gamma_in", gamma_in)
    # print("gamma_out", gamma_out)
    # print("a_elev_in", a_elev_in)
    # print("a_elev_out", a_elev_out)
    # print("A_proj", A_proj)
    # print("v_w_n", v_w_n)

    P_w = 0.5 * v_w_n ** 3 * rho  # Wind Power

    # TODO: do aerodynamic parameters vary with projected area?

    power = P_w * A_proj * (
            eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
            (F_in / eff_in) * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) * (
                    (gamma_out * gamma_in) / (gamma_out + gamma_in))

    return power

def objective(x):
    P_ref = baseline_power
    A_ref = baseline_area
    #return 1/np.sqrt(power_output(x)**2)
    #return -power_output(x)
    # # TODO: fix multi-objective
    return (P_ref / (1 + power_output(x))) + x[0] / A_ref