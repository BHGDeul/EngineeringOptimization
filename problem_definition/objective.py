from constants.constants import rho, F_in, F_out, eff_out, eff_in, baseline_power, baseline_area
import numpy as np

#bounds = [(0, 30), (0, np.inf), (0, 1), (0, np.pi/2), (10 * np.pi/180, 90 * np.pi/180), (0, 10)]
upper_bounds = [1, 2.34, 0.29, 1.570796, 0.17453292519943295, 10.]

def power_output(x):
    gamma_in, gamma_out = x[1], x[2]
    a_elev_in = x[3]
    a_elev_out = x[4]
    A_proj = x[0]
    v_w_n = x[5]

    P_w = 0.5 * v_w_n ** 3 * rho  # Wind Power

    # TODO: do aerodynamic parameters vary with projected area?

    power = P_w * A_proj * (
            eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
            (F_in * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) / eff_in) * (
                    (gamma_out * gamma_in) / (gamma_out + gamma_in))

    return power

def objective(x):
    P_ref = baseline_power
    A_ref = baseline_area
    #return 1/np.sqrt(power_output(x)**2)
    #return -power_output(x)
    # # TODO: fix multi-objective
    return (P_ref / (1 + power_output(x))) + x[0] / A_ref