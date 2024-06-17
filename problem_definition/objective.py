from constants.constants import *

def power_output(x):
    gamma_in, gamma_out = x[1], x[2]
    a_elev_in = x[3]
    a_elev_out = x[4]
    A_proj = x[0]

    P_w = 0.5 * x[5] ** 3 * data['rho']  # Wind Power

    # TODO: do aerodynamic parameters vary with projected area?

    power = P_w * A_proj * (
        eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
        (F_in * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) / eff_in) * (
        (gamma_out * gamma_in) / (gamma_out + gamma_in))

    return power

def objective(x):
    P_ref = baseline_power
    A_ref = baseline_area
    return -1 * power_output(x) / P_ref + x[0] / A_ref