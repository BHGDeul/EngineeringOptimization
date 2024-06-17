from constants import data

def power_output(x):
    gamma_in, gamma_out = x[1], x[2]
    A_proj = x[0]

    P_w = 0.5 * x[5] ** 3 * data['rho']  # Wind Power

    # TODO: do aerodynamic parameters vary with projected area?

    power = P_w * A_proj * (
        data['eff_out'] * data['F_out'] * (1 - gamma_out) ** 2 -
        (data['F_in'] * (1 + gamma_in) ** 2) / data['eff_in']) * ((gamma_out * gamma_in) / (gamma_out + gamma_in))
    return power

def objective(x):
    P_ref = data['baseline_power']
    A_ref = data['baseline_area']
    return -1 * power_output(x) / P_ref + x[0] / A_ref