import numpy as np
from old_constants.input import *

# assume Power and Tether Limit are reached at the same time

def calculate_power_limit_reel_speeds(step):
    F_in = 0.07
    F_out = 5.4

    # call optimal nominal gamma
    # gamma_out_n,gamma_in_n = calculate_opt_gamma_nominal() #nominal reel out from optimal_gamma.py
    gamma_out_n = 0.3
    gamma_in_n = 1.65

    # create lists to have running variables
    v_w = np.linspace(v_w_n, 2.5 * v_w_n, step)
    gamma_in = np.linspace(gamma_in_n, 0.25, step)
    gamma_in_max_f_c = np.zeros(step)
    gamma_out_max_f_c = np.zeros(step)
    f_c_mu = np.zeros(step)
    P_c = np.zeros(step)

    # create counters and set to zero
    ci = 0
    cj = 0

    # loop to go through all the different gamma_in for all v_w
    # also calculates gamma_out
    for i in v_w:
        mu = i / v_w_n
        P_w = 0.5 * rho * i ** 3
        F_out_mu = (F_out / (mu ** 2)) * (((1 - gamma_out_n) ** 2) / ((1 - (gamma_out_n / mu)) ** 2))
        for j in gamma_in:
            f_c_mu[cj] = ((1 / (mu ** 2)) * (1 - gamma_out_n) ** 2 - (F_in / F_out) * (1 + j) ** 2) * (
                        (gamma_out_n * j) / (gamma_out_n + mu * j))
            cj += 1
        max_f_c = np.amax(f_c_mu)
        a = np.where(f_c_mu == max_f_c)
        gamma_in_max_f_c[ci] = gamma_in[a]
        gamma_out_max_f_c[ci] = gamma_out_n / mu
        P_c[ci] = P_w * max_f_c
        ci += 1
        cj = 0
    # print ("gamma_in=",gamma_in_max_f_c,"gamma_out",gamma_out_max_f_c,"wind_speed=",v_w)
    return gamma_in_max_f_c, gamma_out_max_f_c, v_w, P_c


def calculate_power_before_v_w_n(step):
    gamma_out_n, gamma_in_n = calculate_opt_gamma_nominal()

    F_in = 0.07
    F_out = 5.4

    v_w_before_v_n = np.linspace(0, v_w_n, step)
    power_before_v_w_n = np.zeros(step)

    ci = 0
    for i in v_w_before_v_n:
        P_w = 0.5 * rho * i ** 3
        power_before_v_w_n[ci] = (1 / F_out) * P_w * (
                    F_out * (1 - gamma_out_n) ** 2 - (F_in * (1 + gamma_in_n) ** 2)) * (
                                             (gamma_out_n * gamma_in_n) / (gamma_out_n + gamma_in_n))
        ci += 1
    return power_before_v_w_n, v_w_before_v_n
