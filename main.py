import numpy as np
from scipy.optimize import minimize

from constants.constants import *


# Vector of design variables
x = np.zeros(6)
# A_proj    0
# gamma_in  1
# gamma_out 2
# elev_in   3
# elev_out  4
# v_wn      5

# Initial guess
initial_design = x

initial_design[0] = 15
initial_design[1] = 0.5
initial_design[2] = 0.5
initial_design[3] = 30 * np.pi / 180
initial_design[4] = 70 * np.pi / 180
initial_design[5] = 5

bounds = [(0, 17), (0, 1), (0, 1), (0, np.pi/2), (0, np.pi/2), (0, 10)]

# Initial guess
initial_guess = [0.5, 0.5]
# Perform the optimization
result = minimize(objective, initial_guess, bounds=bounds, method='L-BFGS-B')

# Extract the optimal values
optimal_gamma_out, optimal_gamma_in = result.x
max_electrical_power = -result.fun

# Update the optimal values
gamma_out_n = optimal_gamma_out
gamma_in_n = optimal_gamma_in
P_elec_opt_gamma = max_electrical_power
print(max_electrical_power)



# Define a function to calculate electrical power
def electrical_power(gamma_out, gamma_in):
    power = P_w * A_proj * (
        eff_out * F_out * (1 - gamma_out) ** 2 -
        (F_in * (1 + gamma_in) ** 2) / eff_in) * ((gamma_out * gamma_in) / (gamma_out + gamma_in))
    return power
