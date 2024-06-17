import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from constants.constants import *

# Intermediate calculations
F_out = CL_out ** 3 / CD_out ** 2
F_in = CD_in
P_w = 0.5 * v_w_n ** 3 * rho  # Wind Power

# Define gamma range
lim = max_reel_speed / v_w_n
resolution = 100
gamma_in_range = np.linspace(0.01, lim, resolution)
gamma_out_range = np.linspace(0.01, 1, resolution)

# Set empty arrays
power_array_e = np.zeros((resolution, resolution))

# Direct Calculation
for cj, gamma_out in enumerate(gamma_out_range):
    for ci, gamma_in in enumerate(gamma_in_range):
        power_array_e[ci][cj] = P_w * A_proj * (
            eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
            (F_in * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) / eff_in) * (
            (gamma_out * gamma_in) / (gamma_out + gamma_in))

# Find maximal electrical power with Brute Force Method
P_elec_opt_gamma = np.amax(power_array_e)
(a, b) = np.where(power_array_e == P_elec_opt_gamma)

gamma_out_n = gamma_out_range[b][0]
gamma_in_n = gamma_in_range[a][0]

print(gamma_out_n, gamma_in_n, P_elec_opt_gamma)

# Optimization
# def objective(gamma):
#     gamma_out, gamma_in = gamma
#     power = P_w * A_proj * (
#         eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
#         (F_in * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) / eff_in) * (
#         (gamma_out * gamma_in) / (gamma_out + gamma_in))
#     return -power  # Negative because we want to maximize the power

# Define bounds for gamma_out and gamma_in
bounds = [(0.01, 1), (0.01, lim)]
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

# Plot electrical power array for gammas
plt.figure(figsize=(12, 10))
contour = plt.contourf(gamma_out_range, gamma_in_range, power_array_e, levels=200, cmap='viridis')
plt.colorbar(contour, label='Electrical Power (W)')
plt.xlabel('Gamma Out')
plt.ylabel('Gamma In')
plt.title('Electrical Power Distribution')

# Maximum point
plt.scatter(gamma_out_n, gamma_in_n, color='red', label='Max Power', zorder=5)

# Annotate the points
plt.annotate(f'Max\n({gamma_out_n:.2f}, {gamma_in_n:.2f})', xy=(gamma_out_n, gamma_in_n),
             xytext=(gamma_out_n + 0.1, gamma_in_n + 0.1),
             arrowprops=dict(facecolor='red', shrink=0.05))
plt.legend()

plt.show()

# Define a function to calculate electrical power
def electrical_power(gamma_out, gamma_in):
    power = P_w * A_proj * (
        eff_out * F_out * (1 - gamma_out) ** 2 -
        (F_in * (1 + gamma_in) ** 2) / eff_in) * ((gamma_out * gamma_in) / (gamma_out + gamma_in))
    return power

gamma_plots = False
if gamma_plots:
    # Gamma ranges for plotting
    gamma_out_range = np.linspace(0.01, 1, 100)
    gamma_in_range = np.linspace(0.01, lim, 100)

    # Electrical power vs gamma_out for a fixed gamma_in
    plt.figure(figsize=(6, 4))
    gamma_in_fixed = 0.5
    power_values_out = [electrical_power(g, gamma_in_fixed) for g in gamma_out_range]
    plt.plot(gamma_out_range, power_values_out, label=f'Gamma_in = {gamma_in_fixed}')
    plt.xlabel('Gamma Out')
    plt.ylabel('Electrical Power (W)')
    plt.title('Electrical Power vs Gamma Out for fixed Gamma In')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Electrical power vs gamma_in for a fixed gamma_out
    plt.figure(figsize=(12, 6))
    gamma_out_fixed = 0.5
    power_values_in = [electrical_power(gamma_out_fixed, g) for g in gamma_in_range]
    plt.plot(gamma_in_range, power_values_in, label=f'Gamma_out = {gamma_out_fixed}')
    plt.xlabel('Gamma In')
    plt.ylabel('Electrical Power (W)')
    plt.title('Electrical Power vs Gamma In for fixed Gamma Out')
    plt.legend()
    plt.grid(True)
    plt.show()
