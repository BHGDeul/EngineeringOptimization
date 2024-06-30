import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, NonlinearConstraint

# data initialization
v_w_n = 10
CL_out = 1.06
A_proj = 16.65
rho = 1.18
CD_out = 0.15
CD_in = 0.10
eff_in = 0.639
eff_out = 0.652
max_reel_speed = 25
a_elev_out = 30 * np.pi / 180
a_elev_in = 70 * np.pi / 180
T_max_in = 10300  # for visual purposes, real value is also close to 10300
T_max_out = 10300

# intermediate calculations
F_out = CL_out ** 3 / CD_out ** 2
F_in = CD_in
P_w = 0.5 * v_w_n ** 3 * rho  # wind power

# gamma range
lim = 5  # max_reel_speed / v_w_n
resolution = 100
gamma_in_range = np.linspace(0.01, lim, resolution)
gamma_out_range = np.linspace(0.01, 1, resolution)

# initialize empty arrays
power_array_e = np.zeros((resolution, resolution))
t_in_array = np.zeros((resolution, resolution))
t_out_array = np.zeros((resolution, resolution))

# calculate power and boundaries
for cj, gamma_out in enumerate(gamma_out_range):
    for ci, gamma_in in enumerate(gamma_in_range):
        power_array_e[ci][cj] = P_w * A_proj * (
                eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
                (F_in * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) / eff_in) * (
                                        (gamma_out * gamma_in) / (gamma_out + gamma_in))

        t_out_array[ci][cj] = 0.5 * rho * v_w_n ** 2 * A_proj * (np.cos(a_elev_out) - gamma_out) ** 2 * F_out
        t_in_array[ci][cj] = 0.5 * rho * v_w_n ** 2 * A_proj * (1 + 2 * gamma_in * np.cos(a_elev_in) + gamma_in ** 2) * F_in

# optimization
def objective(gamma):
    gamma_out, gamma_in = gamma
    power = P_w * A_proj * (
            eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
            (F_in * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) / eff_in) * (
                    (gamma_out * gamma_in) / (gamma_out + gamma_in))
    return -power  # negative because we want to maximize the power

# bounds for gamma_out and gamma_in
bounds = [(0.01, 1), (0.01, lim)]
# initial guess
initial_guess = [0.5, 0.5]

# perform the optimization
result = minimize(objective, initial_guess, bounds=bounds, method='SLSQP')

def plot_result(result):
    # extract optimal values
    optimal_gamma_out, optimal_gamma_in = result.x
    max_electrical_power = -result.fun

    # electrical power array for gammas
    plt.figure(figsize=(6, 4))
    contour = plt.contourf(gamma_out_range, gamma_in_range, power_array_e, levels=200, cmap='viridis')
    plt.colorbar(contour, label='electrical power (W)')
    plt.xlabel('Gamma Out')
    plt.ylabel('Gamma In')
    plt.title('Electrical Power Distribution')

    # max
    plt.scatter(optimal_gamma_out, optimal_gamma_in, color='red', label='Max Power', zorder=5)

    # annotate
    plt.annotate(f'Max\n({optimal_gamma_out:.2f}, {optimal_gamma_in:.2f}), {round(max_electrical_power)} W', xy=(optimal_gamma_out, optimal_gamma_in),
                 xytext=(optimal_gamma_out + 0.1, optimal_gamma_in + 0.1),
                 arrowprops=dict(facecolor='red', shrink=0.05))

    # max 2.5
    plt.axhline(y=2.5, color='red', linestyle='--', linewidth=2)
    plt.text(gamma_out_range[-1] - 0.1, 2.8, 'Unfeasible Region', color='red', fontsize=20, ha='right')

    plt.legend()
    plt.show()

plot_result(result)

def plot_result_with_boundaries(result):
    # extract optimal values
    optimal_gamma_out, optimal_gamma_in = result.x
    max_electrical_power = -result.fun

    # electrical power array for gammas
    plt.figure(figsize=(6, 4))
    contour = plt.contourf(gamma_out_range, gamma_in_range, power_array_e, levels=200, cmap='viridis')
    plt.colorbar(contour, label='electrical power (W)')
    plt.xlabel('Gamma Out')
    plt.ylabel('Gamma In')
    plt.title('Electrical Power Distribution with Boundaries')

    # max
    plt.scatter(optimal_gamma_out, optimal_gamma_in, color='red', label='Max Power', zorder=5)

    # annotate
    plt.annotate(f'Max\n({optimal_gamma_out:.2f}, {optimal_gamma_in:.2f}), {round(max_electrical_power)} W', xy=(optimal_gamma_out, optimal_gamma_in),
                 xytext=(optimal_gamma_out + 0.1, optimal_gamma_in + 0.1),
                 arrowprops=dict(facecolor='red', shrink=0.05))

    # boundary lines
    plt.contourf(gamma_out_range, gamma_in_range, t_in_array, levels=[T_max_in, max(t_in_array.max(),T_max_in+1)], colors='red',
                 alpha=0.6, hatches=['///'])
    plt.contourf(gamma_out_range, gamma_in_range, t_out_array, levels=[T_max_out, max(t_out_array.max(), T_max_out+1)], colors='blue',
                 alpha=0.3, hatches=['\\\\'])

    plt.legend()
    plt.show()

plot_result_with_boundaries(result)

# constraints
def constraint_traction(gamma):
    gamma_out, gamma_in = gamma
    T_out_n = 0.5 * rho * v_w_n ** 2 * A_proj * (np.cos(a_elev_out) - gamma_out) ** 2 * F_out
    T_in_n = 0.5 * rho * v_w_n ** 2 * A_proj * (1 + 2 * gamma_in * np.cos(a_elev_in) + gamma_in ** 2) * F_in
    return T_out_n / T_max_out - 1, T_in_n / T_max_in - 1

def constraint_reel_speed(gamma):
    gamma_out, gamma_in = gamma
    return gamma_out, gamma_in / (max_reel_speed / v_w_n) - 1

traction_constraint = NonlinearConstraint(constraint_traction, [-1, -1], [0, 0])
reel_speed_constraint = NonlinearConstraint(constraint_reel_speed, [0, -1], [1, 0])

# initial guess
initial_guess = [0.5, 0.5]

# perform the optimization with constraints
result = minimize(objective, initial_guess, bounds=None, method='SLSQP',
                  constraints=[reel_speed_constraint, traction_constraint],
                  options={'disp': True})
plot_result_with_boundaries(result)

# define a function to calculate electrical power
def electrical_power(gamma_out, gamma_in):
    power = P_w * A_proj * (
            eff_out * F_out * (np.cos(a_elev_out) - gamma_out) ** 2 -
            (F_in * (gamma_in ** 2 + 2 * np.cos(a_elev_in) * gamma_in + 1)) / eff_in) * (
                    (gamma_out * gamma_in) / (gamma_out + gamma_in))
    return power

gamma_plots = True
if gamma_plots:
    # gamma ranges for plotting
    gamma_out_range = np.linspace(0.01, 1, 100)
    gamma_in_range = np.linspace(0.01, lim, 100)

    # electrical power vs gamma_out for a fixed gamma_in
    plt.figure(figsize=(6, 4))
    gamma_in_fixed = 1.94
    power_values_out = [electrical_power(g, gamma_in_fixed) for g in gamma_out_range]
    plt.plot(gamma_out_range, power_values_out, label=f'Gamma_in = {gamma_in_fixed}')
    plt.xlabel('Gamma Out')
    plt.ylabel('Electrical Power (W)')
    plt.title('Electrical Power vs Gamma Out for fixed Gamma In')
    plt.legend()
    plt.grid(True)
    plt.show()

    # electrical power vs gamma_in for a fixed gamma_out
    plt.figure(figsize=(12, 6))
    gamma_out_fixed = 0.25
    power_values_in = [electrical_power(gamma_out_fixed, g) for g in gamma_in_range]
    plt.plot(gamma_in_range, power_values_in, label=f'Gamma_out = {gamma_out_fixed}')
    plt.xlabel('Gamma In')
    plt.ylabel('Electrical Power (W)')
    plt.title('Electrical Power vs Gamma In for fixed Gamma Out')
    plt.legend()
    plt.grid(True)
    plt.show()

    # traction vs gamma_in for a fixed gamma_out
    plt.figure(figsize=(12, 6))
    plt.plot(gamma_out_range, constraint_traction((gamma_out_range, gamma_in_range))[0])
    plt.xlabel('Gamma Out')
    plt.ylabel('Traction Out')
    plt.title('Traction Out vs Gamma Out')
    plt.grid(True)
    plt.show()

    # traction vs gamma_in for a fixed gamma_out
    plt.figure(figsize=(12, 6))
    plt.plot(gamma_in_range, constraint_traction((gamma_out_range, gamma_in_range))[1])
    plt.xlabel('Gamma In')
    plt.ylabel('Traction In')
    plt.title('Traction In vs Gamma In')
    plt.grid(True)
    plt.show()

    # reel speeds
    plt.figure(figsize=(12, 6))
    plt.plot(gamma_in_range, constraint_reel_speed((gamma_out_range, gamma_in_range))[1])
    plt.xlabel('Gamma In')
    plt.ylabel('Reel Speed Constraint')
    plt.grid(True)
    plt.show()
