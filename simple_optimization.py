from scipy.optimize import minimize, NonlinearConstraint

from problem_definition.objective import objective
from problem_definition.constraints import *

from constants import data

import numpy as np

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
initial_design[5] = 10

bounds = [(0, 17), (0, 1), (0, 1), (0, np.pi/2), (0, np.pi/2), (0, 10)]

traction_constraint = NonlinearConstraint(constraint_traction, [0, 0], [10E3, 10E3]) # This gives an inequality constraint
reel_speed_constraint = NonlinearConstraint(constraint_reel_speed, [-1], [0]) # This gives an inequality constraint

constraints = {traction_constraint, reel_speed_constraint}

method = 'SLSQP'

"""
‘Nelder-Mead’ (see here)
‘Powell’ (see here)
‘CG’ (see here)
‘BFGS’ (see here)
‘Newton-CG’ (see here)
‘L-BFGS-B’ (see here)
‘TNC’ (see here)
‘COBYLA’ (see here)
‘SLSQP’ (see here)
‘trust-constr’(see here)
‘dogleg’ (see here)
‘trust-ncg’ (see here)
‘trust-exact’ (see here)
‘trust-krylov’ (see here)
"""

# Perform the optimization
result = minimize(objective, initial_design, bounds=bounds, method=method,
                  constraints=constraints,#{traction_constraint, reel_speed_constraint},
                  options={'disp': True})

# Extract the optimal values
from problem_definition.objective import power_output
max_electrical_power = power_output(result.x)

if __name__ == '__main__':
    print(result.message)
    print(result.success)
    print('optimal gamma_in, gamma_out', result.x)
    print('Optimal Power output', max_electrical_power)

    output_file = f'output/{method}_minimize_output.txt'
    with open(output_file, 'w') as f:
        f.write("Optimization Result:\n")
        f.write(f"Success: {result.success}\n")
        f.write(f"Status: {result.status}\n")
        f.write(f"Message: {result.message}\n")
        f.write(f"Number of Iterations: {result.nit}\n")
        f.write(f"Function Value: {result.fun}\n")
        f.write(f"Optimal Parameters: {result.x}\n")


    plotting = False

    if plotting:
        import numpy as np
        import matplotlib.pyplot as plt

        # Define a function to calculate electrical power
        def electrical_power(gamma_out, gamma_in, data):
            power = data['P_w'] * data['A_proj'] * (
                data['eff_out'] * data['F_out'] * (1 - gamma_out) ** 2 -
                (data['F_in'] * (1 + gamma_in) ** 2) / data['eff_in']) * ((gamma_out * gamma_in) / (gamma_out + gamma_in))
            return power

        # gamma ranges for plotting
        gamma_out_range = np.linspace(0.01, 1, 100)
        gamma_in_range = np.linspace(0.01, lim, 100)

        # electrical power vs gamma_out for a fixed gamma_in
        plt.figure(figsize=(6, 6))
        gamma_in_fixed = 0.5
        power_values_out = [electrical_power(g, gamma_in_fixed, data) for g in gamma_out_range]
        plt.plot(gamma_out_range, power_values_out, label=f'Gamma_in = {gamma_in_fixed}')
        plt.xlabel('Gamma Out')
        plt.ylabel('Electrical Power (W)')
        plt.title('Electrical Power vs Gamma Out for fixed Gamma In')
        plt.legend()
        plt.grid(True)
        plt.show()


        #electrical power vs gamma_in for a fixed gamma_out
        plt.figure(figsize=(6, 6))
        gamma_out_fixed = 0.5
        power_values_in = [electrical_power(gamma_out_fixed, g, data) for g in gamma_in_range]
        plt.plot(gamma_in_range, power_values_in, label=f'Gamma_out = {gamma_out_fixed}')
        plt.xlabel('Gamma In')
        plt.ylabel('Electrical Power (W)')
        plt.title('Electrical Power vs Gamma In for fixed Gamma Out')
        plt.legend()
        plt.grid(True)
        plt.show()

        #electrical power vs gamma_in for a fixed gamma_out
        plt.figure(figsize=(6, 6))
        # power_values_in = [constraint_traction((gamma_in_range,gamma_out_range), g, data) for g in gamma_in_range]
        plt.plot(gamma_out_range, constraint_traction((gamma_in_range,gamma_out_range))[0])
        plt.xlabel('Gamma Out')
        plt.ylabel('Traction Out')
        plt.title('Traction Out vs Gamma Out')
        plt.grid(True)
        plt.show()

        #electrical power vs gamma_in for a fixed gamma_out
        plt.figure(figsize=(6, 6))
        # power_values_in = [constraint_traction((gamma_in_range,gamma_out_range), g, data) for g in gamma_in_range]
        plt.plot(gamma_in_range, constraint_traction((gamma_in_range,gamma_out_range))[1])
        plt.xlabel('Gamma In')
        plt.ylabel('Traction In')
        plt.title('Traction In vs Gamma In')
        plt.grid(True)
        plt.show()

        #reel speeds
        plt.figure(figsize=(6, 6))
        # power_values_in = [constraint_traction((gamma_in_range,gamma_out_range), g, data) for g in gamma_in_range]
        plt.plot(gamma_in_range, constraint_reel_speed((gamma_in_range,gamma_out_range))[1])
        plt.xlabel('Gamma In')
        plt.ylabel('Reel speed constraint')
        plt.title('Traction In vs Gamma In')
        plt.grid(True)
        plt.show()