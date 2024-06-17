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
initial_design[5] = 5

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


def callbackF(Xi):
    print(f"{list(np.around(Xi, 5))} \t\t {objective(Xi)}")


# Perform the optimization
result = minimize(objective, initial_design, bounds=bounds, method=method,
                  constraints=constraints,#{traction_constraint, reel_speed_constraint},
                  options={'disp': True},
                  callback=callbackF)

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