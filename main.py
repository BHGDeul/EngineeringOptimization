import numpy as np
from scipy.optimize import minimize, NonlinearConstraint, LinearConstraint

from constants.constants import *

from problem_definition.objective import power_output, objective
from problem_definition.constraints import constraint_traction, constraint_reel_speed

# Initial guess
initial_design = np.zeros(6)

initial_design[0] = 10      # A_proj    0
initial_design[1] = 0.5     # gamma_in  1
initial_design[2] = 0.5     # gamma_out 2
initial_design[3] = 50 * np.pi / 180 # elev_in   3
initial_design[4] = 30 * np.pi / 180 # elev_out  4
initial_design[5] = 5       # v_w_n      5

bounds = [(0, 17), (0, 2.5), (0, 1), (0, 70 * np.pi/180), (0, np.pi/2), (0, 10)]

# TODO: make traction a linear constraint
cons = [NonlinearConstraintLinearConstraint(constraint_traction, lb=(-1, -1), ub=(0, 0)),
        NonlinearConstraint(constraint_reel_speed, -1, 0)]

method = "SLSQP"

# Perform the optimization
result = minimize(objective, initial_design, bounds=bounds, method=method,
                  options={'disp': True},
                  constraints=cons)

# power output for optimal design
optimal_power = power_output(result.x)

print(method)
print(result.x)