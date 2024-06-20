import numpy as np
import pandas as pd
from scipy.optimize import minimize, NonlinearConstraint, LinearConstraint

from constants.constants import *

from problem_definition.objective import power_output, objective
from problem_definition.constraints import (constraint_traction,
                                            A, ub, lb,
                                            constraint_reel_speed,
                                            constraint_elev_out,
                                            constraint_power,
                                            constraint_gamma_out,
                                            constraint_area,
                                            constraint_wind_speed,
                                            constraint_elev_in)

index = ['A_proj', 'gamma_in', 'gamma_out', 'elev_in', 'elev_out', 'v_w_n']

# Initial guess
initial_design = np.zeros(6)

initial_design[0] = 10      # A_proj    0
initial_design[1] = 0.5     # gamma_in  1
initial_design[2] = 0.5     # gamma_out 2
initial_design[3] = 10 * np.pi / 180 # elev_in   3
initial_design[4] = 10 * np.pi / 180 # elev_out  4
initial_design[5] = 1       # v_w_n      5

bounds = [(0, 30), (0, np.inf), (0, 1), (0, np.pi/2), (10 * np.pi/180, 90 * np.pi/180), (0, 10)]

# TODO: make traction a linear constraint
cons = [NonlinearConstraint(constraint_traction, lb=(-np.inf, -np.inf), ub=(0, 0)),
        NonlinearConstraint(constraint_reel_speed, lb=-np.inf, ub=0),
        NonlinearConstraint(constraint_elev_out, lb=-0.85714, ub=0),
        NonlinearConstraint(constraint_power, lb=-np.inf, ub=0),
        NonlinearConstraint(constraint_gamma_out, lb=-np.inf, ub=0),
        NonlinearConstraint(constraint_area, lb=-np.inf, ub=0),
        NonlinearConstraint(constraint_wind_speed, lb=-np.inf, ub=0),
        NonlinearConstraint(constraint_elev_in, lb=-np.inf, ub=0)] # equality constraint, power output must always be positive
        #NonlinearConstraint(constraint_elev_out_lower, lb=-8, ub=0)]#,
        # # LinearConstraint(A, lb, ub, keep_feasible=True)]

method = "SLSQP" # unconstrained

# Perform the optimization
result = minimize(objective, initial_design, bounds=None, method=method,
                  options={'disp': True},
                  constraints=cons)

# power output for optimal design
optimal_power = power_output(result.x)

print(method)
print(pd.DataFrame(index=index, data=result.x, columns=['value']))

print(constraint_traction(result.x))
print(constraint_area(result.x))
print(constraint_reel_speed(result.x))
print(constraint_elev_out(result.x))
print(constraint_power(result.x))
print(constraint_wind_speed(result.x))
print(constraint_elev_in(result.x))
print(constraint_gamma_out(result.x))