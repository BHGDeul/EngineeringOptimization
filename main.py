import pandas as pd
from scipy.optimize import minimize, NonlinearConstraint

from constants.constants import *

from problem_definition.objective import power_output, objective
from problem_definition.constraints import (constraint_traction,
                                            constraint_reel_speed,
                                            constraint_power)

from plotting.monte_carlo import create_boxplot

index = ['A_proj', 'gamma_in', 'gamma_out', 'elev_in', 'elev_out', 'v_w_n']

import random

# Initial guess
initial_design = np.zeros(6)

initial_design[0] = 10      # A_proj    0
initial_design[1] = 0.5     # gamma_in  1
initial_design[2] = 0.5     # gamma_out 2
initial_design[3] = 50 * np.pi / 180 # elev_in   3
initial_design[4] = 30 * np.pi / 180 # elev_out  4
initial_design[5] = 10       # v_w_n      5
# initial_design[6] = 5        # aoa_max_out
# initial_design[7] = 5        # aoa_max_in


bounds = [(1, 30), (0.1, np.inf), (0.1, 1), (0, 70 * np.pi/180), (30 * np.pi/180, 40 * np.pi/180), (0, 10)]#, (8, 16), (-1, 16)]

# TODO: make traction a linear constraint
cons = [NonlinearConstraint(constraint_traction, lb=(-np.inf, -np.inf), ub=(0, 0)),
        NonlinearConstraint(constraint_reel_speed, lb=-np.inf, ub=0),
        NonlinearConstraint(constraint_power, lb=-np.inf, ub=0)]

method = "COBYLA" # unconstrained

cycles = 1
store = np.zeros((6, cycles))
objective_values = np.zeros((cycles,))

for cycle in range(cycles):
    # initial_design[0] = random.randint(1,30)  # A_proj    0
    # initial_design[1] = random.randint(1, 25) / 10  # gamma_in  1
    # initial_design[2] = random.randint(1, 10) / 10  # gamma_out 2
    # initial_design[3] = random.randint(0, 90) * np.pi / 180  # elev_in   3

    design = initial_design.copy()

    # Perform the optimization
    result = minimize(objective, design, bounds=bounds, method=method,
                      options={'disp': True},
                      constraints=cons,
                      )

    # power output for optimal design
    optimal_power = power_output(result.x)

    print(method)
    print(pd.DataFrame(index=index, data=result.x, columns=['value']))
    print('\nPower output', round(power_output((result.x)),2), 'Watt')
    # print(cycle)
    # store[:, cycle] = result.x
    # objective_values[cycle] = objective(result.x)
#
# create_boxplot(store, objective_values)


#
print('Traction constraint', constraint_traction(result.x))
print('reel speed constraint', constraint_reel_speed(result.x))
# print('elevation out constraint', constraint_elev_out(result.x))
print('power constraint', constraint_power(result.x))
# print('wind speed constraint', constraint_wind_speed(result.x))
# print('elevation in constraint', constraint_elev_in(result.x))
# print('gamma out constraint', constraint_gamma_out(result.x))