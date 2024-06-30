import pandas as pd
from scipy.optimize import minimize, NonlinearConstraint

from constants.constants import *

from problem_definition.objective import power_output, objective, aero_module
from problem_definition.constraints import (constraint_traction,
                                            constraint_reel_speed,
                                            constraint_power)

from plotting.monte_carlo import create_boxplot

index = ['A_proj',
         'gamma_in',
         'gamma_out',
         'elev_in',
         'elev_out',
         'v_w_n',
         #, 'aoa_in'] uncomment for aeromodule implementation
]

import random

# Initial guess
initial_design = np.zeros(6)

initial_design[0] = 10      # A_proj    0
initial_design[1] = 1.9     # gamma_in  1
initial_design[2] = 0.2     # gamma_out 2
initial_design[3] = 50 * np.pi / 180 # elev_in   3
initial_design[4] = 30 * np.pi / 180 # elev_out  4
initial_design[5] = 10       # v_w_n      5

## uncomment for aeromodule implementation
# initial_design[6] = -0.4        # aoa_in


bounds = [(0.01, 30),
          (0.1, 2.5),
          (0.1, 1),
          (0, 70 * np.pi/180),
          (30 * np.pi/180, 40 * np.pi/180),
          (0, 10),
          # (-1,5) # uncomment for aeromodule implementation
          ]

# List of non-linear constraints
cons = [NonlinearConstraint(constraint_traction, lb=(-np.inf, -np.inf), ub=(0, 0)),
        NonlinearConstraint(constraint_reel_speed, lb=-np.inf, ub=0),
        NonlinearConstraint(constraint_power, lb=-1., ub=0)]

# Scipy.optimize.minimize method selected
method = "COBYLA"

# number of optimization cycles to conduct
cycles = 1

# matrices below only required for monte-carlo analysis
store = np.zeros((6, cycles))
objective_values = np.zeros((cycles,))

for cycle in range(cycles):
    ## Uncomment below when performing monte-carlo analysis.
    # initial_design[0] = random.randint(1,30)  # A_proj    0
    # initial_design[1] = random.randint(1, 25) / 10  # gamma_in  1
    # initial_design[2] = random.randint(1, 10) / 10  # gamma_out 2
    # initial_design[3] = random.randint(0, 70) * np.pi / 180  # elev_in   3
    # initial_design[4] = random.randint(0, 90) * np.pi / 180  # elev_in   3

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

## uncomment line 86-90 when performing monte-carlo analysis.
#     print(cycle)
#     store[:, cycle] = result.x
#     objective_values[cycle] = objective(result.x)
#
# create_boxplot(store, objective_values)


## print statements to verify compliance with constraints.
print('Traction constraint', constraint_traction(result.x))
print('reel speed constraint', constraint_reel_speed(result.x))
print('power constraint', constraint_power(result.x))