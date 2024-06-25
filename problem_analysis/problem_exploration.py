import matplotlib.pyplot as plt
import numpy as np
from problem_definition.objective import objective, power_output
from problem_definition.constraints import constraint_traction

import seaborn as sns
sns.set()

initial_design = np.zeros(6)

initial_design[0] = 10      # A_proj    0
initial_design[1] = 0.5     # gamma_in  1
initial_design[2] = 0.5     # gamma_out 2
initial_design[3] = 50 * np.pi / 180 # elev_in   3
initial_design[4] = 30 * np.pi / 180 # elev_out  4
initial_design[5] = 5       # v_w_n      5
# initial_design[6] = 5        # aoa_max_out
# initial_design[7] = 5        # aoa_max_in


bounds = [(1, 30), (0.1, 2.5), (0.1, 1), (0, 70 * np.pi/180), (0 * np.pi/180, 40 * np.pi/180), (0, 10)]#, (8, 16), (-1, 16)]

index = ['A_proj', 'gamma_in', 'gamma_out', 'elev_in', 'elev_out', 'v_w_n']

variable_range = np.zeros((6, 100))

variable_range[0] = np.linspace(bounds[0][0], bounds[0][1], 100)
variable_range[1] = np.linspace(bounds[1][0], bounds[1][1], 100)
variable_range[2] = np.linspace(bounds[2][0], bounds[2][1], 100)
variable_range[3] = np.linspace(bounds[3][0], bounds[3][1], 100)
variable_range[4] = np.linspace(bounds[4][0], bounds[4][1], 100)
variable_range[5] = np.linspace(bounds[5][0], bounds[5][1], 100)

gamma_plots = True
if gamma_plots:

    fig, ax = plt.subplots(2,3, figsize = (6.2,4))

    output_matrix = np.zeros((6,100))
    for design_variable in range(6):
        # Electrical power vs gamma_out for a fixed gamma_in

        design = initial_design.copy()

        output = []

        for x in variable_range[design_variable, :]:
            design[design_variable] = x

            if design_variable == 1 or design_variable == 3:
                traction = 1
            else:
                traction = 0

            #output.append(constraint_traction(design)[traction])
            output.append(objective(design))

        output_matrix[design_variable, :] = np.array(output)

    ax[0,0].plot(variable_range[0, :], output_matrix[0,:], label=f'{index[0]}')
    ax[0,0].set_xlabel(f'{index[0]}')
    ax[0,0].set_ylabel('Objective function \n value')

    ax[0,1].plot(variable_range[1, :], output_matrix[1,:], label=f'{index[1]}')
    ax[0,1].set_xlabel(f'{index[1]}')

    ax[0,2].plot(variable_range[2, :], output_matrix[2,:], label=f'{index[2]}')
    ax[0,2].set_xlabel(f'{index[2]}')

    ax[1,0].plot(variable_range[3, :], output_matrix[3,:], label=f'{index[3]}')
    ax[1,0].set_xlabel(f'{index[3]}')
    ax[1, 0].set_ylabel('Objective function \n value')

    ax[1,1].plot(variable_range[4, :], output_matrix[4,:], label=f'{index[4]}')
    ax[1,1].set_xlabel(f'{index[4]}')

    ax[1,2].plot(variable_range[5, :], output_matrix[5,:], label=f'{index[5]}')
    ax[1,2].set_xlabel(f'{index[5]}')

    plt.tight_layout()
    plt.show()
