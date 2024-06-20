import numpy as np

# active constraints:
# constraint_power
# constraint_elev_out

from problem_definition.constraints import (
                                            constraint_reel_speed,
                                            constraint_elev_out,
                                            constraint_power,
                                            constraint_gamma_out,
                                            constraint_wind_speed,
                                            constraint_elev_in)

from problem_definition.objective import objective

constraint_list = [constraint_reel_speed,
                    constraint_elev_out,
                    constraint_power,
                    constraint_gamma_out,
                    constraint_wind_speed,
                    constraint_elev_in]

dG = np.zeros((6,6))
df = np.zeros((6,1))

upper_bounds = [1, 2.34, 0.29, 1.570796, 0.17453292519943295, 10.]

h = 1e-8

for constraint in range(len(constraint_list)):
    for design_variable in range(len(upper_bounds)):
        pertubation = upper_bounds.copy()
        pertubation[design_variable] = upper_bounds[design_variable] + h
        ffd = (constraint_list[constraint](pertubation) - constraint_list[constraint](upper_bounds)) / h

        dG[constraint, design_variable] = ffd

        df[design_variable] = (objective(pertubation) - objective(upper_bounds)) / h

print(np.dot(np.linalg.inv(dG), df))
