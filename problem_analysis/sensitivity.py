import numpy as np
from constants.constants import rho, eff_out, eff_in, F_out, F_in
from problem_definition.objective import power_output

# consider power_output the response of the system

#minimizer = np.array([22.44414079,  2.34574528,  0.29021158,  1.57079633,  0.17453642, 10.])
from main import initial_design as minimizer


pertubations = np.logspace(-15,-1,30)

sensitivities = np.zeros((len(minimizer), len(pertubations)))
log_sensitivities = np.zeros((len(minimizer), len(pertubations)))

j = 0
for pertubation in pertubations:
    i = 0
    for design_variable in range(len(minimizer)):

        perturbed_design = minimizer.copy()
        perturbed_design[design_variable] = minimizer[design_variable] + pertubation

        ffd = (power_output(perturbed_design) - power_output(minimizer)) / pertubation
        sensitivities[i,j] = ffd


        # print log sensitivity
        log_sens = minimizer[design_variable] / power_output(minimizer) * ffd
        print(log_sens)
        log_sensitivities[i, j] = log_sens

        i+=1
    print('\n')
    j += 1

# analytical sensitivities
x = minimizer.copy()
P_w = 0.5 * x[5] ** 3 * rho
b = P_w * x[0]
c = F_in * (x[1]**2 + 2 * np.cos(x[3]) * x[1] + 1) / eff_in
d = F_out * eff_out
e = x[1]

dudx2 = (
    e * b * (
    2 * d * x[2]**3 + 3 * e * d * x[2]**2 + e * d * (np.cos(x[4]))**2 - 2 * d * (x[2])**2 * np.cos(x[4]) \
    - 4 * e * d * x[2] * np.cos(x[4]) - e * c
)/ (x[2] + x[1])**2
)

index = ['A_proj', 'gamma_in', 'gamma_out', 'elev_in', 'elev_out', 'v_w_n']

import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,2, figsize = (8,4))

design_variable = 2
ax[0].plot(pertubations, sensitivities[design_variable], label='ffd')
ax[0].set_xscale('log')
# plt.xscale('log')
ax[0].set_title(f'Comparison of FFD approximation\n and analytical value of sensitivity \n for {index[design_variable]}')
ax[0].hlines(y=dudx2, xmin=pertubations[0],xmax=pertubations[-1], colors='red', label='analytical')
ax[0].set_xlabel('Pertubation size')
ax[0].set_ylabel(f'du/d{index[design_variable]}')

for design_variable in range(6):
    ax[1].plot(pertubations, sensitivities[design_variable], label=f'{index[design_variable]}')
    ax[1].set_xscale('log')
    ax[1].set_title(f'FFD sensitivity of the power output for \n different pertubations in design_variables')
    ax[1].set_ylabel(f'du/dx_i')
    ax[1].set_xlabel('Pertubation size')

plt.legend()
plt.tight_layout()
plt.grid()
plt.savefig(f'../output/plots/sensitivity_gamma_out.png')
plt.show()
# break

