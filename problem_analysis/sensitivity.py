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

import matplotlib.pyplot as plt
plt.plot(pertubations, sensitivities[2], color='blue', label='ffd')
plt.xscale('log')
plt.title('FFD sensitivity of the power output for \n different pertubations in gamma_out')
plt.hlines(y=dudx2, xmin=pertubations[0],xmax=pertubations[-1], colors='red', label='analytical')
plt.xlabel('Pertubation size')
plt.ylabel('du/dgamma_out')
plt.legend()
plt.grid()
plt.savefig('../output/plots/sensitivity_dudx2.png')
plt.show()