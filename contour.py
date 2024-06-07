import matplotlib.pyplot as plt
import numpy as np


def contour_plot():
    # electrical power array for gammas
    plt.figure(figsize=(12, 10))
    contour = plt.contourf(plot_gamma_data['gamma_out'], plot_gamma_data['gamma_in'], plot_gamma_data['power_array_e'],
                           levels=200, cmap='viridis')
    plt.colorbar(contour, label='Electrical Power (W)')
    plt.xlabel('Gamma Out')
    plt.ylabel('Gamma In')
    plt.title('Electrical Power Distribution')

    # maximum point
    plt.scatter(data['gamma_out_n'], data['gamma_in_n'], color='red', label='Max Power', zorder=5)

    # Annotate the points
    plt.annotate(f'Max\n({data["gamma_out_n"]:.2f}, {data["gamma_in_n"]:.2f})',
                 xy=(data['gamma_out_n'], data['gamma_in_n']),
                 xytext=(data['gamma_out_n'] + 0.1, data['gamma_in_n'] + 0.1),
                 arrowprops=dict(facecolor='red', shrink=0.05))

    plt.legend()

    plt.show()

plot_gamma_data = {}
# Define gamma range ##
"""Prohibits reel-in speed from exceeding max reeling speed"""
# should also be a cosntraint
lim = data['max_reel_speed'] / data['v_w_n']
# if data['max_reel_speed'] <= 2 * data['v_w_n']:
#     lim = data['max_reel_speed'] / data['v_w_n']
# else:
#     lim = 2.0
resolution = 100
plot_gamma_data['gamma_in'] = np.linspace(0.01, lim, resolution)
plot_gamma_data['gamma_out'] = np.linspace(0.01, 1, resolution)


""""Original: Brute force method"""
## Set empty arrays ##
plot_gamma_data['power_array_m'] = np.zeros((resolution, resolution))
plot_gamma_data['power_array_e'] = np.zeros((resolution, resolution))
## Initiate counters ##
ci = 0
cj = 0
for cj, gamma_out in enumerate(plot_gamma_data['gamma_out'] ):
    for ci, gamma_in in enumerate(plot_gamma_data['gamma_in'] ):
        # plot_gamma_data['power_array_m'][cj][ci] = data['P_w'] * data['A_proj'] * (
        #         data['F_out'] * (1 - j) ** 2 - (data['F_in'] * (1 + i) ** 2)) * ((j * i) / (j + i))

        plot_gamma_data['power_array_e'][ci][cj] = data['P_w'] * data['A_proj'] * (
                data['eff_out'] * data['F_out'] * (1 - gamma_out) ** 2 - (data['F_in'] * \
                (1 + gamma_in) ** 2) / data['eff_in']) * ((gamma_out * gamma_in) / (gamma_out + gamma_in))

        ci += 1

    cj += 1
    ci = 0

contour_plot()