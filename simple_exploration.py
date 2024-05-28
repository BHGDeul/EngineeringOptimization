from constants.input import get_initial_data
import numpy as np
import matplotlib.pyplot as plt

# data = get_initial_data()
data = {}
data['v_w_adj'] = np.linspace(0.5, 25, 100)
data['A_proj_list'] = np.linspace(7, 25, 25)
data['F_out_list'] = np.linspace(10, 150, 20)  # np.linspace(10, 100, 90)
data['v_w_n'] = 10  # np.linspace(5,15,20) #10 #10.44#np.linspace(5,20,50)
data['CL_out'] = 1.06  # np.linspace(0.6,1.5,50)
data['A_proj'] = 16.65  # 12.302 #10.04#15.681 #9.723#9.34#21.54#24.19 #8.18 #24.19 #20.79 #21.57 #np.linspace(15,35,50) #19.8 #
data['T_out_target'] = 10300
data['T_out_max'] = 10300
data['rho'] = 1.18
data['lc'] = 250
data['CD_out'] = 0.15
# data['CL_in'] = 0.1
data['CD_in'] = 0.10
data['eff_in'] = 0.639  # 0.652 #were switched
data['eff_out'] = 0.652  # 0.639
# data['gamma_out_n'] = 0.43
# data['gamma_in_n'] = 1.78
data['P_avg_e_req'] = 20000  # Nominal electrical power (W)
data['max_reel_speed'] = 25  # m/s

data['a_elev_out'] = 30 * np.pi / 180
data['a_elev_in'] = 70 * np.pi / 180

# data['SF_supercap'] = 1.2  # Safety factor supercap
# data['diameter_drum'] = 0.46586
# data['drum_circum'] = data['diameter_drum'] * np.pi
# data['SF_force'] = 1.35
# data['rpm_min'] = 1500
# data['rpm_n'] = 1525
# data['rpm_max'] = 3000
# data['rpm_motor'] = 750
# data['Generator_lim'] = 50000

## Intermediate calculations ##

data['F_out'] = data['CL_out'] ** 3 / data['CD_out'] ** 2
data['F_in'] = data['CD_in']
data['P_w'] = 0.5 * data['v_w_n'] ** 3 * data['rho']  # Wind Power


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
plot_gamma_data['gamma_in'] = np.linspace(2, lim, resolution)
plot_gamma_data['gamma_out'] = np.linspace(0.01, 1, resolution)

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

    ## Find maximal mechanical power  ##

# data['max_power_m'] = np.amax(plot_gamma_data['power_array_m'])
data['P_elec_opt_gamma'] = np.amax(plot_gamma_data['power_array_e'])
(a, b) = np.where(plot_gamma_data['power_array_e'] == data['P_elec_opt_gamma'])

data['gamma_out_n'] = plot_gamma_data['gamma_out'][b][0]
data['gamma_in_n'] = plot_gamma_data['gamma_in'][a][0]

print(data['gamma_out_n'], data['gamma_in_n'], data['P_elec_opt_gamma'])



# meshgrid for plotting
# gamma_out_mesh, gamma_in_mesh = np.meshgrid(plot_gamma_data['gamma_out'], plot_gamma_data['gamma_in'])

# electrical power array for gammas
plt.figure(figsize=(12, 10))
contour = plt.contourf(plot_gamma_data['gamma_out'], plot_gamma_data['gamma_in'], plot_gamma_data['power_array_e'], levels=200, cmap='viridis')
plt.colorbar(contour, label='Electrical Power (W)')
plt.xlabel('Gamma Out')
plt.ylabel('Gamma In')
plt.title('Electrical Power Distribution')

# true maximum point
plt.scatter(data['gamma_out_n'], data['gamma_in_n'], color='red', label='True Max Power', zorder=5)

# Annotate the points
plt.annotate(f'True Max\n({data["gamma_out_n"]:.2f}, {data["gamma_in_n"]:.2f})', xy=(data['gamma_out_n'], data['gamma_in_n']),
             xytext=(data['gamma_out_n'] + 0.1, data['gamma_in_n'] + 0.1),
             arrowprops=dict(facecolor='red', shrink=0.05))

plt.legend()

# plt.show()
