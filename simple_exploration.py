import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, NonlinearConstraint

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

# TODO: optinally, we could change these efficiencies to stochastic variables!
data['eff_in'] = 0.639  # 0.652 #were switched
data['eff_out'] = 0.652  # 0.639

data['P_avg_e_req'] = 20000  # Nominal electrical power (W)
data['max_reel_speed'] = 25  # m/s

data['a_elev_out'] = 30 * np.pi / 180
data['a_elev_in'] = 70 * np.pi / 180

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

"""Brute Force Solution"""
data['P_elec_opt_gamma'] = np.amax(plot_gamma_data['power_array_e'])
(a, b) = np.where(plot_gamma_data['power_array_e'] == data['P_elec_opt_gamma'])

bf_y_out =  plot_gamma_data['gamma_out'][a][0]
bf_y_in = plot_gamma_data['gamma_in'][b][0]
print("Brute Force gamma in, gamma out =", bf_y_in, bf_y_out, "max_power",data['P_elec_opt_gamma']
      )


"""%%%%% Optimisation using scipy.optimize.minimize %%%%%"""

"""Objective function"""
# Objective function to minimize (negative of electrical power)
def objective(gamma):
    gamma_out, gamma_in = gamma
    power = data['P_w'] * data['A_proj'] * (
        data['eff_out'] * data['F_out'] * (1 - gamma_out) ** 2 -
        (data['F_in'] * ( gamma_in** 2 + 1 + 2) / data['eff_in']) * ((gamma_out * gamma_in) / (gamma_out + gamma_in)))
    return -power  # Negative because we want to maximize the power


"""Constraints"""
def constraint_traction(gamma):
    # TODO: rewrite in negative-null form
    # Assume gamma_in and gamma_out are the nominal values
    gamma_in, gamma_out = gamma

    T_out_n = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (1 - gamma_out) ** 2 * \
                   data['F_out']
    #
    # # @Winter, what is a_elev_out?

    T_in_n = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (1 - gamma_in) ** 2 * \
                   data['F_in']

    return T_out_n, T_in_n

def constraint_reel_speed(gamma):
    gamma_in, gamma_out = gamma
    return gamma_in, gamma_out / data['max_reel_speed'] / data['v_w_n'] - 1

traction_constraint = NonlinearConstraint(constraint_traction, [0, 0], [10E3, 10E3]) # This gives an inequality constraint
reel_speed_constraint = NonlinearConstraint(constraint_reel_speed, [0, -1], [1, 0]) # This gives an inequality constraint

# Initial guess
initial_guess = [0.5, 0.5]

# Perform the optimization
# TODO: find out how SLSQP algorithm works, only scipy algorithm that supports constraints
result = minimize(objective, initial_guess, bounds=None, method='SLSQP',
                  constraints=[traction_constraint, reel_speed_constraint],
                  options={'disp': True})

# Extract the optimal values
optimal_gamma_out, optimal_gamma_in = result.x
max_electrical_power = -result.fun

# Update the data dictionary with optimal values
data['gamma_out_n'] = optimal_gamma_out

data['gamma_in_n'] = optimal_gamma_in
data['P_elec_opt_gamma'] = max_electrical_power

if __name__ == '__main__':
    print(result.message)
    print(result.success)
    print('optimal gamma_in, gamma_out', result.x)
    print('Optimal Power output', max_electrical_power)

    plotting = True

    if plotting:
        def contour_plot():
            # electrical power array for gammas
            plt.figure(figsize=(12, 10))
            contour = plt.contourf(plot_gamma_data['gamma_out'], plot_gamma_data['gamma_in'], plot_gamma_data['power_array_e'], levels=200, cmap='viridis')
            plt.colorbar(contour, label='Electrical Power (W)')
            plt.xlabel('Gamma Out')
            plt.ylabel('Gamma In')
            plt.title('Electrical Power Distribution')

            # maximum point
            plt.scatter(data['gamma_out_n'], data['gamma_in_n'], color='red', label='Max Power', zorder=5)

            # Annotate the points
            plt.annotate(f'Max\n({data["gamma_out_n"]:.2f}, {data["gamma_in_n"]:.2f})', xy=(data['gamma_out_n'], data['gamma_in_n']),
                         xytext=(data['gamma_out_n'] + 0.1, data['gamma_in_n'] + 0.1),
                         arrowprops=dict(facecolor='red', shrink=0.05))

            plt.legend()

            plt.show()

        # contour_plot()

        import numpy as np
        import matplotlib.pyplot as plt

        # Define a function to calculate electrical power
        def electrical_power(gamma_out, gamma_in, data):
            power = data['P_w'] * data['A_proj'] * (
                data['eff_out'] * data['F_out'] * (1 - gamma_out) ** 2 -
                (data['F_in'] * (1 + gamma_in) ** 2) / data['eff_in']) * ((gamma_out * gamma_in) / (gamma_out + gamma_in))
            return power

        # gamma ranges for plotting
        gamma_out_range = np.linspace(0.01, 1, 100)
        gamma_in_range = np.linspace(0.01, lim, 100)

        # electrical power vs gamma_out for a fixed gamma_in
        plt.figure(figsize=(12, 6))
        gamma_in_fixed = 0.5
        power_values_out = [electrical_power(g, gamma_in_fixed, data) for g in gamma_out_range]
        plt.plot(gamma_out_range, power_values_out, label=f'Gamma_in = {gamma_in_fixed}')
        plt.xlabel('Gamma Out')
        plt.ylabel('Electrical Power (W)')
        plt.title('Electrical Power vs Gamma Out for fixed Gamma In')
        plt.legend()
        plt.grid(True)
        plt.show()


        #electrical power vs gamma_in for a fixed gamma_out
        plt.figure(figsize=(12, 6))
        gamma_out_fixed = 0.5
        power_values_in = [electrical_power(gamma_out_fixed, g, data) for g in gamma_in_range]
        plt.plot(gamma_in_range, power_values_in, label=f'Gamma_out = {gamma_out_fixed}')
        plt.xlabel('Gamma In')
        plt.ylabel('Electrical Power (W)')
        plt.title('Electrical Power vs Gamma In for fixed Gamma Out')
        plt.legend()
        plt.grid(True)
        plt.show()

        #electrical power vs gamma_in for a fixed gamma_out
        plt.figure(figsize=(12, 6))
        # power_values_in = [constraint_traction((gamma_in_range,gamma_out_range), g, data) for g in gamma_in_range]
        plt.plot(gamma_out_range, constraint_traction((gamma_in_range,gamma_out_range))[0])
        plt.xlabel('Gamma Out')
        plt.ylabel('Traction Out')
        plt.title('Traction Out vs Gamma Out')
        plt.grid(True)
        plt.show()

        #electrical power vs gamma_in for a fixed gamma_out
        plt.figure(figsize=(12, 6))
        # power_values_in = [constraint_traction((gamma_in_range,gamma_out_range), g, data) for g in gamma_in_range]
        plt.plot(gamma_in_range, constraint_traction((gamma_in_range,gamma_out_range))[1])
        plt.xlabel('Gamma In')
        plt.ylabel('Traction In')
        plt.title('Traction In vs Gamma In')
        plt.grid(True)
        plt.show()

        #reel speeds
        plt.figure(figsize=(12, 6))
        # power_values_in = [constraint_traction((gamma_in_range,gamma_out_range), g, data) for g in gamma_in_range]
        plt.plot(gamma_in_range, constraint_reel_speed((gamma_in_range,gamma_out_range))[1])
        plt.xlabel('Gamma In')
        plt.ylabel('Reel speed constraint')
        plt.grid(True)
        plt.show()