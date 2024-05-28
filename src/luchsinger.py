from InputV2 import *
import numpy as np
import matplotlib.pyplot as plt


# Luchsinger.luchsingermodel.
### This function finds the optimal gamma for nominal flight conditions defined in Input.py and plots this on a heat map ###
def calculate_opt_gamma_nominal(data):
    plot_gamma_data = {}
    # Define gamma range ##
    # Prohibits reel-in speed from exceeding max reeling speed #
    if data['max_reel_speed'] <= 2 * data['v_w_n']:
        lim = data['max_reel_speed'] / data['v_w_n']
    else:
        lim = 2
    lim = 2.0
    plot_gamma_data['gamma_in'] = np.linspace(0.01, lim, 100)
    plot_gamma_data['gamma_out'] = np.linspace(0.01, 1, 100)

    ## Set empty arrays ##
    plot_gamma_data['power_array_m'] = np.zeros((100, 100))
    plot_gamma_data['power_array_e'] = np.zeros((100, 100))
    ## Initiate counters ##
    ci = 0
    cj = 0
    for j in plot_gamma_data['gamma_out']:
        for i in plot_gamma_data['gamma_in']:
            plot_gamma_data['power_array_m'][cj][ci] = data['P_w'] * data['A_proj'] * (
                        data['F_out'] * (1 - j) ** 2 - (data['F_in'] * (1 + i) ** 2)) * ((j * i) / (j + i))
            plot_gamma_data['power_array_e'][cj][ci] = data['P_w'] * data['A_proj'] * (
                        data['eff_out'] * data['F_out'] * (1 - j) ** 2 - (data['F_in'] * (1 + i) ** 2) / data[
                    'eff_in']) * ((j * i) / (j + i))
            ci += 1

        cj += 1
        ci = 0

        ## Find maximal mechanical power  ##

    # data['max_power_m'] = np.amax(plot_gamma_data['power_array_m'])
    data['P_elec_opt_gamma'] = np.amax(plot_gamma_data['power_array_e'])
    (a, b) = np.where(plot_gamma_data['power_array_e'] == data['P_elec_opt_gamma'])

    data['gamma_out_n'] = plot_gamma_data['gamma_out'][a][0]
    data['gamma_in_n'] = plot_gamma_data['gamma_in'][b][0]

    return data


# def calculate_opt_gamma_nominal_elev(data):
#     plot_gamma_data = {}
#     # Define gamma range ##
#     # Prohibits reel-in speed from exceeding max reeling speed #
#     if data['max_reel_speed'] <= 2 * data['v_w_n']:
#         lim = data['max_reel_speed'] / data['v_w_n']
#     else:
#         lim = 2.5
#
#     plot_gamma_data['gamma_in'] = np.linspace(0.01, lim, 100)
#     plot_gamma_data['gamma_out'] = np.linspace(0.01, 1, 100)
#
#
#     ## Set empty arrays ##
#     plot_gamma_data['power_array_m'] = np.zeros((100, 100))
#     plot_gamma_data['power_array_e'] = np.zeros((100, 100))
#     ## Initiate counters ##
#     ci = 0
#     cj = 0
#     for j in plot_gamma_data['gamma_out']:
#         for i in plot_gamma_data['gamma_in']:
#             # plot_gamma_data['power_array_m'][cj][ci] = data['P_w']*data['A_proj']*(data['F_out']*(np.cos(data['a_elev_out'])-j)**2-(data['F_in']*(i**2+2*np.cos(data['a_elev_in'])*i+1)))*((j*i)/(j+i))
#             plot_gamma_data['power_array_e'][cj][ci] = data['A_proj'] * (
#                         data['eff_out'] * data['F_out'] * (np.cos(data['a_elev_out']) - j) ** 2 - (
#                             data['F_in'] * (i ** 2 + 2 * np.cos(data['a_elev_in']) * i + 1)) / data['eff_in']) * (
#                                                                    (j * i) / (j + i))
#             ci += 1
#         # data['P_w']*
#         cj += 1
#         ci = 0
#
#         ## Find maximal mechanical power  ##
#
#     # data['max_power_m'] = np.amax(plot_gamma_data['power_array_m'])
#     data['P_elec_opt_gamma'] = np.amax(plot_gamma_data['power_array_e'])
#     # (a,b) = np.where(power_array_m == max_power_m)
#     (a, b) = np.where(plot_gamma_data['power_array_e'] == data['P_elec_opt_gamma'])
#
#     data['gamma_out_n'] = plot_gamma_data['gamma_out'][a][0]
#     data['gamma_in_n'] = plot_gamma_data['gamma_in'][b][0]
#
#     return data, plot_gamma_data


### This function finds the optimal gamma in when a fixed gamma out is set ###

def calculate_opt_gamma_in(data):
    plot_gamma_data = {}
    ## Define gamma range ##
    # Prohibits reel-in speed from exceeding max reeling speed #
    if data['max_reel_speed'] <= 2 * data['v_w_n']:
        lim = data['max_reel_speed'] / data['v_w_n']
    else:
        lim = 2.5

    plot_gamma_data['gamma_in'] = np.linspace(0.01, lim, 100)

    ## Set empty arrays ##
    plot_gamma_data['power_array_m'] = np.zeros(100)
    plot_gamma_data['power_array_e'] = np.zeros(100)
    ## Initiate counters ##
    ci = 0
    cj = 0

    for i in plot_gamma_data['gamma_in']:
        plot_gamma_data['power_array_e'][ci] = data['P_w'] * data['A_proj'] * (
                    data['eff_out'] * data['F_out'] * (np.cos(data['a_elev_out']) - data['gamma_out_n']) ** 2 - (
                        data['F_in'] * (i ** 2 + 2 * np.cos(data['a_elev_in']) * i + 1)) / data['eff_in']) * (
                                                           (data['gamma_out_n'] * i) / (data['gamma_out_n'] + i))
        ci += 1

    ## Find maximal mechanical power  ##

    # data['max_power_m'] = np.amax(plot_gamma_data['power_array_m'])
    data['P_elec_opt_gamma_in'] = np.amax(plot_gamma_data['power_array_e'])
    # (a,b) = np.where(power_array_m == max_power_m)
    b = np.where(plot_gamma_data['power_array_e'] == data['P_elec_opt_gamma_in'])

    # data['gamma_out_n'] = plot_gamma_data['gamma_out'][a][0]
    data['gamma_in_n'] = plot_gamma_data['gamma_in'][b][0]

    return data

### This function calculates the traction forces for nominal flight conditions ###

def calculate_nominal_tractionF(data):
    data['T_out_n'] = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (1 - data['gamma_out_n']) ** 2 * data[
        'F_out']
    data['T_out_elev_n'] = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (
                np.cos(data['a_elev_out']) - data['gamma_out_n']) ** 2 * data['F_out']

    data['T_in_n'] = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (1 + data['gamma_in_n']) ** 2 * data[
        'F_in']
    data['T_in_elev_n'] = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (
                1 + 2 * data['gamma_in_n'] * np.cos(data['a_elev_in']) + data['gamma_in_n'] ** 2) * data['F_in']

    return data


### This function calculates the powers during all phases of the power cycle ###

def calculate_nominal_powers(data):
    data['P_out'] = data['T_out_n'] * data['gamma_out_n'] * data['v_w_n']
    data['P_out_e'] = data['P_out'] * data['eff_out']

    data['P_out_elev'] = data['T_out_elev_n'] * data['gamma_out_n'] * data['v_w_n']
    data['P_out_e_elev'] = data['P_out_elev'] * data['eff_out']

    data['P_in'] = data['T_in_n'] * data['gamma_in_n'] * data['v_w_n']
    data['P_in_e'] = data['P_in'] / data['eff_in']

    data['P_in_elev'] = data['T_in_elev_n'] * data['gamma_in_n'] * data['v_w_n']
    data['P_in_e_elev'] = data['P_in_elev'] / data['eff_in']

    data['P_avg_mech'] = data['P_w'] * data['A_proj'] * (
                data['F_out'] * (1 - data['gamma_out_n']) ** 2 - (data['F_in'] * (data['gamma_in_n'] + 1) ** 2)) * (
                                     (data['gamma_out_n'] * data['gamma_in_n']) / (
                                         data['gamma_out_n'] + data['gamma_in_n']))
    data['P_avg_elec'] = data['P_w'] * data['A_proj'] * (
                data['eff_out'] * data['F_out'] * (1 - data['gamma_out_n']) ** 2 - (
                    data['F_in'] * (data['gamma_in_n'] + 1) ** 2) / data['eff_in']) * (
                                     (data['gamma_out_n'] * data['gamma_in_n']) / (
                                         data['gamma_out_n'] + data['gamma_in_n']))

    data['P_avg_mech_elev'] = data['P_w'] * data['A_proj'] * (
                data['F_out'] * (np.cos(data['a_elev_out']) - data['gamma_out_n']) ** 2 - (data['F_in'] * (
                    data['gamma_in_n'] ** 2 + 2 * np.cos(data['a_elev_in']) * data['gamma_in_n'] + 1))) * (
                                          (data['gamma_out_n'] * data['gamma_in_n']) / (
                                              data['gamma_out_n'] + data['gamma_in_n']))
    data['P_avg_elec_elev'] = data['P_w'] * data['A_proj'] * (
                data['eff_out'] * data['F_out'] * (np.cos(data['a_elev_out']) - data['gamma_out_n']) ** 2 - (
                    data['F_in'] * (data['gamma_in_n'] ** 2 + 2 * np.cos(data['a_elev_in']) * data['gamma_in_n'] + 1)) /
                data['eff_in']) * ((data['gamma_out_n'] * data['gamma_in_n']) / (
                data['gamma_out_n'] + data['gamma_in_n']))

    data['max_cycle_power_elev'] = data['P_w'] * data['A_proj'] * data['F_out'] * (
        np.cos(data['a_elev_out'])) ** 3 * 4 / 27

    ## Sanity check, code verification ##
    data['P_avg_mech_verif'] = data['P_out'] * (data['gamma_in_n']) / (data['gamma_in_n'] + data['gamma_out_n']) - data[
        'P_in'] * data['gamma_out_n'] / (data['gamma_in_n'] + data['gamma_out_n'])
    data['P_avg_elec_verif'] = data['P_out_e'] * (data['gamma_in_n']) / (data['gamma_in_n'] + data['gamma_out_n']) - \
                               data['P_in_e'] * data['gamma_out_n'] / (data['gamma_in_n'] + data['gamma_out_n'])

    data['P_avg_mech_elev_verif'] = data['P_out_elev'] * (data['gamma_in_n']) / (
                data['gamma_in_n'] + data['gamma_out_n']) - data['P_in_elev'] * data['gamma_out_n'] / (
                                                data['gamma_in_n'] + data['gamma_out_n'])
    data['P_avg_elec_elev_verif'] = data['P_out_e_elev'] * (data['gamma_in_n']) / (
                data['gamma_in_n'] + data['gamma_out_n']) - data['P_in_e_elev'] * data['gamma_out_n'] / (
                                                data['gamma_in_n'] + data['gamma_out_n'])

    return data


### This function updates the projected area based on the average power requirement ###

# def calculate_updated_projected_area(data):
#     data['A_proj_u'] = data['P_avg_e_req'] / data['P_w'] / ((data['eff_out'] * data['F_out'] * (
#                 np.cos(data['a_elev_out']) - data['gamma_out_n']) ** 2 - (data['F_in'] * (
#                 data['gamma_in_n'] ** 2 + 2 * np.cos(data['a_elev_in']) * data['gamma_in_n'] + 1)) / data['eff_in']) * (
#                                                                         (data['gamma_out_n'] * data['gamma_in_n']) / (
#                                                                             data['gamma_out_n'] + data['gamma_in_n'])))
#     # data['A_proj']= data['A_proj_u']
#     return data


### This function calculates the cycle times ###

def calculate_cycle_param(data):
    data['cycle_time'] = (data['lc'] / data['v_w_n']) * (
                (data['gamma_out_n'] + data['gamma_in_n']) / (data['gamma_out_n'] * data['gamma_in_n']))
    data['t_in'] = data['cycle_time'] * data['gamma_out_n'] / (data['gamma_in_n'] + data['gamma_out_n'])
    data['t_out'] = data['cycle_time'] * data['gamma_in_n'] / (data['gamma_in_n'] + data['gamma_out_n'])

    return data


### This function looks for gamma outs to lower the traction force, and directly updates the area accordingly
def evaluate_tether_force(data):
    data['gamma_out_n_init'] = data['gamma_out_n']
    data['gamma_in_n_init'] = data['gamma_in_n']
    data['A_proj_init'] = data['A_proj']

    TF_an = {'gamma_out': [], 'force': [], 'area': []}
    while data['gamma_out_n'] < .6:
        data['gamma_out_n'] += 0.01
        data = calculate_opt_gamma_in(data)
        data = calculate_updated_projected_area(data)
        data['A_proj'] = data['A_proj_u']

        data = calculate_nominal_tractionF(data)

        TF_an['gamma_out'].append(data['gamma_out_n'])
        TF_an['force'].append(data['T_out_elev_n'])
        TF_an['area'].append(data['A_proj'])

    data['A_proj'] = data['A_proj_init']
    data['gamma_out_n'] = np.cos(data['a_elev_out']) - np.sqrt(
        data['T_out_target'] * 2 / (data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * data['F_out']))
    data = calculate_opt_gamma_in(data)
    data = calculate_updated_projected_area(data)
    if data['A_proj'] < data['A_proj_u']:
        print('ok')

    else:
        print('not ok')
    print(data['A_proj'])
    print(data['A_proj_u'])
    data['A_proj'] = data['A_proj_u']
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)

    return data

### Do Analysis ###
def run_nominal_analysis(data):
    ip = 1
    # ip = int(input('Enter 1 if you want to take the tether elevation into account for finding the optimal reeling speeds, 0 for ignoring it: '))
    if ip == 0:
        data = calculate_opt_gamma_nominal(data)[0]
        data_plot = calculate_opt_gamma_nominal(data)[1]
    elif ip == 1:
        data = calculate_opt_gamma_nominal_elev(data)[0]
        data_plot = calculate_opt_gamma_nominal_elev(data)[1]
    else:
        print('Enter a valid number!')
        quit()
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)
    data = calculate_updated_projected_area(data)
    if abs(data['A_proj_u'] - data['A_proj']) > 0.01:
        print('The projected area of the kite should be iterated on!')
    else:
        print('The area of the kite is optimal for the required power output.')
    # data['A_proj'] = data['A_proj_u']
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)
    data = calculate_cycle_param(data)
    data = calculate_apparent_speed(data)
    data = size_supercap(data)

    plot_gamma_power(data_plot)

    # Write to file #
    file = open("Luchsinger\data.txt", "w")
    for key, value in data.items():
        file.write('%s:%s\n' % (key, value))
    file.close()
    print('The extended results of the analysis can be found in the data file added to the directory.')
    # ip2 = input('The tetherforce is now', data['T_out_elev_n'], 'enter 1 if you want to analyse how to lower it, else enter 0: ')
    ip2 = 1
    if ip2 == 0:
        quit()
    elif ip2 == 1:
        run_TF_anal(data)

    else:
        print('Enter a valid number!')
        quit()

    return data


def run_TF_anal(data):
    data = evaluate_tether_force(data)
    data = calculate_cycle_param(data)
    data = calculate_apparent_speed(data)
    data = size_supercap(data)
    data = size_generator(data)

    # Write to file #
    file = open("Luchsinger\data_TF_an.txt", "w")
    for key, value in data.items():
        file.write('%s:%s\n' % (key, value))
    file.close()
    print('The extended results of the analysis can be found in the data file added to the directory.')

data = run_nominal_analysis(get_initial_data())







