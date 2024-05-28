import numpy as np

import constants.input


### This function finds the optimal gamma for nominal flight conditions defined in Input.py and plots this on a heat map ###
def calculate_opt_gamma_nominal(data):

    plot_gamma_data = {}
    ## Define gamma range ##
    # Prohibits reel-in speed from exceeding max reeling speed #
    if data['max_reel_speed'] <= 2* data['v_w_n']:
        lim = data['max_reel_speed'] / data['v_w_n']
    else:
        lim = 2

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
    (a, b) = np.where(plot_gamma_data['power_array_e'] == data['max_power_e'])

    data['gamma_out_n'] = plot_gamma_data['gamma_out'][a][0]
    data['gamma_in_n'] = plot_gamma_data['gamma_in'][b][0]

    return data


### This function calculates the traction forces for nominal flight conditions ###

def calculate_nominal_tractionF(data):
    data['T_out_n'] = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (1 - data['gamma_out_n']) ** 2 * data[
        'F_out']
    data['T_out_n_angle'] = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (
                np.cos(data['a_elev_out']) - data['gamma_out_n']) ** 2 * data['F_out']
    data['T_in_n'] = 0.5 * data['rho'] * data['v_w_n'] ** 2 * data['A_proj'] * (1 + data['gamma_in_n']) ** 2 * data[
        'F_in']

    return data


def calculate_nominal_powers(data):
    data['P_out'] = data['T_out_n'] * data['gamma_out_n'] * data['v_w_n']
    data['P_out_e'] = data['P_out'] * data['eff_out']

    data['P_in'] = data['T_in_n'] * data['gamma_in_n'] * data['v_w_n']
    data['P_in_e'] = data['P_in'] / data['eff_in']

    ## Sanity check, code verification ##
    # P_avg_mech = P_out*(gamma_in)/(gamma_in + gamma_out) - P_in*data['gamma_out_n']/(data['gamma_in_n'] + data['gamma_out_n'])
    # P_avg_elec = P_out_e*(gamma_in)/(gamma_in + gamma_out) - P_in_e*data['gamma_out_n']/(data['gamma_in_n'] + data['gamma_out_n'])
    return data


def calculate_updated_projected_area(data):
    data['A_proj_u'] = data['P_avg_e_n'] / data['P_w'] / ((data['eff_out'] * data['F_out'] * (
                1 - data['gamma_out_n']) ** 2 - (data['F_in'] * (1 + data['gamma_in_n']) ** 2) / data['eff_in']) * (
                                                                      (data['gamma_out_n'] * data['gamma_in_n']) / (
                                                                          data['gamma_out_n'] + data['gamma_in_n'])))
    print(data['A_proj_u'])
    return data


def calculate_cycle_param(data):
    data['cycle_time'] = (data['lc'] / data['v_w_n']) * (
                (data['gamma_out_n'] + data['gamma_in_n']) / (data['gamma_out_n'] * data['gamma_in_n']))
    data['t_in'] = data['cycle_time'] * data['gamma_out_n'] / (data['gamma_in_n'] + data['gamma_out_n'])
    data['t_out'] = data['cycle_time'] * data['gamma_in_n'] / (data['gamma_in_n'] + data['gamma_out_n'])

    return data


# def evaluate_tether_force():


def run_nominal_analysis(data):
    data = calculate_opt_gamma_nominal(data)[0]
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)
    data = calculate_updated_projected_area(data)
    data = calculate_cycle_param(data)
    plot_gamma_power()
    # Write to file #
    file = open("data.txt", "w")
    for key, value in data.items():
        file.write('%s:%s\n' % (key, value))
    file.close()
    print('These are the results for the nominal case defined in the input file:')
    print('')
    return data

if __name__ == '__main__':
    data = constants.input.get_initial_data()
    calculate_opt_gamma_nominal(data)