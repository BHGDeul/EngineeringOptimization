# Optimisation of the Main System Parameters
import numpy as np

from Drum_and_tether_design import structures_calculation
from Luchsinger.luchsingermodel.Nominal_power_cycle import run_nominal_analysis, sensitivity_analysis
from Luchsinger.luchsingermodel.InputV2 import get_initial_data
from Aero.aero_main_function import main_aero_function
from Anchoring_design import anchoring_info
from LaunchingEquipment import launching_equipment_info

A_proj = 16.65
Points = 1000000
Kite_segments = 12
N_split = 5
AoA_range_out = np.arange(8, 15.5, 0.5)
AoA_range_in = np.arange(-1, 3.5, 0.5)
TAS = 32.5
area_diff = 0.001


def iteration_aero_power(area_diff, A_proj, TAS):
    Points = 1000000
    Kite_segments = 12
    N_split = 5
    AoA_range_out = np.arange(8, 15.5, 0.5)
    AoA_range_in = np.arange(-1, 3.5, 0.5)

    A_proj_last = 0
    c = 0

    while abs(A_proj_last - A_proj) > area_diff:
        A_proj_t = A_proj_last
        A_proj_last = A_proj
        if c > 0:
            A_proj = (A_proj_t + A_proj) / 2

        c += 1
        print('Hey Bradda, the optimisation has started. Leggo! This is cycle: ', c)
        ## AERO ##
        print('Aero started')
        CL_average_out, CD_average_out, CL3_CD2_average_out, CD_average_in, A_proj, Strut_area_av, flat_area, flat_area_span, chords = main_aero_function(
            A_proj, Points, Kite_segments, N_split, AoA_range_out, AoA_range_in, TAS, Print=False)
        # CL_average_in, CD_average_in, CL3_CD2_average_in, A_proj, Strut_area_av, flat_area, flat_area_span, chords = main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_in, TAS, Print=False)
        print('Aero part finished')
        ## POWER ##
        data = get_initial_data()
        data['F_out'] = CL3_CD2_average_out
        data['F_in'] = CD_average_in
        data['CL_out'] = CL_average_out
        data['CD_out'] = CD_average_out
        # data['CL_in'] = CL_average_in
        data['CD_in'] = CD_average_in
        # data['']
        data = run_nominal_analysis(data)
        data['Strut_area_av'] = Strut_area_av
        data['flat_area'] = flat_area
        data['flat_area_span'] = flat_area_span
        data['chords'] = chords
        # for next iteration
        A_proj = data['A_proj']
        TAS = data['v_a_out']

        print('Power part finished, the old area was: ', A_proj_last, 'The new one is:', A_proj)
    print('Optimisation completed!')
    # Write to file #
    file = open("data_optim.txt", "w")
    for key, value in data.items():
        file.write('%s:%s,\n' % (key, value))
    file.close()
    print('The extended results of the analysis can be found in the data file added to the directory.')

    return data


data = iteration_aero_power(area_diff, A_proj, TAS)


# Open Data
def import_data(file_name):
    file = open(file_name)
    data = {}
    entries = ''
    for line in file:
        entries = entries + line.replace("\n", "")
    entries = entries[:-1]
    entries = entries.split(',')

    for entry in entries:
        key, value = entry.split(':')
        if '[' in value:
            print(value)
            value = value.replace('[', '')
            value = value.replace(']', '')
            value = value.split()
            data[key] = [float(x) for x in value]
        else:
            data[key] = float(value)
    return data


data = import_data("data_optim.txt")
# T_F_out = data['T_out_elev_n']
# gamma_out = data['gamma_out_n']
# gamma_in = data['gamma_in_n']


# DEFINITIONS
len_drum = 1.10
extra_len = 75
nom_load = data['T_out_elev_n']
saf_fac = 2
kite_mass_margin = 1.05
t_out = data['t_out']
t_in = data['t_in']
Strut_area_av = data['Strut_area_av']

Kite_mass_ALUULA, Tether_diameter, Tether_mass, Load, D_drum = structures_calculation(A_proj, Strut_area_av, len_drum,
                                                                                      data['a_elev_out'], extra_len,
                                                                                      nom_load, saf_fac,
                                                                                      kite_mass_margin, data['t_out'],
                                                                                      data['t_in'])

# Inputs
Riv_w = 3692
Safety_f = 1.2
Kite_f = nom_load * Safety_f

ref_kin_fric_coeff = 0.35
force_duration = 0.1

min_kite_angle = 25
max_kite_angle = 35
road_angle = 2

req_fric_coeff, sliding_d_mm = anchoring_info(Kite_f, Riv_w, road_angle, min_kite_angle, max_kite_angle,
                                              ref_kin_fric_coeff, force_duration)
vw_ground = 10.  # m/s
chord = max(data['chords'])  # m
F_clamp = launching_equipment_info(vw_ground, req_fric_coeff, data['flat_area'], chord, Kite_mass_ALUULA)