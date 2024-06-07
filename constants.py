data = {}
data['CL_out'] = 1.06

# TODO: make data['A_proj'] a design variable
data['T_out_target'] = 10300
data['T_out_max'] = 10300
data['rho'] = 1.18
data['lc'] = 250
data['CD_out'] = 0.15
# data['CL_in'] = 0.1
data['CD_in'] = 0.10

# TODO: optimally, we could change these efficiencies to stochastic variables!
data['eff_in'] = 0.639  # 0.652 #were switched
data['eff_out'] = 0.652  # 0.639

data['P_avg_e_req'] = 20000  # Nominal electrical power (W)
data['max_reel_speed'] = 25  # m/s

## Intermediate calculations ##
data['F_out'] = data['CL_out'] ** 3 / data['CD_out'] ** 2
data['F_in'] = data['CD_in']

# baseline design values
data['baseline_power'] = 20000
data['baseline_area'] = 16.65

