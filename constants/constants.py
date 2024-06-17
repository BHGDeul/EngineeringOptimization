# Data initialization
v_w_adj = np.linspace(0.5, 25, 100)
A_proj_list = np.linspace(7, 25, 25)
F_out_list = np.linspace(10, 150, 20)
v_w_n = 10
CL_out = 1.06
A_proj = 16.65
T_out_target = 10300
T_out_max = 10300
rho = 1.18
lc = 250
CD_out = 0.15
CD_in = 0.10
eff_in = 0.639
eff_out = 0.652
P_avg_e_req = 20000
max_reel_speed = 25
a_elev_out = 30 * np.pi / 180
a_elev_in = 70 * np.pi / 180

baseline_power = 20E3
baseline_area = 16.65

# Intermediate calculations
F_out = CL_out ** 3 / CD_out ** 2
F_in = CD_in
P_w = 0.5 * v_w_n ** 3 * rho  # Wind Power