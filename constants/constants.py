import numpy as np

# Data initialization
v_w_adj = np.linspace(0.5, 25, 100)
A_proj_list = np.linspace(7, 25, 25)
F_out_list = np.linspace(10, 150, 20)

CL_out = 1.06
max_traction = 10E3
rho = 1.18
lc = 250
CD_out = 0.15
CD_in = 0.10
eff_in = 0.639
eff_out = 0.652
P_avg_e_req = 20000
max_reel_speed = 2.5
# a_elev_out = 30 * np.pi / 180
# a_elev_in = 70 * np.pi / 180

baseline_power = 40E3
baseline_area = 30#16.65

# Intermediate calculations
F_out = CL_out ** 3 / CD_out ** 2
F_in = CD_in