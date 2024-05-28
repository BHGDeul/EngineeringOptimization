from setup import Kite

def calculate_opt_gamma_nominal(self):
    ## Define gamma range ##
    # Prohibits reel-in speed from exceeding max reeling speed #
    # TODO: can we add this as constraint?
    if self.max_reel_speed <= 2 * self.v_w_n:
        lim = self.max_reel_speed / self.v_w_n
    else:
        lim = 2

    # gamma_in = np.linspace(0.01, lim, 100)
    # gamma_out = np.linspace(0.01, 1, 100)
    #
    # ## Set empty arrays ##
    # power_array_m = np.zeros((100, 100))
    # power_array_e = np.zeros((100, 100))
    # ## Initiate counters ##
    # ci = 0
    # cj = 0
    # for j in gamma_out:
    #     for i in gamma_in:
    #         power_array_m[cj][ci] = self.wind_power() * self.A_proj * (
    #                 self.F_out() * (1 - j) ** 2 - (self.F_in() * (1 + i) ** 2)) * ((j * i) / (j + i))
    #
    #         power_array_e[cj][ci] = self.wind_power() * self.A_proj * (
    #                 self.eff_out * self.F_out() * (1 - j) ** 2 - (self.F_in() * (1 + i) ** 2) /
    #                 self.eff_in) * ((j * i) / (j + i))
    #         ci += 1
    #
    #     cj += 1
    #     ci = 0

    gamma_in = self.wind_power()

    # TODO: here we have a maximization
    self.P_elec_opt_gamma = np.amax(power_array_e)
    (a, b) = np.where(power_array_e == self.max_power_e)

    self.gamma_out_n = gamma_out[a][0]
    self.gamma_in_n = gamma_in[b][0]

    return self