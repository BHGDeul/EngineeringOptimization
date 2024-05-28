"""
Define here in (pseudo) code the structure of the model
"""

import numpy as np

class Kite():
    def __init__(self):
        # all constant design parameters of the kite
        self.v_w_n = 10  # np.linspace(5,20,50)

        # TODO: make A_proj a design variable
        self.A_proj = 20.79  # 21.57 #np.linspace(15,35,50) #19.8 #

        self.rho = 1.18
        self.lc = 250
        self.CL_out = 1.0
        self.CD_out = 0.2
        self.CL_in = 0.14
        self.CD_in = 0.07

        # can be made variables, keep as is for simplification
        self.eff_in = 0.579
        self.eff_out = 0.591

        self.P_avg_e_n = 20000  # Nominal electrical power (W)
        self.max_reel_speed = 25  # m/s

    def F_out(self):
        return self.CL_out ** 3 / self.CD_out ** 2

    def F_in(self):
        return self.CD_in

    def wind_power(self):
        P_w = 0.5 * self.v_w_n ** 3 * self.rho  # Wind Power
        return P_w

    def power_m(self, gamma_in, gamma_out):
        # TODO: @Winter, is this the power output of the full cycle based on gamma_in and gamma_out?
        # nice multivariable function
        # function taken from find_gamma.py
        # efficiencies are excluded
        power_m = (self.wind_power() * self.A_proj *
                   (self.eff_out * self.F_out()* (1-gamma_out)**2 - (self.F_in() * (1 + gamma_in)**2) / self.eff_in) *
                   ((gamma_out * gamma_in) / (gamma_out + gamma_in))
                   )
        return power_m

    def calculate_nominal_tractionF(self, gamma_in, gamma_out):
        # Assume gamma_in and gamma_out are the nominal values

        # @Winter, is T_out_n the nominal value or the value in a normal direction?
        self.T_out_n = 0.5 * self.rho * self.v_w_n ** 2 * self.A_proj * (1 - gamma_out) ** 2 * \
                          self.F_out()

        # @Winter, what is a_elev_out?
        self.T_out_n_angle = 0.5 * self.rho * self.v_w_n ** 2 * self.A_proj * (
                np.cos(self.a_elev_out) - gamma_out) ** 2 * self.F_out()

        self.T_in_n = 0.5 * self.rho * self.v_w_n ** 2 * self.A_proj * (1 + gamma_in) ** 2 * \
                         self.F_in()

        return self

    def calculate_nominal_powers(self, gamma_in, gamma_out):
        self.P_out = self.T_out_n * gamma_out * self.v_w_n * self.eff_out

        self.P_in = self.T_in_n * gamma_in * self.v_w_n / self.eff_in

        return self


if __name__ == '__main__':
    obj = Kite()