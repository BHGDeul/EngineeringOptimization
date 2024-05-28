# -*- coding: utf-8 -*-
"""
Created on Mon May 23 11:46:02 2022

@author: Joris
"""
from math import sin,radians,pi,ceil
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from tabulate import tabulate
from scipy.optimize import curve_fit


#DEFINITIONS
kite_area_in = 20.16
avg_strut_in = 0.6867
len_drum_in = 1.10
angle_in = 30
extra_len_in = 75
nom_load_in = 10329.316
saf_fac_in = 2
kite_mass_margin = 1.05
t_out = 103.5
t_in = 11.5


def structures_calculation(kite_area_in, avg_strut_in, len_drum_in, angle_in, extra_len_in, nom_load_in, saf_fac_in, kite_mass_margin, t_out, t_in):
    #(double _winding = boolean True if double winded, angle = kite angle, oper_range = operative range in list, tether_d = tether diameter, saf_marg = safety margin, dD_dT = drum diameter of tether diameter, extra_len = extra length of tether )
    def drum_sizing(double_winding , angle, oper_range, tether_d, saf_marg, dD_dT, extra_len):
        # general calculations for drum sizing
        height = sin(radians(angle))
        max_tet_len = np.array(oper_range)/height

        #specific parts for if double winded or not
        if double_winding:
            drum_diam_oper = (tether_d)*dD_dT
            drum_circ_oper = pi*drum_diam_oper

            drum_diam_start = drum_diam_oper+ 2*tether_d
            drum_circ_start = pi*drum_diam_start
            oper_len = max_tet_len[1] - max_tet_len[0] +extra_len
            print(oper_len, max_tet_len[0])

            winds_oper = ceil(oper_len/drum_circ_oper)
            winds_start = ceil(max_tet_len[0]/drum_circ_start)
            len_drum_oper = (saf_marg + tether_d)*winds_oper
            len_drum_start = (saf_marg + tether_d)*winds_start

            len_drum = 1.10 #Fixed width instead of max(len_drum_oper,len_drum_start)
            if len_drum_oper< len_drum_start:
                print("Operation drum length is too short",len_drum_oper-len_drum_start, "m.")
            elif len_drum_oper == len_drum_start:
                print("Operation drum length is equal to starting phase")
            else:
                print(f'{"Operation drum length is longer than starting phase : "}{len_drum_oper-len_drum_start}{" m longer"}')
        else:
            winds = ceil((max_tet_len[1]+extra_len)/(pi*drum_diam_oper))
            len_drum = (saf_marg + tether_d)*winds

        return len_drum, drum_diam_oper, max_tet_len[0], max_tet_len[1]

    def drum_sizing_fixed_len(double_winding, angle, oper_range, tether_d, saf_marg, len_drum, extra_len):
        # general calculations for drum sizing
        height = sin(angle)
        max_tet_len = np.array(oper_range) / height
        # specific parts for if double winded or not
        if double_winding:
            oper_len = max_tet_len[1] - max_tet_len[0] + extra_len
            winds = ceil(len_drum/(saf_marg + tether_d))
            if max_tet_len[0] > (oper_len):
                drum_diam_start = (max_tet_len[0]/winds)/pi
                drum_diam = drum_diam_start
            if max_tet_len[0] < (oper_len):
                drum_diam_oper = (oper_len/winds)/pi + 2 * tether_d
                drum_diam = drum_diam_oper
        dD_dT = drum_diam/tether_d

        return drum_diam, dD_dT, max_tet_len[0], max_tet_len[1]

    #str_value = breaking strength, lst_str = list of breaking strengths of tether, lst_diam = list of diameters of tether, lst_weight = list of weights of tether to interpolate
    def interpolate_t_str_diam(str_value, lst_diam, lst_str, lst_weight):
        #interpolates both str diam and diam weight to continuously estimate both diameter and weight based on discrete sizes.
        f_str_diam = CubicSpline(lst_str,lst_diam)
        f_diam_weight = CubicSpline(lst_diam,lst_weight)
        diam = f_str_diam(str_value)
        weight = f_diam_weight(diam)
        return f_str_diam, f_diam_weight, diam, weight

    def comp_tether(ar,str_value,winner_ind):
        lst = [["hps"],["bio"],["sierra"]]
        for i in range(len(ar)):
            f_s_d, f_d_w, d , w = interpolate_t_str_diam(str_value, ar[i][0],ar[i][1],ar[i][2])
            lst[i].append(d)
            lst[i].append(w)
        header = ["Tether", "Diameter (mm)", "Weight per 100 m (kg/100m)"]
        print(tabulate(lst,headers = header))
        return lst[winner_ind],lst

    def total_SSL_factors(reel_out_perc):
        f_L = 100/reel_out_perc      #Load factor
        f_T = 1.5                   #Day/night temperature fluctuation
        f_S = 1.5                   #Seasonal temperature fluctuation
        return f_L*f_T*f_S          #Total 'safety' factor

    def tension(nom_load, diam):
        tens = nom_load/(((diam/2)**2)*pi)
        return tens

    def get_equation(x,y):
        degree = 2
        coefs, res, _, _, _ = np.polyfit(x,y,degree, full = True)
        ffit = np.poly1d(coefs)
        plt.plot(x, ffit(x))
        return ffit

    def canopy_weight(area, mat_weight):
        lstweights = []
        for i in range(len(area)):
            weight = area[i]*mat_weight
            lstweights.append(weight)
        return lstweights

    def strut_LE_perc(lstweights):
        totalweight = []
        for i in range(len(lstm2)):
            weight = (lstAL[i]-lstweights[i])/lstAL[i]
            totalweight.append(weight)
        return totalweight

    def extrapolation_perc(value):
        popt, pcov = curve_fit(lambda fx,a,b: a*fx**-b,  lstm2,strut_LE_perc(canopy_weight(lstm2, DA_weight)))
        sizes = np.linspace(min(lstm2), 35, 1000)
        percentage = popt[0]*value**-popt[1]
        power_plot = popt[0]*sizes**-popt[1]
        #plt.scatter(lstm2,strut_LE_perc(canopy_weight(lstm2, DA_weight)), label='actual data')
        #plt.plot(sizes, power_plot, label='smooth-power-fit', color = 'b')
        #plt.xlabel("Kite area m$^2$", fontsize = 16)
        #plt.ylabel("Air frame fractional weight -", fontsize = 16)
        return percentage

    def total_weight(area, mat_weight, percentage):
        extra_strut = avg_strut_A*AL_weight*2.5
        kite_weight = canopy_weight(area, mat_weight)/(1-percentage)
        return kite_weight[0]+extra_strut

    #INPUTS

    ar_tether_char = [[[2,3,4,5,6,8,10], [500,1000,1700,3000,3800,7500,10000],[0.3,0.6,0.9,1.6,2.1,4.2,5.5]]
    ,[[2,3,4,5,6,8,10],[450,900,1400,2400,3190,5590,7980],[0.24,0.48,0.84,1.43,1.9,3.3,4.8]]
    ,[[2,3,4,5,6,8,10],[380,1045,1620,2860,3665,7055,9530],[0.2,0.5,0.8,1.5,2.5,3.6,5.3]]]

    saf_fac = saf_fac_in
    nom_load= nom_load_in #10192 #10084.273 #11500 #12000 #15308.556 #14000.11188 #16219.636
    tether_char,tether_list = comp_tether(ar_tether_char,nom_load/9.81*saf_fac,2)

    # get_equation(ar_tether_char[0][0], ar_tether_char[0][1])

    # Drum dimensions
    angle = angle_in #30                          #deg
    oper_range = [165,290]
    tether_d = tether_char[1]/1000
    saf_marg = 0.0005                   #m
    dD_dT = 105                         #drum/tether
    extra_len = extra_len_in                    #m

    # Extrapolate kite weight
    lstm2 = [4,5,6,7,8,9,10,12,14]
    lstm2 = [7,8,9,10,11,12]
    #lstAL = [1.3,1.4,1.5,1.7,1.8,1.9,2,2.3,2.5]
    lbs = 0.4535
    lstAL = [4.65*lbs,5*lbs,5.35*lbs,5.65*lbs,6.1*lbs,6.5*lbs]
    AL_weight = 0.084                   # kg/m^2
    DA_weight = 0.054                   # kg/m^2
    kite_area = [kite_area_in]#[20.16]                 # m^2
    reel_out_perc = 100 * (t_out) / (t_out + t_in) #89                  #Update when possible
    avg_strut_A = avg_strut_in #0.6867                # m^2

    diam = tether_list[2][1]*10**(-3) #is this also hardcoded?

    # len_drum, d_drum, tet_len_reel_out_start, tet_len_reel_out_end = drum_sizing(True, angle, oper_range, tether_d, saf_marg, dD_dT, extra_len)

    d_drum, dD_dT, tet_len_reel_out_start, tet_len_reel_out_end  = drum_sizing_fixed_len(True, angle, oper_range, tether_d, saf_marg, len_drum_in, extra_len)
    print('dD_dT =', dD_dT, '-' )
    #Bridle diameters
    #brid_forces = [3710.24, 2497.78, 1643.51, 1522.16, 1187.60, 1072.19, 810.32]
    #brid_lst = []
    #for i in range(len(brid_forces)):
    #    bridle_char, bridle_list = comp_tether(ar_tether_char,(brid_forces[i]/9.81)*saf_fac,2)
    #    brid_lst.append(bridle_char[1]/1000)
    #print(brid_lst)

    #OUTPUTS
    # Len_drum = len_drum_in
    D_drum = d_drum
    Total_factor = total_SSL_factors(reel_out_perc)
    Load = round(tension(nom_load, diam)/10**6, 3)
    Tether_diameter= (diam)
    Tether_density = 1/((diam/2)**2*pi)/100
    Tether_mass = (tether_list[2][2])*6.55
    Tether_volume = (tet_len_reel_out_end+extra_len)*(diam/2)**2*pi
    #ALUULA weight = get_equation(lstm2, lstAL)(30)}{" kg"}')
    #Dacron weight = get_equation(lstm2, lstDA)(30)}{" kg"}')

    Kite_weight_Teijin = total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))
    Kite_weight_ALUULA = kite_mass_margin *(total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))+ canopy_weight(kite_area, AL_weight)[0]-canopy_weight(kite_area, DA_weight)[0])
    Canopy_weight = canopy_weight(kite_area, AL_weight)[0]
    Airframe_mass = total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))+canopy_weight(kite_area, AL_weight)[0]-\
                    canopy_weight(kite_area, DA_weight)[0]-canopy_weight(kite_area, AL_weight)[0]
    Valve_reinforcement_mass = (((total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))+ canopy_weight(kite_area, AL_weight)[0]-canopy_weight(kite_area, DA_weight)[0]-canopy_weight(kite_area, AL_weight)[0])/AL_weight)-25)*AL_weight

    #Print OUTPUTS
    print(f'{"Len_drum = "}{len_drum_in}{" m"}')
    print(f'{"D_drum = "}{d_drum}{" m"}')
    print(f'{"Total factor = "}{total_SSL_factors(reel_out_perc)}{" [-]"}')
    print(f'{"Load = "}{round(tension(nom_load, diam)/10**6, 3)}{" MPa"}')  #Read fig 33.16, multiply value from graph by factor
    print(f'{"Tether diameter= "}{(diam)}{"mm"}')
    print(f'{"Tether density = "}{1/((diam/2)**2*pi)/100}{" kg/m^3"}')
    print(f'{"Tether mass = "}{(tether_list[2][2])*6.55}{" kg"}')
    print(f'{"Tether volume = "}{(tet_len_reel_out_end+extra_len)*(diam/2)**2*pi}{" m^3"}')
    #print(f'{"ALUULA weight = "}{get_equation(lstm2, lstAL)(30)}{" kg"}')
    #print(f'{"Dacron weight = "}{get_equation(lstm2, lstDA)(30)}{" kg"}')

    print(f'{"Kite weight with Teijin canopy = "}{total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))}{" kg"}')
    print(f'{"Kite weight with ALUULA canopy = "}{kite_mass_margin *(total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))+ canopy_weight(kite_area, AL_weight)[0]-canopy_weight(kite_area, DA_weight)[0])}{" kg"}')
    print(f'{"Canopy weight = "}{canopy_weight(kite_area, AL_weight)[0]}{" kg"}{","}{((canopy_weight(kite_area, AL_weight)[0])/(total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))+ canopy_weight(kite_area, AL_weight)[0]-canopy_weight(kite_area, DA_weight)[0]))*100}{"%"}')
    print(f'{"Airframe mass = "}{total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))+ canopy_weight(kite_area, AL_weight)[0]-canopy_weight(kite_area, DA_weight)[0]-canopy_weight(kite_area, AL_weight)[0]}{" kg"}{","}{(1-((canopy_weight(kite_area, AL_weight)[0])/(total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))+ canopy_weight(kite_area, AL_weight)[0]-canopy_weight(kite_area, DA_weight)[0])))*100}{"%"}')
    print(f'{"Valve and reinforcement mass = "}{(((total_weight(kite_area, DA_weight, extrapolation_perc(kite_area))+ canopy_weight(kite_area, AL_weight)[0]-canopy_weight(kite_area, DA_weight)[0]-canopy_weight(kite_area, AL_weight)[0])/AL_weight)-25)*AL_weight}{" kg"}')

    """The kite weight is determined using a datasheet of a Rise A-Series kite. This kite was used since it has an aspect ratio of 5 and 7 struts, 
    thus it is easier to approximate a 8 strut system with the same aspect ratio. The flat kite area is used as input, to extrapolate the percentage of weight that the struts take. This can be seen in ... .
    The extrapolation is based on a Teijin D2 54 gsm canopy and a full ALUULA airframe. To make the system fully ALUULA, the weight of the canopy is subtracted and the weight of an ALUULA canopy is added.
    Using an averaged value of the surface of all struts (0.687 m\textsuperscript{2}), a single strut can be added to the system to approximate the weight. This weight is multiplied by a safety factor of 2, because the strut is thicker than the canopy. A 5% margin on the weight will be taken, to keep possible inaccuracies of the used data into account.
    The total weight and the individual contributions of the canopy and the airframe can be seen in ... ."""
    return Kite_weight_ALUULA, Tether_diameter, Tether_mass, Load, D_drum

if __name__ == "__main__":
    structures_calculation(kite_area_in, avg_strut_in, len_drum_in, angle_in, extra_len_in, nom_load_in, saf_fac_in, kite_mass_margin, t_out, t_in)