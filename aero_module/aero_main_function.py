# -*- coding: utf-8 -*-
import numpy as np

def main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_out, AoA_range_in, V_wind, Print=False):
    import os
    import numpy as np
    import aero_module.functions_VS_LLT as VSM
    import aero_module.coordinates_generation as CGEN

    script_dir = os.path.dirname(__file__)

    ''' Code iteration part for polar graphs generation '''

    #%%  Input DATA

    tether = True   # Set parameter for tether + KCU drag estimation
    N_split = N_split

    # Definition of the geometry
    Name1 = 'Kite15'
    Atot = A_proj  # Projected area
    Segments = Kite_segments  # Number of kite segments
    Points = Points  # Number of points for coord function
    Plotting = False
    Bridles = False

    # CAD = VSM.get_Kite15_coords()  # Geometry nodes locations
    coords, chords, MAC, arclength, flat_area, flat_area_span = CGEN.Generate_Kite15_coords(Atot, Segments, Points, Plotting, Bridles, Print)
    # coords = VSM.struct2aero_geometry(CAD)  # Change geometry to definition
    N = int(len(coords) / 2)  # Number of sections defined
    dt = 0.0045   # m
    lt = 580        # m
    CDc = 1.1       # -

    # LE thickness at each section [m]
    t = 0.08*MAC*np.ones(len(chords))

    # Only for constant thickness
    Chord_sum = np.sum(chords)
    chord_diam = 0.75*t[0]
    Tube_Area1 = chords * np.pi * chord_diam + 2 * np.pi * (chord_diam / 2) ** 2
    Strut_area_av = np.average(Tube_Area1)
    if Print:
        print('Chords = ', chords)
        print('m2\nAverage strut area = ', Strut_area_av, 'm2')
        print('Thickness = ', t[0], 'm\n ------------------')

    # Camber for each section (ct in my case)
    # If you want it to vary do the same as the thickness
    k = 0.08

    #   Model and program specifics
    ring_geo = '5fil' # System of vorticity defined with 5 or 3 filaments per wing section
    model = 'VSM'     # Choose between Vortex Step Method (VSM) or Lifting Line Method (LLT)
    # Convergence criteria
    conv_crit = {
        'Niterations': 100,
        'error' : 1e-3,
        'Conv_weight': 0.03
        }


    aoa_0 = AoA_range_out
    aoa_1 = AoA_range_in
    CD_x = []
    CL_x = []

    for angle in aoa_0:

        # Wind speed vector definition
        Umag = V_wind               # Wind speed magnitude
        aoa = angle*np.pi/180      # Angle of attack
        sideslip = 0/180*np.pi  # Sideslip angle
        Uinf = np.array([np.cos(aoa)*np.cos(sideslip),np.sin(sideslip),np.sin(aoa)])*Umag # Wind speed vector
        N = int(len(coords) / 2)  # Number of sections defined
        # Plot thickness distribution
        indexes = np.empty(N)
        for i in range(N):
            indexes[i]=coords[2*i,1]
        # fig = plt.figure(figsize = (6,5))
        # plt.plot(indexes,t,'ro',label='CAD Drawing')


        # Number of splits per panel
        N_split = N_split
        # Refine structrural mesh into more panels
        coord = VSM.refine_LEI_mesh(coords, N-1, N_split)
        # Define system of vorticity
        controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_LEI(
            coord, Uinf, int((N-1)*N_split+1),ring_geo,model)
        N = int(len(coord)/2)           # Number of section after refining the mesh


        #%% Airfoil Coefficients

        # Definition of the thickness distribution for the refined mesh
        thicc = np.array([])
        for i in range(Segments-1):
            temp = np.linspace(t[i],t[i+1],N_split+1)
            temp1 = []
            for a in range(len(temp)-1):
                temp1.append((temp[a] +temp[a+1])/2)
            thicc = np.append(thicc,temp1)

        # Definition of airfoil coefficients
        # Based on Breukels (2011) correlation model
        aoas = np.arange(-20,21,1)
        data_airf = np.empty((len(aoas),4,N-1))
        t_c = np.empty(N-1)
        for i in range(N-1):
            for j in range(len(aoas)):
                t_c[i] = thicc[i]/controlpoints[i]['chord']
                alpha = aoas[j]
                Cl,Cd,Cm = VSM.LEI_airf_coeff(t_c[i], k, alpha)
                data_airf[j,0,i] = alpha
                data_airf[j,1,i] = Cl
                data_airf[j,2,i] = Cd
                data_airf[j,3,i] = Cm

        #%% SOLVER
        # Define system of vorticity
        controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_LEI(
            coord, Uinf, N,ring_geo,model)
        # Solve for Gamma
        Fmag, Gamma,aero_coeffs = VSM.solve_lifting_line_system_matrix_approach_semiinfinite(
            ringvec, controlpoints, rings,Uinf,data_airf,conv_crit,model)
        #%OUTPUT Results
        F_rel,F_gl,Ltot,Dtot,CL,CD,Lift = VSM.output_results(Fmag,aero_coeffs,ringvec,Uinf,controlpoints,Atot)

        print(angle, CL,CD)

        CD_x.append(CD)
        CL_x.append(CL)

    CL_x = np.array(CL_x)
    CD_x = np.array(CD_x)

    if tether == True:
        CDt = 1/4 * dt*lt/Atot *CDc
        CD_x = 1.1*(np.array(CD_x)+CDt)

    LD = CL_x/CD_x
    Power = CL_x**3 / CD_x**2

    CL_average = np.average(CL_x)
    CD_average = np.average(CD_x)
    CL3_CD2_average = np.average(Power)

    if Print:
        print(f'CD Tether = {CDt}\n'
              f'CD KCU = {1.1*(np.array(CD_x)+CDt)}')
        print(f'Max CL = {np.max(CL_x)}')
        print(f'Max LD = {np.max(LD)}')
        print(f'Max L3/D2 = {np.max(Power)}')
        print(f'Average LD = {np.average(LD)}')
        print(f'Average L3D2 = {CL3_CD2_average}')
        print('Average CD = ', CD_average)
        print('Average CL = ', CL_average)
        print('Camille L3D2 = ', np.average(CL_x)**3/np.average(CD_x)**2)

    CD_in = []
    CL_in = []

    for angle in aoa_1:

        # Wind speed vector definition
        Umag = V_wind               # Wind speed magnitude
        aoa = angle*np.pi/180      # Angle of attack
        sideslip = 0/180*np.pi  # Sideslip angle
        Uinf = np.array([np.cos(aoa)*np.cos(sideslip),np.sin(sideslip),np.sin(aoa)])*Umag # Wind speed vector
        N = int(len(coords) / 2)  # Number of sections defined
        # Plot thickness distribution
        indexes = np.empty(N)
        for i in range(N):
            indexes[i]=coords[2*i,1]
        # fig = plt.figure(figsize = (6,5))
        # plt.plot(indexes,t,'ro',label='CAD Drawing')


        # Number of splits per panel
        N_split = N_split
        # Refine structrural mesh into more panels
        coord = VSM.refine_LEI_mesh(coords, N-1, N_split)
        # Define system of vorticity
        controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_LEI(
            coord, Uinf, int((N-1)*N_split+1),ring_geo,model)
        N = int(len(coord)/2)           # Number of section after refining the mesh


        #%% Airfoil Coefficients

        # Definition of the thickness distribution for the refined mesh
        thicc = np.array([])
        for i in range(Segments-1):
            temp = np.linspace(t[i],t[i+1],N_split+1)
            temp1 = []
            for a in range(len(temp)-1):
                temp1.append((temp[a] +temp[a+1])/2)
            thicc = np.append(thicc,temp1)

        # Definition of airfoil coefficients
        # Based on Breukels (2011) correlation model
        aoas = np.arange(-20,21,1)
        data_airf = np.empty((len(aoas),4,N-1))
        t_c = np.empty(N-1)
        for i in range(N-1):
            for j in range(len(aoas)):
                t_c[i] = thicc[i]/controlpoints[i]['chord']
                alpha = aoas[j]
                Cl,Cd,Cm = VSM.LEI_airf_coeff(t_c[i], k, alpha)
                data_airf[j,0,i] = alpha
                data_airf[j,1,i] = Cl
                data_airf[j,2,i] = Cd
                data_airf[j,3,i] = Cm

        #%% SOLVER
        # Define system of vorticity
        controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_LEI(
            coord, Uinf, N,ring_geo,model)
        # Solve for Gamma
        Fmag, Gamma,aero_coeffs = VSM.solve_lifting_line_system_matrix_approach_semiinfinite(
            ringvec, controlpoints, rings,Uinf,data_airf,conv_crit,model)
        #%OUTPUT Results
        F_rel,F_gl,Ltot,Dtot,CL,CD,Lift = VSM.output_results(Fmag,aero_coeffs,ringvec,Uinf,controlpoints,Atot)

        print(angle, CL,CD)

        CD_in.append(CD)
        CL_in.append(CL)

    if tether == True:
        CDt = 1 / 4 * dt * lt / Atot * CDc
        CD_in = 1.1 * (np.array(CD_in) + CDt)

    CD_in = np.average(np.array(CD_in))
    print('Aero Analysis Completed')
    return CL_average, CD_average, CL3_CD2_average, CD_in, A_proj, Strut_area_av, flat_area, flat_area_span, chords

# a,b,c,d,e,f,g,h,i = main_aero_function(14.39, 1000000, 12, 5, np.arange(8, 15.5, 0.5), np.arange(-1, 3.5, 0.5), 32.5, False)
# print(a, b, c, d)

if __name__ == '__main__':
    A_proj = 17
    Points = 10000
    Kite_segments = 12
    N_split = 5
    AoA_range_out = np.arange(8, 15.5, 0.5)
    AoA_range_in = np.arange(-1, 3.5, 0.5)
    V_wind = 32.5

    CL_average, CD_average, CL3_CD2_average, CD_in, A_proj, Strut_area_av, flat_area, flat_area_span, chords = \
        main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_out, AoA_range_in, V_wind, Print=True)