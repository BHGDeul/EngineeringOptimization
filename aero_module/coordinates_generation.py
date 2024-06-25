import numpy as np
import matplotlib.pyplot as plt

def Generate_Kite15_coords(Projected_area,segments, Points_no, Plotting, Bridles, Print):
    ''' a is an ellipse's semi-major axis
        b is an ellipse's semi-minor axis '''

    Top_ellipse_ratio = 3
    Front_ellipse_ratio = 2
    AR_Projected = 1/((1/(2*Top_ellipse_ratio) * np.pi + 4/(3*Top_ellipse_ratio))/4)    # Based on 0.4 Taper ratio

    # Inputs
    Points_N = Points_no
    segments = segments
    # AR = 6                  # -
    Projected_area = Projected_area     # m^2
    Taper = 0.4             # -

    # First Calculations
    a_proj = np.sqrt(AR_Projected*Projected_area)/2
    b_proj = a_proj/Top_ellipse_ratio
    # print(b_proj*2/3)
    c_tip = 2/3 * b_proj
    c_root = c_tip/Taper

    y = []
    y2 = []

    # Front view Kite
    x = np.linspace(-1,1,Points_N)
    x = a_proj*x
    a = a_proj*1
    b = a_proj*1/Front_ellipse_ratio

    projected_span = np.max(x)-np.min(x)

    # Generate semi ellipse
    for i in x:
        if np.abs(i) <= x[-1]:
            y.append(1/Front_ellipse_ratio * np.sqrt(a**2-i**2))
            y2.append(-1/Front_ellipse_ratio * np.sqrt(a**2-i**2))
        else:
            y.append(0)

    # Introduce empty lists
    LE_x_segmented1 = []
    z_segmented = []

    arclength = 0

    for i in np.arange(1, Points_N, 1):
        arclength += np.sqrt(((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2))

    # Subdivide based on arclength
    arclist = []
    LE_x_segmented1.append(a)
    z_segmented.append(0)
    for i in range(0, int(segments), 1):
        arclist.append(arclength - i*arclength/(segments-1))

    for length in arclist:
        arc = 0
        for i in range(1, len(x), 1):
            if arc < length:
                arc += np.sqrt(((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2))
            else:
                LE_x_segmented1.append(x[i])
                z_segmented.append(y[i])
                break

    # Hardcoded Bridle System
    # Old Values
    # bridle_x_coords = [-projected_span/2 + 0.28*(projected_span/2), -projected_span/2 + 0.52*(projected_span/2),
    #                    -projected_span/2 + 0.45*(projected_span/2)]
    # bridle_z_coords = [b-1.18*b, b-0.93*b, b-1.62*b]

    # New Values
    bridle_x_coords = [-projected_span / 2 + 0.3 * (projected_span / 2),
                       -projected_span / 2 + 0.60 * (projected_span / 2),
                       -projected_span / 2 + 0.54 * (projected_span / 2),
                       0]
    bridle_z_coords = [b - 1.13 * b, b - 0.85 * b, b - 1.62 * b, -10]

    strut6_x_loc = -np.array(LE_x_segmented1)
    # strut6_x_loc = np.array([-4.71710426, -3.88796858, -2.41709853, -0.81508341,  0.81505511,  2.41707023, 3.88794028,  4.71709483])
    strut_z_loc = 1/Front_ellipse_ratio * np.sqrt(a**2-strut6_x_loc**2)

    # Bridle line matching
    under_bridles_x = np.array([bridle_x_coords[2], bridle_x_coords[0], strut6_x_loc[0], strut6_x_loc[1], bridle_x_coords[1], strut6_x_loc[2], strut6_x_loc[3], bridle_x_coords[3]])
    under_bridles_z = np.array([bridle_z_coords[2], bridle_z_coords[0], strut_z_loc[0], strut_z_loc[1], bridle_z_coords[1], strut_z_loc[2], strut_z_loc[3], bridle_z_coords[3]])
    #index_combos = np.array([[0, 1], [1, 2], [1, 3], [0, 4], [4, 5], [4, 6], [0, 7]])
    # line F3, line F7, line F6, line F2, line F5, line F4, line F1
    index_combos = np.array([[0, 7], [0, 4], [0, 1], [4, 6], [4, 5], [1, 3], [1, 2]])
    #line F1, line F2, line F3, line F4, line F5, line F6, line F7

    if Bridles and Plotting:
        plt.figure()
        list_bridle_length = []
        for i in range(0, len(index_combos)):
            x_plot = np.array([under_bridles_x[index_combos[i, 0]], under_bridles_x[index_combos[i, 1]]])
            z_plot = [under_bridles_z[index_combos[i, 0]], under_bridles_z[index_combos[i, 1]]]
            bridle_length = 2 * np.sqrt((x_plot[1] - x_plot[0]) ** 2 + (z_plot[1] - z_plot[0]) ** 2)
            list_bridle_length.append(bridle_length)
            print(f'bridle_length = {i+1, bridle_length}')
            # if Plotting:
            plt.plot(x_plot, z_plot, marker='x', c='r', linestyle='-')
            plt.plot(-x_plot, z_plot, marker='x', c='r', linestyle='-')
            plt.scatter(x_plot, z_plot, marker='x', c='r')
            plt.scatter(-x_plot, z_plot, marker='x', c='r')
            if i == 1:
                plt.plot(-x_plot, z_plot, c='r', linestyle='-', label='Bridle lines')
                plt.scatter(-x_plot, z_plot, marker='x', c='r', label='Bridle Bifurcations')
        plt.scatter(strut6_x_loc, strut_z_loc, c='r', marker='o', label='Strut locations')
        plt.ylim(-3, 2.5)
        print(f'Bridle lenght = {bridle_length}')
        print(list_bridle_length)
    if Plotting:
        plt.figure()
        plt.plot(x, y, label='Front View')
        # plt.scatter(bridle_x_coords, bridle_z_coords, c='r', marker='x')
        # plt.scatter(bridle_x_coords_sym, bridle_z_coords, c='r', marker='x')
        plt.plot(LE_x_segmented1, z_segmented, label='LE Front View Segmented')
        plt.xlabel('y [m]')
        plt.ylabel('z [m]')
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')

    span_scaling = projected_span/arclength
    flat_area_span = 1/span_scaling * projected_span
    flat_area = Projected_area/span_scaling

    # Projected top view Kite
    y = []
    y2 = []

    '''# TOP view Kite'''
    x = np.linspace(-1, 1, Points_N)
    x = a_proj * x
    a = a_proj
    b = a_proj * 1/Top_ellipse_ratio
    # MAC = 4 * b / (3 * np.pi)+c_tip
    # # print('MAC = ', MAC)
    # # print('b = ', b)
    # x_mac = 4*a/(3*np.pi)

    LE_x_segmented = []
    y_segmented1 = []

    for i in x:
        if np.abs(i) <= x[-1]:
            y.append(b/a * np.sqrt(a**2-i**2))
            y2.append(-b/a * np.sqrt(a**2-i**2))
        else:
            y.append(0)
    sum = 0
    for i in range(0, len(x) - 1):
        sum += (x[i] - x[i + 1]) * (y2[i + 1])
    MAC_1 = sum / (x[-1] - x[0]) + c_tip

    arclength1 = 0

    for i in np.arange(1, Points_N, 1):
        arclength1 += np.sqrt(((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2))

    # Subdivide based on arclength
    # arclist = []
    LE_x_segmented.append(a)
    y_segmented1.append(0)
    # for i in range(0, segments, 1):
    #     arclist.append(arclength1 - i*arclength1/(segments-1))

    y_segmented1 = (b/a * np.sqrt(a**2-np.array(LE_x_segmented1)**2))

    # for length in arclist:
    #     arc = 0
    #     for i in range(1,len(x), 1):
    #         if arc < length:
    #             arc += np.sqrt(((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2))
    #         else:
    #             LE_x_segmented.append(x[i])
    #             y_segmented1.append(y[i])
    #             break

    LE_y_segmented = -y_segmented1 - c_tip
    y_segmented = np.zeros(int(segments))
    y_segmented = np.array(y_segmented)
    z_segmented = np.array(z_segmented)


    # Convert to Uri compatible coordinates

    coord = np.empty((2*len(LE_y_segmented), 3))
    for i in range(len(LE_y_segmented)):
        coord[2*i, 0] = LE_y_segmented[i]
        coord[2*i, 1] = LE_x_segmented1[-i-1]
        coord[2*i, 2] = z_segmented[i]

        coord[2*i + 1, 0] = y_segmented[i]
        coord[2*i + 1, 1] = LE_x_segmented1[-i-1]
        coord[2*i + 1, 2] = z_segmented[i]
    chords = -np.array(LE_y_segmented)

    if Print:
        print("# ----- Kite Data Printing ----- #\nProjected AR = ", AR_Projected)
        print('Projected span = ', 2 * a_proj)
        print('Tip chord = ', c_tip)
        print('Root chord = ', c_root)
        print('arclength', arclength)
        print("Span Scaling = ", span_scaling)
        print('MAC = ', MAC_1)
        print('Flat area = ', flat_area)
        print('Flat span = ', flat_area_span)
        print('AR_flat = ', flat_area_span ** 2 / flat_area)
    print('\nCoordinate Generation Finished\n# ------------------------------ #')

    # print(coord)
    if Plotting:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(coord[::2, 0], coord[::2, 1], coord[::2, 2], marker='x', label='Leading Edge')
        ax.plot(coord[1::2, 0], coord[1::2, 1], coord[1::2, 2], marker='o', label='Trailing Edge')
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')() for a in 'xyz')])
        plt.legend()

        fig = plt.figure()
        plt.plot(LE_x_segmented1, y_segmented, linestyle='-', marker='x', label='Planform View Segmented', color='r')
        plt.plot(LE_x_segmented1, -LE_y_segmented, linestyle='-', marker='x', color='r')
        plt.plot([LE_x_segmented1[0], LE_x_segmented[0]], [y_segmented[0], -LE_y_segmented[0]], linestyle='-', marker='x', color='r')
        plt.plot([LE_x_segmented1[-1], LE_x_segmented1[-1]], [y_segmented[-1], -LE_y_segmented[-1]], linestyle='-', marker='x', color='r')
        plt.plot(x, y + c_tip, '--', label='Ellipsoid LE')
        plt.plot([2.55, 2.55], [0, MAC_1], linestyle='dotted', marker='o', label='MAC')
        plt.xlabel('y [m]')
        plt.ylabel('x [m]')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend()

        plt.show()
    return coord, chords, MAC_1, arclength, flat_area, flat_area_span


if __name__ == "__main__":
    A_proj = 16.65
    # A_proj = 21.54
    segments = 8
    plotting = False
    Bridles = False
    Print = False
    Points = 100000
    flat_area, flat_area_span, coord, chords, MAC, arclength = Generate_Kite15_coords(A_proj, segments, Points, plotting, Bridles, Print)