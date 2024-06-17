import matplotlib.pyplot as plt
import numpy as np

from problem_definition.constraints import constraint_traction

A_range = np.arange(10,17,0.5)

# Vector of design variables
x = np.zeros((6, len(A_range)))
# A_proj    0
# gamma_in  1
# gamma_out 2
# elev_in   3
# elev_out  4
# v_wn      5

# Initial guess
initial_design = x

initial_design[0,:] = A_range
initial_design[1,:] = 0.5
initial_design[2,:] = 0.5
initial_design[3,:] = 30 * np.pi / 180
initial_design[4,:] = 70 * np.pi / 180
initial_design[5,:] = 10

traction_list = []

reel_direction = ['out', 'in']

direction = 1

for design in initial_design.T:
    # print(design)
    traction_list.append(constraint_traction(design)[direction])

plt.figure(figsize=(6, 6))
plt.plot(A_range, traction_list)
plt.xlabel('Projected Area (m^2)')
plt.ylabel('Traction Force (N)')
plt.title(f'Traction when reeling {reel_direction[direction]} for varying projected area')
plt.grid(True)
plt.show()

plt.savefig(f'../output/plots/reel_{reel_direction[direction]}.png')