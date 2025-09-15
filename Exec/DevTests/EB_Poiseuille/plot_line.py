import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Load the CSV file
csv_path = "data.csv"
df = pd.read_csv(csv_path)

dz = 0.1/32
z0 = 5.75 * dz # EB location
U_in = 1.0
mu = 1.0
H = 0.1 - z0
h = H/2

u_approx = df["x_velocity"]
zcoord_org = df["Points:2"]

# Analytical solution

u_exact = np.zeros(len(zcoord_org))
zcoord = np.zeros(len(zcoord_org))
for i in range(len(zcoord)):
    zcoord[i] = zcoord_org[i] - z0
    z = zcoord[i]
    if z < 0:
        u_exact[i] = 0
    else:
        u_exact[i] = 1.5 * U_in * (1 - ((z-h)/h)**2)


# Create a figure and axis

w = 1.5
fs = 18

plt.rc('text', usetex=True) # Activate LaTex rendering
plt.rc('font', family='serif',serif='times')

plt.rcParams.update({'font.family': 'serif'}) # set the font
plt.rcParams.update({'font.size': fs}) # set fontsize
plt.rcParams.update({'font.style': 'normal'})
plt.rcParams.update({'font.weight': 'normal'})

font_label = {'family': 'serif', 'size': fs}
font_tick = font_manager.FontProperties(family='serif', style='normal',
                                        size=fs, weight='normal', stretch='normal')

# Create the plot
plt.figure(figsize=(5, 5))

plt.plot(u_exact, zcoord, marker='none', linestyle='-', color='k', markersize=6, label='Exact', clip_on=True)
plt.plot(u_approx, zcoord, marker='none', linestyle='-', color='tab:blue', markersize=6, label='Approximate', clip_on=True)

# Set custom tick locations
plt.xticks(np.arange(0, 2, 0.2), fontproperties=font_tick)
plt.yticks(np.arange(0, 0.2, 0.01), fontproperties=font_tick)

plt.xlim(0, 1.6)
plt.ylim(0, 0.1)

plt.xlabel(r'Horizontal velocity', labelpad=5, fontdict=font_label)
plt.ylabel(r'Vertical coordinate', labelpad=5, fontdict=font_label)

legend = plt.legend(loc='upper right', prop=font_tick, borderpad=0.0)
legend.get_frame().set_linewidth(0.8)       # Optional: set border width
legend.get_frame().set_edgecolor('black')   # Optional: set border color
legend.get_frame().set_facecolor('white')   # Set white background
legend.get_frame().set_boxstyle('square')   # Set square corners (not rounded)

plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()

plt.savefig("Figure_Velocity_Profile.png", bbox_inches='tight', format='png', dpi=300)
