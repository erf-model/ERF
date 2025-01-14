
# Soonpil Kang (kang18@llnl.gov)
# This script creates the terrain input file for WitchOfAgnesi test

import numpy as np
import matplotlib.pyplot as plt

#------------------------------------
# Inputs
#------------------------------------

nx = 64
ny = 32

x_lo = -72000.
x_hi = 72000.
y_lo = 0.
y_hi = 1.

case = "WitchOfAgnesi"
# case = "Cylinder"
# case = "Hemisphere"
# case = "Gaussian"

#------------------------------------
# Create x and y arrays
#------------------------------------

npoints_x = nx+1
npoints_y = ny+1

x_arr = np.linspace(x_lo, x_hi, npoints_x)
y_arr = np.linspace(y_lo, y_hi, npoints_y)

#------------------------------------
# Create z-height
#------------------------------------

z_arr = np.zeros((npoints_x,npoints_y),dtype=float)

if case == "WitchOfAgnesi":

    a    = 0.5
    num  = 8 * a * a * a #8 * a * a * a;
    xcen = 0.5 * (x_lo + x_hi);
    for i in range(npoints_x):
        x = x_arr[i] - xcen
        for j in range(npoints_y):
            z_arr[i,j] = num / (x*x + 4 * a * a)

elif case == "Cylinder":

    a    = 0.5
    xcen = 0.5 * (x_lo + x_hi);
    for i in range(npoints_x):
        x = x_arr[i] - xcen
        for j in range(npoints_y):
            if abs(x)<a:
                z_arr[i,j] = np.sqrt(a**2 - x**2)
            else:
                z_arr[i,j] = 0.0

elif case == "Hemisphere":

    a    = 0.5
    xcen = 0.5 * (x_lo + x_hi);
    ycen = 0.5 * (y_lo + y_hi);
    for i in range(npoints_x):
        x = x_arr[i] - xcen
        for j in range(npoints_y):
            y = y_arr[j] - ycen
            if np.sqrt(x**2+y**2)<a:
                z_arr[i,j] = np.sqrt(a**2 - x**2 - y**2)
            else:
                z_arr[i,j] = 0.0

elif case == "Gaussian":

    a    = 0.5
    sigma_x2 = (0.25*x_hi)**2
    sigma_y2 = (0.1*y_hi)**2
    xcen = 0.25 * x_hi
    ycen = 0.75 * y_hi
    for i in range(npoints_x):
        x = x_arr[i] - xcen
        for j in range(npoints_y):
            y = y_arr[j] - ycen
            z_arr[i,j] = a * np.exp( -( x**2/2./sigma_x2 + y**2/2./sigma_y2 ) )

else:
    print("Wrong case.")
    exit()


#####################################
# Write input_terrain
#####################################

outfile = open("input_terrain.dat","w")
outfile.write( str(npoints_x) + '\n' )
outfile.write( str(npoints_y) + '\n' )
for i in range(npoints_x):
    outfile.write( str(x_arr[i]) + '\n' )
for i in range(npoints_y):
    outfile.write( str(y_arr[i]) + '\n' )
for j in range(npoints_y):
    for i in range(npoints_x):
        outfile.write( str(z_arr[i,j]) + '\n' )

#####################################
# Surface Plot
# Need to adjust the following:
#     - vmin
#    - vmax
#     - ax.set_xlim, ax.set_ylim, ax.set_zlim
#####################################

x_arr, y_arr = np.meshgrid(x_arr, y_arr)
# print(z_arr.shape)
z_reshaped = np.transpose(z_arr)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
surf = ax.plot_surface(x_arr, y_arr, z_reshaped, cmap='viridis', vmin=0, vmax=1.0)

# Add a color bar
fig.colorbar(surf)

# Set labels
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('height')

# ax.set_xlim(-80000, 80000)
# ax.set_ylim(-0.2, 1.2)
# ax.set_zlim(-0.1, 1.5)

fig.savefig('surface_plot.png')  # Saves as a PNG file

