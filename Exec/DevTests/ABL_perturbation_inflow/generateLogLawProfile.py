import numpy as np

# Given parameters
Re_tau = 395
kinematic_viscosity = 0.001
height_start = 1e-2
height_end = 10
#height_end = 1.0
num_points = 256
kappa = 0.41  # Von Kármán constant
B = 5.2

# Calculate friction velocity (u_tau) based on Re_tau, kinematic viscosity, and height_end
u_tau = (Re_tau * kinematic_viscosity) / height_end
print(f"u_tau = {u_tau}")

# Generate heights (y)
heights = np.linspace(height_start, height_end, num_points)

# Calculate u values using smooth wall log law of the wall formula
u = u_tau * (1. / kappa * np.log(heights * u_tau / kinematic_viscosity) + B)
print(f"u = {u}")

# Create v values (set all v values to 0)
v = np.zeros_like(heights)
w = np.zeros_like(heights)
mixing_ratio = np.zeros_like(heights)
potential_temp = np.full_like(heights, 300.0)

# Prepare data to write to file
data_sounding = np.column_stack((heights, potential_temp, mixing_ratio, u, v))
data_sponge = np.column_stack((heights, u, v))
data_inflow = np.column_stack((heights, u, v, w))

# Add the initial rows
data_sounding = np.vstack([[0.0, 300.0, 0.0, 0.0, 0.0], data_sounding])
data_sponge = np.vstack([[0.0, 0.0, 0.0], data_sponge])
data_inflow = np.vstack([[0.0, 0.0, 0.0, 0.0], data_inflow])

# Save data to file
filename = f"input_ReTau{Re_tau}Ana"
filetype = ".txt"

np.savetxt(filename + "_sounding" + filetype, data_sounding, fmt="%.8e", delimiter=' ')
print(f"Data saved to {filename}_sounding{filetype}.")
np.savetxt(filename + "_sponge" + filetype, data_sponge, fmt="%.8e", delimiter=' ')
print(f"Data saved to {filename}_sponge{filetype}.")
np.savetxt(filename + "_inflow" + filetype, data_inflow, fmt="%.8e", delimiter=' ')
print(f"Data saved to {filename}_inflow{filetype}.")
