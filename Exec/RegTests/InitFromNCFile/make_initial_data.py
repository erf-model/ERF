###################################################
#
# A demonstration of initial condition generation
# for export to a netCDF file, and processing by
# ERF.
#
# This file and method of data generation are based
# on the original python notebook written by
# Timothy Sliwinski at CIRA/CSU/NOAA GSL.
#
###################################################


import numpy as np
import netCDF4 as nc
import math

###################################################
# Constants used in calculations.
# Cf. ERF/Source/ERF_Constants.H.
###################################################

# Constants
PI           = 3.14159265358979323846264338327950288
PIoTwo       = PI/2.0

# Physical Constants
R_d          = 287.0    # dry air constant for dry air [J/(kg-K)]
R_v          = 461.505  # water vapor constant for water vapor [J/(kg-K)]
Cp_d         = 1004.5   # We have set this so that with qv=0 we get identically gamma = 1.4
Cp_v         = 1859.0
Cp_l         = 4200.0

L_v          = 2.5e6    # latent heat of vaporization (J / kg)

p_0          = 1.0e5    # reference surface pressure [Pa]
Gamma        = 1.4      # c_p / c_v [-]
KAPPA        = 0.41     # von Karman constant
CONST_GRAV   = 9.81

# PROBLEM PARAMETERS
p_inf = 1e5  # reference pressure [Pa]
T_inf = 300.0  # reference temperature [K]
M_inf = 1.1952286093343936  # freestream Mach number [-]
alpha = 0.7853981633974483  # inflow angle, 0 --> x-aligned [rad]
beta = 1.1088514254079065  # non-dimensional max perturbation strength [-]
gamma = Gamma  # gamma = Gamma by default
R = 1.0  # characteristic length scale for grid [m]
# R = 2.0
sigma = 1.0  # Gaussian standard deviation [-]
xc = 0
yc = 0
inv_gm1 = 1.0 / (gamma - 1.0)
rho_0 = p_inf / (R_d * T_inf)
a_inf = np.sqrt(gamma * R_d * T_inf)

# Derived constants
rdOcp = R_d / Cp_d

###################################################
# Problem-specific function for calculating
# certain space-dependent values in the
# isentropic vortex
###################################################

def erf_vortex_Gaussian(x, y, xc, yc, R, beta, sigma):
    r2 = ((x-xc) * (x-xc) + (y-yc) * (y-yc)) / (R * R)
    return beta * np.exp(-r2 / (2. * sigma * sigma))


###################################################
# Problem-independent calculation of density given
# potential temperature and pressure.
# Cf. ERF/Source/Utils/ERF_EOS.H
###################################################

def getRhoThetagivenP(p, qv=0.0):
    return np.pow(p * np.pow(p_0, Gamma - 1), iGamma) * iR_d / (1.0 + R_v / R_d * qv)


###################################################
# Data for setting up the grid. These values are
# specific to the problem grid size defined in
# inputs file; cf. the field `amr.n_cell` in
# the local file `inputs`.
###################################################

# Grid shape

# Number of cells in each direction. These are associated with
# cell-centered, conserved quantities.
Nx_cell = 48
Ny_cell = 48
Nz_cell = 4

# The number of faces required for staggered grids in each
# direction. These are used for the velocity components.
Nx_face = Nx_cell + 1
Ny_face = Ny_cell + 1
Nz_face = Nz_cell + 1

# Problem geometry data. Cf. `geometry.prob_lo` and
# `geometry.prob_hi` in `inputs`.
prob_lo = np.array([-12, -12, -1])
prob_hi = np.array([12, 12, 1])

# Problem grid data.
n_cell = np.array([Nx_cell, Ny_cell, Nz_cell])

# Cell size in each direction
dx = (prob_hi - prob_lo) / n_cell


###################################################
# Populate data. Here we use numpy arrays to
# represent our discretized domain and store
# point-wise values. Note that, for clarity,
# we will calculate and store these values for
# time t = 0 in the order x, y, z.
# These will have to rearranged before export
# to netCDF, as ERF is expecting the format z, y, x
# for the spatial grid.
###################################################

# Cell centered, conserved quantities
Rho = np.ndarray(n_cell, np.float64)
RhoTheta = np.ndarray(n_cell, np.float64)
RhoScalar = np.ndarray(n_cell, np.float64)

for i in range(Rho.shape[0]):
    for j in range(Rho.shape[1]):
        for k in range(Rho.shape[2]):
            x = prob_lo[0] + (i + 0.5) * dx[0]
            y = prob_lo[1] + (j + 0.5) * dx[1]
            Omg = erf_vortex_Gaussian(x, y, xc, yc, R, beta, sigma)
            deltaT = -(gamma - 1.0) / (2.0 * sigma * sigma) * Omg * Omg
            rho_norm = (1.0 + deltaT) ** inv_gm1
            Rho[i, j, k] = rho_norm * rho_0

            T = (1.0 + deltaT) * T_inf
            p = rho_norm**Gamma / Gamma * rho_0 * a_inf * a_inf
            RhoTheta[i, j, k] = T * (p_0 / p)**rdOcp

            r2d_xy = math.sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc))
            RhoScalar[i, j, k] = 0.25 * (
                1.0 + math.cos(math.pi * min(r2d_xy, R) / R)
            ) / Rho[i, j, k]

# Staggered quantities
# x-velocity
x_vel = np.ndarray((Nx_face, Ny_cell, Nz_cell), np.float64)

for i in range(x_vel.shape[0]):
    for j in range(x_vel.shape[1]):
        for k in range(x_vel.shape[2]):
            x = prob_lo[0] + i * dx[0]
            y = prob_lo[1] + (j + 0.5) * dx[1]
            Omg = erf_vortex_Gaussian(x, y, xc, yc, R, beta, sigma)
            x_vel[i, j, k] = (M_inf * math.cos(alpha) - (y - yc) / R * Omg) * a_inf

# y-velocity
y_vel = np.ndarray((Nx_cell, Ny_face, Nz_cell), np.float64)

for i in range(y_vel.shape[0]):
    for j in range(y_vel.shape[1]):
        for k in range(y_vel.shape[2]):
            x = prob_lo[0] + (i + 0.5) * dx[0]
            y = prob_lo[1] + j * dx[1]
            Omg = erf_vortex_Gaussian(x, y, xc, yc, R, beta, sigma)
            y_vel[i, j, k] = (M_inf * math.sin(alpha) + (x - xc) / R * Omg) * a_inf

# z-velocity
z_vel = np.zeros((Nx_cell, Ny_cell, Nz_face), np.float64)


###################################################
# Populate netCDF file with the quantities
# calculated above. See documentation for the
# python packages netCDF4.
###################################################

# Init file variable
outfile = nc.Dataset("initial_data.nc", "w")

# Set up array dimensions for exported file
time_dim = outfile.createDimension("time", None)
dateStrLen_dim = outfile.createDimension("DateStrLen", 19)
bottom_top_dim = outfile.createDimension("BottomTop", Nz_cell)
bottom_top_stag_dim = outfile.createDimension("BottomTopStag", Nz_face)
south_north_dim = outfile.createDimension("SouthNorth", Ny_cell)
south_north_stag_dim = outfile.createDimension("SouthNorthStag", Ny_face)
west_east_dim = outfile.createDimension("WestEast", Nx_cell)
west_east_stag_dim = outfile.createDimension("WestEastStag", Nx_face)


dims2dcolumn = ("time", "BottomTop")

dims3dhplane = ("time", "SouthNorth", "WestEast")
dims3dhplane_ustag = ("time", "SouthNorth", "WestEastStag")
dims3dhplane_vstag = ("time", "SouthNorthStag", "WestEast")

dims4d = ("time", "BottomTop", "SouthNorth", "WestEast")
dims4d_ustag = ("time", "BottomTop", "SouthNorth", "WestEastStag")
dims4d_vstag = ("time", "BottomTop", "SouthNorthStag", "WestEast")
dims4d_wstag = ("time", "BottomTopStag", "SouthNorth", "WestEast")

# Global Attributes required by ERF
outfile.DESCRIPTION = "Python generated input file for initialization of ERF Isentropic Vortex problem. For use with file inputs_advecting on 48x48x4 grid."
outfile.SIMULATION_START_DATE = "0001-01-01_00:00:00"
outfile.DX = dx[0]
outfile.DY = dx[1]
# need to use setncattr for these because attribute name has dash ("-") in name, and python won't allow dashes for variable names
outfile.setncattr("WEST-EAST_GRID_DIMENSION", int(Nx_face))  # based on staggered grid
outfile.setncattr("SOUTH-NORTH_GRID_DIMENSION", int(Ny_face))  # based on staggered grid

# Times variable (1 single time for initialization)
times_var = outfile.createVariable("Times", "S1", ("time", "DateStrLen"))

# Follow the naming conventions for the variables that
# ERF is expecting.
Rho_var = outfile.createVariable("RHO", np.float64, dims4d)
RhoTheta_var = outfile.createVariable("T", np.float64, dims4d)
RhoScalar_var = outfile.createVariable("SCAL", np.float64, dims4d)

# Change order from (x, y, z) to (z, y, x).
Rho = np.swapaxes(Rho, 0, 2)
RhoTheta = np.swapaxes(RhoTheta, 0, 2)
RhoScalar = np.swapaxes(RhoScalar, 0, 2)

# Populate NetCDF variables for time-step 0.
Rho_var[0, :, :, :] = Rho
RhoTheta_var[0, :, :, :] = RhoTheta
RhoScalar_var[0, :, :, :] = RhoScalar

uwind_var = outfile.createVariable("U", np.float64, dims4d_ustag)
vwind_var = outfile.createVariable("V", np.float64, dims4d_vstag)
wwind_var = outfile.createVariable("W", np.float64, dims4d_wstag)

# Change order from (x, y, z) to (z, y, x).
x_vel = np.swapaxes(x_vel, 0, 2)
y_vel = np.swapaxes(y_vel, 0, 2)
z_vel = np.swapaxes(z_vel, 0, 2)

# Populate NetCDF variables for time-step 0.
uwind_var[0, :, :, :] = x_vel
vwind_var[0, :, :, :] = y_vel
wwind_var[0, :, :, :] = z_vel

outfile.close()
