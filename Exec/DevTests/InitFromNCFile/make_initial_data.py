import numpy as np
import netCDF4 as nc
import math

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


# Function(s) for computed values
def erf_vortex_Gaussian(x, y, xc, yc, R, beta, sigma):
    r2 = ((x-xc) * (x-xc) + (y-yc) * (y-yc)) / (R * R)
    return beta * np.exp(-r2 / (2. * sigma * sigma))


def getRhoThetagivenP(p, qv=0.0):
    return np.pow(p * np.pow(p_0, Gamma - 1), iGamma) * iR_d / (1.0 + R_v / R_d * qv)


# Grid shape
Nx_cell = 48
Ny_cell = 48
Nz_cell = 4

Nx_face = Nx_cell + 1
Ny_face = Ny_cell + 1
Nz_face = Nz_cell + 1

prob_lo = np.array([-12, -12, -1])
prob_hi = np.array([12, 12, 1])
n_cell = np.array([Nx_cell, Ny_cell, Nz_cell])

dx = (prob_hi - prob_lo) / n_cell

# Cell center quantities
Rho_comp = np.ndarray(n_cell, np.float64)
RhoTheta_comp = np.ndarray(n_cell, np.float64)
RhoScalar_comp = np.ndarray(n_cell, np.float64)

rdOcp = R_d / Cp_d





for i in range(Rho_comp.shape[0]):
    for j in range(Rho_comp.shape[1]):
        for k in range(Rho_comp.shape[2]):
            x = prob_lo[0] + (i + 0.5) * dx[0]
            y = prob_lo[1] + (j + 0.5) * dx[1]
            Omg = erf_vortex_Gaussian(x, y, xc, yc, R, beta, sigma)
            deltaT = -(gamma - 1.0) / (2.0 * sigma * sigma) * Omg * Omg
            rho_norm = (1.0 + deltaT) ** inv_gm1
            Rho_comp[i, j, k] = rho_norm * rho_0

            T = (1.0 + deltaT) * T_inf
            p = rho_norm**Gamma / Gamma * rho_0 * a_inf * a_inf
            # rho_theta = rho_0 * rho_norm * T * (p_0 / p) ** rdOcp
            # RhoTheta_comp[i, j, k] = np.abs((rho_theta - getRhoThetagivenP(p_hse)))
            # RhoTheta_comp[i, j, k] = rho_0 * rho_norm * T * (p_0 / p)**rdOcp / Rho_comp[i, j, k]
            RhoTheta_comp[i, j, k] = T * (p_0 / p)**rdOcp


            r2d_xy = math.sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc))
            RhoScalar_comp[i, j, k] = 0.25 * (
                1.0 + math.cos(math.pi * min(r2d_xy, R) / R)
            ) / Rho_comp[i, j, k]


# TBD: These are only set if use_moisture == true
# RhoQ1_comp
# RhoQ2_comp

# x-velocity
x_vel_pert = np.ndarray((Nx_face, Ny_cell, Nz_cell), np.float64)
for i in range(x_vel_pert.shape[0]):
    for j in range(x_vel_pert.shape[1]):
        for k in range(x_vel_pert.shape[2]):
            x = prob_lo[0] + i * dx[0]
            y = prob_lo[1] + (j + 0.5) * dx[1]
            Omg = erf_vortex_Gaussian(x, y, xc, yc, R, beta, sigma)
            x_vel_pert[i, j, k] = (M_inf * math.cos(alpha) - (y - yc) / R * Omg) * a_inf



# y-velocity
y_vel_pert = np.ndarray((Nx_cell, Ny_face, Nz_cell), np.float64)
for i in range(y_vel_pert.shape[0]):
    for j in range(y_vel_pert.shape[1]):
        for k in range(y_vel_pert.shape[2]):
            x = prob_lo[0] + (i + 0.5) * dx[0]
            y = prob_lo[1] + j * dx[1]
            Omg = erf_vortex_Gaussian(x, y, xc, yc, R, beta, sigma)
            y_vel_pert[i, j, k] = (M_inf * math.sin(alpha) + (x - xc) / R * Omg) * a_inf

# z-velocity
z_vel_pert = np.zeros((Nx_cell, Ny_cell, Nz_face), np.float64)

outfile = nc.Dataset("initial_data.nc", "w")


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

uwind_var = outfile.createVariable("U", np.float64, dims4d_ustag)
vwind_var = outfile.createVariable("V", np.float64, dims4d_vstag)
wwind_var = outfile.createVariable("W", np.float64, dims4d_wstag)

x_vel_pert = np.swapaxes(x_vel_pert, 0, 2)
y_vel_pert = np.swapaxes(y_vel_pert, 0, 2)
z_vel_pert = np.swapaxes(z_vel_pert, 0, 2)

Rho_comp_var = outfile.createVariable("RHO", np.float64, dims4d)
RhoTheta_comp_var = outfile.createVariable("T", np.float64, dims4d)
RhoScalar_comp_var = outfile.createVariable("SCAL", np.float64, dims4d)


uwind_var[0, :, :, :] = x_vel_pert
vwind_var[0, :, :, :] = y_vel_pert
wwind_var[0, :, :, :] = z_vel_pert

Rho_comp = np.swapaxes(Rho_comp, 0, 2)
RhoTheta_comp = np.swapaxes(RhoTheta_comp, 0, 2)
RhoScalar_comp = np.swapaxes(RhoScalar_comp, 0, 2)

Rho_comp_var[0, :, :, :] = Rho_comp
RhoTheta_comp_var[0, :, :, :] = RhoTheta_comp
RhoScalar_comp_var[0, :, :, :] = RhoScalar_comp

outfile.close()
