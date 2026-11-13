"""Turbulence closure, surface layer, forcing and damping, initial state."""

from ..enums import enum_values
from ..options import PBL

SOUNDING_FILE = "input_sounding"
ROUGHNESS_FILE = "roughness.txt"
SFC_FILE = "sfc.txt"
GEO_TABLE_FILE = "geo_wind_table.txt"

_PBL_NAME = {"myj": "MYJ", "mynn25": "MYNN25", "mynnedmf": "MYNNEDMF", "ysu": "YSU",
             "ysunew": "YSUNew", "mrf": "MRF"}


def turbulence(deck, s) -> None:
    sec = deck.section("TURBULENCE AND DIFFUSION")
    sec.note("One closure per level: an LES type, or a RANS type, or a PBL scheme. A PBL scheme may be")
    sec.note("combined with Smagorinsky2D for the horizontal mixing only.")
    t = s.turbulence

    les = {"smagorinsky": "Smagorinsky", "smagorinsky2d": "Smagorinsky2D", "deardorff": "Deardorff"}.get(t, "None")
    sec.choice("erf.les_type", les, enum_values("LESType", ["None", "Smagorinsky", "Smagorinsky2D", "Deardorff"]),
               describe={"Smagorinsky": "needs erf.Cs",
                         "Smagorinsky2D": "horizontal mixing only; needs erf.Cs",
                         "Deardorff": "prognostic subgrid TKE; erf.Ck, erf.Ce"})
    if les in ("Smagorinsky", "Smagorinsky2D"):
        sec.set("erf.Cs", 0.17, "Smagorinsky constant; the default 0 stops the run")
        if les == "Smagorinsky2D":
            sec.alt("erf.smag2d_kh_cap", 0.0, "[m2/s] upper bound of the horizontal diffusivity; 0 = off")
    elif les == "Deardorff":
        sec.set("erf.Ck", 0.1, "eddy-viscosity coefficient (default)")
        sec.set("erf.Ce", 0.93, "dissipation coefficient (default)")
        sec.alt("erf.Ce_wall", 0.0, "extra dissipation in the first cell (default 0)")
        sec.alt("erf.theta_ref", s.theta, "[K] reference for the stability length; default: local theta")
    if t in PBL:
        sec.alt("erf.Cs", 0.25, "with erf.les_type = Smagorinsky2D")
    if les != "None" or t in PBL or t == "keqn":
        sec.alt("erf.Pr_t", 0.3333, "turbulent Prandtl number (default 1/3)")
        sec.alt("erf.Sc_t", 1.0, "turbulent Schmidt number (default 1)")
    sec.blank()

    if t == "keqn":
        sec.set("erf.rans_type", "kEqn", "one-equation k model with a diagnosed length scale")
        sec.alt("erf.Cmu0", 0.5562, "(default)")
        sec.alt("erf.Cb", 0.35, "(default)")
        sec.alt("erf.max_geom_lscale", 30.0, "[m] largest geometric length scale (default)")
        sec.alt("erf.wall_dist_type", "poisson", "or terrain_height")
        sec.alt("erf.rans_realizability", 0.0, "0 = off, else in [0.01, 10]")
    else:
        sec.alt("erf.rans_type", "kEqn", "RANS; cannot be combined with an LES type")
    sec.blank()

    pbl_values = [v for v in enum_values("PBLType", ["None", "MYJ", "MYNN25", "MYNNEDMF", "YSU", "YSUNew", "MRF"])
                  if "SHOC" not in v]
    if t in PBL:
        sec.choice("erf.pbl_type", _PBL_NAME[t], pbl_values,
                   describe={"None": "", "MYJ": "TKE scheme", "MYNN25": "TKE scheme, level 2.5",
                             "MYNNEDMF": "MYNN with mass flux", "YSU": "non-local K profile; stable conditions only",
                             "YSUNew": "YSU, revised", "MRF": "non-local K profile"})
        if t in ("mynn25", "mynnedmf"):
            sec.alt("erf.pbl_mynn_config", "NN09", "or Chen2021")
            sec.alt("erf.advect_tke", "true", "(default)")
            sec.alt("erf.diffuse_tke_3D", "true", "(default)")
        elif t in ("ysu", "ysunew"):
            sec.alt("erf.pbl_ysu_coriolis_freq", 1.0e-4, "[1/s] (default)")
            sec.alt("erf.pbl_ysu_use_consistent_coriolis", "true", "take it from erf.latitude; needs erf.use_coriolis")
            sec.alt("erf.pbl_ysu_land_Ribcr", 0.25, "critical bulk Richardson number over land (default)")
        elif t == "mrf":
            sec.alt("erf.pbl_mrf_coriolis_freq", 1.0e-4, "[1/s] (default)")
            sec.alt("erf.pbl_mrf_Ribcr", 0.5, "critical bulk Richardson number (default)")
            sec.alt("erf.pbl_mrf_const_b", 7.8, "(default)")
            sec.alt("erf.pbl_mrf_sf", 0.1, "surface-layer fraction (default)")
        sec.note("The SHOC values of erf.pbl_type need their own build options and are not listed.")
    else:
        sec.note("erf.pbl_type (default None): " + ", ".join(v for v in pbl_values if v != "None")
                 + ". All but MYJ need zlo.type = surface_layer.")
        sec.alt("erf.pbl_type", "MYNN25")
    sec.blank()

    if t == "dns":
        sec.set("erf.molec_diff_type", "ConstantAlpha", "constant kinematic coefficients")
        sec.alt("erf.molec_diff_type", "Constant", "constant dynamic coefficients")
        nu = 75.0 if s.setup == "density_current" else 1.5e-5
        sec.user("erf.dynamic_viscosity", nu, "[kg/(m s)]; nu = dynamic_viscosity / rho0_trans")
        sec.user("erf.alpha_T", nu, "[m2/s] thermal diffusivity")
        sec.alt("erf.alpha_C", nu, "[m2/s] scalar diffusivity (default 0)")
        sec.alt("erf.rho0_trans", 1.0, "[kg/m3] reference density (default 1)")
    else:
        sec.set("erf.molec_diff_type", "None", "no molecular diffusion")
        sec.alt("erf.molec_diff_type", "ConstantAlpha", "with erf.dynamic_viscosity, erf.alpha_T, erf.alpha_C")
    if s.has_tke:
        sec.alt("erf.init_tke_from_ustar", "true", "initial TKE from the surface friction velocity")



def heat_is_file(s) -> bool:
    return s.zlo_most and s.surface_heat == "file"


def surface_layer(deck, s) -> None:
    sec = deck.section("SURFACE LAYER (zlo.type = surface_layer)")
    on = s.zlo_most
    if not on:
        sec.note("Not active: zlo is a wall. To use it, set zlo.type = surface_layer and a turbulence closure.")
    flux_values = ["moeng", "rico", "bulk_coeff", "custom"]
    describe = {"moeng": "Monin-Obukhov similarity (default)",
                "custom": "needs erf.most.ustar, tstar, qstar; not with a PBL scheme",
                "bulk_coeff": "needs erf.most.Cd, Ch, Cq, pot_surf_temp, surf_moist",
                "rico": "needs erf.most.ustar, tstar, qstar; not with a PBL scheme"}
    if on:
        sec.choice("erf.surface_layer.flux_type", "moeng", flux_values, describe=describe)
    else:
        sec.alt("erf.surface_layer.flux_type", "moeng", describe["moeng"])
    sec.blank()
    if on and s.roughness == "file":
        sec.user("erf.most.roughness_file_name", ROUGHNESS_FILE,
                 "rows 'x y z0' [m], one per grid node; replace the dummy file")
        sec.alt("erf.most.z0", s.z0, "[m] one value instead of the file")
    else:
        sec.block(on, [("erf.most.z0", s.z0, "[m] roughness length (default 0.1); must be below zref")])
        sec.alt("erf.most.roughness_file_name", ROUGHNESS_FILE, "rows 'x y z0' [m], one per grid node")
    sec.alt("erf.most.zref", 0.5 * s.dz_first, "[m] reference height (default: first cell centre, this value)")
    sec.alt("erf.most.roughness_type_sea", "charnock", "over water (erf.is_land = 0): charnock, coare3.0, donelan, "
                                                       "modified_charnock, wave_coupled, constant")
    if heat_is_file(s):
        pass
    else:
        sec.alt("erf.is_land", 1, "1 land (default), 0 water")
    sec.blank()
    sec.note("Surface heating: give a temperature (pot_surf_temp) or a flux (surf_temp_flux), not both.")
    sec.note("With neither, the surface is adiabatic.")
    heat = s.surface_heat if on else "adiabatic"
    if heat == "flux":
        sec.user("erf.most.surf_temp_flux", s.heat_flux, "[K m/s] kinematic heat flux, positive upward")
    else:
        sec.alt("erf.most.surf_temp_flux", s.heat_flux, "[K m/s] kinematic heat flux, positive upward")
    if heat == "temp":
        sec.user("erf.most.pot_surf_temp", s.theta, "[K]; must equal the surface theta of the sounding")
    else:
        sec.alt("erf.most.pot_surf_temp", s.theta, "[K]; must equal the surface theta of the sounding")
    sec.alt("erf.most.surf_heating_rate", -0.25, "[K/h] with pot_surf_temp only")
    if heat == "file":
        sec.set("erf.most.use_sfc_sst", "true", "surface temperature from the time table; applied over water only")
        sec.user("erf.most.sfc_file", SFC_FILE, "header line, then 'day sst[K]'; replace the dummy file")
        sec.set("erf.is_land", 0, "water surface: ERF applies the table only where the surface is not land")
        sec.note("Over water the roughness follows erf.most.roughness_type_sea, not erf.most.z0.")
    else:
        sec.alt("erf.most.use_sfc_sst", "true", "surface temperature from a time table")
        sec.alt("erf.most.sfc_file", SFC_FILE, "header line, then 'day sst[K]' (or 'day sst H LE TAU')")
    sec.alt("erf.most.use_sfc_fluxes", "true", "fluxes from the 5-column table; not with a PBL scheme")
    sec.alt("erf.most.stable_functions", "businger_dyer", "default; or qnse (needs a prescribed heat flux)")


def forcing(deck, s) -> None:
    sec = deck.section("FORCING")
    if s.coriolis:
        sec.set("erf.use_coriolis", "true")
        sec.set("erf.latitude", s.latitude, "[deg] (default 90)")
        sec.set("erf.rotational_time_period", 86164.1, "[s] one sidereal day (default 86400)")
        sec.alt("erf.coriolis_3d", "false", "drop the cos(latitude) terms (default true)")
        if s.real:
            sec.alt("erf.variable_coriolis", "true", "Coriolis parameter from the file's latitudes (default false)")
    else:
        sec.set("erf.use_coriolis", "false", "(default)")
        sec.alt("erf.latitude", s.latitude, "[deg] with use_coriolis (default 90)")
        sec.alt("erf.rotational_time_period", 86164.1, "[s] with use_coriolis (default 86400)")
    sec.blank()
    driver = {"none": "None", "pressure_gradient": "PressureGradient"}.get(s.forcing, "GeostrophicWind")
    sec.choice("erf.abl_driver_type", driver,
               enum_values("ABLDriverType", ["None", "PressureGradient", "GeostrophicWind"]),
               describe={"None": "no large-scale driver",
                         "PressureGradient": "constant pressure gradient erf.abl_pressure_grad",
                         "GeostrophicWind": "needs erf.use_coriolis"})
    if s.forcing == "pressure_gradient":
        sec.user("erf.abl_pressure_grad", (-0.0003, 0.0, 0.0),
                 "[Pa/m] dp/dx dp/dy dp/dz; a negative x value accelerates the flow toward +x")
    else:
        sec.alt("erf.abl_pressure_grad", (-0.0003, 0.0, 0.0), "[Pa/m] with PressureGradient")
    if s.forcing == "geostrophic":
        sec.user("erf.abl_geo_wind", (s.wind[0], s.wind[1], 0.0), "[m/s] geostrophic wind")
        sec.alt("erf.abl_geo_wind_table", GEO_TABLE_FILE, "rows 'z Ug Vg'; replaces erf.abl_geo_wind")
    elif s.forcing == "geo_table":
        sec.user("erf.abl_geo_wind_table", GEO_TABLE_FILE,
                 "rows 'z Ug Vg', z increasing over the domain depth; replace the dummy file")
        sec.alt("erf.abl_geo_wind", (s.wind[0], s.wind[1], 0.0), "[m/s] one value instead of the table")
    else:
        sec.alt("erf.abl_geo_wind", (s.wind[0], s.wind[1], 0.0), "[m/s] with GeostrophicWind")
    sec.blank()
    lz = s.hi[2] - s.lo[2]
    on = s.damping == "rayleigh"
    sec.note("Rayleigh damping below the domain top%s." % ("" if on else " (not active)"))
    # an explicit relaxation: rate x time step must stay well below 2; the rate is capped at 0.25 / step
    rate_cap = float("%.2g" % (0.24 / s.dt_est)) if s.dt_est > 0.0 else 0.2
    rayleigh = min(0.2, rate_cap)
    sec.block(on, [("erf.rayleigh_damp_W", "true", "damp vertical velocity"),
                   ("erf.rayleigh_dampcoef", rayleigh,
                    "[1/s] (default 0.2); explicit: rate x time step is kept at or below 0.25 (step about %g s)" % s.dt_est),
                   ("erf.rayleigh_zdamp", 0.25 * lz, "[m] depth of the layer (default 500)")])
    sec.alt("erf.rayleigh_damp_U", "true", "also damp u toward prob.rayleigh_U_0 (default false)")
    sec.alt("erf.rayleigh_damp_V", "true", "(default false)")
    sec.alt("erf.rayleigh_damp_T", "true", "(default false)")
    sec.alt("erf.rayleigh_damping_type", "SlowExplicit", "default; or FastExplicit, FastImplicit")
    sec.blank()
    on = s.damping == "sponge"
    sec.note("Sponge zone below the domain top%s; the same keys exist for xlo, xhi, ylo, yhi, zlo." % ("" if on else " (not active)"))
    sec.block(on, [("erf.sponge_type", "Standard", "relax toward the target values below (default); or Input_Sponge"),
                   ("erf.sponge_strength", min(1.0, rate_cap),
                    "[1/s] the default 0 turns the sponge off; explicit: strength x time step is kept at or below 0.25"),
                   ("erf.use_zhi_sponge_damping", "true", ""),
                   ("erf.zhi_sponge_start", s.hi[2] - 0.25 * lz, "[m] height where the sponge starts; required with the flag"),
                   ("erf.sponge_x_velocity", s.wind[0], "[m/s] target wind in the sponge (default 0)"),
                   ("erf.sponge_y_velocity", s.wind[1], "[m/s]")])
    sec.alt("erf.sponge_density", 1.16, "[kg/m3] target density; default: the base state")
    sec.alt("erf.input_sponge_file", "input_sponge_file.txt", "rows 'z u v' with sponge_type = Input_Sponge")


def initial(deck, s) -> None:
    sec = deck.section("INITIAL STATE")
    if s.real:
        _real(sec, s)
        return
    init = {"sounding": "Input_Sounding", "uniform": "Uniform", "isentropic": "Isentropic"}[s.init]
    sec.choice("erf.init_type", init, ["Input_Sounding", "Uniform", "Isentropic"],
               describe={"Input_Sounding": "profiles of theta, qv, u, v from a text file; needs gravity",
                         "Uniform": "constant density and theta (prob.rho_0, prob.T_0); no gravity",
                         "Isentropic": "constant theta = prob.T_0 in hydrostatic balance; needs gravity"})
    others = [v for v in enum_values("InitType", []) if v not in ("None", "Input_Sounding", "Uniform", "Isentropic")]
    if others:
        sec.note("Other values: " + ", ".join(others) + ".")
    if s.init == "sounding":
        sec.user("erf.input_sounding_file", SOUNDING_FILE,
                 "line 1 'p_surf[hPa] theta[K] qv[g/kg]', then 'z theta qv u v'; replace the dummy file")
        sec.choice("erf.sounding_type", "Ideal", enum_values("SoundingType", ["ConstantDensity", "Ideal", "Isentropic", "DryIsentropic"]),
                   describe={"Ideal": "hydrostatic integration of the sounding (default)"})
    else:
        sec.alt("erf.input_sounding_file", SOUNDING_FILE, "with init_type = Input_Sounding")
    sec.blank()
    sec.note("Parameters of erf.prob_name = \"%s\"." % s.prob_name)
    if s.setup == "abl":
        _abl(sec, s)
    elif s.setup == "bubble":
        _bubble(sec, s)
    elif s.setup == "hill":
        _hill(sec, s)
    else:
        _density_current(sec, s)


WRFINPUT_FILE = "wrfinput_d01"
WRFBDY_FILE = "wrfbdy_d01"
METEM_FILES = ("met_em.d01.FIRST_TIME.nc", "met_em.d01.SECOND_TIME.nc")


def _real(sec, s) -> None:
    sec.note("Real-data initialisation needs a build with NetCDF (ERF_ENABLE_NETCDF).")
    if s.setup == "wrfinput":
        sec.set("erf.init_type", "WRFInput", "initial state, terrain and vertical levels from a WRF wrfinput file")
        sec.alt("erf.init_type", "Metgrid", "from WPS met_em files instead")
        sec.user("erf.nc_init_file_0", WRFINPUT_FILE, "the wrfinput file of level 0; not generated, supply it")
        sec.user("erf.nc_bdy_file", WRFBDY_FILE, "the wrfbdy file; without it the lateral boundaries are not real-data")
        sec.user("erf.real_width", 5, "cells of the lateral relaxation zone; at most the bdy_width of the wrfbdy file")
        sec.alt("erf.use_real_bcs", "false", "turn the real-data lateral boundaries off (default true with nc_bdy_file)")
        sec.alt("erf.nc_low_file", "wrflowinp_d01", "time-varying SST from a wrflowinp file")
        sec.alt("erf.bdy_nudge_factor", 10.0, "strength of the lateral relaxation (default)")
        sec.alt("erf.use_wrf_bdy_density", "true", "density on the boundary from the file; compressible only")
    else:
        sec.set("erf.init_type", "Metgrid", "initial and boundary data from WPS met_em files")
        sec.alt("erf.init_type", "WRFInput", "from wrfinput and wrfbdy files instead")
        sec.user("erf.nc_init_file_0", ['"%s"' % f for f in METEM_FILES],
                 "two or more met_em files in time order, evenly spaced; not generated, supply them")
        sec.user("erf.real_width", 5, "cells of the lateral relaxation zone")
        sec.alt("erf.metgrid_order", 2, "order of the vertical interpolation, 1 to 9 (default)")
        sec.alt("erf.metgrid_proximity", 500.0, "(default)")
        sec.alt("erf.metgrid_use_sfc", "true", "use the surface level of the met_em data (default)")
        sec.alt("erf.metgrid_interp_theta", "false", "interpolate theta instead of temperature (default false)")
        sec.alt("erf.write_erfbdy", "true", "save the boundary data ERF builds, for later runs (default true)")
    sec.blank()
    sec.note("Parameters of erf.prob_name = \"%s\": initial subgrid TKE only." % s.prob_name)
    sec.block(s.has_tke, [("prob.KE_0", 0.1, "[m2/s2] initial subgrid TKE (default 0)")])
    sec.alt("prob.KE_decay_height", 1000.0, "[m] TKE decays to zero at this height (default: no decay)")


def _abl(sec, s) -> None:
    if s.init == "uniform":
        sec.user("prob.rho_0", 1.16, "[kg/m3] uniform density")
        sec.user("prob.T_0", s.theta, "[K] uniform potential temperature")
        sec.user("prob.U_0", s.wind[0], "[m/s] initial wind")
        sec.user("prob.V_0", s.wind[1], "[m/s]")
    elif s.init == "isentropic":
        sec.set("prob.T_0", s.theta, "[K] potential temperature of the base state (default 300)")
        sec.user("prob.U_0", s.wind[0], "[m/s] initial wind; the isentropic base state is at rest")
        sec.user("prob.V_0", s.wind[1], "[m/s]")
    else:
        sec.note("Mean wind and temperature come from the sounding; prob.U_0 and prob.V_0 are added to it.")
        sec.alt("prob.U_0", 0.0, "[m/s] added to the sounding wind")
    if s.forcing == "none" and s.lateral == "periodic":
        sec.note("No driver is set (see FORCING): surface drag slows the wind through the run.")
    resolved = s.turbulence in ("smagorinsky", "smagorinsky2d", "deardorff", "dns", "none")
    sec.note("Random perturbations below pert_ref_height start the resolved turbulence%s."
             % ("" if resolved else " (not needed with a RANS or PBL closure)"))
    sec.block(resolved, [("prob.U_0_Pert_Mag", 0.5, "[m/s] amplitude of the u perturbation"),
                         ("prob.V_0_Pert_Mag", 0.5, "[m/s]"),
                         ("prob.W_0_Pert_Mag", 0.0, "[m/s]"),
                         ("prob.pert_ref_height", 100.0, "[m] (default)")])
    sec.alt("prob.T_0_Pert_Mag", 0.1, "[K] theta perturbation")
    sec.block(s.has_tke, [("prob.KE_0", 0.1, "[m2/s2] initial subgrid TKE (default)")])
    sec.alt("prob.KE_decay_height", 250.0, "[m] TKE decays to zero at this height (default: no decay)")


def _hill(sec, s) -> None:
    if s.init == "sounding":
        sec.note("The wind comes from the sounding; prob.U_0 and prob.V_0 are added to it.")
        sec.alt("prob.U_0", 0.0, "[m/s] added to the sounding wind")
    else:
        if s.init == "uniform":
            sec.user("prob.rho_0", 1.16, "[kg/m3] uniform density")
        sec.set("prob.T_0", s.theta, "[K] background potential temperature")
        sec.user("prob.U_0", s.wind[0], "[m/s] initial wind")
        sec.user("prob.V_0", s.wind[1], "[m/s]")


def _bubble(sec, s) -> None:
    lo, hi = s.lo, s.hi
    sec.set("prob.T_0", s.theta, "[K] background potential temperature")
    sec.user("prob.T_pert", 2.0, "[K] amplitude; negative for a cold bubble")
    sec.set("prob.T_pert_is_airtemp", "false", "false: T_pert is a theta perturbation (default true)")
    sec.user("prob.x_c", 0.5 * (lo[0] + hi[0]), "[m] bubble centre")
    sec.user("prob.y_c", 0.5 * (lo[1] + hi[1]), "[m]")
    sec.user("prob.z_c", lo[2] + 0.2 * (hi[2] - lo[2]), "[m]")
    r = 0.2 * (hi[2] - lo[2])
    sec.user("prob.x_r", r, "[m] bubble radius")
    sec.user("prob.y_r", r, "[m]")
    sec.user("prob.z_r", r, "[m]")
    if s.dycore == "anelastic":
        sec.set("prob.perturb_rho", "false", "perturb theta, not density (default true); anelastic levels keep density fixed")
    else:
        sec.alt("prob.perturb_rho", "false", "perturb theta instead of density (default true)")
    if s.at_rest:
        sec.alt("prob.U_0", 0.0, "[m/s] uniform background wind (default 0)")
    else:
        sec.user("prob.U_0", s.wind[0], "[m/s] uniform background wind")
        sec.user("prob.V_0", s.wind[1], "[m/s]")
    sec.alt("prob.do_moist_bubble", "true", "needs a moisture model (default false)")


def _density_current(sec, s) -> None:
    sec.set("prob.T_0", s.theta, "[K] background potential temperature")
    sec.alt("prob.T_pert", -15.0, "[K] amplitude of the cold bubble (default)")
    sec.alt("prob.x_c", 0.0, "[m] bubble centre (default)")
    sec.alt("prob.z_c", 3000.0, "[m] (default)")
    sec.alt("prob.x_r", 4000.0, "[m] bubble radius (default)")
    sec.alt("prob.z_r", 2000.0, "[m] (default)")
