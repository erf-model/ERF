"""Moisture, radiation, land surface, forest canopy, wind farm."""

from ..enums import enum_values
from ..options import MOISTURE

FOREST_FILE = "forest.txt"
WINDFARM_LOC_FILE = "windturbines_loc_x_y.txt"
WINDFARM_SPEC_FILE = "windturbines_spec.tbl"
RRTMGP_FILES = [("erf.rrtmgp_coeffs_sw", "rrtmgp-data-sw-g224-2018-12-04.nc"),
                ("erf.rrtmgp_coeffs_lw", "rrtmgp-data-lw-g256-2018-12-04.nc"),
                ("erf.rrtmgp_cloud_optics_sw", "rrtmgp-cloud-optics-coeffs-sw.nc"),
                ("erf.rrtmgp_cloud_optics_lw", "rrtmgp-cloud-optics-coeffs-lw.nc")]

_RAIN = ("kessler", "sam", "sam_noice", "morrison", "morrison_noice", "wsm6", "wdm6")
_ICE = ("sam", "morrison", "wsm6", "wdm6")


def moist_plot_vars(s):
    if not s.moist:
        return []
    names = ["qv", "qc"]
    if s.moisture in _RAIN:
        names.append("qrain")
    if s.moisture in _ICE:
        names += ["qi", "qsnow", "qgraup"]
    return names


def moisture(deck, s) -> None:
    sec = deck.section("MOISTURE")
    chosen = MOISTURE[s.moisture]
    listed = [v for v in enum_values("MoistureType", sorted(set(MOISTURE.values())))
              if v in MOISTURE.values()]
    sec.choice("erf.moisture_model", chosen, listed,
               describe={"None": "dry", "Kessler": "warm rain: qv, qc, qrain",
                         "Kessler_NoRain": "qv, qc", "SatAdj": "saturation adjustment: qv, qc",
                         "MoistNoCondensation": "qv and qc carried, no phase change",
                         "SAM": "qv, qc, qi, rain, snow, graupel; uniform dz only",
                         "SAM_NoIce": "uniform dz only", "SAM_NoPrecip_NoIce": "uniform dz only",
                         "Morrison": "two-moment, with ice", "Morrison_NoIce": "two-moment, warm",
                         "WSM6": "single-moment, 6 classes", "WDM6": "double-moment warm rain, 6 classes"})
    sec.note("SuperDroplets and SBM need their own build options and inputs and are not listed.")
    if s.moist:
        # the default WENO3 is not available with terrain_type = EB (ERF_DataStruct.H:1817-1840)
        adv = "Upwind_3rd" if s.method == "eb" else "WENO3"
        sec.set("erf.moistscal_horiz_adv_type", adv, "moisture advection" + (" (default)" if adv == "WENO3" else "; WENO is not available with EB"))
        sec.set("erf.moistscal_vert_adv_type", adv, "")
        if s.moisture.startswith("morrison"):
            sec.alt("erf.morrison_ndcnst", 250.0, "[1/cm3] droplet number concentration (default)")
        if s.zlo_most:
            sec.alt("erf.most.surf_moist_flux", 0.0, "[kg/kg m/s] surface moisture flux (default 0)")
            sec.alt("erf.most.surf_moist", 0.01, "[kg/kg] surface mixing ratio instead of a flux")
        if s.init == "sounding":
            sec.note("The initial qv is the third column of the sounding file [g/kg]; the dummy file is dry.")
    else:
        sec.alt("erf.moistscal_horiz_adv_type", "WENO3", "moisture advection (default)")


def radiation(deck, s) -> None:
    sec = deck.section("RADIATION")
    chosen = {"none": "None", "twostream": "TwoStream", "rrtmgp": "RRTMGP"}[s.radiation]
    sec.choice("erf.radiation_model", chosen, ["None", "TwoStream", "RRTMGP"],
               describe={"None": "no radiation (default)",
                         "TwoStream": "grey two-stream scheme; no extra build option or data files",
                         "RRTMGP": "needs a build with ERF_ENABLE_RRTMGP and four data files"})
    on = s.radiation != "none"
    if on:
        sec.user("erf.rad_t_sfc", s.theta, "[K] surface temperature where no land model or surface layer supplies one; required")
    else:
        sec.alt("erf.rad_t_sfc", s.theta, "[K] required with a radiation model")
    sec.alt("erf.rad_cons_lat", 39.80986, "[deg] latitude for the solar geometry of an idealised case (default)")
    sec.alt("erf.rad_cons_lon", -98.555183, "[deg] longitude (default)")
    if on and not s.real:
        sec.user("erf.fixed_solar_zenith_angle", 0.5, "cosine of a fixed zenith angle, at most 1; or set start_datetime")
        sec.alt("start_datetime", '"2025-06-21 18:00:00"', "UTC start of the run; the sun then follows the calendar")
    else:
        sec.alt("erf.fixed_solar_zenith_angle", 0.5, "cosine of a fixed zenith angle, at most 1; default: from date and position")
    if s.radiation == "twostream":
        sec.set("erf.radiation.sw_enabled", "true", "(default)")
        sec.set("erf.radiation.lw_enabled", "true", "(default)")
        sec.alt("erf.radiation.tau_per_layer", 0.05, "shortwave optical depth per layer (default)")
        sec.alt("erf.radiation.surface_albedo_sw", 0.3, "(default)")
        sec.alt("erf.radiation.surface_emissivity_lw", 0.99, "(default)")
        sec.alt("erf.radiation.seb_enable", "true", "surface energy balance from the radiative fluxes (default false)")
        sec.note("Exec/CanonicalTests/Radiation/inputs_twostream_master lists every erf.radiation.* key.")
    elif s.radiation == "rrtmgp":
        sec.set("erf.rad_freq_in_steps", 10, "steps between radiation calls (default 1)")
        sec.user("erf.rrtmgp_file_path", ".", "directory holding the four RRTMGP data files")
        for key, name in RRTMGP_FILES:
            sec.set(key, name, "(default name)")
        sec.alt("erf.co2vmr", 388.717e-6, "CO2 volume mixing ratio (default); also o3vmr, n2ovmr, covmr, ch4vmr, o2vmr, n2vmr")


def land_surface(deck, s) -> None:
    sec = deck.section("LAND SURFACE")
    chosen = {"none": "None", "slm": "SLM", "noahmp": "NOAHMP"}[s.lsm]
    sec.choice("erf.land_surface_model", chosen, ["None", "SLM", "NOAHMP"],
               describe={"None": "no land model (default)", "SLM": "simple land model; keys under slm.",
                         "NOAHMP": "needs a build with ERF_ENABLE_NOAHMP, namelist.erf and NoahmpTable.TBL in the run directory"})
    if s.lsm == "slm":
        if s.real:
            sec.note("Real-data case: the soil layers come from the NetCDF file; leave slm.nsoil and slm.soil_dz unset.")
            sec.alt("slm.nsoil", 4, "must equal the file's layer count if given")
            sec.alt("slm.soil_dz", (0.1, 0.3, 0.6, 1.0), "[m] must equal the file's layers if given")
        else:
            sec.set("slm.nsoil", 4, "soil layers (at least 2)")
            sec.user("slm.soil_dz", (0.1, 0.3, 0.6, 1.0), "[m] layer thicknesses, top first; required for idealised cases")
        soil = [("slm.clay0", 13.0, "[%] clay content; one value, or one per layer"),
                ("slm.sand0", 17.0, "[%] sand content"),
                ("slm.sw0", 0.6, "soil wetness fraction, 0 to 1"),
                ("slm.st0", s.theta, "[K] initial soil temperature")]
        for key, value, comment in soil:
            if s.real:
                sec.alt(key, value, comment + "; idealised cases only")
            else:
                sec.user(key, value, comment)
        sec.alt("slm.landtype0", 16, "uniform land-use type (default)")
        sec.alt("slm.LAI0", 0.0, "uniform leaf area index (default)")
        sec.alt("slm.use_parameter_file", "true", "vegetation and soil parameters from slm.parameter_file (default false)")
        sec.alt("slm.parameter_file", "NoahmpTable.TBL", "(default name)")
    elif s.lsm == "noahmp":
        sec.alt("erf.lsm_nsoil", 4, "soil layers; must match the Noah-MP namelist (default)")
        sec.note("Noah-MP reads namelist.erf (&NOAHLSM_OFFLINE, ERF_SETUP_FILE_01..03) and NoahmpTable.TBL from the")
        sec.note("run directory. Templates: Exec/RegTests/NoahMP_Ideal/. Without a radiation model it runs with zero")
        sec.note("radiative forcing.")
    else:
        sec.alt("slm.soil_dz", (0.1, 0.3, 0.6, 1.0), "[m] with SLM: soil layer thicknesses; required for idealised cases")


def canopy(deck, s) -> None:
    sec = deck.section("FOREST CANOPY")
    if s.canopy == "file":
        sec.user("erf.forest_file", FOREST_FILE, "rows 'type xc yc height diameter cd lai laimax'; replace the dummy file")
    else:
        sec.alt("erf.forest_file", FOREST_FILE, "tree patches: rows 'type xc yc height diameter cd lai laimax'")
    sec.note("Gridded alternative (NetCDF build): erf.forest_lai_file, erf.forest_height_file and erf.forest_cd_file")
    sec.note("or erf.forest_cd; the two forms cannot be mixed.")
    sec.alt("erf.forest_canopy_tke", "true", "canopy source of TKE; only with erf.rans_type = kEqn (default false)")


def windfarm(deck, s) -> None:
    sec = deck.section("WIND FARM")
    chosen = {"none": "None", "fitch": "Fitch", "ewp": "EWP", "simplead": "SimpleAD"}[s.windfarm]
    sec.note("Needs a build with ERF_ENABLE_WINDFARM.")
    sec.choice("erf.windfarm_type", chosen, ["None", "Fitch", "EWP", "SimpleAD", "GeneralAD"],
               describe={"None": "(default)", "Fitch": "drag and TKE source per cell", "EWP": "explicit wake parameterisation",
                         "SimpleAD": "actuator disk with uniform thrust",
                         "GeneralAD": "blade-element actuator disk; also needs blade and airfoil tables"})
    on = s.windfarm != "none"
    sec.block(on, [("erf.windfarm_loc_type", "x_y", "turbine positions in metres; or lat_lon with erf.windfarm_x_shift, _y_shift")])
    if on:
        sec.user("erf.windfarm_loc_table", WINDFARM_LOC_FILE, "rows 'x y' [m], one per turbine; replace the dummy file")
        sec.user("erf.windfarm_spec_table", WINDFARM_SPEC_FILE, "turbine geometry and power curve; replace the dummy file")
    else:
        sec.alt("erf.windfarm_loc_table", WINDFARM_LOC_FILE, "rows 'x y' [m], one per turbine")
        sec.alt("erf.windfarm_spec_table", WINDFARM_SPEC_FILE, "turbine geometry and power curve")
    if s.windfarm == "simplead":
        sec.user("erf.sampling_distance_by_D", 2.5, "upstream sampling distance in rotor diameters; required")
        sec.user("erf.turb_disk_angle_from_x", 90.0, "[deg] disk orientation from the x axis; required")
    else:
        sec.alt("erf.sampling_distance_by_D", 2.5, "SimpleAD: upstream sampling distance in rotor diameters")
        sec.alt("erf.turb_disk_angle_from_x", 90.0, "[deg] SimpleAD and GeneralAD: disk orientation from the x axis")
