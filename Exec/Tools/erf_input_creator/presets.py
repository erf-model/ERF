"""Named presets: one name in place of a list of options.

A preset is a set of option values (the names are the command-line options
with ``-`` written as ``_``).  ``--preset NAME`` applies them; an option given
on the command line, or in a ``--case-file``, replaces the preset's value.
Every preset is a case of the self-test, so each one is generated, checked
and, when an executable is given, run.
"""

from typing import Dict

PRESETS: Dict[str, dict] = {
    # ---- boundary layer over a flat surface ------------------------------------------------
    "abl_neutral_les": dict(
        about="neutral boundary layer, Smagorinsky LES, geostrophic wind",
        setup="abl", turbulence="smagorinsky", surface="most", coriolis=True, forcing="geostrophic",
        damping="rayleigh", output="plt,chk,logs,profiles"),
    "abl_convective_les": dict(
        about="convective boundary layer, Deardorff LES, surface heat flux 0.1 K m/s",
        setup="abl", turbulence="deardorff", surface="most", surface_heat="flux", heat_flux=0.1,
        coriolis=True, forcing="geostrophic", damping="rayleigh", output="plt,chk,logs,profiles"),
    "abl_stable_les": dict(
        about="stable boundary layer, Deardorff LES, surface cooling -0.01 K m/s",
        setup="abl", turbulence="deardorff", surface="most", surface_heat="flux", heat_flux=-0.01,
        coriolis=True, forcing="geostrophic", damping="rayleigh", output="plt,chk,logs,profiles"),
    "abl_pbl": dict(
        about="boundary layer with the MYNN level-2.5 PBL scheme on a coarse grid",
        setup="abl", turbulence="mynn25", surface="most", coriolis=True, forcing="geostrophic",
        domain=[16000.0, 16000.0, 4000.0], nx=32, ny=32, nz=40, output="plt,plt2d,chk,logs,profiles"),
    "abl_rans": dict(
        about="boundary layer with the one-equation k RANS closure",
        setup="abl", turbulence="keqn", surface="most", forcing="pressure_gradient",
        output="plt,chk,logs,profiles"),
    "abl_anelastic_les": dict(
        about="neutral boundary layer LES with the anelastic equation set",
        setup="abl", dycore="anelastic", turbulence="smagorinsky", surface="most", coriolis=True,
        forcing="geostrophic", output="plt,chk,logs,profiles"),
    "abl_inflow_outflow": dict(
        about="boundary layer with inflow at xlo and outflow at xhi",
        setup="abl", turbulence="smagorinsky", surface="most", lateral="inflow_outflow",
        output="plt,chk,logs,plane"),
    # ---- terrain, buildings, canopy, wind farm ---------------------------------------------
    "hill_fitted": dict(
        about="flow over an analytic hill on a terrain-fitted mesh",
        setup="hill", terrain_method="fitted", output="plt,chk,logs"),
    "hill_immersed": dict(
        about="flow over an analytic hill by immersed forcing, anelastic",
        setup="hill", terrain_method="immersed", dycore="anelastic", surface="most", output="plt,chk,logs"),
    "hill_eb": dict(
        about="flow over an analytic ridge as an embedded boundary, anelastic",
        setup="hill", terrain_method="eb", dycore="anelastic", output="plt,chk,logs"),
    "terrain_site_les": dict(
        about="LES over terrain from a table, stretched fitted mesh, roughness map",
        setup="abl", terrain="file", terrain_method="fitted", dz0=5.0, roughness="file", turbulence="smagorinsky",
        surface="most", coriolis=True, forcing="geostrophic", damping="rayleigh",
        output="plt,plt2d,chk,logs,profiles,stations"),
    "terrain_site_nested": dict(
        about="terrain from a table by immersed forcing with one refined level",
        setup="abl", terrain="file", terrain_method="immersed", dycore="anelastic", surface="most", amr=1,
        output="plt,chk,logs,stations"),
    "urban_les": dict(
        about="LES over buildings from a table (immersed forcing), anelastic",
        setup="abl", buildings="file", dycore="anelastic", turbulence="smagorinsky", surface="most",
        output="plt,chk,logs,plane"),
    "canopy_rans": dict(
        about="forest canopy drag with the k RANS closure",
        setup="abl", canopy="file", turbulence="keqn", surface="most", forcing="pressure_gradient",
        output="plt,chk,logs,profiles"),
    "windfarm_fitch": dict(
        about="Fitch wind-farm parameterisation with the MYNN PBL scheme (needs a wind-farm build)",
        setup="abl", windfarm="fitch", turbulence="mynn25", surface="most", coriolis=True, forcing="geostrophic",
        domain=[16000.0, 16000.0, 4000.0], nx=32, ny=32, nz=40, output="plt,chk,logs"),
    # ---- buoyant and moist idealised cases -------------------------------------------------
    "warm_bubble": dict(
        about="dry warm bubble in a neutral atmosphere",
        setup="bubble", output="plt,chk,logs"),
    "moist_bubble": dict(
        about="warm bubble with Kessler microphysics",
        setup="bubble", moisture="kessler", output="plt,chk,logs"),
    "density_current": dict(
        about="Straka density current, compressible",
        setup="density_current", output="plt,chk,logs"),
    "moist_abl_radiation": dict(
        about="moist boundary layer with two-stream radiation and the simple land model",
        setup="abl", turbulence="smagorinsky", surface="most", moisture="kessler", radiation="twostream",
        lsm="slm", output="plt,plt2d,chk,logs,profiles"),
    # ---- real-data cases ---------------------------------------------------------------------
    "real_wrf": dict(
        about="real-data case from wrfinput and wrfbdy: MYNN, Kessler",
        setup="wrfinput", output="plt,plt2d,chk,logs"),
    "real_wrf_full": dict(
        about="real-data case from wrfinput and wrfbdy: MRF, Morrison, RRTMGP, Noah-MP",
        setup="wrfinput", turbulence="mrf", moisture="morrison", radiation="rrtmgp", lsm="noahmp",
        output="plt,plt2d,chk,logs,stations"),
    "real_metgrid": dict(
        about="real-data case from WPS met_em files: MYNN, Kessler",
        setup="metgrid", output="plt,plt2d,chk,logs"),
}


def preset_lines():
    """One line per preset for --list."""
    width = max(len(name) for name in PRESETS)
    return ["    %-*s  %s" % (width, name, values["about"]) for name, values in PRESETS.items()]
