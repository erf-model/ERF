"""Phenomena the advisor knows, with a best-guess set-up for each.

Every number here is a best guess from common mesoscale-to-LES downscaling practice, not a
value taken from a paper.  Each entry gives:

  target_dx     finest horizontal cell size [m]
  outer         smallest outer (level-0) domain, x by y [km]
  finest        size of the finest refined box around the point of interest, x by y [km]
  fetch         extra length added upwind of the finest box, for turbulence to develop [km]
  top           model top [km]
  dz0           first-cell height [m]
  dz_top        largest cell height allowed at the top [m]
  hours         run length, spin-up included [h]
  zi            typical boundary-layer depth [m], for the gray-zone check
  umax          largest wind speed expected [m/s], for the time step
  moisture      erf.moisture_model
  notes         what the outer domain must contain, and the main risks
  fire          whether --fire is on by default (an erf.fire.* block in the deck); any
                phenomenon can take --fire, and --no-fire turns it off
"""

from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class Phenomenon:
    key: str
    title: str
    target_dx: float
    outer: Tuple[float, float]
    finest: Tuple[float, float]
    fetch: float
    top: float
    dz0: float
    dz_top: float
    hours: float
    zi: float
    umax: float
    moisture: str
    notes: str
    fire: bool = False


# Target cell sizes below 1000 m are 900 m divided by powers of 3 (300, 100, 33.3, 3.7 m), so
# that the default ratio-3 chain from a 900 m base reaches them exactly.  hurricane_landfall
# (1000 m) runs on its base level alone.
PHENOMENA = {p.key: p for p in [
    Phenomenon("hurricane_landfall", "Hurricane or tropical storm landfall",
               1000.0, (1000.0, 1000.0), (400.0, 400.0), 0.0, 20.0, 50.0, 600.0, 24.0, 1000.0, 70.0,
               "Morrison",
               "The storm and its rain bands must stay inside the outer domain for the whole run. "
               "HRRR covers only CONUS and its coastal waters; a storm offshore may need GFS instead."),
    Phenomenon("deep_convection", "Deep convection, squall line, supercell",
               300.0, (400.0, 400.0), (150.0, 150.0), 0.0, 20.0, 40.0, 500.0, 12.0, 1500.0, 50.0,
               "Morrison",
               "Updrafts are only partly resolved at 300 m. The storm must not reach the outer "
               "relaxation zone during the run."),
    Phenomenon("lake_effect_snow", "Lake-effect snow band",
               300.0, (400.0, 300.0), (150.0, 100.0), 0.0, 12.0, 30.0, 500.0, 24.0, 1000.0, 40.0,
               "Morrison",
               "The whole lake fetch along the wind must be inside the outer domain. HRRR lake "
               "surface temperature is coarse."),
    Phenomenon("sea_breeze", "Sea or lake breeze",
               300.0, (200.0, 200.0), (80.0, 80.0), 0.0, 6.0, 20.0, 400.0, 12.0, 1000.0, 20.0,
               "Kessler",
               "The coastline and about 50 km on each side of it must be inside the outer domain. "
               "Check the land mask and the sea-surface temperature."),
    Phenomenon("downslope_windstorm", "Downslope windstorm, gap wind (Santa Ana, Diablo, chinook)",
               100.0, (300.0, 200.0), (40.0, 40.0), 10.0, 20.0, 20.0, 600.0, 12.0, 1000.0, 60.0,
               "Kessler",
               "The whole ridge and the air upstream of it must be inside the outer domain. Lee "
               "slopes are steep; a deep top with damping keeps waves from reflecting."),
    Phenomenon("mountain_waves", "Mountain waves aloft (aviation turbulence)",
               300.0, (400.0, 300.0), (150.0, 100.0), 0.0, 25.0, 50.0, 600.0, 12.0, 1000.0, 70.0,
               "Kessler",
               "Waves need a deep domain and a damping layer that absorbs them."),
    Phenomenon("haboob_dust", "Haboob, dust storm",
               300.0, (300.0, 300.0), (120.0, 120.0), 0.0, 12.0, 20.0, 500.0, 12.0, 2000.0, 40.0,
               "Morrison",
               "The parent storm must be inside the outer domain; its timing in HRRR sets the "
               "timing of the outflow."),
    Phenomenon("valley_cold_pool", "Valley cold pool, drainage flow",
               900.0 / 27, (80.0, 80.0), (15.0, 15.0), 0.0, 8.0, 5.0, 400.0, 18.0, 300.0, 20.0,
               "Kessler",
               "The steepest-terrain case: valley walls are often too steep for a terrain-following "
               "grid at 33 m. Start the run before sunset."),
    Phenomenon("convective_bl_les", "Convective boundary layer, LES over flat land",
               100.0, (60.0, 60.0), (30.0, 30.0), 15.0, 6.0, 20.0, 400.0, 10.0, 1500.0, 30.0,
               "Kessler",
               "Resolved turbulence needs an upwind fetch to develop; the level in the gray zone "
               "is a transit level and should stay small."),
    Phenomenon("wildfire", "Wildfire (atmosphere around a fire)",
               900.0 / 27, (150.0, 150.0), (15.0, 15.0), 10.0, 10.0, 10.0, 500.0, 12.0, 1500.0, 40.0,
               "Kessler",
               "The finest box covers the fire perimeter plus its run. A Santa Ana or Diablo event "
               "also needs the outer domain of downslope_windstorm. The deck gets an erf.fire.* block; "
               "ignition and fuel keys belong to the fire engine and are left to the user.",
               fire=True),
    Phenomenon("wind_farm", "Wind farm wakes (actuator disk)",
               900.0 / 27, (100.0, 100.0), (25.0, 20.0), 15.0, 6.0, 10.0, 400.0, 12.0, 1000.0, 30.0,
               "Kessler",
               "The finest box covers the farm, its wake and an upwind fetch. Wind-farm keys need a "
               "build with ERF_ENABLE_WINDFARM and are not written."),
    Phenomenon("hill_flow", "Flow over an isolated hill or ridge (Askervein type)",
               900.0 / 27, (80.0, 80.0), (8.0, 8.0), 5.0, 6.0, 3.0, 400.0, 6.0, 1000.0, 30.0,
               "Kessler",
               "Slopes at fine resolution decide whether the terrain-following grid holds."),
    Phenomenon("downburst_gust", "Downburst gusts (for example on power lines)",
               900.0 / 27, (150.0, 150.0), (15.0, 15.0), 0.0, 15.0, 5.0, 500.0, 6.0, 1500.0, 50.0,
               "Morrison",
               "HRRR cannot place a downburst exactly in space and time; treat the run as one "
               "sample of possible outflows."),
    Phenomenon("urban_canopy", "Urban canopy, streets and buildings",
               900.0 / 243, (100.0, 100.0), (3.0, 3.0), 3.0, 4.0, 2.0, 400.0, 3.0, 1000.0, 20.0,
               "Kessler",
               "Five refined levels; a precursor or cycling run may be cheaper. Building keys are "
               "not written."),
    Phenomenon("fog_stratus", "Fog, low stratus",
               300.0, (150.0, 150.0), (60.0, 60.0), 0.0, 4.0, 5.0, 300.0, 18.0, 300.0, 20.0,
               "Morrison",
               "Vertical resolution near the ground matters more than horizontal resolution."),
]}


def names():
    return sorted(PHENOMENA)


def table_lines():
    """One line per phenomenon for --list."""
    out = ["%-20s %-7s %-13s %-11s %-6s %-5s %s" % ("phenomenon", "dx [m]", "outer [km]",
                                                     "finest [km]", "top", "fire", "title")]
    for key in names():
        p = PHENOMENA[key]
        out.append("%-20s %-7.4g %-13s %-11s %-6s %-5s %s" % (
            key, p.target_dx, "%gx%g" % p.outer, "%gx%g" % p.finest, "%g km" % p.top,
            "yes" if p.fire else "-", p.title))
    return out
