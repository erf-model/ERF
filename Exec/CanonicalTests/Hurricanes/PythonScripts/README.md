# Post-processing hurricane simulations

The script `PlotsForHurricaneTracking.py` visualizes hurricane simulation data by plotting storm tracks on a geographic map and generating intensity (maximum wind speed) and minimum pressure time-series plots. It supports comparing multiple data sources, such as **ERF** (Energy Research and Forecasting), **WRF** (Weather Research and Forecasting), and **Observed (Actual)** data.

## Required data

When running the hurricane simulations:

1. The data for hurricane tracks are written in the `Output_StormTracker/xy` directory.
2. The data for the intensity (i.e., maximum wind speed vs. time) are written into the `Output_StormTracker/maxvel` directory.
3. The data for the minimum pressure (i.e., minimum pressure vs. time) are written into the `Output_StormTracker/minpressure` directory.
4. The data is written at the same times as the plot files.

## Features

- **Map Visualization**: Plots storm trajectories using Cartopy with coastlines and state borders.
- **Intensity Comparison**: Generates time-series plots of maximum wind speeds.
- **Minimum Pressure Comparison**: Generates time-series plots of minimum pressure.
- **Unit Conversion**: Automatically converts ERF velocity data (assumed km/hr) to knots (via `/ 1.852`) for standardized comparison.

## Prerequisites

Ensure you have the following Python libraries installed:

- `numpy`, `matplotlib`, `scipy`
- `cartopy` (for geographic mapping)

## Command Line Options

| Option | Required | Description |
|:---|:---:|:---|
| `--area` | **Yes** | Geographic bounding box for the map. Format: `lon_min,lon_max,lat_min,lat_max` |
| `--erf_track` | No | Path to text file with ERF simulation track (Lon, Lat columns). |
| `--actual_track` | No | Path to text file with actual observed track (Lon, Lat columns). |
| `--outfile_track` | No | Filename for the saved map image (Default: `map.png`). |
| `--erf_maxvel` | No | Path to ERF velocity data (Time, Max Wind Speed columns). |
| `--actual_maxvel` | No | Path to actual observed velocity data (Time, Max Wind Speed columns). |
| `--wrf_maxvel` | No | Path to WRF simulation velocity data (Time, Max Wind Speed columns). |
| `--outfile_maxvel` | No | **Required if using maxvel flags.** Filename for the intensity plot. |
| `--erf_minpressure` | No | Path to ERF minimum pressure data (Time, Minimum Pressure columns). |
| `--actual_minpressure` | No | Path to actual observed minimum pressure data (Time, Minimum Pressure columns). |
| `--wrf_minpressure` | No | Path to WRF simulation minimum pressure data (Time, Minimum Pressure columns). |
| `--outfile_minpressure` | No | **Required if using minpressure flags.** Filename for the minimum pressure comparison plot. |

## Usage Examples

### 1. Basic Track and Intensity Plot

For example, to plot the tracks and generate a comparison of maximum wind speeds, the command will look like:

```bash
python plot_hurricane.py \
    --area=-95,-75,25,40 \
    --erf_track=erf_track.txt \
    --erf_maxvel=erf_v.txt \
    --wrf_maxvel=wrf_v.txt \
    --actual_maxvel=obs_v.txt \
    --outfile_track=map.png \
    --outfile_maxvel=intensity_comparison.png
```

### 2. Track, Intensity, and Minimum Pressure Comparison

To generate hurricane tracks, maximum wind speed comparisons, and minimum pressure comparisons:

```bash
python plot_hurricane.py \
    --area=-95,-75,25,40 \
    --erf_track=erf_track.txt \
    --erf_maxvel=erf_v.txt \
    --wrf_maxvel=wrf_v.txt \
    --actual_maxvel=obs_v.txt \
    --erf_minpressure=erf_p.txt \
    --wrf_minpressure=wrf_p.txt \
    --actual_minpressure=obs_p.txt \
    --outfile_track=map.png \
    --outfile_maxvel=intensity_comparison.png \
    --outfile_minpressure=pressure_comparison.png
```

This command generates:

- A hurricane track map (`map.png`)
- A maximum wind speed comparison plot (`intensity_comparison.png`)
- A minimum pressure comparison plot (`pressure_comparison.png`)
