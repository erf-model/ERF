# Post-processing hurricane simulations 

The script `PlotsForHurricaneTracking.py` visualizes hurricane simulation data by plotting storm tracks on a geographic map and generating intensity (maximum wind speed) time-series plots. It supports comparing multiple data sources, such as **ERF** (Energy Research and Forecasting), **WRF** (Weather Research and Forecasting), and **Observed (Actual)** data.

## Required data
When running the hurricane simulations
1. The data for hurricane tracks are written in the `Output_StormTracker/xy` directory. 
2. The data for the intensity (ie. max speed vs time) are written into the `Output_StormTracker/maxvel` directory.
The data is written at the same times as the plot files. 

## Features
- **Map Visualization**: Plots storm trajectories using Cartopy with coastlines and state borders.
- **Intensity Comparison**: Generates time-series plots of wind speeds.
- **Unit Conversion**: Automatically converts ERF velocity data (assumed km/hr) to knots (via `/ 1.852`) for standardized comparison.

## Prerequisites
Ensure you have the following Python libraries installed:
- `numpy`, `matplotlib`, `scipy`
- `cartopy` (for geographic mapping)

## Command Line Options

| Option | Required | Description |
|:--- |:---:|:--- |
| `--area` | **Yes** | Geographic bounding box for the map. Format: `lon_min,lon_max,lat_min,lat_max` |
| `--erf_track` | No | Path to text file with ERF simulation track (Lon, Lat columns). |
| `--actual_track` | No | Path to text file with actual observed track (Lon, Lat columns). |
| `--outfile_track` | No | Filename for the saved map image (Default: `map.png`). |
| `--erf_maxvel` | No | Path to ERF velocity data (Time, Max Wind Speed columns). |
| `--actual_maxvel` | No | Path to actual observed velocity data (Time, Max Wind Speed columns). |
| `--wrf_maxvel` | No | Path to WRF simulation velocity data (Time, Max Wind Speed columns). |
| `--outfile_maxvel` | No | **Required if using maxvel flags.** Filename for the intensity plot. |


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


