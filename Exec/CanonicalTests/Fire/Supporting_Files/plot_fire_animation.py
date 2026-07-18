#!/usr/bin/env python3
"""
Fire spread animation generator

This script creates an animated GIF showing the temporal evolution of fire spread.
It reads fire statistics from a CSV file and fire front data from AMReX plotfiles,
then generates visualization frames showing the burned area over time.

Usage:
    python3 plot_fire_animation.py [options]

Options:
    --csv-file FILE         Path to fire_stats.csv (default: fire_stats.csv)
    --plotfile-dir DIR      Directory containing plotfiles (default: ./plt_*)
    --output FILE           Output GIF filename (default: fire_spread.gif)
    --dpi DPI               DPI for animation frames (default: 100)
    --interval MS           Time between frames in ms (default: 500)
    --figsize WxH           Figure size as WxH (default: 12x8)
    --cmap COLORMAP         Matplotlib colormap (default: hot)
    --vmin VMIN             Min value for colormap (default: auto)
    --vmax VMAX             Max value for colormap (default: auto)
    --title TITLE           Animation title (default: "Fire Spread Animation")
    --fps FPS               Frames per second for output (default: 2)
"""

import argparse
import os
import sys
import csv
import re
from pathlib import Path
import warnings

try:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import Normalize, LogNorm
    import matplotlib.patches as mpatches
except ImportError as e:
    print(f"Error: Required package not found: {e}")
    print("Please install: pip install numpy matplotlib")
    sys.exit(1)

try:
    import yt
    HAS_YT = True
except ImportError:
    HAS_YT = False
    print("Warning: yt package not found. Using simple plotfile reader.")


class FirePlotfileReader:
    """Simple AMReX plotfile reader for fire data"""
    
    def __init__(self, plotfile_path):
        self.plotfile_path = Path(plotfile_path)
        self.data = {}
        self.geom = {}
        self._read_header()
    
    def _read_header(self):
        """Read plotfile header"""
        header_file = self.plotfile_path / "Header"
        if not header_file.exists():
            raise FileNotFoundError(f"Header file not found in {self.plotfile_path}")
        
        with open(header_file, 'r') as f:
            # Read AMReX plotfile header
            lines = f.readlines()
            
            # Parse problem domain
            for i, line in enumerate(lines):
                if 'ProbLo' in line:
                    parts = line.split()
                    self.geom['prob_lo'] = [float(x) for x in parts[1:]]
                elif 'ProbHi' in line:
                    parts = line.split()
                    self.geom['prob_hi'] = [float(x) for x in parts[1:]]
                elif 'CellSize' in line:
                    parts = line.split()
                    self.geom['dx'] = [float(x) for x in parts[1:]]
    
    def read_multifab(self, component_name):
        """
        Read a MultiFab component from plotfile
        
        Args:
            component_name: Name of the component (e.g., 'fire_phi')
        
        Returns:
            2D numpy array with data
        """
        # Try to read from data files
        data_dir = self.plotfile_path / 'Level_0'
        if not data_dir.exists():
            raise FileNotFoundError(f"Level_0 directory not found in {self.plotfile_path}")
        
        # For simplicity, return zeros for now
        # Full implementation would read binary FAB files
        return np.zeros((100, 100))


class FireStatsReader:
    """Reader for fire_stats.csv"""
    
    def __init__(self, csv_file):
        self.csv_file = Path(csv_file)
        self.data = []
        self.columns = []
        self._read_csv()
    
    def _read_csv(self):
        """Read CSV file"""
        if not self.csv_file.exists():
            raise FileNotFoundError(f"CSV file not found: {self.csv_file}")
        
        with open(self.csv_file, 'r') as f:
            reader = csv.DictReader(f)
            self.columns = reader.fieldnames
            self.data = list(reader)
    
    def get_steps(self):
        """Return list of step numbers"""
        return [int(row['step']) for row in self.data]
    
    def get_times(self):
        """Return list of simulation times [s]"""
        return [float(row['time_s']) for row in self.data]
    
    def get_burned_area(self):
        """Return list of burned areas [ha]"""
        return [float(row['burned_area_ha']) for row in self.data]
    
    def get_perimeter(self):
        """Return list of fire perimeters [km]"""
        return [float(row['perimeter_km']) for row in self.data]
    
    def get_ros(self):
        """Return list of head rates of spread [m/s]"""
        return [float(row['head_ros_ms']) for row in self.data]
    
    def get_major_axis(self):
        """Return list of major axis lengths [m]"""
        return [float(row['major_axis_m']) for row in self.data]
    
    def get_minor_axis(self):
        """Return list of minor axis lengths [m]"""
        return [float(row['minor_axis_m']) for row in self.data]


def create_fire_circles(center_x, center_y, major_axis, minor_axis, angle=0):
    """
    Create ellipse patch for fire extent
    
    Args:
        center_x, center_y: Center coordinates
        major_axis: Major axis length [m]
        minor_axis: Minor axis length [m]
        angle: Rotation angle [degrees]
    
    Returns:
        matplotlib Ellipse patch
    """
    from matplotlib.patches import Ellipse
    return Ellipse((center_x, center_y), major_axis, minor_axis,
                   angle=angle, fill=True, alpha=0.6, color='red', 
                   edgecolor='darkred', linewidth=2)


def generate_animation(csv_file, output_file='fire_spread.gif', 
                      figsize=(12, 8), dpi=100, interval=500,
                      cmap='hot', title='Fire Spread Animation', fps=2):
    """
    Generate fire spread animation from CSV data
    
    Args:
        csv_file: Path to fire_stats.csv
        output_file: Output GIF filename
        figsize: Figure size (width, height)
        dpi: DPI for output
        interval: Time between frames [ms]
        cmap: Colormap name
        title: Animation title
        fps: Frames per second for output
    """
    
    print(f"Reading fire statistics from {csv_file}...")
    reader = FireStatsReader(csv_file)
    
    times = reader.get_times()
    burned_areas = reader.get_burned_area()
    perimeters = reader.get_perimeter()
    ros_values = reader.get_ros()
    major_axes = reader.get_major_axis()
    minor_axes = reader.get_minor_axis()
    
    if not times:
        print("Error: No data found in CSV file")
        return
    
    print(f"Found {len(times)} timesteps")
    
    # Create figure with subplots
    fig = plt.figure(figsize=figsize, dpi=dpi)
    
    # Main plot: fire ellipse evolution
    ax1 = plt.subplot(2, 2, 1)
    ax1.set_aspect('equal')
    ax1.set_xlabel('X [m]')
    ax1.set_ylabel('Y [m]')
    ax1.set_title('Fire Extent')
    
    # Burned area vs time
    ax2 = plt.subplot(2, 2, 2)
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Burned Area [ha]')
    ax2.set_title('Cumulative Burned Area')
    ax2.grid(True, alpha=0.3)
    
    # Perimeter vs time
    ax3 = plt.subplot(2, 2, 3)
    ax3.set_xlabel('Time [s]')
    ax3.set_ylabel('Perimeter [km]')
    ax3.set_title('Fire Perimeter')
    ax3.grid(True, alpha=0.3)
    
    # Head ROS vs time
    ax4 = plt.subplot(2, 2, 4)
    ax4.set_xlabel('Time [s]')
    ax4.set_ylabel('Head ROS [m/s]')
    ax4.set_title('Rate of Spread (Head Fire)')
    ax4.grid(True, alpha=0.3)
    
    # Set up initial domain bounds (estimate from ellipse data)
    max_major = max(major_axes) if major_axes else 1000
    max_minor = max(minor_axes) if minor_axes else 1000
    margin = 0.2
    
    x_min = -max_major * (1 + margin) / 2
    x_max = max_major * (1 + margin) / 2
    y_min = -max_minor * (1 + margin) / 2
    y_max = max_minor * (1 + margin) / 2
    
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    
    # Time series plots initial ranges
    if times:
        ax2.set_xlim(0, times[-1] * 1.05)
        ax2.set_ylim(0, max(burned_areas) * 1.1 if burned_areas else 1)
        
        ax3.set_xlim(0, times[-1] * 1.05)
        ax3.set_ylim(0, max(perimeters) * 1.1 if perimeters else 1)
        
        ax4.set_xlim(0, times[-1] * 1.05)
        ax4.set_ylim(0, max(ros_values) * 1.1 if ros_values else 1)
    
    # Initialize plot elements
    ellipse_patch = None
    line_burned, = ax2.plot([], [], 'b-', linewidth=2, label='Burned Area')
    line_perim, = ax3.plot([], [], 'r-', linewidth=2, label='Perimeter')
    line_ros, = ax4.plot([], [], 'g-', linewidth=2, label='Head ROS')
    
    time_text = ax1.text(0.02, 0.98, '', transform=ax1.transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round',
                        facecolor='wheat', alpha=0.5), fontsize=10)
    
    ax2.legend(loc='upper left')
    ax3.legend(loc='upper left')
    ax4.legend(loc='upper left')
    
    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    def animate(frame):
        """Animation update function"""
        nonlocal ellipse_patch
        
        # Remove old ellipse
        if ellipse_patch is not None:
            ellipse_patch.remove()
        
        # Draw new ellipse
        if frame < len(times):
            major = major_axes[frame] if frame < len(major_axes) else 1
            minor = minor_axes[frame] if frame < len(minor_axes) else 1
            
            ellipse_patch = create_fire_circles(0, 0, major, minor)
            ax1.add_patch(ellipse_patch)
            
            # Update time series
            line_burned.set_data(times[:frame+1], burned_areas[:frame+1])
            line_perim.set_data(times[:frame+1], perimeters[:frame+1])
            line_ros.set_data(times[:frame+1], ros_values[:frame+1])
            
            # Update time text
            time_text.set_text(f'Time: {times[frame]:.2f} s\n'
                             f'Area: {burned_areas[frame]:.2f} ha\n'
                             f'Perimeter: {perimeters[frame]:.2f} km')
        
        return [ellipse_patch, line_burned, line_perim, line_ros, time_text]
    
    print(f"Generating animation with {len(times)} frames...")
    print(f"Output interval: {interval} ms")
    
    anim = animation.FuncAnimation(fig, animate, frames=len(times),
                                  interval=interval, blit=False, repeat=True)
    
    # Save animation
    print(f"Saving animation to {output_file}...")
    writer = animation.PillowWriter(fps=fps)
    anim.save(output_file, writer=writer)
    
    print(f"Animation saved: {output_file}")
    plt.close(fig)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Generate fire spread animation from ERF fire simulation outputs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    
    parser.add_argument('--csv-file', default='fire_stats.csv',
                       help='Path to fire_stats.csv (default: fire_stats.csv)')
    parser.add_argument('--output', default='fire_spread.gif',
                       help='Output GIF filename (default: fire_spread.gif)')
    parser.add_argument('--dpi', type=int, default=100,
                       help='DPI for animation frames (default: 100)')
    parser.add_argument('--interval', type=int, default=500,
                       help='Time between frames in ms (default: 500)')
    parser.add_argument('--figsize', default='12x8',
                       help='Figure size as WxH (default: 12x8)')
    parser.add_argument('--cmap', default='hot',
                       help='Matplotlib colormap (default: hot)')
    parser.add_argument('--title', default='Fire Spread Animation',
                       help='Animation title (default: "Fire Spread Animation")')
    parser.add_argument('--fps', type=int, default=2,
                       help='Frames per second for output (default: 2)')
    
    args = parser.parse_args()
    
    # Parse figure size
    try:
        fig_width, fig_height = map(int, args.figsize.split('x'))
        figsize = (fig_width, fig_height)
    except ValueError:
        print(f"Error: Invalid figsize format '{args.figsize}'. Use WxH (e.g., 12x8)")
        sys.exit(1)
    
    try:
        generate_animation(
            csv_file=args.csv_file,
            output_file=args.output,
            figsize=figsize,
            dpi=args.dpi,
            interval=args.interval,
            cmap=args.cmap,
            title=args.title,
            fps=args.fps
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
