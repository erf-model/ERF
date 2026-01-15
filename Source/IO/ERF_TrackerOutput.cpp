#include <iomanip>

#include "ERF.H"
#include "ERF_EOS.H"
#include <filesystem>
namespace fs = std::filesystem;

using namespace amrex;

std::string
ERF::MakeVTKFilename(int nstep) {
    // Ensure output directory exists
    const std::string dir = "Output_StormTracker";
    if (!fs::exists(dir)) {
        fs::create_directory(dir);
    }

    std::ostringstream oss;
    oss << dir << "/storm_track_" << std::setw(7) << std::setfill('0') << nstep << ".vtk";
    return oss.str();
}

std::string
ERF::MakeVTKFilename_TrackerCircle(int nstep) {
    // Ensure output directory exists
    const std::string dir = "Output_StormTracker/tracker_circle";
    if (!fs::exists(dir)) {
        fs::create_directories(dir);
    }

    // Construct filename with zero-padded step
    std::ostringstream oss;
    oss << dir << "/storm_tracker_circle_" << std::setw(7) << std::setfill('0') << nstep << ".vtk";
    return oss.str();
}

std::string
ERF::MakeVTKFilename_EyeTracker_xy(int nstep) {
    // Ensure output directory exists
    const std::string dir = "Output_StormTracker/xy";
    if (!fs::exists(dir)) {
        fs::create_directories(dir);
    }

    // Construct filename with zero-padded step
    std::ostringstream oss;
    oss << dir << "/storm_track_xy_" << std::setw(7) << std::setfill('0') << nstep << ".vtk";
    return oss.str();
}

std::string
ERF::MakeFilename_EyeTracker_latlon(int nstep) {
    // Ensure output directory exists
    const std::string dir = "Output_StormTracker/latlon";
    if (!fs::exists(dir)) {
        fs::create_directories(dir);
    }

    // Construct filename with zero-padded step
    std::ostringstream oss;
    oss << dir << "/storm_track_latlon" << std::setw(7) << std::setfill('0') << nstep << ".txt";
    return oss.str();
}

std::string
ERF::MakeFilename_EyeTracker_maxvel(int nstep) {
    // Ensure output directory exists
    const std::string dir = "Output_StormTracker/maxvel";
    if (!fs::exists(dir)) {
        fs::create_directories(dir);
    }

    // Construct filename with zero-padded step
    std::ostringstream oss;
    oss << dir << "/storm_maxvel_" << std::setw(7) << std::setfill('0') << nstep << ".txt";
    return oss.str();
}

void
ERF::WriteVTKPolyline(const std::string& filename,
                      Vector<std::array<Real, 2>>& points_xy)
{
    std::ofstream vtkfile(filename);
    if (!vtkfile.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    int num_points = points_xy.size();
    if (num_points == 0) {
        vtkfile << "# vtk DataFile Version 3.0\n";
        vtkfile << "Storm Track\n";
        vtkfile << "ASCII\n";
        vtkfile << "DATASET POLYDATA\n";
        vtkfile << "POINTS " << num_points << " float\n";
        vtkfile.close();
        return;
    }
    if (num_points < 2) {
        points_xy.push_back(points_xy[0]);
    }
    num_points = points_xy.size();

    vtkfile << "# vtk DataFile Version 3.0\n";
    vtkfile << "Storm Track\n";
    vtkfile << "ASCII\n";
    vtkfile << "DATASET POLYDATA\n";

    // Write points (Z=0 assumed)
    vtkfile << "POINTS " << num_points << " float\n";
    for (const auto& pt : points_xy) {
        vtkfile << pt[0] << " " << pt[1] << " 10000.0\n";
    }

    // Write polyline connectivity
    vtkfile << "LINES 1 " << num_points + 1 << "\n";
    vtkfile << num_points << " ";
    for (int i = 0; i < num_points; ++i) {
        vtkfile << i << " ";
    }
    vtkfile << "\n";

    vtkfile.close();
}
