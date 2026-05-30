#ifndef ERF_FORESTDRAG_H_
#define ERF_FORESTDRAG_H_

#include <memory>
#include <AMReX_MultiFab.H>

/*
  ForestDrag flow physics adapted from:
  Lalic & Mihailovic (2004)
  https://doi.org/amrex::Real(10.1175)/1520-0450(2004)043<0641:AERDLD>two.CO;2
 */
class ForestDrag
{
public:

    // Constructor for discrete patch mode (text file)
    explicit ForestDrag (std::string forestfile);

    // Constructor for gridded NetCDF mode with Cd from file
    ForestDrag (std::string lai_file,
                std::string height_file,
                std::string cd_file,
                std::string vegtype_file = "",
                int tree_type = 1,
                amrex::Real laimax = 0.8);

    // Constructor for gridded NetCDF mode with constant Cd
    ForestDrag (std::string lai_file,
                std::string height_file,
                amrex::Real cd_const,
                std::string vegtype_file = "",
                int tree_type = 1,
                amrex::Real laimax = 0.8);

    ~ForestDrag () = default;

    void
    define_drag_field (const amrex::BoxArray& ba,
                       const amrex::DistributionMapping& dm,
                       amrex::Geometry& geom,
                       amrex::MultiFab* z_phys_cc,
                       amrex::MultiFab* z_phys_nd);

    amrex::MultiFab*
    get_drag_field () { return m_forest_drag.get(); }

    amrex::MultiFab*
    get_frontal_area () { return m_frontal_area.get(); }

    amrex::MultiFab*
    get_cumLAI () { return m_cumLAI.get(); }

private:
    // Mode flag: true for gridded NetCDF, false for discrete patches
    bool m_use_gridded_data{false};

    // Constant Cd mode (skips Cd file; used when forest_cd = value is specified)
    bool m_use_const_cd{false};
    amrex::Real m_cd_const{0.0};

    // Discrete patch data (for text file mode)
    amrex::Vector<amrex::Real> m_type_forest;
    amrex::Vector<amrex::Real> m_x_forest;
    amrex::Vector<amrex::Real> m_y_forest;
    amrex::Vector<amrex::Real> m_height_forest;
    amrex::Vector<amrex::Real> m_diameter_forest;
    amrex::Vector<amrex::Real> m_cd_forest;
    amrex::Vector<amrex::Real> m_lai_forest;
    amrex::Vector<amrex::Real> m_laimax_forest;

    // Gridded data (for NetCDF mode) - stored on CPU
    amrex::Vector<amrex::Real> m_gridded_lai;
    amrex::Vector<amrex::Real> m_gridded_height;
    amrex::Vector<amrex::Real> m_gridded_cd;
    amrex::Vector<int> m_gridded_vegtype;

    // Grid metadata for gridded data
    int m_grid_nx{0};
    int m_grid_ny{0};
    amrex::Real m_grid_dx{0.0};
    amrex::Real m_grid_dy{0.0};
    amrex::Real m_grid_xmin{0.0};
    amrex::Real m_grid_ymin{0.0};

    // Tree type and LAImax for gridded data (applies to entire domain)
    int m_tree_type{1};
    amrex::Real m_laimax{0.8};

    // Output drag field: Cd * a(z)
    std::unique_ptr<amrex::MultiFab> m_forest_drag;

    // Frontal area density a(z) (without Cd factor)
    std::unique_ptr<amrex::MultiFab> m_frontal_area;

    // Cumulative LAI from canopy top down to each cell
    std::unique_ptr<amrex::MultiFab> m_cumLAI;

    // Helper: read 2D gridded data from a NetCDF file
    void read_netcdf_file (const std::string& filename,
                           const std::string& varname,
                           amrex::Vector<amrex::Real>& data,
                           int& nx, int& ny,
                           amrex::Real& dx, amrex::Real& dy,
                           amrex::Real& xmin, amrex::Real& ymin);

    amrex::Real interpolate_gridded_data (const amrex::Vector<amrex::Real>& data,
                                          amrex::Real x, amrex::Real y) const;
};
#endif
