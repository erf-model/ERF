#ifndef ERF_FORESTDRAG_H_
#define ERF_FORESTDRAG_H_

#include <memory>
#include <AMReX_MultiFab.H>

/**
 * Build canopy leaf-area-density and momentum-drag fields.
 *
 * The vertical leaf-area-density profiles follow Lalic and Mihailovic (2004),
 * https://doi.org/10.1175/1520-0450(2004)043%3C0641:AERDLD%3E2.0.CO;2.
 * Canopies may be supplied as discrete patches in a text file or as gridded
 * LAI, height, and drag-coefficient NetCDF fields.
 */
class ForestDrag
{
public:

    /** Construct a discrete-patch canopy from a text file. */
    explicit ForestDrag (std::string forestfile);

    /**
     * Construct a gridded canopy with drag coefficient read from NetCDF.
     * @param lai_file NetCDF file containing the horizontal ``LAI`` field.
     * @param height_file NetCDF file containing canopy ``height``.
     * @param cd_file NetCDF file containing the ``cd`` field.
     * @param tree_type LAD profile selector (1 uniform, 2 Lalic--Mihailovic).
     * @param laimax Fractional height of maximum LAD for profile type 2.
     */
    ForestDrag (std::string lai_file,
                std::string height_file,
                std::string cd_file,
                int tree_type = 1,
                amrex::Real laimax = 0.8);

    /**
     * Construct a gridded canopy with a spatially constant drag coefficient.
     * @param lai_file NetCDF file containing the horizontal ``LAI`` field.
     * @param height_file NetCDF file containing canopy ``height``.
     * @param cd_const Spatially constant drag coefficient.
     * @param tree_type LAD profile selector (1 uniform, 2 Lalic--Mihailovic).
     * @param laimax Fractional height of maximum LAD for profile type 2.
     */
    ForestDrag (std::string lai_file,
                std::string height_file,
                amrex::Real cd_const,
                int tree_type = 1,
                amrex::Real laimax = 0.8);

    ~ForestDrag () = default;

    /** Allocate and fill drag and, optionally, frontal-area-density fields.
     *
     * @param need_frontal_area Allocate the frontal-area-density field, which
     *        is only consumed by the fixed-leaf-temperature path.
     * @param level Level index, used only in diagnostic messages.
     *
     * Neither argument is defaulted on purpose: every caller must state
     * whether the frontal-area field is needed, so that a new call site cannot
     * silently allocate and fill it when forest biophysics is off, nor report
     * a placeholder level in the coverage warning.
     */
    void
    define_drag_field (const amrex::BoxArray& ba,
                       const amrex::DistributionMapping& dm,
                       amrex::Geometry& geom,
                       amrex::MultiFab* z_phys_cc,
                       amrex::MultiFab* z_phys_nd,
                       bool need_frontal_area,
                       int level);

    /** Return the cell-centered product of drag coefficient and LAD. */
    amrex::MultiFab*
    get_drag_field () { return m_forest_drag.get(); }

    /** Return cell-centered leaf/frontal area density without drag coefficient. */
    amrex::MultiFab*
    get_frontal_area () { return m_frontal_area.get(); }

private:
    // Mode flag: true for gridded NetCDF, false for discrete patches
    bool m_use_gridded_data{false};

    // Frontal-area storage is only needed by the fixed-leaf-temperature path.
    bool m_need_frontal_area{true};

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
    // Dimensionless integral of the type-2 Lalic--Mihailovic profile.  It is
    // independent of horizontal position and canopy height, so compute it
    // once on the host rather than once per canopy cell.
    amrex::Real m_lad_normalization{1.0};

    // Output drag field: Cd * a(z)
    std::unique_ptr<amrex::MultiFab> m_forest_drag;

    // Frontal area density a(z) (without Cd factor)
    std::unique_ptr<amrex::MultiFab> m_frontal_area;

    // Helper: read gridded data from a NetCDF file
    void read_netcdf_file (const std::string& filename,
                           const std::string& varname,
                           amrex::Vector<amrex::Real>& data,
                           int& nx, int& ny,
                           amrex::Real& dx, amrex::Real& dy,
                           amrex::Real& xmin, amrex::Real& ymin);

};
#endif
