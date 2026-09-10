#include <ERF_ForestDrag.H>
#ifdef ERF_USE_NETCDF
#include <ERF_NCInterface.H>
#endif
#include <ERF_Constants.H>
#include <ERF_ForestUtils.H>
#include <ERF_GridUtils.H>
#include <AMReX_Reduce.H>

using namespace amrex;

namespace {

constexpr int lad_quadrature_points = 100;

Real
compute_lad_normalization (const Real laimax)
{
    const std::string error = erf_forest_utils::validate_laimax(
        laimax, "erf.forest_laimax");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error.empty(), error.c_str());

    const Real inv_n = Real(1.0) / Real(lad_quadrature_points);
    Real profile_sum = Real(0.0);
    for (int kk = 0; kk < lad_quadrature_points; ++kk) {
        const Real eta = Real(kk) * inv_n;
        const Real ratio = (Real(1.0) - laimax) / (Real(1.0) - eta);
        if (eta < laimax) {
            profile_sum += amrex::Math::powi<6>(ratio) *
                           std::exp(Real(6.0) * (Real(1.0) - ratio));
        } else {
            profile_sum += std::sqrt(ratio) *
                           std::exp(Real(0.5) * (Real(1.0) - ratio));
        }
    }

    const Real normalization = profile_sum * inv_n;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(normalization) && normalization > Real(0.0),
        "Lalic--Mihailovic LAD normalization must be finite and positive");
    return normalization;
}

void
warn_if_forest_grid_does_not_cover_targets (const MultiFab& target_field,
                                            const Geometry& geom,
                                            Real grid_xmin, Real grid_ymin,
                                            Real grid_dx, Real grid_dy,
                                            int grid_nx, int grid_ny,
                                            int level)
{
    ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Long, Long> reduce_data(reduce_op);
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const int source_xmax = grid_nx - 1;
    const int source_ymax = grid_ny - 1;

    // The coverage test is purely horizontal, so count each target column
    // once: reduce over a single k plane rather than the full 3-D box.  Tiles
    // partition the valid box, so restricting to the one tile that contains
    // k_ref visits every (i,j) column exactly once.  Note also that the box
    // must come from tilebox(), not validbox(): with TilingIfNotGPU() true on
    // CPU the iterator yields one entry per tile while validbox() always
    // returns the whole FAB, which would re-count the entire box per tile.
    const int k_ref = geom.Domain().smallEnd(2);

    for (MFIter mfi(target_field, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box box = mfi.tilebox();
        if (box.smallEnd(2) > k_ref || box.bigEnd(2) < k_ref) { continue; }
        box.setSmall(2, k_ref);
        box.setBig(2, k_ref);
        reduce_op.eval(box, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) -> GpuTuple<Long, Long> {
                const Real x = prob_lo[0] + (static_cast<Real>(i) + Real(0.5)) * dx[0];
                const Real y = prob_lo[1] + (static_cast<Real>(j) + Real(0.5)) * dx[1];
                const auto x_stencil = erf_grid_utils::uniform_interpolation_stencil(
                    x, grid_xmin, grid_dx, grid_nx);
                const auto y_stencil = erf_grid_utils::uniform_interpolation_stencil(
                    y, grid_ymin, grid_dy, grid_ny);
                const Long outside = (x_stencil.inside && y_stencil.inside) ? Long(0) : Long(1);
                return {Long(1), outside};
            });
    }

    const auto local = reduce_data.value();
    Long total_targets = amrex::get<0>(local);
    Long outside_targets = amrex::get<1>(local);
    ParallelDescriptor::ReduceLongSum(total_targets);
    ParallelDescriptor::ReduceLongSum(outside_targets);
    if (ParallelDescriptor::IOProcessor() && outside_targets > 0) {
        const Real source_xhi = grid_xmin + static_cast<Real>(source_xmax) * grid_dx;
        const Real source_yhi = grid_ymin + static_cast<Real>(source_ymax) * grid_dy;
        const Real target_xlo = geom.ProbLo(0) + Real(0.5) * geom.CellSize(0);
        const Real target_ylo = geom.ProbLo(1) + Real(0.5) * geom.CellSize(1);
        const Real target_xhi = geom.ProbHi(0) - Real(0.5) * geom.CellSize(0);
        const Real target_yhi = geom.ProbHi(1) - Real(0.5) * geom.CellSize(1);
        const Real percentage = total_targets > 0
            ? Real(100.0) * static_cast<Real>(outside_targets) /
              static_cast<Real>(total_targets) : Real(0.0);
        Print() << "WARNING: Forest source grid does not cover ERF target cell centers"
                << " at level " << level << ": outside " << outside_targets << " of "
                << total_targets << " (" << percentage << "%). Source extent x=["
                << grid_xmin << ", " << source_xhi << "] y=[" << grid_ymin << ", "
                << source_yhi << "]; target extent x=[" << target_xlo << ", "
                << target_xhi << "] y=[" << target_ylo << ", " << target_yhi << "].\n";
    }
}

} // namespace

/*
  Constructor to get the forest parameters:
  TreeType xc, yc, height, diameter, cd, lai, laimax
*/
ForestDrag::ForestDrag (std::string forestfile)
{
    std::ifstream file(forestfile, std::ios::in);
    if (!file.good()) {
        Abort("Cannot find forest file: " + forestfile);
    }
    // TreeType xc yc height diameter cd lai laimax
    Real value1, value2, value3, value4, value5, value6, value7, value8;
    int row_number = 0;
    while (file >> value1 >> value2 >> value3 >> value4 >> value5 >> value6 >>
           value7 >> value8) {
        ++row_number;
        if (value1 == Real(2.0)) {
            const std::string error = erf_forest_utils::validate_laimax(
                value8, "Forest file '" + forestfile + "' row " +
                std::to_string(row_number) + " column 8 (laimax)");
            if (!error.empty()) {
                Abort(error);
            }
        }
        m_type_forest.push_back(value1);
        m_x_forest.push_back(value2);
        m_y_forest.push_back(value3);
        m_height_forest.push_back(value4);
        m_diameter_forest.push_back(value5);
        m_cd_forest.push_back(value6);
        m_lai_forest.push_back(value7);
        m_laimax_forest.push_back(value8);
    }
    file.close();
}

/*
  NetCDF Constructor — gridded LAI/height from files, constant Cd
*/
ForestDrag::ForestDrag (std::string lai_file,
                        std::string height_file,
                        Real cd_const,
                        int tree_type,
                        Real laimax)
{
    m_use_gridded_data = true;
    m_use_const_cd     = true;
    m_cd_const         = cd_const;
    m_tree_type        = tree_type;
    m_laimax           = laimax;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        tree_type == 1 || tree_type == 2,
        "forest_tree_type must be 1 or 2");
    if (tree_type == 2) {
        const std::string error = erf_forest_utils::validate_laimax(
            laimax, "erf.forest_laimax");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error.empty(), error.c_str());
        m_lad_normalization = compute_lad_normalization(laimax);
    }

#ifdef ERF_USE_NETCDF
    read_netcdf_file(lai_file, "LAI", m_gridded_lai,
                     m_grid_nx, m_grid_ny,
                     m_grid_dx, m_grid_dy,
                     m_grid_xmin, m_grid_ymin);

    int nx_check, ny_check;
    Real dx_check, dy_check, xmin_check, ymin_check;
    read_netcdf_file(height_file, "height", m_gridded_height,
                     nx_check, ny_check, dx_check, dy_check,
                     xmin_check, ymin_check);

    const erf_grid_utils::UniformGridMetadata lai_grid{
        m_grid_nx, m_grid_ny, m_grid_dx, m_grid_dy, m_grid_xmin, m_grid_ymin};
    const erf_grid_utils::UniformGridMetadata height_grid{
        nx_check, ny_check, dx_check, dy_check, xmin_check, ymin_check};
    const std::string height_grid_error = erf_grid_utils::validate_matching_grid(
        lai_grid, height_grid, "forest field 'LAI' in '" + lai_file + "'",
        "forest field 'height' in '" + height_file + "'");
    if (!height_grid_error.empty()) {
        Abort(height_grid_error);
    }

    Print() << "ForestDrag: Gridded LAI/height + constant Cd=" << m_cd_const << "\n"
            << "  Grid size: " << m_grid_nx << " x " << m_grid_ny << "\n"
            << "  Grid spacing: dx=" << m_grid_dx << ", dy=" << m_grid_dy << "\n"
            << "  Tree type: " << m_tree_type << "  LAImax: " << m_laimax << "\n";

#else
    amrex::ignore_unused(lai_file, height_file);
    Abort("ERF must be compiled with NetCDF support to use gridded forest data");
#endif
}

/*
  NetCDF Constructor — gridded LAI/height/Cd all from files
*/
ForestDrag::ForestDrag (std::string lai_file,
                        std::string height_file,
                        std::string cd_file,
                        int tree_type,
                        Real laimax)
{
    m_use_gridded_data = true;
    m_tree_type = tree_type;
    m_laimax = laimax;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        tree_type == 1 || tree_type == 2,
        "forest_tree_type must be 1 or 2");
    if (tree_type == 2) {
        const std::string error = erf_forest_utils::validate_laimax(
            laimax, "erf.forest_laimax");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error.empty(), error.c_str());
        m_lad_normalization = compute_lad_normalization(laimax);
    }

#ifdef ERF_USE_NETCDF
    read_netcdf_file(lai_file, "LAI", m_gridded_lai,
                     m_grid_nx, m_grid_ny,
                     m_grid_dx, m_grid_dy,
                     m_grid_xmin, m_grid_ymin);

    int nx_check, ny_check;
    Real dx_check, dy_check, xmin_check, ymin_check;
    read_netcdf_file(height_file, "height", m_gridded_height,
                     nx_check, ny_check, dx_check, dy_check,
                     xmin_check, ymin_check);

    const erf_grid_utils::UniformGridMetadata lai_grid{
        m_grid_nx, m_grid_ny, m_grid_dx, m_grid_dy, m_grid_xmin, m_grid_ymin};
    const erf_grid_utils::UniformGridMetadata height_grid{
        nx_check, ny_check, dx_check, dy_check, xmin_check, ymin_check};
    std::string grid_error = erf_grid_utils::validate_matching_grid(
        lai_grid, height_grid, "forest field 'LAI' in '" + lai_file + "'",
        "forest field 'height' in '" + height_file + "'");
    if (!grid_error.empty()) {
        Abort(grid_error);
    }

    read_netcdf_file(cd_file, "cd", m_gridded_cd,
                     nx_check, ny_check, dx_check, dy_check,
                     xmin_check, ymin_check);

    const erf_grid_utils::UniformGridMetadata cd_grid{
        nx_check, ny_check, dx_check, dy_check, xmin_check, ymin_check};
    grid_error = erf_grid_utils::validate_matching_grid(
        lai_grid, cd_grid, "forest field 'LAI' in '" + lai_file + "'",
        "forest field 'cd' in '" + cd_file + "'");
    if (!grid_error.empty()) {
        Abort(grid_error);
    }

    Print() << "ForestDrag: Successfully read gridded forest data\n"
            << "  Grid size: " << m_grid_nx << " x " << m_grid_ny << "\n"
            << "  Grid spacing: dx=" << m_grid_dx << ", dy=" << m_grid_dy << "\n"
            << "  Domain: x=[" << m_grid_xmin << ", "
            << (m_grid_xmin + (m_grid_nx-1)*m_grid_dx) << "], y=["
            << m_grid_ymin << ", "
            << (m_grid_ymin + (m_grid_ny-1)*m_grid_dy) << "]\n"
            << "  Tree type: " << m_tree_type << "\n"
            << "  LAImax: " << m_laimax << "\n";

#else
    amrex::ignore_unused(lai_file, height_file, cd_file);
    Abort("ERF must be compiled with NetCDF support to use gridded forest data");
#endif
}

void
ForestDrag::define_drag_field (const BoxArray& ba,
                               const DistributionMapping& dm,
                               Geometry& geom,
                               MultiFab* z_phys_cc,
                               MultiFab* z_phys_nd,
                               bool need_frontal_area,
                               int level)
{
    m_need_frontal_area = need_frontal_area;
    // Geometry params
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();

    bool all_boxes_touch_bottom = true;
    for (int i = 0; i < ba.size(); i++) {
        if (ba[i].smallEnd(2) != geom.Domain().smallEnd(2)) {
            all_boxes_touch_bottom = false;
        }
    }
    AMREX_ALWAYS_ASSERT(all_boxes_touch_bottom);

    // Allocate the forest drag MF and frontal area MF
    // NOTE: 1 ghost cell for averaging to faces
    m_forest_drag.reset();
    m_forest_drag = std::make_unique<MultiFab>(ba,dm,1,1);
    m_forest_drag->setVal(zero);

    m_frontal_area.reset();
    if (m_need_frontal_area) {
        m_frontal_area = std::make_unique<MultiFab>(ba,dm,1,1);
        m_frontal_area->setVal(zero);
    }

    // Copy namespace-scope constants into automatic-storage values before
    // entering GPU lambdas.  NVCC does not make these host-side constexpr
    // variables available as device symbols automatically.
    const Real zero_d   = zero;
    const Real one_d    = one;
    const Real myhalf_d = myhalf;
    const Real fourth_d = fourth;
    const bool store_frontal_area = m_need_frontal_area;

    if (m_use_gridded_data) {
        // =====================================================================
        // Gridded NetCDF mode: interpolate from gridded LAI/height/cd data
        // =====================================================================

        const Real* lai_data_h    = m_gridded_lai.data();
        const Real* height_data_h = m_gridded_height.data();
        const Real* cd_data_h     = m_use_const_cd ? nullptr : m_gridded_cd.data();

        const int grid_nx   = m_grid_nx;
        const int grid_ny   = m_grid_ny;
        const Real grid_dx  = m_grid_dx;
        const Real grid_dy  = m_grid_dy;
        const Real grid_xmin = m_grid_xmin;
        const Real grid_ymin = m_grid_ymin;

        warn_if_forest_grid_does_not_cover_targets(
            *m_forest_drag, geom, grid_xmin, grid_ymin, grid_dx, grid_dy,
            grid_nx, grid_ny, level);

        const int grid_size = m_grid_nx * m_grid_ny;
        Gpu::DeviceVector<Real> lai_data_d(grid_size);
        Gpu::DeviceVector<Real> height_data_d(grid_size);
        Gpu::DeviceVector<Real> cd_data_d;
        if (!m_use_const_cd) {
            cd_data_d.resize(grid_size);
        }

        Gpu::copy(Gpu::hostToDevice, lai_data_h, lai_data_h + grid_size,
                  lai_data_d.begin());
        Gpu::copy(Gpu::hostToDevice, height_data_h, height_data_h + grid_size,
                  height_data_d.begin());
        if (!m_use_const_cd) {
            Gpu::copy(Gpu::hostToDevice, cd_data_h, cd_data_h + grid_size,
                      cd_data_d.begin());
        }

        const Real* lai_data    = lai_data_d.data();
        const Real* height_data = height_data_d.data();
        const Real* cd_data     = m_use_const_cd ? nullptr : cd_data_d.data();
        const Real cd_const     = m_cd_const;
        const bool use_const_cd = m_use_const_cd;

        const int tree_type = m_tree_type;
        const Real laimax   = m_laimax;
        const Real lad_normalization = m_lad_normalization;

        for (MFIter mfi(*m_forest_drag); mfi.isValid(); ++mfi) {
            Box gtbx = mfi.growntilebox();
            const Array4<Real>& levelDrag   = m_forest_drag->array(mfi);
            const Array4<Real> frontalArea = store_frontal_area
                ? m_frontal_area->array(mfi) : Array4<Real>{};
            const Array4<const Real>& z_cc  = z_phys_cc->const_array(mfi);
            const Array4<const Real>& z_nd  = z_phys_nd->const_array(mfi);

            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real x = prob_lo[0] + (i + myhalf_d) * dx[0];
                const Real y = prob_lo[1] + (j + myhalf_d) * dx[1];

                const Real z_sfc = fourth_d * (z_nd(i, j, 0) + z_nd(i + 1, j, 0) +
                                               z_nd(i, j + 1, 0) + z_nd(i + 1, j + 1, 0));
                const Real z = amrex::max((z_cc(i, j, k) - z_sfc), zero_d);

                const auto x_stencil = erf_grid_utils::uniform_interpolation_stencil(
                    x, grid_xmin, grid_dx, grid_nx);
                const auto y_stencil = erf_grid_utils::uniform_interpolation_stencil(
                    y, grid_ymin, grid_dy, grid_ny);

                if (x_stencil.inside && y_stencil.inside) {
                    const int ii_c = x_stencil.lower;
                    const int jj_c = y_stencil.lower;
                    const Real wx = x_stencil.weight;
                    const Real wy = y_stencil.weight;

                    // Bilinear interpolation of LAI
                    Real lai00 = lai_data[jj_c * grid_nx + ii_c];
                    Real lai10 = lai_data[jj_c * grid_nx + (ii_c + 1)];
                    Real lai01 = lai_data[(jj_c + 1) * grid_nx + ii_c];
                    Real lai11 = lai_data[(jj_c + 1) * grid_nx + (ii_c + 1)];
                    Real lai_interp = (lai00 * (one_d - wx) + lai10 * wx) * (one_d - wy) +
                                      (lai01 * (one_d - wx) + lai11 * wx) * wy;

                    // Bilinear interpolation of height
                    Real h00 = height_data[jj_c * grid_nx + ii_c];
                    Real h10 = height_data[jj_c * grid_nx + (ii_c + 1)];
                    Real h01 = height_data[(jj_c + 1) * grid_nx + ii_c];
                    Real h11 = height_data[(jj_c + 1) * grid_nx + (ii_c + 1)];
                    Real height_interp = (h00 * (one_d - wx) + h10 * wx) * (one_d - wy) +
                                         (h01 * (one_d - wx) + h11 * wx) * wy;

                    Real cd_interp = cd_const;
                    if (!use_const_cd) {
                        // Bilinear interpolation of Cd from file.
                        Real cd00 = cd_data[jj_c * grid_nx + ii_c];
                        Real cd10 = cd_data[jj_c * grid_nx + (ii_c + 1)];
                        Real cd01 = cd_data[(jj_c + 1) * grid_nx + ii_c];
                        Real cd11 = cd_data[(jj_c + 1) * grid_nx + (ii_c + 1)];
                        cd_interp = (cd00 * (one_d - wx) + cd10 * wx) * (one_d - wy) +
                                    (cd01 * (one_d - wx) + cd11 * wx) * wy;
                    }

                    // Compute drag if within canopy
                    if (z < height_interp && height_interp > zero_d && lai_interp > zero_d) {
                        Real af;
                        Real factor = one_d;

                        if (tree_type == 1) {
                            af = lai_interp / height_interp;
                        } else {
                            const Real treeZm = laimax * height_interp;
                            af = lai_interp / (height_interp * lad_normalization);

                            Real ratio = (height_interp - treeZm) / (height_interp - z);
                            if (z < treeZm) {
                                factor = amrex::Math::powi<6>(ratio) *
                                         std::exp(Real(6.0) * (one_d - ratio));
                            } else {
                                factor = std::sqrt(ratio) *
                                         std::exp(myhalf_d * (one_d - ratio));
                            }
                        }

                        levelDrag(i, j, k)   = cd_interp * af * factor;
                        if (store_frontal_area) {
                            frontalArea(i, j, k) = af * factor;
                        }
                    }
                }
            });
        } // mfi

        // Synchronize before DeviceVectors go out of scope
        Gpu::synchronize();

    } else {
        // =====================================================================
        // Discrete patch mode: original implementation
        // =====================================================================

        for (unsigned ii = 0; ii < m_x_forest.size(); ++ii) {
            Real af;
            Real treeZm = zero_d;
            int  tf      = int(m_type_forest[ii]);
            Real hf      = m_height_forest[ii];
            Real xf      = m_x_forest[ii];
            Real yf      = m_y_forest[ii];
            Real df      = m_diameter_forest[ii];
            Real cdf     = m_cd_forest[ii];
            Real laif    = m_lai_forest[ii];
            Real laimaxf = m_laimax_forest[ii];

            if (tf == 1) {
                af = laif / hf;
            } else {
                treeZm = laimaxf * hf;
                af = laif / (hf * compute_lad_normalization(laimaxf));
            }

            for (MFIter mfi(*m_forest_drag); mfi.isValid(); ++mfi) {
                Box gtbx = mfi.growntilebox();
                const Array4<Real>& levelDrag   = m_forest_drag->array(mfi);
                const Array4<Real> frontalArea = store_frontal_area
                    ? m_frontal_area->array(mfi) : Array4<Real>{};
                const Array4<const Real>& z_cc  = z_phys_cc->const_array(mfi);
                const Array4<const Real>& z_nd  = z_phys_nd->const_array(mfi);

                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    const Real x = prob_lo[0] + (i + myhalf_d) * dx[0];
                    const Real y = prob_lo[1] + (j + myhalf_d) * dx[1];

                    const Real z_sfc = fourth_d * (z_nd(i,j  ,0) + z_nd(i+1,j  ,0)
                                              + z_nd(i,j+1,0) + z_nd(i+1,j+1,0));
                    const Real z = std::max((z_cc(i,j,k)-z_sfc), amrex::Real(0));

                    const Real radius = std::sqrt((x - xf) * (x - xf) +
                                                  (y - yf) * (y - yf));

                    Real factor = one_d;
                    if ((z <= hf) && (radius <= (myhalf_d * df))) {
                        if (tf == 2) {
                            Real ratio = (hf - treeZm) / (hf - z);
                            if (z < treeZm) {
                                factor = amrex::Math::powi<6>(ratio) *
                                         std::exp(Real(6.0) * (one_d - ratio));
                            } else if (z <= hf) {
                                factor = std::sqrt(ratio) *
                                         std::exp(myhalf_d * (one_d - ratio));
                            }
                        }
                        levelDrag(i, j, k)   = cdf * af * factor;
                        if (store_frontal_area) {
                            frontalArea(i, j, k) = af * factor;
                        }
                    }
                });
            } // mfi
        } // ii (forest patch)

    } // end else (discrete patch mode)

    // Fillboundary for periodic ghost cell copy
    m_forest_drag->FillBoundary(geom.periodicity());
    if (m_frontal_area != nullptr) {
        m_frontal_area->FillBoundary(geom.periodicity());
    }

} // define_drag_field

#ifdef ERF_USE_NETCDF
/*
  Helper: read a 2D variable from a NetCDF file and return grid metadata.
  Data is broadcast to all MPI ranks.
*/
void
ForestDrag::read_netcdf_file (const std::string& filename,
                              const std::string& varname,
                              Vector<Real>& data,
                              int& nx, int& ny,
                              Real& dx, Real& dy,
                              Real& xmin, Real& ymin)
{
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(filename, NC_NOWRITE);

        if (!ncf.has_var(varname)) {
            Abort("Variable '" + varname + "' not found in file: " + filename);
        }

        auto var  = ncf.var(varname);
        auto dims = var.shape();

        const auto dim_names = var.dimnames();
        const auto is_y_dimension = [] (const std::string& name) {
            return name == "y" || name == "south_north";
        };
        const auto is_x_dimension = [] (const std::string& name) {
            return name == "x" || name == "west_east";
        };
        const auto is_time_dimension = [] (const std::string& name) {
            return name == "time" || name == "Time";
        };

        if (dims.size() == 2) {
            if (dim_names.size() != 2 || !is_y_dimension(dim_names[0]) ||
                !is_x_dimension(dim_names[1])) {
                Abort("Forest field '" + varname + "' in '" + filename +
                      "' must use dimensions (y,x) or (south_north,west_east)");
            }
            ny = dims[0];
            nx = dims[1];
        } else if (dims.size() == 3) {
            if (dim_names.size() != 3 || !is_time_dimension(dim_names[0]) ||
                !is_y_dimension(dim_names[1]) || !is_x_dimension(dim_names[2])) {
                Abort("Forest field '" + varname + "' in '" + filename +
                      "' must use dimensions (time,y,x) with time leading");
            }
            if (dims[0] < 1) {
                Abort("Forest field '" + varname + "' in '" + filename +
                      "' has no time records");
            }
            ny = dims[1];
            nx = dims[2];
        } else {
            Abort("Forest field '" + varname + "' in '" + filename +
                  "' must be 2D or 3D with a leading time dimension");
        }

        if (nx < 2 || ny < 2) {
            Abort("Forest field '" + varname + "' in '" + filename +
                  "' requires at least two x and y coordinates for bilinear interpolation");
        }

        const bool has_x = ncf.has_var("x");
        const bool has_y = ncf.has_var("y");
        const bool has_lon = ncf.has_var("lon");
        const bool has_lat = ncf.has_var("lat");
        const std::string coordinate_source =
            "Forest field '" + varname + "' in '" + filename + "'";
        const std::string coordinate_kind_error =
            erf_forest_utils::validate_cartesian_coordinates(
                has_x, has_y, has_lon, has_lat, coordinate_source);
        if (!coordinate_kind_error.empty()) {
            Abort(coordinate_kind_error);
        }
        if (has_x && has_y) {
            const std::string x_name = "x";
            const std::string y_name = "y";
            auto x_var = ncf.var(x_name);
            auto y_var = ncf.var(y_name);

            const auto x_shape = x_var.shape();
            const auto y_shape = y_var.shape();
            if (x_shape.size() != 1 || x_shape[0] != static_cast<std::size_t>(nx)) {
                Abort("Forest field '" + varname + "' in '" + filename + "' has " +
                      x_name + " coordinate shape inconsistent with its x dimension");
            }
            if (y_shape.size() != 1 || y_shape[0] != static_cast<std::size_t>(ny)) {
                Abort("Forest field '" + varname + "' in '" + filename + "' has " +
                      y_name + " coordinate shape inconsistent with its y dimension");
            }

            Vector<Real> x_coords(nx);
            Vector<Real> y_coords(ny);
            x_var.get(x_coords.data());
            y_var.get(y_coords.data());

            const std::string field_description =
                "Forest field '" + varname + "' in '" + filename + "'";
            std::string coordinate_error = erf_grid_utils::validate_uniform_axis(
                x_coords, nx, x_name, field_description, xmin, dx);
            if (coordinate_error.empty()) {
                coordinate_error = erf_grid_utils::validate_uniform_axis(
                    y_coords, ny, y_name, field_description, ymin, dy);
            }
            if (!coordinate_error.empty()) {
                Abort(coordinate_error);
            }

        } else {
            Abort(coordinate_source + " must provide projected Cartesian x and y coordinate variables");
        }

        data.resize(nx * ny);

        if (dims.size() == 2) {
            var.get(data.data());
        } else {
            std::vector<size_t> start = {0, 0, 0};
            std::vector<size_t> count = {1, static_cast<size_t>(ny),
                                             static_cast<size_t>(nx)};
            var.get(data.data(), start, count);
        }

        ncf.close();

        Print() << "ForestDrag: Read '" << varname << "' from " << filename << "\n"
                << "  Dimensions: " << nx << " x " << ny << "\n"
                << "  Grid: xmin=" << xmin << ", ymin=" << ymin
                << ", dx=" << dx << ", dy=" << dy << "\n";
    }

    ParallelDescriptor::Bcast(&nx,   1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ny,   1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&dx,   1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&dy,   1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&xmin, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ymin, 1, ParallelDescriptor::IOProcessorNumber());

    if (!ParallelDescriptor::IOProcessor()) {
        data.resize(nx * ny);
    }
    ParallelDescriptor::Bcast(data.data(), data.size(),
                              ParallelDescriptor::IOProcessorNumber());
}
#endif
