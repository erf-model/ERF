#include <ERF_ForestDrag.H>
#ifdef ERF_USE_NETCDF
#include <ERF_NCInterface.H>
#endif
#include <ERF_Constants.H>
#include <ERF_GridUtils.H>

using namespace amrex;

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
    while (file >> value1 >> value2 >> value3 >> value4 >> value5 >> value6 >>
           value7 >> value8) {
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

    m_gridded_cd.assign(m_grid_nx * m_grid_ny, m_cd_const);

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
                               MultiFab* z_phys_nd)
{
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
    m_frontal_area = std::make_unique<MultiFab>(ba,dm,1,1);
    m_frontal_area->setVal(zero);

    // Copy namespace-scope constants into automatic-storage values before
    // entering GPU lambdas.  NVCC does not make these host-side constexpr
    // variables available as device symbols automatically.
    const Real zero_d   = zero;
    const Real one_d    = one;
    const Real myhalf_d = myhalf;
    const Real fourth_d = fourth;

    if (m_use_gridded_data) {
        // =====================================================================
        // Gridded NetCDF mode: interpolate from gridded LAI/height/cd data
        // =====================================================================

        const Real* lai_data_h    = m_gridded_lai.data();
        const Real* height_data_h = m_gridded_height.data();
        const Real* cd_data_h     = m_gridded_cd.data();

        const int grid_nx   = m_grid_nx;
        const int grid_ny   = m_grid_ny;
        const Real grid_dx  = m_grid_dx;
        const Real grid_dy  = m_grid_dy;
        const Real grid_xmin = m_grid_xmin;
        const Real grid_ymin = m_grid_ymin;

        const int grid_size = m_grid_nx * m_grid_ny;
        Gpu::DeviceVector<Real> lai_data_d(grid_size);
        Gpu::DeviceVector<Real> height_data_d(grid_size);
        Gpu::DeviceVector<Real> cd_data_d(grid_size);

        Gpu::copy(Gpu::hostToDevice, lai_data_h, lai_data_h + grid_size,
                  lai_data_d.begin());
        Gpu::copy(Gpu::hostToDevice, height_data_h, height_data_h + grid_size,
                  height_data_d.begin());
        Gpu::copy(Gpu::hostToDevice, cd_data_h, cd_data_h + grid_size,
                  cd_data_d.begin());

        const Real* lai_data    = lai_data_d.data();
        const Real* height_data = height_data_d.data();
        const Real* cd_data     = cd_data_d.data();

        const int tree_type = m_tree_type;
        const Real laimax   = m_laimax;

        for (MFIter mfi(*m_forest_drag); mfi.isValid(); ++mfi) {
            Box gtbx = mfi.growntilebox();
            const Array4<Real>& levelDrag   = m_forest_drag->array(mfi);
            const Array4<Real>& frontalArea = m_frontal_area->array(mfi);
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

                    // Bilinear interpolation of Cd
                    Real cd00 = cd_data[jj_c * grid_nx + ii_c];
                    Real cd10 = cd_data[jj_c * grid_nx + (ii_c + 1)];
                    Real cd01 = cd_data[(jj_c + 1) * grid_nx + ii_c];
                    Real cd11 = cd_data[(jj_c + 1) * grid_nx + (ii_c + 1)];
                    Real cd_interp = (cd00 * (one_d - wx) + cd10 * wx) * (one_d - wy) +
                                     (cd01 * (one_d - wx) + cd11 * wx) * wy;

                    // Compute drag if within canopy
                    if (z < height_interp && height_interp > zero_d && lai_interp > zero_d) {
                        Real af;
                        Real factor = one_d;

                        if (tree_type == 1) {
                            af = lai_interp / height_interp;
                        } else {
                            Real treeZm = laimax * height_interp;

                            const int nk = 100;
                            const Real dz = height_interp / Real(nk);
                            Real expFun = zero_d;
                            Real ztree  = zero_d;

                            for (int kk = 0; kk < nk; ++kk) {
                                Real ratio = (height_interp - treeZm) / (height_interp - ztree);
                                if (ztree < treeZm) {
                                    expFun += amrex::Math::powi<6>(ratio) *
                                              std::exp(Real(6.0) * (one_d - ratio));
                                } else {
                                    expFun += std::sqrt(ratio) *
                                              std::exp(myhalf_d * (one_d - ratio));
                                }
                                ztree += dz;
                            }
                            af = lai_interp / (expFun * dz);

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
                        frontalArea(i, j, k) = af * factor;
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
                int nk      = 100;
                Real ztree  = 0;
                Real expFun = zero_d;
                const Real dz = hf / Real(nk);
                treeZm = laimaxf * hf;
                for (int k(0); k<nk; ++k) {
                    Real ratio = (hf - treeZm) / (hf - ztree);
                    if (ztree < treeZm) {
                        expFun += amrex::Math::powi<6>(ratio) *
                    std::exp(Real(6.0) * (one_d - ratio));
                    } else {
                        expFun += std::sqrt(ratio) *
                                  std::exp(myhalf_d * (one_d - ratio));
                    }
                    ztree += dz;
                }
                af = laif / (expFun * dz);
            }

            for (MFIter mfi(*m_forest_drag); mfi.isValid(); ++mfi) {
                Box gtbx = mfi.growntilebox();
                const Array4<Real>& levelDrag   = m_forest_drag->array(mfi);
                const Array4<Real>& frontalArea = m_frontal_area->array(mfi);
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
                        frontalArea(i, j, k) = af * factor;
                    }
                });
            } // mfi
        } // ii (forest patch)

    } // end else (discrete patch mode)

    // Fillboundary for periodic ghost cell copy
    m_forest_drag->FillBoundary(geom.periodicity());
    m_frontal_area->FillBoundary(geom.periodicity());

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
        if (has_x != has_y) {
            Abort("Forest field '" + varname + "' in '" + filename +
                  "' must provide both x and y coordinate variables");
        }
        if (!has_x && has_lon != has_lat) {
            Abort("Forest field '" + varname + "' in '" + filename +
                  "' must provide both lon and lat coordinate variables");
        }

        if ((has_x && has_y) || (has_lon && has_lat)) {
            const std::string x_name = has_x ? "x" : "lon";
            const std::string y_name = has_y ? "y" : "lat";
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

            if (!has_x) {
                Print() << "Warning: Using lon/lat coordinates as domain coordinates.\n";
            }
        } else {
            Abort("Forest field '" + varname + "' in '" + filename +
                  "' must provide x/y (or lon/lat) coordinate variables; "
                  "unit-spacing fallback is not supported");
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
