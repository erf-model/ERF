#include <ERF_ForestDrag.H>
#ifdef ERF_USE_NETCDF
#include <ERF_NCInterface.H>
#endif
#include <ERF_Constants.H>

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
                        std::string vegtype_file,
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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nx_check == m_grid_nx && ny_check == m_grid_ny &&
        std::abs(dx_check - m_grid_dx) < 1e-6 &&
        std::abs(dy_check - m_grid_dy) < 1e-6,
        "Grid dimensions must match across LAI and height files");

    if (!vegtype_file.empty()) {
        m_gridded_vegtype.resize(m_grid_nx * m_grid_ny);
        Print() << "Warning: vegetation type reading not yet implemented\n";
    }

    m_gridded_cd.assign(m_grid_nx * m_grid_ny, m_cd_const);

    Print() << "ForestDrag: Gridded LAI/height + constant Cd=" << m_cd_const << "\n"
            << "  Grid size: " << m_grid_nx << " x " << m_grid_ny << "\n"
            << "  Grid spacing: dx=" << m_grid_dx << ", dy=" << m_grid_dy << "\n"
            << "  Tree type: " << m_tree_type << "  LAImax: " << m_laimax << "\n";

#else
    Abort("ERF must be compiled with NetCDF support to use gridded forest data");
#endif
}

/*
  NetCDF Constructor — gridded LAI/height/Cd all from files
*/
ForestDrag::ForestDrag (std::string lai_file,
                        std::string height_file,
                        std::string cd_file,
                        std::string vegtype_file,
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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nx_check == m_grid_nx && ny_check == m_grid_ny &&
        std::abs(dx_check - m_grid_dx) < 1e-6 &&
        std::abs(dy_check - m_grid_dy) < 1e-6,
        "Grid dimensions must match across all forest input files");

    read_netcdf_file(cd_file, "cd", m_gridded_cd,
                     nx_check, ny_check, dx_check, dy_check,
                     xmin_check, ymin_check);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nx_check == m_grid_nx && ny_check == m_grid_ny,
        "Grid dimensions must match across all forest input files");

    if (!vegtype_file.empty()) {
        m_gridded_vegtype.resize(m_grid_nx * m_grid_ny);
        Print() << "Warning: vegetation type reading not yet implemented\n";
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
        if (ba[i].smallEnd(2) != geom.ProbLo(2)) {
            all_boxes_touch_bottom = false;
        }
    }
    AMREX_ALWAYS_ASSERT(all_boxes_touch_bottom);

    // Allocate the forest drag MF and frontal area MF
    // NOTE: 1 ghost cell for averaging to faces
    m_forest_drag.reset();
    m_forest_drag = std::make_unique<MultiFab>(ba,dm,1,1);
    m_forest_drag->setVal(0.);

    m_frontal_area.reset();
    m_frontal_area = std::make_unique<MultiFab>(ba,dm,1,1);
    m_frontal_area->setVal(0.);

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
                const Real x = prob_lo[0] + (i + 0.5) * dx[0];
                const Real y = prob_lo[1] + (j + 0.5) * dx[1];

                const Real z_sfc = 0.25 * (z_nd(i, j, 0) + z_nd(i + 1, j, 0) +
                                           z_nd(i, j + 1, 0) + z_nd(i + 1, j + 1, 0));
                const Real z = amrex::max((z_cc(i, j, k) - z_sfc), 0.0);

                Real fx = (x - grid_xmin) / grid_dx;
                Real fy = (y - grid_ymin) / grid_dy;

                int ii = static_cast<int>(amrex::Math::floor(fx));
                int jj = static_cast<int>(amrex::Math::floor(fy));

                if (ii >= 0 && ii <= grid_nx - 1 && jj >= 0 && jj <= grid_ny - 1) {
                    int ii_c = amrex::min(ii, grid_nx - 2);
                    int jj_c = amrex::min(jj, grid_ny - 2);

                    Real wx = fx - ii_c;
                    Real wy = fy - jj_c;

                    // Bilinear interpolation of LAI
                    Real lai00 = lai_data[jj_c * grid_nx + ii_c];
                    Real lai10 = lai_data[jj_c * grid_nx + (ii_c + 1)];
                    Real lai01 = lai_data[(jj_c + 1) * grid_nx + ii_c];
                    Real lai11 = lai_data[(jj_c + 1) * grid_nx + (ii_c + 1)];
                    Real lai_interp = (lai00 * (1.0 - wx) + lai10 * wx) * (1.0 - wy) +
                                      (lai01 * (1.0 - wx) + lai11 * wx) * wy;

                    // Bilinear interpolation of height
                    Real h00 = height_data[jj_c * grid_nx + ii_c];
                    Real h10 = height_data[jj_c * grid_nx + (ii_c + 1)];
                    Real h01 = height_data[(jj_c + 1) * grid_nx + ii_c];
                    Real h11 = height_data[(jj_c + 1) * grid_nx + (ii_c + 1)];
                    Real height_interp = (h00 * (1.0 - wx) + h10 * wx) * (1.0 - wy) +
                                         (h01 * (1.0 - wx) + h11 * wx) * wy;

                    // Bilinear interpolation of Cd
                    Real cd00 = cd_data[jj_c * grid_nx + ii_c];
                    Real cd10 = cd_data[jj_c * grid_nx + (ii_c + 1)];
                    Real cd01 = cd_data[(jj_c + 1) * grid_nx + ii_c];
                    Real cd11 = cd_data[(jj_c + 1) * grid_nx + (ii_c + 1)];
                    Real cd_interp = (cd00 * (1.0 - wx) + cd10 * wx) * (1.0 - wy) +
                                     (cd01 * (1.0 - wx) + cd11 * wx) * wy;

                    // Compute drag if within canopy
                    if (z < height_interp && height_interp > 0.0 && lai_interp > 0.0) {
                        Real af;
                        Real factor = 1.0;

                        if (tree_type == 1) {
                            af = lai_interp / height_interp;
                        } else {
                            Real treeZm = laimax * height_interp;

                            const int nk = 100;
                            const Real dz = height_interp / Real(nk);
                            Real expFun = 0.0;
                            Real ztree  = 0.0;

                            for (int kk = 0; kk < nk; ++kk) {
                                Real ratio = (height_interp - treeZm) / (height_interp - ztree);
                                if (ztree < treeZm) {
                                    expFun += std::pow(ratio, 6.0) *
                                              std::exp(6.0 * (1.0 - ratio));
                                } else {
                                    expFun += std::pow(ratio, 0.5) *
                                              std::exp(0.5 * (1.0 - ratio));
                                }
                                ztree += dz;
                            }
                            af = lai_interp / (expFun * dz);

                            Real ratio = (height_interp - treeZm) / (height_interp - z);
                            if (z < treeZm) {
                                factor = std::pow(ratio, 6.0) *
                                         std::exp(6.0 * (1.0 - ratio));
                            } else {
                                factor = std::pow(ratio, 0.5) *
                                         std::exp(0.5 * (1.0 - ratio));
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
            Real treeZm = zero;
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
                Real expFun = 0;
                const Real dz = hf / Real(nk);
                treeZm = laimaxf * hf;
                for (int k(0); k<nk; ++k) {
                    Real ratio = (hf - treeZm) / (hf - ztree);
                    if (ztree < treeZm) {
                        expFun += amrex::Math::powi<6>(ratio) *
                                  std::exp(6 * (1 - ratio));
                    } else {
                        expFun += std::pow(ratio, myhalf) *
                                  std::exp(myhalf * (1 - ratio));
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
                    const Real x = prob_lo[0] + (i + myhalf) * dx[0];
                    const Real y = prob_lo[1] + (j + myhalf) * dx[1];

                    const Real z_sfc = fourth * (z_nd(i,j  ,0) + z_nd(i+1,j  ,0)
                                              + z_nd(i,j+1,0) + z_nd(i+1,j+1,0));
                    const Real z = std::max((z_cc(i,j,k)-z_sfc), amrex::Real(0));

                    const Real radius = std::sqrt((x - xf) * (x - xf) +
                                                  (y - yf) * (y - yf));

                    Real factor = 1;
                    if ((z <= hf) && (radius <= (myhalf * df))) {
                        if (tf == 2) {
                            Real ratio = (hf - treeZm) / (hf - z);
                            if (z < treeZm) {
                                factor = amrex::Math::powi<6>(ratio) *
                                         std::exp(Real(6.0) * (one - ratio));
                            } else if (z <= hf) {
                                factor = std::pow(ratio, myhalf) *
                                         std::exp(myhalf * (one - ratio));
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

    // =========================================================================
    // Compute cumulative LAI from canopy top downward (for Beer's law radiation)
    //   cumLAI(k) = integral_{z_k}^{h_c} a(z') dz'
    // =========================================================================
    m_cumLAI.reset();
    m_cumLAI = std::make_unique<MultiFab>(ba, dm, 1, 1);
    m_cumLAI->setVal(0.);

    Gpu::synchronize();

    for (MFIter mfi(*m_cumLAI, false); mfi.isValid(); ++mfi) {
        const Box& vbx     = mfi.validbox();
        const int  klo     = vbx.smallEnd(2);
        const int  khi_box = vbx.bigEnd(2);

        Box bx2d(IntVect(vbx.smallEnd(0), vbx.smallEnd(1), 0),
                 IntVect(vbx.bigEnd(0),   vbx.bigEnd(1),   0));

        const auto& fa  = m_frontal_area->const_array(mfi);
        const auto& cum = m_cumLAI->array(mfi);
        const Real dz_cell = dx[2];

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            Real accum = 0.0;
            for (int k = khi_box; k >= klo; --k) {
                cum(i, j, k) = accum;
                accum += fa(i, j, k) * dz_cell;
            }
        });
    }
    m_cumLAI->FillBoundary(geom.periodicity());

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

        if (dims.size() == 2) {
            ny = dims[0];
            nx = dims[1];
        } else if (dims.size() == 3) {
            ny = dims[1];
            nx = dims[2];
        } else {
            Abort("Expected 2D or 3D (with time) data in NetCDF file");
        }

        if (ncf.has_var("x") && ncf.has_var("y")) {
            auto x_var = ncf.var("x");
            auto y_var = ncf.var("y");

            Vector<Real> x_coords(nx);
            Vector<Real> y_coords(ny);
            x_var.get(x_coords.data());
            y_var.get(y_coords.data());

            xmin = x_coords[0];
            ymin = y_coords[0];
            dx = (nx > 1) ? (x_coords[1] - x_coords[0]) : 1.0;
            dy = (ny > 1) ? (y_coords[1] - y_coords[0]) : 1.0;

        } else if (ncf.has_var("lon") && ncf.has_var("lat")) {
            auto lon_var = ncf.var("lon");
            auto lat_var = ncf.var("lat");

            Vector<Real> lon_coords(nx);
            Vector<Real> lat_coords(ny);
            lon_var.get(lon_coords.data());
            lat_var.get(lat_coords.data());

            xmin = lon_coords[0];
            ymin = lat_coords[0];
            dx = (nx > 1) ? (lon_coords[1] - lon_coords[0]) : 1.0;
            dy = (ny > 1) ? (lat_coords[1] - lat_coords[0]) : 1.0;

            Print() << "Warning: Using lon/lat coordinates as domain coordinates.\n";
        } else {
            Print() << "Warning: No coordinate variables found. Assuming unit spacing at (0,0).\n";
            xmin = 0.0;
            ymin = 0.0;
            dx = 1.0;
            dy = 1.0;
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

Real
ForestDrag::interpolate_gridded_data (const Vector<Real>& data,
                                      Real x, Real y) const
{
    Real fx = (x - m_grid_xmin) / m_grid_dx;
    Real fy = (y - m_grid_ymin) / m_grid_dy;

    int i = static_cast<int>(std::floor(fx));
    int j = static_cast<int>(std::floor(fy));

    if (i < 0 || i >= m_grid_nx - 1 || j < 0 || j >= m_grid_ny - 1) {
        return 0.0;
    }

    Real wx = fx - i;
    Real wy = fy - j;

    Real val00 = data[j * m_grid_nx + i];
    Real val10 = data[j * m_grid_nx + (i + 1)];
    Real val01 = data[(j + 1) * m_grid_nx + i];
    Real val11 = data[(j + 1) * m_grid_nx + (i + 1)];

    return (val00 * (1.0 - wx) + val10 * wx) * (1.0 - wy) +
           (val01 * (1.0 - wx) + val11 * wx) * wy;
}
