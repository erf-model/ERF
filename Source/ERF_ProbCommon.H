#ifndef ERF_PROBCOMMON_H_
#define ERF_PROBCOMMON_H_

#include <algorithm>
#include <cmath>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include "ERF_DataStruct.H"
#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_HSEUtils.H"
#include "ERF_GridUtils.H"
#include "ERF_TileNoZ.H"

#ifdef ERF_USE_NETCDF
#include <ERF_NCInterface.H>
#endif

struct ProbParmDefaults {
    amrex::Real rho_0 = one;
    amrex::Real T_0 = amrex::Real(300.0);
};

template <typename T>
T erf_get_single_value (std::istream& is, int n)
{
  T r = 0;
  std::string line;
  if (std::getline(is, line)) {
    std::stringstream ss(line);
    ss >> r;
    if (ss.fail()) {
      amrex::Abort("Failed to read");
    }
    if (ss >> r) {
      amrex::Print() << "Trying to read line " << n << " in this file " << std::endl;
      amrex::Abort("Wrong format: more than one number in this line");
    }
  } else {
    amrex::Abort("Unable to read this line");
  }
  return r;
}

void
init_my_custom_terrain (const amrex::Geometry& geom,
                        amrex::FArrayBox& terrain_fab,
                        const double& time);

/**
 * Class to hold problem-specific routines
*/
class ProblemBase
{
public:

    /**
     * Virtual destructor to avoid data leakage with derived class
     */
    virtual ~ProblemBase () = default;

    /**
     * Function to initialize the hydrostatic reference density
     *
     * @param[out] rho_hse hydrostatic reference density
     * @param[in] z_phys_nd height coordinate at nodes
     * @param[in] z_phys_cc height coordinate at cell centers
     * @param[in] geom container for geometric information
    */
    virtual void
    erf_init_dens_hse_dry (amrex::MultiFab& /*rho_hse*/,
                           std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                           std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/,
                           amrex::Geometry const& /*geom*/,
                           const amrex::Vector<amrex::Real>& /*stretched_dz_h*/,
                           bool /*is_constant_dz*/, bool /*is_stretched_dz*/)
    {
        amrex::Print() << "Hydrostatically balanced density was NOT set"
            << " -- an appropriate init_type should probably have been specified"
            << " (e.g., input_sounding, WRFInput, or Metgrid)"
            << std::endl;
        amrex::Error("Should never call this version of erf_init_dens_hse_dry for "+name()+" problem");
    }

    virtual void
    erf_init_const_dens_hse(amrex::MultiFab& /*rho_hse*/)
    {
        amrex::Error("Should never call this version of erf_init_const_dens_hse for "+name()+" problem");
    }

    virtual void
    erf_init_const_dens_and_th_hse (amrex::MultiFab& /*rho_hse*/, amrex::MultiFab& /* p_hse*/,
                                   amrex::MultiFab& /* pi_hse*/, amrex::MultiFab& /*th_hse*/,
                                   amrex::MultiFab& /* qv_hse*/, amrex::Real /*l_rdOcp*/)
    {
        amrex::Error("Should never call this version of erf_init_const_dens_and_th_hse for "+name()+" problem");
    }

    virtual void
    erf_init_const_dens_and_linear_th_hse (amrex::MultiFab& /*rho_hse*/, amrex::MultiFab& /* p_hse*/,
                                           amrex::MultiFab& /* pi_hse*/, amrex::MultiFab& /*th_hse*/,
                                           amrex::MultiFab& /* qv_hse*/, amrex::Real /*l_rdOcp*/,
                                           std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/)
    {
        amrex::Error("Should never call this version of erf_init_const_dens_and_linear_th_hse for "+name()+" problem");
    }

    virtual void
    erf_init_dens_hse_moist (amrex::MultiFab& /*rho_hse*/,
                             std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                             amrex::Geometry const& /*geom*/)
    {

    }

    /**
     * Function to perform custom initialization of a test problem
     *
     * @param[in]  bx cell-centered box on which to initialize scalars
     * @param[out] state cell-centered variables to be filled in this routine
     * @param[out] r_hse       hydrostatic reference density
     * @param[out] p_hse       hydrostatic reference pressure
     * @param[in] z_nd         height coordinate at nodes
     * @param[in] z_cc         height coordinate at cell centers
     * @param[in] mf_m         map factor on cell centers
     * @param[in] mf_u         map factor on x-faces
     * @param[in] mf_v         map factor on y-faces
     * @param[in] sc           SolverChoice structure that carries parameters
     * @param[in] lev          AMR level
    */
    virtual void
    init_custom_pert (const amrex::Box&  /*bx*/,
                      amrex::Array4<amrex::Real const> const& /*state*/,
                      amrex::Array4<amrex::Real      > const& /*state_pert*/,
                      amrex::Array4<amrex::Real      > const& /*r_hse*/,
                      amrex::Array4<amrex::Real      > const& /*p_hse*/,
                      amrex::Array4<amrex::Real const> const& /*z_nd*/,
                      amrex::Array4<amrex::Real const> const& /*z_cc*/,
                      amrex::GeometryData const& /*geomdata*/,
                      amrex::Array4<amrex::Real const> const& /*mf_m*/,
                      const SolverChoice& /*sc*/, const int /*lev*/)
    {
        amrex::Print() << "No perturbation to background fields supplied for "
            << name() << " problem" << std::endl;
    }

    /**
     * Function to perform custom initialization of the velocities
     *
     * @param[in] xbx box on which to initialize x_vel_pert
     * @param[in] ybx box on which to initialize y_vel_pert
     * @param[in] zbx box on which to initialize z_vel_pert
     * @param[out] x_vel_pert  x-component of velocity perturbation to be filled in this routine
     * @param[out] y_vel_pert  y-component of velocity perturbation to be filled in this routine
     * @param[out] z_vel_pert  z-component of velocity perturbation to be filled in this routine
     * @param[in] z_nd         height coordinate at nodes
     * @param[in] mf_m         map factor on cell centers
     * @param[in] mf_u         map factor on x-faces
     * @param[in] mf_v         map factor on y-faces
     * @param[in] sc           SolverChoice structure that carries parameters
     * @param[in] lev          AMR level
    */
    virtual void
    init_custom_pert_vels (const amrex::Box& /*xbx*/,
                           const amrex::Box& /*ybx*/,
                           const amrex::Box& /*zbx*/,
                           amrex::Array4<amrex::Real      > const& /*x_vel_pert*/,
                           amrex::Array4<amrex::Real      > const& /*y_vel_pert*/,
                           amrex::Array4<amrex::Real      > const& /*z_vel_pert*/,
                           amrex::Array4<amrex::Real const> const& /*z_nd*/,
                           amrex::GeometryData const& /*geomdata*/,
                           amrex::Array4<amrex::Real const> const& /*mf_u*/,
                           amrex::Array4<amrex::Real const> const& /*mf_v*/,
                           const SolverChoice& /*sc*/, const int /*lev*/)
    {
        amrex::Print() << "No perturbation velocities supplied for " << name() << " problem" << std::endl;
    }

    /**
     * Function to update user-specified temperature source terms that can
     * vary with time and height.
     *
     * @param[in]  time             current time
     * @param[out] rhotheta_source  forcing profile
     * @param[in]  geom             container for geometric information
     * @param[in]  z_phys_cc        height coordinate at cell centers
    */
    virtual void
    update_rhotheta_sources (const double& /*time*/,
                             amrex::MultiFab* src,
                             const amrex::Geometry& /*geom*/,
                             std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/)
    {
        if (src->empty()) return;

        amrex::Warning("Temperature forcing not defined for "+name()+" problem");
        for ( amrex::MFIter mfi(*src, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const auto &box = mfi.tilebox();
            const amrex::Array4<amrex::Real>& src_arr = src->array(mfi);
            // src is a spatial function if erf.spatial_rhotheta_forcing = true
            // otherwise, qsrc_arr is defined only in over Z (box x and y dimensions are 1)
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // const amrex::Real z_cc = prob_lo[2] + (k+myhalf)* dx[2];
                // set RHS term of RhoTheta equation based on time, z_cc here
                src_arr(i, j, k) = zero;
            });
        }
    }

    /**
     * Function to update user-specified moisture source terms that can
     * vary with time and height.
     *
     * @param[in]  time             current time
     * @param[out] rhoqt_source     moisture forcing profile
     * @param[in]  geom             container for geometric information
     * @param[in]  z_phys_cc        height coordinate at cell centers
    */
    virtual void
    update_rhoqt_sources (const double& /*time*/,
                             amrex::MultiFab* qsrc,
                             const amrex::Geometry& /*geom*/,
                             std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/)
    {
        if (qsrc->empty()) return;

        amrex::Warning("Moisture forcing not defined for "+name()+" problem");
        for ( amrex::MFIter mfi(*qsrc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const auto &box = mfi.tilebox();
            const amrex::Array4<amrex::Real>& qsrc_arr = qsrc->array(mfi);
            // src is a spatial function if erf.spatial_moisture_forcing = true
            // otherwise, qsrc_arr is defined only in over Z (box x and y dimensions are 1)
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // const amrex::Real z_cc = prob_lo[2] + (k+myhalf)* dx[2];
                // set RHS term of RhoQ1 equation based on time, z_cc here
                qsrc_arr(i, j, k) = zero;
            });
        }
    }

    /**
     * Function to update the vertical velocity profile, used to add subsidence
     * source terms for x-mom, y-mom, rho*theta, rho*Q1, and rho*Q2.
     *
     * TODO: Currently, this is only called by InitData, so there is no time
     * dependence.
     *
     * @param[in]  time             current time
     * @param[out] wbar             w vel forcing profile
     * @param[in]  geom             container for geometric information
     * @param[in]  z_phys_cc        height coordinate at cell centers
    */
    virtual void
    update_w_subsidence (const double& /*time*/,
                         amrex::Vector<amrex::Real>& wbar,
                         amrex::Gpu::DeviceVector<amrex::Real>& d_wbar,
                         const amrex::MultiFab& /*state*/,
                         const amrex::Geometry& geom,
                         std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/)
    {
        if (wbar.empty()) return;

        amrex::Warning("Moisture forcing not defined for "+name()+" problem");

        const int khi              = geom.Domain().bigEnd()[2];
        // const amrex::Real* prob_lo = geom.ProbLo();
        // const auto dx              = geom.CellSize();
        for (int k = 0; k <= khi; k++)
        {
            // const amrex::Real z_cc = prob_lo[2] + (k+myhalf)* dx[2];
            // set vertical velocity profile based on time, z_cc here
            wbar[k] = zero;
        }

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
    }

    /**
     * Function to update user-specified geostrophic wind profile.
     *
     * @param[in]  time             current time
     * @param[out] u_geos           geostrophic wind profile
     * @param[out] v_geos           geostrophic wind profile
     * @param[in]  geom             container for geometric information
     * @param[in]  z_phys_cc        height coordinate at cell centers
    */
    virtual void
    update_geostrophic_profile (const double& /*time*/,
                         amrex::Vector<amrex::Real>& u_geos,
                         amrex::Gpu::DeviceVector<amrex::Real>& d_u_geos,
                         amrex::Vector<amrex::Real>& v_geos,
                         amrex::Gpu::DeviceVector<amrex::Real>& d_v_geos,
                         const amrex::Geometry& geom,
                         std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/)
    {
        if (u_geos.empty()) return;

        amrex::Warning("Geostrophic wind profile not defined for "+name()+" problem");

        const int khi              = geom.Domain().bigEnd()[2];
        // const amrex::Real* prob_lo = geom.ProbLo();
        // const auto dx              = geom.CellSize();
        for (int k = 0; k <= khi; k++)
        {
            // const amrex::Real z_cc = prob_lo[2] + (k+myhalf)* dx[2];
            // set RHS term of RhoTheta equation based on time, z_cc here
            u_geos[k] = zero;
            v_geos[k] = zero;
        }

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, u_geos.begin(), u_geos.end(), d_u_geos.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, v_geos.begin(), v_geos.end(), d_v_geos.begin());
    }

    /**
     * Function to perform custom initialization of buildings
     *
     * @param[in] geom container for geometric information
     * @param[out] z_phys_nd height coordinate at nodes
     * @param[in] time current time
    */
    void
    init_buildings_surface (const amrex::Geometry& geom,
                          amrex::FArrayBox& buildings_fab,
                          const double& time)
    {
        // Check if a valid text file exists for the buildings
        std::string fname;
        amrex::ParmParse pp("erf");
        auto valid_fname = pp.queryAdd("buildings_file_name",fname);

        if (valid_fname) {
            read_custom_terrain(fname,false,geom,buildings_fab,time);
        } else {
            init_my_custom_terrain(geom, buildings_fab, time);
        }
    }


    /**
     * Function to perform custom initialization of terrain
     *
     * @param[in] geom container for geometric information
     * @param[out] z_phys_nd height coordinate at nodes
     * @param[in] time current time
    */
    void
    init_terrain_surface (const amrex::Geometry& geom,
                          amrex::FArrayBox& terrain_fab,
                          const double& time)
    {
        amrex::ParmParse pp("erf");

        // NetCDF terrain file takes priority
        std::string fname_nc;
        auto valid_fname_nc = pp.query("terrain_file_name_nc", fname_nc);

        std::string fname, fname_usgs;
        auto valid_fname      = pp.query("terrain_file_name",      fname);
        auto valid_fname_USGS = pp.query("terrain_file_name_USGS", fname_usgs);

        if (valid_fname_nc) {
#ifdef ERF_USE_NETCDF
            read_terrain_netcdf(fname_nc, geom, terrain_fab, time);
#else
            amrex::Abort("ERF must be compiled with NetCDF support to use terrain_file_name_nc");
#endif
        } else if (valid_fname) {
            read_custom_terrain(fname, false, geom, terrain_fab, time);

        } else if (valid_fname_USGS) {
            read_custom_terrain(fname_usgs, true, geom, terrain_fab, time);

        } else {
            init_my_custom_terrain (geom, terrain_fab, time);
        }
    }

#ifdef ERF_USE_NETCDF
    void
    read_terrain_netcdf (const std::string& fname,
                         const amrex::Geometry& geom,
                         amrex::FArrayBox& terrain_fab,
                         const amrex::Real& /*time*/)
    {
        amrex::Print() << "Reading terrain from NetCDF file: " << fname << std::endl;

        amrex::Vector<amrex::Real> m_xterrain, m_yterrain, m_zterrain;
        int nx = 0, ny = 0;
        // Broadcast this as a trivially broadcastable integer before the
        // interpolation kernel is launched.
        int mass_point_terrain = 0;

        if (amrex::ParallelDescriptor::IOProcessor()) {
            auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);

            // Try common terrain variable names
            std::string height_varname = "height";
            if (!ncf.has_var("height")) {
                if      (ncf.has_var("z"))       height_varname = "z";
                else if (ncf.has_var("terrain"))  height_varname = "terrain";
                else if (ncf.has_var("HGT_M"))    height_varname = "HGT_M";
                else if (ncf.has_var("AGL"))       height_varname = "AGL";
                else
                    amrex::Abort("No terrain height variable found in NetCDF file. "
                                 "Expected: 'height', 'z', 'terrain', 'HGT_M', or 'AGL'");
            }
            amrex::Print() << "  Using terrain variable: " << height_varname << std::endl;

            auto var  = ncf.var(height_varname);
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

            bool has_time_dimension = false;
            if (dims.size() == 2) {
                if (dim_names.size() != 2 || !is_y_dimension(dim_names[0]) ||
                    !is_x_dimension(dim_names[1])) {
                    amrex::Abort("NetCDF terrain field '" + height_varname + "' in '" +
                                 fname + "' must use dimensions (y,x) or "
                                 "(south_north,west_east)");
                }
                ny = dims[0];
                nx = dims[1];
            } else if (dims.size() == 3) {
                if (dim_names.size() != 3 || !is_time_dimension(dim_names[0]) ||
                    !is_y_dimension(dim_names[1]) || !is_x_dimension(dim_names[2])) {
                    amrex::Abort("NetCDF terrain field '" + height_varname + "' in '" +
                                 fname + "' must use dimensions (time,y,x) with time leading; "
                                 "accepted layouts are (y,x) and (time,y,x)");
                }
                if (dims[0] < 1) {
                    amrex::Abort("NetCDF terrain field '" + height_varname + "' in '" +
                                 fname + "' has no time records");
                }
                has_time_dimension = true;
                ny = dims[1];
                nx = dims[2];
            } else {
                amrex::Abort("NetCDF terrain field '" + height_varname + "' in '" + fname +
                             "' has unsupported rank " + std::to_string(dims.size()) +
                             "; accepted layouts are (y,x) and (time,y,x)");
            }

            if (nx < 2 || ny < 2) {
                amrex::Abort("NetCDF terrain field '" + height_varname + "' in '" +
                             fname + "' requires at least two x and y coordinates "
                             "for bilinear interpolation");
            }

            const bool has_x = ncf.has_var("x");
            const bool has_y = ncf.has_var("y");
            if (has_x != has_y) {
                amrex::Abort("NetCDF terrain file '" + fname +
                             "' must provide both one-dimensional 'x' and 'y' coordinates");
            }

            const std::string terrain_description =
                "NetCDF terrain field '" + height_varname + "' in '" + fname + "'";
            m_zterrain.resize(nx * ny);
            if (!has_x) {
                if (height_varname != "HGT_M" || !has_time_dimension ||
                    dim_names != std::vector<std::string>{"Time", "south_north", "west_east"}) {
                    amrex::Abort("NetCDF terrain field '" + height_varname + "' in '" + fname +
                                 "' requires one-dimensional 'x' and 'y' coordinates; a genuine "
                                 "WPS file must provide HGT_M(Time,south_north,west_east)");
                }

                const auto require_real_attribute = [&] (const std::string& name) {
                    if (!ncf.has_attr(name)) {
                        amrex::Abort("WPS terrain file '" + fname +
                                     "' is missing required global attribute '" + name + "'");
                    }
                    std::vector<double> values;
                    ncf.get_attr(name, values);
                    if (values.size() != 1 || !std::isfinite(values[0])) {
                        amrex::Abort("WPS terrain file '" + fname +
                                     "' must have one finite global attribute '" + name + "'");
                    }
                    return values[0];
                };
                const auto require_int_attribute = [&] (const std::string& name) {
                    if (!ncf.has_attr(name)) {
                        amrex::Abort("WPS terrain file '" + fname +
                                     "' is missing required global attribute '" + name + "'");
                    }
                    std::vector<int> values;
                    ncf.get_attr(name, values);
                    if (values.size() != 1) {
                        amrex::Abort("WPS terrain file '" + fname +
                                     "' must have one global attribute '" + name + "'");
                    }
                    return values[0];
                };

                const int map_proj = require_int_attribute("MAP_PROJ");
                // Projection 6 (latitude/longitude) is deliberately excluded.
                // WPS writes DX/DY in degrees for that projection, whereas
                // everything below consumes them as metres -- both the
                // geometry-consistency check and the m_xterrain/m_yterrain
                // axes built from ProbLo + (i+0.5)*DX.  Accepting it would
                // either abort with a misleading "grid does not match"
                // message, or, if prob_hi were specified in degrees to get
                // past that check, silently place metre-valued terrain
                // heights on a degree-spaced horizontal grid and corrupt the
                // terrain-fitted metrics with no warning.
                if (map_proj != 1 && map_proj != 2 && map_proj != 3) {
                    std::string reason =
                        (map_proj == 6)
                        ? "; latitude/longitude WPS grids are not supported because their "
                          "DX/DY are in degrees rather than metres"
                        : "; supported projections are 1 (Lambert), 2 (polar), "
                          "and 3 (Mercator)";
                    amrex::Abort("WPS terrain file '" + fname + "' uses unsupported MAP_PROJ=" +
                                 std::to_string(map_proj) + reason);
                }
                const double cen_lat = require_real_attribute("CEN_LAT");
                const double cen_lon = require_real_attribute("CEN_LON");
                const double stand_lon = require_real_attribute("STAND_LON");
                const double wps_dx = require_real_attribute("DX");
                const double wps_dy = require_real_attribute("DY");
                if (wps_dx <= 0.0 || wps_dy <= 0.0 || std::abs(cen_lat) > 90.0 ||
                    std::abs(cen_lon) > 360.0 || std::abs(stand_lon) > 360.0) {
                    amrex::Abort("WPS terrain file '" + fname +
                                 "' has invalid projection or grid metadata");
                }
                // TRUELAT1 is required by all three remaining projections.
                const double truelat1 = require_real_attribute("TRUELAT1");
                if (std::abs(truelat1) > 90.0) {
                    amrex::Abort("WPS terrain file '" + fname +
                                 "' has invalid TRUELAT1");
                }
                if (map_proj == 1) {
                    const double truelat2 = require_real_attribute("TRUELAT2");
                    if (std::abs(truelat2) > 90.0) {
                        amrex::Abort("WPS terrain file '" + fname +
                                     "' has invalid TRUELAT2");
                    }
                }

                const int wps_nx = require_int_attribute("WEST-EAST_GRID_DIMENSION");
                const int wps_ny = require_int_attribute("SOUTH-NORTH_GRID_DIMENSION");
                if (wps_nx != nx + 1 || wps_ny != ny + 1) {
                    amrex::Abort("WPS terrain file '" + fname +
                                 "' has grid dimensions inconsistent with HGT_M: expected " +
                                 std::to_string(nx + 1) + " x " + std::to_string(ny + 1));
                }

                const auto validate_wps_coordinate = [&] (const std::string& name) {
                    if (!ncf.has_var(name)) {
                        amrex::Abort("WPS terrain file '" + fname +
                                     "' is missing required coordinate variable '" + name + "'");
                    }
                    const auto coordinate = ncf.var(name);
                    if (coordinate.shape() != dims ||
                        coordinate.dimnames() != dim_names) {
                        amrex::Abort("WPS terrain coordinate '" + name + "' in '" + fname +
                                     "' must have dimensions (Time,south_north,west_east) "
                                     "matching HGT_M");
                    }
                    amrex::Vector<amrex::Real> values(nx * ny);
                    const std::vector<size_t> start{0, 0, 0};
                    const std::vector<size_t> count{
                        1, static_cast<size_t>(ny), static_cast<size_t>(nx)};
                    coordinate.get(values.data(), start, count);
                    if (!std::all_of(values.begin(), values.end(),
                                     [] (amrex::Real value) { return std::isfinite(value); })) {
                        amrex::Abort("WPS terrain coordinate '" + name + "' in '" + fname +
                                     "' contains nonfinite values");
                    }
                };
                validate_wps_coordinate("XLAT_M");
                validate_wps_coordinate("XLONG_M");

                const auto dx = geom.CellSizeArray();
                const auto close_enough = [] (amrex::Real lhs, double rhs) {
                    const amrex::Real scale = std::max(amrex::Real(1.0),
                                                        std::abs(static_cast<amrex::Real>(rhs)));
                    return std::abs(lhs - static_cast<amrex::Real>(rhs)) <=
                           amrex::Real(1.0e-6) * scale;
                };
                if (geom.Domain().length(0) != nx || geom.Domain().length(1) != ny ||
                    !close_enough(dx[0], wps_dx) || !close_enough(dx[1], wps_dy) ||
                    !close_enough(geom.ProbHi(0) - geom.ProbLo(0), wps_dx * nx) ||
                    !close_enough(geom.ProbHi(1) - geom.ProbLo(1), wps_dy * ny)) {
                    amrex::Abort("WPS terrain grid does not match the ERF Geometry; "
                                 "HGT_M dimensions, DX/DY, and the ERF domain extent must agree");
                }

                // WPS HGT_M/XLAT_M/XLONG_M are mass-point fields.  Their
                // logical grid is anchored at the ERF lower corner, with the
                // first mass point half a cell from that corner.  The
                // explicit clamp below is the documented one-sided boundary
                // conversion from mass points to ERF terrain nodes.
                m_xterrain.resize(nx);
                m_yterrain.resize(ny);
                for (int i = 0; i < nx; ++i) {
                    m_xterrain[i] = geom.ProbLo(0) +
                                    (static_cast<amrex::Real>(i) + amrex::Real(0.5)) *
                                    static_cast<amrex::Real>(wps_dx);
                }
                for (int j = 0; j < ny; ++j) {
                    m_yterrain[j] = geom.ProbLo(1) +
                                    (static_cast<amrex::Real>(j) + amrex::Real(0.5)) *
                                    static_cast<amrex::Real>(wps_dy);
                }
                const std::vector<size_t> start{0, 0, 0};
                const std::vector<size_t> count{
                    1, static_cast<size_t>(ny), static_cast<size_t>(nx)};
                var.get(m_zterrain.data(), start, count);
                mass_point_terrain = 1;
            } else {
                const auto x_shape = ncf.var("x").shape();
                const auto y_shape = ncf.var("y").shape();
                if (x_shape.size() != 1 || x_shape[0] != static_cast<std::size_t>(nx)) {
                    amrex::Abort("NetCDF terrain x coordinate in '" + fname +
                                 "' must be one-dimensional with length " + std::to_string(nx));
                }
                if (y_shape.size() != 1 || y_shape[0] != static_cast<std::size_t>(ny)) {
                    amrex::Abort("NetCDF terrain y coordinate in '" + fname +
                                 "' must be one-dimensional with length " + std::to_string(ny));
                }

                m_xterrain.resize(nx);
                m_yterrain.resize(ny);
                ncf.var("x").get(m_xterrain.data());
                ncf.var("y").get(m_yterrain.data());

                amrex::Real x_origin, x_spacing, y_origin, y_spacing;
                std::string coordinate_error = erf_grid_utils::validate_uniform_axis(
                    m_xterrain, nx, "x", terrain_description, x_origin, x_spacing);
                if (coordinate_error.empty()) {
                    coordinate_error = erf_grid_utils::validate_uniform_axis(
                        m_yterrain, ny, "y", terrain_description, y_origin, y_spacing);
                }
                if (!coordinate_error.empty()) {
                    amrex::Abort(coordinate_error);
                }

                if (has_time_dimension) {
                    const std::vector<size_t> start{0, 0, 0};
                    const std::vector<size_t> count{
                        1, static_cast<size_t>(ny), static_cast<size_t>(nx)};
                    var.get(m_zterrain.data(), start, count);
                } else {
                    var.get(m_zterrain.data());
                }
            }
            if (!std::all_of(m_zterrain.begin(), m_zterrain.end(),
                             [] (amrex::Real value) { return std::isfinite(value); })) {
                amrex::Abort(terrain_description + " contains nonfinite height values");
            }
            ncf.close();

            amrex::Print() << "  Grid size: " << nx << " x " << ny << "\n"
                           << "  X range: [" << m_xterrain[0] << ", " << m_xterrain[nx-1] << "]\n"
                           << "  Y range: [" << m_yterrain[0] << ", " << m_yterrain[ny-1] << "]\n"
                           << "  Height range: ["
                           << *std::min_element(m_zterrain.begin(), m_zterrain.end()) << ", "
                           << *std::max_element(m_zterrain.begin(), m_zterrain.end()) << "]\n";
        }

        // Broadcast to all MPI ranks
        amrex::ParallelDescriptor::Bcast(&nx, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&ny, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&mass_point_terrain, 1,
                                         amrex::ParallelDescriptor::IOProcessorNumber());

        if (!amrex::ParallelDescriptor::IOProcessor()) {
            m_xterrain.resize(nx);
            m_yterrain.resize(ny);
            m_zterrain.resize(nx * ny);
        }

        amrex::ParallelDescriptor::Bcast(m_xterrain.data(), nx,      amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(m_yterrain.data(), ny,      amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(m_zterrain.data(), nx * ny, amrex::ParallelDescriptor::IOProcessorNumber());

        // Copy to GPU
        amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(nx), d_yterrain(ny), d_zterrain(nx * ny);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_xterrain.begin(), m_xterrain.end(), d_xterrain.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_yterrain.begin(), m_yterrain.end(), d_yterrain.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_zterrain.begin(), m_zterrain.end(), d_zterrain.begin());

        amrex::Real* d_xt = d_xterrain.data();
        amrex::Real* d_yt = d_yterrain.data();
        amrex::Real* d_zt = d_zterrain.data();

        auto dx_arr  = geom.CellSizeArray();
        auto prob_lo = geom.ProbLoArray();

        int ilo = geom.Domain().smallEnd(0);
        int jlo = geom.Domain().smallEnd(1);
        int klo = geom.Domain().smallEnd(2);
        int ihi = geom.Domain().bigEnd(0) + 1;
        int jhi = geom.Domain().bigEnd(1) + 1;

        amrex::Box zbx = terrain_fab.box();
        amrex::Array4<amrex::Real> const& z_arr = terrain_fab.array();
        const amrex::Real terrain_dx = m_xterrain[1] - m_xterrain[0];
        const amrex::Real terrain_dy = m_yterrain[1] - m_yterrain[0];

        if (!mass_point_terrain && !zbx.isEmpty()) {
            amrex::Box target_box = zbx & amrex::convert(geom.Domain(), zbx.ixType());
            if (!target_box.isEmpty()) {
                target_box.setSmall(2, zbx.smallEnd(2));
                target_box.setBig(2, zbx.smallEnd(2));
                amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
                amrex::ReduceData<amrex::Long, amrex::Long> reduce_data(reduce_op);
                reduce_op.eval(target_box, reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) -> amrex::GpuTuple<amrex::Long, amrex::Long> {
                        const amrex::Real x = prob_lo[0] + static_cast<amrex::Real>(i) * dx_arr[0];
                        const amrex::Real y = prob_lo[1] + static_cast<amrex::Real>(j) * dx_arr[1];
                        const auto x_stencil = erf_grid_utils::uniform_interpolation_stencil(
                            x, d_xt[0], terrain_dx, nx);
                        const auto y_stencil = erf_grid_utils::uniform_interpolation_stencil(
                            y, d_yt[0], terrain_dy, ny);
                        const amrex::Long outside =
                            (x_stencil.inside && y_stencil.inside) ? amrex::Long(0) : amrex::Long(1);
                        return {amrex::Long(1), outside};
                    });
                // No parallel reduction here: init_terrain_surface is called
                // on every rank with its own private full-domain FArrayBox,
                // not a distributed MultiFab, so every rank has already
                // counted the whole domain.  Summing across ranks would
                // simply multiply both counts by the rank count.
                const auto local = reduce_data.value();
                const amrex::Long total_targets = amrex::get<0>(local);
                const amrex::Long outside_targets = amrex::get<1>(local);
                if (amrex::ParallelDescriptor::IOProcessor() && outside_targets > 0) {
                    const amrex::Real source_xhi = m_xterrain[nx-1];
                    const amrex::Real source_yhi = m_yterrain[ny-1];
                    const amrex::Real target_xlo = geom.ProbLo(0);
                    const amrex::Real target_ylo = geom.ProbLo(1);
                    const amrex::Real target_xhi = geom.ProbHi(0);
                    const amrex::Real target_yhi = geom.ProbHi(1);
                    const amrex::Real percentage = total_targets > 0
                        ? amrex::Real(100.0) * static_cast<amrex::Real>(outside_targets) /
                          static_cast<amrex::Real>(total_targets) : amrex::Real(0.0);
                    amrex::Print() << "WARNING: Terrain source grid does not cover ERF terrain nodes: outside "
                                   << outside_targets << " of " << total_targets << " ("
                                   << percentage << "%). Source extent x=[" << m_xterrain[0]
                                   << ", " << source_xhi << "] y=[" << m_yterrain[0] << ", "
                                   << source_yhi << "]; target extent x=[" << target_xlo << ", "
                                   << target_xhi << "] y=[" << target_ylo << ", " << target_yhi
                                   << "]. Outside nodes are set to zero.\n";
                }
            }
        }

        amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
        {
            int ii = amrex::min(amrex::max(i, ilo), ihi);
            int jj = amrex::min(amrex::max(j, jlo), jhi);

            amrex::Real x = prob_lo[0] + static_cast<amrex::Real>(ii) * dx_arr[0];
            amrex::Real y = prob_lo[1] + static_cast<amrex::Real>(jj) * dx_arr[1];

            const amrex::Real dx_t = d_xt[1] - d_xt[0];
            const amrex::Real dy_t = d_yt[1] - d_yt[0];
            if (mass_point_terrain) {
                // A WPS mass grid does not contain the ERF terrain nodes at
                // the four outer edges.  Hold the nearest mass value there;
                // this is intentional, finite one-sided staggering rather
                // than accidental extrapolation.
                x = amrex::min(amrex::max(x, d_xt[0]), d_xt[nx-1]);
                y = amrex::min(amrex::max(y, d_yt[0]), d_yt[ny-1]);
            }
            const auto x_stencil = erf_grid_utils::uniform_interpolation_stencil(
                x, d_xt[0], dx_t, nx);
            const auto y_stencil = erf_grid_utils::uniform_interpolation_stencil(
                y, d_yt[0], dy_t, ny);

            if (x_stencil.inside && y_stencil.inside) {
                const int i_t = x_stencil.lower;
                const int j_t = y_stencil.lower;
                const amrex::Real wx = x_stencil.weight;
                const amrex::Real wy = y_stencil.weight;

                // Data is stored [y, x] in NetCDF
                amrex::Real z00 = d_zt[ j_t      * nx +  i_t     ];
                amrex::Real z10 = d_zt[ j_t      * nx + (i_t + 1)];
                amrex::Real z01 = d_zt[(j_t + 1) * nx +  i_t     ];
                amrex::Real z11 = d_zt[(j_t + 1) * nx + (i_t + 1)];

                z_arr(i, j, klo) = (amrex::Real(1.0) - wx) * (amrex::Real(1.0) - wy) * z00 +
                                         wx * (amrex::Real(1.0) - wy) * z10 +
                                   (amrex::Real(1.0) - wx) * wy * z01 +
                                         wx * wy * z11;
            } else {
                z_arr(i, j, klo) = amrex::Real(0.0);
            }
        });

        // The interpolation kernel dereferences the device vectors above;
        // synchronize before their scope ends and their storage is released.
        amrex::Gpu::synchronize();
    }
#endif

    void
    read_custom_terrain(const std::string& fname,
                         const bool is_usgs,
                         const amrex::Geometry& geom,
                         amrex::FArrayBox& terrain_fab,
                         const double& /*time*/)
    {
        amrex::Vector<amrex::Real> m_xterrain,m_yterrain,m_zterrain;

        int nx = 0; int ny = 0;

        if (amrex::ParallelDescriptor::IOProcessor()) {
            // Read terrain file
            amrex::Print()<<"Reading terrain file: "<< fname<< std::endl;
            std::ifstream file(fname);

            if (!file.is_open()) {
                amrex::Abort("Error: Could not open the file " + fname+ "\n");
            }

            // Check if file is empty
            if (file.peek() == std::ifstream::traits_type::eof()) {
                amrex::Abort("Error: The file  " + fname+ " is empty.\n");
            }

            amrex::Real value1,value2,value3;

            if (is_usgs) {
                amrex::Real lat_min, lon_min;

                file >> lon_min >> lat_min;
                 if(std::fabs(lon_min) > amrex::Real(180.0)) {
                    amrex::Error("The value of longitude for entry in the first line in " + fname
                                 + " should not exceed amrex::Real(180.) It is " + std::to_string(lon_min));
                }
                 if(std::fabs(lat_min) > amrex::Real(90.0)) {
                    amrex::Error("The value of latitude for entry in the first line in " + fname
                                 + " should not exceed amrex::Real(90.) It is " + std::to_string(lat_min));
                }

                file >> nx >> ny;

                int counter = 0;
                while (file >> value1 >> value2 >> value3) {
                    m_xterrain.push_back(value1);
                    if(counter%nx==0) {
                        m_yterrain.push_back(value2);
                    }
                    m_zterrain.push_back(value3);
                    counter += 1;
                }
                AMREX_ASSERT(m_xterrain.size() == static_cast<long int>(nx*ny));
                AMREX_ASSERT(m_yterrain.size() == static_cast<long int>(ny));
                AMREX_ASSERT(m_zterrain.size() == static_cast<long int>(nx*ny));

            } else {
                int cnt = 1;
                nx = erf_get_single_value<int>(file,cnt); cnt++;
                ny = erf_get_single_value<int>(file,cnt); cnt++;
                amrex::Print()<<"Expecting " << nx << " values of x,    " <<
                                                ny << " values of y, and " <<
                                             nx*ny << " values of z" << std::endl;
                AMREX_ALWAYS_ASSERT(nx > 0);
                AMREX_ALWAYS_ASSERT(ny > 0);
                m_xterrain.resize(nx);
                m_yterrain.resize(ny);
                m_zterrain.resize(nx * ny);
                for (int n = 0; n < nx; n++) {
                    m_xterrain[n] = erf_get_single_value<amrex::Real>(file,cnt);
                    cnt++;
                }
                for (int n = 0; n < ny; n++) {
                    m_yterrain[n] = erf_get_single_value<amrex::Real>(file,cnt);
                    cnt++;
                }
                for (int n = 0; n < nx * ny; n++) {
                    m_zterrain[n] = erf_get_single_value<amrex::Real>(file,cnt);
                    cnt++;
                }
            }

            // Close the file after reading
            file.close();
        }

        amrex::ParallelDescriptor::Bcast(&nx,1,amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&ny,1,amrex::ParallelDescriptor::IOProcessorNumber());

        AMREX_ALWAYS_ASSERT(nx > 0 && ny > 0);
        if (is_usgs) {
            // The usgs branch below always uses a full 2x2 stencil
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nx > 1 && ny > 1,
                "USGS terrain file must have at least two points in each direction");
        }

        int nz = nx * ny;
        int nx_vals = is_usgs ? (nx * ny) : nx;

        m_xterrain.resize(nx_vals);
        m_yterrain.resize(ny);
        m_zterrain.resize(nz);

        amrex::ParallelDescriptor::Bcast(m_xterrain.data(),nx_vals,amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(m_yterrain.data(),ny,amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(m_zterrain.data(),nz,amrex::ParallelDescriptor::IOProcessorNumber());

        // Copy data to the GPU
        amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(nx_vals),d_yterrain(ny),d_zterrain(nz);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_xterrain.begin(), m_xterrain.end(), d_xterrain.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_yterrain.begin(), m_yterrain.end(), d_yterrain.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_zterrain.begin(), m_zterrain.end(), d_zterrain.begin());

        amrex::Real* d_xt = d_xterrain.data();
        amrex::Real* d_yt = d_yterrain.data();
        amrex::Real* d_zt = d_zterrain.data();

        auto dx = geom.CellSizeArray();
        auto ProbLoArr = geom.ProbLoArray();

        int ilo = geom.Domain().smallEnd(0);
        int jlo = geom.Domain().smallEnd(1);
        int klo = geom.Domain().smallEnd(2);
        int ihi = geom.Domain().bigEnd(0) + 1;
        int jhi = geom.Domain().bigEnd(1) + 1;

        amrex::Box zbx = terrain_fab.box();
        amrex::Array4<amrex::Real> const& z_arr = terrain_fab.array();

        amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
        {
            // Clip indices for ghost cells
            int ii = amrex::min(amrex::max(i,ilo),ihi);
            int jj = amrex::min(amrex::max(j,jlo),jhi);

            // Location of nodes
            amrex::Real x = ProbLoArr[0] + ii * dx[0] + amrex::Real(1e-3);
            amrex::Real y = ProbLoArr[1] + jj * dx[1] + amrex::Real(1e-3);

            int ind11, ind12, ind21, ind22;
            amrex::Real x1, x2, y1, y2;

            int iindex_terrain=-1;
            int jindex_terrain=-1;
            // *******************************************************************
            // NOTE: usgs-format is contiguous in x
            // *******************************************************************
            if (is_usgs) {

                // NOTE: the searches below must distinguish "the node is below the first
                //       table value" (bracket index -1) from "no table value exceeds the
                //       node" (search never triggered), so we can not use the initial
                //       value of -1 as the not-found flag
                bool jfound = false;
                for (int it = 0; it < ny && !jfound; it++) {
                    if (d_yt[it] > y) {
                        jindex_terrain = it-1; jfound = true;
                    }
                }
                if (!jfound) {
                    jindex_terrain = ny-1;
                }
                // Keep the 2x2 stencil inside the table
                jindex_terrain = amrex::min(amrex::max(jindex_terrain,0), ny-2);

                int gstart = (jindex_terrain  )*nx;
                int gend   = (jindex_terrain+1)*nx-1;
                bool ifound = false;
                for (int it = gstart; it <= gend && !ifound; it++) {
                    if (d_xt[it] > x) {
                        iindex_terrain = it-gstart-1; ifound = true;
                    }
                }
                if (!ifound) {
                    iindex_terrain = nx-1;
                }
                // Keep the 2x2 stencil inside the table
                iindex_terrain = amrex::min(amrex::max(iindex_terrain,0), nx-2);

                // Define the four values to interpolate between
                ind11 = jindex_terrain*nx + iindex_terrain; // (x1,y1)
                ind12 = ind11+nx;                           // (x1,y2)
                ind21 = ind11+1;                            // (x2,y1)
                ind22 = ind12+1;                            // (x2,y2)

                x1 = d_xt[ind11];
                x2 = d_xt[ind21];
                y1 = d_yt[jindex_terrain];
                y2 = d_yt[jindex_terrain+1];

                amrex::Real denom = (x2-x1)*(y2-y1);

                if (denom == amrex::Real(0.0)) {
                    // Degenerate cell (repeated x or y in the table) -- use the corner value
                    z_arr(i,j,klo) = d_zt[ind11];
                } else {
                    amrex::Real  w_11 = (x2-x)*(y2-y)/denom; // (x1,y1)
                    amrex::Real  w_12 = (x2-x)*(y-y1)/denom; // (x1,y2)
                    amrex::Real  w_21 = (x-x1)*(y2-y)/denom; // (x2,y1)
                    amrex::Real  w_22 = (x-x1)*(y-y1)/denom; // (x2,y2)

                    // Do bilinear interpolation
                    z_arr(i,j,klo) = w_11*d_zt[ind11] + w_12*d_zt[ind12] + w_21*d_zt[ind21] + w_22*d_zt[ind22];
                }

            } else {

                bool jfound = false;
                for (int it = 0; it < ny && !jfound; it++) {
                    if (d_yt[it] > y) {
                        jindex_terrain = it-1; jfound = true;
                    }
                }
                if (!jfound) {
                    // No table value exceeds y -- fall back to the last row
                    jindex_terrain = ny-1;
                }
                jindex_terrain = amrex::max(jindex_terrain,0);

                bool ifound = false;
                for (int it = 0; it < nx && !ifound; it++) {
                    if (d_xt[it] > x) {
                        iindex_terrain = it-1; ifound = true;
                    }
                }
                if (!ifound) {
                    // No table value exceeds x -- fall back to the last column
                    iindex_terrain = nx-1;
                }
                iindex_terrain = amrex::max(iindex_terrain,0);

                // NOTE: these must be clamped since the fallbacks above can leave
                //       iindex_terrain = nx-1 (jindex_terrain = ny-1)
                int ip1 = amrex::min(iindex_terrain+1,nx-1);
                int jp1 = amrex::min(jindex_terrain+1,ny-1);

                // Define the four values to interpolate between
                x1 = d_xt[iindex_terrain];
                x2 = d_xt[ip1];
                y1 = d_yt[jindex_terrain];
                y2 = d_yt[jp1];

#if 1
                // *******************************************************************
                // NOTE: this format is contiguous in y to match the AMR-Wind read
                // *******************************************************************
                ind11 = iindex_terrain * ny + jindex_terrain;  // (x1,y1)
                ind21 = ip1            * ny + jindex_terrain;  // (x2,y1)

                ind12 = iindex_terrain * ny + jp1;             // (x1,y2)
                ind22 = ip1            * ny + jp1;             // (x2,y2)
#else
                // *******************************************************************
                // NOTE: this format is contiguous in x as an alternative
                // *******************************************************************

                ind11 = jindex_terrain * nx + iindex_terrain;  // (x1,y1)
                ind12 = jp1            * nx + iindex_terrain;  // (x1,y2)

                ind21 = jindex_terrain * nx + ip1;             // (x2,y1)
                ind22 = jp1            * nx + ip1;             // (x2,y2)
#endif

                // Only interpolate in a direction if we have a non-degenerate interval there
                bool interp_x = (ip1 != iindex_terrain) && (x2 != x1);
                bool interp_y = (jp1 != jindex_terrain) && (y2 != y1);

                if (!interp_x && !interp_y)
                {
                    z_arr(i,j,klo) = d_zt[ind11];
                }
                else if (interp_x && !interp_y)
                {
                    amrex::Real  w_11 = (x2-x); // (x1,y1)
                    amrex::Real  w_21 = (x-x1); // (x2,y1)
                    amrex::Real denom = (x2-x1);
                    z_arr(i,j,klo) = (w_11*d_zt[ind11] + w_21*d_zt[ind21])/denom;
                }
                else if (!interp_x && interp_y)
                {
                    amrex::Real  w_11 = (y2-y); // (x1,y1)
                    amrex::Real  w_12 = (y-y1); // (x1,y2)
                    amrex::Real denom = (y2-y1);
                    z_arr(i,j,klo) = (w_11*d_zt[ind11] + w_12*d_zt[ind12])/denom;
                }
                else
                {
                    amrex::Real  w_11 = (x2-x)*(y2-y); // (x1,y1)
                    amrex::Real  w_21 = (x-x1)*(y2-y); // (x2,y1)
                    amrex::Real  w_12 = (x2-x)*(y-y1); // (x1,y2)
                    amrex::Real  w_22 = (x-x1)*(y-y1); // (x2,y2)
                    amrex::Real denom = (x2-x1)*(y2-y1);
                    z_arr(i,j,klo) = (w_11*d_zt[ind11] + w_12*d_zt[ind12] + w_21*d_zt[ind21] + w_22*d_zt[ind22]) / denom;
                }
            } // usgs?
        });
    }

    /**
     * Function to perform custom initialization of terrain
     *
     * Note: Terrain functionality can also be used to provide grid stretching.
     *
     * @param[in] geom container for geometric information
     * @param[out] terrain_mf height coordinate at nodes
     * @param[in] time current time
    */
    virtual void
    init_custom_terrain (const amrex::Geometry& geom,
                         amrex::FArrayBox& terrain_fab,
                         const double& time)
    {
        std::string custom_terrain_type = "None";
        amrex::ParmParse pp_prob("prob"); pp_prob.queryAdd("custom_terrain_type",custom_terrain_type);

        if (custom_terrain_type != "None")
        {
            amrex::Print() << "Calling custom terrain initialization" << std::endl;
            init_my_custom_terrain(geom,terrain_fab,time);

        } else {
            amrex::Print() << "Initializing flat terrain" << std::endl;
            terrain_fab.template setVal<amrex::RunOn::Device>(0);

            if (SolverChoice::mesh_type == MeshType::VariableDz) {
                SolverChoice::set_mesh_type(MeshType::StretchedDz);
                amrex::Print() << "Resetting mesh type to StretchedDz" << std::endl;
            }
        }
    }

#ifdef ERF_USE_TERRAIN_VELOCITY
    virtual amrex::Real compute_terrain_velocity(const double /*time*/)
    {
        amrex::Error("Should never call compute_terrain_velocity for "+name()+" problem");
    }
#endif

    /**
     * Function to define the quantities needed to impose Rayleigh damping
     *
     * @param[out] rayleigh_ptrs = {strength of Rayleigh damping, reference values for xvel/yvel/zvel/theta used to define Rayleigh damping}
     * @param[in] geom container for geometric information
     * @param[in] z_phys_cc height coordinate at cell centers
    */
    virtual void
    erf_init_rayleigh (amrex::Vector<amrex::Vector<amrex::Real> >& /*rayleigh_ptrs*/,
                       amrex::Geometry const& /*geom*/,
                       std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                       amrex::Real /*zdamp*/)
    {
        // Default which does no harm
    }

protected:
    // Struct to store problem parameters
    ProbParmDefaults base_parms;

    /**
     * Function to update default base parameters, currently only used for
     * init_type == InitType::Uniform
     */
    void init_base_parms (amrex::Real rho_0, amrex::Real T_0) {
        base_parms.rho_0 = rho_0;
        base_parms.T_0 = T_0;
    }

    std::string name() {
        std::string prob_name;
        amrex::ParmParse pp("erf");
        pp.get("prob_name",prob_name); return prob_name;
    }
};


/**
 * Function to init the physical bounds of the domain
 * and instantiate a Problem derived from ProblemBase
*/
extern std::unique_ptr<ProblemBase> amrex_probinit (const amrex_real* problo,
                                                    const amrex_real* probhi) AMREX_ATTRIBUTE_WEAK;
#endif
