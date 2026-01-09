#ifndef ERF_PROBCOMMON_H_
#define ERF_PROBCOMMON_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>

#include "ERF_DataStruct.H"
#include "ERF_EOS.H"
#include "ERF_HSEUtils.H"
#include "ERF_TileNoZ.H"

struct ProbParmDefaults {
    amrex::Real rho_0 = 1.0;
    amrex::Real T_0 = 300.0;
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
    erf_init_dens_hse (amrex::MultiFab& /*rho_hse*/,
                       std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                       std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/,
                       amrex::Geometry const& /*geom*/)
    {
        amrex::Print() << "Hydrostatically balanced density was NOT set"
            << " -- an appropriate init_type should probably have been specified"
            << " (e.g., input_sounding, WRFInput, or Metgrid)"
            << std::endl;
        amrex::Error("Should never call erf_init_dens_hse for "+name()+" problem");
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
     * @param[in] xbx box on which to initialize x-component of velocity
     * @param[in] ybx box on which to initialize y-component of velocity
     * @param[in] zbx box on which to initialize z-component of velocity
     * @param[out] state cell-centered variables to be filled in this routine
     * @param[out] x_velocity  x-component of velocity to be filled in this routine
     * @param[out] y_velocity  y-component of velocity to be filled in this routine
     * @param[out] z_velocity  z-component of velocity to be filled in this routine
     * @param[out] r_hse       hydrostatic reference density
     * @param[out] p_hse       hydrostatic reference pressure
     * @param[in] z_nd         height coordinate at nodes
     * @param[in] z_cc         height coordinate at cell centers
     * @param[in] qv           water vapor
     * @param[in] qc           cloud water
     * @param[in] qi           cloud ice
     * @param[in] mf_m         map factor on cell centers
     * @param[in] mf_u         map factor on x-faces
     * @param[in] mf_v         map factor on y-faces
     * @param[in] sc           SolverChoice structure that carries parameters
    */
    virtual void
    init_custom_pert (const amrex::Box&  /*bx*/,
                      const amrex::Box& /*xbx*/,
                      const amrex::Box& /*ybx*/,
                      const amrex::Box& /*zbx*/,
                      amrex::Array4<amrex::Real const> const& /*state*/,
                      amrex::Array4<amrex::Real      > const& /*state_pert*/,
                      amrex::Array4<amrex::Real      > const& /*x_vel_pert*/,
                      amrex::Array4<amrex::Real      > const& /*y_vel_pert*/,
                      amrex::Array4<amrex::Real      > const& /*z_vel_pert*/,
                      amrex::Array4<amrex::Real      > const& /*r_hse*/,
                      amrex::Array4<amrex::Real      > const& /*p_hse*/,
                      amrex::Array4<amrex::Real const> const& /*z_nd*/,
                      amrex::Array4<amrex::Real const> const& /*z_cc*/,
                      amrex::GeometryData const& /*geomdata*/,
                      amrex::Array4<amrex::Real const> const& /*mf_m*/,
                      amrex::Array4<amrex::Real const> const& /*mf_u*/,
                      amrex::Array4<amrex::Real const> const& /*mf_v*/,
                      const SolverChoice& /*sc*/,
                      const int /*lev*/)
    {
        amrex::Print() << "No perturbation to background fields supplied for "
            << name() << " problem" << std::endl;
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
    update_rhotheta_sources (const amrex::Real& /*time*/,
                             amrex::Vector<amrex::Real>& src,
                             amrex::Gpu::DeviceVector<amrex::Real>& d_src,
                             const amrex::Geometry& geom,
                             std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/)
    {
        if (src.empty()) return;

        amrex::Warning("Temperature forcing not defined for "+name()+" problem");

        const int khi              = geom.Domain().bigEnd()[2];
        // const amrex::Real* prob_lo = geom.ProbLo();
        // const auto dx              = geom.CellSize();
        for (int k = 0; k <= khi; k++)
        {
            // const amrex::Real z_cc = prob_lo[2] + (k+0.5)* dx[2];
            // set RHS term of RhoTheta equation based on time, z_cc here
            src[k] = 0.0;
        }

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, src.begin(), src.end(), d_src.begin());
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
    update_rhoqt_sources (const amrex::Real& /*time*/,
                             amrex::Vector<amrex::Real>& qsrc,
                             amrex::Gpu::DeviceVector<amrex::Real>& d_qsrc,
                             const amrex::Geometry& geom,
                             std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/)
    {
        if (qsrc.empty()) return;

        amrex::Warning("Moisture forcing not defined for "+name()+" problem");

        const int khi              = geom.Domain().bigEnd()[2];
        // const amrex::Real* prob_lo = geom.ProbLo();
        // const auto dx              = geom.CellSize();
        for (int k = 0; k <= khi; k++)
        {
            // const amrex::Real z_cc = prob_lo[2] + (k+0.5)* dx[2];
            // set RHS term of RhoQ1 equation based on time, z_cc here
            qsrc[k] = 0.0;
        }

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, qsrc.begin(), qsrc.end(), d_qsrc.begin());
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
    update_w_subsidence (const amrex::Real& /*time*/,
                         amrex::Vector<amrex::Real>& wbar,
                         amrex::Gpu::DeviceVector<amrex::Real>& d_wbar,
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
            // const amrex::Real z_cc = prob_lo[2] + (k+0.5)* dx[2];
            // set vertical velocity profile based on time, z_cc here
            wbar[k] = 0.0;
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
    update_geostrophic_profile (const amrex::Real& /*time*/,
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
            // const amrex::Real z_cc = prob_lo[2] + (k+0.5)* dx[2];
            // set RHS term of RhoTheta equation based on time, z_cc here
            u_geos[k] = 0.0;
            v_geos[k] = 0.0;
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
                          const amrex::Real& time)
    {
        // Check if a valid text file exists for the buildings
        std::string fname;
        amrex::ParmParse pp("erf");
        auto valid_fname = pp.query("buildings_file_name",fname);

        if (valid_fname) {
            read_custom_terrain(fname,false,geom,buildings_fab,time);
        } else {
            init_custom_terrain(geom, buildings_fab, time);
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
                          const amrex::Real& time)
    {
        // Check if a valid text file exists for the terrain
        std::string fname, fname_usgs;
        amrex::ParmParse pp("erf");
        auto valid_fname      = pp.query("terrain_file_name",fname);
        auto valid_fname_USGS = pp.query("terrain_file_name_USGS",fname_usgs);

        bool is_usgs;

        if (valid_fname) {
            is_usgs = false;
            read_custom_terrain(fname,is_usgs,geom,terrain_fab,time);

        } else if (valid_fname_USGS) {
            is_usgs = true;
            read_custom_terrain(fname_usgs,is_usgs,geom,terrain_fab,time);

        } else {
            init_custom_terrain (geom, terrain_fab, time);
        }
    }

    void
    read_custom_terrain(const std::string& fname,
                         const bool is_usgs,
                         const amrex::Geometry& geom,
                         amrex::FArrayBox& terrain_fab,
                         const amrex::Real& /*time*/)
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
                 if(std::fabs(lon_min) > 180.0) {
                    amrex::Error("The value of longitude for entry in the first line in " + fname
                                 + " should not exceed 180. It is " + std::to_string(lon_min));
                }
                 if(std::fabs(lat_min) > 90.0) {
                    amrex::Error("The value of latitude for entry in the first line in " + fname
                                 + " should not exceed 90. It is " + std::to_string(lat_min));
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

        int nz = nx * ny;

        m_xterrain.resize(nx);
        m_yterrain.resize(ny);
        m_zterrain.resize(nz);

        amrex::ParallelDescriptor::Bcast(m_xterrain.data(),nx,amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(m_yterrain.data(),ny,amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(m_zterrain.data(),nz,amrex::ParallelDescriptor::IOProcessorNumber());

        // Copy data to the GPU
        amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(nx),d_yterrain(ny),d_zterrain(nz);
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
            amrex::Real x = ProbLoArr[0] + ii * dx[0] + 1e-3;
            amrex::Real y = ProbLoArr[1] + jj * dx[1] + 1e-3;

            int ind11, ind12, ind21, ind22;
            amrex::Real x1, x2, y1, y2;

            int iindex_terrain=-1;
            int jindex_terrain=-1;
            // *******************************************************************
            // NOTE: usgs-format is contiguous in x
            // *******************************************************************
            if (is_usgs) {

                for (int it = 0; it < ny && jindex_terrain == -1; it++) {
                    if (d_yt[it] > y) {
                        jindex_terrain = it-1;
                    }
                }
                if (jindex_terrain == -1) {
                    jindex_terrain = ny-1;
                }

                int gstart = (jindex_terrain  )*nx;
                int gend   = (jindex_terrain+1)*nx-1;
                for (int it = gstart; it <= gend && iindex_terrain == -1; it++) {
                    if (d_xt[it] > x) {
                        iindex_terrain = it-gstart;
                    }
                }

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
                amrex::Real  w_11 = (x2-x)*(y2-y)/denom; // (x1,y1)
                amrex::Real  w_12 = (x2-x)*(y-y1)/denom; // (x1,y2)
                amrex::Real  w_21 = (x-x1)*(y2-y)/denom; // (x2,y1)
                amrex::Real  w_22 = (x-x1)*(y-y1)/denom; // (x2,y2)

                // Do bilinear interpolation
                z_arr(i,j,klo) = w_11*d_zt[ind11] + w_12*d_zt[ind12] + w_21*d_zt[ind21] + w_22*d_zt[ind22];

            } else {

                for (int it = 0; it < ny && jindex_terrain == -1; it++) {
                    if (d_yt[it] > y) {
                        jindex_terrain = it-1;
                    }
                }
                if (jindex_terrain == -1) {
                    jindex_terrain = ny-1;
                }

                for (int it = 0; it < nx && iindex_terrain == -1; it++) {
                    if (d_xt[it] > x) {
                        iindex_terrain = it-1;
                    }
                }
                if (iindex_terrain == -1) {
                    iindex_terrain = nx-1;
                }

                // Define the four values to interpolate between
                x1 = d_xt[iindex_terrain];
                x2 = d_xt[iindex_terrain+1];
                y1 = d_yt[jindex_terrain];
                y2 = d_yt[jindex_terrain+1];

#if 1
                // *******************************************************************
                // NOTE: this format is contiguous in y to match the AMR-Wind read
                // *******************************************************************
                ind11 =           iindex_terrain         * ny + jindex_terrain;                   // (x1,y1)
                ind21 =  std::min(iindex_terrain+1,nx-1) * ny + jindex_terrain;                   // (x2,y1)

                ind12 =           iindex_terrain         * ny + std::min(jindex_terrain+1,ny-1);  // (x1,y2)
                ind22  = std::min(iindex_terrain+1,nx-1) * ny + std::min(jindex_terrain+1,ny-1);  // (x1,y2)
#else
                // *******************************************************************
                // NOTE: this format is contiguous in x as an alternative
                // *******************************************************************

                ind11 =           jindex_terrain         * nx + iindex_terrain;                   // (x1,y1)
                ind12 =  std::min(jindex_terrain+1,ny-1) * nx + iindex_terrain;                   // (x1,y2)

                ind21 =           jindex_terrain         * nx + std::min(iindex_terrain+1,nx-1);  // (x2,y1)
                ind22  = std::min(jindex_terrain+1,ny-1) * nx + std::min(iindex_terrain+1,nx-1);  // (x2,y2)
#endif

                if (iindex_terrain == nx-1 && jindex_terrain == ny-1)
                {
                    z_arr(i,j,klo) = d_zt[ind11];
                }
                else if (iindex_terrain != nx-1 && jindex_terrain == ny-1)
                {
                    amrex::Real  w_11 = (x2-x); // (x1,y1)
                    amrex::Real  w_21 = (x-x1); // (x2,y1)
                    amrex::Real denom = (x2-x1);
                    z_arr(i,j,klo) = (w_11*d_zt[ind11] + w_21*d_zt[ind21])/denom;
                }
                else if (iindex_terrain == nx-1 && jindex_terrain != ny-1)
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
    init_custom_terrain (const amrex::Geometry& /*geom*/,
                         amrex::FArrayBox& terrain_fab,
                         const amrex::Real& /*time*/)
    {
        amrex::Print() << "Initializing flat terrain" << std::endl;

        if (SolverChoice::mesh_type == MeshType::VariableDz) {
            SolverChoice::set_mesh_type(MeshType::StretchedDz);
            amrex::Print() << "Resetting mesh type to StretchedDz" << std::endl;
        }

        terrain_fab.template setVal<amrex::RunOn::Device>(0.0);
    }

#ifdef ERF_USE_TERRAIN_VELOCITY
    virtual amrex::Real compute_terrain_velocity(const amrex::Real /*time*/)
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
        amrex::Error("Should never call erf_init_rayleigh for "+name()+" problem");
    }

    /**
     * Function to set uniform background density and potential temperature fields
    */
    void
    init_uniform (const amrex::Box& bx, amrex::Array4<amrex::Real> const& state)
    {
        amrex::Real rho_0 = base_parms.rho_0;
        amrex::Real T_0 = base_parms.T_0;
        amrex::Print() << "Initializing uniform fields"
            << " rho=" << rho_0 << " theta=" << T_0
            << std::endl;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            state(i, j, k, Rho_comp) = rho_0;
            state(i, j, k, RhoTheta_comp) = rho_0 * T_0;
        });
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

    // descriptor for problem definition
    virtual std::string name () = 0;
};


/**
 * Function to init the physical bounds of the domain
 * and instantiate a Problem derived from ProblemBase
*/
extern std::unique_ptr<ProblemBase> amrex_probinit (const amrex_real* problo,
                                                    const amrex_real* probhi) AMREX_ATTRIBUTE_WEAK;
#endif
