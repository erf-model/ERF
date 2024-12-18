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
            << " (e.g., input_sounding, ideal, real, or metgrid)"
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
                      const SolverChoice& /*sc*/)
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
                             std::unique_ptr<amrex::MultiFab>& z_phys_cc)
    {
        if (src.empty()) return;

        if (z_phys_cc) {
            // use_terrain=1
            amrex::Error("Temperature forcing not defined for "+name()+" problem with terrain");
        } else {
            const int khi              = geom.Domain().bigEnd()[2];
            // const amrex::Real* prob_lo = geom.ProbLo();
            // const auto dx              = geom.CellSize();
            for (int k = 0; k <= khi; k++)
            {
                // const amrex::Real z_cc = prob_lo[2] + (k+0.5)* dx[2];
                // set RHS term of RhoTheta equation based on time, z_cc here
                src[k] = 0.0;
            }
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
                             std::unique_ptr<amrex::MultiFab>& z_phys_cc)
    {
        if (qsrc.empty()) return;

        if (z_phys_cc) {
            // use_terrain=1
            amrex::Error("Moisture forcing not defined for "+name()+" problem with terrain");
        } else {
            const int khi              = geom.Domain().bigEnd()[2];
            // const amrex::Real* prob_lo = geom.ProbLo();
            // const auto dx              = geom.CellSize();
            for (int k = 0; k <= khi; k++)
            {
                // const amrex::Real z_cc = prob_lo[2] + (k+0.5)* dx[2];
                // set RHS term of RhoQ1 equation based on time, z_cc here
                qsrc[k] = 0.0;
            }
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
                         std::unique_ptr<amrex::MultiFab>& z_phys_cc)
    {
        if (wbar.empty()) return;

        if (z_phys_cc) {
            // use_terrain=1
            amrex::Error("Moisture forcing not defined for "+name()+" problem with terrain");
        } else {
            const int khi              = geom.Domain().bigEnd()[2];
            // const amrex::Real* prob_lo = geom.ProbLo();
            // const auto dx              = geom.CellSize();
            for (int k = 0; k <= khi; k++)
            {
                // const amrex::Real z_cc = prob_lo[2] + (k+0.5)* dx[2];
                // set vertical velocity profile based on time, z_cc here
                wbar[k] = 0.0;
            }
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
                         std::unique_ptr<amrex::MultiFab>& z_phys_cc)
    {
        if (u_geos.empty()) return;

        if (z_phys_cc) {
            // use_terrain=1
            amrex::Error("Geostrophic wind profile not defined for "+name()+" problem with terrain");
        } else {
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
        }

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, u_geos.begin(), u_geos.end(), d_u_geos.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, v_geos.begin(), v_geos.end(), d_v_geos.begin());
    }

    /**
     * Function to perform custom initialization of terrain
     *
     * Note: Terrain functionality can also be used to provide grid stretching.
     *
     * @param[in] geom container for geometric information
     * @param[out] z_phys_nd height coordinate at nodes
     * @param[in] time current time
    */
    virtual void
    init_custom_terrain (const amrex::Geometry& geom,
                         amrex::MultiFab& z_phys_nd,
                         const amrex::Real& time)
    {
        // Check if a valid text file exists for the terrain
        std::string fname, fname_usgs, ftype;
        amrex::ParmParse pp("erf");
        auto valid_fname = pp.query("terrain_file_name",fname);
        auto valid_fname_USGS = pp.query("terrain_file_name_USGS",fname_usgs);
        if (valid_fname) {
            this->read_custom_terrain(fname,geom,z_phys_nd,time);
        } else if(valid_fname_USGS) {
            this->read_custom_terrain_USGS(fname_usgs,geom,z_phys_nd,time);
        } else {
            // Note that this only sets the terrain value at the ground IF k=0 is in the box
            amrex::Print() << "Initializing flat terrain at z=0" << std::endl;

            // Number of ghost cells
            int ngrow = z_phys_nd.nGrow();

            // Bottom of domain
            int k0 = 0;

            for ( amrex::MFIter mfi(z_phys_nd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                amrex::Box zbx = mfi.nodaltilebox(2);
                if (zbx.smallEnd(2) <= 0) {
                    // Grown box with no z range
                    amrex::Box xybx = mfi.growntilebox(ngrow);
                    xybx.setRange(2,0);

                    amrex::Array4<amrex::Real> const& z_arr = z_phys_nd.array(mfi);

                    ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
                        z_arr(i,j,k0) = 0.0;
                    });
                }
            }
            SolverChoice::set_flat_terrain_flag();
        }
    }


    /**
    * Function to perform custom initialization of terrain
    *
    * This version takes a single FArrayBox instead of a MultiFab
    *
    */
    virtual void
    init_custom_terrain (const amrex::Geometry& /*geom*/,
                         amrex::FArrayBox& z_phys_nd,
                         const amrex::Real& /*time*/)
    {
        // Note that this only sets the terrain value at the ground IF k=0 is in the box
        amrex::Print() << "Initializing flat terrain at z=0" << std::endl;

        // Bottom of domain
        int k0 = 0;

        // Grown box with no z range
        amrex::Box bx = z_phys_nd.box();
        amrex::Array4<amrex::Real> const& z_arr = z_phys_nd.array();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            z_arr(i,j,k0) = 0.0;
        });

        SolverChoice::set_flat_terrain_flag();
    }

    /**
     * Function to perform custom initialization of terrain from an input file
     *
     * Note: Terrain functionality can also be used to provide grid stretching.
     *
     * @param[in] geom container for geometric information
     * @param[out] z_phys_nd height coordinate at nodes
     * @param[in] time current time
    */
    virtual void
    read_custom_terrain (const std::string& fname,
                         const amrex::Geometry& geom,
                         amrex::MultiFab& z_phys_nd,
                         const amrex::Real& /*time*/)
    {
        // Read terrain file on the host
        amrex::Print()<<"Reading terrain file: "<< fname << std::endl;
        std::ifstream file(fname);
        amrex::Gpu::HostVector<amrex::Real> m_xterrain,m_yterrain,m_zterrain;
        amrex::Real value1,value2,value3;
        while(file>>value1>>value2>>value3){
            m_xterrain.push_back(value1);
            m_yterrain.push_back(value2);
            m_zterrain.push_back(value3);
        }
        file.close();

        // Copy data to the GPU
        int nnode = m_xterrain.size();
        amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(nnode),d_yterrain(nnode),d_zterrain(nnode);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_xterrain.begin(), m_xterrain.end(), d_xterrain.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_yterrain.begin(), m_yterrain.end(), d_yterrain.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_zterrain.begin(), m_zterrain.end(), d_zterrain.begin());
        amrex::Real* d_xt = d_xterrain.data();
        amrex::Real* d_yt = d_yterrain.data();
        amrex::Real* d_zt = d_zterrain.data();

        // Populate z_phys data
        amrex::Real tol = 1.0e-4;
        int ngrow = z_phys_nd.nGrow();
        auto dx = geom.CellSizeArray();
        auto ProbLoArr = geom.ProbLoArray();
        int ilo = geom.Domain().smallEnd(0);
        int jlo = geom.Domain().smallEnd(1);
        int klo = geom.Domain().smallEnd(2);
        int ihi = geom.Domain().bigEnd(0) + 1;
        int jhi = geom.Domain().bigEnd(1) + 1;
        for (amrex::MFIter mfi(z_phys_nd,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Grown box with no z range
            amrex::Box xybx = mfi.growntilebox(ngrow);
            xybx.setRange(2,0);

            amrex::Array4<amrex::Real> const& z_arr = z_phys_nd.array(mfi);
            amrex::ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,ilo),ihi);
                int jj = amrex::min(amrex::max(j,jlo),jhi);

                // Location of nodes
                amrex::Real x = ProbLoArr[0]  + ii  * dx[0];
                amrex::Real y = ProbLoArr[1]  + jj  * dx[1];
                int inode = ii + jj * (ihi-ilo+2); // stride is Nx+1
                if (std::sqrt(std::pow(x-d_xt[inode],2)+std::pow(y-d_yt[inode],2)) < tol) {
                    z_arr(i,j,klo) = d_zt[inode];
                } else {
                    // Unexpected list order, do brute force search
                    amrex::Real zloc = 0.0;
                    bool found = false;
                    for (int n=0; n<nnode; ++n)
                    {
                        amrex::Real delta=std::sqrt(std::pow(x-d_xt[n],2)+std::pow(y-d_yt[n],2));
                        if (delta < tol) {
                            found = true;
                            zloc = d_zt[n];
                            break;
                        }
                    }
                    AMREX_ASSERT_WITH_MESSAGE(found, "Location read from terrain file does not match the grid!");
                    amrex::ignore_unused(found);
                    z_arr(i,j,klo) = zloc;
                }
            });
        }
    }

    virtual void
    read_custom_terrain_USGS (const std::string& fname_usgs,
                               const amrex::Geometry& geom,
                              amrex::MultiFab& z_phys_nd,
                               const amrex::Real& /*time*/)
    {
        // Read terrain file on the host
        amrex::Print()<<"Reading terrain file: "<< fname_usgs << std::endl;
        std::ifstream file(fname_usgs);

        if (!file.is_open()) {
            amrex::Abort("Error: Could not open the file " + fname_usgs + "\n");
        }

        // Check if file is empty
        if (file.peek() == std::ifstream::traits_type::eof()) {
            amrex::Abort("Error: The file  " + fname_usgs + " is empty.\n");
        }

        amrex::Gpu::HostVector<amrex::Real> m_xterrain,m_yterrain,m_zterrain;
        amrex::Real value1,value2,value3;
        amrex::Real lat_min, lon_min;
        int nlons, nlats;

        file >> lon_min >> lat_min;
         if(std::fabs(lon_min) > 180.0) {
            amrex::Error("The value of longitude for entry in the first line in " + fname_usgs
                         + " should not exceed 180. It is " + std::to_string(lon_min));
        }
         if(std::fabs(lat_min) > 90.0) {
            amrex::Error("The value of latitude for entry in the first line in " + fname_usgs
                         + " should not exceed 90. It is " + std::to_string(lat_min));
        }

        file >> nlons >> nlats;

        int counter = 0;
        while(file >> value1 >> value2 >> value3){
            m_xterrain.push_back(value1);
            if(counter%nlons==0) {
                m_yterrain.push_back(value2);
            }
            m_zterrain.push_back(value3);
            counter += 1;
        }
        file.close();

        AMREX_ASSERT(m_xterrain.size() == static_cast<long unsigned int>(nlons*nlats));
        AMREX_ASSERT(m_yterrain.size() == static_cast<long unsigned int>(nlats));
        AMREX_ASSERT(m_zterrain.size() == static_cast<long unsigned int>(nlons*nlats));

        // Shift the terrain so that there is some flat region for the
        // inflow to come in

        amrex::ParmParse pp("erf");
        amrex::Real x_shift, y_shift;
        pp.query("windfarm_x_shift",x_shift);
        pp.query("windfarm_x_shift",y_shift);


        for(auto& v : m_xterrain){
            v += x_shift;
        }
        for(auto& v : m_yterrain){
            v += y_shift;
        }

        if (m_xterrain.empty()) {
            amrex::Abort("Error: m_xterrain is empty!\n");
        }

        auto index = std::max_element(m_xterrain.begin(), m_xterrain.end());
        amrex::Real xmax_terrain = *index;
        index = std::max_element(m_yterrain.begin(), m_yterrain.end());
        amrex::Real ymax_terrain = *index;

        index = std::min_element(m_zterrain.begin(), m_zterrain.end());
        amrex::Real zmin_terrain = *index;

        for(auto& v : m_zterrain){
            v -= zmin_terrain;
        }

        // Copy data to the GPU
        int nnode = nlons*nlats;
        amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(nnode),d_yterrain(nlats),d_zterrain(nnode);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_xterrain.begin(), m_xterrain.end(), d_xterrain.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_yterrain.begin(), m_yterrain.end(), d_yterrain.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_zterrain.begin(), m_zterrain.end(), d_zterrain.begin());
        amrex::Real* d_xt = d_xterrain.data();
        amrex::Real* d_yt = d_yterrain.data();
        amrex::Real* d_zt = d_zterrain.data();

        // Populate z_phys data
        int ngrow = z_phys_nd.nGrow();
        auto dx = geom.CellSizeArray();
        auto ProbLoArr = geom.ProbLoArray();
        int ilo = geom.Domain().smallEnd(0);
        int jlo = geom.Domain().smallEnd(1);
        int klo = geom.Domain().smallEnd(2);
        int ihi = geom.Domain().bigEnd(0) + 1;
        int jhi = geom.Domain().bigEnd(1) + 1;

        for (amrex::MFIter mfi(z_phys_nd,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Grown box with no z range
            amrex::Box xybx = mfi.growntilebox(ngrow);
            xybx.setRange(2,0);

            amrex::Array4<amrex::Real> const& z_arr = z_phys_nd.array(mfi);
            amrex::ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,ilo),ihi);
                int jj = amrex::min(amrex::max(j,jlo),jhi);

                // Location of nodes
                amrex::Real x = ProbLoArr[0]  + ii  * dx[0] + 1e-3;
                amrex::Real y = ProbLoArr[1]  + jj  * dx[1] + 1e-3;

                if(x < x_shift or x > xmax_terrain or
                   y < y_shift or y > ymax_terrain) {
                    z_arr(i,j,klo) = 0.0;
                } else {
                    int jindex_terrain = -1;
                    for (int it = 0; it < nlats && jindex_terrain == -1; it++) {
                        if (d_yt[it] > y) {
                            jindex_terrain = it-1;
                        }
                    }
                    //if(jindex_terrain == -1){
                    //    printf("Could not find the jindex for a y value of %0.15g\n",y);
                    //}
                    // Now the index to check for x goes from jindex_terrain*nlons to
                    // (jindex_terrain+1)*nlons-1
                    int iindex_terrain=-1;
                    int gstart = jindex_terrain*nlons;
                    int gend = (jindex_terrain+1)*nlons-1;
                    for (int it = gstart; it <= gend && iindex_terrain == -1; it++) {
                        if (d_xt[it] > x) {
                            iindex_terrain = it-gstart;
                        }
                    }
                    //if(iindex_terrain == -1){
                    //    printf("Could not find the iindex for a x value of %0.15g for a y index of %d and y value of %0.15g\n",x, jindex_terrain,y);
                    //}

                    // Now do an averaging
                    int ind11, ind12, ind21, ind22;
                    ind11 = jindex_terrain*nlons + iindex_terrain;
                    ind12 = ind11+nlons;
                    ind21 = ind11+1;
                    ind22 = ind12+1;
                    amrex::Real x1 = d_xt[ind11];
                    amrex::Real x2 = d_xt[ind21];
                    amrex::Real y1 = d_yt[jindex_terrain];
                    amrex::Real y2 = d_yt[jindex_terrain+1];
                    // Do bilinear interpolation
                    amrex::Real denom = (x2-x1)*(y2-y1);
                    amrex::Real w_11 = (x2-x)*(y2-y)/denom;
                    amrex::Real w_12 = (x2-x)*(y-y1)/denom;
                    amrex::Real w_21 = (x-x1)*(y2-y)/denom;
                    amrex::Real w_22 = (x-x1)*(y-y1)/denom;

                    z_arr(i,j,klo) = w_11*d_zt[ind11] + w_12*d_zt[ind12] + w_21*d_zt[ind21] + w_22*d_zt[ind22];

                    //z_arr(i,j,klo) = 0.25*(d_zt[ind1] + d_zt[ind2] + d_zt[ind3] + d_zt[ind4]);
                    // Smoothen out the edges using Gaussian
                    /*amrex::Real sigma = 600.0;
                    amrex::Real fac = 2.0*sigma*sigma;
                    z_arr(i,j,klo) = z_arr(i,j,klo)*(1.0-exp(-(x-x_shift)*(x-x_shift)/fac))*
                                                    (1.0-exp(-(x-xmax_terrain)*(x-xmax_terrain)/fac))*
                                                    (1.0-exp(-(y-y_shift)*(y-y_shift)/fac))*
                                                    (1.0-exp(-(y-ymax_terrain)*(y-ymax_terrain)/fac)); */
                    //std::cout << "Values of z is " << d_zt[ind1] << " " << d_zt[ind2] << " " << d_zt[ind3] << " " << d_zt[ind4] << "\n";
                }
            });
        }
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
     * Function to set uniform background density and pressure fields
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
