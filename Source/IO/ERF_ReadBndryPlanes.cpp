#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include <AMReX_PlotFileUtil.H>
#include "ERF_ReadBndryPlanes.H"
#include "IndexDefines.H"
#include "AMReX_MultiFabUtil.H"
#include "EOS.H"

//! Return closest index (from lower) of value in vector
AMREX_FORCE_INLINE int
closest_index(const amrex::Vector<amrex::Real>& vec, const amrex::Real value)
{
    auto const it = std::upper_bound(vec.begin(), vec.end(), value);
    AMREX_ALWAYS_ASSERT(it != vec.end());

    const int idx = std::distance(vec.begin(), it);
    return std::max(idx - 1, 0);
}

//! Return offset vector
AMREX_FORCE_INLINE amrex::IntVect offset(const int face_dir, const int normal)
{
    amrex::IntVect offset(amrex::IntVect::TheDimensionVector(normal));
    if (face_dir == 1) {
        for (auto& o : offset) {
            o *= -1;
        }
    }
    return offset;
}

void ReadBndryPlanes::define_level_data(int lev)
{
    amrex::Print() << "ReadBndryPlanes::define_level_data" << std::endl;
    // *********************************************************
    // Allocate space for all of the boundary planes we may need
    // *********************************************************
    int ncomp = BCVars::NumTypes;
    const amrex::Box& domain = m_geom.Domain();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if (ori.coordDir() < 2) {

            m_data_n[ori]      = std::make_unique<PlaneVector>();
            m_data_np1[ori]    = std::make_unique<PlaneVector>();
            m_data_np2[ori]    = std::make_unique<PlaneVector>();
            m_data_interp[ori] = std::make_unique<PlaneVector>();

            const auto& lo = domain.loVect();
            const auto& hi = domain.hiVect();

            amrex::IntVect plo(lo);
            amrex::IntVect phi(hi);
            const int normal = ori.coordDir();
            plo[normal] = ori.isHigh() ? hi[normal] + 1 : -1;
            phi[normal] = ori.isHigh() ? hi[normal] + 1 : -1;
            const amrex::Box pbx(plo, phi);
            m_data_n[ori]->push_back(amrex::FArrayBox(pbx, ncomp));
            m_data_np1[ori]->push_back(amrex::FArrayBox(pbx, ncomp));
            m_data_np2[ori]->push_back(amrex::FArrayBox(pbx, ncomp));
            m_data_interp[ori]->push_back(amrex::FArrayBox(pbx, ncomp));
        }
    }
}

amrex::Vector<std::unique_ptr<PlaneVector>>&
ReadBndryPlanes::interp_in_time(const amrex::Real& time)
{
    AMREX_ALWAYS_ASSERT(m_tn <= time && time <= m_tnp2);

    //amrex::Print() << "interp_in_time at time " << time << " given " << m_tn << " " << m_tnp1 << " " << m_tnp2 << std::endl;
    //amrex::Print() << "m_tinterp " << m_tinterp << std::endl;

    if (time == m_tinterp) {
        // We have already interpolated to this time
        return m_data_interp;

    } else {

        // We must now interpolate to a new time
        m_tinterp = time;

        if (time < m_tnp1) {
            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                auto ori = oit();
                if (ori.coordDir() < 2) {
                    const int nlevels = m_data_n[ori]->size();
                    for (int lev = 0; lev < nlevels; ++lev) {
                        const auto& datn   = (*m_data_n[ori])[lev];
                        const auto& datnp1 = (*m_data_np1[ori])[lev];
                        auto& dati = (*m_data_interp[ori])[lev];
                        dati.linInterp<amrex::RunOn::Device>(
                            datn, 0, datnp1, 0, m_tn, m_tnp1, m_tinterp, datn.box(), 0, dati.nComp());
                    }
                }
            }
        } else {
            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                auto ori = oit();
                if (ori.coordDir() < 2) {
                    const int nlevels = m_data_n[ori]->size();
                    for (int lev = 0; lev < nlevels; ++lev) {
                        const auto& datnp1 = (*m_data_np1[ori])[lev];
                        const auto& datnp2 = (*m_data_np2[ori])[lev];
                        auto& dati = (*m_data_interp[ori])[lev];
                        dati.linInterp<amrex::RunOn::Device>(
                            datnp1, 0, datnp2, 0, m_tnp1, m_tnp2, m_tinterp, datnp1.box(), 0,
                            dati.nComp());
                    }
                }
            }
        }
    }
    return m_data_interp;
}

ReadBndryPlanes::ReadBndryPlanes(const amrex::Geometry& geom): m_geom(geom)
{
    amrex::ParmParse pp("erf");

    last_file_read = -1;

    m_tinterp = -1.;

    // What folder will the time series of planes be read from
    pp.get("bndry_file", m_filename);

    is_velocity_read     = 0;
    is_density_read      = 0;
    is_temperature_read  = 0;
    is_theta_read        = 0;
    is_scalar_read       = 0;
    is_qv_read           = 0;
    is_qc_read           = 0;
    is_KE_read           = 0;
    is_QKE_read           = 0;

    if (pp.contains("bndry_input_var_names"))
    {
        int num_vars = pp.countval("bndry_input_var_names");
        m_var_names.resize(num_vars);
        pp.queryarr("bndry_input_var_names",m_var_names,0,num_vars);
        for (int i = 0; i < m_var_names.size(); i++) {
            if (m_var_names[i] == "velocity")     is_velocity_read = 1;
            if (m_var_names[i] == "density")      is_density_read = 1;
            if (m_var_names[i] == "temperature")  is_temperature_read = 1;
            if (m_var_names[i] == "theta")        is_theta_read = 1;
            if (m_var_names[i] == "scalar")       is_scalar_read = 1;
            if (m_var_names[i] == "qv")           is_qv_read = 1;
            if (m_var_names[i] == "qc")           is_qc_read = 1;
            if (m_var_names[i] == "KE")           is_KE_read = 1;
            if (m_var_names[i] == "QKE")          is_QKE_read = 1;
        }
    }

    // time.dat will be in the same folder as the time series of data
    m_time_file = m_filename + "/time.dat";

    // each pointer (at at given time) has 6 components, one for each orientation
    // TODO: we really only need 4 not 6
    int size = 2*AMREX_SPACEDIM;
    m_data_n.resize(size);
    m_data_np1.resize(size);
    m_data_np2.resize(size);
    m_data_interp.resize(size);
}

void ReadBndryPlanes::read_time_file()
{
    BL_PROFILE("ERF::ReadBndryPlanes::read_time_file");

    // *********************************************************
    // Read the time.data file and store the timesteps and times
    // *********************************************************
    int time_file_length = 0;

    if (amrex::ParallelDescriptor::IOProcessor()) {

        std::string line;
        std::ifstream time_file(m_time_file);
        if (!time_file.good()) {
            amrex::Abort("Cannot find time file: " + m_time_file);
        }
        while (std::getline(time_file, line)) {
            ++time_file_length;
        }

        time_file.close();
    }

    amrex::ParallelDescriptor::Bcast(
        &time_file_length, 1,
        amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());

    m_in_times.resize(time_file_length);
    m_in_timesteps.resize(time_file_length);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream time_file(m_time_file);
        for (int i = 0; i < time_file_length; ++i) {
            time_file >> m_in_timesteps[i] >> m_in_times[i];
        }
        // Sanity check that there are no duplicates or mis-orderings
        for (int i = 1; i < time_file_length; ++i) {
            if (m_in_timesteps[i] <= m_in_timesteps[i-1])
                amrex::Error("Bad timestep in time.dat file");
            if (m_in_times[i] <= m_in_times[i-1])
                amrex::Error("Bad time in time.dat file");
        }
        time_file.close();
    }

    amrex::ParallelDescriptor::Bcast(
        m_in_timesteps.data(), time_file_length,
        amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());

    amrex::ParallelDescriptor::Bcast(
        m_in_times.data(), time_file_length,
        amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());

    // Allocate data we will need -- for now just at one level
    int lev = 0;
    define_level_data(lev);
    amrex::Print() << "Successfully read time file and allocated data" << std::endl;
}

void ReadBndryPlanes::read_input_files(amrex::Real time, amrex::Real dt,
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals)
{
    BL_PROFILE("ERF::ReadBndryPlanes::read_input_files");

    // Assert that both the current time and the next time are within the bounds
    // of the data that we can read
    AMREX_ALWAYS_ASSERT((m_in_times[0] <= time) && (time <= m_in_times.back()));
    AMREX_ALWAYS_ASSERT((m_in_times[0] <= time+dt) && (time+dt <= m_in_times.back()));

    int ncomp = 1;

    const amrex::Box& domain = m_geom.Domain();
    amrex::BoxArray ba(domain);
    amrex::DistributionMapping dm{ba};
    amrex::BndryRegister bndryn(ba, dm, m_in_rad, m_out_rad, m_extent_rad, ncomp);
    bndryn.setVal(1.0e13);

    // The first time we enter this routine we read the first three files
    if (last_file_read == -1)
    {
        int idx_init = 0;
        read_file(idx_init,m_data_n,m_bc_extdir_vals);
        read_file(idx_init,m_data_interp,m_bc_extdir_vals); // We want to start with this filled
        m_tn = m_in_times[idx_init];

        idx_init = 1;
        read_file(idx_init,m_data_np1,m_bc_extdir_vals);
        m_tnp1 = m_in_times[idx_init];

        idx_init = 2;
        read_file(idx_init,m_data_np2,m_bc_extdir_vals);
        m_tnp2 = m_in_times[idx_init];

        last_file_read = idx_init;
    }

    // Compute the index such that time falls between times[idx] and times[idx+1]
    const int idx = closest_index(m_in_times, time);

    // Now we need to read another file
    if (idx >= last_file_read-1 && last_file_read != m_in_times.size()-1) {
        int new_read = last_file_read+1;

        // We need to change which data the pointers point to before we read in the new data
        // This doesn't actually move the data, just swaps the pointers
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            std::swap(m_data_n[ori]  ,m_data_np1[ori]);
            std::swap(m_data_np1[ori],m_data_np2[ori]);
        }

        // Set the times corresponding to the post-swap pointers
        m_tn   = m_tnp1;
        m_tnp1 = m_tnp2;
        m_tnp2 = m_in_times[new_read];

        read_file(new_read,m_data_np2,m_bc_extdir_vals);
        last_file_read = new_read;
    }

    AMREX_ASSERT(time    >= m_tn && time    <= m_tnp2);
    AMREX_ASSERT(time+dt >= m_tn && time+dt <= m_tnp2);
}

void ReadBndryPlanes::read_file(const int idx, amrex::Vector<std::unique_ptr<PlaneVector>>& data_to_fill,
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals)
{
    const int t_step = m_in_timesteps[idx];
    const std::string chkname1 = m_filename + amrex::Concatenate("/bndry_output", t_step);

    const std::string level_prefix = "Level_";
    const int lev = 0;

    const amrex::Box& domain = m_geom.Domain();
    amrex::BoxArray ba(domain);
    amrex::DistributionMapping dm{ba};

    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>,
                                                 AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

    for (int i = 0; i < BCVars::NumTypes; i++)
    {
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[i][ori];
        }
    }

    int n_for_density = -1;
    for (int i = 0; i < m_var_names.size(); i++)
    {
       if (m_var_names[i] == "density") n_for_density = i;
    }

    // We need to initialize all the components because we may not fill all of them from files,
    //    but the loop in the interpolate routine goes over all the components anyway
    int ncomp = BCVars::NumTypes;
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if (ori.coordDir() < 2) {
            amrex::FArrayBox& d = (*data_to_fill[ori])[lev];
            const auto& bx = d.box();
            amrex::Array4<amrex::Real> d_arr    = d.array();
            amrex::ParallelFor(
                bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                d_arr(i,j,k,n) = 0.;
            });
        }
    }

    for (int ivar = 0; ivar < m_var_names.size(); ivar++)
    {
        std::string var_name = m_var_names[ivar];

        std::string filename1 = amrex::MultiFabFileFullPrefix(lev, chkname1, level_prefix, var_name);

        int ncomp;

        if (var_name == "velocity") {
            ncomp = AMREX_SPACEDIM;
        } else {
            ncomp = 1;
        }

        int n_offset;
        if (var_name == "density")     n_offset = BCVars::Rho_bc_comp;
        if (var_name == "theta")       n_offset = BCVars::RhoTheta_bc_comp;
        if (var_name == "temperature") n_offset = BCVars::RhoTheta_bc_comp;
        if (var_name == "scalar")      n_offset = BCVars::RhoScalar_bc_comp;
        if (var_name == "qv")          n_offset = BCVars::RhoQv_bc_comp;
        if (var_name == "qc")          n_offset = BCVars::RhoQc_bc_comp;
        if (var_name == "velocity")    n_offset = BCVars::xvel_bc;

        // amrex::Print() << "Reading " << chkname1 << " for variable " << var_name << " with n_offset == " << n_offset << std::endl;

        amrex::BndryRegister bndry(ba, dm, m_in_rad, m_out_rad, m_extent_rad, ncomp);
        bndry.setVal(1.0e13);

        // *********************************************************
        // Read in the BndryReg for all non-z faces
        // *********************************************************
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
          auto ori = oit();
          if (ori.coordDir() < 2) {

            std::string facename1 = amrex::Concatenate(filename1 + '_', ori, 1);
            bndry[ori].read(facename1);

            const int normal = ori.coordDir();
            const amrex::IntVect v_offset = offset(ori.faceDir(), normal);

            const auto& bbx = (*data_to_fill[ori])[lev].box();

            // *********************************************************
            // Copy from the BndryReg into a MultiFab then use copyTo
            //     to write from the MultiFab to a single FAB for each face
            // *********************************************************
            amrex::MultiFab bndryMF(
                bndry[ori].boxArray(), bndry[ori].DistributionMap(),
                ncomp, 0, amrex::MFInfo());

            for (amrex::MFIter mfi(bndryMF); mfi.isValid(); ++mfi) {

                const auto& vbx = mfi.validbox();
                const auto& bndry_read_arr = bndry[ori].array(mfi);
                const auto& bndry_mf_arr   = bndryMF.array(mfi);

                const auto& bx = bbx & vbx;
                if (bx.isEmpty()) {
                    continue;
                }

                // We average the two cell-centered data points in the normal direction
                //    to define a Dirichlet value on the face itself.

                // This is the scalars -- they all get multiplied by rho, and in the case of
                //   reading in temperature, we must convert to theta first
                if (n_for_density >= 0) {
                  if (var_name == "temperature") {
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                             amrex::Real R1 =  bndry_read_arr(i, j, k, n_for_density);
                             amrex::Real R2 =  bndry_read_arr(i+v_offset[0],j+v_offset[1],k+v_offset[2],n_for_density);
                             amrex::Real T1 =  bndry_read_arr(i, j, k, 0);
                             amrex::Real T2 =  bndry_read_arr(i+v_offset[0],j+v_offset[1],k+v_offset[2],0);
                             amrex::Real Th1 = getThgivenRandT(R1,T1);
                             amrex::Real Th2 = getThgivenRandT(R2,T2);
                             bndry_mf_arr(i, j, k, 0) = 0.5 * (R1*Th1 + R2*Th2);
                        });
                  } else if (var_name == "scalar" || var_name == "qv" || var_name == "qc") {
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                             amrex::Real R1 =  bndry_read_arr(i, j, k, n_for_density);
                             amrex::Real R2 =  bndry_read_arr(i+v_offset[0],j+v_offset[1],k+v_offset[2],n_for_density);
                             bndry_mf_arr(i, j, k, 0) = 0.5 *
                                  ( R1 * bndry_read_arr(i, j, k, 0) +
                                    R2 * bndry_read_arr(i+v_offset[0],j+v_offset[1],k+v_offset[2], 0));
                        });
                   } else if (var_name == "density") {
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                             bndry_mf_arr(i, j, k, 0) = 0.5 *
                                  ( bndry_read_arr(i, j, k, 0) +
                                    bndry_read_arr(i+v_offset[0],j+v_offset[1],k+v_offset[2], 0));
                        });
                   }
                } else if (!ingested_density()) {
                  if (var_name == "temperature") {
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                             amrex::Real R1  = l_bc_extdir_vals_d[BCVars::Rho_bc_comp][ori];
                             amrex::Real R2  = l_bc_extdir_vals_d[BCVars::Rho_bc_comp][ori];
                             amrex::Real T1  = bndry_read_arr(i, j, k, 0);
                             amrex::Real T2  = bndry_read_arr(i+v_offset[0],j+v_offset[1],k+v_offset[2], 0);
                             amrex::Real Th1 = getThgivenRandT(R1,T1);
                             amrex::Real Th2 = getThgivenRandT(R2,T2);
                             bndry_mf_arr(i, j, k, 0) = 0.5 * (R1*Th1 + R2*Th2);
                        });
                  } else if (var_name == "scalar" || var_name == "qv" || var_name == "qc") {
                      amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                             amrex::Real R1  = l_bc_extdir_vals_d[BCVars::Rho_bc_comp][ori];
                             amrex::Real R2  = l_bc_extdir_vals_d[BCVars::Rho_bc_comp][ori];
                             bndry_mf_arr(i, j, k, 0) = 0.5 *
                                (R1 * bndry_read_arr(i, j, k, 0) +
                                 R2 * bndry_read_arr(i+v_offset[0],j+v_offset[1],k+v_offset[2], 0));
                        });
                  }
                }

                // This is velocity
                if (var_name == "velocity") {
                    amrex::ParallelFor(
                        bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                                bndry_mf_arr(i, j, k, n) = 0.5 *
                                  (bndry_read_arr(i, j, k, n) +
                                   bndry_read_arr(i+v_offset[0],j+v_offset[1],k+v_offset[2], n));
                        });
                }

            } // mfi
            bndryMF.copyTo((*data_to_fill[ori])[lev], 0, n_offset, ncomp);
          } // coordDir < 2
        } // ori
    } // var_name
}
