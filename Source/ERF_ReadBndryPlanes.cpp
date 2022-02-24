#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include <AMReX_PlotFileUtil.H>
#include "ERF_ReadBndryPlanes.H"
#include "AMReX_MultiFabUtil.H"

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

void ReadBndryPlanes::define_level_data(
    const amrex::Orientation ori, const amrex::Box& bx, const size_t nc)
{
    m_data_n[ori]->push_back(amrex::FArrayBox(bx, nc));
    m_data_np1[ori]->push_back(amrex::FArrayBox(bx, nc));
    m_data_np2[ori]->push_back(amrex::FArrayBox(bx, nc));
    m_data_interp[ori]->push_back(amrex::FArrayBox(bx, nc));
}

void ReadBndryPlanes::read_data_native(amrex::FArrayBox& xlo_plane, amrex::FArrayBox& ylo_plane,
                                 amrex::FArrayBox& xhi_plane, amrex::FArrayBox& yhi_plane)
{
#if 0
    const amrex::OrientationIter oit,
    amrex::BndryRegister& bndry_n,
    amrex::BndryRegister& bndry_np1,
    const int lev,
    const amrex::MultiFab* fld,
    const amrex::Real time,
    const amrex::Vector<amrex::Real>& times)
    const size_t nc = fld->nComp();
    const int nstart = 0;  //  = m_components[fld->id()];

    const int idx = closest_index(times, time);
    const int idxp1 = idx + 1;

    m_tn   = times[idx];
    m_tnp1 = times[idx+1];
    m_tnp2 = times[idx+2];

    auto ori = oit();

    //AMREX_ALWAYS_ASSERT(((m_tn <= time) && (time <= m_tnp1)));
    //AMREX_ALWAYS_ASSERT(fld->num_comp() == bndry_n[ori].nComp());
    //AMREX_ASSERT(bndry_n[ori].boxArray() == bndry_np1[ori].boxArray());

    const int normal = ori.coordDir();
    const auto& bbx = (*m_data_n[ori])[lev].box();
    const amrex::IntVect v_offset = offset(ori.faceDir(), normal);

    amrex::MultiFab bndry(
        bndry_n[ori].boxArray(), bndry_n[ori].DistributionMap(),
        bndry_n[ori].nComp(), 0, amrex::MFInfo());

    for (amrex::MFIter mfi(bndry); mfi.isValid(); ++mfi) {

        const auto& vbx = mfi.validbox();
        const auto& bndry_n_arr = bndry_n[ori].array(mfi);
        const auto& bndry_arr = bndry.array(mfi);

        const auto& bx = bbx & vbx;
        if (bx.isEmpty()) {
            continue;
        }

        amrex::ParallelFor(
            bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                bndry_arr(i, j, k, n) =
                    0.5 *
                    (bndry_n_arr(i, j, k, n) +
                     bndry_n_arr(
                         i + v_offset[0], j + v_offset[1], k + v_offset[2], n));
            });
    }

    bndry.copyTo((*m_data_n[ori])[lev], 0, nstart, nc);

    for (amrex::MFIter mfi(bndry); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& bndry_np1_arr = bndry_np1[ori].array(mfi);
        const auto& bndry_arr = bndry.array(mfi);

        const auto& bx = bbx & vbx;
        if (bx.isEmpty()) {
            continue;
        }

        amrex::ParallelFor(
            bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                bndry_arr(i, j, k, n) =
                    0.5 *
                    (bndry_np1_arr(i, j, k, n) +
                     bndry_np1_arr(
                         i + v_offset[0], j + v_offset[1], k + v_offset[2], n));
            });
    }

    bndry.copyTo((*m_data_np1[ori])[lev], 0, nstart, nc);
#endif
}

void ReadBndryPlanes::interpolate(const amrex::Real time)
{
    m_tinterp = time;
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();

        const int nlevels = m_data_n[ori]->size();
        for (int lev = 0; lev < nlevels; ++lev) {

            const auto& datn = (*m_data_n[ori])[lev];
            const auto& datnp1 = (*m_data_np1[ori])[lev];
            auto& dati = (*m_data_interp[ori])[lev];
            dati.linInterp<amrex::RunOn::Device>(
                datn, 0, datnp1, 0, m_tn, m_tnp1, m_tinterp, datn.box(), 0,
                dati.nComp());
        }
    }
}

ReadBndryPlanes::ReadBndryPlanes(const amrex::Geometry& geom): m_geom(geom)
{
    amrex::ParmParse pp("erf");

    last_file_read = -1;

    // What folder will the time series of planes be read from
    pp.get("bndry_file", m_filename);

    // time.dat will be in the same folder as the time series of data
    m_time_file = m_filename + "/time.dat";
}

void ReadBndryPlanes::read_header()
{
    BL_PROFILE("ERF::ReadBndryPlanes::read_header");

    // each pointer (at at given time) has 6 components, one for each orientation
    // TODO: we really only need 4 not 6
    int size = 2*AMREX_SPACEDIM;
    m_data_n.resize(size);
    m_data_np1.resize(size);
    m_data_np2.resize(size);
    m_data_interp.resize(size);

    // *********************************************************
    // Read the time.data file and store the timesteps and times
    // *********************************************************
    int time_file_length = 0;

    if (amrex::ParallelDescriptor::IOProcessor()) {

        std::string line;
        amrex::Print() <<" TRYING TO READ TIME FILE " << m_time_file << std::endl;
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


    amrex::Print() << "First time in time.dat: " << m_in_times[0] << std::endl;
    amrex::Print() << "Last  time in time.dat: " << m_in_times.back() << std::endl;

    // *********************************************************
    // Allocate space for all of the boundary planes we may need
    // *********************************************************
    const int lev = 0;
    int ncomp = 1;
    const amrex::Box& domain = m_geom.Domain();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();

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
        define_level_data(ori, pbx, ncomp);
    }
}

void ReadBndryPlanes::read_input_files(amrex::Real time, amrex::Real dt)
{
    BL_PROFILE("ERF::ReadBndryPlanes::read_input_files");

    amrex::Print() << "TIME " << time << " " << std::endl;

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
        read_file(idx_init);

        amrex::Print() << "Reading file with tstep " << m_in_timesteps[idx_init] << std::endl;

        idx_init = 1;
        read_file(idx_init);

        amrex::Print() << "Reading file with tstep " << m_in_timesteps[idx_init] << std::endl;

        idx_init = 2;
        read_file(idx_init);

        amrex::Print() << "Reading file with tstep " << m_in_timesteps[idx_init] << std::endl;

        last_file_read = idx_init;
    }

    // Compute the index such that time falls between times[idx] and times[idx+1]
    const int idx = closest_index(m_in_times, time);

    // These are the times that we need in order to ensure we have the necessary data
    //   for this timestep
    amrex::Real tn   = m_in_times[idx];
    amrex::Real tnp1 = m_in_times[idx+1];
    amrex::Real tnp2 = m_in_times[idx+2];

    amrex::Print() << " IDX IS " << idx << std::endl;
    amrex::Print() << " TIMES GIVEN IDX " << tn << " " << tnp1 << " " << tnp2 << std::endl;
    amrex::Print() << " idx / last_file_read: " << idx << " " << last_file_read << std::endl;

    // Now we need to read another file
    if (idx >= last_file_read-1 && last_file_read != m_in_times.size()-1) {
        int new_read = last_file_read+1;
        amrex::Print() << "Reading file with idx / tstep " << new_read << " " << m_in_timesteps[new_read] << std::endl;
        read_file(new_read);
        last_file_read = new_read;
        amrex::Print() << "Now setting last_file_read to " << new_read << std::endl;
    } else {
        amrex::Print() << " idx < last_file_read: " << idx << " " << last_file_read << std::endl;
    }

#if 0
    // Go ahead and interpolate in time if we have already read the files we need
    if ((tn() <= time) && (time+dt <= tnp2()))
    {
        interpolate(time);
        return;
    } else

        amrex::Print() << "TIME " << time << std::endl;
        amrex::Print() << "MTIME " << m_in_timesteps[0] << " " << m_in_timesteps[1] << std::endl;
        const int index = closest_index(m_in_times, time);
        const int t_step1 = m_in_timesteps[index];
        //const int t_step2 = t_step1 + 1;
        const int t_step2 = m_in_timesteps[index+1];

        AMREX_ALWAYS_ASSERT(
            (m_in_times[index] <= time) && (time <= m_in_times[index + 1]));
#endif

    //interpolate(time);
}

void ReadBndryPlanes::read_file(const int idx)
{
#if 0
    const int t_step = m_in_timesteps[idx];
    const std::string chkname1 = m_filename + amrex::Concatenate("/bndry_output", t_step);

    const std::string level_prefix = "Level_";
    const int lev = 0;

    const amrex::Box& domain = m_geom.Domain();
    amrex::BoxArray ba(domain);
    amrex::DistributionMapping dm{ba};
    amrex::BndryRegister bndry(ba, dm, m_in_rad, m_out_rad, m_extent_rad, ncomp);
    bndry.setVal(1.0e13);

    std::string state_name = "temperature";
    std::string filename1 = amrex::MultiFabFileFullPrefix(lev, chkname1, level_prefix, state_name);

    int nstart = 0;
    int ncomp  = 1;

    // *********************************************************
    // Read in the BndryReg for all non-z faces
    // *********************************************************
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if (ori.coordDir() < 2) {

            std::string facename1 = amrex::Concatenate(filename1 + '_', ori, 1);
            bndry[ori].read(facename1);
        }
    }

    // *********************************************************
    // Copy from the BndryReg into a MultiFab then use copyTo
    //     to write from the MultiFab to a single FAB for each face
    // *********************************************************
    amrex::MultiFab bndryMF(
        bndry_n[ori].boxArray(), bndry_n[ori].DistributionMap(),
        bndry_n[ori].nComp(), 0, amrex::MFInfo());

    for (amrex::MFIter mfi(bndry); mfi.isValid(); ++mfi) {

        const auto& vbx = mfi.validbox();
        const auto& bndry_n_arr = bndry_n[ori].array(mfi);
        const auto& bndry_arr = bndryMF.array(mfi);

        const auto& bx = bbx & vbx;
        if (bx.isEmpty()) {
            continue;
        }

        amrex::ParallelFor(
            bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                bndry_arr(i, j, k, n) =
                    0.5 *
                    (bndry_n_arr(i, j, k, n) +
                     bndry_n_arr(
                         i + v_offset[0], j + v_offset[1], k + v_offset[2], n));
            });
    }
    bndry.copyTo((*m_data_n[ori])[lev], 0, nstart, nc);
#endif
}
