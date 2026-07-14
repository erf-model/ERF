/*
 * Driver-owned I/O for the ERF Noah-MP interface: the per-step land plotfile
 * and checkpoint/restart of the full Noah-MP prognostic state (plus the precip
 * snapshots).
 */

#include <string>
#include <memory>

#include <AMReX_VisMF.H>
#include <AMReX_Utility.H>

#include <ERF_NOAHMP.H>

using namespace amrex;

// Per-step NetCDF land output; cadence owned by the caller. Fans out over the local
// boxes.
void
NOAHMP::Plot_Landfile (const int& nstep)
{
#ifdef AMREX_USE_MPI
    // Collective nf90_create/close runs on a land-owning sub-communicator so
    // land-free ranks never block land-owning ones.
    const int color = noahmpio_vect.empty() ? MPI_UNDEFINED : 1;
    MPI_Comm land_comm = MPI_COMM_NULL;
    MPI_Comm_split(ParallelDescriptor::Communicator(), color,
                   ParallelDescriptor::MyProc(), &land_comm);
    if (land_comm == MPI_COMM_NULL) return; // land-free rank
    const int land_comm_f = static_cast<int>(MPI_Comm_c2f(land_comm));
    for (NoahmpIO_type &noahmpio : noahmpio_vect) {
        const int saved = noahmpio.comm;
        noahmpio.comm   = land_comm_f;
        noahmpio.WriteLand(nstep);
        noahmpio.comm   = saved;
    }
    MPI_Comm_free(&land_comm);
#else
    for (NoahmpIO_type &noahmpio : noahmpio_vect) {
        noahmpio.WriteLand(nstep);
    }
#endif
}

// Write/read the full NoahMP prognostic state to/from a NetCDF restart dir so a
// restart reproduces a cold-start trajectory bitwise (#3255).
void
NOAHMP::Write_Lsm_Restart (const std::string& dir) const
{
#ifdef AMREX_USE_MPI
    const int color = noahmpio_vect.empty() ? MPI_UNDEFINED : 1;
    MPI_Comm land_comm = MPI_COMM_NULL;
    MPI_Comm_split(ParallelDescriptor::Communicator(), color,
                   ParallelDescriptor::MyProc(), &land_comm);
    if (land_comm == MPI_COMM_NULL) return; // land-free rank
    const int land_comm_f = static_cast<int>(MPI_Comm_c2f(land_comm));
    for (const NoahmpIO_type& noahmpio : noahmpio_vect) {
        NoahmpIO_type& io = const_cast<NoahmpIO_type&>(noahmpio);
        const int saved = io.comm;
        io.comm   = land_comm_f;
        io.WriteRestart(dir);
        io.comm   = saved;
    }
    MPI_Comm_free(&land_comm);
#else
    for (const NoahmpIO_type& noahmpio : noahmpio_vect) {
        const_cast<NoahmpIO_type&>(noahmpio).WriteRestart(dir);
    }
#endif

    // Also save the precip snapshots so the next land call after restart differences
    // against the same baseline (no precip dropped or double-counted).
    if (m_lev < int(m_precip_accum_prev.size())) {
        for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
            if (m_precip_accum_prev[m_lev][s]) {
                amrex::VisMF::Write(*m_precip_accum_prev[m_lev][s],
                                    m_precip_snapshot_name(dir, m_lev, s));
            }
        }
    }
}

void
NOAHMP::Read_Lsm_Restart (const std::string& dir)
{
#ifdef AMREX_USE_MPI
    const int color = noahmpio_vect.empty() ? MPI_UNDEFINED : 1;
    MPI_Comm land_comm = MPI_COMM_NULL;
    MPI_Comm_split(ParallelDescriptor::Communicator(), color,
                   ParallelDescriptor::MyProc(), &land_comm);
    if (land_comm == MPI_COMM_NULL) return; // land-free rank
    const int land_comm_f = static_cast<int>(MPI_Comm_c2f(land_comm));
    for (NoahmpIO_type& noahmpio : noahmpio_vect) {
        const int saved = noahmpio.comm;
        noahmpio.comm   = land_comm_f;
        noahmpio.ReadRestart(dir);
        noahmpio.comm   = saved;
    }
    MPI_Comm_free(&land_comm);
#else
    for (NoahmpIO_type& noahmpio : noahmpio_vect) {
        noahmpio.ReadRestart(dir);
    }
#endif

    // Read saved precip snapshots into a staging store; the first Advance copies them
    // into the live snapshots. Legacy checkpoints lack the files -> Advance cold-seeds.
    for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
        const std::string name = m_precip_snapshot_name(dir, m_lev, s);
        if (amrex::FileExists(name + "_H")) {
            auto mf = std::make_unique<amrex::MultiFab>();
            amrex::VisMF::Read(*mf, name);
            m_precip_accum_restored[s] = std::move(mf);
        }
    }
}
