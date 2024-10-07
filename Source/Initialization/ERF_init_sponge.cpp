/**
 * \file ERF_init_sponge.cpp
 */
#include <ERF.H>
#include <ERF_ParFunctions.H>
#include <ERF_Interpolation_1D.H>

using namespace amrex;

/**
 * High level wrapper for sponge x and y velocities
 * level data from input sponge data.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::input_sponge (int lev)
{
    // We only want to read the file once
    if (lev == 0) {
        if (input_sponge_data.input_sponge_file.empty())
            Error("input_sounding file name must be provided via input");

        // this will interpolate the input profiles to the nominal height levels
        // (ranging from 0 to the domain top)
        input_sponge_data.read_from_file(geom[lev], zlevels_stag[lev]);
    }
}

/**
 * Initialization function for host and device vectors
 * used to store the effects of sponge Damping.
 */
void
ERF::initSponge ()
{
    h_sponge_ptrs.resize(max_level+1);
    d_sponge_ptrs.resize(max_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        // These have 2 components: ubar, vbar
        h_sponge_ptrs[lev].resize(Sponge::nvars_sponge);
        d_sponge_ptrs[lev].resize(Sponge::nvars_sponge);

        const int zlen_sponge = geom[lev].Domain().length(2);

        // Allocate space for these 1D vectors
        for (int n = 0; n < Sponge::nvars_sponge; n++) {
            h_sponge_ptrs[lev][n].resize(zlen_sponge, 0.0_rt);
            d_sponge_ptrs[lev][n].resize(zlen_sponge, 0.0_rt);
        }

    }
}

/**
 * Sets the sponge damping averaged quantities from an
 * externally supplied input sponge data file.
 *
 * @param[in] restarting Boolean parameter that indicates whether
                         we are currently restarting from a checkpoint file.
 */
void
ERF::setSpongeRefFromSounding (bool restarting)
{
    // If we are restarting then we haven't read the input_sponge file yet
    //    so we need to read it here
    // TODO: should we store this information in the checkpoint file instead?
    if (restarting) {
        input_sponge_data.read_from_file(geom[0], zlevels_stag[0]);
    }

    const Real* z_inp_sponge     = input_sponge_data.z_inp_sponge.dataPtr();
    const Real* U_inp_sponge     = input_sponge_data.U_inp_sponge.dataPtr();
    const Real* V_inp_sponge     = input_sponge_data.V_inp_sponge.dataPtr();
    const int   inp_sponge_size  = input_sponge_data.size();

    for (int lev = 0; lev <= finest_level; lev++)
    {
        const int khi = geom[lev].Domain().bigEnd()[2];
        Vector<Real> zcc(khi+1);

        if (z_phys_cc[lev]) {
            // use_terrain=1
            // calculate the damping strength based on the max height at each k
            reduce_to_max_per_height(zcc, z_phys_cc[lev]);
        } else {
            const auto *const prob_lo = geom[lev].ProbLo();
            const auto *const dx = geom[lev].CellSize();
            for (int k = 0; k <= khi; k++)
            {
                zcc[k] = prob_lo[2] + (k+0.5) * dx[2];
            }
        }

        for (int k = 0; k <= khi; k++)
        {
            h_sponge_ptrs[lev][Sponge::ubar_sponge][k] = interpolate_1d(z_inp_sponge, U_inp_sponge, zcc[k], inp_sponge_size);
            h_sponge_ptrs[lev][Sponge::vbar_sponge][k] = interpolate_1d(z_inp_sponge, V_inp_sponge, zcc[k], inp_sponge_size);
        }

        // Copy from host version to device version
        Gpu::copy(Gpu::hostToDevice, h_sponge_ptrs[lev][Sponge::ubar_sponge].begin(), h_sponge_ptrs[lev][Sponge::ubar_sponge].end(),
                         d_sponge_ptrs[lev][Sponge::ubar_sponge].begin());
        Gpu::copy(Gpu::hostToDevice, h_sponge_ptrs[lev][Sponge::vbar_sponge].begin(), h_sponge_ptrs[lev][Sponge::vbar_sponge].end(),
                         d_sponge_ptrs[lev][Sponge::vbar_sponge].begin());
    }
}
