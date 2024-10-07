/**
 * \file ERF_init_rayleigh.cpp
 */
#include <ERF.H>

using namespace amrex;

/**
 * Initialization function for host and device vectors
 * used to store averaged quantities when calculating
 * the effects of Rayleigh Damping.
 */
void
ERF::initRayleigh ()
{
    const int khi = geom[0].Domain().bigEnd(2);
    solverChoice.rayleigh_ztop = (solverChoice.use_terrain) ? zlevels_stag[0][khi+1] : geom[0].ProbHi(2);

    h_rayleigh_ptrs.resize(max_level+1);
    d_rayleigh_ptrs.resize(max_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        // These have 4 components: ubar, vbar, wbar, thetabar
        h_rayleigh_ptrs[lev].resize(Rayleigh::nvars);
        d_rayleigh_ptrs[lev].resize(Rayleigh::nvars);

        const int zlen_rayleigh = geom[lev].Domain().length(2);

        // Allocate space for these 1D vectors
        for (int n = 0; n < Rayleigh::nvars; n++) {
            h_rayleigh_ptrs[lev][n].resize(zlen_rayleigh, 0.0_rt);
            d_rayleigh_ptrs[lev][n].resize(zlen_rayleigh, 0.0_rt);
        }

        // Init the host vectors
        prob->erf_init_rayleigh(h_rayleigh_ptrs[lev], geom[lev], z_phys_nd[lev], solverChoice.rayleigh_zdamp);

        // Copy from host vectors to device vectors
        for (int n = 0; n < Rayleigh::nvars; n++) {
            Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][n].begin(), h_rayleigh_ptrs[lev][n].end(),
                      d_rayleigh_ptrs[lev][n].begin());
        }
    }
}

/**
 * Sets the Rayleigh Damping averaged quantities from an
 * externally supplied input sounding data file.
 *
 * @param[in] restarting Boolean parameter that indicates whether
                         we are currently restarting from a checkpoint file.
 */
void
ERF::setRayleighRefFromSounding (bool restarting)
{
    // If we are restarting then we haven't read the input_sounding file yet
    //    so we need to read it here
    // TODO: should we store this information in the checkpoint file instead?
    if (restarting) {
        input_sounding_data.resize_arrays();
        for (int n = 0; n < input_sounding_data.n_sounding_files; n++) {
            input_sounding_data.read_from_file(geom[0], zlevels_stag[0], n);
        }
    }

    const Real* z_inp_sound     = input_sounding_data.z_inp_sound[0].dataPtr();
    const Real* U_inp_sound     = input_sounding_data.U_inp_sound[0].dataPtr();
    const Real* V_inp_sound     = input_sounding_data.V_inp_sound[0].dataPtr();
    const Real* theta_inp_sound = input_sounding_data.theta_inp_sound[0].dataPtr();
    const int   inp_sound_size  = input_sounding_data.size(0);

    int refine_fac{1};
    for (int lev = 0; lev <= finest_level; lev++)
    {
        const int klo = geom[lev].Domain().smallEnd(2);
        const int khi = geom[lev].Domain().bigEnd(2);
        const int Nz = khi - klo + 1;

        Vector<Real> zcc(Nz);
        Vector<Real> zlevels_sub(zlevels_stag[0].begin()+klo/refine_fac,
                                 zlevels_stag[0].begin()+khi/refine_fac+2);
        expand_and_interpolate_1d(zcc, zlevels_sub, refine_fac, true);
#if 0
        amrex::AllPrint() << "lev="<<lev<<" : (refine_fac="<<refine_fac<<",klo="<<klo<<",khi="<<khi<<") ";
        for (int k = 0; k < zlevels_sub.size(); k++) { amrex::AllPrint() << zlevels_sub[k] << " "; }
        amrex::AllPrint() << " --> ";
        for (int k = 0; k < Nz; k++) { amrex::AllPrint() << zcc[k] << " "; }
        amrex::AllPrint() << std::endl;
#endif

        for (int k = 0; k < Nz; k++)
        {
            h_rayleigh_ptrs[lev][Rayleigh::ubar][k]     = interpolate_1d(z_inp_sound, U_inp_sound, zcc[k], inp_sound_size);
            h_rayleigh_ptrs[lev][Rayleigh::vbar][k]     = interpolate_1d(z_inp_sound, V_inp_sound, zcc[k], inp_sound_size);
            h_rayleigh_ptrs[lev][Rayleigh::wbar][k]     = Real(0.0);
            h_rayleigh_ptrs[lev][Rayleigh::thetabar][k] = interpolate_1d(z_inp_sound, theta_inp_sound, zcc[k], inp_sound_size);
        }

        // Copy from host version to device version
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][Rayleigh::ubar].begin(), h_rayleigh_ptrs[lev][Rayleigh::ubar].end(),
                         d_rayleigh_ptrs[lev][Rayleigh::ubar].begin());
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][Rayleigh::vbar].begin(), h_rayleigh_ptrs[lev][Rayleigh::vbar].end(),
                         d_rayleigh_ptrs[lev][Rayleigh::vbar].begin());
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][Rayleigh::wbar].begin(), h_rayleigh_ptrs[lev][Rayleigh::wbar].end(),
                         d_rayleigh_ptrs[lev][Rayleigh::wbar].begin());
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][Rayleigh::thetabar].begin(), h_rayleigh_ptrs[lev][Rayleigh::thetabar].end(),
                         d_rayleigh_ptrs[lev][Rayleigh::thetabar].begin());

        refine_fac *= ref_ratio[lev][2];
    }
}
