/**
 * \file ERF_InitRayleigh.cpp
 */
#include <ERF.H>

using namespace amrex;

/**
 * Initialization function for host and device vectors
 * used to store averaged quantities when calculating
 * the effects of Rayleigh Damping.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::initRayleigh_at_level (const int& lev)
{
    const int khi = geom[0].Domain().bigEnd(2);
    solverChoice.dampingChoice.rayleigh_ztop = (solverChoice.terrain_type == TerrainType::None) ?  geom[0].ProbHi(2) : zlevels_stag[0][khi+1];

    // These have 4 components: ubar, vbar, wbar, thetabar
    h_rayleigh_ptrs[lev].resize(Rayleigh::nvars);
    d_rayleigh_ptrs[lev].resize(Rayleigh::nvars);

    const int zlen_rayleigh = geom[lev].Domain().length(2);

    // Allocate space for these 1D vectors
    for (int n = 0; n < Rayleigh::nvars; n++) {
        h_rayleigh_ptrs[lev][n].resize(zlen_rayleigh, 0.0_rt);
        d_rayleigh_ptrs[lev][n].resize(zlen_rayleigh, 0.0_rt);
    }

    h_sinesq_ptrs[lev].resize(zlen_rayleigh);
    d_sinesq_ptrs[lev].resize(zlen_rayleigh);

    h_sinesq_stag_ptrs[lev].resize(zlen_rayleigh+1);
    d_sinesq_stag_ptrs[lev].resize(zlen_rayleigh+1);

    Real ztop     = solverChoice.dampingChoice.rayleigh_ztop;
    Real zdamp    = solverChoice.dampingChoice.rayleigh_zdamp;

    for (int k = 0; k < zlen_rayleigh; k++) {
        Real z = myhalf * (zlevels_stag[lev][k] + zlevels_stag[lev][k+1]);
        if (z > (ztop - zdamp)) {
            Real zfrac = one - (ztop - z) / zdamp;
            Real s = std::sin(PIoTwo*zfrac);
            h_sinesq_ptrs[lev][k] = s*s;
        } else {
            h_sinesq_ptrs[lev][k] = zero;
        }
    }

    for (int k = 0; k < zlen_rayleigh+1; k++) {
        Real z = zlevels_stag[lev][k];
        if (z > (ztop - zdamp)) {
            Real zfrac = one - (ztop - z) / zdamp;
            Real s = std::sin(PIoTwo*zfrac);
            h_sinesq_stag_ptrs[lev][k] = s*s;
        } else {
            h_sinesq_stag_ptrs[lev][k] = zero;
        }
    }

    // Init the host vectors for the reference states
    prob->erf_init_rayleigh(h_rayleigh_ptrs[lev], geom[lev], z_phys_nd[lev],
                            solverChoice.dampingChoice.rayleigh_zdamp);

    // With an input sounding, ubar, vbar and thetabar come from it (wbar is zero); the
    // damping coefficient set by erf_init_rayleigh is still used.  Doing it here, rather
    // than once at start-up, keeps the sounding on a level a regrid remakes or creates.
    // Every caller (InitData_post, and RemakeLevel / MakeNewLevelFromCoarse, which run only
    // from regrid after InitData_post) comes after the sounding is loaded, so a missing
    // sounding is an error rather than a level to leave with the problem defaults.
    if (solverChoice.init_type == InitType::Input_Sounding) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!input_sounding_data.z_inp_sound.empty() &&
                                         !input_sounding_data.z_inp_sound[0].empty(),
            "initRayleigh_at_level: the input sounding must be loaded before the Rayleigh "
            "reference profiles are set");
        setRayleighRefFromSounding_at_level(lev);
    }

    // Copy from host vectors to device vectors
    for (int n = 0; n < Rayleigh::nvars; n++) {
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][n].begin(), h_rayleigh_ptrs[lev][n].end(),
                  d_rayleigh_ptrs[lev][n].begin());
    }
    Gpu::copy(Gpu::hostToDevice, h_sinesq_ptrs[lev].begin(), h_sinesq_ptrs[lev].end(), d_sinesq_ptrs[lev].begin());
    Gpu::copy(Gpu::hostToDevice, h_sinesq_stag_ptrs[lev].begin(), h_sinesq_stag_ptrs[lev].end(), d_sinesq_stag_ptrs[lev].begin());
}

/**
 * Sets the Rayleigh Damping reference profiles at one level from the input sounding
 * (sounding 0), at the cell-center heights of that level.  The sounding is the one
 * sampled at this level's own cell centers (sounding_at_level), as for the level's
 * initial state, so that the damping targets are the profile the level starts from.
 * Called from initRayleigh_at_level, before the profiles are copied to the device.  The
 * sounding is already loaded: at start-up by init_only, on a restart by
 * MakeNewLevelFromScratch.
 *
 * @param[in] lev Level whose reference profiles to set
 */
void
ERF::setRayleighRefFromSounding_at_level (int lev)
{
    const InputSoundingData& snd = sounding_at_level(lev);
    const Real* z_inp_sound     = snd.z_inp_sound[0].dataPtr();
    const Real* U_inp_sound     = snd.U_inp_sound[0].dataPtr();
    const Real* V_inp_sound     = snd.V_inp_sound[0].dataPtr();
    const Real* theta_inp_sound = snd.theta_inp_sound[0].dataPtr();
    const int   inp_sound_size  = snd.size(0);

    const int Nz = geom[lev].Domain().length(2);
    for (int k = 0; k < Nz; k++)
    {
        const Real zcc = myhalf * (zlevels_stag[lev][k] + zlevels_stag[lev][k+1]);
        h_rayleigh_ptrs[lev][Rayleigh::ubar][k]     = interpolate_1d(z_inp_sound, U_inp_sound, zcc, inp_sound_size);
        h_rayleigh_ptrs[lev][Rayleigh::vbar][k]     = interpolate_1d(z_inp_sound, V_inp_sound, zcc, inp_sound_size);
        h_rayleigh_ptrs[lev][Rayleigh::wbar][k]     = zero;
        h_rayleigh_ptrs[lev][Rayleigh::thetabar][k] = interpolate_1d(z_inp_sound, theta_inp_sound, zcc, inp_sound_size);
    }
}
