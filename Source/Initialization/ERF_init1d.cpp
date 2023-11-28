/**
 * \file ERF_init1d.cpp
 */
#include <EOS.H>
#include <ERF.H>
#include <TileNoZ.H>
#include <prob_common.H>

using namespace amrex;

/**
 * Initialization function for host and device vectors
 * used to store averaged quantities when calculating
 * the effects of Rayleigh Damping.
 */
void
ERF::initRayleigh ()
{
    AMREX_ALWAYS_ASSERT(solverChoice.use_rayleigh_damping);

    h_rayleigh_tau.resize(max_level+1, amrex::Vector<Real>(0));
    h_rayleigh_ubar.resize(max_level+1, amrex::Vector<Real>(0));
    h_rayleigh_vbar.resize(max_level+1, amrex::Vector<Real>(0));
    h_rayleigh_wbar.resize(max_level+1, amrex::Vector<Real>(0));
    h_rayleigh_thetabar.resize(max_level+1, amrex::Vector<Real>(0));
    d_rayleigh_tau.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    d_rayleigh_ubar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    d_rayleigh_vbar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    d_rayleigh_wbar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    d_rayleigh_thetabar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));

    for (int lev = 0; lev <= finest_level; lev++)
    {
        const int zlen_rayleigh = geom[lev].Domain().length(2);
        h_rayleigh_tau[lev].resize(zlen_rayleigh, 0.0_rt);
        d_rayleigh_tau[lev].resize(zlen_rayleigh, 0.0_rt);
        h_rayleigh_ubar[lev].resize(zlen_rayleigh, 0.0_rt);
        d_rayleigh_ubar[lev].resize(zlen_rayleigh, 0.0_rt);
        h_rayleigh_vbar[lev].resize(zlen_rayleigh, 0.0_rt);
        d_rayleigh_vbar[lev].resize(zlen_rayleigh, 0.0_rt);
        h_rayleigh_wbar[lev].resize(zlen_rayleigh, 0.0_rt);
        d_rayleigh_wbar[lev].resize(zlen_rayleigh, 0.0_rt);
        h_rayleigh_thetabar[lev].resize(zlen_rayleigh, 0.0_rt);
        d_rayleigh_thetabar[lev].resize(zlen_rayleigh, 0.0_rt);

        prob->erf_init_rayleigh(h_rayleigh_tau[lev], h_rayleigh_ubar[lev], h_rayleigh_vbar[lev],
                          h_rayleigh_wbar[lev], h_rayleigh_thetabar[lev], geom[lev]);

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_tau[lev].begin(), h_rayleigh_tau[lev].end(),
                         d_rayleigh_tau[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_ubar[lev].begin(), h_rayleigh_ubar[lev].end(),
                         d_rayleigh_ubar[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_vbar[lev].begin(), h_rayleigh_vbar[lev].end(),
                         d_rayleigh_vbar[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_wbar[lev].begin(), h_rayleigh_wbar[lev].end(),
                         d_rayleigh_wbar[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_thetabar[lev].begin(), h_rayleigh_thetabar[lev].end(),
                         d_rayleigh_thetabar[lev].begin());
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
        input_sounding_data.read_from_file(input_sounding_file, geom[0]);
    }

    const Real* z_inp_sound     = input_sounding_data.z_inp_sound.dataPtr();
    const Real* U_inp_sound     = input_sounding_data.U_inp_sound.dataPtr();
    const Real* V_inp_sound     = input_sounding_data.V_inp_sound.dataPtr();
    const Real* theta_inp_sound = input_sounding_data.theta_inp_sound.dataPtr();
    const int   inp_sound_size  = input_sounding_data.size();

    for (int lev = 0; lev <= finest_level; lev++)
    {
        const int khi = geom[lev].Domain().bigEnd()[2];
        const auto *const prob_lo = geom[lev].ProbLo();
        const auto *const dx = geom[lev].CellSize();

        for (int k = 0; k <= khi; k++)
        {
            const Real z = prob_lo[2] + (k + 0.5) * dx[2];
            h_rayleigh_ubar[lev][k]     = interpolate_1d(z_inp_sound, U_inp_sound, z, inp_sound_size);
            h_rayleigh_vbar[lev][k]     = interpolate_1d(z_inp_sound, V_inp_sound, z, inp_sound_size);
            h_rayleigh_wbar[lev][k]     = 0.0;
            h_rayleigh_thetabar[lev][k] = interpolate_1d(z_inp_sound, theta_inp_sound, z, inp_sound_size);
            if (h_rayleigh_tau[lev][k] > 0) {
                amrex::Print() << z << ":" << " tau=" << h_rayleigh_tau[lev][k];
                if (solverChoice.rayleigh_damp_U) amrex::Print() << " ubar=" << h_rayleigh_ubar[lev][k];
                if (solverChoice.rayleigh_damp_V) amrex::Print() << " vbar=" << h_rayleigh_vbar[lev][k];
                if (solverChoice.rayleigh_damp_W) amrex::Print() << " wbar=" << h_rayleigh_wbar[lev][k];
                if (solverChoice.rayleigh_damp_T) amrex::Print() << " thetabar=" << h_rayleigh_thetabar[lev][k];
                amrex::Print() << std::endl;
            }
        }

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_ubar[lev].begin(), h_rayleigh_ubar[lev].end(),
                         d_rayleigh_ubar[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_vbar[lev].begin(), h_rayleigh_vbar[lev].end(),
                         d_rayleigh_vbar[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_wbar[lev].begin(), h_rayleigh_wbar[lev].end(),
                         d_rayleigh_wbar[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_thetabar[lev].begin(), h_rayleigh_thetabar[lev].end(),
                         d_rayleigh_thetabar[lev].begin());
    }
}

/**
 * Initialize density and pressure base state in
 * hydrostatic equilibrium.
 */
void
ERF::initHSE (int lev)
{
    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    // Initial r_hse may or may not be in HSE -- defined in prob.cpp
    if(solverChoice.use_moist_background){
        prob->erf_init_dens_hse_moist(r_hse, z_phys_nd[lev], z_phys_cc[lev], geom[lev]);
    }else{
        prob->erf_init_dens_hse(r_hse, z_phys_nd[lev], z_phys_cc[lev], geom[lev]);
    }
    // This integrates up through column to update p_hse, pi_hse;
    // r_hse is not const b/c FillBoundary is called at the end for r_hse and p_hse
    erf_enforce_hse(lev, r_hse, p_hse, pi_hse, z_phys_cc[lev], z_phys_nd[lev]);

}

void
ERF::initHSE ()
{
    AMREX_ALWAYS_ASSERT(!init_sounding_ideal);
    for (int lev = 0; lev <= finest_level; lev++)
    {
        initHSE(lev);
    }
}

/**
 * Enforces hydrostatic equilibrium when using terrain.
 *
 * @param[in]  lev  Integer specifying the current level
 * @param[out] dens MultiFab storing base state density
 * @param[out] pres MultiFab storing base state pressure
 * @param[out] pi   MultiFab storing base state Exner function
 * @param[in]  z_cc Pointer to MultiFab storing cell centered z-coordinates
 * @param[in]  z_nd Pointer to MultiFab storing node centered z-coordinates
 */
void
ERF::erf_enforce_hse (int lev,
                      MultiFab& dens, MultiFab& pres, MultiFab& pi,
                      std::unique_ptr<MultiFab>& z_cc,
                      std::unique_ptr<MultiFab>& z_nd)
{
    amrex::Real l_gravity = solverChoice.gravity;
    bool l_use_terrain = solverChoice.use_terrain;

    const auto geomdata = geom[lev].data();
    const Real dz = geomdata.CellSize(2);
    int nz = geom[lev].Domain().length(2);

    const Box& domain = geom[lev].Domain();

    for ( MFIter mfi(dens, TileNoZ()); mfi.isValid(); ++mfi )
    {
        // Create a flat box with same horizontal extent but only one cell in vertical
        const Box& tbz = mfi.nodaltilebox(2);
        amrex::Box b2d = tbz; // Copy constructor
        b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
        b2d.setRange(2,0);

        // We integrate to the first cell (and below) by using rho in this cell
        // If gravity == 0 this is constant pressure
        // If gravity != 0, hence this is a wall, this gives gp0 = dens[0] * gravity
        // (dens_hse*gravity would also be dens[0]*gravity because we use foextrap for rho at k = -1)
        // Note ng_pres_hse = 1

       // We start by assuming pressure on the ground is p_0 (in ERF_Constants.H)
       // Note that gravity is positive

        Array4<Real>  rho_arr = dens.array(mfi);
        Array4<Real> pres_arr = pres.array(mfi);
        Array4<Real>   pi_arr =   pi.array(mfi);
        Array4<Real> zcc_arr;
        Array4<Real> znd_arr;
        if (l_use_terrain) {
           zcc_arr = z_cc->array(mfi);
           znd_arr = z_nd->array(mfi);
        }

        const Real rdOcp = solverChoice.rdOcp;

        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            int k0  = 0;
            // Physical height of the terrain at cell center
            Real hz;
            if (l_use_terrain) {
                hz = .125 * ( znd_arr(i,j,0) + znd_arr(i+1,j,0) + znd_arr(i,j+1,0) + znd_arr(i+1,j+1,0)
                             +znd_arr(i,j,1) + znd_arr(i+1,j,1) + znd_arr(i,j+1,1) + znd_arr(i+1,j+1,1) );
            } else {
                hz = 0.5*dz;
            }

            // Set value at surface from Newton iteration for rho
            pres_arr(i,j,k0  ) = p_0 - hz * rho_arr(i,j,k0) * l_gravity;
            pi_arr(i,j,k0  ) = getExnergivenP(pres_arr(i,j,k0  ), rdOcp);

            // Set ghost cell with dz and rho at boundary
            pres_arr(i,j,k0-1) = p_0 + hz * rho_arr(i,j,k0) * l_gravity;
            pi_arr(i,j,k0-1) = getExnergivenP(pres_arr(i,j,k0-1), rdOcp);

            Real dens_interp;
            if (l_use_terrain) {
                for (int k = 1; k <= nz; k++) {
                    Real dz_loc = (zcc_arr(i,j,k) - zcc_arr(i,j,k-1));
                    dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                    pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz_loc * dens_interp * l_gravity;
                    pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
                }
            } else {
                for (int k = 1; k <= nz; k++) {
                    dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                    pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz * dens_interp * l_gravity;
                    pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
                }
            }
        });

        int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0);
        int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1);

        if (pres[mfi].box().smallEnd(0) < domlo_x)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(0,domlo_x-1);
            bx.setBig(0,domlo_x-1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(domlo_x,j,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().bigEnd(0) > domhi_x)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(0,domhi_x+1);
            bx.setBig(0,domhi_x+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(domhi_x,j,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().smallEnd(1) < domlo_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domlo_y-1);
            bx.setBig(1,domlo_y-1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domlo_y,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().bigEnd(1) > domhi_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domhi_y+1);
            bx.setBig(1,domhi_y+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domhi_y,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }
    }

    dens.FillBoundary(geom[lev].periodicity());
    pres.FillBoundary(geom[lev].periodicity());
}
