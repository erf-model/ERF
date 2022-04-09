#include <ERF.H>
#include <prob_common.H>

void
ERF::initRayleigh()
{
    AMREX_ALWAYS_ASSERT(solverChoice.use_rayleigh_damping);

    h_rayleigh_tau.resize(max_level+1, amrex::Vector<Real>(0));
    h_rayleigh_ubar.resize(max_level+1, amrex::Vector<Real>(0));
    h_rayleigh_vbar.resize(max_level+1, amrex::Vector<Real>(0));
    h_rayleigh_thetabar.resize(max_level+1, amrex::Vector<Real>(0));
    d_rayleigh_tau.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    d_rayleigh_ubar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    d_rayleigh_vbar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
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
        h_rayleigh_thetabar[lev].resize(zlen_rayleigh, 0.0_rt);
        d_rayleigh_thetabar[lev].resize(zlen_rayleigh, 0.0_rt);

        erf_init_rayleigh(h_rayleigh_tau[lev], h_rayleigh_ubar[lev], h_rayleigh_vbar[lev],
                          h_rayleigh_thetabar[lev], geom[lev]);

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_tau[lev].begin(), h_rayleigh_tau[lev].end(),
                         d_rayleigh_tau[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_ubar[lev].begin(), h_rayleigh_ubar[lev].end(),
                         d_rayleigh_ubar[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_vbar[lev].begin(), h_rayleigh_vbar[lev].end(),
                         d_rayleigh_vbar[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_rayleigh_thetabar[lev].begin(), h_rayleigh_thetabar[lev].end(),
                         d_rayleigh_thetabar[lev].begin());
    }
}

void
ERF::initHSE()
{
#ifdef ERF_USE_TERRAIN
    for (int lev = 0; lev <= finest_level; lev++)
    {
        pres_hse[lev].define(grids[lev],dmap[lev],1,0);
        dens_hse[lev].define(grids[lev],dmap[lev],1,0);

        for ( MFIter mfi(pres_hse[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Create a flat box with same horizontal extent but only one cell in vertical
            const Box& tbz = mfi.nodaltilebox(2);
            amrex::Box b2d = tbz; // Copy constructor
            int klo = tbz.smallEnd(2);
            int khi = tbz.bigEnd(2);
            b2d.setRange(2,0);

            Array4<Real> rho_arr = dens_hse[lev].array(mfi);
            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                rho_arr(i,j,k) = 1.0;
            });
        }
        erf_enforce_hse(lev, dens_hse[lev], pres_hse[lev]);
    }
#else
    //
    // Setup Base State Arrays
    //
    h_dens_hse.resize(max_level+1, amrex::Vector<Real>(0));
    h_pres_hse.resize(max_level+1, amrex::Vector<Real>(0));
    d_dens_hse.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    d_pres_hse.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));

    for (int lev = 0; lev <= finest_level; lev++)
    {
        const int zlen_dens = geom[lev].Domain().length(2) + 2*ng_dens_hse;
        h_dens_hse[lev].resize(zlen_dens, 0.0_rt);
        d_dens_hse[lev].resize(zlen_dens, 0.0_rt);

        const int zlen_pres = geom[lev].Domain().length(2) + 2*ng_pres_hse;
        h_pres_hse[lev].resize(zlen_pres, p_0);
        d_pres_hse[lev].resize(zlen_pres, p_0);

        Real* hptr_dens = h_dens_hse[lev].data() + ng_dens_hse;

        erf_init_dens_hse(hptr_dens,geom[lev],ng_dens_hse);

        erf_enforce_hse(lev,h_dens_hse[lev],h_pres_hse[lev]);

        // Copy from host version to device version
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_dens_hse[lev].begin(), h_dens_hse[lev].end(),
                         d_dens_hse[lev].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_pres_hse[lev].begin(), h_pres_hse[lev].end(),
                     d_pres_hse[lev].begin());
    }
#endif
}

#ifdef ERF_USE_TERRAIN
void
ERF::erf_enforce_hse(int lev,
                     MultiFab& dens,
                     MultiFab& pres)
{
    amrex::Real l_gravity = solverChoice.gravity;

    const auto geomdata = geom[lev].data();
    const Real dz = geomdata.CellSize(2);
    int nz = geom[lev].Domain().length(2);

    for ( MFIter mfi(dens, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Create a flat box with same horizontal extent but only one cell in vertical
        const Box& tbz = mfi.nodaltilebox(2);
        amrex::Box b2d = tbz; // Copy constructor
        int klo = tbz.smallEnd(2);
        int khi = tbz.bigEnd(2);
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
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            pres_arr(i,j,-1) = p_0 + (0.5*dz) * rho_arr(i,j,0) * l_gravity;
            pres_arr(i,j, 0) = p_0 - (0.5*dz) * rho_arr(i,j,0) * l_gravity;
            for (int k = 1; k <= nz; k++) {
                Real dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz * dens_interp * l_gravity;
            }
        });
    }
}
#else
void
ERF::erf_enforce_hse(int lev,
                     amrex::Vector<amrex::Real>& dens,
                     amrex::Vector<amrex::Real>& pres)
{
    AMREX_ALWAYS_ASSERT(dens.size() == pres.size());

    amrex::Real l_gravity = solverChoice.gravity;

    const auto geomdata = geom[lev].data();
    const Real dz = geomdata.CellSize(2);
    int nz = geom[lev].Domain().length(2);

    // We start by assuming pressure on the ground is p_0 (in ERF_Constants.H)
    // Note that gravity is positive

    Real* hptr_dens = h_dens_hse[lev].data() + ng_dens_hse;
    Real* hptr_pres = h_pres_hse[lev].data() + ng_pres_hse;

    Real dens_interp;

    // We integrate to the first cell (and below) by using rho in this cell
    // If gravity == 0 this is constant pressure
    // If gravity != 0, hence this is a wall, this gives gp0 = dens[0] * gravity
    // (dens_hse*gravity would also be dens[0]*gravity because we use foextrap for rho at k = -1)
    // Note ng_pres_hse = 1
    hptr_pres[-1] = p_0 + (0.5*dz) * dens[0] * l_gravity;
    hptr_pres[ 0] = p_0 - (0.5*dz) * dens[0] * l_gravity;
    //amrex::Print() << "erf_enforce_hse: p[-1] = " << hptr_pres[-1] << " (ghost)" << std::endl;
    //amrex::Print() << "erf_enforce_hse: p[ 0] = " << hptr_pres[ 0] << std::endl;

    for (int k = 1; k < nz+ng_pres_hse; k++)
    {
       dens_interp = 0.5*(hptr_dens[k] + hptr_dens[k-1]);
       hptr_pres[k] = hptr_pres[k-1] - dz * dens_interp * l_gravity;
    }
}
#endif
