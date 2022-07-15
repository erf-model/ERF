#include <ERF.H>
#include <EOS.H>
#include <prob_common.H>

using namespace amrex;

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
    if (solverChoice.use_terrain) {
        amrex::Print() << "Using terrain in initHSE" << std::endl;
        for (int lev = 0; lev <= finest_level; lev++)
        {
            MultiFab* z_nd = nullptr;
            MultiFab* z_cc = nullptr;
            //if (solverChoice.use_terrain) {
            z_nd = &z_phys_nd[lev];
            z_cc = &z_phys_cc[lev];
            //}
            
            amrex::Print() << "I'm in initHSE() before erf_init_dens_hse()!" << std::endl;
            erf_init_dens_hse(dens_hse[lev], z_nd, z_cc, geom[lev], solverChoice.use_terrain);
            amrex::Print() << "I'm in initHSE() after erf_init_dens_hse()!" << std::endl;
            
            erf_enforce_hse(lev, dens_hse[lev], pres_hse[lev], z_phys_cc[lev]);
        }
    }
    else {
        amrex::Print() << "NOT using terrain in initHSE" << std::endl;

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
    }
}

/*
void
ERF::erf_enforce_hse(int lev, MultiFab& dens, MultiFab& pres)
{
    amrex::Real l_gravity = solverChoice.gravity;

    int nz = geom[lev].Domain().length(2);

    const Box& domain = geom[lev].Domain();

    bool l_use_terrain = solverChoice.use_terrain;

    const auto geomdata = geom[lev].data();
    const Real dz = geomdata.CellSize(2);

    for ( MFIter mfi(dens, TilingIfNotGPU()); mfi.isValid(); ++mfi )
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
        Array4<Real const> zcc_arr  = l_use_terrain ? z_phys_cc[lev].const_array(mfi) : Array4<Real const>{};
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            int k0  = 0;
            // Physical height of the terrain at cell center
            Real hz;
            if (l_use_terrain) 
                hz = zcc_arr(i,j,k0);
            else
                hz = k0 * dz;

            // Set value at surface from Newton iteration for rho
            pres_arr(i,j,k0  ) = p_0 - hz * rho_arr(i,j,k0) * l_gravity;

            // Set ghost cell with dz and rho at boundary
            pres_arr(i,j,k0-1) = p_0 + hz * rho_arr(i,j,k0) * l_gravity;

            Real dens_interp;
            for (int k = 1; k <= nz; k++) {
                Real dz_loc;
                if (l_use_terrain) {
                    dz_loc = (zcc_arr(i,j,k) - zcc_arr(i,j,k-1));
                } else {
                    dz_loc = dz;
                }
                dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz_loc * dens_interp * l_gravity;
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
            });
        }

        if (pres[mfi].box().bigEnd(0) > domhi_x)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(0,domhi_x+1);
            bx.setBig(0,domhi_x+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(domhi_x,j,k);
            });
        }

        if (pres[mfi].box().smallEnd(1) < domlo_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domlo_y-1);
            bx.setBig(1,domlo_y-1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domlo_y,k);
            });
        }

        if (pres[mfi].box().bigEnd(1) > domhi_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domhi_y+1);
            bx.setBig(1,domhi_y+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domhi_y,k);
            });
        }
    }
    dens.FillBoundary(geom[lev].periodicity());
    pres.FillBoundary(geom[lev].periodicity());
}
*/
void
ERF::erf_enforce_hse(int lev,
                     MultiFab& dens, MultiFab& pres,
                     MultiFab& z_cc)
{
    amrex::Real l_gravity = solverChoice.gravity;

    int nz = geom[lev].Domain().length(2);

    const Box& domain = geom[lev].Domain();

    for ( MFIter mfi(dens, TilingIfNotGPU()); mfi.isValid(); ++mfi )
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
        Array4<Real> zcc_arr  = z_cc.array(mfi);
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            int k0  = 0;
            // Physical height of the terrain at cell center
            Real hz = zcc_arr(i,j,k0);

            // Set value at surface from Newton iteration for rho
            pres_arr(i,j,k0  ) = p_0 - hz * rho_arr(i,j,k0) * l_gravity;

            // Set ghost cell with dz and rho at boundary
            pres_arr(i,j,k0-1) = p_0 + hz * rho_arr(i,j,k0) * l_gravity;

            Real dens_interp;
            for (int k = 1; k <= nz; k++) {
                Real dz_loc = (zcc_arr(i,j,k) - zcc_arr(i,j,k-1));
                dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz_loc * dens_interp * l_gravity;
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
            });
        }

        if (pres[mfi].box().bigEnd(0) > domhi_x)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(0,domhi_x+1);
            bx.setBig(0,domhi_x+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(domhi_x,j,k);
            });
        }

        if (pres[mfi].box().smallEnd(1) < domlo_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domlo_y-1);
            bx.setBig(1,domlo_y-1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domlo_y,k);
            });
        }

        if (pres[mfi].box().bigEnd(1) > domhi_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domhi_y+1);
            bx.setBig(1,domhi_y+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domhi_y,k);
            });
        }
    }
    dens.FillBoundary(geom[lev].periodicity());
    pres.FillBoundary(geom[lev].periodicity());
}

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
