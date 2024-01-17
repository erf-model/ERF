#include <SLM.H>

using namespace amrex;

/* Initialize lsm data structures */
void
SLM::Init (const MultiFab& cons_in,
           const Geometry& geom,
           const Real& dt)
{
    m_dt = dt;
    m_geom = geom;

    Box domain = geom.Domain();
    khi_lsm    = domain.smallEnd(2) - 1;

    LsmVarMap.resize(m_lsm_size);
    LsmVarMap = {LsmVar_SLM::theta};

    for (auto ivar = 0; ivar < LsmVar_SLM::NumVars; ++ivar) {
        IntVect ng(0,0,1);
        BoxArray ba = cons_in.boxArray();
        DistributionMapping dm = cons_in.DistributionMap();

        // NOTE: All boxes in ba extend from zlo to zhi, so this transform is valid.
        //       If that were to change, the dm and new ba are no longer valid and
        //       direct copying between lsm data/flux vars cannot be done in a parfor.
        BoxList bl_lsm = ba.boxList();
        for (auto& b : bl_lsm) {
            b.setBig(2,khi_lsm);              // First point below the surface
            b.setSmall(2,khi_lsm - m_nz + 1); // Last point below the surface
        }
        BoxArray ba_lsm(std::move(bl_lsm));

        // State vars are CC
        Real theta_0 = m_theta_dir;
        lsm_fab_vars[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_vars[ivar]->setVal(theta_0);

        // Fluxes are nodal in z
        lsm_fab_flux[ivar] = std::make_shared<MultiFab>(convert(ba_lsm, IntVect(0,0,1)), dm, 1, IntVect(0,0,0));
        lsm_fab_flux[ivar]->setVal(0.);
    }
}

/* Extrapolate surface temperature and store in ghost cell */
void
SLM::ComputeTsurf ()
{
    // Expose for GPU copy
    int khi = khi_lsm;

    for ( MFIter mfi(*(lsm_fab_vars[LsmVar_SLM::theta])); mfi.isValid(); ++mfi) {
        auto box2d = mfi.tilebox(); box2d.makeSlab(2,khi);

        auto theta_array = lsm_fab_vars[LsmVar_SLM::theta]->array(mfi);

        ParallelFor( box2d, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            theta_array(i,j,khi+1) = 1.5*theta_array(i,j,khi) - 0.5*theta_array(i,j,khi-1);
        });
    }
}

/* Compute the diffusive fluxes */
void
SLM::ComputeFluxes ()
{
    // Expose for GPU copy
    int khi = khi_lsm;
    Real Dsoil = m_d_soil;
    Real dzInv = 1.0/dz;

    for ( MFIter mfi(*(lsm_fab_flux[LsmVar_SLM::theta])); mfi.isValid(); ++mfi) {
        auto box3d = mfi.tilebox();

        // Do not overwrite the flux at the top (comes from MOST BC)
        if (box3d.bigEnd(2) == khi+1) box3d.setBig(2,khi);

        auto theta_array = lsm_fab_vars[LsmVar_SLM::theta]->array(mfi);
        auto theta_flux  = lsm_fab_flux[LsmVar_SLM::theta]->array(mfi);

        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            theta_flux(i,j,k) = Dsoil * ( theta_array(i,j,k) - theta_array(i,j,k-1) ) * dzInv;
        });
    }
}

/* Advance the solution with a simple explicit update */
void
SLM::AdvanceSLM ()
{
    // Expose for GPU copy
    Real dt = m_dt;
    Real dzInv = 1.0/dz;

    for ( MFIter mfi(*(lsm_fab_vars[LsmVar_SLM::theta])); mfi.isValid(); ++mfi) {
        auto box3d = mfi.tilebox();

        auto theta_array = lsm_fab_vars[LsmVar_SLM::theta]->array(mfi);
        auto theta_flux  = lsm_fab_flux[LsmVar_SLM::theta]->array(mfi);

        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            theta_array(i,j,k) += dt * ( theta_flux(i,j,k+1) - theta_flux(i,j,k) ) * dzInv;
        });
    }
}
