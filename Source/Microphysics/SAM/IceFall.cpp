#include <AMReX_ParReduce.H>
#include "SAM.H"
#include "TileNoZ.H"

using namespace amrex;

/**
 * Sedimentation of cloud ice (A32)
 */
void SAM::IceFall () {

    Real dz   = m_geom.CellSize(2);
    Real dtn  = dt;
    Real coef = dtn/dz;

    auto domain = m_geom.Domain();
    int k_lo = domain.smallEnd(2);
    int k_hi = domain.bigEnd(2);

    auto qcl   = mic_fab_vars[MicVar::qcl];
    auto qci   = mic_fab_vars[MicVar::qci];
    auto qn    = mic_fab_vars[MicVar::qn];
    auto qt    = mic_fab_vars[MicVar::qt];
    auto rho   = mic_fab_vars[MicVar::rho];
    auto tabs  = mic_fab_vars[MicVar::tabs];

    MultiFab fz;
    IntVect  ng = qcl->nGrowVect();
    BoxArray ba = qcl->boxArray();
    DistributionMapping dm = qcl->DistributionMap();
    fz.define(convert(ba, IntVect(0,0,1)), dm, 1, ng);
    fz.setVal(0.);

    for (MFIter mfi(fz, TileNoZ()); mfi.isValid(); ++mfi) {
        auto qci_array = qci->array(mfi);
        auto rho_array = rho->array(mfi);
        auto fz_array  = fz.array(mfi);

        const auto& box3d  = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real rho_avg, qci_avg;
            if (k==k_lo) {
                rho_avg = rho_array(i,j,k);
                qci_avg = qci_array(i,j,k);
            } else if (k==k_hi+1) {
                rho_avg = rho_array(i,j,k-1);
                qci_avg = qci_array(i,j,k-1);
            } else {
                rho_avg = 0.5*(rho_array(i,j,k-1) + rho_array(i,j,k));
                qci_avg = 0.5*(qci_array(i,j,k-1) + qci_array(i,j,k));
            }
            Real vt_ice = min( 0.4 , 8.66 * pow( (max(0.,qci_avg)+1.e-10) , 0.24) );

            // NOTE: Fz is the sedimentation flux from the advective operator.
            //       In the terrain-following coordinate system, the z-deriv in
            //       the divergence uses the normal velocity (Omega). However,
            //       there are no u/v components to the sedimentation velocity.
            //       Therefore, we simply end up with a division by detJ when
            //       evaluating the source term: dJinv * (flux_hi - flux_lo) * dzinv.
            fz_array(i,j,k) = rho_avg*vt_ice*qci_avg;
        });
    }

    for (MFIter mfi(*qci, TileNoZ()); mfi.isValid(); ++mfi) {
        auto qci_array   = qci->array(mfi);
        auto qn_array    = qn->array(mfi);
        auto qt_array    = qt->array(mfi);
        auto rho_array   = rho->array(mfi);
        auto fz_array    = fz.array(mfi);

        const auto dJ_array = (m_detJ_cc) ? m_detJ_cc->const_array(mfi) : Array4<const Real>{};

        const auto& box3d  = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Jacobian determinant
            Real dJinv = (dJ_array) ? 1.0/dJ_array(i,j,k) : 1.0;

            //==================================================
            // Cloud ice sedimentation (A32)
            //==================================================
            Real dqi  = dJinv * (1.0/rho_array(i,j,k)) * ( fz_array(i,j,k+1) - fz_array(i,j,k) ) * coef;
            dqi = std::max(-qci_array(i,j,k), dqi);

            // Add this increment to both non-precipitating and total water.
            qci_array(i,j,k) += dqi;
             qn_array(i,j,k) += dqi;
             qt_array(i,j,k) += dqi;

            // NOTE: Sedimentation does not affect the potential temperature,
            //       but it does affect the liquid/ice static energy.
            //       No source to Theta occurs here.
        });
    }
}

