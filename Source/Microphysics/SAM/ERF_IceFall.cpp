#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include "ERF_SAM.H"
#include "ERF_SAMUtils.H"
#include "ERF_TileNoZ.H"
#include <cmath>

using namespace amrex;

/**
 * Sedimentation of cloud ice (A32)
 */
void SAM::IceFall (const SolverChoice& sc) {

    if (sam_is_no_ice(sc.moisture_type))
      return;

    Real dtn  = dt;
    Real coef = dtn/m_dzmin;

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
            const Real rho_km1 = (k == k_lo) ? rho_array(i,j,k) : rho_array(i,j,k-1);
            const Real rho_k = (k == k_hi + 1) ? rho_array(i,j,k-1) : rho_array(i,j,k);
            const Real qci_km1 = (k == k_lo) ? qci_array(i,j,k) : qci_array(i,j,k-1);
            const Real qci_k = (k == k_hi + 1) ? qci_array(i,j,k-1) : qci_array(i,j,k);
            const SAMFaceState face_state = sam_face_average_state(k, k_lo, k_hi,
                                                                   rho_km1, rho_k,
                                                                   zero, zero,
                                                                   qci_km1, qci_k,
                                                                   zero, zero);
            const Real vt_ice = sam_cloud_ice_terminal_velocity(face_state.qci_avg);

            // NOTE: Fz is the sedimentation flux from the advective operator.
            //       In the terrain-following coordinate system, the z-deriv in
            //       the divergence uses the normal velocity (Omega). However,
            //       there are no u/v components to the sedimentation velocity.
            //       Therefore, we simply end up with a division by detJ when
            //       evaluating the source term: dJinv * (flux_hi - flux_lo) * dzinv.
            fz_array(i,j,k) = face_state.rho_avg * vt_ice * face_state.qci_avg;
        });
    }

    // Compute number of substeps from maximum terminal velocity
    Real wt_max;
    int n_substep;
    auto const& ma_fz_arr = fz.const_arrays();
    GpuTuple<Real> max = ParReduce(TypeList<ReduceOpMax>{},
                                   TypeList<Real>{},
                                   fz, IntVect(0),
                         [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                         -> GpuTuple<Real>
                         {
                             return { ma_fz_arr[box_no](i,j,k) };
                         });
    wt_max = get<0>(max);
    // ParReduce over a MultiFab gives the local rank maximum here. IceFall
    // needs one global substep count across ranks, matching PrecipFall.
    ParallelDescriptor::ReduceRealMax(wt_max);
    wt_max += std::numeric_limits<Real>::epsilon();
    n_substep = sam_substep_count_from_reduced_flux(wt_max, dtn, m_dzmin);
    AMREX_ALWAYS_ASSERT(n_substep >= 1);
    coef /= Real(n_substep);
    dtn  /= Real(n_substep);

    // Substep the vertical advection
    for (int nsub(0); nsub<n_substep; ++nsub) {
        for (MFIter mfi(*qci, TileNoZ()); mfi.isValid(); ++mfi) {
            auto qci_array   = qci->array(mfi);
            auto qn_array    = qn->array(mfi);
            auto qt_array    = qt->array(mfi);
            auto rho_array   = rho->array(mfi);
            auto fz_array    = fz.array(mfi);

            const auto dJ_array = (m_detJ_cc) ? m_detJ_cc->const_array(mfi) : Array4<const Real>{};

            const auto& tbx  = mfi.tilebox();
            const auto& tbz = mfi.tilebox(IntVect(0,0,1),IntVect(0));

            // Update vertical flux every substep
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real rho_km1 = (k == k_lo) ? rho_array(i,j,k) : rho_array(i,j,k-1);
                const Real rho_k = (k == k_hi + 1) ? rho_array(i,j,k-1) : rho_array(i,j,k);
                const Real qci_km1 = (k == k_lo) ? qci_array(i,j,k) : qci_array(i,j,k-1);
                const Real qci_k = (k == k_hi + 1) ? qci_array(i,j,k-1) : qci_array(i,j,k);
                const SAMFaceState face_state = sam_face_average_state(k, k_lo, k_hi,
                                                                       rho_km1, rho_k,
                                                                       zero, zero,
                                                                       qci_km1, qci_k,
                                                                       zero, zero);
                const Real vt_ice = sam_cloud_ice_terminal_velocity(face_state.qci_avg);

                // NOTE: Fz is the sedimentation flux from the advective operator.
                //       In the terrain-following coordinate system, the z-deriv in
                //       the divergence uses the normal velocity (Omega). However,
                //       there are no u/v components to the sedimentation velocity.
                //       Therefore, we simply end up with a division by detJ when
                //       evaluating the source term: dJinv * (flux_hi - flux_lo) * dzinv.
                fz_array(i,j,k) = face_state.rho_avg * vt_ice * face_state.qci_avg;
            });

            // Update precip every substep
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Jacobian determinant
                Real dJinv = (dJ_array) ? one/dJ_array(i,j,k) : one;

                //==================================================
                // Cloud ice sedimentation (A32)
                //==================================================
                Real dqi  = sam_sedimentation_tendency(fz_array(i,j,k+1), fz_array(i,j,k),
                                                       rho_array(i,j,k), dJinv, coef);
                dqi = std::max(-qci_array(i,j,k), dqi);

                // Add this increment to both non-precipitating and total water.
                qci_array(i,j,k) += dqi;
                 qn_array(i,j,k) += dqi;
                 qt_array(i,j,k) += dqi;

                // NOTE: Sedimentation does not affect the potential temperature,
                //       but it does affect the liquid/ice static energy.
                //       No source to Theta occurs here.
            });
        } // mfi
    } // nsub
}

