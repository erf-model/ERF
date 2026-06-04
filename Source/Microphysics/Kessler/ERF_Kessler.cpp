#include <ERF_EOS.H>
#include <ERF_TileNoZ.H>

#include "ERF_DataStruct.H"
#include "ERF_Kessler.H"
#include "ERF_KesslerUtils.H"

using namespace amrex;

/**
 * Compute Precipitation-related Microphysics quantities.
 */

void Kessler::AdvanceKessler (const SolverChoice &solverChoice)
{
    bool do_cond = m_do_cond;
    auto tabs    = mic_fab_vars[MicVar_Kess::tabs];
    auto domain  = m_geom.Domain();
    int k_lo = domain.smallEnd(2);
    int k_hi = domain.bigEnd(2);
    if (solverChoice.moisture_type == MoistureType::Kessler)
    {
        MultiFab fz;
        auto ba    = tabs->boxArray();
        auto dm    = tabs->DistributionMap();
        fz.define(convert(ba, IntVect(0,0,1)), dm, 1, 0); // No ghost cells

        Real dtn  = dt;
        Real coef = dtn/m_dzmin;

        for ( MFIter mfi(*tabs,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto qv_array    = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
            auto qc_array    = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);
            auto qp_array    = mic_fab_vars[MicVar_Kess::qp]->array(mfi);
            auto qt_array    = mic_fab_vars[MicVar_Kess::qt]->array(mfi);
            auto tabs_array  = mic_fab_vars[MicVar_Kess::tabs]->array(mfi);
            auto pres_array  = mic_fab_vars[MicVar_Kess::pres]->array(mfi);
            auto theta_array = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
            auto rho_array   = mic_fab_vars[MicVar_Kess::rho]->array(mfi);

            auto tbx = mfi.tilebox();

            Real d_fac_cond = m_fac_cond;

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                qv_array(i,j,k) = std::max(Real(0), qv_array(i,j,k));
                qc_array(i,j,k) = std::max(Real(0), qc_array(i,j,k));
                qp_array(i,j,k) = std::max(Real(0), qp_array(i,j,k));

                Real qsat_local, dtqsat_local;
                Real pressure = pres_array(i,j,k);
                // Kessler stores pressure in mbar / hPa for the qsat helpers.
                erf_qsatw(tabs_array(i,j,k), pressure, qsat_local);
                erf_dtqsatw(tabs_array(i,j,k), pressure, dtqsat_local);

                if (qsat_local <= Real(0)) {
                    amrex::Warning("qsat computed as non-positive; setting to Real(0)!");
                    qsat_local = Real(0);
                }

                const KesslerSourceTerms source_terms = kessler_warm_rain_sources(
                    qv_array(i,j,k), qc_array(i,j,k), qp_array(i,j,k), rho_array(i,j,k),
                    pressure, qsat_local, dtqsat_local, dtn, do_cond, d_fac_cond);

                qv_array(i,j,k) += -source_terms.dq_vapor_to_cloud
                                 +  source_terms.dq_cloud_to_vapor
                                 +  source_terms.dq_rain_to_vapor;
                qc_array(i,j,k) +=  source_terms.dq_vapor_to_cloud
                                 -  source_terms.dq_cloud_to_vapor
                                 -  source_terms.dq_cloud_to_rain;
                qp_array(i,j,k) +=  source_terms.dq_cloud_to_rain
                                 -  source_terms.dq_rain_to_vapor;

                Real theta_over_T = theta_array(i,j,k)/tabs_array(i,j,k);
                theta_array(i,j,k) += theta_over_T * d_fac_cond
                    * (source_terms.dq_vapor_to_cloud
                       - source_terms.dq_cloud_to_vapor
                       - source_terms.dq_rain_to_vapor);

                qv_array(i,j,k) = std::max(Real(0), qv_array(i,j,k));
                qc_array(i,j,k) = std::max(Real(0), qc_array(i,j,k));
                qp_array(i,j,k) = std::max(Real(0), qp_array(i,j,k));

                qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);
            });
        }

        for ( MFIter mfi(fz, TilingIfNotGPU()); mfi.isValid(); ++mfi ){
            auto rho_array = mic_fab_vars[MicVar_Kess::rho]->array(mfi);
            auto qp_array  = mic_fab_vars[MicVar_Kess::qp]->array(mfi);

            auto fz_array  = fz.array(mfi);
            const Box& tbz = mfi.tilebox();

            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real rho_km1 = (k == k_lo) ? rho_array(i,j,k) : rho_array(i,j,k-1);
                const Real rho_k = (k == k_hi+1) ? rho_array(i,j,k-1) : rho_array(i,j,k);
                const Real qp_km1 = (k == k_lo) ? qp_array(i,j,k) : qp_array(i,j,k-1);
                const Real qp_k = (k == k_hi+1) ? qp_array(i,j,k-1) : qp_array(i,j,k);
                const KesslerFaceState face_state =
                    kessler_face_state(k, k_hi, rho_km1, rho_k, qp_km1, qp_k);

                const Real terminal_velocity = kessler_terminal_velocity(face_state.rho, face_state.qp);
                fz_array(i,j,k) = kessler_precip_flux(face_state.rho, terminal_velocity, face_state.qp);
            });
        }

        auto const& ma_rho_arr = mic_fab_vars[MicVar_Kess::rho]->const_arrays();
        auto const& ma_qp_arr = mic_fab_vars[MicVar_Kess::qp]->const_arrays();
        // The sedimentation CFL uses fall speed, not precipitating mass flux.
        // fz stores rho * Vt * qp for the flux divergence below, so the reduction
        // recomputes Vt from the same face state used for the flux.
        GpuTuple<Real> max_terminal_velocity = ParReduce(TypeList<ReduceOpMax>{},
                                                         TypeList<Real>{},
                                                         fz, IntVect(0),
                                                         [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                                         -> GpuTuple<Real>
                                                         {
                                                             const auto& rho_arr = ma_rho_arr[box_no];
                                                             const auto& qp_arr = ma_qp_arr[box_no];
                                                             const Real rho_km1 = (k == k_lo) ? rho_arr(i,j,k) : rho_arr(i,j,k-1);
                                                             const Real rho_k = (k == k_hi+1) ? rho_arr(i,j,k-1) : rho_arr(i,j,k);
                                                             const Real qp_km1 = (k == k_lo) ? qp_arr(i,j,k) : qp_arr(i,j,k-1);
                                                             const Real qp_k = (k == k_hi+1) ? qp_arr(i,j,k-1) : qp_arr(i,j,k);
                                                             const KesslerFaceState face_state =
                                                                 kessler_face_state(k, k_hi, rho_km1, rho_k, qp_km1, qp_k);
                                                             return { kessler_terminal_velocity(face_state.rho, face_state.qp) };
                                                         });
        int n_substep = kessler_num_sedimentation_substeps(get<0>(max_terminal_velocity),
                                                           dt, m_dzmin);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_substep >= 1,
                                         "Kessler: Number of precipitation substeps must be greater than 0!");
        coef /= Real(n_substep);
        dtn  /= Real(n_substep);

        for (int nsub(0); nsub<n_substep; ++nsub) {
            for ( MFIter mfi(*tabs, TilingIfNotGPU()); mfi.isValid(); ++mfi ){
                auto rho_array = mic_fab_vars[MicVar_Kess::rho]->array(mfi);
                auto qp_array  = mic_fab_vars[MicVar_Kess::qp]->array(mfi);
                auto rain_accum_array = mic_fab_vars[MicVar_Kess::rain_accum]->array(mfi);
                auto fz_array  = fz.array(mfi);

                const auto dJ_array = (m_detJ_cc) ? m_detJ_cc->const_array(mfi) : Array4<const Real>{};

                const Box& tbx = mfi.tilebox();
                const Box& tbz = mfi.tilebox(IntVect(0,0,1),IntVect(0));

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    const Real rho_km1 = (k == k_lo) ? rho_array(i,j,k) : rho_array(i,j,k-1);
                    const Real rho_k = (k == k_hi+1) ? rho_array(i,j,k-1) : rho_array(i,j,k);
                    const Real qp_km1 = (k == k_lo) ? qp_array(i,j,k) : qp_array(i,j,k-1);
                    const Real qp_k = (k == k_hi+1) ? qp_array(i,j,k-1) : qp_array(i,j,k);
                    const int donor_k = kessler_face_donor_k(k, k_hi);
                    const KesslerFaceState face_state =
                        kessler_face_state(k, k_hi, rho_km1, rho_k, qp_km1, qp_k);

                    const Real terminal_velocity = kessler_terminal_velocity(face_state.rho, face_state.qp);
                    const Real donor_rho = rho_array(i,j,donor_k);
                    const Real donor_qp = std::max(Real(0), qp_array(i,j,donor_k));
                    const Real donor_detJ = (dJ_array) ? dJ_array(i,j,donor_k) : Real(1);
                    // The face flux uses face rho and donor qp. The limiter uses donor rho,
                    // donor qp, and donor detJ because it caps how much rain can leave the donor
                    // cell in this substep.
                    // Cap outgoing flux by the donor cell's detJ-weighted available rain water:
                    // F * dt / dz <= rho_donor * qp_donor * detJ_donor.
                    // This keeps compressed cells from losing more rain than they contain.
                    const Real max_flux = donor_rho * donor_qp * donor_detJ / coef;
                    fz_array(i,j,k) = amrex::min(
                        kessler_precip_flux(face_state.rho, terminal_velocity, face_state.qp), max_flux);

                    if(k==k_lo){
                        // Surface accumulation stores the bottom-face precipitation mass per area
                        // increment converted to liquid-water depth.
                        rain_accum_array(i,j,k) = rain_accum_array(i,j,k)
                                                + kessler_rain_accumulation_increment(fz_array(i,j,k) * dtn);
                    }
                });

                ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    Real dJinv = (dJ_array) ? Real(1)/dJ_array(i,j,k) : Real(1);

                    // Threshold local face-flux copies only. Neighboring cells share fz
                    // faces, so the cell update must not mutate fz while applying the
                    // small-value cutoff.
                    Real f_hi = fz_array(i,j,k+1);
                    Real f_lo = fz_array(i,j,k  );

                    if (kessler_is_small_sedimentation_value(f_hi)) {
                        f_hi = Real(0);
                    }
                    if (kessler_is_small_sedimentation_value(f_lo)) {
                        f_lo = Real(0);
                    }
                    Real dq_sed = kessler_sedimentation_tendency(
                        f_hi, f_lo, rho_array(i,j,k), dJinv, coef);
                    if (kessler_is_small_sedimentation_value(dq_sed)) {
                        dq_sed = Real(0);
                    }

                    qp_array(i,j,k) +=  dq_sed;
                    qp_array(i,j,k)  = std::max(Real(0), qp_array(i,j,k));
                });
            }
        }
    }

    if (solverChoice.moisture_type == MoistureType::Kessler_NoRain) {
        if (!do_cond) { return; }
        for ( MFIter mfi(*tabs,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto qv_array    = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
            auto qc_array    = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);
            auto qt_array    = mic_fab_vars[MicVar_Kess::qt]->array(mfi);
            auto tabs_array  = mic_fab_vars[MicVar_Kess::tabs]->array(mfi);
            auto theta_array = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
            auto pres_array  = mic_fab_vars[MicVar_Kess::pres]->array(mfi);

            auto tbx = mfi.tilebox();

            Real d_fac_cond = m_fac_cond;

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                qc_array(i,j,k) = std::max(Real(0), qc_array(i,j,k));

                Real qsat, dtqsat;

                Real pressure = pres_array(i,j,k);
                // Kessler stores pressure in mbar / hPa for the qsat helpers.
                erf_qsatw(tabs_array(i,j,k), pressure, qsat);
                erf_dtqsatw(tabs_array(i,j,k), pressure, dtqsat);

                const KesslerSaturationAdjustment saturation_adjustment =
                    kessler_saturation_adjustment(qv_array(i,j,k), qc_array(i,j,k), qsat, dtqsat,
                                                  do_cond, d_fac_cond);

                qv_array(i,j,k) += -saturation_adjustment.dq_vapor_to_cloud
                                 +  saturation_adjustment.dq_cloud_to_vapor;
                qc_array(i,j,k) +=  saturation_adjustment.dq_vapor_to_cloud
                                 -  saturation_adjustment.dq_cloud_to_vapor;

                Real theta_over_T = theta_array(i,j,k)/tabs_array(i,j,k);
                theta_array(i,j,k) += theta_over_T * d_fac_cond
                    * (saturation_adjustment.dq_vapor_to_cloud - saturation_adjustment.dq_cloud_to_vapor);

                qv_array(i,j,k) = std::max(Real(0), qv_array(i,j,k));
                qc_array(i,j,k) = std::max(Real(0), qc_array(i,j,k));

                qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);
            });
        }
    }
}
