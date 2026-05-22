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
    int i_lo = domain.smallEnd(0);
    int i_hi = domain.bigEnd(0);
    int j_lo = domain.smallEnd(1);
    int j_hi = domain.bigEnd(1);
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
            if (tbx.smallEnd(0) == i_lo) { tbx.growLo(0,-m_real_width); }
            if (tbx.bigEnd(0)   == i_hi) { tbx.growHi(0,-m_real_width); }
            if (tbx.smallEnd(1) == j_lo) { tbx.growLo(1,-m_real_width); }
            if (tbx.bigEnd(1)   == j_hi) { tbx.growHi(1,-m_real_width); }

            Real d_fac_cond = m_fac_cond;

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                qv_array(i,j,k) = std::max(Real(0), qv_array(i,j,k));
                qc_array(i,j,k) = std::max(Real(0), qc_array(i,j,k));
                qp_array(i,j,k) = std::max(Real(0), qp_array(i,j,k));

                Real qsat_local, dtqsat_local;
                Real pressure = pres_array(i,j,k);
                erf_qsatw(tabs_array(i,j,k), pressure, qsat_local);
                erf_dtqsatw(tabs_array(i,j,k), pressure, dtqsat_local);

                if (qsat_local <= Real(0)) {
                    amrex::Warning("qsat computed as non-positive; setting to Real(0)!");
                    qsat_local = Real(0);
                }

                const KesslerSourceTerms source_terms = kessler_warm_rain_sources(
                    qv_array(i,j,k), qc_array(i,j,k), qp_array(i,j,k), rho_array(i,j,k),
                    pressure, qsat_local, dtqsat_local, dtn, do_cond);

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
                    kessler_face_state(k, k_lo, k_hi, rho_km1, rho_k, qp_km1, qp_k);

                const Real terminal_velocity = kessler_terminal_velocity(face_state.rho, face_state.qp);
                fz_array(i,j,k) = kessler_precip_flux(face_state.rho, terminal_velocity, face_state.qp);
            });
        }

        auto const& ma_fz_arr = fz.const_arrays();
        GpuTuple<Real> max = ParReduce(TypeList<ReduceOpMax>{},
                                       TypeList<Real>{},
                                       fz, IntVect(0),
                                       [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                       -> GpuTuple<Real>
                                       {
                                           return { ma_fz_arr[box_no](i,j,k) };
                                       });
        int n_substep = kessler_num_sedimentation_substeps(get<0>(max), dt, m_dzmin);
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
                    const KesslerFaceState face_state =
                        kessler_face_state(k, k_lo, k_hi, rho_km1, rho_k, qp_km1, qp_k);

                    const Real terminal_velocity = kessler_terminal_velocity(face_state.rho, face_state.qp);
                    fz_array(i,j,k) = kessler_precip_flux(face_state.rho, terminal_velocity, face_state.qp);

                    if(k==k_lo){
                        rain_accum_array(i,j,k) = rain_accum_array(i,j,k)
                                                + face_state.rho * face_state.qp * terminal_velocity * dtn / Real(1000.0) * Real(1000.0);
                    }
                });

                ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    Real dJinv = (dJ_array) ? Real(1)/dJ_array(i,j,k) : Real(1);

                    if (kessler_is_small_sedimentation_value(fz_array(i,j,k+1))) {
                        fz_array(i,j,k+1) = Real(0);
                    }
                    if (kessler_is_small_sedimentation_value(fz_array(i,j,k  ))) {
                        fz_array(i,j,k  ) = Real(0);
                    }
                    Real dq_sed = kessler_sedimentation_tendency(
                        fz_array(i,j,k+1), fz_array(i,j,k), rho_array(i,j,k), dJinv, coef);
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
            if (tbx.smallEnd(0) == i_lo) { tbx.growLo(0,-m_real_width); }
            if (tbx.bigEnd(0)   == i_hi) { tbx.growHi(0,-m_real_width); }
            if (tbx.smallEnd(1) == j_lo) { tbx.growLo(1,-m_real_width); }
            if (tbx.bigEnd(1)   == j_hi) { tbx.growHi(1,-m_real_width); }

            Real d_fac_cond = m_fac_cond;

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                qc_array(i,j,k) = std::max(Real(0), qc_array(i,j,k));

                Real qsat, dtqsat;

                Real pressure = pres_array(i,j,k);
                erf_qsatw(tabs_array(i,j,k), pressure, qsat);
                erf_dtqsatw(tabs_array(i,j,k), pressure, dtqsat);

                const KesslerSaturationAdjustment saturation_adjustment =
                    kessler_saturation_adjustment(qv_array(i,j,k), qc_array(i,j,k), qsat, dtqsat, do_cond);

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
