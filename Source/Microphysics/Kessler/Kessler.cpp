#include <EOS.H>
#include <TileNoZ.H>
#include "Kessler.H"
#include "DataStruct.H"

using namespace amrex;

/**
 * Compute Precipitation-related Microphysics quantities.
 */
void Kessler::AdvanceKessler (const SolverChoice &solverChoice)
{
    if (solverChoice.moisture_type == MoistureType::Kessler){
        auto dz = m_geom.CellSize(2);
        auto domain = m_geom.Domain();
        int k_lo = domain.smallEnd(2);
        int k_hi = domain.bigEnd(2);

        MultiFab fz;
        auto ba    = mic_fab_vars[MicVar_Kess::tabs]->boxArray();
        auto dm    = mic_fab_vars[MicVar_Kess::tabs]->DistributionMap();
        fz.define(convert(ba, IntVect(0,0,1)), dm, 1, 0); // No ghost cells

        Real dtn = dt;

        for ( MFIter mfi(fz, TilingIfNotGPU()); mfi.isValid(); ++mfi ){
            auto rho_array = mic_fab_vars[MicVar_Kess::rho]->array(mfi);
            auto qp_array  = mic_fab_vars[MicVar_Kess::qp]->array(mfi);
            auto rain_accum_array = mic_fab_vars[MicVar_Kess::rain_accum]->array(mfi);

            auto fz_array  = fz.array(mfi);
            const Box& tbz = mfi.tilebox();

            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real rho_avg, qp_avg;

                if (k==k_lo) {
                    rho_avg = rho_array(i,j,k);
                     qp_avg = qp_array(i,j,k);
                } else if (k==k_hi+1) {
                    rho_avg = rho_array(i,j,k-1);
                     qp_avg = qp_array(i,j,k-1);
                } else {
                    rho_avg = 0.5*(rho_array(i,j,k-1) + rho_array(i,j,k)); // Convert to g/cm^3
                     qp_avg = 0.5*(qp_array(i,j,k-1)  + qp_array(i,j,k));
                }

                qp_avg = std::max(0.0, qp_avg);

                // NOTE: rho is in units of g/cm^3 for the Kessler equations but we transport kg/m^3
                Real rho_cgs = rho_avg * 0.001;

                // TODO: This should come from rho_hse at k=0!
                Real rho_hse = 1.16;

                // NOTE: Converted from cm/s -> m/s
                Real V_terminal = 36.34*std::pow(rho_cgs*qp_avg, 0.1346)*std::pow(rho_avg/rho_hse, -0.5); // in m/s

                // NOTE: Fz is the sedimentation flux from the advective operator.
                //       In the terrain-following coordinate system, the z-deriv in
                //       the divergence uses the normal velocity (Omega). However,
                //       there are no u/v components to the sedimentation velocity.
                //       Therefore, we simply end up with a division by detJ when
                //       evaluating the source term: dJinv * (flux_hi - flux_lo) * dzinv.
                fz_array(i,j,k) = rho_avg*V_terminal*qp_avg;

                if(k==k_lo){
                    rain_accum_array(i,j,k) = rain_accum_array(i,j,k) + rho_avg*qp_avg*V_terminal*dtn/1000.0*1000.0; // Divide by rho_water and convert to mm
                }
            });
        }


        for ( MFIter mfi(*mic_fab_vars[MicVar_Kess::tabs],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto qv_array    = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
            auto qc_array    = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);
            auto qp_array    = mic_fab_vars[MicVar_Kess::qp]->array(mfi);
            auto qt_array    = mic_fab_vars[MicVar_Kess::qt]->array(mfi);
            auto tabs_array  = mic_fab_vars[MicVar_Kess::tabs]->array(mfi);
            auto pres_array  = mic_fab_vars[MicVar_Kess::pres]->array(mfi);
            auto theta_array = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
            auto rho_array   = mic_fab_vars[MicVar_Kess::rho]->array(mfi);

            const auto dJ_array = (m_detJ_cc) ? m_detJ_cc->const_array(mfi) : Array4<const Real>{};

            const auto& box3d = mfi.tilebox();

            auto fz_array  = fz.array(mfi);

            // Expose for GPU
            Real d_fac_cond = m_fac_cond;

            ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Jacobian determinant
                Real dJinv = (dJ_array) ? 1.0/dJ_array(i,j,k) : 1.0;

                qv_array(i,j,k) = std::max(0.0, qv_array(i,j,k));
                qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));
                qp_array(i,j,k) = std::max(0.0, qp_array(i,j,k));
                qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);

                //------- Autoconversion/accretion
                Real qsat, dqsat;
                Real qcc, autor, accrr;
                Real delta_qv, delta_qc;
                Real dq_clwater_to_rain, dq_rain_to_vapor;
                Real dq_clwater_to_vapor, dq_vapor_to_clwater;

                // Pressure alread in mbar
                Real tabs = tabs_array(i,j,k);
                Real pres = pres_array(i,j,k);
                erf_qsatw(tabs, pres, qsat);
                erf_dtqsatw(tabs, pres, dqsat);

                // Rho cgs units
                Real rho_cgs = rho_array(i,j,k) * 0.001;

                // Zero source terms to begin
                dq_clwater_to_rain  = 0.0;
                dq_rain_to_vapor    = 0.0;
                dq_vapor_to_clwater = 0.0;
                dq_clwater_to_vapor = 0.0;

                // Factor for cloud analytical solution (\partial_{T} qvs * Lv/Cp)
                Real fac = dqsat * d_fac_cond;

                // Vapor condenses to water (exothermic)
                if(qv_array(i,j,k) > qsat){
                    dq_vapor_to_clwater = std::min(qv_array(i,j,k), (qv_array(i,j,k)-qsat)/(1.0 + fac));
                }

                // Water evaporated to vapor (endothermic)
                if(qv_array(i,j,k) < qsat && qc_array(i,j,k) > 0.0){
                    dq_clwater_to_vapor = std::min(qc_array(i,j,k), (qsat - qv_array(i,j,k))/(1.0 + fac));
                }

                // Evaporation of rain to vapor
                if(qp_array(i,j,k) > 0.0 && qv_array(i,j,k) < qsat) {
                    Real C = 1.6 + 124.9*std::pow(rho_cgs*qp_array(i,j,k),0.2046);
                    dq_rain_to_vapor = (1.0/rho_cgs)*(1.0 - qv_array(i,j,k)/qsat)*C*std::pow(rho_cgs*qp_array(i,j,k),0.525) /
                                       (5.4e5 + 2.55e6/(pres*qsat))*dtn;
                    dq_rain_to_vapor = std::min(qp_array(i,j,k), dq_rain_to_vapor);
                }

                // Accretion and autoconversion (cloud to rain)
                if (qc_array(i,j,k) > 0.0) {
                    qcc = qc_array(i,j,k);

                    autor = 0.0;
                    if (qcc > qcw0) {
                        autor = 0.001;
                    }

                    accrr = 2.2 * std::pow(qp_array(i,j,k) , 0.875);
                    dq_clwater_to_rain = dtn *(accrr*qcc + autor*(qcc - qcw0));
                    dq_clwater_to_rain = std::min(qc_array(i,j,k), dq_clwater_to_rain);
                }

                if(std::fabs(fz_array(i,j,k+1)) < 1e-14) fz_array(i,j,k+1) = 0.0;
                if(std::fabs(fz_array(i,j,k  )) < 1e-14) fz_array(i,j,k  ) = 0.0;
                Real dq_sed = dtn * dJinv * (1.0/rho_array(i,j,k)) * (fz_array(i,j,k+1) - fz_array(i,j,k))/dz;
                if(std::fabs(dq_sed) < 1e-14) dq_sed = 0.0;

                qv_array(i,j,k) = qv_array(i,j,k) - dq_vapor_to_clwater + dq_clwater_to_vapor + dq_rain_to_vapor;
                qc_array(i,j,k) = qc_array(i,j,k) + dq_vapor_to_clwater - dq_clwater_to_vapor - dq_clwater_to_rain;
                qp_array(i,j,k) = qp_array(i,j,k) + dq_sed + dq_clwater_to_rain - dq_rain_to_vapor;
                qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);

                theta_array(i,j,k) = theta_array(i,j,k) +
                                     theta_array(i,j,k)/tabs_array(i,j,k)*d_fac_cond*
                                     (dq_vapor_to_clwater - dq_clwater_to_vapor - dq_rain_to_vapor);
                tabs_array(i,j,k)  = getTgivenRandRTh(rho_array(i,j,k),
                                                      rho_array(i,j,k)*theta_array(i,j,k),
                                                      qv_array(i,j,k));
                pres_array(i,j,k)  = getPgivenRTh(rho_array(i,j,k)*theta_array(i,j,k), qv_array(i,j,k)) * 0.01;
            });
        }
    }


    if(solverChoice.moisture_type == MoistureType::Kessler_NoRain){

        // get the temperature, dentisy, theta, qt and qc from input
        for ( MFIter mfi(*mic_fab_vars[MicVar_Kess::tabs],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto qv_array    = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
            auto qc_array    = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);
            auto qt_array    = mic_fab_vars[MicVar_Kess::qt]->array(mfi);
            auto tabs_array  = mic_fab_vars[MicVar_Kess::tabs]->array(mfi);
            auto theta_array = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
            auto pres_array  = mic_fab_vars[MicVar_Kess::pres]->array(mfi);
            auto rho_array   = mic_fab_vars[MicVar_Kess::rho]->array(mfi);

            const auto& box3d = mfi.tilebox();

            // Expose for GPU
            Real d_fac_cond = m_fac_cond;

            ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                qv_array(i,j,k) = std::max(0.0, qv_array(i,j,k));
                qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));
                qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);

                //------- Autoconversion/accretion
                Real qsat, dqsat;
                Real dq_clwater_to_vapor, dq_vapor_to_clwater;
                Real delta_qv, delta_qc;

                // Pressure alread in mbar
                Real tabs = tabs_array(i,j,k);
                Real pres = pres_array(i,j,k);
                erf_qsatw(tabs, pres, qsat);
                erf_dtqsatw(tabs, pres, dqsat);

                // Zero source terms to begin
                dq_vapor_to_clwater = 0.0;
                dq_clwater_to_vapor = 0.0;

                // Factor for cloud analytical solution (\partial_{T} qvs * Lv/Cp)
                Real fac = dqsat * d_fac_cond;

                // Vapor condenses to water (exothermic)
                if(qv_array(i,j,k) > qsat){
                    dq_vapor_to_clwater = std::min(qv_array(i,j,k), (qv_array(i,j,k)-qsat)/(1.0 + fac));
                }

                // Water evaporated to vapor (endothermic)
                if(qv_array(i,j,k) < qsat && qc_array(i,j,k) > 0.0){
                    dq_clwater_to_vapor = std::min(qc_array(i,j,k), (qsat - qv_array(i,j,k))/(1.0 + fac));
                }

                qv_array(i,j,k) = qv_array(i,j,k) - dq_vapor_to_clwater + dq_clwater_to_vapor;
                qc_array(i,j,k) = qc_array(i,j,k) + dq_vapor_to_clwater - dq_clwater_to_vapor;
                qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);

                theta_array(i,j,k) = theta_array(i,j,k) +
                                     theta_array(i,j,k)/tabs_array(i,j,k)*d_fac_cond*
                                     (dq_vapor_to_clwater - dq_clwater_to_vapor);
                tabs_array(i,j,k)  = getTgivenRandRTh(rho_array(i,j,k),
                                                      rho_array(i,j,k)*theta_array(i,j,k),
                                                      qv_array(i,j,k));
                pres_array(i,j,k)  = getPgivenRTh(rho_array(i,j,k)*theta_array(i,j,k), qv_array(i,j,k)) * 0.01;
            });
        }
    }
}
