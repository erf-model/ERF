#include <ERF_EOS.H>
#include <ERF_TileNoZ.H>
#include "ERF_Kessler.H"
#include "ERF_DataStruct.H"

using namespace amrex;

/**
 * Compute Precipitation-related Microphysics quantities.
 */
void Kessler::AdvanceKessler (const SolverChoice &solverChoice)
{
    auto tabs  = mic_fab_vars[MicVar_Kess::tabs];
    if (solverChoice.moisture_type == MoistureType::Kessler){
        auto dz = m_geom.CellSize(2);
        auto domain = m_geom.Domain();
        int k_lo = domain.smallEnd(2);
        int k_hi = domain.bigEnd(2);

        MultiFab fz;
        auto ba    = tabs->boxArray();
        auto dm    = tabs->DistributionMap();
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
                    qp_avg  = qp_array(i,j,k);
                } else if (k==k_hi+1) {
                    rho_avg = rho_array(i,j,k-1);
                    qp_avg  = qp_array(i,j,k-1);
                } else {
                    rho_avg = 0.5*(rho_array(i,j,k-1) + rho_array(i,j,k)); // Convert to g/cm^3
                    qp_avg = 0.5*(qp_array(i,j,k-1)  + qp_array(i,j,k));
                }

                qp_avg = std::max(0.0, qp_avg);

                Real V_terminal = 36.34*std::pow(rho_avg*0.001*qp_avg, 0.1346)*std::pow(rho_avg/1.16, -0.5); // in m/s

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

                /*if(k==0){
                  fz_array(i,j,k) = 0;
                  }*/
            });
        }

        for ( MFIter mfi(*tabs,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
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

                //------- Autoconversion/accretion
                Real qcc, auto_r, accrr;
                Real qsat, dtqsat;
                Real dq_clwater_to_rain, dq_rain_to_vapor, dq_clwater_to_vapor, dq_vapor_to_clwater;

                Real pressure = pres_array(i,j,k);
                erf_qsatw(tabs_array(i,j,k), pressure, qsat);
                erf_dtqsatw(tabs_array(i,j,k), pressure, dtqsat);

                if (qsat <= 0.0) {
                    amrex::Warning("qsat computed as non-positive; setting to 0.!");
                    qsat = 0.0;
                }

                // If there is precipitating water (i.e. rain), and the cell is not saturated
                // then the rain water can evaporate leading to extraction of latent heat, hence
                // reducing temperature and creating negative buoyancy

                dq_clwater_to_rain  = 0.0;
                dq_rain_to_vapor    = 0.0;
                dq_vapor_to_clwater = 0.0;
                dq_clwater_to_vapor = 0.0;

                //Real fac = qsat*4093.0*L_v/(Cp_d*std::pow(tabs_array(i,j,k)-36.0,2));
                //Real fac = qsat*L_v*L_v/(Cp_d*R_v*tabs_array(i,j,k)*tabs_array(i,j,k));
                Real fac = 1.0 + (L_v/Cp_d)*dtqsat;

                // If water vapor content exceeds saturation value, then vapor condenses to water and latent heat is released, increasing temperature
                if (qv_array(i,j,k) > qsat) {
                    dq_vapor_to_clwater = std::min(qv_array(i,j,k), (qv_array(i,j,k)-qsat)/(1.0 + fac));
                }

                // If water vapor is less than the saturated value, then the cloud water can evaporate,
                // leading to evaporative cooling and reducing temperature
                if (qv_array(i,j,k) < qsat && qc_array(i,j,k) > 0.0) {
                    dq_clwater_to_vapor = std::min(qc_array(i,j,k), (qsat - qv_array(i,j,k))/(1.0 + fac));
                }

                if (qp_array(i,j,k) > 0.0 && qv_array(i,j,k) < qsat) {
                    Real C = 1.6 + 124.9*std::pow(0.001*rho_array(i,j,k)*qp_array(i,j,k),0.2046);
                    dq_rain_to_vapor = 1.0/(0.001*rho_array(i,j,k))*(1.0 - qv_array(i,j,k)/qsat)*C*std::pow(0.001*rho_array(i,j,k)*qp_array(i,j,k),0.525)/
                        (5.4e5 + 2.55e6/(pressure*qsat))*dtn;
                    // The negative sign is to make this variable (vapor formed from evaporation)
                    // a positive quantity (as qv/qs < 1)
                    dq_rain_to_vapor = std::min({qp_array(i,j,k), dq_rain_to_vapor});

                    // Removing latent heat due to evaporation from rain water to water vapor, reduces the (potential) temperature
                }

                // If there is cloud water present then do accretion and autoconversion to rain
                if (qc_array(i,j,k) > 0.0) {
                    qcc = qc_array(i,j,k);

                    auto_r = 0.0;
                    if (qcc > qcw0) {
                        auto_r = alphaelq;
                    }

                    accrr = 0.0;
                    accrr = 2.2 * std::pow(qp_array(i,j,k) , 0.875);
                    dq_clwater_to_rain = dtn *(accrr*qcc + auto_r*(qcc - qcw0));

                    // If the amount of change is more than the amount of qc present, then dq = qc
                    dq_clwater_to_rain = std::min(dq_clwater_to_rain, qc_array(i,j,k));
                }

                if(std::fabs(fz_array(i,j,k+1)) < 1e-14) fz_array(i,j,k+1) = 0.0;
                if(std::fabs(fz_array(i,j,k  )) < 1e-14) fz_array(i,j,k  ) = 0.0;
                Real dq_sed = dtn * dJinv * (1.0/rho_array(i,j,k)) * (fz_array(i,j,k+1) - fz_array(i,j,k))/dz;
                if(std::fabs(dq_sed) < 1e-14) dq_sed = 0.0;

                qv_array(i,j,k) += -dq_vapor_to_clwater + dq_clwater_to_vapor + dq_rain_to_vapor;
                qc_array(i,j,k) +=  dq_vapor_to_clwater - dq_clwater_to_vapor - dq_clwater_to_rain;
                qp_array(i,j,k) +=  dq_sed + dq_clwater_to_rain - dq_rain_to_vapor;

                Real theta_over_T = theta_array(i,j,k)/tabs_array(i,j,k);
                theta_array(i,j,k) += theta_over_T * d_fac_cond * (dq_vapor_to_clwater - dq_clwater_to_vapor - dq_rain_to_vapor);

                qv_array(i,j,k) = std::max(0.0, qv_array(i,j,k));
                qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));
                qp_array(i,j,k) = std::max(0.0, qp_array(i,j,k));

                qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);
            });
        }
    }

    if (solverChoice.moisture_type == MoistureType::Kessler_NoRain){

        // get the temperature, dentisy, theta, qt and qc from input
        for ( MFIter mfi(*tabs,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto qv_array    = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
            auto qc_array    = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);
            auto qt_array    = mic_fab_vars[MicVar_Kess::qt]->array(mfi);
            auto tabs_array  = mic_fab_vars[MicVar_Kess::tabs]->array(mfi);
            auto theta_array = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
            auto pres_array  = mic_fab_vars[MicVar_Kess::pres]->array(mfi);

            const auto& box3d = mfi.tilebox();

            // Expose for GPU
            Real d_fac_cond = m_fac_cond;

            ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));

                //------- Autoconversion/accretion
                Real qsat, dtqsat;
                Real dq_clwater_to_vapor, dq_vapor_to_clwater;

                Real pressure = pres_array(i,j,k);
                erf_qsatw(tabs_array(i,j,k), pressure, qsat);
                erf_dtqsatw(tabs_array(i,j,k), pressure, dtqsat);

                // If there is precipitating water (i.e. rain), and the cell is not saturated
                // then the rain water can evaporate leading to extraction of latent heat, hence
                // reducing temperature and creating negative buoyancy

                dq_vapor_to_clwater = 0.0;
                dq_clwater_to_vapor = 0.0;

                //Real fac = qsat*4093.0*L_v/(Cp_d*std::pow(tabs_array(i,j,k)-36.0,2));
                //Real fac = qsat*L_v*L_v/(Cp_d*R_v*tabs_array(i,j,k)*tabs_array(i,j,k));
                Real fac = 1.0 + (L_v/Cp_d)*dtqsat;

                // If water vapor content exceeds saturation value, then vapor condenses to water and latent heat is released, increasing temperature
                if (qv_array(i,j,k) > qsat){
                    dq_vapor_to_clwater = std::min(qv_array(i,j,k), (qv_array(i,j,k)-qsat)/(1.0 + fac));
                }
                // If water vapor is less than the saturated value, then the cloud water can evaporate, leading to evaporative cooling and
                // reducing temperature
                if (qv_array(i,j,k) < qsat && qc_array(i,j,k) > 0.0){
                    dq_clwater_to_vapor = std::min(qc_array(i,j,k), (qsat - qv_array(i,j,k))/(1.0 + fac));
                }

                qv_array(i,j,k) += -dq_vapor_to_clwater + dq_clwater_to_vapor;
                qc_array(i,j,k) +=  dq_vapor_to_clwater - dq_clwater_to_vapor;

                Real theta_over_T = theta_array(i,j,k)/tabs_array(i,j,k);

                theta_array(i,j,k) += theta_over_T * d_fac_cond * (dq_vapor_to_clwater - dq_clwater_to_vapor);

                qv_array(i,j,k) = std::max(0.0, qv_array(i,j,k));
                qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));

                qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);
            });
        }
    }
}
