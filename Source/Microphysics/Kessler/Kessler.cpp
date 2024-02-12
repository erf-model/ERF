/*
 * this file is modified from precip_proc from samxx
 */
#include "Microphysics.H"

#include <EOS.H>
#include <TileNoZ.H>

using namespace amrex;

/**
 * Compute Precipitation-related Microphysics quantities.
 */
void Kessler::AdvanceKessler ()
{
    auto qv   = mic_fab_vars[MicVar_Kess::qv];
    auto qc   = mic_fab_vars[MicVar_Kess::qcl];
    auto qp   = mic_fab_vars[MicVar_Kess::qp];
    auto tabs = mic_fab_vars[MicVar_Kess::tabs];
    auto pres = mic_fab_vars[MicVar_Kess::pres];
    auto theta = mic_fab_vars[MicVar_Kess::theta];
    auto rho   = mic_fab_vars[MicVar_Kess::rho];

    auto dz = m_geom.CellSize(2);
    auto domain = m_geom.Domain();
    int k_lo = domain.smallEnd(2);
    int k_hi = domain.bigEnd(2);

    MultiFab fz;
    auto ba    = tabs->boxArray();
    auto dm    = tabs->DistributionMap();
    fz.define(convert(ba, IntVect(0,0,1)), dm, 1, 0); // No ghost cells

    for ( MFIter mfi(fz, TilingIfNotGPU()); mfi.isValid(); ++mfi ){
        auto rho_array = mic_fab_vars[MicVar_Kess::rho]->array(mfi);
        auto qp_array  = mic_fab_vars[MicVar_Kess::qp]->array(mfi);
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

            fz_array(i,j,k) = rho_avg*V_terminal*qp_avg;

            /*if(k==0){
              fz_array(i,j,k) = 0;
              }*/
        });
    }

    Real dtn = dt;

    for ( MFIter mfi(*tabs,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto qv_array   = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
        auto qc_array   = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);
        auto qp_array   = mic_fab_vars[MicVar_Kess::qp]->array(mfi);
        auto qt_array   = mic_fab_vars[MicVar_Kess::qt]->array(mfi);
        auto tabs_array = mic_fab_vars[MicVar_Kess::tabs]->array(mfi);
        auto pres_array = mic_fab_vars[MicVar_Kess::pres]->array(mfi);
        auto theta_array = theta->array(mfi);
        auto rho_array   = mic_fab_vars[MicVar_Kess::rho]->array(mfi);

        const auto& box3d = mfi.tilebox();

        auto fz_array  = fz.array(mfi);

        // Expose for GPU
        Real d_fac_cond = m_fac_cond;

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            qv_array(i,j,k) = std::max(0.0, qv_array(i,j,k));
            qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));
            qp_array(i,j,k) = std::max(0.0, qp_array(i,j,k));

            //------- Autoconversion/accretion
            Real qcc, autor, accrr, dq_clwater_to_rain, dq_rain_to_vapor, dq_clwater_to_vapor, dq_vapor_to_clwater, qsat;

            Real pressure = pres_array(i,j,k);
            erf_qsatw(tabs_array(i,j,k), pressure, qsat);

            // If there is precipitating water (i.e. rain), and the cell is not saturated
            // then the rain water can evaporate leading to extraction of latent heat, hence
            // reducing temperature and creating negative buoyancy

            dq_clwater_to_rain  = 0.0;
            dq_rain_to_vapor    = 0.0;
            dq_vapor_to_clwater = 0.0;
            dq_clwater_to_vapor = 0.0;


            Real fac = qsat*4093.0*L_v/(Cp_d*std::pow(tabs_array(i,j,k)-36.0,2));
            //Real fac = qsat*L_v*L_v/(Cp_d*R_v*tabs_array(i,j,k)*tabs_array(i,j,k));

            // If water vapor content exceeds saturation value, then vapor condenses to waterm and latent heat is released, increasing temperature
            if(qv_array(i,j,k) > qsat){
                dq_vapor_to_clwater = std::min(qv_array(i,j,k), (qv_array(i,j,k)-qsat)/(1.0 + fac));
            }


            // If water vapor is less than the satruated value, then the cloud water can evaporate, leading to evaporative cooling and
            // reducing temperature
            if(qv_array(i,j,k) < qsat and qc_array(i,j,k) > 0.0){
                dq_clwater_to_vapor = std::min(qc_array(i,j,k), (qsat - qv_array(i,j,k))/(1.0 + fac));
            }

            if(qp_array(i,j,k) > 0.0 && qv_array(i,j,k) < qsat) {
                Real C = 1.6 + 124.9*std::pow(0.001*rho_array(i,j,k)*qp_array(i,j,k),0.2046);
                dq_rain_to_vapor = 1.0/(0.001*rho_array(i,j,k))*(1.0 - qv_array(i,j,k)/qsat)*C*std::pow(0.001*rho_array(i,j,k)*qp_array(i,j,k),0.525)/
                    (5.4e5 + 2.55e6/(pressure*qsat))*dtn;
                // The negative sign is to make this variable (vapor formed from evaporation)
                // a poistive quantity (as qv/qs < 1)
                dq_rain_to_vapor = std::min({qp_array(i,j,k), dq_rain_to_vapor});

                // Removing latent heat due to evaporation from rain water to water vapor, reduces the (potential) temperature
            }

            // If there is cloud water present then do accretion and autoconversion to rain

            if (qc_array(i,j,k) > 0.0) {
                qcc = qc_array(i,j,k);

                autor = 0.0;
                if (qcc > qcw0) {
                    autor = 0.001;
                }

                accrr = 0.0;
                accrr = 2.2 * std::pow(qp_array(i,j,k) , 0.875);
                dq_clwater_to_rain = dtn *(accrr*qcc + autor*(qcc - 0.001));

                // If the amount of change is more than the amount of qc present, then dq = qc
                dq_clwater_to_rain = std::min(dq_clwater_to_rain, qc_array(i,j,k));
            }

            if(std::fabs(fz_array(i,j,k+1)) < 1e-14) fz_array(i,j,k+1) = 0.0;
            if(std::fabs(fz_array(i,j,k)) < 1e-14) fz_array(i,j,k) = 0.0;
            Real dq_sed = 1.0/rho_array(i,j,k)*(fz_array(i,j,k+1) - fz_array(i,j,k))/dz*dtn;
            if(std::fabs(dq_sed) < 1e-14)dq_sed = 0.0;
            //dq_sed = 0.0;

            qv_array(i,j,k) = qv_array(i,j,k) - dq_vapor_to_clwater + dq_clwater_to_vapor + dq_rain_to_vapor;
            qc_array(i,j,k) = qc_array(i,j,k) + dq_vapor_to_clwater - dq_clwater_to_vapor - dq_clwater_to_rain;
            qp_array(i,j,k) = qp_array(i,j,k) + dq_sed + dq_clwater_to_rain - dq_rain_to_vapor;

            theta_array(i,j,k) = theta_array(i,j,k) + theta_array(i,j,k)/tabs_array(i,j,k)*d_fac_cond*(dq_vapor_to_clwater - dq_clwater_to_vapor - dq_rain_to_vapor);

            qv_array(i,j,k) = std::max(0.0, qv_array(i,j,k));
            qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));
            qp_array(i,j,k) = std::max(0.0, qp_array(i,j,k));

            qt_array(i,j,k) = qv_array(i,j,k) + qc_array(i,j,k);
        });
    }
}
