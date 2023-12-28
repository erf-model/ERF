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
void FastEddy::AdvanceFE ()
{
    auto tabs  = mic_fab_vars[MicVar_FE::tabs];

    // get the temperature, dentisy, theta, qt and qc from input
    for ( MFIter mfi(*tabs,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto qv_array    = mic_fab_vars[MicVar_FE::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_FE::qc]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_FE::tabs]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_FE::theta]->array(mfi);
        auto rho_array   = mic_fab_vars[MicVar_FE::rho]->array(mfi);

        const auto& box3d = mfi.tilebox();

        // Expose for GPU
        Real d_fac_cond = m_fac_cond;

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {

            qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));

            //------- Autoconversion/accretion
            Real dq_clwater_to_vapor, dq_vapor_to_clwater, qsat;

            Real pressure = getPgivenRTh(rho_array(i,j,k)*theta_array(i,j,k),qv_array(i,j,k))/100.0;
            erf_qsatw(tabs_array(i,j,k), pressure, qsat);

            // If there is precipitating water (i.e. rain), and the cell is not saturated
            // then the rain water can evaporate leading to extraction of latent heat, hence
            // reducing temperature and creating negative buoyancy

            dq_vapor_to_clwater = 0.0;
            dq_clwater_to_vapor = 0.0;

            //Real fac = qsat*4093.0*L_v/(Cp_d*std::pow(tabs_array(i,j,k)-36.0,2));
            Real fac = qsat*L_v*L_v/(Cp_d*R_v*tabs_array(i,j,k)*tabs_array(i,j,k));

            // If water vapor content exceeds saturation value, then vapor condenses to waterm and latent heat is released, increasing temperature
            if(qv_array(i,j,k) > qsat){
                dq_vapor_to_clwater = std::min(qv_array(i,j,k), (qv_array(i,j,k)-qsat)/(1.0 + fac));
            }
            // If water vapor is less than the satruated value, then the cloud water can evaporate, leading to evaporative cooling and
            // reducing temperature
            if(qv_array(i,j,k) < qsat and qc_array(i,j,k) > 0.0){
                dq_clwater_to_vapor = std::min(qc_array(i,j,k), (qsat - qv_array(i,j,k))/(1.0 + fac));
            }

            qv_array(i,j,k) = qv_array(i,j,k) - dq_vapor_to_clwater + dq_clwater_to_vapor;
            qc_array(i,j,k) = qc_array(i,j,k) + dq_vapor_to_clwater - dq_clwater_to_vapor;

            theta_array(i,j,k) = theta_array(i,j,k) + theta_array(i,j,k)/tabs_array(i,j,k)*d_fac_cond*(dq_vapor_to_clwater - dq_clwater_to_vapor);

            qv_array(i,j,k) = std::max(0.0, qv_array(i,j,k));
            qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));

        });
    }
}
