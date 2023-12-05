#include <ERF.H>
#include <EOS.H>
#include <TileNoZ.H>

using namespace amrex;

/**
 * Function for computing the source terms for cloud vapor, cloud water and heat due to condensation.
 * This routine is only called when using the warm moisture only formulation.
 *
 * @param[out] source the source terms for cloud vapor, cloud water and heat computed here
 * @param[in]  S current solution
 * @param[in]  tau_cond condensation process parameter
 * @param[in]  c_p    Specific heat
 */

#if defined(ERF_USE_WARM_NO_PRECIP)
void
ERF::condensation_source (MultiFab& source, MultiFab& S, Real tau_cond, Real c_p)
{
    for ( MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();
        auto const&   s_arr = S.const_array(mfi);
        auto        src_arr = source.array(mfi);

        // Coefficients and formula from Flatau et al (1992)
        // "Polynomial Fits to Saturation Vapor Pressure", J. Applied Meteorology, 31(12), 1507-1513, 1992.
        // DOI: https://doi.org/10.1175/1520-0450(1992)031<1507:pftsvp>2.0.co;2
        // as used in Munoz-Esparza et al (2022), JAMES
        // "The FastEddy Resident-GPU Acclerated Large-Eddy Simulation Framework:
        //   Moist Dynamics Extension, Validation and Sensitivities of Modeling Non-Precipitating Shallow Cumulus Clouds"
        // DOI: https://doi.org/10.1029/2021MS002904

        Real g0 = -0.29912729e4;
        Real g1 = -0.60170128e4;
        Real g2 = 0.1887643854e2;
        Real g3 = -0.28354721e-1;
        Real g4 =  0.17838301e-4;
        Real g5 = -0.84150417e-9;
        Real g6 =  0.44412543e-12;
        Real g7 =  0.2858487e1;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real      rho_v = s_arr(i,j,k,RhoQv_comp);
            Real      rho_c = s_arr(i,j,k,RhoQc_comp);

            if ( (rho_v + rho_c) > 0.) {
                Real      rho_d = s_arr(i,j,k,Rho_comp);
                Real rhotheta_d = s_arr(i,j,k,RhoTheta_comp);
                Real    theta_d = rhotheta_d / rho_d;

                Real T   = getTgivenRandRTh(rho_d, rhotheta_d);
                Real Tsq = getTgivenRandRTh(rho_d, rhotheta_d);

                Real ln_pvs = ( g0 + (g1 + (g2 + g7 * std::log(T) + (g3 + (g4 + (g5 + g6*T) * T) * T) * T) * T) ) / Tsq;
                Real p_vs   = std::exp(ln_pvs); // saturation vapor pressure computed using 8th-order polynomial fitting
                Real rho_vs = p_vs / (R_d * T);
                Real q_vs   = rho_vs / rho_d;

                Real denom = 1.0 + (L_v * L_v * q_vs)/(c_p * R_v * Tsq);
                Real f_cond_star = (rho_v - rho_vs) * denom;

                Real f_lim  = rho_c;
                Real f_cond = std::max(f_cond_star, -f_lim) / tau_cond*0.0;

                src_arr(i,j,k,RhoQv_comp)    = -f_cond;
                src_arr(i,j,k,RhoQc_comp)    =  f_cond;

                src_arr(i,j,k,RhoTheta_comp) =  theta_d * L_v / (T * c_p) * f_cond;
            } // net water > 0
        });

    } // mfi
}
#endif
