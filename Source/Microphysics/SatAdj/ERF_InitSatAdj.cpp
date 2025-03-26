#include "ERF_SatAdj.H"

using namespace amrex;


/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 * @param[in] qc_in Cloud variables input
 * @param[in,out] qv_in Vapor variables input
 * @param[in] qi_in Ice variables input
 * @param[in] grids The boxes on which we will evolve the solution
 * @param[in] geom Geometry associated with these MultiFabs and grids
 * @param[in] dt_advance Timestep for the advance
 */
void SatAdj::Init (const MultiFab& cons_in,
                   const BoxArray& /*grids*/,
                   const Geometry& geom,
                   const Real& dt_advance,
                   std::unique_ptr<MultiFab>& /*z_phys_nd*/,
                   std::unique_ptr<MultiFab>& /*detJ_cc*/)
{
    dt = dt_advance;
    m_geom = geom;

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar_SatAdj::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.);
    }
}

/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 */
void SatAdj::Copy_State_to_Micro (const MultiFab& cons_in)
{
    // Get the temperature, density, theta, qt and qp from input
    for (MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& tbx = mfi.tilebox();

        auto states_array = cons_in.array(mfi);

        auto qv_array    = mic_fab_vars[MicVar_SatAdj::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_SatAdj::qc]->array(mfi);

        auto rho_array   = mic_fab_vars[MicVar_SatAdj::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_SatAdj::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_SatAdj::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_SatAdj::pres]->array(mfi);

        // Get pressure, theta, temperature, density, and qt, qp
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            theta_array(i,j,k) = states_array(i,j,k,RhoTheta_comp)/states_array(i,j,k,Rho_comp);
            qv_array(i,j,k)    = states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp);
            qc_array(i,j,k)    = states_array(i,j,k,RhoQ2_comp)/states_array(i,j,k,Rho_comp);

            tabs_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),
                                                  states_array(i,j,k,RhoTheta_comp),
                                                  qv_array(i,j,k));

            // Pressure in [mbar] for qsat evaluation
            pres_array(i,j,k)  = getPgivenRTh(states_array(i,j,k,RhoTheta_comp), qv_array(i,j,k)) * 0.01;
        });
    }
}

