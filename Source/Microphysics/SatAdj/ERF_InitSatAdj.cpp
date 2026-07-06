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
 * Copies conserved ERF state into SatAdj's internal specific variables.
 * rho is copied but not modified by SatAdj.
 * theta, qv, and qc are formed by dividing conserved densities by rho.
 * tabs and pressure are diagnosed from the same conserved state snapshot.
 * pressure is converted from Pa to mbar/hPa for saturation helper calls.
 *
 * @param[in] cons_in Conserved variables input
 */
void SatAdj::Copy_State_to_Micro (const MultiFab& cons_in)
{
    for (MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& tbx = mfi.tilebox();

        auto states_array = cons_in.array(mfi);

        auto qv_array    = mic_fab_vars[MicVar_SatAdj::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_SatAdj::qc]->array(mfi);

        auto rho_array   = mic_fab_vars[MicVar_SatAdj::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_SatAdj::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_SatAdj::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_SatAdj::pres]->array(mfi);

        // Use local scalars so temperature and pressure are diagnosed from the same
        // conserved state values copied into the microphysics arrays. This avoids
        // hidden read-after-write dependencies inside the kernel body.
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real rho      = states_array(i,j,k,Rho_comp);
            const Real rhoTheta = states_array(i,j,k,RhoTheta_comp);
            const Real qv       = states_array(i,j,k,RhoQ1_comp) / rho;
            const Real qc       = states_array(i,j,k,RhoQ2_comp) / rho;
            const Real theta    = rhoTheta / rho;

            rho_array(i,j,k)   = rho;
            theta_array(i,j,k) = theta;
            qv_array(i,j,k)    = qv;
            qc_array(i,j,k)    = qc;
            tabs_array(i,j,k)  = getTgivenRandRTh(rho, rhoTheta, qv);

            // Pressure is stored in mbar/hPa because erf_qsatw expects that unit.
            pres_array(i,j,k)  = getPgivenRTh(rhoTheta, qv) * Real(0.01);
        });
    }
}

