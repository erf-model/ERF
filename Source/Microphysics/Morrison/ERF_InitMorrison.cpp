#include <AMReX_GpuContainers.H>
#include "ERF_Morrison.H"
#include "ERF_IndexDefines.H"
#include "ERF_PlaneAverage.H"
#include "ERF_EOS.H"
#include "ERF_TileNoZ.H"

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
void
Morrison::Init (const MultiFab& cons_in,
                const BoxArray& /*grids*/,
                const Geometry& geom,
                const Real& dt_advance,
                std::unique_ptr<MultiFab>& z_phys_nd,
                std::unique_ptr<MultiFab>& detJ_cc)
{
    [[maybe_unused]] amrex::Real dt     = dt_advance;
    m_geom = geom;

    m_z_phys_nd = z_phys_nd.get();
    m_detJ_cc   = detJ_cc.get();

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar_Morr::rain_accum, MicVar_Morr::snow_accum, MicVar_Morr::graup_accum};

#if defined(ERF_USE_MORR_FORT) && defined(AMREX_USE_GPU)
    Arena* Arena_Used = The_Managed_Arena();
#else
    Arena* Arena_Used = The_Arena();
#endif

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar_Morr::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect(),
                                                        MFInfo().SetArena(Arena_Used));
        mic_fab_vars[ivar]->setVal(0.);
    }

#ifdef ERF_USE_MORR_FORT
    bool use_cpp;
    amrex::ParmParse pp("erf");
    pp.query("use_morr_cpp_answer", use_cpp);

    if (!use_cpp) {
        MoistureType moisture_type;
        pp.query_enum_case_insensitive("moisture_model",moisture_type);
        int morr_noice     = (moisture_type == MoistureType::Morrison_NoIce);
        int morr_rimed_ice = 0; // This is used to set something called "ihail"
        morr_two_moment_init_c(morr_rimed_ice, morr_noice);
    }
#endif
}


/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 */
void
Morrison::Copy_State_to_Micro (const MultiFab& cons_in)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.growntilebox();

        auto states_array = cons_in.array(mfi);

        // Non-precipitating
        auto qv_array    = mic_fab_vars[MicVar_Morr::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_Morr::qcl]->array(mfi);
        auto qi_array    = mic_fab_vars[MicVar_Morr::qci]->array(mfi);
        auto qn_array    = mic_fab_vars[MicVar_Morr::qn]->array(mfi);
        auto qt_array    = mic_fab_vars[MicVar_Morr::qt]->array(mfi);

        // Precipitating
        auto qpr_array   = mic_fab_vars[MicVar_Morr::qpr]->array(mfi);
        auto qps_array   = mic_fab_vars[MicVar_Morr::qps]->array(mfi);
        auto qpg_array   = mic_fab_vars[MicVar_Morr::qpg]->array(mfi);
        auto qp_array    = mic_fab_vars[MicVar_Morr::qp]->array(mfi);

        auto nc_array    = mic_fab_vars[MicVar_Morr::nc]->array(mfi);
        auto ni_array    = mic_fab_vars[MicVar_Morr::ni]->array(mfi);
        auto nr_array    = mic_fab_vars[MicVar_Morr::nr]->array(mfi);
        auto ns_array    = mic_fab_vars[MicVar_Morr::ns]->array(mfi);
        auto ng_array    = mic_fab_vars[MicVar_Morr::ng]->array(mfi);

        auto rho_array   = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_Morr::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_Morr::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_Morr::pres]->array(mfi);

        // Get pressure, theta, temperature, density, and qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            theta_array(i,j,k) = states_array(i,j,k,RhoTheta_comp)/states_array(i,j,k,Rho_comp);

            qv_array(i,j,k)    = std::max(Real(0),states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp));
            qc_array(i,j,k)    = std::max(Real(0),states_array(i,j,k,RhoQ2_comp)/states_array(i,j,k,Rho_comp));
            qi_array(i,j,k)    = std::max(Real(0),states_array(i,j,k,RhoQ3_comp)/states_array(i,j,k,Rho_comp));
            qn_array(i,j,k)    = qc_array(i,j,k) + qi_array(i,j,k);
            qt_array(i,j,k)    = qv_array(i,j,k) + qn_array(i,j,k);

            qpr_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ4_comp)/states_array(i,j,k,Rho_comp));
            qps_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ5_comp)/states_array(i,j,k,Rho_comp));
            qpg_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ6_comp)/states_array(i,j,k,Rho_comp));

             qp_array(i,j,k)   = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

            nc_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ7_comp) /states_array(i,j,k,Rho_comp));
            ni_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ8_comp) /states_array(i,j,k,Rho_comp));
            nr_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ9_comp) /states_array(i,j,k,Rho_comp));
            ns_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ10_comp)/states_array(i,j,k,Rho_comp));
            ng_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ11_comp)/states_array(i,j,k,Rho_comp));

            tabs_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),
                                                  states_array(i,j,k,RhoTheta_comp),
                                                  qv_array(i,j,k));

            // NOTE: the Morrison Fortran version uses Pa not hPa so we don't divideby 100!
            pres_array(i,j,k)  = getPgivenRTh(states_array(i,j,k,RhoTheta_comp), qv_array(i,j,k)); //  * Real(0.01);
        });
    }
}

