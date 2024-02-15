#include <AMReX_GpuContainers.H>
#include "Kessler.H"
#include "IndexDefines.H"
#include "PlaneAverage.H"
#include "EOS.H"
#include "TileNoZ.H"

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
void Kessler::Init (const MultiFab& cons_in,
                    const BoxArray& grids,
                    const Geometry& geom,
                    const Real& dt_advance)
{
    dt = dt_advance;
    m_geom = geom;
    m_gtoe = grids;

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar_Kess::qt, MicVar_Kess::qv, MicVar_Kess::qcl, MicVar_Kess::qp};

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar_Kess::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.);
    }

    // Set class data members
    for ( MFIter mfi(cons_in, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        const auto& lo = amrex::lbound(box3d);
        const auto& hi = amrex::ubound(box3d);

        nlev = box3d.length(2);
        zlo  = lo.z;
        zhi  = hi.z;
    }
}

/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 */
void Kessler::Copy_State_to_Micro (const MultiFab& cons_in)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto states_array = cons_in.array(mfi);

        auto qv_array    = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);
        auto qp_array    = mic_fab_vars[MicVar_Kess::qp]->array(mfi);
        auto qt_array    = mic_fab_vars[MicVar_Kess::qt]->array(mfi);

        auto rho_array   = mic_fab_vars[MicVar_Kess::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_Kess::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_Kess::pres]->array(mfi);

        // Get pressure, theta, temperature, density, and qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            theta_array(i,j,k) = states_array(i,j,k,RhoTheta_comp)/states_array(i,j,k,Rho_comp);
            qv_array(i,j,k)    = states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp);
            qc_array(i,j,k)    = states_array(i,j,k,RhoQ2_comp)/states_array(i,j,k,Rho_comp);
            qp_array(i,j,k)    = states_array(i,j,k,RhoQ3_comp)/states_array(i,j,k,Rho_comp);
            qt_array(i,j,k)    = qv_array(i,j,k) + qc_array(i,j,k);

            tabs_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),
                                                  states_array(i,j,k,RhoTheta_comp),
                                                  qv_array(i,j,k));
            pres_array(i,j,k)  = getPgivenRTh(states_array(i,j,k,RhoTheta_comp), qv_array(i,j,k))/100.;
        });
    }
}

