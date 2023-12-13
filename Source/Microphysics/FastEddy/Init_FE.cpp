
#include <AMReX_GpuContainers.H>
#include "Microphysics.H"
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
void FastEddy::Init (const MultiFab& cons_in, MultiFab& qmoist,
                     const BoxArray& grids,
                     const Geometry& geom,
                     const Real& dt_advance)
{
    dt = dt_advance;
    m_geom = geom;
    m_gtoe = grids;

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar_FE::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(), 1, cons_in.nGrowVect());
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
void FastEddy::Copy_State_to_Micro (const MultiFab& cons_in)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto states_array = cons_in.array(mfi);

        auto qv_array     = mic_fab_vars[MicVar_FE::qv]->array(mfi);
        auto qc_array     = mic_fab_vars[MicVar_FE::qc]->array(mfi);
        auto rho_array    = mic_fab_vars[MicVar_FE::rho]->array(mfi);
        auto theta_array  = mic_fab_vars[MicVar_FE::theta]->array(mfi);
        auto temp_array   = mic_fab_vars[MicVar_FE::tabs]->array(mfi);

        // Get pressure, theta, temperature, density, and qt, qp
        amrex::ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            theta_array(i,j,k) = states_array(i,j,k,RhoTheta_comp)/states_array(i,j,k,Rho_comp);
            qv_array(i,j,k)    = states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp);
            qc_array(i,j,k)    = states_array(i,j,k,RhoQ2_comp)/states_array(i,j,k,Rho_comp);
            temp_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),states_array(i,j,k,RhoTheta_comp), qv_array(i,j,k));
        });
    }
}



#if 0
#include <AMReX_GpuContainers.H>
#include "Microphysics.H"
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
void FastEddy::Init (const MultiFab& cons_in, MultiFab& qmoist,
                     const BoxArray& grids,
                     const Geometry& geom,
                     const Real& dt_advance)
{
    m_geom = geom;
    m_gtoe = grids;

    dt = dt_advance;

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar_FE::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(), 1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.);
    }

    // We must initialize to zero since we now need boundary values for the call to getP and we need all values filled
    // The ghost cells of these arrays aren't filled in the boundary condition calls for the state

    qmoist.setVal(0.);

    for ( MFIter mfi(cons_in, TileNoZ()); mfi.isValid(); ++mfi) {

        const auto& box3d = mfi.tilebox();

        const auto& lo = amrex::lbound(box3d);
        const auto& hi = amrex::ubound(box3d);

        nlev = box3d.length(2);
        zlo  = lo.z;
        zhi  = hi.z;
    }

    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in, false); mfi.isValid(); ++mfi) {
        auto states_array = cons_in.array(mfi);

        auto qv_array     = mic_fab_vars[MicVar_FE::qv]->array(mfi);
        auto qc_array     = mic_fab_vars[MicVar_FE::qc]->array(mfi);
        auto rho_array    = mic_fab_vars[MicVar_FE::rho]->array(mfi);
        auto theta_array  = mic_fab_vars[MicVar_FE::theta]->array(mfi);
        auto temp_array   = mic_fab_vars[MicVar_FE::tabs]->array(mfi);

        const auto& box3d = mfi.tilebox();

        // Get pressure, theta, temperature, density, and qt, qp
        amrex::ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            theta_array(i,j,k) = states_array(i,j,k,RhoTheta_comp)/states_array(i,j,k,Rho_comp);
            qv_array(i,j,k)    = states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp);
            qc_array(i,j,k)    = states_array(i,j,k,RhoQ2_comp)/states_array(i,j,k,Rho_comp);
            temp_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),states_array(i,j,k,RhoTheta_comp), qv_array(i,j,k));
        });
    }
}
#endif
