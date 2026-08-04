#include <AMReX_GpuContainers.H>
#include "ERF_Kessler.H"
#include "ERF_IndexDefines.H"
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
void Kessler::Init (const MultiFab& cons_in,
                    const BoxArray& /*grids*/,
                    const Geometry& geom,
                    const Real& dt_advance,
                    std::unique_ptr<MultiFab>& z_phys_nd,
                    std::unique_ptr<MultiFab>& detJ_cc)
{
    dt = dt_advance;
    m_geom = geom;

    m_z_phys_nd = z_phys_nd.get();
    m_detJ_cc   = detJ_cc.get();

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar_Kess::rain_accum};

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar_Kess::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.);
    }

    // Set class data members
    for ( MFIter mfi(cons_in, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        const auto& lo = lbound(box3d);
        const auto& hi = ubound(box3d);

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
    Copy_State_to_Micro(cons_in, nullptr);
}

void Kessler::Copy_State_to_Micro (const MultiFab& cons_in,
                                   const MultiFab* base_state)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto states_array = cons_in.array(mfi);

        auto qt_array    = mic_fab_vars[MicVar_Kess::qt]->array(mfi);
        auto qv_array    = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);

        auto qp_array    = mic_fab_vars[MicVar_Kess::qp]->array(mfi);

        auto rho_array   = mic_fab_vars[MicVar_Kess::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_Kess::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_Kess::pres]->array(mfi);
        const auto base_array = base_state ? base_state->const_array(mfi) : Array4<Real const>{};
        const bool use_anelastic_reference_pressure = (base_state != nullptr);
        const Real rdOcp = m_rdOcp;

        // The shared diagnosis returns Pa. Kessler stores pressure in mbar / hPa
        // for the qsat helper path, so convert once at copy-in.
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real rho = states_array(i,j,k,Rho_comp);
            const Real rho_theta = states_array(i,j,k,RhoTheta_comp);
            const Real qv = states_array(i,j,k,RhoQ1_comp) / rho;

            rho_array(i,j,k)   = rho;
            theta_array(i,j,k) = rho_theta / rho;
            qv_array(i,j,k)    = qv;
            qc_array(i,j,k)    = states_array(i,j,k,RhoQ2_comp) / rho;
            qp_array(i,j,k)    = states_array(i,j,k,RhoQ3_comp) / rho;
            qt_array(i,j,k)    = qv_array(i,j,k) + qc_array(i,j,k);

            const Real p0 = use_anelastic_reference_pressure
                ? base_array(i,j,k,BaseState::p0_comp) : Real(0.0);
            const MicrophysicsThermoState thermo = diagnose_microphysics_thermo_state(
                rho, rho_theta, qv, rdOcp, use_anelastic_reference_pressure, p0);
            tabs_array(i,j,k) = thermo.temperature;
            pres_array(i,j,k) = thermo.pressure_pa * Real(0.01);
        });
    }
}

void Kessler::Update_Micro_Vars (MultiFab& cons_in,
                                 const MultiFab* base_state)
{
    Copy_State_to_Micro(cons_in, base_state);
}
