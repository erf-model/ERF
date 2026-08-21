#include "ERF_WDM6.H"
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <iomanip>

using namespace amrex;

void
WDM6::Copy_Micro_to_State(MultiFab& cons)
{
    // Conservative update of all fields
    for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto& box3d = mfi.tilebox();
        auto states = cons.array(mfi);

        const auto& rho   = mic_fab_vars[MicVar_WDM6::rho]->array(mfi);
        const auto& theta = mic_fab_vars[MicVar_WDM6::theta]->array(mfi);
        const auto& qv    = mic_fab_vars[MicVar_WDM6::qv]->array(mfi);
        const auto& qc    = mic_fab_vars[MicVar_WDM6::qc]->array(mfi);
        const auto& qi    = mic_fab_vars[MicVar_WDM6::qi]->array(mfi);
        const auto& qr    = mic_fab_vars[MicVar_WDM6::qr]->array(mfi);
        const auto& qs    = mic_fab_vars[MicVar_WDM6::qs]->array(mfi);
        const auto& qg    = mic_fab_vars[MicVar_WDM6::qg]->array(mfi);
        const auto& nn    = mic_fab_vars[MicVar_WDM6::nn]->array(mfi);
        const auto& nc    = mic_fab_vars[MicVar_WDM6::nc]->array(mfi);
        const auto& nr    = mic_fab_vars[MicVar_WDM6::nr]->array(mfi);

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            states(i,j,k,RhoTheta_comp) = rho(i,j,k) * theta(i,j,k);

            states(i,j,k,RhoQ1_comp) = rho(i,j,k) * amrex::max(Real(0), qv(i,j,k));
            states(i,j,k,RhoQ2_comp) = rho(i,j,k) * amrex::max(Real(0), qc(i,j,k));
            states(i,j,k,RhoQ3_comp) = rho(i,j,k) * amrex::max(Real(0), qi(i,j,k));
            states(i,j,k,RhoQ4_comp) = rho(i,j,k) * amrex::max(Real(0), qr(i,j,k));
            states(i,j,k,RhoQ5_comp) = rho(i,j,k) * amrex::max(Real(0), qs(i,j,k));
            states(i,j,k,RhoQ6_comp) = rho(i,j,k) * amrex::max(Real(0), qg(i,j,k));

            // Number concentrations
            states(i,j,k,RhoQ7_comp) = rho(i,j,k) * amrex::max(Real(0), nc(i,j,k));
            states(i,j,k,RhoQ8_comp) = rho(i,j,k) * amrex::max(Real(0), nn(i,j,k));
            states(i,j,k,RhoQ9_comp) = rho(i,j,k) * amrex::max(Real(0), nr(i,j,k));
        });
    }

    // No nn bookkeeping needed here. Copy_State_to_Micro decides per cell from
    // the state itself, so once RhoQ8 has been written by the copyback above it
    // is read back like any other variable.

    cons.FillBoundary(m_geom.periodicity());
}
