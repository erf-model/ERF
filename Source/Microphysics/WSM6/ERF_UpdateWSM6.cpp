#include "ERF_WSM6.H"
#include "ERF_EOS.H"

using namespace amrex;

void
WSM6::Copy_Micro_to_State(MultiFab& cons)
{
    AMREX_ALWAYS_ASSERT(!m_use_anelastic_reference_pressure || m_rdOcp > Real(0));

    for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();
        auto states = cons.array(mfi);

        auto rho = mic_fab_vars[MicVar_WSM6::rho]->array(mfi);
        auto theta = mic_fab_vars[MicVar_WSM6::theta]->array(mfi);
        auto tabs = mic_fab_vars[MicVar_WSM6::tabs]->array(mfi);
        auto pres = mic_fab_vars[MicVar_WSM6::pres]->array(mfi);
        const bool use_anelastic_reference_pressure = m_use_anelastic_reference_pressure;
        const Real rdOcp = m_rdOcp;

        auto qv = mic_fab_vars[MicVar_WSM6::qv]->array(mfi);
        auto qc = mic_fab_vars[MicVar_WSM6::qc]->array(mfi);
        auto qi = mic_fab_vars[MicVar_WSM6::qi]->array(mfi);
        auto qr = mic_fab_vars[MicVar_WSM6::qr]->array(mfi);
        auto qs = mic_fab_vars[MicVar_WSM6::qs]->array(mfi);
        auto qg = mic_fab_vars[MicVar_WSM6::qg]->array(mfi);

        ParallelFor(box3d, [=]
                    AMREX_GPU_DEVICE(int i, int j, int k) {
            wsm6_copy_micro_to_state_cell(
                states, rho, theta, tabs, pres, qv, qc, qi, qr, qs, qg,
                use_anelastic_reference_pressure, rdOcp, i, j, k);
        });
    }

    cons.FillBoundary(m_geom.periodicity());
}
