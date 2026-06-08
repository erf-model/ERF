#include "ERF_WSM6.H"

using namespace amrex;

void
WSM6::Copy_Micro_to_State(MultiFab& cons)
{
    for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();
        auto states = cons.array(mfi);

        auto rho = mic_fab_vars[MicVar_WSM6::rho]->array(mfi);
        auto theta = mic_fab_vars[MicVar_WSM6::theta]->array(mfi);

        auto qv = mic_fab_vars[MicVar_WSM6::qv]->array(mfi);
        auto qc = mic_fab_vars[MicVar_WSM6::qc]->array(mfi);
        auto qi = mic_fab_vars[MicVar_WSM6::qi]->array(mfi);
        auto qr = mic_fab_vars[MicVar_WSM6::qr]->array(mfi);
        auto qs = mic_fab_vars[MicVar_WSM6::qs]->array(mfi);
        auto qg = mic_fab_vars[MicVar_WSM6::qg]->array(mfi);

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            states(i,j,k,RhoTheta_comp) = rho(i,j,k) * theta(i,j,k);
            states(i,j,k,RhoQ1_comp) = rho(i,j,k) * amrex::max(Real(0), qv(i,j,k));
            states(i,j,k,RhoQ2_comp) = rho(i,j,k) * amrex::max(Real(0), qc(i,j,k));
            states(i,j,k,RhoQ3_comp) = rho(i,j,k) * amrex::max(Real(0), qi(i,j,k));
            states(i,j,k,RhoQ4_comp) = rho(i,j,k) * amrex::max(Real(0), qr(i,j,k));
            states(i,j,k,RhoQ5_comp) = rho(i,j,k) * amrex::max(Real(0), qs(i,j,k));
            states(i,j,k,RhoQ6_comp) = rho(i,j,k) * amrex::max(Real(0), qg(i,j,k));
        });
    }

    cons.FillBoundary(m_geom.periodicity());
}
