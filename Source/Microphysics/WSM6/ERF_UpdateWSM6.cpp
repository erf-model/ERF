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
            theta(i,j,k) = use_anelastic_reference_pressure
                ? wsm6_theta_from_temperature_pressure(tabs(i,j,k), pres(i,j,k), rdOcp)
                : getThgivenRandT(rho(i,j,k), tabs(i,j,k), rdOcp, qv(i,j,k));
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
