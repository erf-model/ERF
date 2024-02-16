#include "SAM.H"

void SAM::Diagnose ()
{
    auto qt   = mic_fab_vars[MicVar::qt];
    auto qp   = mic_fab_vars[MicVar::qp];
    auto qv   = mic_fab_vars[MicVar::qv];
    auto qn   = mic_fab_vars[MicVar::qn];
    auto qcl  = mic_fab_vars[MicVar::qcl];
    auto qci  = mic_fab_vars[MicVar::qci];
    auto qpr  = mic_fab_vars[MicVar::qpr];
    auto qps  = mic_fab_vars[MicVar::qps];
    auto qpg  = mic_fab_vars[MicVar::qpg];
    auto tabs = mic_fab_vars[MicVar::tabs];

    for ( amrex::MFIter mfi(*tabs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto qt_array    = qt->array(mfi);
        auto qp_array    = qp->array(mfi);
        auto qv_array    = qv->array(mfi);
        auto qn_array    = qn->array(mfi);
        auto qcl_array   = qcl->array(mfi);
        auto qci_array   = qci->array(mfi);
        auto qpr_array   = qpr->array(mfi);
        auto qps_array   = qps->array(mfi);
        auto qpg_array   = qpg->array(mfi);
        auto tabs_array  = tabs->array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::Real omn  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tbgmin)*a_bg));
            amrex::Real omp  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
            amrex::Real omg  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tgrmin)*a_gr));

            qt_array(i,j,k)  = qv_array(i,j,k) + qn_array(i,j,k);
            qcl_array(i,j,k) = qn_array(i,j,k)*omn;
            qci_array(i,j,k) = qn_array(i,j,k)*(1.0-omn);

            qpr_array(i,j,k) = qp_array(i,j,k)*omp;
            qps_array(i,j,k) = qp_array(i,j,k)*(1.0-omp)*(1.0-omg);
            qpg_array(i,j,k) = qp_array(i,j,k)*(1.0-omp)*omg;
        });
    }
}

