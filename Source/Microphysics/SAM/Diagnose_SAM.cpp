#include "SAM.H"
#include "EOS.H"

/**
 * Computes diagnostic quantities like cloud ice/liquid and precipitation ice/liquid
 * from the existing Microphysics variables.
 */
void SAM::Diagnose () {
    /*
    for ( amrex::MFIter mfi(*tabs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto qt_array    = mic_fab_vars[MicVar::qt]->array(mfi);
        auto qv_array    = mic_fab_vars[MicVar::qv]->array(mfi);
        auto qn_array    = mic_fab_vars[MicVar::qn]->array(mfi);
        auto qcl_array   = mic_fab_vars[MicVar::qcl]->array(mfi);
        auto qci_array   = mic_fab_vars[MicVar::qci]->array(mfi);

        auto qp_array    = mic_fab_vars[MicVar::qp]->array(mfi);
        auto qpr_array   = mic_fab_vars[MicVar::qpr]->array(mfi);
        auto qps_array   = mic_fab_vars[MicVar::qps]->array(mfi);
        auto qpg_array   = mic_fab_vars[MicVar::qpg]->array(mfi);

        auto tabs_array  = tabs->array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::Real omn  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tbgmin)*a_bg));
            amrex::Real omp  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
            amrex::Real omg  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tgrmin)*a_gr));

            // NOTE: Cloud_SAM.cpp iterates for phase change. The result
            //       of those iterations is stored in qn. Therfore,
            qcl_array(i,j,k) = qn_array(i,j,k)*omn;
            qci_array(i,j,k) = qn_array(i,j,k)*(1.0-omn);

            qpr_array(i,j,k) = qp_array(i,j,k)*omp;
            qps_array(i,j,k) = qp_array(i,j,k)*(1.0-omp)*(1.0-omg);
            qpg_array(i,j,k) = qp_array(i,j,k)*(1.0-omp)*omg;
        });
    }
    */
}

