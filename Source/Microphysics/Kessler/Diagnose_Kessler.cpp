#include "Microphysics.H"
#include "EOS.H"

/**
 * Computes diagnostic quantities like cloud ice/liquid and precipitation ice/liquid
 * from the existing Microphysics variables.
 */
void Kessler::Diagnose ()
{
    auto qt   = mic_fab_vars[MicVar_Kess::qt];
    auto qp   = mic_fab_vars[MicVar_Kess::qp];
    auto qv   = mic_fab_vars[MicVar_Kess::qv];
    auto qn   = mic_fab_vars[MicVar_Kess::qn];
    auto qcl  = mic_fab_vars[MicVar_Kess::qcl];
    auto qci  = mic_fab_vars[MicVar_Kess::qci];
    auto qpl  = mic_fab_vars[MicVar_Kess::qpl];
    auto qpi  = mic_fab_vars[MicVar_Kess::qpi];
    auto tabs = mic_fab_vars[MicVar_Kess::tabs];

    for ( amrex::MFIter mfi(*tabs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto qt_array    = qt->array(mfi);
        auto qp_array    = qp->array(mfi);
        auto qv_array    = qv->array(mfi);
        auto qn_array    = qn->array(mfi);
        auto qcl_array   = qcl->array(mfi);
        auto qci_array   = qci->array(mfi);
        auto qpl_array   = qpl->array(mfi);
        auto qpi_array   = qpi->array(mfi);

        const auto& box3d = mfi.tilebox();

        amrex::ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            qt_array(i,j,k)  = qv_array(i,j,k) + qn_array(i,j,k);
            amrex::Real omn  = 1.0;
            qcl_array(i,j,k) = qn_array(i,j,k)*omn;
            qci_array(i,j,k) = qn_array(i,j,k)*(1.0-omn);
            amrex::Real omp  = 1.0;
            qpl_array(i,j,k) = qp_array(i,j,k)*omp;
            qpi_array(i,j,k) = qp_array(i,j,k)*(1.0-omp);
        });
    }
}

