#include "Microphysics.H"

/**
 * Computes diagnostic quantities like cloud ice/liquid and precipitation ice/liquid
 * from the existing Microphysics variables.
 */
void Microphysics::Diagnose () {

  auto qt   = mic_fab_vars[MicVar::qt];
  auto qp   = mic_fab_vars[MicVar::qp];
  auto qv   = mic_fab_vars[MicVar::qv];
  auto qn   = mic_fab_vars[MicVar::qn];
  auto qcl  = mic_fab_vars[MicVar::qcl];
  auto qci  = mic_fab_vars[MicVar::qci];
  auto qpl  = mic_fab_vars[MicVar::qpl];
  auto qpi  = mic_fab_vars[MicVar::qpi];
  auto tabs = mic_fab_vars[MicVar::tabs];

  for ( amrex::MFIter mfi(*tabs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     auto qt_array    = qt->array(mfi);
     auto qp_array    = qp->array(mfi);
     auto qv_array    = qv->array(mfi);
     auto qn_array    = qn->array(mfi);
     auto qcl_array   = qcl->array(mfi);
     auto qci_array   = qci->array(mfi);
     auto qpl_array   = qpl->array(mfi);
     auto qpi_array   = qpi->array(mfi);
     auto tabs_array  = tabs->array(mfi);

     const auto& box3d = mfi.tilebox();

     amrex::ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
       qv_array(i,j,k)  = qt_array(i,j,k) - qn_array(i,j,k);
       amrex::Real omn  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tbgmin)*a_bg));
       qcl_array(i,j,k) = qn_array(i,j,k)*omn;
       qci_array(i,j,k) = qn_array(i,j,k)*(1.0-omn);
       amrex::Real omp  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
       qpl_array(i,j,k) = qp_array(i,j,k)*omp;
       qpi_array(i,j,k) = qp_array(i,j,k)*(1.0-omp);
     });
  }
}

/**
 * Wrapper for the PrecipFall, Cloud, Precipitation, and Diagnostics routines.
 */
void Microphysics::Proc () {

    MicroPrecipFall();

    if (docloud) {
        Cloud();
        if (doprecip) {
            Precip();
        }
    }

    Diagnose();
}
