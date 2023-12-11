#include "Microphysics.H"

/**
 * Computes diagnostic quantities like cloud ice/liquid and precipitation ice/liquid
 * from the existing Microphysics variables.
 */
void FastEddy::Diagnose () {

  auto qt   = mic_fab_vars[MicVar_FE::qt];
  auto qp   = mic_fab_vars[MicVar_FE::qp];
  auto qv   = mic_fab_vars[MicVar_FE::qv];
  auto qn   = mic_fab_vars[MicVar_FE::qn];
  auto qcl  = mic_fab_vars[MicVar_FE::qcl];
  auto qci  = mic_fab_vars[MicVar_FE::qci];
  auto qpl  = mic_fab_vars[MicVar_FE::qpl];
  auto qpi  = mic_fab_vars[MicVar_FE::qpi];
  auto tabs = mic_fab_vars[MicVar_FE::tabs];

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

     amrex::ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
       qv_array(i,j,k)  = std::max(0.0,qt_array(i,j,k) - qn_array(i,j,k));
       amrex::Real omn  = 1.0;
       qcl_array(i,j,k) = qn_array(i,j,k)*omn;
       qci_array(i,j,k) = qn_array(i,j,k)*(1.0-omn);
       amrex::Real omp  = 1.0;;
       qpl_array(i,j,k) = qp_array(i,j,k)*omp;
       qpi_array(i,j,k) = qp_array(i,j,k)*(1.0-omp);
     });
  }
}
