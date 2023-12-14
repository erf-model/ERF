#include "Microphysics.H"
#include "EOS.H"

/**
 * Computes diagnostic quantities like cloud ice/liquid and precipitation ice/liquid
 * from the existing Microphysics variables.
 */
void SAM::Diagnose () {

  auto qt   = mic_fab_vars[MicVar::qt];
  auto qp   = mic_fab_vars[MicVar::qp];
  auto qv   = mic_fab_vars[MicVar::qv];
  auto qn   = mic_fab_vars[MicVar::qn];
  auto qcl  = mic_fab_vars[MicVar::qcl];
  auto qci  = mic_fab_vars[MicVar::qci];
  auto qpl  = mic_fab_vars[MicVar::qpl];
  auto qpi  = mic_fab_vars[MicVar::qpi];
  auto qg  = mic_fab_vars[MicVar::qg];
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
     auto qg_array    = qg->array(mfi);
     auto tabs_array  = tabs->array(mfi);

     const auto& box3d = mfi.tilebox();

     amrex::ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         qv_array(i,j,k)  = qt_array(i,j,k) - qn_array(i,j,k);
         amrex::Real omn  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tbgmin)*a_bg));
         qcl_array(i,j,k) = qn_array(i,j,k)*omn;
         qci_array(i,j,k) = qn_array(i,j,k)*(1.0-omn);
         amrex::Real omp  = std::max(0.0, std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
         qpl_array(i,j,k) = qp_array(i,j,k)*omp;
         qpi_array(i,j,k) = qp_array(i,j,k)*(1.0-omp);

         // TODO: If omp above is bound by [0, 1], shouldn't this always be 0?
         // Graupel == precip total - rain - snow (but must be >= 0)
         qg_array(i,j,k)  = std::max(0.0, qp_array(i,j,k)-qpl_array(i,j,k)-qpi_array(i,j,k));
     });
  }
}


/**
 * Computes diagnostic quantities like temperature and pressure
 * from the existing Microphysics variables.
 */
void SAM::Compute_Tabs_and_Pres ()
{
    auto rho   = mic_fab_vars[MicVar::rho];
    auto theta = mic_fab_vars[MicVar::theta];
    auto qv    = mic_fab_vars[MicVar::qv];
    auto tabs  = mic_fab_vars[MicVar::tabs];
    auto pres  = mic_fab_vars[MicVar::pres];

    for ( amrex::MFIter mfi(*tabs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto rho_array   = rho->array(mfi);
        auto theta_array = theta->array(mfi);
        auto qv_array    = qv->array(mfi);
        auto tabs_array  = tabs->array(mfi);
        auto pres_array  = pres->array(mfi);

        amrex::ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::Real RhoTheta = rho_array(i,j,k)*theta_array(i,j,k);
            tabs_array(i,j,k)  = getTgivenRandRTh(rho_array(i,j,k), RhoTheta, qv_array(i,j,k));
            pres_array(i,j,k)  = getPgivenRTh(RhoTheta, qv_array(i,j,k))/100.;
        });
    }
}
