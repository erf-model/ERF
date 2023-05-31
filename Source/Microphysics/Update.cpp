
#include "Microphysics.H"
#include "IndexDefines.H"
#include "TileNoZ.H"

/**
 * Updates conserved and microphysics variables in the provided MultiFabs from
 * the internal MultiFabs that store Microphysics module data.
 *
 * @param[out] cons_in Conserved variables
 * @param[out] qv_in Water vapor component
 * @param[out] qc_in Cloud liquid component
 * @param[out] qi_in Cloud ice component
 * @param[out] qrain_in Precipitation liquid component
 * @param[out] qsnow_in Precipitation ice component
 * @param[out] qgraup_in Graup component
 */
void Microphysics::Update(amrex::MultiFab& cons_in,
                          amrex::MultiFab& qv_in,
                          amrex::MultiFab& qc_in,
                          amrex::MultiFab& qi_in,
                          amrex::MultiFab& qrain_in,
                          amrex::MultiFab& qsnow_in,
                          amrex::MultiFab& qgraup_in)
{
  // copy multifab data to qc, qv, and qi
  amrex::MultiFab::Copy(qv_in    , *mic_fab_vars[MicVar::qv],  0, 0, 1, mic_fab_vars[MicVar::qv]->nGrowVect());
  amrex::MultiFab::Copy(qc_in    , *mic_fab_vars[MicVar::qcl], 0, 0, 1, mic_fab_vars[MicVar::qcl]->nGrowVect());
  amrex::MultiFab::Copy(qi_in    , *mic_fab_vars[MicVar::qci], 0, 0, 1, mic_fab_vars[MicVar::qci]->nGrowVect());
  amrex::MultiFab::Copy(qrain_in , *mic_fab_vars[MicVar::qpl], 0, 0, 1, mic_fab_vars[MicVar::qpl]->nGrowVect());
  amrex::MultiFab::Copy(qsnow_in , *mic_fab_vars[MicVar::qpi], 0, 0, 1, mic_fab_vars[MicVar::qpi]->nGrowVect());
  amrex::MultiFab::Copy(qgraup_in, *mic_fab_vars[MicVar::qpi], 0, 0, 1, mic_fab_vars[MicVar::qci]->nGrowVect());

  // get the temperature, density, theta, qt and qp from input
  for ( amrex::MFIter mfi(cons_in,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     auto states_array = cons_in.array(mfi);
     auto rho_array    = mic_fab_vars[MicVar::rho]->array(mfi);
     auto theta_array  = mic_fab_vars[MicVar::theta]->array(mfi);
     auto qt_array     = mic_fab_vars[MicVar::qt]->array(mfi);
     auto qp_array     = mic_fab_vars[MicVar::qp]->array(mfi);
     auto qpl_array    = mic_fab_vars[MicVar::qpl]->array(mfi);
     auto qpi_array    = mic_fab_vars[MicVar::qpi]->array(mfi);
     auto qgraup_array = qgraup_in.array(mfi);

     const auto& box3d = mfi.tilebox();

     // get potential total density, temperature, qt, qp
     amrex::ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
       states_array(i,j,k,Rho_comp)      = rho_array(i,j,k);
       states_array(i,j,k,RhoTheta_comp) = rho_array(i,j,k)*theta_array(i,j,k);
       states_array(i,j,k,RhoQt_comp)    = rho_array(i,j,k)*qt_array(i,j,k);
       states_array(i,j,k,RhoQp_comp)    = rho_array(i,j,k)*qp_array(i,j,k);
       qgraup_array(i,j,k)               = std::max(0.0, qp_array(i,j,k)-qpl_array(i,j,k)-qpi_array(i,j,k)); // negative is unphysical
     });
  }

  // Fill interior ghost cells and periodic boundaries
  cons_in.FillBoundary(m_geom.periodicity());

  qv_in.FillBoundary(m_geom.periodicity());
  qc_in.FillBoundary(m_geom.periodicity());
  qi_in.FillBoundary(m_geom.periodicity());
  qrain_in. FillBoundary(m_geom.periodicity());
  qsnow_in. FillBoundary(m_geom.periodicity());
  qgraup_in.FillBoundary(m_geom.periodicity());
}


