
#include "Microphysics.H"
#include "IndexDefines.H"
#include "EOS.H"

void Microphysics::Update(amrex::MultiFab& cons_in,
                          amrex::MultiFab& qv_in,
                          amrex::MultiFab& qc_in,
                          amrex::MultiFab& qi_in)
{
  // copy multifab data to qc, qv, and qi
  amrex::MultiFab::Copy(qv_in, *mic_fab_vars[MicVar::qv],  0, 0, 1, mic_fab_vars[MicVar::qv]->nGrowVect());
  amrex::MultiFab::Copy(qc_in, *mic_fab_vars[MicVar::qcl], 0, 0, 1, mic_fab_vars[MicVar::qcl]->nGrowVect());
  amrex::MultiFab::Copy(qi_in, *mic_fab_vars[MicVar::qci], 0, 0, 1, mic_fab_vars[MicVar::qci]->nGrowVect());

  // get the temperature, dentisy, theta, qt and qp from input
  for ( amrex::MFIter mfi(cons_in,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     auto states_array = cons_in.array(mfi);
     auto rho_array    = mic_fab_vars[MicVar::rho]->array(mfi);
     auto theta_array  = mic_fab_vars[MicVar::theta]->array(mfi);
     auto qt_array     = mic_fab_vars[MicVar::qt]->array(mfi);
     auto qp_array     = mic_fab_vars[MicVar::qp]->array(mfi);

     const auto& box3d = mfi.tilebox();

     // get potential total density, temperature, qt, qp
     amrex::ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
       states_array(i,j,k,Rho_comp)      = rho_array(i,j,k);
       states_array(i,j,k,RhoTheta_comp) = rho_array(i,j,k)*theta_array(i,j,k);
       states_array(i,j,k,RhoQt_comp)    = rho_array(i,j,k)*qt_array(i,j,k);
       states_array(i,j,k,RhoQp_comp)    = rho_array(i,j,k)*qp_array(i,j,k);
     });
  }

  // fill the boundary
  cons_in.FillBoundary(m_geom.periodicity());
  qv_in.FillBoundary(m_geom.periodicity());
  qc_in.FillBoundary(m_geom.periodicity());
  qi_in.FillBoundary(m_geom.periodicity());
}


