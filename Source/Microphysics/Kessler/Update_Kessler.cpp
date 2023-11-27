
#include "Microphysics.H"
#include "IndexDefines.H"
#include "TileNoZ.H"

/**
 * Updates conserved and microphysics variables in the provided MultiFabs from
 * the internal MultiFabs that store Microphysics module data.
 *
 * @param[out] cons Conserved variables
 * @param[out] qmoist: qv, qc, qi, qr, qs, qg
 */
void Kessler::Update (amrex::MultiFab& cons,
                  amrex::MultiFab& qmoist)
{
  // copy multifab data to qc, qv, and qi
  amrex::MultiFab::Copy(qmoist, *mic_fab_vars[MicVar_Kess::qv],  0, 0, 1, mic_fab_vars[MicVar_Kess::qv]->nGrowVect());  // vapor
  amrex::MultiFab::Copy(qmoist, *mic_fab_vars[MicVar_Kess::qcl], 0, 1, 1, mic_fab_vars[MicVar_Kess::qcl]->nGrowVect()); // cloud water
  amrex::MultiFab::Copy(qmoist, *mic_fab_vars[MicVar_Kess::qci], 0, 2, 1, mic_fab_vars[MicVar_Kess::qci]->nGrowVect()); // cloud ice
  amrex::MultiFab::Copy(qmoist, *mic_fab_vars[MicVar_Kess::qpl], 0, 3, 1, mic_fab_vars[MicVar_Kess::qpl]->nGrowVect()); // rain
  amrex::MultiFab::Copy(qmoist, *mic_fab_vars[MicVar_Kess::qpi], 0, 4, 1, mic_fab_vars[MicVar_Kess::qpi]->nGrowVect()); // snow

  // Don't need to copy this since it is filled below
  // amrex::MultiFab::Copy(qmoist, *mic_fab_vars[MicVar_Kess::qpi], 0, 5, 1, mic_fab_vars[MicVar_Kess::qci]->nGrowVect()); // graupel

  amrex::MultiFab qgraup_mf(qmoist, amrex::make_alias, 5, 1);

  // Get the temperature, density, theta, qt and qp from input
  for ( amrex::MFIter mfi(cons,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     auto states_arr = cons.array(mfi);

     auto rho_arr    = mic_fab_vars[MicVar_Kess::rho]->array(mfi);
     auto theta_arr  = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
     auto qt_arr     = mic_fab_vars[MicVar_Kess::qt]->array(mfi);
     auto qp_arr     = mic_fab_vars[MicVar_Kess::qp]->array(mfi);
     auto qpl_arr    = mic_fab_vars[MicVar_Kess::qpl]->array(mfi);
     auto qpi_arr    = mic_fab_vars[MicVar_Kess::qpi]->array(mfi);

     auto qgraup_arr= qgraup_mf.array(mfi);

     const auto& box3d = mfi.tilebox();

     // get potential total density, temperature, qt, qp
     amrex::ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
       states_arr(i,j,k,Rho_comp)      = rho_arr(i,j,k);
       states_arr(i,j,k,RhoTheta_comp) = rho_arr(i,j,k)*theta_arr(i,j,k);
       states_arr(i,j,k,RhoQt_comp)    = rho_arr(i,j,k)*qt_arr(i,j,k);
       states_arr(i,j,k,RhoQp_comp)    = rho_arr(i,j,k)*qp_arr(i,j,k);

       // Graupel == precip total - rain - snow (but must be >= 0)
       qgraup_arr(i,j,k)  = 0.0;// 
     });
  }

  // Fill interior ghost cells and periodic boundaries
  cons.FillBoundary(m_geom.periodicity());
  qmoist.FillBoundary(m_geom.periodicity());
}


