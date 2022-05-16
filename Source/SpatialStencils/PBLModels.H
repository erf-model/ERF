#ifndef _PBLMODELS_H_
#define _PBLMODELS_H_

#include "ABLMost.H"
#include "DirectionSelector.H"

inline
void ComputeTurbulentViscosityPBL(const amrex::MultiFab& xvel,
                                  const amrex::MultiFab& yvel,
                                  const amrex::MultiFab& zvel,
                                  const amrex::MultiFab& cons_in,
                                  amrex::MultiFab& eddyViscosity,
                                  const amrex::Geometry& geom,
                                  const SolverChoice& solverChoice,
                                  std::unique_ptr<ABLMost>& most)
{
  // MYNN Level 2.5 PBL Model
  if (solverChoice.pbl_type == PBLType::MYNN25) {

    const amrex::Real A1 = solverChoice.pbl_A1;
    const amrex::Real A2 = solverChoice.pbl_A2;
    //const amrex::Real B1 = solverChoice.pbl_B1;
    const amrex::Real B2 = solverChoice.pbl_B2;
    const amrex::Real C1 = solverChoice.pbl_C1;
    const amrex::Real C2 = solverChoice.pbl_C2;
    const amrex::Real C3 = solverChoice.pbl_C3;
    //const amrex::Real C4 = solverChoice.pbl_C4;
    const amrex::Real C5 = solverChoice.pbl_C5;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const amrex::Box &bx = mfi.growntilebox(1);
      const amrex::Array4<amrex::Real const > &cell_data = cons_in.array(mfi);
      const amrex::Array4<amrex::Real> &K_turb = eddyViscosity.array(mfi);
      const amrex::Array4<amrex::Real const> &uvel = xvel.array(mfi);
      const amrex::Array4<amrex::Real const> &vvel = yvel.array(mfi);

      // Compute some quantities that are constant in each column
      // Sbox is shrunk to only include the interior of the domain in the vertical direction to compute integrals
      // NOTE: Here we requite that sbx covers the entire vertical domain
      const amrex::Box &dbx = geom.Domain();
      amrex::Box sbx(bx.smallEnd(), bx.bigEnd());
      sbx.grow(2,-1);
      AMREX_ALWAYS_ASSERT(sbx.smallEnd(2) == dbx.smallEnd(2) && sbx.bigEnd(2) == dbx.bigEnd(2));

      const amrex::GeometryData gdata = geom.data();

      const amrex::Box xybx = PerpendicularBox<ZDir>(bx, amrex::IntVect{0,0,0});
      amrex::FArrayBox qintegral(xybx,2);
      qintegral.setVal<amrex::RunOn::Device>(0.0);
      amrex::FArrayBox qturb(bx,1);
      const amrex::Array4<amrex::Real> qint = qintegral.array();
      const amrex::Array4<amrex::Real> qvel= qturb.array();

      amrex::ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), bx,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept
      {
          const amrex::Real Zval = gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
          const amrex::Real rho = cell_data(i,j,k,Rho_comp);
          qvel(i,j,k) = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / rho);
          // We will divide by qvel later
          AMREX_ASSERT_WITH_MESSAGE(qvel(i,j,k) > 0.0, "QKE must have a positive value");
          if (sbx.contains(i,j,k)) {
              amrex::Gpu::deviceReduceSum(&qint(i,j,0,0), Zval*qvel(i,j,k), handler);
              amrex::Gpu::deviceReduceSum(&qint(i,j,0,1), qvel(i,j,k), handler);
          }
      });

      amrex::Real dz_inv = geom.InvCellSize(2);
      int izmin = geom.Domain().smallEnd(2);
      int izmax = geom.Domain().bigEnd(2);
      const amrex::Real l_obukhov = most->obukhov_len;
      const amrex::Real surface_heat_flux = most->surf_temp_flux;
      const amrex::Real theta0 = most->theta_mean; //(TODO: IS THIS ACTUALLY RHOTHETA)
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
          // Compute some partial derivatives that we will need (1st order at domain boundary)
          // U and V derivatives are interpolated to account for staggered grid
          amrex::Real dthetadz, dudz, dvdz;
          if (k == izmax) {
              dthetadz = (cell_data(i,j,k,RhoTheta_comp)/cell_data(i,j,k,Rho_comp) -
                          cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp))*dz_inv;
              dudz = 0.5*(uvel(i,j,k) - uvel(i,j,k-1) + uvel(i+1,j,k) - uvel(i+1,j,k-1))*dz_inv;
              dvdz = 0.5*(vvel(i,j,k) - vvel(i,j,k-1) + vvel(i,j+1,k) - vvel(i,j+1,k-1))*dz_inv;
          } else if (k == izmin){
              dthetadz = (cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp) -
                          cell_data(i,j,k,RhoTheta_comp)/cell_data(i,j,k,Rho_comp))*dz_inv;
              dudz = 0.5*(uvel(i,j,k+1) - uvel(i,j,k) + uvel(i+1,j,k+1) - uvel(i+1,j,k))*dz_inv;
              dvdz = 0.5*(vvel(i,j,k+1) - vvel(i,j,k) + vvel(i,j+1,k+1) - vvel(i,j+1,k))*dz_inv;
          } else {
              dthetadz = 0.5*(cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp) -
                              cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp))*dz_inv;
              dudz = 0.25*(uvel(i,j,k+1) - uvel(i,j,k-1) + uvel(i+1,j,k+1) - uvel(i+1,j,k-1))*dz_inv;
              dvdz = 0.25*(vvel(i,j,k+1) - vvel(i,j,k-1) + vvel(i,j+1,k+1) - vvel(i,j+1,k-1))*dz_inv;
          }

          // First Length Scale
          AMREX_ASSERT(l_obukhov != 0);
          const amrex::Real zval = gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
          const amrex::Real zeta = zval/l_obukhov;
          amrex::Real l_S;
          if (zeta >= 1.0) {
              l_S = KAPPA*zval/3.7;
          } else if (zeta >= 0) {
              l_S = KAPPA*zval/(1+2.7*zeta);
          } else {
              l_S = KAPPA*zval*std::pow(1.0 - 100.0 * zeta, 0.2);
          }

          // Second Length Scale
          amrex::Real l_T;
          if (qint(i,j,0,1) > 0.0) {
              l_T = 0.23*qint(i,j,0,0)/qint(i,j,0,1);
          } else {
              l_T = std::numeric_limits<amrex::Real>::max();
          }

          // Third Length Scale
          amrex::Real l_B;
          if (dthetadz > 0) {
              amrex::Real N_brunt_vaisala = CONST_GRAV/theta0 * std::sqrt(dthetadz);
              if (zeta < 0) {
                  amrex::Real qc = CONST_GRAV/theta0 * surface_heat_flux * l_T;
                  qc = std::pow(qc,1.0/3.0);
                  l_B = (1.0 + 5.0*std::sqrt(qc/(N_brunt_vaisala * l_T))) * qvel(i,j,k)/N_brunt_vaisala;
              } else {
                  l_B = qvel(i,j,k) / N_brunt_vaisala;
              }
          } else {
              l_B = std::numeric_limits<amrex::Real>::max();
          }

          // Overall Length Scale
          amrex::Real l_comb = 1.0 / (1.0/l_S + 1.0/l_T + 1.0/l_B);

          // Compute non-dimensional parameters
          amrex::Real l2_over_q2 = l_comb*l_comb/(qvel(i,j,k)*qvel(i,j,k));
          amrex::Real GM = l2_over_q2 * (dudz*dudz + dvdz*dvdz);
          amrex::Real GH = -l2_over_q2 / theta0 * dthetadz;
          amrex::Real E1 = 1.0 + 6.0*A1*A1*GM - 9.0*A1*A2*(1.0-C2)*GH;
          amrex::Real E2 = -3.0*A1*(4.0*A1 + 3.0*A2*(1.0-C5))*(1.0-C2)*GH;
          amrex::Real E3 = 6.0*A2*A1*GM;
          amrex::Real E4 = 1.0 - 12.0*A2*A1*(1.0-C2)*GH -3.0*A2*B2*(1.0-C3)*GH;
          amrex::Real R1 = A1*(1.0-3.0*C1);

          amrex::Real SM = (A2*E2 - R1*E4)/(E2*E3 - E1*E4);
          amrex::Real SH = (R1*E3 - A2*E1)/(E2*E3 - E1*E4);
          amrex::Real SQ = 3.0 * SM;

          // Finally, compute the eddy viscosity/diffusivities
          const amrex::Real rho = cell_data(i,j,k,Rho_comp);
          K_turb(i,j,k,EddyDiff::Mom_v)   = rho * l_comb * qvel(i,j,k) * SM;
          K_turb(i,j,k,EddyDiff::Theta_v) = rho * l_comb * qvel(i,j,k) * SH;
          K_turb(i,j,k,EddyDiff::QKE_v)   = rho * l_comb * qvel(i,j,k) * 3.0 * SQ;

          K_turb(i,j,k,EddyDiff::PBL_lengthscale) = l_comb;
          // TODO: How should this be done for other components (scalars, moisture)
      });
    }
  }
}

AMREX_GPU_DEVICE
inline
amrex::Real ComputeQKESourceTerms(int i, int j, int k,
                           const amrex::Array4<const amrex::Real>& uvel,
                           const amrex::Array4<const amrex::Real>& vvel,
                           const amrex::Array4<const amrex::Real>& cell_data,
                           const amrex::Array4<const amrex::Real>& cell_prim,
                           const amrex::Array4<amrex::Real>& K_turb,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                           const amrex::Box& domain,
                           const SolverChoice& solverChoice,
                           const amrex::Real theta_mean)
{
  // Compute some relevant derivatives
  amrex::Real dz_inv = cellSizeInv[2];
  int izmin = domain.smallEnd(2);
  int izmax = domain.bigEnd(2);
  amrex::Real dthetadz, dudz, dvdz;
  amrex::Real source_term = 0.0;
  if (k == izmax) {
    dthetadz = (cell_prim(i,j,k,PrimTheta_comp) - cell_prim(i,j,k-1,PrimTheta_comp))*dz_inv;
    dudz = 0.5*(uvel(i,j,k) - uvel(i,j,k-1) + uvel(i+1,j,k) - uvel(i+1,j,k-1))*dz_inv;
    dvdz = 0.5*(vvel(i,j,k) - vvel(i,j,k-1) + vvel(i,j+1,k) - vvel(i,j+1,k-1))*dz_inv;
  } else if (k == izmin){
    dthetadz = (cell_prim(i,j,k+1,PrimTheta_comp) - cell_prim(i,j,k,PrimTheta_comp))*dz_inv;
    dudz = 0.5*(uvel(i,j,k+1) - uvel(i,j,k) + uvel(i+1,j,k+1) - uvel(i+1,j,k))*dz_inv;
    dvdz = 0.5*(vvel(i,j,k+1) - vvel(i,j,k) + vvel(i,j+1,k+1) - vvel(i,j+1,k))*dz_inv;
  } else {
    dthetadz = 0.5*(cell_prim(i,j,k+1,PrimTheta_comp) - cell_prim(i,j,k-1,PrimTheta_comp))*dz_inv;
    dudz = 0.25*(uvel(i,j,k+1) - uvel(i,j,k-1) + uvel(i+1,j,k+1) - uvel(i+1,j,k-1))*dz_inv;
    dvdz = 0.25*(vvel(i,j,k+1) - vvel(i,j,k-1) + vvel(i,j+1,k+1) - vvel(i,j+1,k-1))*dz_inv;
  }

  // Bouyancy
  source_term += 2*CONST_GRAV/theta_mean*K_turb(i,j,k,EddyDiff::Theta_v)*dthetadz;

  // Production
  source_term += K_turb(i,j,k,EddyDiff::Mom_v) * (dudz*dudz + dvdz*dvdz);

  // Dissipation
  amrex::Real qke = cell_prim(i,j,k,PrimQKE_comp);
  if (std::abs(qke) > 0.0) {
    source_term += 2.0 * cell_data(i,j,k,Rho_comp) * std::pow(qke,1.5) /
      (solverChoice.pbl_B1 * K_turb(i,j,k,EddyDiff::PBL_lengthscale));
  }

  return source_term;
}
#endif
