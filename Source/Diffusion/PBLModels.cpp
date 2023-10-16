#include "ABLMost.H"
#include "DirectionSelector.H"
#include "Diffusion.H"
#include "ERF_Constants.H"
#include "TurbStruct.H"

/**
 * Function to compute turbulent viscosity with PBL.
 *
 * @param[in] xvel velocity in x-dir
 * @param[in] yvel velocity in y-dir
 * @param[in] cons_in cell center conserved quantities
 * @param[out] eddyViscosity holds turbulent viscosity
 * @param[in] geom problem geometry
 * @param[in] turbChoice container with turbulence parameters
 * @param[in] most pointer to Monin-Obukhov class if instantiated
 */
void
ComputeTurbulentViscosityPBL (const amrex::MultiFab& xvel,
                              const amrex::MultiFab& yvel,
                              const amrex::MultiFab& cons_in,
                              const amrex::MultiFab& cons_old,
                              amrex::MultiFab& eddyViscosity,
                              const amrex::Geometry& geom,
                              const TurbChoice& turbChoice,
                              std::unique_ptr<ABLMost>& most,
                              const amrex::BCRec* bc_ptr,
                              bool /*vert_only*/)
{
  // MYNN Level 2.5 PBL Model
  if (turbChoice.pbl_type == PBLType::MYNN25) {

    const amrex::Real A1 = turbChoice.pbl_A1;
    const amrex::Real A2 = turbChoice.pbl_A2;
    const amrex::Real B1 = turbChoice.pbl_B1;
    const amrex::Real B2 = turbChoice.pbl_B2;
    const amrex::Real C1 = turbChoice.pbl_C1;
    const amrex::Real C2 = turbChoice.pbl_C2;
    const amrex::Real C3 = turbChoice.pbl_C3;
    //const amrex::Real C4 = turbChoice.pbl_C4;
    const amrex::Real C5 = turbChoice.pbl_C5;

    // Dirichlet flags to switch derivative stencil
    bool c_ext_dir_on_zlo = ( (bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir) );
    bool c_ext_dir_on_zhi = ( (bc_ptr[BCVars::cons_bc].lo(5) == ERFBCType::ext_dir) );
    bool u_ext_dir_on_zlo = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) );
    bool u_ext_dir_on_zhi = ( (bc_ptr[BCVars::xvel_bc].lo(5) == ERFBCType::ext_dir) );
    bool v_ext_dir_on_zlo = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) );
    bool v_ext_dir_on_zhi = ( (bc_ptr[BCVars::yvel_bc].lo(5) == ERFBCType::ext_dir) );

    // Epsilon
    amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const amrex::Box &bx = mfi.growntilebox(1);
      const amrex::Array4<amrex::Real const > &cell_data     = cons_in.array(mfi);
      const amrex::Array4<amrex::Real const > &cell_data_old = cons_old.array(mfi);
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
      amrex::FArrayBox qturb(bx,1); amrex::FArrayBox qturb_old(bx,1);
      const amrex::Array4<amrex::Real> qint = qintegral.array();
      const amrex::Array4<amrex::Real> qvel = qturb.array();
      const amrex::Array4<amrex::Real> qvel_old = qturb_old.array();

      amrex::ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), bx,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept
      {
          qvel(i,j,k)     = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp));
          qvel_old(i,j,k) = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp) + eps);
          AMREX_ASSERT_WITH_MESSAGE(qvel(i,j,k) > 0.0, "QKE must have a positive value");
          AMREX_ASSERT_WITH_MESSAGE(qvel_old(i,j,k) > 0.0, "Old QKE must have a positive value");

          const amrex::Real Zval = gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
          if (sbx.contains(i,j,k)) {
              amrex::Gpu::deviceReduceSum(&qint(i,j,0,0), Zval*qvel(i,j,k), handler);
              amrex::Gpu::deviceReduceSum(&qint(i,j,0,1), qvel(i,j,k), handler);
          }
      });

      amrex::Real dz_inv = geom.InvCellSize(2);
      int izmin = geom.Domain().smallEnd(2);
      int izmax = geom.Domain().bigEnd(2);

      // Spatially varying MOST
      amrex::Real d_kappa   = KAPPA;
      amrex::Real d_gravity = CONST_GRAV;

      const auto& t_mean_mf = most->get_mac_avg(0,2); // TODO: IS THIS ACTUALLY RHOTHETA
      const auto& u_star_mf = most->get_u_star(0);    // Use coarsest level
      const auto& t_star_mf = most->get_t_star(0);    // Use coarsest level

      const auto& tm_arr     = t_mean_mf->array(mfi); // TODO: IS THIS ACTUALLY RHOTHETA
      const auto& u_star_arr = u_star_mf->array(mfi);
      const auto& t_star_arr = t_star_mf->array(mfi);

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
          // Compute some partial derivatives that we will need (second order)
          // U and V derivatives are interpolated to account for staggered grid
          amrex::Real dthetadz, dudz, dvdz;
          if ( k==izmax && c_ext_dir_on_zhi ) {
              dthetadz = (1.0/3.0)*(-cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp)
                             - 3.0 * cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp)
                             + 4.0 * cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp) )*dz_inv;
          } else if ( k==izmin && c_ext_dir_on_zlo ) {
              dthetadz = (1.0/3.0)*( cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp)
                             + 3.0 * cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp)
                             - 4.0 * cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp) )*dz_inv;
          } else {
              dthetadz = 0.5*(cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp)
                            - cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp))*dz_inv;
          }

          if ( k==izmax && u_ext_dir_on_zhi ) {
              dudz = (1.0/6.0)*( (-uvel(i  ,j,k-1) - 3.0 * uvel(i  ,j,k  ) + 4.0 * uvel(i  ,j,k+1))
                               + (-uvel(i+1,j,k-1) - 3.0 * uvel(i+1,j,k  ) + 4.0 * uvel(i+1,j,k+1)) )*dz_inv;
          } else if ( k==izmin && u_ext_dir_on_zlo ) {
              dudz = (1.0/6.0)*( (uvel(i  ,j,k+1) + 3.0 * uvel(i  ,j,k  ) - 4.0 * uvel(i  ,j,k-1))
                               + (uvel(i+1,j,k+1) + 3.0 * uvel(i+1,j,k  ) - 4.0 * uvel(i+1,j,k-1)) )*dz_inv;
          } else {
              dudz = 0.25*(uvel(i,j,k+1) - uvel(i,j,k-1) + uvel(i+1,j,k+1) - uvel(i+1,j,k-1))*dz_inv;
          }

          if ( k==izmax && v_ext_dir_on_zhi ) {
              dvdz = (1.0/6.0)*( (-vvel(i,j  ,k-1) - 3.0 * vvel(i,j  ,k  ) + 4.0 * vvel(i,j  ,k+1))
                               + (-vvel(i,j+1,k-1) - 3.0 * vvel(i,j+1,k  ) + 4.0 * vvel(i,j+1,k+1)) )*dz_inv;
          } else if ( k==izmin && v_ext_dir_on_zlo ) {
              dvdz = (1.0/6.0)*( (vvel(i,j  ,k+1) + 3.0 * vvel(i,j  ,k  ) - 4.0 * vvel(i,j  ,k-1))
                               + (vvel(i,j+1,k+1) + 3.0 * vvel(i,j+1,k  ) - 4.0 * vvel(i,j+1,k-1)) )*dz_inv;
          } else {
              dvdz = 0.25*(vvel(i,j,k+1) - vvel(i,j,k-1) + vvel(i,j+1,k+1) - vvel(i,j+1,k-1))*dz_inv;
          }

          // Spatially varying MOST
          amrex::Real surface_heat_flux = -u_star_arr(i,j,0) * t_star_arr(i,j,0);
          amrex::Real theta0            = tm_arr(i,j,0); // TODO: IS THIS ACTUALLY RHOTHETA
          amrex::Real l_obukhov;
          if (std::abs(surface_heat_flux) > eps) {
              l_obukhov = ( theta0 * u_star_arr(i,j,0) * u_star_arr(i,j,0) ) /
                          ( d_kappa * d_gravity * t_star_arr(i,j,0) );
          } else {
              l_obukhov = std::numeric_limits<amrex::Real>::max();
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
              amrex::Real N_brunt_vaisala = std::sqrt(CONST_GRAV/theta0 * dthetadz);
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

          // NOTE: Level 2 limiting from balance of production and dissipation.
          //       K_turb has a setval of 0.0 when the MF is created (NOT EACH STEP).
          //       We do this inline to avoid storing qe^2 at each cell.
          amrex::Real shearProd = dudz*dudz + dvdz*dvdz;
          amrex::Real buoyProd  = -(CONST_GRAV/theta0) * dthetadz;
          amrex::Real lSM       = K_turb(i,j,k,EddyDiff::Mom_v)   / (qvel_old(i,j,k) + eps);
          amrex::Real lSH       = K_turb(i,j,k,EddyDiff::Theta_v) / (qvel_old(i,j,k) + eps);
          amrex::Real qe2       = B1 * l_comb * ( lSM * shearProd + lSH * buoyProd );
          amrex::Real qe        = (qe2 < 0.0) ? 0.0 : std::sqrt(qe2);

          // Compute non-dimensional parameters
          amrex::Real one_m_alpha  = (qvel(i,j,k) > qe) ? 1.0 : qvel(i,j,k) / (qe + eps);
          amrex::Real one_m_alpha2 = one_m_alpha * one_m_alpha;
          amrex::Real l2_over_q2   = l_comb*l_comb/(qvel(i,j,k)*qvel(i,j,k));
          amrex::Real GM = l2_over_q2 * shearProd;
          amrex::Real GH = l2_over_q2 * buoyProd;
          amrex::Real E1 = 1.0 + one_m_alpha2 * ( 6.0*A1*A1*GM - 9.0*A1*A2*(1.0-C2)*GH );
          amrex::Real E2 = one_m_alpha2 * ( -3.0*A1*(4.0*A1 + 3.0*A2*(1.0-C5))*(1.0-C2)*GH );
          amrex::Real E3 = one_m_alpha2 * ( 6.0*A2*A1*GM );
          amrex::Real E4 = 1.0 + one_m_alpha2 * ( -12.0*A2*A1*(1.0-C2)*GH - 3.0*A2*B2*(1.0-C3)*GH );
          amrex::Real R1 = one_m_alpha * ( A1*(1.0-3.0*C1) );
          amrex::Real R2 = one_m_alpha * A2;

          amrex::Real SM = (R2*E2 - R1*E4)/(E2*E3 - E1*E4);
          amrex::Real SH = (R1*E3 - A2*E1)/(E2*E3 - E1*E4);
          amrex::Real SQ = 3.0 * SM;

          // Finally, compute the eddy viscosity/diffusivities
          const amrex::Real rho = cell_data(i,j,k,Rho_comp);
          K_turb(i,j,k,EddyDiff::Mom_v)   = rho * l_comb * qvel(i,j,k) * SM * 0.5; // 0.5 for mu_turb
          K_turb(i,j,k,EddyDiff::Theta_v) = rho * l_comb * qvel(i,j,k) * SH;
          K_turb(i,j,k,EddyDiff::QKE_v)   = rho * l_comb * qvel(i,j,k) * SQ;

          K_turb(i,j,k,EddyDiff::PBL_lengthscale) = l_comb;
          // TODO: How should this be done for other components (scalars, moisture)
      });
    }
  } else if (turbChoice.pbl_type == PBLType::YSU) {
      amrex::Error("YSU Model not implemented yet");
  }
}
