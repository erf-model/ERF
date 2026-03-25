#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_SolveTridiag.H"
#include "ERF_GetRhoAlpha.H"
#include "ERF_GetRhoAlphaForFaces.H"

using namespace amrex;

/**
 * Function for computing the implicit contribution to the vertical diffusion
 * of theta, with terrain.
 *
 * @param[in   ] bx cell-centered box to loop over
 * @param[in   ] domain box of the whole domain
 * @param[in   ] dt time step
 * @param[in   ] bc_neumann_vals values of derivatives if bc_type == Neumann
 * @param[inout] cell_data conserved cell-centered rho, rho theta
 * @param[in   ] z_nd nodal array of z
 * @param[in   ] detJ Jacobian determinant
 * @param[in   ] cellSizeInv inverse cell size array
 * @param[inout] hfx_z heat flux in z-dir
 * @param[in   ] mu_turb turbulent viscosity
 * @param[in   ] solverChoice container of parameters
 * @param[in   ] bc_ptr container with boundary conditions
 * @param[in   ] use_SurfLayer whether we have turned on subgrid diffusion
 * @param[in   ] implicit_fac if 1 then fully implicit; if 0 then fully explicit
 */
void
ImplicitDiffForStateLU_T (const Box& bx,
                          const Box& domain,
                          const int level,
                          const Real dt,
                          const GpuArray<Real, AMREX_SPACEDIM*2>& bc_neumann_vals,
                          const Array4<      Real>& cell_data,
                          const Array4<const Real>& z_nd,
                          const Array4<const Real>& detJ,
                          const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                          const Array4<const Real>& hfx_z,
                          const Array4<const Real>& mu_turb,
                          const SolverChoice &solverChoice,
                          const BCRec* bc_ptr,
                          const bool use_SurfLayer,
                          const Real implicit_fac)
{
    BL_PROFILE_VAR("ImplicitDiffForState_T()",ImplicitDiffForState_T);

    // setup quantities for getRhoAlpha()
#include "ERF_SetupVertDiff.H"
    const int         n = RhoTheta_comp;
    const int qty_index = RhoTheta_comp;
    const int prim_index = qty_index - 1;
    const int prim_scal_index = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ? PrimScalar_comp : prim_index;

    // Box bounds
    int ilo = bx.smallEnd(0);
    int ihi = bx.bigEnd(0);
    int jlo = bx.smallEnd(1);
    int jhi = bx.bigEnd(1);
    int klo = bx.smallEnd(2);
    int khi = bx.bigEnd(2);
    amrex::ignore_unused(ilo, ihi, jlo, jhi);

    // Temporary FABs for tridiagonal solve (allocated on column)
    //   A[k] * x[k-1] + B[k] * x[k] + C[k+1] = RHS[k]

    // With LU decomposition, M * x = r is written as L * U * x = r with U * x = rho
    // We then first have L * rho = r and U * x = rho
    amrex::FArrayBox RHS_fab, soln_fab, coeffG_fab;
           RHS_fab.resize(bx,1, amrex::The_Async_Arena());
          soln_fab.resize(bx,1, amrex::The_Async_Arena());
        coeffG_fab.resize(bx,1, amrex::The_Async_Arena());
    auto const& RHS_a        =        RHS_fab.array();
    auto const& soln_a       =       soln_fab.array();
    auto const& coeffG_a     =     coeffG_fab.array();

    Real dz_inv = cellSizeInv[2];

    int bc_comp = qty_index;
    bool foextrap_on_zlo = (bc_ptr[bc_comp].lo(2) == ERFBCType::foextrap);
    bool foextrap_on_zhi = (bc_ptr[bc_comp].hi(2) == ERFBCType::foextrap);
    bool neumann_on_zlo  = (bc_ptr[bc_comp].lo(2) == ERFBCType::neumann);
    bool neumann_on_zhi  = (bc_ptr[bc_comp].hi(2) == ERFBCType::neumann);
    amrex::ignore_unused(foextrap_on_zlo, foextrap_on_zhi);

    AMREX_ASSERT_WITH_MESSAGE(foextrap_on_zlo || neumann_on_zlo || use_SurfLayer,
                              "Unexpected lower BC used with implicit vertical diffusion");
    AMREX_ASSERT_WITH_MESSAGE(foextrap_on_zhi || neumann_on_zhi,
                              "Unexpected upper BC used with implicit vertical diffusion");

    Real Fact = implicit_fac * dt * dz_inv;

#ifdef AMREX_USE_GPU
    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
#else
    for (int j(jlo); j<=jhi; ++j) {
        for (int i(ilo); i<=ihi; ++i) {
#endif
            // Notes: The detJ for divergence was multiplied through.
            //        Therefore, it doesn't show up in the A/B denominator,
            //        but it does modify B and the RHS.

            // Bottom boundary coefficients and RHS for L decomp
            //===================================================
            Real rhoAlpha_lo, rhoAlpha_hi;
            Real met_h_zeta_lo, met_h_zeta_hi;
            Real a_tmp, b_tmp, c_tmp, inv_b2_tmp;
            {
                getRhoAlpha(i, j, klo, rhoAlpha_lo, rhoAlpha_hi,
                            cell_data, mu_turb, d_alpha_eff, d_eddy_diff_idz,
                            prim_index, prim_scal_index, l_consA, l_turb);

                met_h_zeta_lo = Compute_h_zeta_AtKface(i,j,klo  ,cellSizeInv,z_nd);
                met_h_zeta_hi = Compute_h_zeta_AtKface(i,j,klo+1,cellSizeInv,z_nd);

                a_tmp      = 0.;
                c_tmp      = -Fact * rhoAlpha_hi * dz_inv / met_h_zeta_hi;
                b_tmp      = detJ(i,j,klo) * cell_data(i,j,klo,Rho_comp) - a_tmp - c_tmp;
                inv_b2_tmp = 1.;

                RHS_a(i,j,klo) = detJ(i,j,klo) * cell_data(i,j,klo,n); // NOTE: this is rho*theta; solution is theta
                if (use_SurfLayer) {
                    RHS_a(i,j,klo) +=  Fact * hfx_z(i,j,klo); // NOTE: hfx_z = -K*d_z(\theta)
                } else if (neumann_on_zlo) {
                    RHS_a(i,j,klo) += -Fact * rhoAlpha_lo * bc_neumann_vals[2]; // NOTE: N_val = d_z(\theta)
                }

                RHS_a(i,j,klo)    /= b_tmp;         // NOTE: this is now "rho"
                coeffG_a(i,j,klo)  = c_tmp / b_tmp; // NOTE: this is now "gamma"
            }

            // Build the coefficients and RHS for L decomp
            //===================================================
            for (int k(klo+1); k < khi; k++) {
                getRhoAlpha(i, j, k, rhoAlpha_lo, rhoAlpha_hi,
                            cell_data, mu_turb, d_alpha_eff, d_eddy_diff_idz,
                            prim_index, prim_scal_index, l_consA, l_turb);

                met_h_zeta_lo = Compute_h_zeta_AtKface(i,j,k  ,cellSizeInv,z_nd);
                met_h_zeta_hi = Compute_h_zeta_AtKface(i,j,k+1,cellSizeInv,z_nd);

                a_tmp      = -Fact * rhoAlpha_lo * dz_inv / met_h_zeta_lo;
                c_tmp      = -Fact * rhoAlpha_hi * dz_inv / met_h_zeta_hi;
                b_tmp      = detJ(i,j,k) * cell_data(i,j,k,Rho_comp) - a_tmp - c_tmp;
                inv_b2_tmp = 1. / (b_tmp - a_tmp * coeffG_a(i,j,k-1));

                RHS_a(i,j,k)    = detJ(i,j,k) * cell_data(i,j,k,n); // NOTE: this is rho*theta; solution is theta

                RHS_a(i,j,k)    = (RHS_a(i,j,k) - a_tmp * RHS_a(i,j,k-1)) * inv_b2_tmp; // NOTE: This is now "rho"
                coeffG_a(i,j,k) = c_tmp * inv_b2_tmp; // NOTE: this is now "gamma"
            } // k

            // Top boundary coefficients and RHS for L decomp
            //===================================================
            {
                getRhoAlpha(i, j, khi, rhoAlpha_lo, rhoAlpha_hi,
                            cell_data, mu_turb, d_alpha_eff, d_eddy_diff_idz,
                            prim_index, prim_scal_index, l_consA, l_turb);

                met_h_zeta_lo = Compute_h_zeta_AtKface(i,j,khi  ,cellSizeInv,z_nd);
                met_h_zeta_hi = Compute_h_zeta_AtKface(i,j,khi+1,cellSizeInv,z_nd);

                a_tmp      = -Fact * rhoAlpha_lo * dz_inv / met_h_zeta_hi;
                c_tmp      = 0.;
                b_tmp      = detJ(i,j,khi) * cell_data(i,j,khi,Rho_comp) - a_tmp - c_tmp;
                inv_b2_tmp = 1. / (b_tmp - a_tmp * coeffG_a(i,j,khi-1));

                RHS_a(i,j,khi) = detJ(i,j,khi) * cell_data(i,j,khi,n); // NOTE: this is rho*theta; solution is theta
                if (neumann_on_zhi) {
                    RHS_a(i,j,khi) -= -Fact * rhoAlpha_hi * bc_neumann_vals[5]; // NOTE: N_val = d_z(\theta)
                }

                // First solve
                soln_a(i,j,khi) = (RHS_a(i,j,khi) - a_tmp * RHS_a(i,j,khi-1)) * inv_b2_tmp;
            }

            // Back sweep the U decomp solution
            //===================================================
            for (int k(khi-1); k>=klo; --k) {
                soln_a(i,j,k) = RHS_a(i,j,k) - coeffG_a(i,j,k) * soln_a(i,j,k+1);
            }

            // Convert back to rho*theta
            //===================================================
            for (int k(klo); k<=khi; ++k) {
                cell_data(i,j,k,n) = cell_data(i,j,k,Rho_comp) * soln_a(i,j,k);
            }

#ifdef AMREX_USE_GPU
    });
#else
      } // i
    } // j
#endif
}

/**
 * Function for computing the implicit contribution to the vertical diffusion
 * of momentum, over terrain.
 *
 * This function (explicitly instantiated below) handles staggering in x, y, or
 * z through the template parameter, stagdir. NOTE: implicit diffusion of w has
 * remains an experimental feature and has *not* been tested yet with terrain.
 *
 * @param[in   ] bx cell-centered box to loop over
 * @param[in   ] domain box of the whole domain
 * @param[in   ] dt time step
 * @param[in   ] cell_data conserved cell-centered rho
 * @param[inout] face_data conserved momentum
 * @param[in   ] tau_corr stress contribution to momentum that will be corrected by the implicit solve
 * @param[in   ] z_nd nodal array of z
 * @param[in   ] detJ Jacobian determinant
 * @param[in   ] cellSizeInv inverse cell size array
 * @param[in   ] mu_turb turbulent viscosity
 * @param[in   ] solverChoice container of parameters
 * @param[in   ] bc_ptr container with boundary conditions
 * @param[in   ] use_SurfLayer whether we have turned on subgrid diffusion
 * @param[in   ] implicit_fac if 1 then fully implicit; if 0 then fully explicit
 */
template <int stagdir>
void
ImplicitDiffForMomLU_T (const Box& bx,
                        const Box& /*domain*/,
                        const int level,
                        const Real dt,
                        const Array4<const Real>& cell_data,
                        const Array4<      Real>& face_data,
                        const Array4<const Real>& tau,
                        const Array4<const Real>& tau_corr,
                        const Array4<const Real>& z_nd,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const BCRec* bc_ptr,
                        const bool use_SurfLayer,
                        const Real implicit_fac)
{
    BL_PROFILE_VAR("ImplicitDiffForMom_T()",ImplicitDiffForMom_T);

    // setup quantities for getRhoAlphaAtFaces()
    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];
    bool l_consA  = (dc.molec_diff_type == MolecDiffType::ConstantAlpha);
    bool l_turb   = tc.use_kturb;
    Real mu_eff = (l_consA) ? 2.0 * dc.dynamic_viscosity / dc.rho0_trans
                            : 2.0 * dc.dynamic_viscosity;

    // g(S*) coefficient
    // stagdir==0: tau_corr = 0.5 * du/dz * mu_tot
    // stagdir==1: tau_corr = 0.5 * dv/dz * mu_tot
    // stagdir==2: tau_corr =       dw/dz * mu_tot
    constexpr Real gfac = (stagdir == 2) ? 2.0/3.0 : 1.0;

    // offsets used to average to faces
    constexpr int ioff = (stagdir == 0) ? 1 : 0;
    constexpr int joff = (stagdir == 1) ? 1 : 0;

    // Box bounds
    int ilo = bx.smallEnd(0);
    int ihi = bx.bigEnd(0);
    int jlo = bx.smallEnd(1);
    int jhi = bx.bigEnd(1);
    int klo = bx.smallEnd(2);
    int khi = bx.bigEnd(2);
    amrex::ignore_unused(ilo, ihi, jlo, jhi);

    // Temporary FABs for tridiagonal solve (allocated on column)
    //   A[k] * x[k-1] + B[k] * x[k] + C[k+1] = RHS[k]
    amrex::FArrayBox RHS_fab, soln_fab, coeffG_fab;
           RHS_fab.resize(bx,1, amrex::The_Async_Arena());
          soln_fab.resize(bx,1, amrex::The_Async_Arena());
        coeffG_fab.resize(bx,1, amrex::The_Async_Arena());
    auto const& RHS_a        =        RHS_fab.array();
    auto const& soln_a       =       soln_fab.array();
    auto const& coeffG_a     =     coeffG_fab.array();

    Real dz_inv = cellSizeInv[2];

    int bc_comp = BCVars::xvel_bc + stagdir;
    bool ext_dir_on_zlo  = (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir ||
                            bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim);
    bool ext_dir_on_zhi  = (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir ||
                            bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim);
    bool foextrap_on_zhi = (bc_ptr[bc_comp].hi(2) == ERFBCType::foextrap);
    amrex::ignore_unused(foextrap_on_zhi);

    AMREX_ASSERT_WITH_MESSAGE(ext_dir_on_zlo || use_SurfLayer,
                              "Unexpected lower BC used with implicit vertical diffusion");
    AMREX_ASSERT_WITH_MESSAGE(foextrap_on_zhi || ext_dir_on_zhi,
                              "Unexpected upper BC used with implicit vertical diffusion");
    if (stagdir < 2 && (ext_dir_on_zlo || ext_dir_on_zhi)) {
        amrex::Warning("No-slip walls have not been fully tested");
    }

    Real Fact = implicit_fac * dt * dz_inv;

#ifdef AMREX_USE_GPU
    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
#else
    for (int j(jlo); j<=jhi; ++j) {
      for (int i(ilo); i<=ihi; ++i) {
#endif
          // Notes:
          //
          // - In DiffusionSrcForMom (e.g., for x-mom)
          //
          //     Real diffContrib = ...
          //                      + (tau13(i,j,k+1) - tau13(i,j,k)) / dzinv
          //     rho_u_rhs(i,j,k) -= diffContrib;  // note the negative sign
          //
          // - We need to scale the explicit _part_ of `tau13` (for x-mom) by (1 - implicit_fac)
          //   The part that needs to be scaled is stored in `tau_corr`.
          //   E.g., tau13 = 0.5 * (du/dz + dw/dx)
          //         tau13_corr = 0.5 * du/dz
          //
          // - The momentum (`face_data`) was set to `S_old + S_rhs * dt`
          //   prior to including "ERF_Implicit.H". Recall that S_rhs includes
          //   sources from advection and other forcings, not just diffusion.
          //
          // - To correct momentum, we need to subtract `implicit_fac * diffContrib_corr`
          //   from S_rhs to recover `(1 - implicit_fac) * diffContrib_corr`,
          //   where `diffContrib_corr = -d(tau_corr)/dz`. The negative sign
          //   comes from our convention for the RHS diffusion source.
          //
          //   Subtracting a negative gives the += below; multiply by dt to
          //   get the intermediate momentum on the RHS of the tridiagonal
          //   system.
          //
          // - With a surface_layer BC, tau13/23 holds the vertical flux -d_z(k*u_i)
          //   directly. We must use tau at klo (not tau_corr) with SL BCs.
          //
          // - The detJ for divergence was multiplied through.
          //   Therefore, it doesn't show up in the A/B denominator,
          //   but it does modify B and the RHS.
          //
          // - Finally, the terms ~ RHS += (tau_corr_hi - tau_corr_lo) / dz (below)
          //   essentially undo the explicit diffusion update that will be
          //   handled here implicitly.

          // Bottom boundary coefficients and RHS for L decomp
          //===================================================
          Real rhoface, rhoAlpha_lo, rhoAlpha_hi;
          Real detJface, met_h_zeta_lo, met_h_zeta_hi;
          Real a_tmp, b_tmp, c_tmp, inv_b2_tmp;
          {
              detJface = 0.5 * (detJ(i,j,klo) + detJ(i-ioff,j-joff,klo));
              rhoface  = 0.5 * (cell_data(i,j,klo,Rho_comp) + cell_data(i-ioff,j-joff,klo,Rho_comp));
              getRhoAlphaForFaces(i, j, klo, ioff, joff, rhoAlpha_lo, rhoAlpha_hi,
                                  cell_data, mu_turb, mu_eff,
                                  l_consA, l_turb);

              met_h_zeta_lo = 0.5 * ( Compute_h_zeta_AtKface(i     ,j     ,klo  ,cellSizeInv,z_nd)
                                    + Compute_h_zeta_AtKface(i-ioff,j-joff,klo  ,cellSizeInv,z_nd) );
              met_h_zeta_hi = 0.5 * ( Compute_h_zeta_AtKface(i     ,j     ,klo+1,cellSizeInv,z_nd)
                                    + Compute_h_zeta_AtKface(i-ioff,j-joff,klo+1,cellSizeInv,z_nd) );

              a_tmp = 0.;
              c_tmp = -Fact * gfac * rhoAlpha_hi * dz_inv / met_h_zeta_hi;

              RHS_a(i,j,klo) = detJface * face_data(i,j,klo); // NOTE: this is momenta; solution is velocity

              // BCs: Dirichlet (u_i = val), slip wall (w = 0), or surface layer (w = 0)
              if (ext_dir_on_zlo) {
                  RHS_a(i,j,klo) += Fact * gfac * (tau_corr(i,j,klo+1) - tau_corr(i,j,klo));
                  if (stagdir==2) {
                      c_tmp = 0.;
                      RHS_a(i,j,klo) = 0.;
                  } else {
                      a_tmp = -2.0 * Fact * rhoAlpha_lo * dz_inv / met_h_zeta_lo;
                      RHS_a(i,j,klo) += 2.0 * rhoAlpha_lo * face_data(i,j,klo-1) * dz_inv * dz_inv / met_h_zeta_lo;
                  }
              } else if (use_SurfLayer) {
                  // NOTE: tau = -mu*d_z(u_i) w/ SL
                  RHS_a(i,j,klo) += Fact * gfac * (tau_corr(i,j,klo+1) - tau(i,j,klo));
                  RHS_a(i,j,klo) += Fact * tau(i,j,klo);
              }

              b_tmp      = detJface * rhoface - a_tmp - c_tmp;
              inv_b2_tmp = 1.;

              RHS_a(i,j,klo)    /= b_tmp;         // NOTE: this is now "rho"
              coeffG_a(i,j,klo)  = c_tmp / b_tmp; // NOTE: this is now "gamma"
          }

          // Build the coefficients and RHS for L decomp
          //===================================================
          for (int k(klo+1); k < khi; k++) {
              detJface = 0.5 * (detJ(i,j,k) + detJ(i-ioff,j-joff,k));
              rhoface  = 0.5 * (cell_data(i,j,k,Rho_comp) + cell_data(i-ioff,j-joff,k,Rho_comp));
              getRhoAlphaForFaces(i, j, k, ioff, joff, rhoAlpha_lo, rhoAlpha_hi,
                                  cell_data, mu_turb, mu_eff,
                                  l_consA, l_turb);

              met_h_zeta_lo = 0.5 * ( Compute_h_zeta_AtKface(i     ,j     ,k  ,cellSizeInv,z_nd)
                                    + Compute_h_zeta_AtKface(i-ioff,j-joff,k  ,cellSizeInv,z_nd) );
              met_h_zeta_hi = 0.5 * ( Compute_h_zeta_AtKface(i     ,j     ,k+1,cellSizeInv,z_nd)
                                    + Compute_h_zeta_AtKface(i-ioff,j-joff,k+1,cellSizeInv,z_nd) );

              a_tmp      = -Fact * rhoAlpha_lo * dz_inv / met_h_zeta_lo;
              c_tmp      = -Fact * rhoAlpha_hi * dz_inv / met_h_zeta_hi;
              b_tmp      = detJface * rhoface - a_tmp - c_tmp;
              inv_b2_tmp = 1. / (b_tmp - a_tmp * coeffG_a(i,j,k-1));

              RHS_a(i,j,k)    = detJface * face_data(i,j,k); // NOTE: this is momenta; solution is velocity
              RHS_a(i,j,k)   += Fact * gfac * (tau_corr(i,j,k+1) - tau_corr(i,j,k));

              RHS_a(i,j,k)    = (RHS_a(i,j,k) - a_tmp * RHS_a(i,j,k-1)) * inv_b2_tmp; // NOTE: This is now "rho"
              coeffG_a(i,j,k) = c_tmp * inv_b2_tmp; // NOTE: this is now "gamma"
          } // k

          // Top boundary coefficients and RHS for L decomp
          //===================================================
          {
              detJface = 0.5 * (detJ(i,j,khi) + detJ(i-ioff,j-joff,khi));
              rhoface  = 0.5 * (cell_data(i,j,khi,Rho_comp) + cell_data(i-ioff,j-joff,khi,Rho_comp));
              getRhoAlphaForFaces(i, j, khi, ioff, joff, rhoAlpha_lo, rhoAlpha_hi,
                                  cell_data, mu_turb, mu_eff,
                                  l_consA, l_turb);

              met_h_zeta_lo = 0.5 * ( Compute_h_zeta_AtKface(i     ,j     ,khi  ,cellSizeInv,z_nd)
                                    + Compute_h_zeta_AtKface(i-ioff,j-joff,khi  ,cellSizeInv,z_nd) );
              met_h_zeta_hi = 0.5 * ( Compute_h_zeta_AtKface(i     ,j     ,khi+1,cellSizeInv,z_nd)
                                    + Compute_h_zeta_AtKface(i-ioff,j-joff,khi+1,cellSizeInv,z_nd) );

              a_tmp = -Fact * gfac * rhoAlpha_lo * dz_inv / met_h_zeta_lo;
              c_tmp = 0.;

              RHS_a(i,j,khi)  = detJface * face_data(i,j,khi); // NOTE: this is momenta; solution is velocity
              RHS_a(i,j,khi) += Fact * gfac * (tau_corr(i,j,khi+1) - tau_corr(i,j,khi));

              // BCs: Dirichlet (u_i = val), slip wall (w = 0)
              if (ext_dir_on_zhi) {
                  if (stagdir==2) {
                      a_tmp = 0.;
                      RHS_a(i,j,khi) = 0.;
                  } else {
                      c_tmp = -2.0 * Fact * rhoAlpha_hi * dz_inv / met_h_zeta_hi;
                      RHS_a(i,j,khi) += 2.0 * rhoAlpha_hi * face_data(i,j,khi+1) * dz_inv * dz_inv / met_h_zeta_hi;
                  }
              }

              b_tmp      = detJface * rhoface - a_tmp - c_tmp;
              inv_b2_tmp = 1. / (b_tmp - a_tmp * coeffG_a(i,j,khi-1));

              // First solve
              soln_a(i,j,khi) = (RHS_a(i,j,khi) - a_tmp * RHS_a(i,j,khi-1)) * inv_b2_tmp;
          }

          // Back sweep the U decomp solution
          //===================================================
          for (int k(khi-1); k>=klo; --k) {
              soln_a(i,j,k) = RHS_a(i,j,k) - coeffG_a(i,j,k) * soln_a(i,j,k+1);
          }

          // Convert back to momenta
          //===================================================
          for (int k(klo); k<=khi; ++k) {
              rhoface = 0.5 * (cell_data(i,j,k,Rho_comp) + cell_data(i-ioff,j-joff,k,Rho_comp));
              face_data(i,j,k) = rhoface * soln_a(i,j,k);
          }

#ifdef AMREX_USE_GPU
    });
#else
      } // i
    } // j
#endif
}

#define INSTANTIATE_IMPLICIT_DIFF_FOR_MOM_LU(STAGDIR) \
    template void ImplicitDiffForMomLU_T<STAGDIR> ( \
        const Box&, \
        const Box&, \
        const int, \
        const Real, \
        const Array4<const Real>&, \
        const Array4<      Real>&, \
        const Array4<const Real>&, \
        const Array4<const Real>&, \
        const Array4<const Real>&, \
        const Array4<const Real>&, \
        const GpuArray<Real, AMREX_SPACEDIM>&, \
        const Array4<const Real>&, \
        const SolverChoice&, \
        const BCRec*, \
        const bool, \
        const Real);
INSTANTIATE_IMPLICIT_DIFF_FOR_MOM_LU(0)
INSTANTIATE_IMPLICIT_DIFF_FOR_MOM_LU(1)
INSTANTIATE_IMPLICIT_DIFF_FOR_MOM_LU(2)
#undef INSTANTIATE_IMPLICIT_DIFF_FOR_MOM
