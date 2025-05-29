#include <AMReX.H>
#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>

using namespace amrex;

/**
 * Function for computing the momentum RHS for diffusion operator with terrain.
 *
 * @param[in]  bxx nodal x box for x-mom
 * @param[in]  bxy nodal y box for y-mom
 * @param[in]  bxz nodal z box for z-mom
 * @param[out] rho_u_rhs RHS for x-mom
 * @param[out] rho_v_rhs RHS for y-mom
 * @param[out] rho_w_rhs RHS for z-mom
 * @param[in]  tau11 11 stress
 * @param[in]  tau22 22 stress
 * @param[in]  tau33 33 stress
 * @param[in]  tau12 12 stress
 * @param[in]  tau13 13 stress
 * @param[in]  tau21 21 stress
 * @param[in]  tau23 23 stress
 * @param[in]  tau31 31 stress
 * @param[in]  tau32 32 stress
 * @param[in]  detJ Jacobian determinant
 * @param[in]  differChoice container with diffusion parameters
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_m map factor at cell center
 */
void
DiffusionSrcForMom_T (const Box& bxx, const Box& bxy , const Box& bxz,
                      const Array4<Real>& rho_u_rhs  ,
                      const Array4<Real>& rho_v_rhs  ,
                      const Array4<Real>& rho_w_rhs  ,
                      const Array4<const Real>& tau11,
                      const Array4<const Real>& tau22,
                      const Array4<const Real>& tau33,
                      const Array4<const Real>& tau12, const Array4<const Real>& tau21,
                      const Array4<const Real>& tau13, const Array4<const Real>& tau31,
                      const Array4<const Real>& tau23, const Array4<const Real>& tau32,
                      const Array4<const Real>& detJ ,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                      const Array4<const Real>& mf_mx,
                      const Array4<const Real>& mf_ux,
                      const Array4<const Real>& mf_vx,
                      const Array4<const Real>& mf_my,
                      const Array4<const Real>& mf_uy,
                      const Array4<const Real>& mf_vy)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_T()",DiffusionSrcForMom_T);

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];

    ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // Inv Jacobian
        Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

        Real diffContrib  = ( (tau11(i  , j  , k  ) - tau11(i-1, j  ,k  )) * dxinv * mfsq  // Contribution to x-mom eqn from diffusive flux in x-dir
                            + (tau12(i  , j+1, k  ) - tau12(i  , j  ,k  )) * dyinv * mfsq  // Contribution to x-mom eqn from diffusive flux in y-dir
                            + (tau13(i  , j  , k+1) - tau13(i  , j  ,k  )) * dzinv );     // Contribution to x-mom eqn from diffusive flux in z-dir;
        diffContrib      /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k));
        rho_u_rhs(i,j,k) -= diffContrib;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // Inv Jacobian
        Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

        Real diffContrib  = ( (tau21(i+1, j  , k  ) - tau21(i  , j  , k  )) * dxinv * mfsq  // Contribution to y-mom eqn from diffusive flux in x-dir
                            + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  )) * dyinv * mfsq // Contribution to y-mom eqn from diffusive flux in y-dir
                            + (tau23(i  , j  , k+1) - tau23(i  , j  , k  )) * dzinv );     // Contribution to y-mom eqn from diffusive flux in z-dir;
        diffContrib      /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k));
        rho_v_rhs(i,j,k) -= diffContrib;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // Inv Jacobian
        Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

        Real diffContrib  = ( (tau31(i+1, j  , k  ) - tau31(i  , j  , k  )) * dxinv * mfsq  // Contribution to z-mom eqn from diffusive flux in x-dir
                            + (tau32(i  , j+1, k  ) - tau32(i  , j  , k  )) * dyinv * mfsq  // Contribution to z-mom eqn from diffusive flux in y-dir
                            + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1)) * dzinv );     // Contribution to z-mom eqn from diffusive flux in z-dir;
        diffContrib      /= 0.5*(detJ(i,j,k) + detJ(i,j,k-1));
        rho_w_rhs(i,j,k) -= diffContrib;
    });
}


void ImplicitXmomDiffusion (int& /*lev*/,
                            amrex::Real& dt,
                            amrex::Real mu,
                            const amrex::Geometry& /*geom*/,
                            amrex::MultiFab& xmom,
                            std::unique_ptr<amrex::MultiFab>& z_phys_nd)
{
  amrex::Print() << "IMPLICIT DIFFUSION: Starting Tridiagonal solve...";

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for ( amrex::MFIter mfi(xmom); mfi.isValid(); ++mfi)
  {
    // Boxes
    amrex::Box vbx = mfi.validbox();
    amrex::Box gbx = amrex::grow(vbx,2,1);
    amrex::Box b2d = gbx;
    b2d.setRange(2,0);

    // State data
    auto  xm_arr = xmom.array(mfi);
    auto znd_arr = z_phys_nd->array(mfi);

    // Box bounds
    int ilo = gbx.smallEnd(0);
    int ihi = gbx.bigEnd(0);
    int jlo = gbx.smallEnd(1);
    int jhi = gbx.bigEnd(1);
    int klo = gbx.smallEnd(2);
    int khi = gbx.bigEnd(2);

    // Temporary FABs for tridiagonal solve (allocated on column)
    amrex::FArrayBox RHS_fab, soln_fab, coeffA_fab, inv_coeffB_fab, coeffC_fab;
           RHS_fab.resize(gbx,1, amrex::The_Async_Arena());
          soln_fab.resize(gbx,1, amrex::The_Async_Arena());
        coeffA_fab.resize(gbx,1, amrex::The_Async_Arena());
    inv_coeffB_fab.resize(gbx,1, amrex::The_Async_Arena());
        coeffC_fab.resize(gbx,1, amrex::The_Async_Arena());
    auto const& RHS_a        =        RHS_fab.array();
    auto const& soln_a       =       soln_fab.array();
    auto const& coeffA_a     =     coeffA_fab.array();
    auto const& inv_coeffB_a = inv_coeffB_fab.array();
    auto const& coeffC_a     =     coeffC_fab.array();

    for (int lj(jlo); lj<=jhi; ++lj) {
      for (int li(ilo); li<=ihi; ++li) {

        // Dirichlet top/bottom
        RHS_a(li,lj,klo) = 0.;
        RHS_a(li,lj,khi) = 0.;
        coeffA_a(li,lj,klo) = 0.0;
        coeffA_a(li,lj,khi) = 0.0;
        coeffC_a(li,lj,klo) = 0.0;
        coeffC_a(li,lj,khi) = 0.0;
        inv_coeffB_a(li,lj,klo) = 1.0;
        inv_coeffB_a(li,lj,khi) = 1.0;


        // Build the coefficients and RHS
        for (int lk(klo+1); lk<khi; lk++) {
          amrex::Real z0 = (lk==klo+1) ? 0.5 *(znd_arr(li,lj,lk  ) + znd_arr(li,lj+1,lk  )) :
                                         0.25*(znd_arr(li,lj,lk  ) + znd_arr(li,lj+1,lk  )
                                              +znd_arr(li,lj,lk-1) + znd_arr(li,lj+1,lk-1));
          amrex::Real z1 = 0.25*(znd_arr(li,lj,lk  ) + znd_arr(li,lj+1,lk  )
                                +znd_arr(li,lj,lk+1) + znd_arr(li,lj+1,lk+1));
          amrex::Real z2 = (lk==khi-1) ? 0.5 *(znd_arr(li,lj,lk+1) + znd_arr(li,lj+1,lk+1)) :
                                         0.25*(znd_arr(li,lj,lk+1) + znd_arr(li,lj+1,lk+1)
                                              +znd_arr(li,lj,lk+2) + znd_arr(li,lj+1,lk+2));
          amrex::Real dz_lo = z1 - z0;
          amrex::Real dz_hi = z2 - z1;
          amrex::Real dz    = z2 - z0;

          coeffA_a(li,lj,lk)     = -mu * dt / (dz_lo*dz);
          coeffC_a(li,lj,lk)     = -mu * dt / (dz_hi*dz);
          amrex::Real b = 1.0 + (mu * dt) / (dz_lo*dz_hi);
          inv_coeffB_a(li,lj,lk) = 1.0 / (b - coeffA_a(li,lj,lk)*(coeffC_a(li,lj,lk)/b));
          RHS_a(li,lj,lk)        = xm_arr(li,lj,lk);
        }

        // Tridiagonal solve
        soln_a(li,lj,klo) = RHS_a(li,lj,klo) * inv_coeffB_a(li,lj,klo);
        for (int lk(klo+1); lk<=khi; ++lk) {
          soln_a(li,lj,lk) = (RHS_a(li,lj,lk)-coeffA_a(li,lj,lk)*soln_a(li,lj,lk-1)) * inv_coeffB_a(li,lj,lk);
        }
        for (int lk(khi-1); lk>=klo; --lk) {
          soln_a(li,lj,lk) -= ( coeffC_a(li,lj,lk) * inv_coeffB_a(li,lj,lk) ) * soln_a(li,lj,lk+1);
        }

        // Write soln back to state
        for (int lk(klo+1); lk<khi; ++lk) {
          xm_arr(li,lj,lk) = soln_a(li,lj,lk);
        }
      } // li
    } // lj
  } // mfi
  amrex::Print() << "DONE\n";
}
