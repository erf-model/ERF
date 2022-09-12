#include <AMReX.H>
#include <DiffusionSrcForMom_N.H>
#include <IndexDefines.H>

using namespace amrex;

void
DiffusionSrcForMom_N (const Box& bxx, const Box& bxy, const Box& bxz, const Box& domain,
                      const Array4<      Real>& rho_u_rhs, const Array4<      Real>& rho_v_rhs,
                      const Array4<      Real>& rho_w_rhs,
                      const Array4<const Real>& u        , const Array4<const Real>& v,
                      const Array4<const Real>& w        , const Array4<const Real>& K_turb,
                      const Array4<const Real>& cell_data, const Array4<const Real>& er_arr,
                      const SolverChoice& solverChoice   , const BCRec* bc_ptr,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_N()",DiffusionSrcForMom_N);

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {

    const int l_use_terrain   = solverChoice.use_terrain;

    int domhi_z = domain.bigEnd(2);

    AMREX_ALWAYS_ASSERT(!l_use_terrain);

    // *********************************************************************
    // Define diffusive updates in the RHS of {x, y, z}-momentum equations
    // *********************************************************************
#if 1
    amrex::ParallelFor(bxx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        rho_u_rhs(i, j, k) += DiffusionSrcForXMom(i, j, k, u, v, w, cell_data,
                                                  dxInv, K_turb, solverChoice,
                                                  domain, bc_ptr, er_arr);
    });
    amrex::ParallelFor(bxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        rho_v_rhs(i, j, k) += DiffusionSrcForYMom(i, j, k, u, v, w, cell_data,
                                                  dxInv, K_turb, solverChoice,
                                                  domain, bc_ptr, er_arr);
    });
    amrex::ParallelFor(bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        int k_diff = (k == domhi_z+1) ? domhi_z : k;
        Real temp = DiffusionSrcForZMom(i, j, k_diff, u, v, w, cell_data,
                                        dxInv, K_turb, solverChoice,
                                        domain, bc_ptr, er_arr);
        rho_w_rhs(i, j, k) += temp;
    });
# else
    #else
    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];
    Real OneThird = (1./3.);

    amrex::ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        amrex::Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;

        tau11Prev = (u(i  , j, k) - u(i-1, j, k))*dxinv - OneThird*er_arr(i-1, j, k);
        tau11Next = (u(i+1, j, k) - u(i  , j, k))*dxinv - OneThird*er_arr(i  , j, k);

        // SAVE
        tau12Prev = (u(i, j  , k) - u(i, j-1, k))*dyinv + (v(i, j  , k) - v(i-1, j  , k)) * dxinv;
        tau12Next = (u(i, j+1, k) - u(i, j  , k))*dyinv + (v(i, j+1, k) - v(i-1, j+1, k)) * dxinv;

        // SAVE
        tau13Prev = (u(i, j, k  ) - u(i, j, k-1))*dzinv + (w(i, j, k  ) - w(i-1, j, k  )) * dxinv;
        tau13Next = (u(i, j, k+1) - u(i, j, k  ))*dzinv + (w(i, j, k+1) - w(i-1, j, k+1)) * dxinv;

        rho_u_rhs(i, j, k) +=        (tau11Next - tau11Prev) * dxinv  // Contribution to x-mom eqn from diffusive flux in x-dir
                             + 0.5 * (tau12Next - tau12Prev) * dyinv  // Contribution to x-mom eqn from diffusive flux in y-dir
                             + 0.5 * (tau13Next - tau13Prev) * dzinv; // Contribution to x-mom eqn from diffusive flux in z-dir;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        amrex::Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;

        // RE-USE
        tau21Prev = (u(i  , j, k) - u(i  , j-1, k))*dyinv + (v(i  , j, k) - v(i-1, j, k)) * dxinv;
        tau21Next = (u(i+1, j, k) - u(i+1, j-1, k))*dyinv + (v(i+1, j, k) - v(i  , j, k)) * dxinv;

        tau22Prev = (v(i, j  , k) - v(i, j-1, k))*dyinv - OneThird*er_arr(i, j-1, k);
        tau22Next = (v(i, j+1, k) - v(i, j  , k))*dyinv - OneThird*er_arr(i, j  , k);

        // SAVE
        tau23Prev = (v(i, j, k  ) - v(i, j, k-1))*dzinv + (w(i, j, k  ) - w(i, j-1, k  )) * dyinv;
        tau23Next = (v(i, j, k+1) - v(i, j, k  ))*dzinv + (w(i, j, k+1) - w(i, j-1, k+1)) * dyinv;

        rho_v_rhs(i, j, k) +=  0.5 * (tau21Next - tau21Prev) * dxinv  // Contribution to y-mom eqn from diffusive flux in x-dir
                             +       (tau22Next - tau22Prev) * dyinv  // Contribution to y-mom eqn from diffusive flux in y-dir
                             + 0.5 * (tau23Next - tau23Prev) * dzinv; // Contribution to y-mom eqn from diffusive flux in z-dir;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        amrex::Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev;

        // RE-USE
        tau31Prev = (u(i  , j, k) - u(i  , j, k-1))*dzinv + (w(i  , j, k) - w(i-1, j, k)) * dxinv;
        tau31Next = (u(i+1, j, k) - u(i+1, j, k-1))*dzinv + (w(i+1, j, k) - w(i  , j, k)) * dxinv;

        // RE-USE
        tau32Prev = (v(i, j  , k) - v(i, j  , k-1))*dzinv + (w(i, j  , k) - w(i, j-1, k)) * dyinv;
        tau32Next = (v(i, j+1, k) - v(i, j+1, k-1))*dzinv + (w(i, j+1, k) - w(i, j  , k)) * dyinv;

        tau33Prev = (w(i, j, k  ) - w(i, j, k-1))*dzinv - OneThird*er_arr(i, j, k-1);
        tau33Next = (w(i, j, k+1) - w(i, j, k  ))*dzinv - OneThird*er_arr(i, j, k  );

        rho_w_rhs(i, j, k) +=  0.5 * (tau31Next - tau31Prev) * dxinv  // Contribution to z-mom eqn from diffusive flux in x-dir
                             + 0.5 * (tau32Next - tau32Prev) * dyinv  // Contribution to z-mom eqn from diffusive flux in y-dir
                             +       (tau33Next - tau33Prev) * dzinv; // Contribution to z-mom eqn from diffusive flux in z-dir;
    });
#endif
    }
}
