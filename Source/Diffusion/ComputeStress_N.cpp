#include <Diffusion.H>

using namespace amrex;

/**
 * Function for computing the stress with constant viscosity and without terrain.
 *
 * @param[in] bxcc cell center box for tau_ii
 * @param[in] tbxxy nodal xy box for tau_12
 * @param[in] tbxxz nodal xz box for tau_13
 * @param[in] tbxyz nodal yz box for tau_23
 * @param[in] mu_eff constant molecular viscosity
 * @param[in,out] tau11 11 strain -> stress
 * @param[in,out] tau22 22 strain -> stress
 * @param[in,out] tau33 33 strain -> stress
 * @param[in,out] tau12 12 strain -> stress
 * @param[in,out] tau13 13 strain -> stress
 * @param[in,out] tau23 23 strain -> stress
 * @param[in] er_arr expansion rate
 */
void
ComputeStressConsVisc_N(Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Real mu_eff,
                        Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                        Array4<Real>& tau12, Array4<Real>& tau13, Array4<Real>& tau23,
                        const Array4<const Real>& er_arr)
{
    Real OneThird   = (1./3.);

    // Cell centered strains
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        tau11(i,j,k) = -mu_eff * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
        tau22(i,j,k) = -mu_eff * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
        tau33(i,j,k) = -mu_eff * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau12(i,j,k) *= -mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau13(i,j,k) *= -mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau23(i,j,k) *= -mu_eff;
    });

}

/**
 * Function for computing the stress with variable viscosity and without terrain.
 *
 * @param[in] bxcc cell center box for tau_ii
 * @param[in] tbxxy nodal xy box for tau_12
 * @param[in] tbxxz nodal xz box for tau_13
 * @param[in] tbxyz nodal yz box for tau_23
 * @param[in] mu_eff constant molecular viscosity
 * @param[in] mu_turb variable turbulent viscosity
 * @param[in,out] tau11 11 strain -> stress
 * @param[in,out] tau22 22 strain -> stress
 * @param[in,out] tau33 33 strain -> stress
 * @param[in,out] tau12 12 strain -> stress
 * @param[in,out] tau13 13 strain -> stress
 * @param[in,out] tau23 23 strain -> stress
 * @param[in] er_arr expansion rate
 */
void
ComputeStressVarVisc_N(Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Real mu_eff,
                       const Array4<const Real>& mu_turb,
                       Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                       Array4<Real>& tau12, Array4<Real>& tau13, Array4<Real>& tau23,
                       const Array4<const Real>& er_arr)
{
    Real OneThird   = (1./3.);

    // Cell centered strains
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_11 = mu_eff + 2.0 * mu_turb(i, j, k, EddyDiff::Mom_h);
        Real mu_22 = mu_11;
        Real mu_33 = mu_eff + 2.0 * mu_turb(i, j, k, EddyDiff::Mom_v);
        tau11(i,j,k) = -mu_11 * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
        tau22(i,j,k) = -mu_22 * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
        tau33(i,j,k) = -mu_33 * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_bar = 0.25*( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i, j  , k, EddyDiff::Mom_h)
                           + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i, j-1, k, EddyDiff::Mom_h) );
        Real mu_12  = mu_eff + 2.0*mu_bar;
        tau12(i,j,k) *= -mu_12;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_bar = 0.25*( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                           + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
        Real mu_13  = mu_eff + 2.0*mu_bar;
        tau13(i,j,k) *= -mu_13;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_bar = 0.25*( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                           + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
        Real mu_23  = mu_eff + 2.0*mu_bar;
        tau23(i,j,k) *= -mu_23;
    });
}
