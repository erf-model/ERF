#include <ERF_Diffusion.H>

using namespace amrex;

/**
 * Function for computing the stress with constant viscosity for EB.
 *
 * @param[in] bxcc cell center box for tau_ii
 * @param[in] tbxxy nodal xy box for tau_12
 * @param[in] tbxxz nodal xz box for tau_13
 * @param[in] tbxyz nodal yz box for tau_23
 * @param[in] mu_eff constant molecular viscosity
 * @param[in] cell_data to access rho if ConstantAlpha
 * @param[in,out] tau11 11 strain -> stress
 * @param[in,out] tau22 22 strain -> stress
 * @param[in,out] tau33 33 strain -> stress
 * @param[in,out] tau12 12 strain -> stress
 * @param[in,out] tau13 13 strain -> stress
 * @param[in,out] tau23 23 strain -> stress
 * @param[in] er_arr expansion rate
 * @param[in,out] tau13i contribution to stress from du/dz
 * @param[in,out] tau23i contribution to stress from dv/dz
 * @param[in,out] tau33i contribution to stress from dw/dz
 */
void
ComputeStressConsVisc_EB (Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Real mu_eff,
                         const Array4<const Real>& cell_data,
                         Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                         Array4<Real>& tau12, Array4<Real>& tau13, Array4<Real>& tau23,
                         const Array4<const Real>& er_arr,
                         Array4<const Real>& vfrac,
                         Array4<Real>& tau13i,
                         Array4<Real>& tau23i,
                         Array4<Real>& tau33i)
{
    Real OneThird   = (1./3.);

    // NOTE: mu_eff includes factor of 2

    if (cell_data)
    // constant alpha (stored in mu_eff)
    {
        // Cell centered strains
        ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha  = cell_data(i, j, k, Rho_comp) * mu_eff;
            if (tau33i) tau33i(i,j,k) = -rhoAlpha * tau33(i,j,k);
            tau11(i,j,k) = -rhoAlpha * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
            tau22(i,j,k) = -rhoAlpha * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
            tau33(i,j,k) = -rhoAlpha * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
        });

        // Off-diagonal strains
        ParallelFor(tbxxy,tbxxz,tbxyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i-1,j,k) + vfrac(i,j-1,k) + vfrac(i-1,j-1,k);
            Real rho_bar = 0.0;
            if (vol_sum > 1.e-16) {
                rho_bar = ( vfrac(i-1,j,k) * cell_data(i-1, j  , k, Rho_comp)
                            + vfrac(i,j,k) * cell_data(i, j  , k, Rho_comp)
                            + vfrac(i-1,j-1,k) * cell_data(i-1, j-1, k, Rho_comp)
                            + vfrac(i,j-1,k) * cell_data(i, j-1, k, Rho_comp) ) / vol_sum;
            } else {
                rho_bar = 0.25*( cell_data(i-1, j  , k, Rho_comp) + vfrac(i,j,k) * cell_data(i, j  , k, Rho_comp)
                            + cell_data(i-1, j-1, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            }
            tau12(i,j,k) *= -rho_bar * mu_eff;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i-1,j,k) + vfrac(i,j,k-1) + vfrac(i-1,j,k-1);
            Real rho_bar = 0.0;
            if (vol_sum > 1.e-16) {
                rho_bar = ( vfrac(i-1,j,k) * cell_data(i-1, j, k  , Rho_comp)
                            + vfrac(i,j,k) * cell_data(i, j, k  , Rho_comp)
                            + vfrac(i-1,j,k-1) * cell_data(i-1, j, k-1, Rho_comp)
                            + vfrac(i,j,k-1) * cell_data(i, j, k-1, Rho_comp) )/ vol_sum;
            } else {
                rho_bar = 0.25*( cell_data(i-1, j, k  , Rho_comp) + cell_data(i, j, k  , Rho_comp)
                                + cell_data(i-1, j, k-1, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            }
            tau13(i,j,k) *= -rho_bar * mu_eff;

            if (tau13i) tau13i(i,j,k) *= -rho_bar * mu_eff;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i,j-1,k) + vfrac(i,j,k-1) + vfrac(i,j-1,k-1);
            Real rho_bar = 0.0;
            if (vol_sum > 1.e-16) {
                rho_bar = ( vfrac(i,j-1,k) * cell_data(i, j-1, k  , Rho_comp)
                            + vfrac(i,j,k) * cell_data(i, j, k  , Rho_comp)
                            + vfrac(i,j-1,k-1) * cell_data(i, j-1, k-1, Rho_comp)
                            + vfrac(i,j,k-1) * cell_data(i, j, k-1, Rho_comp) ) / vol_sum;
            } else {
                rho_bar = 0.25*( cell_data(i, j-1, k  , Rho_comp) + cell_data(i, j, k  , Rho_comp)
                                + cell_data(i, j-1, k-1, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            }
            tau23(i,j,k) *= -rho_bar * mu_eff;

            if (tau23i) tau23i(i,j,k) *= -rho_bar * mu_eff;
        });
    }
    else
    // constant mu_eff
    {
        // Cell centered strains
        ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (tau33i) tau33i(i,j,k) = -mu_eff * tau33(i,j,k);
            tau11(i,j,k) = -mu_eff * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
            tau22(i,j,k) = -mu_eff * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
            tau33(i,j,k) = -mu_eff * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
        });

        // Off-diagonal strains
        ParallelFor(tbxxy,tbxxz,tbxyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau12(i,j,k) *= -mu_eff;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau13(i,j,k) *= -mu_eff;

            if (tau13i) tau13i(i,j,k) *= -mu_eff;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau23(i,j,k) *= -mu_eff;

            if (tau23i) tau23i(i,j,k) *= -mu_eff;
        });
    }
}