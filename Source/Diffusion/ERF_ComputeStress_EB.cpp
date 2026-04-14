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
    Real OneThird   = (one/three);

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
            Real rho_bar = zero;
            if (vol_sum > Real(1.e-16)) {
                rho_bar = ( vfrac(i-1,j,k) * cell_data(i-1, j  , k, Rho_comp)
                            + vfrac(i,j,k) * cell_data(i, j  , k, Rho_comp)
                            + vfrac(i-1,j-1,k) * cell_data(i-1, j-1, k, Rho_comp)
                            + vfrac(i,j-1,k) * cell_data(i, j-1, k, Rho_comp) ) / vol_sum;
            } else {
                rho_bar = fourth*( cell_data(i-1, j  , k, Rho_comp) + cell_data(i, j  , k, Rho_comp)
                               + cell_data(i-1, j-1, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            }
            tau12(i,j,k) *= -rho_bar * mu_eff;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i-1,j,k) + vfrac(i,j,k-1) + vfrac(i-1,j,k-1);
            Real rho_bar = zero;
            if (vol_sum > Real(1.e-16)) {
                rho_bar = ( vfrac(i-1,j,k) * cell_data(i-1, j, k  , Rho_comp)
                            + vfrac(i,j,k) * cell_data(i, j, k  , Rho_comp)
                            + vfrac(i-1,j,k-1) * cell_data(i-1, j, k-1, Rho_comp)
                            + vfrac(i,j,k-1) * cell_data(i, j, k-1, Rho_comp) )/ vol_sum;
            } else {
                rho_bar = fourth*( cell_data(i-1, j, k  , Rho_comp) + cell_data(i, j, k  , Rho_comp)
                                + cell_data(i-1, j, k-1, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            }
            tau13(i,j,k) *= -rho_bar * mu_eff;

            if (tau13i) tau13i(i,j,k) *= -rho_bar * mu_eff;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i,j-1,k) + vfrac(i,j,k-1) + vfrac(i,j-1,k-1);
            Real rho_bar = zero;
            if (vol_sum > Real(1.e-16)) {
                rho_bar = ( vfrac(i,j-1,k) * cell_data(i, j-1, k  , Rho_comp)
                            + vfrac(i,j,k) * cell_data(i, j, k  , Rho_comp)
                            + vfrac(i,j-1,k-1) * cell_data(i, j-1, k-1, Rho_comp)
                            + vfrac(i,j,k-1) * cell_data(i, j, k-1, Rho_comp) ) / vol_sum;
            } else {
                rho_bar = fourth*( cell_data(i, j-1, k  , Rho_comp) + cell_data(i, j, k  , Rho_comp)
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

/**
 * Function for computing the stress with variable viscosity for EB.
 * rho_bar at off-diagonal locations uses vfrac-weighted averaging as in
 * ComputeStressConsVisc_EB; mu_turb is averaged with the same vfrac weights.
 *
 * @param[in] bxcc cell center box for tau_ii
 * @param[in] tbxxy nodal xy box for tau_12
 * @param[in] tbxxz nodal xz box for tau_13
 * @param[in] tbxyz nodal yz box for tau_23
 * @param[in] mu_eff constant molecular viscosity
 * @param[in] mu_turb turbulent viscosity (cell-centered)
 * @param[in] cell_data to access rho if ConstantAlpha
 * @param[in,out] tau11 11 strain -> stress
 * @param[in,out] tau22 22 strain -> stress
 * @param[in,out] tau33 33 strain -> stress
 * @param[in,out] tau12 12 strain -> stress
 * @param[in,out] tau13 13 strain -> stress
 * @param[in,out] tau23 23 strain -> stress
 * @param[in] er_arr expansion rate
 * @param[in] vfrac volume fractions
 * @param[in,out] tau13i contribution to stress from du/dz
 * @param[in,out] tau23i contribution to stress from dv/dz
 * @param[in,out] tau33i contribution to stress from dw/dz
 */
void
ComputeStressVarVisc_EB (Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Real mu_eff,
                         const Array4<const Real>& mu_turb,
                         const Array4<const Real>& cell_data,
                         Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                         Array4<Real>& tau12, Array4<Real>& tau13, Array4<Real>& tau23,
                         const Array4<const Real>& er_arr,
                         Array4<const Real>& vfrac,
                         Array4<Real>& tau13i,
                         Array4<Real>& tau23i,
                         Array4<Real>& tau33i)
{
    Real OneThird = (one/three);

    // NOTE: mu_eff includes factor of 2

    if (cell_data)
    // constant alpha (stored in mu_eff)
    {
        // Cell centered strains (no EB correction needed — cell-centered quantities)
        ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (vfrac(i,j,k) > zero) {
                Real rhoAlpha = cell_data(i, j, k, Rho_comp) * mu_eff;
                Real mu_11 = rhoAlpha + two * mu_turb(i, j, k, EddyDiff::Mom_h);
                Real mu_22 = mu_11;
                Real mu_33 = rhoAlpha + two * mu_turb(i, j, k, EddyDiff::Mom_v);
                if (tau33i) tau33i(i,j,k) = -mu_33 * tau33(i,j,k);
                tau11(i,j,k) = -mu_11 * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
                tau22(i,j,k) = -mu_22 * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
                tau33(i,j,k) = -mu_33 * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
            } else {
                if (tau33i) tau33i(i,j,k) = zero;
                tau11(i,j,k) = zero;
                tau22(i,j,k) = zero;
                tau33(i,j,k) = zero;
            }
        });

        // Off-diagonal strains: vfrac-weighted rho_bar and mu_bar
        ParallelFor(tbxxy,tbxxz,tbxyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i-1,j,k) + vfrac(i,j-1,k) + vfrac(i-1,j-1,k);
            Real rho_bar, mu_bar;
            if (vol_sum > Real(1.e-16)) {
                rho_bar = ( vfrac(i-1,j  ,k) * cell_data(i-1, j  , k, Rho_comp)
                          + vfrac(i  ,j  ,k) * cell_data(i  , j  , k, Rho_comp)
                          + vfrac(i-1,j-1,k) * cell_data(i-1, j-1, k, Rho_comp)
                          + vfrac(i  ,j-1,k) * cell_data(i  , j-1, k, Rho_comp) ) / vol_sum;
                mu_bar  = ( vfrac(i-1,j  ,k) * mu_turb(i-1, j  , k, EddyDiff::Mom_h)
                          + vfrac(i  ,j  ,k) * mu_turb(i  , j  , k, EddyDiff::Mom_h)
                          + vfrac(i-1,j-1,k) * mu_turb(i-1, j-1, k, EddyDiff::Mom_h)
                          + vfrac(i  ,j-1,k) * mu_turb(i  , j-1, k, EddyDiff::Mom_h) ) / vol_sum;
            } else {
                rho_bar = fourth*( cell_data(i-1, j  , k, Rho_comp) + cell_data(i  , j  , k, Rho_comp)
                               + cell_data(i-1, j-1, k, Rho_comp) + cell_data(i  , j-1, k, Rho_comp) );
                mu_bar  = fourth*( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i  , j  , k, EddyDiff::Mom_h)
                               + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i  , j-1, k, EddyDiff::Mom_h) );
            }
            Real mu_12 = rho_bar*mu_eff + two*mu_bar;
            tau12(i,j,k) *= -mu_12;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i-1,j,k) + vfrac(i,j,k-1) + vfrac(i-1,j,k-1);
            Real rho_bar, mu_bar;
            if (vol_sum > Real(1.e-16)) {
                rho_bar = ( vfrac(i-1,j,k  ) * cell_data(i-1, j, k  , Rho_comp)
                          + vfrac(i  ,j,k  ) * cell_data(i  , j, k  , Rho_comp)
                          + vfrac(i-1,j,k-1) * cell_data(i-1, j, k-1, Rho_comp)
                          + vfrac(i  ,j,k-1) * cell_data(i  , j, k-1, Rho_comp) ) / vol_sum;
                mu_bar  = ( vfrac(i-1,j,k  ) * mu_turb(i-1, j, k  , EddyDiff::Mom_v)
                          + vfrac(i  ,j,k  ) * mu_turb(i  , j, k  , EddyDiff::Mom_v)
                          + vfrac(i-1,j,k-1) * mu_turb(i-1, j, k-1, EddyDiff::Mom_v)
                          + vfrac(i  ,j,k-1) * mu_turb(i  , j, k-1, EddyDiff::Mom_v) ) / vol_sum;
            } else {
                rho_bar = fourth*( cell_data(i-1, j, k  , Rho_comp) + cell_data(i  , j, k  , Rho_comp)
                               + cell_data(i-1, j, k-1, Rho_comp) + cell_data(i  , j, k-1, Rho_comp) );
                mu_bar  = fourth*( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i  , j, k  , EddyDiff::Mom_v)
                               + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i  , j, k-1, EddyDiff::Mom_v) );
            }
            Real mu_13 = rho_bar*mu_eff + two*mu_bar;
            tau13(i,j,k) *= -mu_13;
            if (tau13i) tau13i(i,j,k) *= -mu_13;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i,j-1,k) + vfrac(i,j,k-1) + vfrac(i,j-1,k-1);
            Real rho_bar, mu_bar;
            if (vol_sum > Real(1.e-16)) {
                rho_bar = ( vfrac(i,j-1,k  ) * cell_data(i, j-1, k  , Rho_comp)
                          + vfrac(i,j  ,k  ) * cell_data(i, j  , k  , Rho_comp)
                          + vfrac(i,j-1,k-1) * cell_data(i, j-1, k-1, Rho_comp)
                          + vfrac(i,j  ,k-1) * cell_data(i, j  , k-1, Rho_comp) ) / vol_sum;
                mu_bar  = ( vfrac(i,j-1,k  ) * mu_turb(i, j-1, k  , EddyDiff::Mom_v)
                          + vfrac(i,j  ,k  ) * mu_turb(i, j  , k  , EddyDiff::Mom_v)
                          + vfrac(i,j-1,k-1) * mu_turb(i, j-1, k-1, EddyDiff::Mom_v)
                          + vfrac(i,j  ,k-1) * mu_turb(i, j  , k-1, EddyDiff::Mom_v) ) / vol_sum;
            } else {
                rho_bar = fourth*( cell_data(i, j-1, k  , Rho_comp) + cell_data(i, j  , k  , Rho_comp)
                               + cell_data(i, j-1, k-1, Rho_comp) + cell_data(i, j  , k-1, Rho_comp) );
                mu_bar  = fourth*( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j  , k  , EddyDiff::Mom_v)
                               + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j  , k-1, EddyDiff::Mom_v) );
            }
            Real mu_23 = rho_bar*mu_eff + two*mu_bar;
            tau23(i,j,k) *= -mu_23;
            if (tau23i) tau23i(i,j,k) *= -mu_23;
        });
    }
    else
    // constant mu_eff
    {
        // Cell centered strains
        ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (vfrac(i,j,k) > zero) {
                Real mu_11 = mu_eff + two * mu_turb(i, j, k, EddyDiff::Mom_h);
                Real mu_22 = mu_11;
                Real mu_33 = mu_eff + two * mu_turb(i, j, k, EddyDiff::Mom_v);
                if (tau33i) tau33i(i,j,k) = -mu_33 * tau33(i,j,k);
                tau11(i,j,k) = -mu_11 * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
                tau22(i,j,k) = -mu_22 * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
                tau33(i,j,k) = -mu_33 * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
            } else {
                if (tau33i) tau33i(i,j,k) = zero;
                tau11(i,j,k) = zero;
                tau22(i,j,k) = zero;
                tau33(i,j,k) = zero;
            }
        });

        // Off-diagonal strains: vfrac-weighted mu_bar
        ParallelFor(tbxxy,tbxxz,tbxyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i-1,j,k) + vfrac(i,j-1,k) + vfrac(i-1,j-1,k);
            Real mu_bar;
            if (vol_sum > Real(1.e-16)) {
                mu_bar = ( vfrac(i-1,j  ,k) * mu_turb(i-1, j  , k, EddyDiff::Mom_h)
                         + vfrac(i  ,j  ,k) * mu_turb(i  , j  , k, EddyDiff::Mom_h)
                         + vfrac(i-1,j-1,k) * mu_turb(i-1, j-1, k, EddyDiff::Mom_h)
                         + vfrac(i  ,j-1,k) * mu_turb(i  , j-1, k, EddyDiff::Mom_h) ) / vol_sum;
            } else {
                mu_bar = fourth*( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i  , j  , k, EddyDiff::Mom_h)
                              + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i  , j-1, k, EddyDiff::Mom_h) );
            }
            Real mu_12 = mu_eff + two*mu_bar;
            tau12(i,j,k) *= -mu_12;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i-1,j,k) + vfrac(i,j,k-1) + vfrac(i-1,j,k-1);
            Real mu_bar;
            if (vol_sum > Real(1.e-16)) {
                mu_bar = ( vfrac(i-1,j,k  ) * mu_turb(i-1, j, k  , EddyDiff::Mom_v)
                         + vfrac(i  ,j,k  ) * mu_turb(i  , j, k  , EddyDiff::Mom_v)
                         + vfrac(i-1,j,k-1) * mu_turb(i-1, j, k-1, EddyDiff::Mom_v)
                         + vfrac(i  ,j,k-1) * mu_turb(i  , j, k-1, EddyDiff::Mom_v) ) / vol_sum;
            } else {
                mu_bar = fourth*( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i  , j, k  , EddyDiff::Mom_v)
                              + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i  , j, k-1, EddyDiff::Mom_v) );
            }
            Real mu_13 = mu_eff + two*mu_bar;
            tau13(i,j,k) *= -mu_13;
            if (tau13i) tau13i(i,j,k) *= -mu_13;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real vol_sum = vfrac(i,j,k) + vfrac(i,j-1,k) + vfrac(i,j,k-1) + vfrac(i,j-1,k-1);
            Real mu_bar;
            if (vol_sum > Real(1.e-16)) {
                mu_bar = ( vfrac(i,j-1,k  ) * mu_turb(i, j-1, k  , EddyDiff::Mom_v)
                         + vfrac(i,j  ,k  ) * mu_turb(i, j  , k  , EddyDiff::Mom_v)
                         + vfrac(i,j-1,k-1) * mu_turb(i, j-1, k-1, EddyDiff::Mom_v)
                         + vfrac(i,j  ,k-1) * mu_turb(i, j  , k-1, EddyDiff::Mom_v) ) / vol_sum;
            } else {
                mu_bar = fourth*( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j  , k  , EddyDiff::Mom_v)
                              + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j  , k-1, EddyDiff::Mom_v) );
            }
            Real mu_23 = mu_eff + two*mu_bar;
            tau23(i,j,k) *= -mu_23;
            if (tau23i) tau23i(i,j,k) *= -mu_23;
        });
    }
}
