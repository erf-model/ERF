/**
 * \file ERF_MomentumToVelocity.cpp
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <ERF_Utils.H>

using namespace amrex;

/**
 * Convert momentum to velocity by dividing by density averaged onto faces
 *
 * @param[out] xvel x-component of velocity
 * @param[out] yvel y-component of velocity
 * @param[out] zvel z-component of velocity
 * @param[in] density density at cell centers
 * @param[in] xmom_in x-component of momentum
 * @param[in] ymom_in y-component of momentum
 * @param[in] zmom_in z-component of momentum
 * @param[in] domain  Domain at this level
 * @param[in] domain_bcs_type_h   host vector for domain boundary conditions
 */

void
MomentumToVelocity (MultiFab& xvel, MultiFab& yvel, MultiFab& zvel,
                    const MultiFab& density,
                    const MultiFab& xmom_in, const MultiFab& ymom_in, const MultiFab& zmom_in,
                    const Box& domain,
                    const Vector<BCRec>& domain_bcs_type_h,
                    const MultiFab* c_vfrac // optional
) {
    BL_PROFILE_VAR("MomentumToVelocity()",MomentumToVelocity);

    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // We need velocity in the interior ghost cells (init == real)
        Box bx = mfi.tilebox();

        const Box& tbx = surroundingNodes(bx,0);
        const Box& tby = surroundingNodes(bx,1);
        const Box& tbz = surroundingNodes(bx,2);

        // Conserved variables on cell centers -- we use this for density
        const Array4<const Real>& dens_arr = density.array(mfi);

        // Momentum on faces
        Array4<Real const> const& momx = xmom_in.const_array(mfi);
        Array4<Real const> const& momy = ymom_in.const_array(mfi);
        Array4<Real const> const& momz = zmom_in.const_array(mfi);

        // Velocity on faces
        const Array4<Real>& velx = xvel.array(mfi);
        const Array4<Real>& vely = yvel.array(mfi);
        const Array4<Real>& velz = zvel.array(mfi);

        if (c_vfrac==nullptr) {
            ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                velx(i,j,k) = momx(i,j,k) * two / (dens_arr(i,j,k,Rho_comp) + dens_arr(i-1,j,k,Rho_comp));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                vely(i,j,k) = momy(i,j,k) * two / (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j-1,k,Rho_comp));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                velz(i,j,k) = momz(i,j,k) * two / (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j,k-1,Rho_comp));
            });
        } else {
            // EB
            const Array4<const Real>& c_vfrac_arr = c_vfrac->const_array(mfi);

            ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real vfrac_i = c_vfrac_arr(i,j,k);
                Real vfrac_im1 = c_vfrac_arr(i-1,j,k);
                Real vfrac_sum = vfrac_i + vfrac_im1;
                if (vfrac_sum > zero) {
                    Real rho = (vfrac_i > zero ? vfrac_i * dens_arr(i,j,k,Rho_comp) : zero)
                             + (vfrac_im1 > zero ? vfrac_im1 * dens_arr(i-1,j,k,Rho_comp) : zero);
                    rho /= vfrac_sum;
                    velx(i,j,k) = momx(i,j,k) / rho;
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real vfrac_j = c_vfrac_arr(i,j,k);
                Real vfrac_jm1 = c_vfrac_arr(i,j-1,k);
                Real vfrac_sum = vfrac_j + vfrac_jm1;
                if (vfrac_sum > zero) {
                    Real rho = (vfrac_j > zero ? vfrac_j * dens_arr(i,j,k,Rho_comp) : zero)
                             + (vfrac_jm1 > zero ? vfrac_jm1 * dens_arr(i,j-1,k,Rho_comp) : zero);
                    rho /= vfrac_sum;
                    vely(i,j,k) = momy(i,j,k) / rho;
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real vfrac_k = c_vfrac_arr(i,j,k);
                Real vfrac_km1 = c_vfrac_arr(i,j,k-1);
                Real vfrac_sum = vfrac_k + vfrac_km1;
                if (vfrac_sum > zero) {
                    Real rho = (vfrac_k > zero ? vfrac_k * dens_arr(i,j,k,Rho_comp) : zero)
                             + (vfrac_km1 > zero ? vfrac_km1 * dens_arr(i,j,k-1,Rho_comp) : zero);
                    rho /= vfrac_sum;
                    velz(i,j,k) = momz(i,j,k) / rho;
                }
            });
        }

        if (bx.smallEnd(0) == domain.smallEnd(0)) {
            if (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::ext_dir)
            {
                ParallelFor(makeSlab(tbx,0,domain.smallEnd(0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    velx(i,j,k) = momx(i,j,k) / dens_arr(i-1,j,k,Rho_comp);
                });
            }
            else if (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::ext_dir_upwind)
            {
                ParallelFor(makeSlab(tbx,0,domain.smallEnd(0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (momx(i,j,k) >= zero) {
                        velx(i,j,k) = momx(i,j,k) / dens_arr(i-1,j,k,Rho_comp);
                    }
                });
            }
        }

        if (bx.bigEnd(0) == domain.bigEnd(0)) {
            if (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir)
            {
                ParallelFor(makeSlab(tbx,0,domain.bigEnd(0)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    velx(i,j,k) = momx(i,j,k) / dens_arr(i,j,k,Rho_comp);
                });
            }
            else if (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir_upwind)
            {
                ParallelFor(makeSlab(tbx,0,domain.bigEnd(0)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (momx(i,j,k) <= zero) {
                        velx(i,j,k) = momx(i,j,k) / dens_arr(i,j,k,Rho_comp);
                    }
                });
            }
        }

        if (bx.smallEnd(1) == domain.smallEnd(1)) {
            if (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::ext_dir)
            {
                ParallelFor(makeSlab(tby,1,domain.smallEnd(1)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    vely(i,j,k) = momy(i,j,k) / dens_arr(i,j-1,k,Rho_comp);
                });
            }
            else if (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::ext_dir_upwind)
            {
                ParallelFor(makeSlab(tby,1,domain.smallEnd(1)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (momy(i,j,k) >= zero) {
                        vely(i,j,k) = momy(i,j,k) / dens_arr(i,j-1,k,Rho_comp);
                    }
                });
            }
        }

        if (bx.bigEnd(1) == domain.bigEnd(1)) {
            if (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::ext_dir)
            {
                ParallelFor(makeSlab(tby,1,domain.bigEnd(1)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    vely(i,j,k) = momy(i,j,k) / dens_arr(i,j,k,Rho_comp);
                });
            }
            else if (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::ext_dir_upwind)
            {
                ParallelFor(makeSlab(tby,1,domain.bigEnd(1)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (momy(i,j,k) <= zero) {
                        vely(i,j,k) = momy(i,j,k) / dens_arr(i,j,k,Rho_comp);
                    }
                });
            }
        }
    } // end MFIter
}
