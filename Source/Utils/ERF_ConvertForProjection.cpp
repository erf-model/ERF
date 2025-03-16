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
ConvertForProjection (const MultiFab& den_div, const MultiFab& den_mlt,
                      MultiFab& xmom, MultiFab& ymom, MultiFab& zmom,
                      const Box& domain, const Vector<BCRec>& domain_bcs_type_h)
{
    BL_PROFILE_VAR("ConvertForProjection()",onvertForProjection);

    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(den_div,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // We need velocity in the interior ghost cells (init == real)
        Box bx = mfi.tilebox();

        Box tbx = surroundingNodes(bx,0);
        Box tby = surroundingNodes(bx,1);
        Box tbz = surroundingNodes(bx,2);

        if ( (bx.smallEnd(0) == domain.smallEnd(0)) &&
             ( (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::ext_dir) ||
               (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::ext_dir_upwind) ) )
        {
            tbx.growLo(0,-1);
        }
        if ( (bx.bigEnd(0) == domain.bigEnd(0)) &&
             ( (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir) ||
               (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir_upwind) ) )
        {
            tbx.growHi(0,-1);
        }
        if ( (bx.smallEnd(1) == domain.smallEnd(1)) &&
             ( (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::ext_dir) ||
               (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::ext_dir_upwind) ) )
        {
            tby.growLo(1,-1);
        }
        if ( (bx.bigEnd(1) == domain.bigEnd(1)) &&
             ( (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::ext_dir) ||
               (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::ext_dir_upwind) ) )
        {
            tby.growHi(1,-1);
        }

        // Conserved variables on cell centers -- we use this for density
        const Array4<const Real>& den_div_arr = den_div.const_array(mfi);
        const Array4<const Real>& den_mlt_arr = den_mlt.const_array(mfi);

        // Momentum on faces
        Array4<Real> const& momx = xmom.array(mfi);
        Array4<Real> const& momy = ymom.array(mfi);
        Array4<Real> const& momz = zmom.array(mfi);

        ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momx(i,j,k) *= ( den_mlt_arr(i,j,k,Rho_comp) + den_mlt_arr(i-1,j,k,Rho_comp) )
                         / ( den_div_arr(i,j,k,Rho_comp) + den_div_arr(i-1,j,k,Rho_comp) );
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) *= ( den_mlt_arr(i,j,k,Rho_comp) + den_mlt_arr(i,j-1,k,Rho_comp) )
                         / ( den_div_arr(i,j,k,Rho_comp) + den_div_arr(i,j-1,k,Rho_comp) );
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) *= ( den_mlt_arr(i,j,k,Rho_comp) + den_mlt_arr(i,j,k-1,Rho_comp) )
                         / ( den_div_arr(i,j,k,Rho_comp) + den_div_arr(i,j,k-1,Rho_comp) );
        });

        if (bx.smallEnd(0) == domain.smallEnd(0)) {
            if (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::ext_dir)
            {
                ParallelFor(makeSlab(tbx,0,domain.smallEnd(0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    momx(i,j,k) *= den_mlt_arr(i-1,j,k,Rho_comp) / den_div_arr(i-1,j,k,Rho_comp) ;
                });
            }
            else if (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::ext_dir_upwind)
            {
                ParallelFor(makeSlab(tbx,0,domain.smallEnd(0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (momx(i,j,k) >= 0.) {
                        momx(i,j,k) *= den_mlt_arr(i-1,j,k,Rho_comp) / den_div_arr(i-1,j,k,Rho_comp) ;
                    } else {
                        momx(i,j,k) *= ( den_mlt_arr(i,j,k,Rho_comp) + den_mlt_arr(i-1,j,k,Rho_comp) )
                                     / ( den_div_arr(i,j,k,Rho_comp) + den_div_arr(i-1,j,k,Rho_comp) );
                    }
                });
            }
        }

        if (bx.bigEnd(0) == domain.bigEnd(0)) {
            if (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir)
            {
                ParallelFor(makeSlab(tbx,0,domain.bigEnd(0)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    momx(i,j,k) *= den_mlt_arr(i-1,j,k,Rho_comp) / den_div_arr(i-1,j,k,Rho_comp) ;
                });
            }
            else if (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir_upwind)
            {
                ParallelFor(makeSlab(tbx,0,domain.smallEnd(0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (momx(i,j,k) <= 0.) {
                        momx(i,j,k) *= den_mlt_arr(i-1,j,k,Rho_comp) / den_div_arr(i-1,j,k,Rho_comp) ;
                    } else {
                        momx(i,j,k) *= ( den_mlt_arr(i,j,k,Rho_comp) + den_mlt_arr(i-1,j,k,Rho_comp) )
                                     / ( den_div_arr(i,j,k,Rho_comp) + den_div_arr(i-1,j,k,Rho_comp) );
                    }
                });
            }
        }

        if (bx.smallEnd(1) == domain.smallEnd(1)) {
            if (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::ext_dir)
            {
                ParallelFor(makeSlab(tby,1,domain.smallEnd(1)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    momy(i,j,k) *= den_mlt_arr(i-1,j,k,Rho_comp) / den_div_arr(i-1,j,k,Rho_comp) ;
                });
            }
            else if (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::ext_dir_upwind)
            {
                ParallelFor(makeSlab(tby,1,domain.smallEnd(1)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (momy(i,j,k) >= 0.) {
                        momy(i,j,k) *= den_mlt_arr(i-1,j,k,Rho_comp) / den_div_arr(i-1,j,k,Rho_comp) ;
                    } else {
                        momy(i,j,k) *= ( den_mlt_arr(i,j,k,Rho_comp) + den_mlt_arr(i-1,j,k,Rho_comp) )
                                     / ( den_div_arr(i,j,k,Rho_comp) + den_div_arr(i-1,j,k,Rho_comp) );
                    }
                });
            }
        }

        if (bx.bigEnd(1) == domain.bigEnd(1)) {
            if (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::ext_dir)
            {
                ParallelFor(makeSlab(tby,1,domain.bigEnd(1)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    momy(i,j,k) *= den_mlt_arr(i-1,j,k,Rho_comp) / den_div_arr(i-1,j,k,Rho_comp) ;
                });
            }
            else if (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::ext_dir_upwind)
            {
                ParallelFor(makeSlab(tby,1,domain.bigEnd(1)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (momy(i,j,k) <= 0.) {
                        momy(i,j,k) *= den_mlt_arr(i-1,j,k,Rho_comp) / den_div_arr(i-1,j,k,Rho_comp) ;
                    } else {
                        momy(i,j,k) *= ( den_mlt_arr(i,j,k,Rho_comp) + den_mlt_arr(i-1,j,k,Rho_comp) )
                                     / ( den_div_arr(i,j,k,Rho_comp) + den_div_arr(i-1,j,k,Rho_comp) );
                    }
                });
            }
        }


    } // end MFIter
}
