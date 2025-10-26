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

void compute_influx_outflux(
    Array<MultiFab*, AMREX_SPACEDIM>& vels_vec,
    Array<MultiFab*, AMREX_SPACEDIM>& area_vec,
    const Geometry& geom,
    Real& influx, Real& outflux)
{
    influx = 0.0, outflux = 0.0;

    const Box domain = geom.Domain();
    const auto& domlo = lbound(domain);
    const auto& domhi = ubound(domain);

    // Normal face area (of undistorted mesh)
    const Real* a_dx = geom.CellSize();
    const Real ds_x = a_dx[1];
    const Real ds_y = a_dx[0];

    IntVect ngrow = {0,0,0};

    // X-dir
    auto const&  vel_x = vels_vec[0]->const_arrays();
    auto const& area_x = area_vec[0]->const_arrays();
    influx += ds_x *
        ParReduce(TypeList<ReduceOpSum>{},
                  TypeList<Real>{},
                  *vels_vec[0], ngrow,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
            noexcept -> GpuTuple<Real>
        {
            if ( (i == domlo.x   && vel_x[box_no](i,j,k) > 0.0) ||
                 (i == domhi.x+1 && vel_x[box_no](i,j,k) < 0.0) ) {
                return { std::abs(vel_x[box_no](i,j,k)) * area_x[box_no](i,j,k) };
            } else {
                return { 0. };
            }
        });

    outflux += ds_x *
        ParReduce(TypeList<ReduceOpSum>{},
                  TypeList<Real>{},
                  *vels_vec[0], ngrow,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
            noexcept -> GpuTuple<Real>
        {
            if ( (i == domlo.x   && vel_x[box_no](i,j,k) < 0.0) ||
                 (i == domhi.x+1 && vel_x[box_no](i,j,k) > 0.0) ) {
                return { std::abs(vel_x[box_no](i,j,k)) * area_x[box_no](i,j,k) };
            } else {
                return { 0. };
            }
        });

    // Y-dir
    auto const&  vel_y=  vels_vec[1]->const_arrays();
    auto const& area_y = area_vec[1]->const_arrays();
    influx += ds_y *
        ParReduce(TypeList<ReduceOpSum>{},
                  TypeList<Real>{},
                  *vels_vec[1], ngrow,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
            noexcept -> GpuTuple<Real>
        {
            if ( (j == domlo.y   && vel_y[box_no](i,j,k) > 0.0) ||
                 (j == domhi.y+1 && vel_y[box_no](i,j,k) < 0.0) ) {
                return { std::abs(vel_y[box_no](i,j,k)) * area_y[box_no](i,j,k) };
            } else {
                return { 0. };
            }
        });

    outflux += ds_y *
        ParReduce(TypeList<ReduceOpSum>{},
                  TypeList<Real>{},
                  *vels_vec[1], ngrow,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
            noexcept -> GpuTuple<Real>
        {
            if ( (j == domlo.y   && vel_y[box_no](i,j,k) < 0.0) ||
                 (j == domhi.y+1 && vel_y[box_no](i,j,k) > 0.0) ) {
                return { std::abs(vel_y[box_no](i,j,k)) * area_y[box_no](i,j,k) };
            } else {
                return { 0. };
            }
        });

    ParallelDescriptor::ReduceRealSum(influx);
    ParallelDescriptor::ReduceRealSum(outflux);
}

void correct_outflow(
    const Geometry& geom_lev,
    Array<MultiFab*, AMREX_SPACEDIM>& vels_vec,
    const Box& domain,
    const Real alpha_fcf)
{
    for (OrientationIter oit; oit != nullptr; ++oit) {
        const auto ori = oit();
        const int dir = ori.coordDir();
        const auto oriIsLow = ori.isLow();
        const auto oriIsHigh = ori.isHigh();

        if (dir < 2) {

        // MultiFab for normal velocity
        const auto& vel_mf = vels_vec[dir];

        // Domain extent indices for the velocities
        int dlo = domain.smallEnd(dir);
        int dhi = domain.bigEnd(dir)+1;

        // get BCs for the normal velocity and set the boundary index
        int bndry;
        if (oriIsLow) {
            bndry = dlo;
        } else {
            bndry = dhi;
        }

        // Assume here that all domain boundaries are viewed as direction_dependent!
        if (!geom_lev.isPeriodic(dir))
        {
            for (MFIter mfi(*vel_mf, false); mfi.isValid(); ++mfi) {

                Box box = mfi.validbox();

                // Enter further only if the box boundary is at the domain boundary
                if ((oriIsLow  && (box.smallEnd(dir) == dlo))
                 || (oriIsHigh && (box.bigEnd(dir)   == dhi))) {

                    // create a 2D box normal to dir at the low/high boundary
                    Box box2d(box); box2d.setRange(dir, bndry);

                    auto vel_arr = vel_mf->array(mfi);

                    ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if ((oriIsLow  && vel_arr(i,j,k) < 0)
                         || (oriIsHigh && vel_arr(i,j,k) > 0)) {
                            vel_arr(i,j,k) *= alpha_fcf;
                        }
                    });
                }
            } // mfi
        } // geom
      } // dir < 2
    } // ori
}

void enforceInOutSolvability (int /*lev*/,
    Array<MultiFab*, AMREX_SPACEDIM>& vels_vec,
    Array<MultiFab*, AMREX_SPACEDIM>& area_vec,
    const Geometry& geom)
{
    Real small_vel = 1.e-8;

    Real influx = 0.0, outflux = 0.0;

    const Box domain = geom.Domain();

    Real influx_lev = 0.0, outflux_lev = 0.0;
    compute_influx_outflux(vels_vec, area_vec, geom, influx_lev, outflux_lev);
    influx += influx_lev;
    outflux += outflux_lev;
    amrex::Print() <<" TOTAL INFLUX / OUTFLOW " << influx << " " << outflux << std::endl;

    if ((influx > small_vel) && (outflux < small_vel)) {
        Abort("Cannot enforce solvability, no outflow from the direction dependent boundaries");
    } else if ((influx < small_vel) && (outflux < small_vel)) {
        return; // do nothing
    } else {
        const Real alpha_fcf = influx/outflux;  // flux correction factor
        correct_outflow(geom, vels_vec, domain, alpha_fcf);

        // Just for diagnostic purposes!
        // Real influx_lev = 0.0, outflux_lev = 0.0;
        // compute_influx_outflux(vels_vec, area_vec, geom, influx_lev, outflux_lev);
        // amrex::Print() <<" TOTAL INFLUX / OUTFLOW " << influx_lev << " " << outflux_lev << std::endl;
    }
}
