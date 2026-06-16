/**
 * \file ERF_EnforceConstraintOnBdy.cpp
 */
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>

#include <ERF_Utils.H>

using namespace amrex;

namespace
{
void
scale_bdy_normal_by_rho0 (const MultiFab& rho0,
                          FArrayBox& bdy_fab,
                          const int dir,
                          const Box& domain,
                          const Vector<BCRec>& domain_bcs_type_h,
                          const bool multiply_by_rho0)
{
    AMREX_ALWAYS_ASSERT(bdy_fab.nComp() == 1);

    FArrayBox bdy_tmp(bdy_fab.box(), 1, The_Managed_Arena());
    bdy_tmp.template setVal<RunOn::Device>(zero);

    const BCRec* bc_ptr_h = domain_bcs_type_h.data();
    const int domlo = domain.smallEnd(dir);
    const int domhi = domain.bigEnd(dir) + 1;
    const int bc_lo = bc_ptr_h[BCVars::cons_bc].lo(dir);
    const int bc_hi = bc_ptr_h[BCVars::cons_bc].hi(dir);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rho0,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box bx = mfi.nodaltilebox(dir) & bdy_fab.box();
        if (!bx.ok()) { continue; }

        const Array4<const Real>& rho0_arr = rho0.const_array(mfi);
        const Array4<const Real>& bdy_arr  = bdy_fab.const_array();
        const Array4<Real>&       tmp_arr  = bdy_tmp.array();

        // amrex::Print() << "MULTIPLYING BY RHO0 ON BX " << dir << " " << bx << std::endl;

        if (dir == 0) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real normal_vel = bdy_arr(i,j,k);
                const bool use_xlo =
                    (i == domlo) &&
                    ((bc_lo == ERFBCType::ext_dir) ||
                     (bc_lo == ERFBCType::ext_dir_upwind && normal_vel >= zero));
                const bool use_xhi =
                    (i == domhi) &&
                    ((bc_hi == ERFBCType::ext_dir) ||
                     (bc_hi == ERFBCType::ext_dir_upwind && normal_vel <= zero));

                Real fac;
                if (use_xlo) {
                    fac = rho0_arr(i-1,j,k,Rho_comp);
                } else if (use_xhi) {
                    fac = rho0_arr(i,j,k,Rho_comp);
                } else {
                    fac = Real(0.5) * (rho0_arr(i,j,k,Rho_comp) + rho0_arr(i-1,j,k,Rho_comp));
                }
                tmp_arr(i,j,k) = multiply_by_rho0 ? normal_vel * fac : normal_vel / fac;
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real normal_vel = bdy_arr(i,j,k);
                const bool use_ylo =
                    (j == domlo) &&
                    ((bc_lo == ERFBCType::ext_dir) ||
                     (bc_lo == ERFBCType::ext_dir_upwind && normal_vel >= zero));
                const bool use_yhi =
                    (j == domhi) &&
                    ((bc_hi == ERFBCType::ext_dir) ||
                     (bc_hi == ERFBCType::ext_dir_upwind && normal_vel <= zero));

                Real fac;
                if (use_ylo) {
                    fac = rho0_arr(i,j-1,k,Rho_comp);
                } else if (use_yhi) {
                    fac = rho0_arr(i,j,k,Rho_comp);
                } else {
                    fac = Real(0.5) * (rho0_arr(i,j,k,Rho_comp) + rho0_arr(i,j-1,k,Rho_comp));
                }
                tmp_arr(i,j,k) = multiply_by_rho0 ? normal_vel * fac : normal_vel / fac;
            });
        }
    }

    ParallelAllReduce::Sum(bdy_tmp.dataPtr(), bdy_tmp.size(), ParallelContext::CommunicatorAll());
    bdy_fab.template copy<RunOn::Device>(bdy_tmp, 0, 0, 1);
}

void
correct_bdy_outflow_on_face (FArrayBox& bdy_fab,
                             const int dir,
                             const bool is_low,
                             const Box& domain,
                             const Real alpha_fcf,
                             const int n)
{
    const int bndry = is_low ? domain.smallEnd(dir) + n : domain.bigEnd(dir)+1 - n;
    Box bx = bdy_fab.box();

    if (bx.smallEnd(dir) > bndry || bx.bigEnd(dir) < bndry) {
        return;
    }
    // amrex::Print() << "CORRECTING OUTFLOW VALS BY MULTIPLIER " << alpha_fcf << std::endl;

    bx.setRange(dir, bndry);
    bx.grow(1-dir, -n);

    const Array4<Real>& bdy_arr = bdy_fab.array();

    // amrex::Print() << "OPERATING ON BX " << bx << std::endl;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ((is_low && bdy_arr(i,j,k) < zero) ||
           (!is_low && bdy_arr(i,j,k) > zero)) {
            bdy_arr(i,j,k) *= alpha_fcf;
        }
    });
}
}

void compute_influx_outflux_bdy (FArrayBox& bdy_data_xlo, FArrayBox& bdy_data_xhi,
                                 FArrayBox& bdy_data_ylo, FArrayBox& bdy_data_yhi,
                                 Array<MultiFab*, AMREX_SPACEDIM>& area_vec,
                                 const Geometry& geom,
                                 Real& influx, Real& outflux,
                                 const int n)
{
    BL_PROFILE_VAR("compute_influx_outflux_bdy()", computeInfluxOutfluxBdy);

    influx = zero, outflux = zero;

    const Box domain = geom.Domain();
    const auto& domlo = lbound(domain);
    const auto& domhi = ubound(domain);

    // Normal face area (of undistorted mesh)
    const Real* a_dx = geom.CellSize();
    const Real ds_x = a_dx[1];
    const Real ds_y = a_dx[0];

    IntVect ngrow = {0,0,0};

    // X-dir
    const auto& bdatxlo = bdy_data_xlo.const_array();
    const auto& bdatxhi = bdy_data_xhi.const_array();
    auto const& area_x = area_vec[0]->const_arrays();
    influx += ds_x *
        ParReduce(TypeList<ReduceOpSum>{},
                  TypeList<Real>{},
                  *area_vec[0], ngrow,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
            noexcept -> GpuTuple<Real>
        {
            if ( (i == domlo.x + n) && (j >= domlo.y+n && j<= domhi.y-n) && bdatxlo(i,j,k) > zero) {
                return { std::abs(bdatxlo(i,j,k)) * area_x[box_no](i,j,k) };
            } else if ( (i == domhi.x+1 - n) && (j >= domlo.y+n && j<= domhi.y-n) && bdatxhi(i,j,k) < zero) {
                return { std::abs(bdatxhi(i,j,k)) * area_x[box_no](i,j,k) };
            } else {
                return { zero };
            }
        });

    outflux += ds_x *
        ParReduce(TypeList<ReduceOpSum>{},
                  TypeList<Real>{},
                  *area_vec[0], ngrow,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
            noexcept -> GpuTuple<Real>
        {
            if (i == (domlo.x + n) && (j >= domlo.y+n && j<= domhi.y-n) && bdatxlo(i,j,k) < zero) {
                return { std::abs(bdatxlo(i,j,k)) * area_x[box_no](i,j,k) };
            } else if ( (i == domhi.x+1 - n) && (j >= domlo.y+n && j<= domhi.y-n) && bdatxhi(i,j,k) > zero) {
                return { std::abs(bdatxhi(i,j,k)) * area_x[box_no](i,j,k) };
            } else {
                return { zero };
            }
        });

    // Y-dir
    const auto& bdatylo = bdy_data_ylo.const_array();
    const auto& bdatyhi = bdy_data_yhi.const_array();
    auto const& area_y = area_vec[1]->const_arrays();
    influx += ds_y *
        ParReduce(TypeList<ReduceOpSum>{},
                  TypeList<Real>{},
                  *area_vec[1], ngrow,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
            noexcept -> GpuTuple<Real>
        {
            if ( (j == domlo.y + n) && (i >= domlo.x+n && i<= domhi.x-n) && bdatylo(i,j,k) > zero) {
                return { std::abs(bdatylo(i,j,k)) * area_y[box_no](i,j,k) };
            } else if ( (j == domhi.y+1 - n) && (i >= domlo.x+n && i<= domhi.x-n) && bdatyhi(i,j,k) < zero) {
                return { std::abs(bdatyhi(i,j,k)) * area_y[box_no](i,j,k) };
            } else {
                return { zero };
            }
        });

    outflux += ds_y *
        ParReduce(TypeList<ReduceOpSum>{},
                  TypeList<Real>{},
                  *area_vec[1], ngrow,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
            noexcept -> GpuTuple<Real>
        {
            if ( (j == domlo.y + n) && (i >= domlo.x+n && i<= domhi.x-n) && bdatylo(i,j,k) < zero) {
                return { std::abs(bdatylo(i,j,k)) * area_y[box_no](i,j,k) };
            } else if ( (j == domhi.y+1 - n) && (i >= domlo.x+n && i<= domhi.x-n) && bdatyhi(i,j,k) > zero) {
                return { std::abs(bdatyhi(i,j,k)) * area_y[box_no](i,j,k) };
            } else {
                return { zero };
            }
        });

    ParallelDescriptor::ReduceRealSum(influx);
    ParallelDescriptor::ReduceRealSum(outflux);
}

void enforceInOutSolvability_bdy (const MultiFab& rho0,
                                  FArrayBox& bdy_data_xlo,
                                  FArrayBox& bdy_data_xhi,
                                  FArrayBox& bdy_data_ylo,
                                  FArrayBox& bdy_data_yhi,
                                  Array<MultiFab*, AMREX_SPACEDIM>& area_vec,
                                  const Geometry& geom,
                                  const Vector<BCRec>& domain_bcs_type_h)
{
    BL_PROFILE_VAR("enforceInOutSolvability_bdy()", enforceInOutSolvabilityBdy);

#ifdef AMREX_USE_FLOAT
    Real small_vel = Real(1.e-6);
#else
    Real small_vel = Real(1.e-8);
#endif

    Box domain(geom.Domain());

    bool multiply_by_rho0 = true;
    scale_bdy_normal_by_rho0(rho0, bdy_data_xlo, 0, domain, domain_bcs_type_h, multiply_by_rho0);
    scale_bdy_normal_by_rho0(rho0, bdy_data_xhi, 0, domain, domain_bcs_type_h, multiply_by_rho0);
    scale_bdy_normal_by_rho0(rho0, bdy_data_ylo, 1, domain, domain_bcs_type_h, multiply_by_rho0);
    scale_bdy_normal_by_rho0(rho0, bdy_data_yhi, 1, domain, domain_bcs_type_h, multiply_by_rho0);

    int width =bdy_data_xlo.box().length(0);
    AMREX_ASSERT(width == bdy_data_xhi.box().length(0));
    AMREX_ASSERT(width == bdy_data_ylo.box().length(1));
    AMREX_ASSERT(width == bdy_data_yhi.box().length(1));
    // amrex::Print() << "width in enforce solvability " << width << std::endl;

    for (int n = 0; n < width; n++)
    {
        Real influx = zero, outflux = zero;
        compute_influx_outflux_bdy(bdy_data_xlo, bdy_data_xhi,
                                   bdy_data_ylo, bdy_data_yhi,
                                   area_vec, geom, influx, outflux, n);

        // amrex::Print() << " TOTAL BDY INFLUX / OUTFLOW AT WIDTH " << n << " " << influx << " " << outflux << "\n";

        if ((influx > small_vel) && (outflux < small_vel)) {
            Abort("Cannot enforce solvability on bdy_data, no outflow from the direction dependent boundaries");
        } else if ((influx < small_vel) && (outflux < small_vel)) {
            // do nothing
        } else {
            const Real alpha_fcf = influx/outflux;  // flux correction factor

            correct_bdy_outflow_on_face(bdy_data_xlo, 0, true , domain, alpha_fcf, n);
            correct_bdy_outflow_on_face(bdy_data_xhi, 0, false, domain, alpha_fcf, n);
            correct_bdy_outflow_on_face(bdy_data_ylo, 1, true , domain, alpha_fcf, n);
            correct_bdy_outflow_on_face(bdy_data_yhi, 1, false, domain, alpha_fcf, n);

            // For diagnostic purposes.
            // Real influx_dbg = zero, outflux_dbg = zero;
            // compute_influx_outflux_bdy(bdy_data_xlo, bdy_data_xhi,
            //                            bdy_data_ylo, bdy_data_yhi,
            //                            area_vec, geom, influx_dbg, outflux_dbg, n);
            // amrex::Print() << " TOTAL BDY INFLUX / OUTFLOW AT WIDTH " << n << " " << influx_dbg << " " << outflux_dbg << "\n";
        }
    }

    multiply_by_rho0 = false;
    scale_bdy_normal_by_rho0(rho0, bdy_data_xlo, 0, domain, domain_bcs_type_h, multiply_by_rho0);
    scale_bdy_normal_by_rho0(rho0, bdy_data_xhi, 0, domain, domain_bcs_type_h, multiply_by_rho0);
    scale_bdy_normal_by_rho0(rho0, bdy_data_ylo, 1, domain, domain_bcs_type_h, multiply_by_rho0);
    scale_bdy_normal_by_rho0(rho0, bdy_data_yhi, 1, domain, domain_bcs_type_h, multiply_by_rho0);
}
