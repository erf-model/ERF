#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

/**
 * Impose bc's on the pressure that comes out of the solve
 */
void ERF::ImposeBCsOnPhi (int lev, MultiFab& phi)
{
    BL_PROFILE("ERF::ImposeBCsOnPhi()");

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    phi.FillBoundary(geom[lev].periodicity());

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& pp_arr  = phi.array(mfi);
        Box const& bx    = mfi.tilebox();
        auto const bx_lo = lbound(bx);
        auto const bx_hi = ubound(bx);

        if (bx_lo.x == dom_lo.x) {
            auto bc_type = domain_bc_type[Orientation(0,Orientation::low)];
            if (bc_type == "Outflow" || bc_type == "Open") {
                ParallelFor(makeSlab(bx,0,dom_lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i-1,j,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,0,dom_lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i-1,j,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_hi.x == dom_hi.x) {
            auto bc_type = domain_bc_type[Orientation(0,Orientation::high)];
            if (bc_type == "Outflow" || bc_type == "Open") {
                ParallelFor(makeSlab(bx,0,dom_hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i+1,j,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,0,dom_hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i+1,j,k) = pp_arr(i,j,k);
                });
            }
        }

        if (bx_lo.y == dom_lo.y) {
            auto bc_type = domain_bc_type[Orientation(1,Orientation::low)];
            Box ybx(bx); ybx.grow(0,1); // Grow in x-dir because we have filled that above
            if (bc_type == "Outflow" || bc_type == "Open") {
                ParallelFor(makeSlab(ybx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j-1,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(ybx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j-1,k) = pp_arr(i,j,k);
                });
            }
        }

        if (bx_hi.y == dom_hi.y) {
            auto bc_type = domain_bc_type[Orientation(1,Orientation::high)];
            Box ybx(bx); ybx.grow(0,1); // Grow in x-dir because we have filled that above
            if (bc_type == "Outflow" || bc_type == "Open") {
                ParallelFor(makeSlab(ybx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j+1,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(ybx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j+1,k) = pp_arr(i,j,k);
                });
            }
        }

        Box zbx(bx); zbx.grow(0,1); zbx.grow(1,1); // Grow in x-dir and y-dir because we have filled that above
        if (bx_lo.z == dom_lo.z) {
            ParallelFor(makeSlab(zbx,2,dom_lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k-1) = pp_arr(i,j,k);
            });
        }

        auto zbc_type = domain_bc_type[Orientation(2,Orientation::high)];
        if (bx_hi.z == dom_hi.z) {
            if (zbc_type == "Outflow" || zbc_type == "Open") {
                ParallelFor(makeSlab(zbx,2,dom_hi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j,k+1) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(zbx,2,dom_hi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j,k+1) = pp_arr(i,j,k);
                });
            }
        }
    } // mfi

    // Now overwrite with periodic fill outside domain and fine-fine fill inside
    phi.FillBoundary(geom[lev].periodicity());
}
