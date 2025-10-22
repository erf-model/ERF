#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

/**
 * Impose bc's on the pressure that comes out of the solve
 */
void ERF::ImposeBCsOnPhi (int lev, MultiFab& phi, const Box& subdomain)
{
    BL_PROFILE("ERF::ImposeBCsOnPhi()");

    auto const sub_lo = lbound(subdomain);
    auto const sub_hi = ubound(subdomain);

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    phi.setBndry(1.e25);
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

        auto bc_type_xlo = domain_bc_type[Orientation(0,Orientation::low)];
        auto bc_type_xhi = domain_bc_type[Orientation(0,Orientation::high)];
        auto bc_type_ylo = domain_bc_type[Orientation(1,Orientation::low)];
        auto bc_type_yhi = domain_bc_type[Orientation(1,Orientation::high)];
        auto bc_type_zhi = domain_bc_type[Orientation(2,Orientation::high)];

        if ( (bx_lo.x == dom_lo.x) && (bc_type_xlo == "Outflow" || bc_type_xlo == "Open") && !solverChoice.use_real_bcs) {
            ParallelFor(makeSlab(bx,0,dom_lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i-1,j,k) = -pp_arr(i,j,k);
            });
        } else if (bx_lo.x == sub_lo.x) {
            ParallelFor(makeSlab(bx,0,sub_lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i-1,j,k) = pp_arr(i,j,k);
            });
        }

        if ( (bx_hi.x == dom_hi.x) && (bc_type_xhi == "Outflow" || bc_type_xhi == "Open") && !solverChoice.use_real_bcs) {
            ParallelFor(makeSlab(bx,0,dom_hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i+1,j,k) = -pp_arr(i,j,k);
            });
        } else if (bx_hi.x == sub_hi.x) {
            ParallelFor(makeSlab(bx,0,sub_hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i+1,j,k) = pp_arr(i,j,k);
            });
        }

        if ( (bx_lo.y == dom_lo.y) && (bc_type_ylo == "Outflow" || bc_type_ylo == "Open") && !solverChoice.use_real_bcs) {
            ParallelFor(makeSlab(bx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j-1,k) = -pp_arr(i,j,k);
            });
        } else if (bx_lo.y == sub_lo.y) {
            ParallelFor(makeSlab(bx,1,sub_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j-1,k) = pp_arr(i,j,k);
            });
        }

        if ( (bx_hi.y == dom_hi.y) && (bc_type_yhi == "Outflow" || bc_type_yhi == "Open") && !solverChoice.use_real_bcs) {
            ParallelFor(makeSlab(bx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j+1,k) = -pp_arr(i,j,k);
            });
        } else if (bx_hi.y == sub_hi.y) {
            ParallelFor(makeSlab(bx,1,sub_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j+1,k) = pp_arr(i,j,k);
            });
        }

        // At low z we are always Neumann whether the box touches the bottom boundary or not
        Box zbx(bx); zbx.grow(0,1); zbx.grow(1,1); // Grow in x-dir and y-dir because we have filled that above
        if (bx_lo.z == sub_lo.z) {
            ParallelFor(makeSlab(zbx,2,dom_lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k-1) = pp_arr(i,j,k);
            });
        }

        if ( (bx_hi.z == dom_hi.z) && (bc_type_zhi == "Outflow" || bc_type_zhi == "Open") ) {
            ParallelFor(makeSlab(bx,2,dom_hi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k+1) = -pp_arr(i,j,k);
            });
        } else if (bx_hi.z == sub_hi.z) {
            ParallelFor(makeSlab(bx,2,sub_hi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k+1) = pp_arr(i,j,k);
            });
        }
    } // mfi

    // Now overwrite with periodic fill outside domain and fine-fine fill inside
    phi.FillBoundary(geom[lev].periodicity());
}
