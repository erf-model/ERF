#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>

using namespace amrex;

/**
 * Define the domain boundary conditions for the (optional) Poisson solve
 */

using BCType = LinOpBCType;

/**
 * Solve the Poisson equation using EB_enabled MLMG
 * Note that the level may or may not be level 0.
 */
void ERF::solve_with_EB_mlmg (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi, Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes)
{
    BL_PROFILE("ERF::solve_with_EB_mlmg()");

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    LPInfo info;
    // Allow a hidden direction if the domain is one cell wide in any lateral direction
    if (dom_lo.x == dom_hi.x) {
        info.setHiddenDirection(0);
    } else if (dom_lo.y == dom_hi.y) {
        info.setHiddenDirection(1);
    }

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(rhs[0].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(rhs[0].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);

    // amrex::Print() << "BCLO " << bclo[0] << " " << bclo[1] << " " << bclo[2] << std::endl;
    // amrex::Print() << "BCHI " << bchi[0] << " " << bchi[1] << " " << bchi[2] << std::endl;

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    // ****************************************************************************
    // Multigrid solve
    // ****************************************************************************

    MLEBABecLap mleb (geom_tmp, ba_tmp, dm_tmp, info, {m_factory[lev].get()});

    mleb.setMaxOrder(2);
    mleb.setDomainBC(bclo, bchi);
    mleb.setLevelBC(0, nullptr);

    //
    // This sets A = 0, B = 1 so that
    // the operator A alpha - b del dot beta grad to b
    // becomes  - del dot beta grad
    //
    mleb.setScalars(0.0, 1.0);

    Array<MultiFab,AMREX_SPACEDIM> bcoef;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        bcoef[idim].define(convert(ba_tmp[0],IntVect::TheDimensionVector(idim)),
                            dm_tmp[0], 1, 0, MFInfo(), *m_factory[lev]);
        bcoef[idim].setVal(1.0);
    }
    mleb.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef));

    MLMG mlmg(mleb);

    int max_iter = 100;
    mlmg.setMaxIter(max_iter);
    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);

    mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), reltol, abstol);

    mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

    phi[0].FillBoundary(geom[lev].periodicity());

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& pp_arr  = phi[0].array(mfi);
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
        if (bx_lo.y == dom_lo.y) {
            auto bc_type = domain_bc_type[Orientation(1,Orientation::low)];
            if (bc_type == "Outflow" || bc_type == "Open") {
                ParallelFor(makeSlab(bx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j-1,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j-1,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_lo.z == dom_lo.z) {
            ParallelFor(makeSlab(bx,2,dom_lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k-1) = pp_arr(i,j,k);
            });
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
        if (bx_hi.y == dom_hi.y) {
            auto bc_type = domain_bc_type[Orientation(1,Orientation::high)];
            if (bc_type == "Outflow" || bc_type == "Open") {
                ParallelFor(makeSlab(bx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j+1,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j+1,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_hi.z == dom_hi.z) {
            ParallelFor(makeSlab(bx,2,dom_hi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k+1) = pp_arr(i,j,k);
            });
        }
    } // mfi

    // Now overwrite with periodic fill outside domain and fine-fine fill inside
    phi[0].FillBoundary(geom[lev].periodicity());
}
