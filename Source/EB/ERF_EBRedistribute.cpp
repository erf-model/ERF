
#include <AMReX_Config.H>
#include <AMReX_Geometry.H>

#include <ERF.H>
#include <ERF_EB.H>
#include <ERF_EBRedistribute.H>
#include "AMReX_EB_Redistribution.H"

using namespace amrex;

void
redistribute_term ( int ncomp,
                    const Geometry& geom,
                    MultiFab& result,
                    MultiFab& result_tmp,
                    MultiFab const& state,
                    EBFArrayBoxFactory const& ebfact,
                    BCRec const* bc, // this is bc for the state (needed for SRD slopes)
                    Real const local_dt)
{
    BL_PROFILE_VAR("redistribute_term1", redistribute_term1);
    // ************************************************************************
    // Redistribute result_tmp and pass out result
    // ************************************************************************
    AMREX_ASSERT(result.nComp() == state.nComp());

    result_tmp.FillBoundary(geom.periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flag = flagfab.const_array();

        bool regular = (flagfab.getType(amrex::grow(bx,4)) == FabType::regular);
        bool covered = (flagfab.getType(bx) == FabType::covered);

        Array4<Real> out = result.array(mfi);
        Array4<Real> in  = result_tmp.array(mfi);

        if (!regular && !covered)
        {
            auto const& vfrac = ebfact.getVolFrac().const_array(mfi);
            auto const& ccc   = ebfact.getCentroid().const_array(mfi);

            auto const& apx = ebfact.getAreaFrac()[0]->const_array(mfi);
            auto const& apy = ebfact.getAreaFrac()[1]->const_array(mfi);
            auto const& apz = ebfact.getAreaFrac()[2]->const_array(mfi);

            auto const& fcx = ebfact.getFaceCent()[0]->const_array(mfi);
            auto const& fcy = ebfact.getFaceCent()[1]->const_array(mfi);
            auto const& fcz = ebfact.getFaceCent()[2]->const_array(mfi);

            Box gbx = bx; gbx.grow(3);

            FArrayBox scratch_fab(gbx,ncomp);
            Array4<Real> scratch = scratch_fab.array();
            Elixir eli_scratch = scratch_fab.elixir();

            std::string redistribution_type = "StateRedist";

            // State redist acts on the state.
            Array4<Real const> state_arr = state.const_array(mfi);
            ApplyRedistribution( bx, ncomp, out, in, state_arr,
                                 scratch, flag,
                                 apx, apy, apz, vfrac,
                                 fcx, fcy, fcz, ccc,
                                 bc, geom, local_dt, redistribution_type,
                                 false, 2, 0.5_rt, {});
        }
        else
        {
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                out(i,j,k,n) = in(i,j,k,n);
            });
        }
    } // MFIter
}

void
redistribute_term ( int ncomp,
                    const Geometry& geom,
                    MultiFab& result,
                    MultiFab& result_tmp,
                    MultiFab const& state,
                    eb_aux_ const& ebfact,
                    BCRec const* bc, // this is bc for the state (needed for SRD slopes)
                    Real const local_dt,
                    int const igrid)
{
    BL_PROFILE_VAR("redistribute_term2", redistribute_term2);
    // ************************************************************************
    // Redistribute result_tmp and pass out result
    // ************************************************************************
    AMREX_ASSERT(result.nComp() == state.nComp());

    result_tmp.FillBoundary(geom.periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flag = flagfab.const_array();

        bool regular = (flagfab.getType(amrex::grow(bx,4)) == FabType::regular);
        bool covered = (flagfab.getType(bx) == FabType::covered);

        Array4<Real> out = result.array(mfi);
        Array4<Real> in  = result_tmp.array(mfi);

        if (!regular && !covered)
        {
            auto const& vfrac = ebfact.getVolFrac().const_array(mfi);
            auto const& ccc   = ebfact.getCentroid().const_array(mfi);

            auto const& apx = ebfact.getAreaFrac()[0]->const_array(mfi);
            auto const& apy = ebfact.getAreaFrac()[1]->const_array(mfi);
            auto const& apz = ebfact.getAreaFrac()[2]->const_array(mfi);

            auto const& fcx = ebfact.getFaceCent()[0]->const_array(mfi);
            auto const& fcy = ebfact.getFaceCent()[1]->const_array(mfi);
            auto const& fcz = ebfact.getFaceCent()[2]->const_array(mfi);

            Box bx_cc = bx;
            bx_cc = bx_cc.convert(IntVect::TheZeroVector());

            // Extend box for staggered grids
            if (igrid == IntVars::xmom) bx_cc = bx_cc.setBig(0, bx_cc.bigEnd(0) + 1);
            if (igrid == IntVars::ymom) bx_cc = bx_cc.setBig(1, bx_cc.bigEnd(1) + 1);
            if (igrid == IntVars::zmom) bx_cc = bx_cc.setBig(2, bx_cc.bigEnd(2) + 1);

            Box gbx = bx_cc; gbx.grow(3);

            // Extended geometry domain
            Box domain_grown = geom.Domain();
            domain_grown.grow(igrid-1, 1); // Extend geometry domain by 1 in the staggering direction
            Geometry geom_new(domain_grown, geom.ProbDomain(), geom.Coord(), geom.isPeriodic());

            FArrayBox scratch_fab(gbx,ncomp);
            Array4<Real> scratch = scratch_fab.array();
            Elixir eli_scratch = scratch_fab.elixir();

            // This is scratch space if calling StateRedistribute

            std::string redistribution_type = "StateRedist";

            // State redist acts on the state.
            Array4<Real const> state_arr = state.const_array(mfi);
            ApplyRedistribution( bx_cc, ncomp, out, in, state_arr,
                                 scratch, flag,
                                 apx, apy, apz, vfrac,
                                 fcx, fcy, fcz, ccc,
                                 bc, geom_new, local_dt, redistribution_type,
                                 false, 2, 0.5_rt, {});
        }
        else
        {
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                out(i,j,k,n) = in(i,j,k,n);
            });
        }
    } // MFIter
}
