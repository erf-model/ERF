/**
 * \file ERF_EBRedistribute.cpp
 * \brief Implements EB redistribution for native and auxiliary EB factories.
 */

#include <AMReX_Config.H>
#include <AMReX_Geometry.H>

#include <ERF.H>
#include <ERF_EB.H>
#include <ERF_EBRedistribute.H>
#include "AMReX_EB_Redistribution.H"

using namespace amrex;

/**
 * \brief Apply EB state redistribution to result_tmp and write the redistributed result.
 */
template <typename EBFactType>
void
redistribute_term ( int ncomp,
                    const Geometry& geom,
                    MultiFab& result,
                    MultiFab& result_tmp,
                    MultiFab const& state,
                    EBFactType const& ebfact,
                    BCRec const* bc, // this is bc for the state (needed for SRD slopes)
                    double local_dt_d,
                    int const igrid)
{
    BL_PROFILE_VAR("redistribute_term", redistribute_term);
    // ************************************************************************
    // Redistribute result_tmp and pass out result
    // ************************************************************************
    AMREX_ASSERT(result.nComp() == state.nComp());

    Real local_dt = static_cast<Real>(local_dt_d);

    result_tmp.FillBoundary(geom.periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flag = flagfab.const_array();

        // Redistribution operates on a grown box and needs EB geometry data
        bool is_singlevalued = (flagfab.getType() == FabType::singlevalued);
        bool is_regular = (flagfab.getType(amrex::grow(bx,4)) == FabType::regular);
        bool is_covered = (flagfab.getType(bx) == FabType::covered);

        Array4<Real> out = result.array(mfi);
        Array4<Real> in  = result_tmp.array(mfi);

        // Only apply redistribution if we have actual cut cells
        if (is_singlevalued && !is_regular && !is_covered)
        {
            auto const& vfrac = ebfact.getVolFrac().const_array(mfi);
            auto const& ccc = ebfact.getCentroid().const_array(mfi);
            auto const& apx = ebfact.getAreaFrac()[0]->const_array(mfi);
            auto const& apy = ebfact.getAreaFrac()[1]->const_array(mfi);
            auto const& apz = ebfact.getAreaFrac()[2]->const_array(mfi);
            auto const& fcx = ebfact.getFaceCent()[0]->const_array(mfi);
            auto const& fcy = ebfact.getFaceCent()[1]->const_array(mfi);
            auto const& fcz = ebfact.getFaceCent()[2]->const_array(mfi);

            Box bx_cc = bx;
            Geometry geom_used = geom;

            // Handle staggered grids if igrid is specified
            if (igrid >= 0) {
                bx_cc = bx_cc.convert(IntVect::TheZeroVector());

                // Extend box for staggered grids
                if (igrid == IntVars::xmom) bx_cc = bx_cc.setBig(0, bx_cc.bigEnd(0) + 1);
                if (igrid == IntVars::ymom) bx_cc = bx_cc.setBig(1, bx_cc.bigEnd(1) + 1);
                if (igrid == IntVars::zmom) bx_cc = bx_cc.setBig(2, bx_cc.bigEnd(2) + 1);

                // Extended geometry domain
                Box domain_grown = geom.Domain();
                domain_grown.growHi(igrid-1, 1); // Extend geometry domain by 1 in the staggering direction
                geom_used = Geometry(domain_grown, geom.ProbDomain(), geom.Coord(), geom.isPeriodic());
            }

            Box gbx = bx_cc; gbx.grow(3);

            FArrayBox scratch_fab(gbx,ncomp);
            Array4<Real> scratch = scratch_fab.array();
            Elixir eli_scratch = scratch_fab.elixir();

            std::string redistribution_type = "StateRedist";

            // State redist acts on the state.
            Array4<Real const> state_arr = state.const_array(mfi);
            ApplyRedistribution( bx_cc, ncomp, out, in, state_arr,
                                 scratch, flag,
                                 apx, apy, apz, vfrac,
                                 fcx, fcy, fcz, ccc,
                                 bc, geom_used, local_dt, redistribution_type,
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

//! \cond
// Explicit template instantiations for the types we use
template void redistribute_term(int, const Geometry&, MultiFab&, MultiFab&,
                                MultiFab const&, EBFArrayBoxFactory const&,
                                BCRec const*, double, int const);
template void redistribute_term(int, const Geometry&, MultiFab&, MultiFab&,
                                MultiFab const&, eb_aux_ const&,
                                BCRec const*, double, int const);
//! \endcond
