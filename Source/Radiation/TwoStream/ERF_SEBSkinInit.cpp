#include <limits>

#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_Reduce.H>

#include "ERF_SEBLandMask.H"
#include "ERF_SEBSkinInit.H"

using namespace amrex;

SEBSkinInitReport
seb_init_skin_from_file (MultiFab& t_sfc, const MultiFab& file_skin,
                         const iMultiFab* land_mask, const Real t_lo, const Real t_hi)
{
    if (t_sfc.boxArray() != file_skin.boxArray() ||
        t_sfc.DistributionMap() != file_skin.DistributionMap()) {
        Abort("seb_init_skin_from_file: the file's skin temperature is not on the grids of "
              "the balance's skin temperature");
    }
    seb_assert_land_mask_layout(land_mask, t_sfc, "seb_init_skin_from_file");

    // One pass: the copy and the counts. A refused column is left as it was, so the
    // caller's message is the only effect of a bad value.
    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_op;
    ReduceData<Long, Long, Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    const Real huge = std::numeric_limits<Real>::max();

    for (MFIter mfi(t_sfc, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& skin = t_sfc.array(mfi);
        const Array4<const Real>& file = file_skin.const_array(mfi);
        const SEBLandMaskView mask = seb_land_mask_view(land_mask, mfi);
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            if (mask.water(i, j)) { return {Long(0), Long(0), huge, -huge}; }
            const Real t = file(i, j, k);
            if (!seb_skin_init_usable(t, t_lo, t_hi)) { return {Long(0), Long(1), huge, -huge}; }
            skin(i, j, k) = t;
            return {Long(1), Long(0), t, t};
        });
    }

    const ReduceTuple result = reduce_data.value(reduce_op);
    SEBSkinInitReport report;
    report.n_set = amrex::get<0>(result);
    report.n_refused = amrex::get<1>(result);
    report.t_min = amrex::get<2>(result);
    report.t_max = amrex::get<3>(result);
    return report;
}
