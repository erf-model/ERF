#ifndef ERF_ParFunctions_H
#define ERF_ParFunctions_H

#include <limits>
#include <memory>

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

/**
 * Reduce a multifab to a vector of max values at each height: v[k] is the largest valid value
 * of mf at index k over all boxes of the level.  A box contributes only the k it holds, since
 * the boxes of a level need not span the column (a refined region that does not reach the top,
 * or grids split in z).  A k that no box holds gets std::numeric_limits<Real>::lowest().
 */
AMREX_FORCE_INLINE
void reduce_to_max_per_height (amrex::Vector<amrex::Real>& v,
                               std::unique_ptr<amrex::MultiFab>& mf)
{
    for (int k = 0; k < static_cast<int>(v.size()); k++) {
        amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
        amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (amrex::MFIter mfi(*mf); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.validbox();
            if (k < bx.smallEnd(2) || k > bx.bigEnd(2)) { continue; }
            bx.setSmall(2, k);
            bx.setBig(2, k);
            const amrex::Array4<amrex::Real const> a = mf->const_array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept -> ReduceTuple
            {
                return { a(i,j,kk) };
            });
        }
        // the identity of ReduceOpMax, lowest(), when no box of this rank holds k
        v[k] = amrex::get<0>(reduce_data.value(reduce_op));
    }
    amrex::ParallelDescriptor::ReduceRealMax(v.data(), static_cast<int>(v.size()));
}

#endif /* ERF_ParFunctions.H */
