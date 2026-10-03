#ifndef ERF_ParFunctions_H
#define ERF_ParFunctions_H

#include <limits>
#include <memory>

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

/**
 * Reduce a multifab of heights to a vector of max values at each height: v[k] is the largest
 * valid value of mf at index k over all boxes of the level.  A box contributes only the k it
 * holds, since the boxes of a level need not span the column (a refined region that does not
 * reach the top or the ground, or grids split in z).  A k that no box holds gets the nominal
 * height of the level from zlevels_stag (the face height if mf is nodal in z, the mean of the
 * two faces otherwise), so every entry of v is a height and no caller has to test for a k
 * that is not covered.
 *
 * @param[out] v            max of mf per k
 * @param[in]  mf           heights (z_phys_cc or z_phys_nd of a level)
 * @param[in]  zlevels_stag nominal heights of the z faces of the level
 */
AMREX_FORCE_INLINE
void reduce_to_max_per_height (amrex::Vector<amrex::Real>& v,
                               std::unique_ptr<amrex::MultiFab>& mf,
                               const amrex::Vector<amrex::Real>& zlevels_stag)
{
    const int  nk      = static_cast<int>(v.size());
    const bool nodal_z = mf->ixType().nodeCentered(2);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<int>(zlevels_stag.size()) >= (nodal_z ? nk : nk+1),
        "reduce_to_max_per_height: zlevels_stag is shorter than the column");

    // the identity of ReduceOpMax: the value of a k that no box holds
    constexpr amrex::Real not_covered = std::numeric_limits<amrex::Real>::lowest();

    for (int k = 0; k < nk; k++) {
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
        v[k] = amrex::get<0>(reduce_data.value(reduce_op));
    }
    amrex::ParallelDescriptor::ReduceRealMax(v.data(), nk);

    for (int k = 0; k < nk; k++) {
        if (v[k] == not_covered) {
            v[k] = nodal_z ? zlevels_stag[k]
                           : amrex::Real(0.5) * (zlevels_stag[k] + zlevels_stag[k+1]);
        }
    }
}

#endif /* ERF_ParFunctions.H */
