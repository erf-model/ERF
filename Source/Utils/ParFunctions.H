#ifndef ParFunctions_H
#define ParFunctions_H

/**
 * Reduce a multifab to a vector of values at height levels
 */
AMREX_FORCE_INLINE
void reduce_to_max_per_level (amrex::Vector<amrex::Real>& v,
                              std::unique_ptr<amrex::MultiFab>& mf)
{
    amrex::MultiArray4<const amrex::Real> const& ma = mf->const_arrays();
    for (int k = 0; k < v.size(); k++) {
        v[k] = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpMax>{},
                                amrex::TypeList<amrex::Real>{},
                                *mf,
                                amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
            -> amrex::GpuTuple<amrex::Real>
        {
            return { ma[box_no](i,j,k) };
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(v.data(), v.size());
}

#endif /* ParFunctions_H */
