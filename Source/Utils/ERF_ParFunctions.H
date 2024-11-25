#ifndef ERF_ParFunctions_H
#define ERF_ParFunctions_H

/**
 * Reduce a multifab to a vector of max values at each height
 */
AMREX_FORCE_INLINE
void reduce_to_max_per_height (amrex::Vector<amrex::Real>& v,
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

// Partial Sum over specific components.
// This is an interface for the micro models.
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
comp_partial_sum (const int& i,
                  const int& j,
                  const int& k,
                  const int& nsize,
                  const int* ind_map,
                  const amrex::Array4<const amrex::Real>& data)
{
    amrex::Real partial_sum = 0.0;
    for (int n(0); n<nsize; ++n) {
        int ind = RhoQ1_comp + ind_map[n];
        partial_sum += data(i,j,k,ind);
    }
    return partial_sum;
}
#endif /* ERF_ParFunctions.H */
