#ifndef ERF_SHOC_MOMENTS_H_
#define ERF_SHOC_MOMENTS_H_

#include "ERF_ShocColumnData.H"

class ShocMoments
{
public:
    static void calc_var_or_covar (const ShocColumnData& col,
                                   amrex::Real tunefac,
                                   const amrex::FArrayBox& isotropy_zi,
                                   const amrex::FArrayBox& tkh_zi,
                                   const amrex::FArrayBox& invar1,
                                   const amrex::FArrayBox& invar2,
                                   amrex::FArrayBox& outvar);

    static void calc_vertflux (const ShocColumnData& col,
                               const amrex::FArrayBox& tkh_zi,
                               const amrex::FArrayBox& invar,
                               amrex::FArrayBox& vertflux);

    static void diagnose_second_moments (ShocColumnData& col,
                                         const ShocRuntimeOptions& opts);

    static void clip_third_moments (const ShocColumnData& col,
                                    const amrex::FArrayBox& w_sec_zi,
                                    amrex::FArrayBox& w3);

    static void diagnose_third_moments (ShocColumnData& col,
                                        const ShocRuntimeOptions& opts);

    static void diagnose_moments (ShocColumnData& col,
                                  const ShocRuntimeOptions& opts);
};

#endif
