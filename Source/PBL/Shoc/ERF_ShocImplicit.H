#ifndef ERF_SHOC_IMPLICIT_H_
#define ERF_SHOC_IMPLICIT_H_

#include "ERF_ShocColumnData.H"

class ShocImplicit
{
public:
    static void cache_baseline_state (ShocColumnData& col);

    static void compute_tmpi (const ShocColumnData& col,
                              int ic,
                              amrex::Real dt,
                              const amrex::Vector<amrex::Real>& rho_zi,
                              amrex::Vector<amrex::Real>& tmpi);

    static void compute_dp_inverse (const ShocColumnData& col,
                                    int ic,
                                    amrex::Vector<amrex::Real>& rdp_zt);

    AMREX_GPU_HOST_DEVICE
    static amrex::Real compute_temperature (amrex::Real thetal,
                                            amrex::Real ql,
                                            amrex::Real exner);

    AMREX_GPU_HOST_DEVICE
    static amrex::Real compute_vapor (amrex::Real qw,
                                      amrex::Real ql);

    static void advance_implicit_state (ShocColumnData& col,
                                        const ShocRuntimeOptions& opts,
                                        amrex::Real dt);

    static void finalize_from_pdf (ShocColumnData& col,
                                   const ShocRuntimeOptions& opts,
                                   amrex::Real dt);

    static void update_prognostics (ShocColumnData& col,
                                    const ShocRuntimeOptions& opts,
                                    amrex::Real dt);
};

#endif
