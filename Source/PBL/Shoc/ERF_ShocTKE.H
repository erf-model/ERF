#ifndef ERF_SHOC_TKE_H_
#define ERF_SHOC_TKE_H_

#include "ERF_ShocColumnData.H"

class ShocTKE
{
public:
    static void compute_shear_production (const ShocColumnData& col,
                                          amrex::FArrayBox& sterm_iface);

    static void integrate_column_stability (const ShocColumnData& col,
                                            amrex::Vector<amrex::Real>& brunt_int);

    static void diagnose_tke_and_diffusivities (ShocColumnData& col,
                                                const ShocRuntimeOptions& opts,
                                                amrex::Real dt);
};

#endif
