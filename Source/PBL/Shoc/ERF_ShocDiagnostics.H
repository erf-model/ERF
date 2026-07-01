#ifndef ERF_SHOC_DIAGNOSTICS_H_
#define ERF_SHOC_DIAGNOSTICS_H_

#include "ERF_ShocColumnData.H"
#include "ERF_ShocMoments.H"
#include "ERF_ShocPDF.H"
#include "ERF_ShocStructure.H"
#include "ERF_ShocTKE.H"

class ShocDiagnostics
{
public:
    static void diagnose_pre_implicit (ShocColumnData& col,
                                       const ShocRuntimeOptions& opts,
                                       amrex::Real dx,
                                       amrex::Real dy,
                                       amrex::Real dt);

    static void diagnose_post_implicit (ShocColumnData& col,
                                        const ShocRuntimeOptions& opts,
                                        amrex::Real dt);
};

#endif
