#include "ERF_ShocDiagnostics.H"

#include <AMReX_BLProfiler.H>

void
ShocDiagnostics::diagnose_pre_implicit (ShocColumnData& col,
                                        const ShocRuntimeOptions& opts,
                                        amrex::Real dx,
                                        amrex::Real dy,
                                        amrex::Real dt)
{
    {
        BL_PROFILE("SHOC::advance::structure");
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
        ShocStructure::diagnose_length_and_brunt(col, opts, dx, dy);
    }
    {
        BL_PROFILE("SHOC::advance::tke");
        ShocTKE::diagnose_tke_and_diffusivities(col, opts, dt);
    }
}

void
ShocDiagnostics::diagnose_post_implicit (ShocColumnData& col,
                                         const ShocRuntimeOptions& opts,
                                         amrex::Real dt)
{
    {
        BL_PROFILE("SHOC::advance::moments");
        ShocMoments::diagnose_moments(col, opts);
    }
    {
        BL_PROFILE("SHOC::advance::pdf");
        ShocPDF::diagnose_pdf(col, opts, dt);
    }
}
