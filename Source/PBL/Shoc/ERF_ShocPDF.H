#ifndef ERF_SHOC_PDF_H_
#define ERF_SHOC_PDF_H_

#include "ERF_ShocColumnData.H"

class ShocPDF
{
public:
    static void diagnose_pdf (ShocColumnData& col,
                              const ShocRuntimeOptions& opts,
                              amrex::Real dt);
};

#endif
