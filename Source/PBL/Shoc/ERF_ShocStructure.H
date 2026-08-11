#ifndef ERF_SHOC_STRUCTURE_H_
#define ERF_SHOC_STRUCTURE_H_

#include "ERF_ShocColumnData.H"

class ShocStructure
{
public:
    static void diagnose_surface_layer (ShocColumnData& col);
    static void diagnose_pblh (ShocColumnData& col);
    static void diagnose_length_and_brunt (ShocColumnData& col,
                                           const ShocRuntimeOptions& opts,
                                           amrex::Real dx,
                                           amrex::Real dy);
};

#endif
