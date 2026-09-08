#ifndef ERF_SHOC_PREPROCESS_H_
#define ERF_SHOC_PREPROCESS_H_

#include "ERF_ShocColumnData.H"

#include <AMReX_MFIter.H>

class ShocPreprocess
{
public:
    static void fill_columns (ShocColumnData& col,
                              const amrex::MFIter& mfi,
                              const amrex::MultiFab& cons,
                              const amrex::MultiFab& xvel,
                              const amrex::MultiFab& yvel,
                              const amrex::MultiFab& zvel,
                              const amrex::MultiFab* hfx3,
                              const amrex::MultiFab* qfx3,
                              const amrex::MultiFab* tau13,
                              const amrex::MultiFab* tau23,
                              const amrex::MultiFab& z_phys_nd,
                              const amrex::Geometry& geom,
                              const MoistureComponentIndices& moisture_indices);
};

#endif
