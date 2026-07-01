#include "ERF.H"

using namespace amrex;

void
ERF::compute_native_shoc_tendencies (int lev,
                                     MultiFab* cons,
                                     MultiFab* xvel,
                                     MultiFab* yvel,
                                     MultiFab* zvel,
                                     Real* /*w_subsid*/,
                                     MultiFab* tau13,
                                     MultiFab* tau23,
                                     MultiFab* hfx3,
                                     MultiFab* qfx3,
                                     MultiFab* eddyDiffs,
                                     MultiFab* z_phys_nd_in,
                                     const Real& dt_advance)
{
    AMREX_ALWAYS_ASSERT(native_shoc_driver[lev]);
    native_shoc_driver[lev]->advance(*cons, *xvel, *yvel, *zvel,
                                     tau13, tau23, hfx3, qfx3, eddyDiffs,
                                     *z_phys_nd_in, Geom(lev), dt_advance);
}
