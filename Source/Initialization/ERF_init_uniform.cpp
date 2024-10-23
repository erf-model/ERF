/**
 * \file ERF_init_uniform.cpp
 */

#include <ERF.H>
#include <ERF_prob_common.H>

using namespace amrex;

/**
 * Use problem-specific reference density and temperature to set the
 * background state to a uniform value.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_uniform (int lev)
{
    auto& lev_new = vars_new[lev];
    for (MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box &gbx = mfi.growntilebox(1);
        const auto &cons_arr = lev_new[Vars::cons].array(mfi);
        prob->init_uniform(gbx, cons_arr);
    }
}
