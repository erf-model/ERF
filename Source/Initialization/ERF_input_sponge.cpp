/**
 * \file ERF_input_sponge.cpp
 */

#include <ERF.H>
#include <EOS.H>
#include <ERF_Constants.H>
#include <Utils.H>
#include <prob_common.H>

using namespace amrex;

/**
 * High level wrapper for sponge x and y velocities
 * level data from input sponge data.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::input_sponge (int lev)
{
    // We only want to read the file once
    if (lev == 0) {
        if (input_sponge_file.empty())
            Error("input_sounding file name must be provided via input");

        // this will interpolate the input profiles to the nominal height levels
        // (ranging from 0 to the domain top)
        input_sponge_data.read_from_file(input_sponge_file, geom[lev], zlevels_stag);
    }
}
