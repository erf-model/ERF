
#include<iostream>

#include <AMReX_Print.H>
#include <ERF_NOAH.H>

using namespace amrex;

/* Initialize lsm data structures */
void
NOAH::Init (const MultiFab& cons_in,
            const Geometry& geom,
            const Real& dt)
{ 
    // Initialize Noahmp IO
    amrex::Print() << "Initializing Noahmp IO" << std::endl;

    /*
     * noahmpio.xstart = 1;
     * noahmpio.xend = 4;
     * noahmpio.ystart = 1;
     * noahmpio.yend = 2;
     *
     */

    NoahmpIOVarInitDefault(&noahmpio);
    NoahmpInitMain(&noahmpio);

    amrex::Print() << "Noahmp IO Initialized" << std::endl;
};
