
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

    NoahmpIOVarInit(&noahmp_io, &noahmp_arr);
    NoahmpInitMain(&noahmp_io);

    amrex::Print() << "Noahmp IO Initialized" << std::endl;
};
