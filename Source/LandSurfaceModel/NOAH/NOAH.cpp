#include <NOAH.H>
#include <NoahmpIO.H>

using namespace amrex;

/* Initialize lsm data structures */
void
NOAH::Init (const MultiFab& cons_in,
           const Geometry& geom,
           const Real& dt)
{ 
    // Initialize Noahmp IO
    //NoahmpIOVarInitDefault(&noahmpio);
};
