
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
     * noahmpio.nsoil = 1;
     * noahmpio.nsnow = 3;
     *
     * noahmpio.ids = 1;
     * noahmpio.ide = 1;
     * noahmpio.ims = 1;
     * noahmpio.ime = 1;
     * noahmpio.its = 1;
     * noahmpio.ite = 1;
     *
     * noahmpio.jds = 1;
     * noahmpio.jde = 1;
     * noahmpio.jms = 1;
     * noahmpio.jme = 1;
     * noahmpio.jts = 1;
     * noahmpio.jte = 1;
     *
     * noahmpio.kds = 1;
     * noahmpio.kde = 1;
     * noahmpio.kms = 1;
     * noahmpio.kme = 1;
     * noahmpio.kts = 1;
     * noahmpio.kte = 1;
     */

    NoahmpIOVarInitDefault(&noahmpio);
    NoahmpInitMain(&noahmpio);

    amrex::Print() << "Noahmp IO Initialized" << std::endl;
};
