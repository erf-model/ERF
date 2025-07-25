
#include<iostream>
#include<string>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include <ERF_NOAH.H>

using namespace amrex;

/* Initialize lsm data structures */
void
NOAH::Init (const int& lev,
            const MultiFab& cons_in,
            const Geometry& geom,
            const Real& dt)
{

    m_dt = dt;
    m_geom = geom;

    Box domain = geom.Domain();
    khi_lsm    = domain.smallEnd(2) - 1;

    LsmVarMap.resize(m_lsm_size);
    LsmVarMap = {LsmVar_NOAH::t_sfc, LsmVar_NOAH::emis_sfc,
                 LsmVar_NOAH::alb_dir_vis, LsmVar_NOAH::alb_dir_nir,
                 LsmVar_NOAH::alb_dif_vis, LsmVar_NOAH::alb_dif_nir,
                 LsmVar_NOAH::sw_flux_dn , LsmVar_NOAH::lw_flux_dn };

    LsmVarName.resize(m_lsm_size);
    LsmVarName = {"t_sfc"      , "sfc_emis"   ,
                  "sfc_alb_dir_vis", "sfc_alb_dir_nir",
                  "sfc_alb_dif_vis", "sfc_alb_dif_nir",
                  "sw_flux_dn" , "lw_flux_dn" };

    // NOTE: lsm data is not used for Noahmp, however, the initialization is done
    //       to maintin consistency with IO and Driver interfaces that depend on
    //       this data. We eventually want to tweak those interfaces so we don't
    //       have to allocate lsm_data while using Noahmp lsm.

    // NOTE: All boxes in ba extend from zlo to zhi, so this transform is valid.
    //       If that were to change, the dm and new ba are no longer valid and
    //       direct copying between lsm data/flux vars cannot be done in a parfor.

    // Set 2D box array for lsm data
    IntVect ng(0,0,0);
    BoxArray ba = cons_in.boxArray();
    DistributionMapping dm = cons_in.DistributionMap();
    BoxList bl_lsm = ba.boxList();
    for (auto& b : bl_lsm) {
        b.setRange(2,0)
    }
    BoxArray ba_lsm(std::move(bl_lsm));

    // Set up lsm geometry
    const RealBox& dom_rb = m_geom.ProbDomain();
    const Real*    dom_dx = m_geom.CellSize();
    RealBox lsm_rb = dom_rb;
    Real lsm_dx[AMREX_SPACEDIM] = {AMREX_D_DECL(dom_dx[0],dom_dx[1],m_dz_lsm)};
    Real lsm_z_hi = dom_rb.lo(2);
    Real lsm_z_lo = lsm_z_hi - Real(m_nz_lsm)*lsm_dx[2];
    lsm_rb.setHi(2,lsm_z_hi); lsm_rb.setLo(2,lsm_z_lo);
    m_lsm_geom.define( ba_lsm.minimalBox(), lsm_rb, m_geom.Coord(), m_geom.isPeriodic() );

    // Create the data and fluxes
    for (auto ivar = 0; ivar < LsmVar_NOAH::NumVars; ++ivar) {
        // State vars are CC
        lsm_fab_vars[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_vars[ivar]->setVal(0.0);

        // Fluxes are nodal in z
        lsm_fab_flux[ivar] = std::make_shared<MultiFab>(convert(ba_lsm, IntVect(0,0,1)), dm, 1, IntVect(0,0,0));
        lsm_fab_flux[ivar]->setVal(0.);
    }

    // NOTE: Actual NoahmpIO interface that is relevant for the
    //       implementation of this lsm

    amrex::Print() << "Noah-MP initialization started" << std::endl;

    // Set noahmpio_vect to the size of local blocks (boxes)
    noahmpio_vect.resize(cons_in.local_size(), lev);

    // Iterate over multifab and noahmpio object together. Multifabs is
    // used to extract size of blocks and set bounds for noahmpio objects.
    int idb = 0;
    for (amrex::MFIter mfi(cons_in, false); mfi.isValid(); ++mfi, ++idb) {

        // Get bounds for the tile
        const amrex::Box& bx = mfi.tilebox();

        // Check if tile is at the lower boundary in lower z direction
        if (bx.smallEnd(2) == domain.smallEnd(2)) {

            // Get reference to the noahmpio object
            NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

            // Pass idb context to noahmpio
            noahmpio->blkid = idb;

            // Pass level context to noahmpio
            noahmpio->level = lev;

            // Initialize scalar values
            noahmpio->ScalarInitDefault();

            // Store the rank of process for noahmp
            noahmpio->rank = amrex::ParallelDescriptor::MyProc();

            // Store parallel communicator for noahmp
            noahmpio->comm = MPI_Comm_c2f(amrex::ParallelDescriptor::Communicator());

            // Read namelist.erf file. This file contains
            // noahmpio specific parameters and is read by
            // the Fortran side of the implementation.
            noahmpio->ReadNamelist();

            // Read the headers from the NetCDF land file. This is also
            // implemented on the Fortran side of things currently.
            noahmpio->ReadLandHeader();

            // Extract tile bounds and set them to their corresponding
            // noahmpio variables. At present we will set all the variables
            // corresponding to domain, memory, and tile to the same bounds.
            // This will be changed later if we want to do special memory
            // management for expensive use cases.
            noahmpio->xstart = bx.smallEnd(0);
            noahmpio->xend = bx.bigEnd(0);
            noahmpio->ystart = bx.smallEnd(1);
            noahmpio->yend = bx.bigEnd(1);

            // Domain bounds
            noahmpio->ids = noahmpio->xstart;
            noahmpio->ide = noahmpio->xend;
            noahmpio->jds = noahmpio->ystart;
            noahmpio->jde = noahmpio->yend;
            noahmpio->kds = 1;
            noahmpio->kde = 2;

            // Tile bounds
            noahmpio->its = noahmpio->xstart;
            noahmpio->ite = noahmpio->xend;
            noahmpio->jts = noahmpio->ystart;
            noahmpio->jte = noahmpio->yend;
            noahmpio->kts = 1;
            noahmpio->kte = 2;

            // Memory bounds
            noahmpio->ims = noahmpio->xstart;
            noahmpio->ime = noahmpio->xend;
            noahmpio->jms = noahmpio->ystart;
            noahmpio->jme = noahmpio->yend;
            noahmpio->kms = 1;
            noahmpio->kme = 2;

            // This procedure allocates memory in Fortran for IO variables
            // using bounds that are set above and read from namelist.erf
            // and headers from the NetCDF land file
            noahmpio->VarInitDefault();

            // This reads NoahmpTable.TBL file which is another input file
            // we need to set some IO variables.
            noahmpio->ReadTable();

            // Read and initialize data from the NetCDF land file.
            noahmpio->ReadLandMain();

            // Compute additional initial values that were not supplied
            // by the NetCDF land file.
            noahmpio->InitMain();

            // Write initial plotfile for land with the tag 0
            amrex::Print() << "Noah-MP writing lnd.nc file at lev: " << lev << std::endl;
            noahmpio->WriteLand(0);
        }
  }

  amrex::Print() << "Noah-MP initialization completed" << std::endl;

};

void
NOAH::Advance_With_State (const int& lev,
                          MultiFab& cons_in,
                          MultiFab& xvel_in,
                          MultiFab& yvel_in,
                          MultiFab* hfx3_out,
                          MultiFab* qfx3_out,
                          const amrex::Real& dt,
                          const int& nstep) {

    Box domain = m_geom.Domain();

    amrex::Print () << "Noah-MP driver started at time step: " << nstep+1 << std::endl;

    // Loop over blocks to copy forcing data to Noahmp, drive the land model,
    // and copy data back to ERF Multifabs.
    int idb = 0;
    for (amrex::MFIter mfi(xvel_in, false); mfi.isValid(); ++mfi, ++idb) {

        const amrex::Box& bx = mfi.tilebox();

        // Check if tile is at the lower boundary in lower z direction
        if (bx.smallEnd(2) == domain.smallEnd(2)) {

            NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

            const amrex::Array4<const amrex::Real>& U_PHY = xvel_in.const_array(mfi);
            const amrex::Array4<const amrex::Real>& V_PHY = yvel_in.const_array(mfi);
            const amrex::Array4<const amrex::Real>& QV_TH = cons_in.const_array(mfi);

            amrex::Array4<amrex::Real> SHBXY = hfx3_out->array(mfi);
            amrex::Array4<amrex::Real> EVBXY = qfx3_out->array(mfi);

            // Copy forcing data from ERF to Noahmp.
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
            {
                noahmpio->U_PHY(i,1,j) = U_PHY(i,j,0);
                noahmpio->V_PHY(i,1,j) = V_PHY(i,j,0);
                noahmpio->T_PHY(i,1,j) = QV_TH(i,j,0,RhoTheta_comp)/QV_TH(i,j,0,Rho_comp);
                noahmpio->QV_CURR(i,1,j) = QV_TH(i,j,0,RhoQ1_comp)/QV_TH(i,j,0,Rho_comp);

            });

            // Call the noahmpio driver code. This runs the land model forcing for
            // each object in noahmpio_vect that represent a block in the domain.
            noahmpio->itimestep = nstep+1;
            noahmpio->DriverMain();

            // Copy forcing data from Noahmp to ERF
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
            {
                SHBXY(i,j,0) = noahmpio->SHBXY(i,j);
                EVBXY(i,j,0) = noahmpio->EVBXY(i,j);
            });

        }
    }
    amrex::Print () << "Noah-MP driver completed" << std::endl;
};
