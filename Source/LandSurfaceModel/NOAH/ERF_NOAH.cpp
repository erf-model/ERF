
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
    LsmVarMap = {LsmVar_NOAH::theta};

    LsmVarName.resize(m_lsm_size);
    LsmVarName = {"theta"};

    // NOTE: lsm data is not used for Noahmp, however, the initialization is done
    //       to maintin consistency with IO and Driver interfaces that depend on
    //       this data. We eventually want to tweak those interfaces so we don't
    //       have to allocate lsm_data while using Noahmp lsm.

    // NOTE: All boxes in ba extend from zlo to zhi, so this transform is valid.
    //       If that were to change, the dm and new ba are no longer valid and
    //       direct copying between lsm data/flux vars cannot be done in a parfor.

    // Set box array for lsm data
    IntVect ng(0,0,1);
    BoxArray ba = cons_in.boxArray();
    DistributionMapping dm = cons_in.DistributionMap();
    BoxList bl_lsm = ba.boxList();
    for (auto& b : bl_lsm) {
        b.setBig(2,khi_lsm);                  // First point below the surface
        b.setSmall(2,khi_lsm - m_nz_lsm + 1); // Last point below the surface
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
        Real theta_0 = m_theta_dir;
        lsm_fab_vars[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_vars[ivar]->setVal(theta_0);

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

    int idb = 0;
    for (amrex::MFIter mfi(xvel_in, false); mfi.isValid(); ++mfi, ++idb)
    {
        const amrex::Box& bx = mfi.tilebox();

        // Fix vertical level: work on 2D slice at bottom vertical level (usually 0)
        int klev = domain.smallEnd(2);
        // Construct a 2D box with vertical thickness 1 at klev = 0
        amrex::Box bx_2d({bx.smallEnd(0), bx.smallEnd(1), 0},
                         {bx.bigEnd(0), bx.bigEnd(1), 0});

        if (bx.smallEnd(2) == domain.smallEnd(2)) {

            NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

            // Choose pinned arena on GPU builds, default arena on CPU
            amrex::Arena* pinned_arena = nullptr;
#ifdef AMREX_USE_GPU
            pinned_arena = amrex::The_Pinned_Arena();
#else
            pinned_arena = amrex::The_Arena();
#endif

            // Allocate pinned host FABs for velocities (1 comp each)
            FArrayBox U_host(bx_2d, 1, pinned_arena);
            FArrayBox V_host(bx_2d, 1, pinned_arena);

            // Number of cons_in components used in calculation (adjust to your actual comps)
            constexpr int ncons = 3; // e.g., RhoTheta, RhoQ1, Rho
            FArrayBox cons_host(bx_2d, ncons, pinned_arena);

            // Copy velocity slices (klev fixed) from device to pinned host
#ifdef AMREX_USE_GPU
            {
                // Copy xvel_in comp 0, slice at klev
                auto const& xvel_arr = xvel_in.const_array(mfi);
                auto const& yvel_arr = yvel_in.const_array(mfi);
                auto& U_arr = U_host.array();
                auto& V_arr = V_host.array();

                // Loop on CPU to copy slice at klev, since no direct 2D copy in AMReX API
                // This is host code; amrex::Gpu::copy cannot copy a slice directly, so do manual copy:
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
                {
                   U_arr(i,j,0) = xvel_arr(i,j,klev);
                   V_arr(i,j,0) = yvel_arr(i,j,klev);
                }
            }

            // Copy cons_in components slice at klev
            {
                auto const& cons_arr = cons_in.const_array(mfi);
                auto& c_arr = cons_host.array();

                // For each component copy slice at klev
                for (int comp = 0; comp < ncons; ++comp) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
                    {
                       c_arr(i,j,0,comp) = cons_arr(i,j,klev,comp);
                    }
                }
            }
#else
            // CPU only: just copy boxes directly (all comps, 1 k level)
            U_host.copy(xvel_in[mfi], 0, 0, 1);
            V_host.copy(yvel_in[mfi], 0, 0, 1);
            cons_host.copy(cons_in[mfi], 0, 0, ncons);
#endif

            // Create array accessors for CPU computations
            const auto& U_arr = U_host.array();
            const auto& V_arr = V_host.array();
            const auto& cons_arr = cons_host.array();

            // Aliases for noahmpio CPU arrays (assumed allocated)
            auto& T_cpu  = noahmpio->T_PHY;
            auto& QV_cpu = noahmpio->QV_CURR;
            auto& U_cpu  = noahmpio->U_PHY;
            auto& V_cpu  = noahmpio->V_PHY;

            // Loop over 2D slice and compute temperature and mixing ratio
            for (int j = bx_2d.smallEnd(1); j <= bx_2d.bigEnd(1); ++j) {
                for (int i = bx_2d.smallEnd(0); i <= bx_2d.bigEnd(0); ++i) {
                    U_cpu(i,1,j) = U_arr(i,j,0);
                    V_cpu(i,1,j) = V_arr(i,j,0);

                    amrex::Real rho_theta = cons_arr(i,j,0,RhoTheta_comp);
                    amrex::Real rho_q1 = cons_arr(i,j,0,RhoQ1_comp);
                    amrex::Real rho = cons_arr(i,j,0,Rho_comp);

                    T_cpu(i,1,j) = rho_theta / rho;
                    QV_cpu(i,1,j) = rho_q1 / rho;
                }
            }

            // Run land surface model timestep
            noahmpio->itimestep = nstep + 1;
            noahmpio->DriverMain();

            // Copy land surface outputs back to MultiFabs
#ifdef AMREX_USE_GPU
            // Create pinned FABs for SHBXY and EVBXY output on CPU host
            FArrayBox SHBXY_host(bx_2d, 1, pinned_arena);
            FArrayBox EVBXY_host(bx_2d, 1, pinned_arena);
            auto& SHB_arr = SHBXY_host.array();
            auto& EVB_arr = EVBXY_host.array();
            amrex::Array4<amrex::Real> SHBXY = hfx3_out->array(mfi);
            amrex::Array4<amrex::Real> EVBXY = qfx3_out->array(mfi);

            // Copy from noahmpio CPU arrays into pinned host FABs
            for (int j = bx_2d.smallEnd(1); j <= bx_2d.bigEnd(1); ++j) {
                for (int i = bx_2d.smallEnd(0); i <= bx_2d.bigEnd(0); ++i) {
                    SHB_arr(i,j,0) = noahmpio->SHBXY(i,j);
                    EVB_arr(i,j,0) = noahmpio->EVBXY(i,j);
                }
            }

            // Copy from pinned host FABs back to device MultiFabs
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
            {
                SHBXY(i,j,0) = SHB_arr(i,j);
                EVBXY(i,j,0) = EVB_arr(i,j); 
            } 

#else
            // CPU only: copy arrays directly
            amrex::Array4<amrex::Real> SHBXY = hfx3_out->array(mfi);
            amrex::Array4<amrex::Real> EVBXY = qfx3_out->array(mfi);

            for (int j = bx_2d.smallEnd(1); j <= bx_2d.bigEnd(1); ++j) {
                for (int i = bx_2d.smallEnd(0); i <= bx_2d.bigEnd(0); ++i) {
                    SHBXY(i,j,0) = noahmpio->SHBXY(i,j);
                    EVBXY(i,j,0) = noahmpio->EVBXY(i,j);
                }
            }
#endif

        } // if vertical slice

    } // MFIter

    amrex::Print() << "Noah-MP driver completed" << std::endl;
};
