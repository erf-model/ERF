/*
 * Lifecycle: NOAHMP::Init. Builds the surface (lsm) geometry and coupling
 * MultiFabs, sizes + initializes one NoahmpIO_type per box, and broadcasts the
 * firing parameters. See dev/spec-noahmp-api.md §4.
 */

#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include <ERF_NOAHMP.H>
#include <ERF_Constants.H>
#include <NoahmpFatal.H>

using namespace amrex;

/* Initialize lsm data structures */
void
NOAHMP::Init (const int& lev,
              const MultiFab& cons_in,
              const Geometry& geom,
              const Geometry& geom0,
              Vector<BCRec>& domain_bcs_type,
              IntVect& refRatio,
              const Real& dt,
              Vector<Vector<std::string>>& nc_init_file)
{
    // Install Noah-MP's fatal-error handler once (thread-safe, runs on the first
    // Init across all levels). Noah-MP carries no MPI/AMReX dependency of its own
    // and instead calls NoahmpIO_fatal(); routing that through amrex::Abort makes a
    // fatal error on any rank -- from either the C++ or the Fortran coupling side
    // -- propagate to all ranks via MPI_Abort, rather than deadlocking peers in
    // the next collective. See NoahmpFatal.H.
    static const bool noahmp_fatal_installed = []() {
        NoahmpIO_set_fatal_handler([](const char* msg){
            amrex::Abort(msg ? msg : "Noah-MP fatal error");
        });
        return true;
    }();
    amrex::ignore_unused(noahmp_fatal_installed);

    m_lev   = lev;
    m_dt    = dt;
    m_geom  = geom;
    m_geom0 = geom0;
    m_domain_bcs_type = domain_bcs_type;
    m_refRatio = refRatio;

    Box domain = geom.Domain();
    khi_lsm    = domain.smallEnd(2) - 1;

    // Number of Noah-MP soil layers. Sizing must be known on EVERY rank and level
    // before the (collective) LSM data fabs are built, but the authoritative NSOIL
    // lives in the Fortran namelist.erf, which is only read later, per land box.
    // Resolve it from ParmParse (erf.lsm_nsoil, default 4 = standard config); this
    // is the same value the parent already used to size lsm_data in MakeNewLevel
    // (it calls Lsm_Data_Size() before Init). The namelist NSOIL is asserted to
    // agree below, so a mismatch fails loudly. m_lsm_data_size is set here too.
    m_ensure_nsoil_resolved();

    // The fixed 2D fields are identity-mapped; their names mirror the enum order.
    LsmDataMap.resize(m_lsm_data_size);
    LsmDataName.resize(m_lsm_data_size);
    for (int i(0); i < LsmData_NOAHMP::NumVars; ++i) { LsmDataMap[i] = i; }
    {
        // Generated from the same X-macro registry as the LsmData_NOAHMP enum
        // (ERF_NOAHMP_Fields.H), so the names and the enum order can never drift
        // -- the historical hand-maintained list is gone. See spec-noahmp-reorg.md.
        const std::vector<std::string> fixed_names = {
#define NOAHMP_X(name) #name,
            NOAHMP_LSMDATA_FIELDS(NOAHMP_X)
#undef NOAHMP_X
        };
        AMREX_ALWAYS_ASSERT(int(fixed_names.size()) == LsmData_NOAHMP::NumVars);
        for (int i(0); i < LsmData_NOAHMP::NumVars; ++i) { LsmDataName[i] = fixed_names[i]; }
    }
    // Append the per-layer soil profile: 3 contiguous groups of m_nsoil, layer index
    // 1-based in the field name to match the Noah-MP / WRF SMOIS_k convention.
    {
        const char* group[3] = {"smois", "sh2o", "tslb"};
        for (int g(0); g < 3; ++g) {
            for (int k(0); k < m_nsoil; ++k) {
                int idx = soil_data_idx(g,k);
                LsmDataMap[idx]  = idx;
                LsmDataName[idx] = std::string(group[g]) + "_" + std::to_string(k+1);
            }
        }
    }


    LsmFluxMap.resize(m_lsm_flux_size);
    LsmFluxMap = {LsmFlux_NOAHMP::t_flux         , LsmFlux_NOAHMP::q_flux         ,
                  LsmFlux_NOAHMP::tau13          , LsmFlux_NOAHMP::tau23          };
    LsmFluxName.resize(m_lsm_flux_size);
    LsmFluxName = {"t_flux"         , "q_flux"         ,
                   "tau13"          , "tau23"          };

    ParmParse pp("erf");
    pp.query("plot_int_1" , m_plot_int_1);

    // NOTE: All boxes in ba extend from zlo to zhi, so this transform is valid.
    //       If that were to change, the dm and new ba are no longer valid and
    //       direct copying between lsm data/flux vars cannot be done in a parfor.

    // Set 2D box array for lsm data
    IntVect ng(1,1,0);
    BoxArray ba = cons_in.boxArray();
    DistributionMapping dm = cons_in.DistributionMap();
    BoxList bl_lsm = ba.boxList();
    for (auto& b : bl_lsm) { b.setRange(2,0); }
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

    // Create the data. lsm_fab_data / lsm_lev0_data are runtime-sized (fixed 2D
    // fields + 3*m_nsoil soil-layer fields) rather than a compile-time Array, so
    // the soil profile scales with NSOIL. lsm_lev0_data pointers are populated
    // later by the parent via Lsm_Set_Lev0_Data_Ptr.
    lsm_fab_data.resize(m_lsm_data_size);
    lsm_lev0_data.resize(m_lsm_data_size, nullptr);
    for (auto ivar = 0; ivar < m_lsm_data_size; ++ivar) {
        // State vars are CC
        lsm_fab_data[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_data[ivar]->setVal(lsm_undefined);
    }

    // Create the fluxes
    for (auto ivar = 0; ivar < LsmFlux_NOAHMP::NumVars; ++ivar) {
        // NOTE: Fluxes are CC with ghost cells for averaging
        lsm_fab_flux[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, IntVect(1,1,0));
        lsm_fab_flux[ivar]->setVal(lsm_undefined);
    }

    m_has_nc_file = (!nc_init_file[lev].empty());
    if (m_has_nc_file) {
        Print() << "Noah-MP initialization started" << std::endl;

        // Size noahmpio_vect to the local blocks (boxes). A rank owning no boxes
        // leaves it empty and relies on the class-level m_itimestep/m_dtbl instead.
        if (cons_in.local_size() > 0) {
            noahmpio_vect.resize(cons_in.local_size(), lev);
        }

        // Allocate pinned buffer space for all the boxes
        noahmp_input_tmp.resize(cons_in.local_size());
        noahmp_output_tmp.resize(cons_in.local_size());

        int klo = domain.smallEnd(2);

        // Iterate over multifab and noahmpio object together. Multifabs is
        // used to extract size of blocks and set bounds for noahmpio objects.
        int idb = 0;
        for (MFIter mfi(cons_in); mfi.isValid(); ++mfi, ++idb) {

            // Get bounds for the tile
            Box bx = mfi.tilebox();

            // Check if tile is at the lower boundary in lower z direction
            if (bx.smallEnd(2) != klo) { continue; }

            // Make a slab
            bx.makeSlab(2,klo);

            // Allocate pinned buffers for each box
            // Output buffer carries the fixed 2D outputs plus 3 soil groups of m_nsoil.
            noahmp_input_tmp[idb]  = std::make_unique<FArrayBox>(bx, NoahmpInputComp::NumComps , The_Pinned_Arena());
            noahmp_output_tmp[idb] = std::make_unique<FArrayBox>(bx, NoahmpOutputComp::NumComps + 3*m_nsoil, The_Pinned_Arena());

            // Get reference to the noahmpio object
            NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

            // Pass idb context to noahmpio
            noahmpio->blkid = idb;

            // Pass level context to noahmpio
            noahmpio->level = lev;

            // Initialize scalar values
            noahmpio->ScalarInitDefault();

            // Store the rank of process for noahmp
            noahmpio->rank = ParallelDescriptor::MyProc();

            // Store parallel communicator for noahmp
            noahmpio->comm = MPI_Comm_c2f(ParallelDescriptor::Communicator());

            // Read namelist.erf file. This file contains
            // noahmpio specific parameters and is read by
            // the Fortran side of the implementation.
            noahmpio->ReadNamelist();

            // The LSM data fabs were sized above from erf.lsm_nsoil (default 4).
            // Assert the authoritative namelist NSOIL agrees, so a non-standard soil
            // discretization that forgot to set erf.lsm_nsoil fails loudly here
            // rather than silently truncating / overrunning the soil diagnostics.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(noahmpio->nsoil == m_nsoil,
                "namelist.erf NSOIL does not match erf.lsm_nsoil (default 4); "
                "set erf.lsm_nsoil to the Noah-MP soil-layer count");

            // Read the headers from the NetCDF land file. This is also
            // implemented on the Fortran side of things currently.
            noahmpio->ReadLandHeader();

            // Extract tile bounds and set them to their corresponding
            // noahmpio variables. At present we will set all the variables
            // corresponding to domain, memory, and tile to the same bounds.
            // This will be changed later if we want to do special memory
            // management for expensive use cases.
            noahmpio->xstart = bx.smallEnd(0);
            noahmpio->xend   = bx.bigEnd(0);
            noahmpio->ystart = bx.smallEnd(1);
            noahmpio->yend   = bx.bigEnd(1);

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
            Print() << "Noah-MP writing lnd.nc file at lev: " << lev << std::endl;
            noahmpio->WriteLand(0);
        }

        // Broadcast DTBL and the initial substep counter to every rank so the firing
        // decision in Advance_With_State is identical everywhere. Land-free ranks use
        // sentinels that lose the max-reduction to any real value.
        m_dtbl      = noahmpio_vect.empty() ? std::numeric_limits<Real>::lowest()
                                            : static_cast<Real>(noahmpio_vect[0].DTBL);
        m_itimestep = noahmpio_vect.empty() ? std::numeric_limits<int>::lowest()
                                            : noahmpio_vect[0].itimestep;
        ParallelDescriptor::ReduceRealMax(m_dtbl);
        ParallelDescriptor::ReduceIntMax(m_itimestep);

        // Guard against a degenerate decomposition in which no rank owns a land box.
        AMREX_ALWAYS_ASSERT(m_dtbl > Real(0.0));
        AMREX_ALWAYS_ASSERT(m_dt <= m_dtbl);

        Print() << "Noah-MP initialization completed" << std::endl;
    } // has nc_init_file

};
