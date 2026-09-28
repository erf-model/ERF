/*
 * NOAHMP::Init: builds the surface (lsm) geometry and coupling MultiFabs, sizes
 * and initializes one NoahmpIO_type per box, and broadcasts the firing parameters.
 */

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <limits>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <ERF_NOAHMP.H>
#include <ERF_Constants.H>
#include <NoahmpFatal.H>

using namespace amrex;

namespace {

//
// Check, before the Fortran driver runs, every file it is about to open.
//
// The driver reports a missing file with a Fortran WRITE on rank 0 and then calls
// NoahmpIO_abort(), which reaches amrex::Abort through the handler installed in
// NOAHMP::Init with no message at all. The user sees only "Noah-MP fatal error" -- for
// a missing namelist.erf, a missing NoahmpTable.TBL and a missing land setup file alike
// -- with the one line that says which file it wanted lost with the Fortran unit. So
// look for them here and name the missing one.
//
// The driver reads its setup file from ERF_SETUP_FILE_01/02/03 in namelist.erf for
// levels 0/1/2 (NoahmpReadNamelistMod.F90) and has no fourth.
//
void
noahmp_preflight (int lev)
{
    const std::string namelist = "namelist.erf";
    const std::string table    = "NoahmpTable.TBL";

    if (!amrex::FileExists(namelist)) {
        amrex::Abort("Noah-MP: " + namelist + " was not found in the run directory. The "
                     "Noah-MP driver reads its physics options and the name of its land "
                     "setup file from it; a template is in "
                     "Submodules/Noah-MP/drivers/erf/tests/namelist.erf.");
    }
    if (!amrex::FileExists(table)) {
        amrex::Abort("Noah-MP: " + table + " was not found in the run directory. Copy "
                     "it from Submodules/Noah-MP/parameters/NoahmpTable.TBL.");
    }
    if (lev > 2) {
        amrex::Abort("Noah-MP: the driver reads a land setup file for levels 0-2 only "
                     "(ERF_SETUP_FILE_01..03 in namelist.erf); level " + std::to_string(lev) +
                     " has none. Limit amr.max_level to 2.");
    }

    // Fortran namelists are case-insensitive and take either quote. The key has to open
    // its line, so a commented-out assignment ("! ERF_SETUP_FILE_01 = ...") is skipped;
    // a later assignment overrides an earlier one, as it does in Fortran.
    const std::string key = "ERF_SETUP_FILE_0" + std::to_string(lev + 1);
    const std::regex assign("^\\s*" + key + "\\s*=\\s*['\"]([^'\"]*)['\"]",
                            std::regex::icase);
    std::string setup_file;
    bool found_key = false;
    {
        std::ifstream in(namelist);
        std::string line;
        std::smatch m;
        while (std::getline(in, line)) {
            if (std::regex_search(line, m, assign)) {
                found_key  = true;
                setup_file = m[1].str();
            }
        }
    }
    // Trim, since Fortran pads its character variables with blanks.
    const auto first = setup_file.find_first_not_of(' ');
    const auto last  = setup_file.find_last_not_of(' ');
    setup_file = (first == std::string::npos) ? std::string{}
                                              : setup_file.substr(first, last - first + 1);

    if (!found_key || setup_file.empty()) {
        amrex::Abort("Noah-MP: " + namelist + " does not set " + key + ", the wrfinput-format "
                     "NetCDF file Noah-MP reads its land state from at level " +
                     std::to_string(lev) + " (soil and vegetation type, soil temperature and "
                     "moisture, and the WRF grid attributes).");
    }
    if (!amrex::FileExists(setup_file)) {
        amrex::Abort("Noah-MP: the land setup file '" + setup_file + "' named by " + key +
                     " in " + namelist + " does not exist.");
    }
}

} // namespace

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
    // Install Noah-MP's fatal-error handler once: route NoahmpIO_fatal() through
    // amrex::Abort so a fatal error propagates via MPI_Abort. See NoahmpFatal.H.
    static const bool noahmp_fatal_installed = []() {
        // The driver's own aborts pass no message: it WRITEs the reason to standard output
        // on rank 0 first, and that line is usually lost when the job is killed. NOAHMP::Init
        // checks the files the driver opens beforehand (noahmp_preflight); what reaches
        // this default is a problem inside one of them, so point at the likely ones.
        NoahmpIO_set_fatal_handler([](const char* msg){
            amrex::Abort(msg ? msg :
                "Noah-MP fatal error in the Fortran driver. Its explanation is written to "
                "standard output on rank 0 and is usually lost when the job stops; check for a "
                "malformed namelist.erf (every *_TIMESTEP and *_OPTION is an integer), a "
                "NoahmpTable.TBL that does not match the land-use scheme (MMINLU), or a land "
                "setup file missing one of the fields Noah-MP reads.");
        });
        return true;
    }();
    amrex::ignore_unused(noahmp_fatal_installed);

    // Noah-MP's own physics checks end the run with a bare Fortran STOP -- there are 27
    // across src/ and drivers/erf/, e.g. "Error: Solar radiation budget problem in NoahMP
    // LSM" -- and a STOP exits with status 0. The run then dies part-way, with no final
    // plotfile or checkpoint. Under Open MPI's launcher that is still reported as a
    // failure, because the rank exits without MPI_Finalize; but run directly, without a
    // launcher, or under one that does not enforce that rule, it reports success.
    //
    // STOP leaves through exit(), which runs atexit handlers. A normal ERF run has called
    // amrex::Finalize before it gets there, and amrex::Abort never calls exit() at all, so
    // a process leaving through exit() while AMReX is still initialised has stopped where
    // it should not have. Say so and leave with a failure status. The Fortran is not
    // modified; this only changes what the process reports.
    static const bool noahmp_stop_trap_installed = []() {
        std::atexit([]() {
            if (amrex::Initialized()) {
                std::fputs("ERF: the run ended through exit() while AMReX was still running "
                           "-- typically a Fortran STOP inside Noah-MP, whose message is "
                           "printed above. Reporting failure instead of exit status 0.\n",
                           stderr);
                std::fflush(stderr);
                std::_Exit(EXIT_FAILURE);
            }
        });
        return true;
    }();
    amrex::ignore_unused(noahmp_stop_trap_installed);

    // dt is a placeholder: ERF::make_lsm_at_level passes zero, because the level's step is
    // not known yet when the level is built. Nothing here may depend on it; the step is
    // checked against NOAH_TIMESTEP in Advance_With_State, where it is real.
    amrex::ignore_unused(dt);

    m_lev   = lev;
    m_geom  = geom;
    m_geom0 = geom0;
    m_domain_bcs_type = domain_bcs_type;
    m_refRatio = refRatio;

    Box domain = geom.Domain();

    // Resolve NSOIL from erf.lsm_nsoil before building the collective LSM fabs (same
    // value the parent used via Lsm_Data_Size()); namelist NSOIL asserted below.
    m_ensure_nsoil_resolved();

    // The fixed 2D fields are identity-mapped; their names mirror the enum order.
    LsmDataMap.resize(m_lsm_data_size);
    LsmDataName.resize(m_lsm_data_size);
    for (int i(0); i < LsmData_NOAHMP::NumVars; ++i) { LsmDataMap[i] = i; }
    {
        // Names from the same registry as the enum, so they cannot drift.
        const std::vector<std::string> fixed_names = {
            NOAHMP_LSMDATA_FIELDS(NOAHMP_QUOTE)
        };
        AMREX_ALWAYS_ASSERT(int(fixed_names.size()) == LsmData_NOAHMP::NumVars);
        for (int i(0); i < LsmData_NOAHMP::NumVars; ++i) { LsmDataName[i] = fixed_names[i]; }
    }
    // Per-layer soil profile: 3 groups of m_nsoil, layer index 1-based (WRF SMOIS_k).
    {
        const char* group[m_num_soil_groups] = {"smois", "sh2o", "tslb"};
        for (int g(0); g < m_num_soil_groups; ++g) {
            for (int k(0); k < m_nsoil; ++k) {
                int idx = soil_data_idx(g,k);
                LsmDataMap[idx]  = idx;
                LsmDataName[idx] = std::string(group[g]) + "_" + std::to_string(k+1);
            }
        }
    }

    LsmFluxMap  = {LsmFlux_NOAHMP::t_flux         , LsmFlux_NOAHMP::q_flux         ,
                  LsmFlux_NOAHMP::tau13          , LsmFlux_NOAHMP::tau23          };
    LsmFluxName = {"t_flux"         , "q_flux"         ,
                   "tau13"          , "tau23"          };

    // NOTE: relies on all boxes in ba spanning zlo..zhi; otherwise dm/ba no longer
    //       line up and lsm data/flux vars can't be copied directly in a parfor.

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

    // Create the data (CC), runtime-sized (fixed 2D fields + 3*m_nsoil) so soil scales
    // with NSOIL. lsm_lev0_data pointers are populated later by the parent.
    lsm_fab_data.resize(m_lsm_data_size);
    lsm_lev0_data.resize(m_lsm_data_size, nullptr);
    for (auto ivar = 0; ivar < m_lsm_data_size; ++ivar) {
        lsm_fab_data[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_data[ivar]->setVal(lsm_undefined);
    }

    // Create the fluxes (CC with ghost cells for averaging)
    for (auto ivar = 0; ivar < LsmFlux_NOAHMP::NumVars; ++ivar) {
        lsm_fab_flux[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_flux[ivar]->setVal(lsm_undefined);
    }

    // Level 0 always runs the Noah-MP driver: the driver reads its land state from the
    // file ERF_SETUP_FILE_01 in namelist.erf names, whatever ERF's own initialization is,
    // and there is no coarser level for the other branch (interp_from_lev0) to take it
    // from. This used to follow only from ERF::nc_init_file's static default being
    // {{""}} -- one empty string at level 0, so the vector below is never empty there,
    // with or without erf.nc_init_file_0 -- which is how an idealized run could use
    // Noah-MP at all. Say it outright so it does not rest on that placeholder.
    // A finer level runs the driver only if it has an init file of its own.
    m_has_nc_file = (lev == 0) || (!nc_init_file[lev].empty());
    if (m_has_nc_file) {
        Print() << "Noah-MP initialization started" << std::endl;

        noahmp_preflight(lev);

        // Size noahmpio_vect to the local boxes. A rank owning no boxes leaves it
        // empty and relies on the class-level m_itimestep/m_dtbl instead.
        if (cons_in.local_size() > 0) {
            noahmpio_vect.resize(cons_in.local_size(), lev);
        }

        // Pinned buffer space for all the boxes
        noahmp_input_tmp.resize(cons_in.local_size());
        noahmp_output_tmp.resize(cons_in.local_size());

        int klo = domain.smallEnd(2);

        // Iterate over the multifab and noahmpio objects together, using the
        // multifab to set the per-box bounds on each noahmpio object.
        int idb = 0;
        for (MFIter mfi(cons_in); mfi.isValid(); ++mfi, ++idb) {

            Box bx = mfi.tilebox();

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bx.smallEnd(2) == klo,
                "NoahMP Init: box does not start at klo; z-decomposed grids are unsupported.");

            bx.makeSlab(2,klo);

            // Pinned buffers per box; output carries the 2D outputs + 3 soil groups.
            noahmp_input_tmp[idb]  = std::make_unique<FArrayBox>(bx, NoahmpInputComp::NumComps , The_Pinned_Arena());
            noahmp_output_tmp[idb] = std::make_unique<FArrayBox>(bx, NoahmpOutputComp::NumComps + m_num_soil_groups*m_nsoil, The_Pinned_Arena());

            NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

            noahmpio->blkid = idb;
            noahmpio->level = lev;
            noahmpio->ScalarInitDefault();
            noahmpio->rank = ParallelDescriptor::MyProc();
            noahmpio->comm = MPI_Comm_c2f(ParallelDescriptor::Communicator());

            // namelist.erf holds noahmpio-specific parameters, read Fortran-side.
            noahmpio->ReadNamelist();

            // Assert namelist NSOIL matches erf.lsm_nsoil (used to size the fabs) so a
            // mismatch fails loudly rather than truncating the soil diagnostics.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(noahmpio->nsoil == m_nsoil,
                "namelist.erf NSOIL does not match erf.lsm_nsoil (default 4); "
                "set erf.lsm_nsoil to the Noah-MP soil-layer count");

            // NetCDF land-file headers (also Fortran-side)
            noahmpio->ReadLandHeader();

            // Set domain/memory/tile bounds from the tile. All three are set to the
            // same bounds for now; may change for special memory management later.
            noahmpio->xstart = bx.smallEnd(0);
            noahmpio->xend   = bx.bigEnd(0);
            noahmpio->ystart = bx.smallEnd(1);
            noahmpio->yend   = bx.bigEnd(1);

            // Domain, tile, and memory bounds are all equal for now (single-slab model).
            auto set_grid_bounds = [](NoahmpIO_type* io, int x0, int x1, int y0, int y1) {
                io->ids=io->its=io->ims=x0; io->ide=io->ite=io->ime=x1;
                io->jds=io->jts=io->jms=y0; io->jde=io->jte=io->jme=y1;
                io->kds=io->kts=io->kms=1;  io->kde=io->kte=io->kme=2;
            };
            set_grid_bounds(noahmpio, noahmpio->xstart, noahmpio->xend, noahmpio->ystart, noahmpio->yend);

            // Allocate Fortran IO memory from the bounds above + namelist/header info
            noahmpio->VarInitDefault();

            // NoahmpTable.TBL input
            noahmpio->ReadTable();

            // Read/initialize from the NetCDF land file
            noahmpio->ReadLandMain();

            // Compute initial values not supplied by the land file
            noahmpio->InitMain();
        }

        // Initial land plotfile (tag 0); must run on the land-owning subcomm, not the
        // full comm, or a land-free rank never joins the collective nf90_create.
        Print() << "Noah-MP writing lnd.nc file at lev: " << lev << std::endl;
        with_land_comm([](NoahmpIO_type& noahmpio) { noahmpio.WriteLand(0); });

        // Broadcast DTBL and the initial substep counter so the firing decision is
        // identical on every rank. Land-free ranks use max-reduction-losing sentinels.
        m_dtbl      = noahmpio_vect.empty() ? std::numeric_limits<Real>::lowest()
                                            : static_cast<Real>(noahmpio_vect[0].DTBL);
        m_itimestep = noahmpio_vect.empty() ? std::numeric_limits<int>::lowest()
                                            : noahmpio_vect[0].itimestep;
        ParallelDescriptor::ReduceRealMax(m_dtbl);
        ParallelDescriptor::ReduceIntMax(m_itimestep);

        // Guard against a decomposition in which no rank owns a land box.
        // DTBL is real(NOAH_TIMESTEP) from namelist.erf, and NOAH_TIMESTEP defaults to
        // -9999 in the driver, so this is what an absent NOAH_TIMESTEP looks like. (It is
        // not about land: noahmpio_vect is sized by boxes, land or not, so every rank that
        // owns a box contributes a DTBL to the max.)
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_dtbl > Real(0.0),
            "Noah-MP: the Noah-MP timestep is not positive. Set NOAH_TIMESTEP (in seconds, "
            "an integer) in namelist.erf.");
        // The ERF-step-versus-NOAH_TIMESTEP constraint is enforced in Advance_With_State:
        // the dt this function receives is a placeholder zero (see the top of Init), so a
        // check against it here -- there used to be one -- could never fire.

        Print() << "Noah-MP initialization completed" << std::endl;
    } // has nc_init_file

};
