/*
 * NOAHMP::Init: builds the surface (lsm) geometry and coupling MultiFabs, sizes
 * and initializes one NoahmpIO_type per box, and broadcasts the firing parameters.
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

    m_has_nc_file = (!nc_init_file[lev].empty());
    if (m_has_nc_file) {
        Print() << "Noah-MP initialization started" << std::endl;

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
        AMREX_ALWAYS_ASSERT(m_dtbl > Real(0.0));
        AMREX_ALWAYS_ASSERT(m_dt <= m_dtbl);

        Print() << "Noah-MP initialization completed" << std::endl;
    } // has nc_init_file

};
