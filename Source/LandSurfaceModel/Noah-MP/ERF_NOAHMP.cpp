
#include<iostream>
#include<string>
#include<limits>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include <ERF_NOAHMP.H>
#include <ERF_Constants.H>
#include <ERF_EOS.H>
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
        const std::vector<std::string> fixed_names = {
                   "t_sfc"             , "sfc_emis"          ,
                   "sfc_alb_dir_vis"   , "sfc_alb_dir_nir"   ,
                   "sfc_alb_dif_vis"   , "sfc_alb_dif_nir"   ,
                   "cos_zenith_angle"  , "sw_flux_dn"        ,
                   "sw_flux_dn_dir_vis", "sw_flux_dn_dir_nir",
                   "sw_flux_dn_dif_vis", "sw_flux_dn_dif_nir",
                   "lw_flux_dn"        ,
                   "grdflx"            , "fira"              ,
                   "sav"               , "sag"               ,
                   "albedo"            , "sfcrunoff"         ,
                   "udrunoff"          , "smstav"            ,
                   "smstot"            };
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

void
NOAHMP::Plot_Landfile(const int& nstep)
{
    for (NoahmpIO_type &noahmpio : noahmpio_vect) {
        noahmpio.WriteLand(nstep);
    }
}

void
NOAHMP::Advance_With_State (const int& lev,
                            MultiFab& cons_in,
                            MultiFab& xvel_in,
                            MultiFab& yvel_in,
                            MultiFab* /*hfx3_out*/,
                            MultiFab* /*qfx3_out*/,
                            const SurfacePrecipAccumulationSources& precip_sources,
                            const Real& elapsed_time,
                            const Real& dt,
                            const int& nstep,
                            const bool updated_lev0)
{
    if (!m_has_nc_file) {
        // NOTE: Do not try to interpolate if lev 0 was just updated. Since Noah is
        //       called post-step the fluxes & data will contain lsm_undefined values.
        if (!updated_lev0) {
            Print () << "Noah-MP interpolation at level " << lev << " started at time step: " << nstep+1 << std::endl;
            m_updated = true;
            for (int ivar(0); ivar<m_lsm_data_size; ++ivar) {
                // Interpolate from lev 0 to obtain the lsm data
                InterpFromCoarseLevel(*lsm_fab_data[ivar], lsm_fab_data[ivar]->nGrowVect(),
                                      IntVect(0,0,0), // do NOT fill ghost cells outside the domain
                                      *lsm_lev0_data[ivar], 0, 0, 1,
                                      m_geom0, m_geom,
                                      m_refRatio, &cell_cons_interp,
                                      m_domain_bcs_type, BCVars::cons_bc);
            } // ivar

            // NOTE: The surface layer class wrote into the noah flux data structures
            //       where lsm_undefined values existed. This makes the noah flux
            //       complete on the coarse grid and we can safely interpolate.
            for (int ivar(0); ivar<LsmFlux_NOAHMP::NumVars; ++ivar) {
                // Interpolate from lev 0 to obtain the lsm fluxes
                InterpFromCoarseLevel(*lsm_fab_flux[ivar], lsm_fab_flux[ivar]->nGrowVect(),
                                      IntVect(0,0,0), // do NOT fill ghost cells outside the domain
                                      *lsm_lev0_flux[ivar], 0, 0, 1,
                                      m_geom0, m_geom,
                                      m_refRatio, &cell_cons_interp,
                                      m_domain_bcs_type, BCVars::cons_bc);
            } // ivar
            Print () << "Noah-MP interpolation at level " << lev << " completed" << std::endl;
        } else {
            m_updated = false;
        }
    } else {
        // Verify we need to take another LSM step. Use the class-level counter/dtbl
        // (valid on land-free ranks) so every rank decides identically -- the
        // FillBoundary at the end of this routine is collective.
        Real NOAH_time = static_cast<Real>(m_itimestep-1) * m_dtbl;
        if (elapsed_time < NOAH_time) {
            m_updated = false;
            return;
        }

        // We are updating
        m_updated = true;

        // Advance the counter once per firing, in lockstep on every rank.
        m_itimestep += 1;

        Box domain = m_geom.Domain();

        Print () << "Noah-MP driver at level " << lev << " started at time step: " << nstep+1 << std::endl;

        bool is_moist = (cons_in.nComp() > RhoQ1_comp);

        int klo = domain.smallEnd(2);

        // -------------------------------------------------------------------------
        // Precipitation forcing into Noah-MP (RAINBL / MP_RAINNC / MP_SNOW / MP_GRAUP /
        // SR). Works for ANY microphysics scheme through the typed source interface:
        // precip_sources carries a borrowed const view of each scheme-native cumulative
        // accumulator plus its native->kg/m^2 (== water-equivalent mm) factor. Slots are
        // scheme-specific -- some expose a `total` accumulator with `snow`/`graupel`
        // subsets (Morrison, WSM6), others expose an explicit `rain` accumulator with
        // density-scaled frozen species (SAM), and warm-rain schemes expose `rain` only
        // (Kessler). Empty sources (moisture/precip off) -> RAINBL stays 0, land runs
        // precip-free. The interval precip = (current - snapshot)*factor over the land-
        // call interval (WRF convention; Noah-MP divides by DTBL internally).
        //
        // Collect the five typed slots into a fixed-index array (NoahmpPrecipSlot order)
        // so the snapshot bookkeeping and the device delta kernel can loop generically.
        const SurfacePrecipAccumulationSource slot_src[NoahmpPrecipSlot::NumSlots] = {
            precip_sources.total, precip_sources.rain, precip_sources.snow,
            precip_sources.graupel, precip_sources.hail };
        const MultiFab* slot_accum[NoahmpPrecipSlot::NumSlots];
        Real            slot_factor[NoahmpPrecipSlot::NumSlots];
        bool            slot_has[NoahmpPrecipSlot::NumSlots];
        for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
            slot_accum[s]  = slot_src[s].accum;
            slot_factor[s] = slot_src[s].native_to_kg_m2;
            slot_has[s]    = (slot_src[s].accum != nullptr);
        }
        const bool have_precip = surface_precip_has_any_source(precip_sources);
        if (have_precip) {
            // Lazily allocate the per-level, per-slot "previous cumulative accum" snapshot.
            // On the first call (cold start OR restart) seed each present slot to its
            // current accumulation so the first interval's delta is 0 (no spurious precip
            // spike from a restored or non-zero initial accumulation).
            if (int(m_precip_accum_prev.size()) <= lev) { m_precip_accum_prev.resize(lev+1); }
            auto& prev = m_precip_accum_prev[lev];
            for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
                if (slot_has[s] && prev[s] == nullptr) {
                    // NOTE: allocate with ZERO ghost cells and setVal(0) first so no cell is
                    // ever uninitialized (an uninitialized prev-snapshot produced garbage
                    // ~1e25 deltas -> Noah-MP water-balance abort). We only ever read the
                    // valid region (surface slab), so 0 ghost cells is sufficient and avoids
                    // any ghost/nGrow mismatch in the Copy.
                    prev[s] = std::make_unique<MultiFab>(
                        slot_accum[s]->boxArray(), slot_accum[s]->DistributionMap(), 1, 0);
                    prev[s]->setVal(0.0);
                    MultiFab::Copy(*prev[s], *slot_accum[s], 0, 0, 1, 0);
                }
            }
        }

        // Loop over blocks to copy forcing data to Noahmp, drive the land model,
        // and copy data back to ERF Multifabs.
        int idb = 0;
        for (MFIter mfi(cons_in); mfi.isValid(); ++mfi, ++idb) {

            Box bx  = mfi.tilebox();
            Box gbx = mfi.tilebox(IntVect(0,0,0),IntVect(1,1,0));

            // Check if tile is at the lower boundary in lower z direction
            if (bx.smallEnd(2) != klo) { continue; }

            bx.makeSlab(2,klo);
            gbx.makeSlab(2,klo);

            // For limiting when populating ghost cells
            int i_lo = bx.smallEnd(0); int i_hi = bx.bigEnd(0);
            int j_lo = bx.smallEnd(1); int j_hi = bx.bigEnd(1);

            NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

            const Array4<const Real>& U_PHY  = xvel_in.const_array(mfi);
            const Array4<const Real>& V_PHY  = yvel_in.const_array(mfi);
            const Array4<const Real>& CONS   = cons_in.const_array(mfi);

            // Into NOAH-MP
            const Array4<const Real>& SWDOWN = lsm_fab_data[LsmData_NOAHMP::sw_flux_dn]->const_array(mfi);
            const Array4<const Real>& GLW    = lsm_fab_data[LsmData_NOAHMP::lw_flux_dn]->const_array(mfi);
            const Array4<const Real>& COSZEN = lsm_fab_data[LsmData_NOAHMP::cos_zenith_angle]->const_array(mfi);

            // Out of NOAH-MP
            Array4<Real> TSK           = lsm_fab_data[LsmData_NOAHMP::t_sfc]->array(mfi);
            Array4<Real> EMISS         = lsm_fab_data[LsmData_NOAHMP::sfc_emis]->array(mfi);
            Array4<Real> ALBSFCDIR_VIS = lsm_fab_data[LsmData_NOAHMP::sfc_alb_dir_vis]->array(mfi);
            Array4<Real> ALBSFCDIR_NIR = lsm_fab_data[LsmData_NOAHMP::sfc_alb_dir_nir]->array(mfi);
            Array4<Real> ALBSFCDIF_VIS = lsm_fab_data[LsmData_NOAHMP::sfc_alb_dif_vis]->array(mfi);
            Array4<Real> ALBSFCDIF_NIR = lsm_fab_data[LsmData_NOAHMP::sfc_alb_dif_nir]->array(mfi);
            // Land return-term diagnostic destinations (output only)
            Array4<Real> GRDFLX_o    = lsm_fab_data[LsmData_NOAHMP::grdflx]->array(mfi);
            Array4<Real> FIRA_o      = lsm_fab_data[LsmData_NOAHMP::fira]->array(mfi);
            Array4<Real> SAV_o       = lsm_fab_data[LsmData_NOAHMP::sav]->array(mfi);
            Array4<Real> SAG_o       = lsm_fab_data[LsmData_NOAHMP::sag]->array(mfi);
            Array4<Real> ALBEDO_o    = lsm_fab_data[LsmData_NOAHMP::albedo]->array(mfi);
            Array4<Real> SFCRUNOFF_o = lsm_fab_data[LsmData_NOAHMP::sfcrunoff]->array(mfi);
            Array4<Real> UDRUNOFF_o  = lsm_fab_data[LsmData_NOAHMP::udrunoff]->array(mfi);
            Array4<Real> SMSTAV_o    = lsm_fab_data[LsmData_NOAHMP::smstav]->array(mfi);
            Array4<Real> SMSTOT_o    = lsm_fab_data[LsmData_NOAHMP::smstot]->array(mfi);

            // Per-layer soil output fabs (3 groups x m_nsoil). Gather their Array4s
            // into a device-visible vector so the scatter kernel below can loop over
            // an arbitrary NSOIL. Layout matches soil_data_idx / soil_out_idx:
            // group g in {0:smois,1:sh2o,2:tslb}, layer k in [0,nsoil).
            const int nsoil = m_nsoil;
            const int n_soil_fld = 3*nsoil;
            Gpu::DeviceVector<Array4<Real>> soil_arr_d(n_soil_fld);
            {
                Vector<Array4<Real>> soil_arr_h(n_soil_fld);
                for (int g(0); g < 3; ++g) {
                    for (int k(0); k < nsoil; ++k) {
                        soil_arr_h[g*nsoil+k] = lsm_fab_data[soil_data_idx(g,k)]->array(mfi);
                    }
                }
                Gpu::copyAsync(Gpu::hostToDevice, soil_arr_h.begin(), soil_arr_h.end(), soil_arr_d.begin());
                Gpu::streamSynchronize();
            }
            Array4<Real>* soil_arr = soil_arr_d.data();
            const int soil_out_base = NoahmpOutputComp::NumComps; // soil outputs start here

            // NOTE: Need to expose stresses and get stresses from NOAHMP
            Array4<Real> q_flux_arr    = lsm_fab_flux[LsmFlux_NOAHMP::q_flux]->array(mfi);
            Array4<Real> t_flux_arr    = lsm_fab_flux[LsmFlux_NOAHMP::t_flux]->array(mfi);
            Array4<Real> tau13_arr     = lsm_fab_flux[LsmFlux_NOAHMP::tau13]->array(mfi);
            Array4<Real> tau23_arr     = lsm_fab_flux[LsmFlux_NOAHMP::tau23]->array(mfi);

            // Use The_Pinned_Arena() for host-accessible memory that can be used with GPU
            Array4<Real> noah_input_arr  =  noahmp_input_tmp[idb]->array();
            Array4<Real> noah_output_arr =  noahmp_output_tmp[idb]->array();

            // Per-slot precip accumulation views (current + previous-call snapshot) and
            // their native->kg/m^2 factors, for building the water-equivalent interval
            // precip. These MultiFabs live in the (device) default arena, so they are read
            // ONLY inside the device ParallelFor below and staged into the pinned
            // noah_input_arr — never dereferenced on the host (that segfaults on GPU).
            // have_precip is false for schemes with no precip (RAINBL stays 0). GpuArrays
            // are captured by value into the device kernel; absent slots keep null Array4s
            // that are guarded by slot_present before any dereference.
            const bool hp = have_precip;
            GpuArray<Array4<const Real>, NoahmpPrecipSlot::NumSlots> accum_now, accum_prv;
            GpuArray<Real, NoahmpPrecipSlot::NumSlots> accum_fac;
            GpuArray<int,  NoahmpPrecipSlot::NumSlots> slot_present;
            for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
                slot_present[s] = (hp && slot_has[s]) ? 1 : 0;
                accum_fac[s]    = slot_factor[s];
                if (slot_present[s]) {
                    accum_now[s] = slot_accum[s]->const_array(mfi);
                    accum_prv[s] = m_precip_accum_prev[lev][s]->const_array(mfi);
                }
            }
            const int kklo = klo;

            // Copy forcing data from ERF to Noahmp (device; stage into pinned buffer).
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real qv = (is_moist) ? CONS(i,j,k,RhoQ1_comp)/CONS(i,j,k,Rho_comp) : zero;
                noah_input_arr(i,j,0,NoahmpInputComp::u_phy)   = myhalf*(U_PHY(i,j,k)+U_PHY(i+1,j,k));
                noah_input_arr(i,j,0,NoahmpInputComp::v_phy)   = myhalf*(V_PHY(i,j,k)+V_PHY(i  ,j+1,k));
                noah_input_arr(i,j,0,NoahmpInputComp::t_phy)   = getTgivenRandRTh(CONS(i,j,k,Rho_comp),CONS(i,j,k,RhoTheta_comp),qv);
                noah_input_arr(i,j,0,NoahmpInputComp::qv_curr) = qv;
                noah_input_arr(i,j,0,NoahmpInputComp::p8w)     = getPgivenRTh(CONS(i,j,k,RhoTheta_comp),qv);
                noah_input_arr(i,j,0,NoahmpInputComp::swdown)  = SWDOWN(i,j,0);
                noah_input_arr(i,j,0,NoahmpInputComp::glw)     = GLW(i,j,0);
                noah_input_arr(i,j,0,NoahmpInputComp::coszen)  = COSZEN(i,j,0);

                // Water-equivalent interval precip [mm] for this land call (WRF convention;
                // Noah-MP divides by DTBL -> mm/s). Per slot: d = (current - snapshot) *
                // native_to_kg_m2, so each species is converted to water-equivalent mm
                // BEFORE it is combined (kg/m^2 == mm of water). Accumulators are on the
                // surface slab kklo.
                Real drain = zero, sr = zero, dsnow = zero, dgraup = zero;
                if (hp) {
                    Real dd[NoahmpPrecipSlot::NumSlots];
                    for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
                        dd[s] = slot_present[s]
                              ? amrex::max(zero, (accum_now[s](i,j,kklo) - accum_prv[s](i,j,kklo)) * accum_fac[s])
                              : zero;
                    }
                    // Frozen species -> MP_SNOW / MP_GRAUP (hail folded into graupel;
                    // Noah-MP has no hail input). Matches WRF MP_SNOW / MP_GRAUP.
                    dsnow  = dd[NoahmpPrecipSlot::snow];
                    dgraup = dd[NoahmpPrecipSlot::graupel] + dd[NoahmpPrecipSlot::hail];
                    const Real dfroz = dsnow + dgraup;
                    // Total non-convective precip -> RAINBL / MP_RAINNC. Prefer the scheme's
                    // `total` accumulator; otherwise reconstruct from explicit species
                    // (rain + frozen), mirroring the IO layer's derivation. This is the
                    // correct total for SAM/Kessler (rain slot, no total), where the old
                    // "rain_accum is total" assumption undercounted frozen precip.
                    drain = slot_present[NoahmpPrecipSlot::total]
                          ? dd[NoahmpPrecipSlot::total]
                          : (dd[NoahmpPrecipSlot::rain] + dfroz);
                    // Physical guard: precip accumulated over one land interval cannot
                    // exceed a sane ceiling (lsm_max_precip_interval). Protects Noah-MP's
                    // water-balance check from any corrupt prev-snapshot delta -- a larger
                    // value can only come from noisy IC/BC, never real precip.
                    drain = amrex::min(drain, lsm_max_precip_interval);
                    sr = (drain > zero) ? amrex::min(Real(1.0), dfroz/drain) : zero;
                }
                noah_input_arr(i,j,0,NoahmpInputComp::rainbl)    = drain;
                noah_input_arr(i,j,0,NoahmpInputComp::sr_in)     = sr;
                // WRF-faithful precip breakdown fed straight to NoahMP's MP_* forcing so the
                // rain/snow partition (opt_snf=4) matches WRF; opt_snf=1 ignores MP_*.
                noah_input_arr(i,j,0,NoahmpInputComp::mp_rainnc) = drain;
                noah_input_arr(i,j,0,NoahmpInputComp::mp_snow)   = dsnow;
                noah_input_arr(i,j,0,NoahmpInputComp::mp_graup)  = dgraup;
            });

            // Synchronize to ensure GPU kernel is complete before host access
            Gpu::streamSynchronize();

            // Now on the host, copy the pinned staged data to NoahmpIO arrays
            LoopOnCpu(bx, [&] (int i, int j, int ) noexcept
            {
                noahmpio->U_PHY(i,1,j)   = noah_input_arr(i,j,0,NoahmpInputComp::u_phy);
                noahmpio->V_PHY(i,1,j)   = noah_input_arr(i,j,0,NoahmpInputComp::v_phy);
                noahmpio->T_PHY(i,1,j)   = noah_input_arr(i,j,0,NoahmpInputComp::t_phy);
                noahmpio->QV_CURR(i,1,j) = noah_input_arr(i,j,0,NoahmpInputComp::qv_curr);
                noahmpio->P8W(i,1,j)     = noah_input_arr(i,j,0,NoahmpInputComp::p8w);
                noahmpio->SWDOWN(i,j)    = noah_input_arr(i,j,0,NoahmpInputComp::swdown);
                noahmpio->GLW(i,j)       = noah_input_arr(i,j,0,NoahmpInputComp::glw);
                noahmpio->COSZEN(i,j)    = noah_input_arr(i,j,0,NoahmpInputComp::coszen);
                noahmpio->RAINBL(i,j)    = noah_input_arr(i,j,0,NoahmpInputComp::rainbl);  // [mm]
                noahmpio->SR(i,j)        = noah_input_arr(i,j,0,NoahmpInputComp::sr_in);   // [-]
                noahmpio->MP_RAINNC(i,j) = noah_input_arr(i,j,0,NoahmpInputComp::mp_rainnc); // [mm]
                noahmpio->MP_SNOW(i,j)   = noah_input_arr(i,j,0,NoahmpInputComp::mp_snow);   // [mm]
                noahmpio->MP_GRAUP(i,j)  = noah_input_arr(i,j,0,NoahmpInputComp::mp_graup);  // [mm]
            });

            // Call the noahmpio driver code. This runs the land model forcing for
            // each object in noahmpio_vect that represent a block in the domain.
            // Mirror the authoritative counter into each block.
            noahmpio->itimestep = m_itimestep;
            noahmpio->DriverMain();

            // The soil-layer count was fixed in Init (erf.lsm_nsoil, asserted to match
            // the namelist NSOIL there), so the per-layer reads below are always in range.

            // Copy results from NoahmpIO back to temporary arrays
            LoopOnCpu(bx, [&] (int i, int j, int ) noexcept
            {
                noah_output_arr(i,j,0,NoahmpOutputComp::hfx)           = noahmpio->HFX(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::lh)            = noahmpio->LH(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::tau_ew)        = noahmpio->TAU_EW(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::tau_ns)        = noahmpio->TAU_NS(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::tsk)           = noahmpio->TSK(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::emiss)         = noahmpio->EMISS(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdir_vis) = noahmpio->ALBSFCDIRXY(i,1,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdir_nir) = noahmpio->ALBSFCDIRXY(i,2,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdif_vis) = noahmpio->ALBSFCDIFXY(i,1,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdif_nir) = noahmpio->ALBSFCDIFXY(i,2,j);
                // Land return-term diagnostics (output only). SMSTAV/SMSTOT are not
                // computed by the Noah-MP core in this driver (read back 0).
                noah_output_arr(i,j,0,NoahmpOutputComp::o_grdflx)      = noahmpio->GRDFLX(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::o_fira)        = noahmpio->FIRAXY(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::o_sav)         = noahmpio->SAVXY(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::o_sag)         = noahmpio->SAGXY(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::o_albedo)      = noahmpio->ALBEDO(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::o_sfcrunoff)   = noahmpio->SFCRUNOFF(i,j);
                noah_output_arr(i,j,0,NoahmpOutputComp::o_udrunoff)    = noahmpio->UDRUNOFF(i,j);
                // SMSTAV/SMSTOT are intentionally NOT staged here: this Noah-MP core does
                // not compute them, so they are emitted as lsm_undefined downstream rather
                // than as a misleading physical 0 (see the scatter kernel below).
                // 3D soil-profile state, one 2D component per soil layer, for any NSOIL.
                // Fortran layer index is 1-based; output components are laid out as
                // [NumComps .. +nsoil) smois, [+nsoil .. +2*nsoil) sh2o, [+2*nsoil ..) tslb.
                for (int k(0); k < nsoil; ++k) {
                    noah_output_arr(i,j,0,soil_out_base + 0*nsoil + k) = noahmpio->SMOIS(i,k+1,j);
                    noah_output_arr(i,j,0,soil_out_base + 1*nsoil + k) = noahmpio->SH2O (i,k+1,j);
                    noah_output_arr(i,j,0,soil_out_base + 2*nsoil + k) = noahmpio->TSLB (i,k+1,j);
                }
            });

            // Copy forcing data from Noahmp to ERF
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Limit indices to the valid box. FillBoundary will pick these up below.
                int ii = std::min(std::max(i,i_lo),i_hi);
                int jj = std::min(std::max(j,j_lo),j_hi);

                // SurfaceLayer fluxes at CC.
                // Noah-MP returns the -9999 fill value for cells it does NOT process
                // (sea-ice / open-water points, which still have LANDMASK=1). Applying
                // that as a flux gives -9999/(rho*Cp) ~ -7.6 K*m/s and crashes the
                // lowest cell to ~200 K. Detect the fill and instead write the
                // lsm_undefined sentinel; the surface layer then falls back to
                // the MOST flux for those cells (see ERF_SurfaceLayer.cpp).
                // NOTE: tau13/tau23 are nodal in xz/yz; the 2D MFs have 1 ghost cell
                //       so the surface layer can average them.
                Real hfx_lsm = noah_output_arr(ii,jj,0,NoahmpOutputComp::hfx);
                if (hfx_lsm > Real(-9990.0)) {
                    // Convert Noah-MP surface fluxes to the KINEMATIC form ERF's
                    // surface layer stores for its MOST parameter updates (u*, t*, L).
                    // The LSM value is interchanged with the MOST value <x'w'> in the
                    // same array; ERF's surface-layer/PBL convention is the KINEMATIC
                    // MOST form (no rho):
                    //   t_flux = <theta' w'>   [K m/s]   (ERF_MOSTStress.H, surf_temp_flux)
                    //   q_flux = <qv' w'>      [kg/kg m/s]
                    //   tau13/23 = u*^2 comps  [m2/s2]   (SurfaceLayer.cpp: u*=sqrt(tau))
                    // Noah-MP returns HFX,LH [W/m2] and TAU_EW/NS [N/m2] (rho INCLUDED),
                    // so we divide by rho. rho is dry-air density (~1.14 kg/m3).
                    //
                    // NOTE on Exner: strictly, ERF's energy variable is POTENTIAL
                    // temperature, so the exact conversion of the Noah-MP sensible-heat
                    // flux (built from ACTUAL temperatures) to a theta-flux would carry an
                    // Exner factor:  theta-flux = HFX/(rho*Cp) * (p0/p)^(Rd/Cp).
                    // We deliberately OMIT it here to match WRF, which applies HFX/(rho*Cp)
                    // directly into its potential-temperature PBL equation (WRF
                    // module_pbl_driver.F: b_t = hfx/rho/CP, no Exner). This is an
                    // approximation relative to ERF's exact potential-temperature flux, not
                    // an identity: the error scales with (p0/p)^(Rd/Cp)-1, which is small
                    // near sea level (p ~= p0 -> Exner ~= 1.00-1.01, <~1%) but grows over
                    // low-pressure / elevated terrain (e.g. p ~ 800 hPa -> ~6-7%). Reinstate
                    // by multiplying t_flux by std::pow(p_0/getPgivenRTh(
                    // CONS(ii,jj,k,RhoTheta_comp),
                    // (is_moist)?CONS(ii,jj,k,RhoQ1_comp)/CONS(ii,jj,k,Rho_comp):zero),
                    // R_d/Cp_d) if that bias matters for the case.
                    Real rho_l  = CONS(ii,jj,k,Rho_comp);
                    t_flux_arr(i,j,k) = hfx_lsm/(rho_l*Cp_d);
                    q_flux_arr(i,j,k) = noah_output_arr(ii,jj,0,NoahmpOutputComp::lh)/(rho_l*L_v);
                    tau13_arr(i,j,k)  = noah_output_arr(ii,jj,0,NoahmpOutputComp::tau_ew)/rho_l;
                    tau23_arr(i,j,k)  = noah_output_arr(ii,jj,0,NoahmpOutputComp::tau_ns)/rho_l;
                } else {
                    t_flux_arr(i,j,k) = lsm_undefined;
                    q_flux_arr(i,j,k) = lsm_undefined;
                    tau13_arr(i,j,k)  = lsm_undefined;
                    tau23_arr(i,j,k)  = lsm_undefined;
                }

                // RRTMGP variables
                if (hfx_lsm > Real(-9990.0)) {
                    TSK(i,j,0)           = noah_output_arr(ii,jj,0,NoahmpOutputComp::tsk);
                    EMISS(i,j,0)         = noah_output_arr(ii,jj,0,NoahmpOutputComp::emiss);
                    ALBSFCDIR_VIS(i,j,0) = noah_output_arr(ii,jj,0,NoahmpOutputComp::albsfcdir_vis);
                    ALBSFCDIR_NIR(i,j,0) = noah_output_arr(ii,jj,0,NoahmpOutputComp::albsfcdir_nir);
                    ALBSFCDIF_VIS(i,j,0) = noah_output_arr(ii,jj,0,NoahmpOutputComp::albsfcdif_vis);
                    ALBSFCDIF_NIR(i,j,0) = noah_output_arr(ii,jj,0,NoahmpOutputComp::albsfcdif_nir);
                    // Land return-term diagnostics (output only; not fed to dynamics)
                    GRDFLX_o(i,j,0)      = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_grdflx);
                    FIRA_o(i,j,0)        = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_fira);
                    SAV_o(i,j,0)         = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_sav);
                    SAG_o(i,j,0)         = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_sag);
                    ALBEDO_o(i,j,0)      = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_albedo);
                    SFCRUNOFF_o(i,j,0)   = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_sfcrunoff);
                    UDRUNOFF_o(i,j,0)    = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_udrunoff);
                    // SMSTAV/SMSTOT are not computed by this Noah-MP core, so emit the
                    // lsm_undefined sentinel rather than the uncomputed core value: a
                    // physical-looking 0 would misrepresent them as valid diagnostics.
                    SMSTAV_o(i,j,0)      = lsm_undefined;
                    SMSTOT_o(i,j,0)      = lsm_undefined;
                    // Per-layer soil profile (any NSOIL): copy each output component
                    // into its dedicated 2D fab (soil_arr laid out group-major).
                    for (int s(0); s < n_soil_fld; ++s) {
                        soil_arr[s](i,j,0) = noah_output_arr(ii,jj,0,soil_out_base + s);
                    }
                } else {
                    TSK(i,j,0)           = lsm_undefined;
                    EMISS(i,j,0)         = lsm_undefined;
                    ALBSFCDIR_VIS(i,j,0) = lsm_undefined;
                    ALBSFCDIR_NIR(i,j,0) = lsm_undefined;
                    ALBSFCDIF_VIS(i,j,0) = lsm_undefined;
                    ALBSFCDIF_NIR(i,j,0) = lsm_undefined;
                    GRDFLX_o(i,j,0)      = lsm_undefined;
                    FIRA_o(i,j,0)        = lsm_undefined;
                    SAV_o(i,j,0)         = lsm_undefined;
                    SAG_o(i,j,0)         = lsm_undefined;
                    ALBEDO_o(i,j,0)      = lsm_undefined;
                    SFCRUNOFF_o(i,j,0)   = lsm_undefined;
                    UDRUNOFF_o(i,j,0)    = lsm_undefined;
                    SMSTAV_o(i,j,0)      = lsm_undefined;
                    SMSTOT_o(i,j,0)      = lsm_undefined;
                    // Per-layer soil profile: undefined sentinel on non-land cells.
                    for (int s(0); s < n_soil_fld; ++s) {
                        soil_arr[s](i,j,0) = lsm_undefined;
                    }
                }
            });
            // The scatter kernel reads soil_arr_d (a per-iteration DeviceVector);
            // synchronize before it is destroyed at the end of this MFIter body.
            Gpu::streamSynchronize();
        }

        // Advance the precip "previous cumulative accum" snapshot to the current values
        // now that this land call's RAINBL/SR have been consumed, so the next call
        // differences against this step. Device-safe whole-MultiFab copy (these live in
        // the device arena; never touch them on the host).
        if (have_precip) {
            for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
                if (slot_has[s]) {
                    MultiFab::Copy(*m_precip_accum_prev[lev][s], *slot_accum[s], 0, 0, 1, 0);
                }
            }
        }

        // Fill the ghost cells
        for (auto ivar = 0; ivar < LsmFlux_NOAHMP::NumVars; ++ivar) {
            lsm_fab_flux[ivar]->FillBoundary(m_geom.periodicity());
        }
        Print () << "Noah-MP driver at level " << lev << " completed" << std::endl;
    } // lev == 0
};
