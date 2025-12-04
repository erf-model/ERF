
#include<iostream>
#include<string>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include <ERF_NOAHMP.H>
#include <ERF_Constants.H>
#include <ERF_EOS.H>

using namespace amrex;

/* Initialize lsm data structures */
void
NOAHMP::Init (const int& lev,
              const MultiFab& cons_in,
              const Geometry& geom,
              const Real& dt)
{

    m_dt   = dt;
    m_geom = geom;

    Box domain = geom.Domain();
    khi_lsm    = domain.smallEnd(2) - 1;

    LsmDataMap.resize(m_lsm_data_size);
    LsmDataMap = {LsmData_NOAHMP::t_sfc           , LsmData_NOAHMP::sfc_emis       ,
                  LsmData_NOAHMP::sfc_alb_dir_vis , LsmData_NOAHMP::sfc_alb_dir_nir,
                  LsmData_NOAHMP::sfc_alb_dif_vis , LsmData_NOAHMP::sfc_alb_dif_nir,
                  LsmData_NOAHMP::cos_zenith_angle, LsmData_NOAHMP::sw_flux_dn     ,
                  LsmData_NOAHMP::lw_flux_dn                                        };
    LsmDataName.resize(m_lsm_data_size);
    LsmDataName = {"t_sfc"           , "sfc_emis"        ,
                   "sfc_alb_dir_vis" , "sfc_alb_dir_nir" ,
                   "sfc_alb_dif_vis" , "sfc_alb_dif_nir" ,
                   "cos_zenith_angle", "sw_flux_dn"      ,
                   "lw_flux_dn"      };


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
    for (auto& b : bl_lsm) {
        b.setRange(2,0);
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

    // Create the data
    for (auto ivar = 0; ivar < LsmData_NOAHMP::NumVars; ++ivar) {
        // State vars are CC
        lsm_fab_data[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);

        // NOTE: Radiation steps first so we set values
        //       to reasonable initialization for coupling
        Real val_to_set = 0.0;
        if (ivar == LsmData_NOAHMP::t_sfc) {
            val_to_set = 300.0;
        } else if (ivar == LsmData_NOAHMP::sfc_emis) {
            val_to_set = 0.9;
        } else if ( (ivar>=LsmData_NOAHMP::sfc_alb_dir_vis) &&
                    (ivar<=LsmData_NOAHMP::sfc_alb_dif_nir) ) {
            val_to_set = 0.06;
        } else {
            val_to_set = 0.0;
        }
        lsm_fab_data[ivar]->setVal(val_to_set);
    }

    // Create the fluxes
    for (auto ivar = 0; ivar < LsmFlux_NOAHMP::NumVars; ++ivar) {
        // NOTE: Fluxes are CC with ghost cells for averaging
        lsm_fab_flux[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, IntVect(1,1,0));
        lsm_fab_flux[ivar]->setVal(0.);
    }

    Print() << "Noah-MP initialization started" << std::endl;

    // Set noahmpio_vect to the size of local blocks (boxes)
    noahmpio_vect.resize(cons_in.local_size(), lev);

    // Iterate over multifab and noahmpio object together. Multifabs is
    // used to extract size of blocks and set bounds for noahmpio objects.
    int idb = 0;
    for (MFIter mfi(cons_in, false); mfi.isValid(); ++mfi, ++idb) {

        // Get bounds for the tile
        const Box& bx = mfi.tilebox();

        // Check if tile is at the lower boundary in lower z direction
        if (bx.smallEnd(2) != domain.smallEnd(2)) { continue; }

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

  Print() << "Noah-MP initialization completed" << std::endl;

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
                            const Real& dt,
                            const int& nstep)
{

    Box domain = m_geom.Domain();

    Print () << "Noah-MP driver started at time step: " << nstep+1 << std::endl;

    // Loop over blocks to copy forcing data to Noahmp, drive the land model,
    // and copy data back to ERF Multifabs.
    int idb = 0;
    for (MFIter mfi(cons_in, false); mfi.isValid(); ++mfi, ++idb) {

        Box bx  = mfi.tilebox();
        Box gbx = mfi.tilebox(IntVect(0,0,0),IntVect(1,1,0));

        // Check if tile is at the lower boundary in lower z direction
        if (bx.smallEnd(2) != domain.smallEnd(2)) { continue; }

        bx.makeSlab(2,domain.smallEnd(2));
        gbx.makeSlab(2,domain.smallEnd(2));

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

        // NOTE: Need to expose stresses and get stresses from NOAHMP
        Array4<Real> q_flux_arr    = lsm_fab_flux[LsmFlux_NOAHMP::q_flux]->array(mfi);
        Array4<Real> t_flux_arr    = lsm_fab_flux[LsmFlux_NOAHMP::t_flux]->array(mfi);
        Array4<Real> tau13_arr     = lsm_fab_flux[LsmFlux_NOAHMP::tau13]->array(mfi);
        Array4<Real> tau23_arr     = lsm_fab_flux[LsmFlux_NOAHMP::tau23]->array(mfi);

        // Create temporary BaseFabs that will be accessible on host
        // Use The_Pinned_Arena() for host-accessible memory that can be used with GPU
        FArrayBox tmp_u_phy(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_v_phy(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_t_phy(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_qv_curr(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_p8w(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_swdown(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_glw(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_coszen(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_hfx(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_lh(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_tau_ew(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_tau_ns(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_tsk(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_emiss(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_albsfcdir_vis(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_albsfcdir_nir(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_albsfcdif_vis(bx, 1, The_Pinned_Arena());
        FArrayBox tmp_albsfcdif_nir(bx, 1, The_Pinned_Arena());

        // Get array views
        auto const& tmp_u_phy_arr   = tmp_u_phy.array();
        auto const& tmp_v_phy_arr   = tmp_v_phy.array();
        auto const& tmp_t_phy_arr   = tmp_t_phy.array();
        auto const& tmp_qv_curr_arr = tmp_qv_curr.array();
        auto const& tmp_p8w_arr     = tmp_p8w.array();
        auto const& tmp_swdown_arr  = tmp_swdown.array();
        auto const& tmp_glw_arr     = tmp_glw.array();
        auto const& tmp_coszen_arr  = tmp_coszen.array();

        // Copy forcing data from ERF to Noahmp.
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
        {
            tmp_u_phy_arr(i,j,0)   = 0.5*(U_PHY(i,j,0)+U_PHY(i+1,j  ,0));
            tmp_v_phy_arr(i,j,0)   = 0.5*(V_PHY(i,j,0)+V_PHY(i  ,j+1,0));
            tmp_t_phy_arr(i,j,0)   = getTgivenRandRTh(CONS(i,j,0,Rho_comp),CONS(i,j,0,RhoTheta_comp));
            tmp_qv_curr_arr(i,j,0) = CONS(i,j,0,RhoQ1_comp)/CONS(i,j,0,Rho_comp);
            tmp_p8w_arr(i,j,0)     = getPgivenRTh(CONS(i,j,0,RhoTheta_comp));
            tmp_swdown_arr(i,j,0)  = SWDOWN(i,j,0);
            tmp_glw_arr(i,j,0)     = GLW(i,j,0);
            tmp_coszen_arr(i,j,0)  = COSZEN(i,j,0);
        });

        // Synchronize to ensure GPU kernel is complete before host access
        Gpu::streamSynchronize();

        // Now on the host, copy data to NoahmpIO arrays
        // Use LoopOnCpu for CPU-side parallelization
        const auto& h_u_arr      = tmp_u_phy.const_array();
        const auto& h_v_arr      = tmp_v_phy.const_array();
        const auto& h_t_arr      = tmp_t_phy.const_array();
        const auto& h_qv_arr     = tmp_qv_curr.const_array();
        const auto& h_p8w_arr    = tmp_p8w.const_array();
        const auto& h_swdown_arr = tmp_swdown.const_array();
        const auto& h_glw_arr    = tmp_glw.const_array();
        const auto& h_coszen_arr = tmp_coszen.const_array();

        LoopOnCpu(bx, [&] (int i, int j, int ) noexcept
        {
            noahmpio->U_PHY(i,1,j)   = h_u_arr(i,j,0);
            noahmpio->V_PHY(i,1,j)   = h_v_arr(i,j,0);
            noahmpio->T_PHY(i,1,j)   = h_t_arr(i,j,0);
            noahmpio->QV_CURR(i,1,j) = h_qv_arr(i,j,0);
            noahmpio->P8W(i,1,j)     = h_p8w_arr(i,j,0);
            noahmpio->SWDOWN(i,j)    = h_swdown_arr(i,j,0);
            noahmpio->GLW(i,j)       = h_glw_arr(i,j,0);
            noahmpio->COSZEN(i,j)    = h_coszen_arr(i,j,0);
        });

        // Call the noahmpio driver code. This runs the land model forcing for
        // each object in noahmpio_vect that represent a block in the domain.
        noahmpio->itimestep = nstep+1;
        noahmpio->DriverMain();

        // Copy results from NoahmpIO back to temporary arrays
        auto h_hfx_arr           = tmp_hfx.array();
        auto h_lh_arr            = tmp_lh.array();
        auto h_tau_ew_arr        = tmp_tau_ew.array();
        auto h_tau_ns_arr        = tmp_tau_ns.array();
        auto h_tsk_arr           = tmp_tsk.array();
        auto h_emiss_arr         = tmp_emiss.array();
        auto h_albsfcdir_vis_arr = tmp_albsfcdir_vis.array();
        auto h_albsfcdir_nir_arr = tmp_albsfcdir_nir.array();
        auto h_albsfcdif_vis_arr = tmp_albsfcdif_vis.array();
        auto h_albsfcdif_nir_arr = tmp_albsfcdif_nir.array();

        LoopOnCpu(bx, [&] (int i, int j, int ) noexcept
        {
            h_hfx_arr(i,j,0)           = noahmpio->HFX(i,j);
            h_lh_arr(i,j,0)            = noahmpio->LH(i,j);
            h_tau_ew_arr(i,j,0)        = noahmpio->TAU_EW(i,j);
            h_tau_ns_arr(i,j,0)        = noahmpio->TAU_NS(i,j);
            h_tsk_arr(i,j,0)           = noahmpio->TSK(i,j);
            h_emiss_arr(i,j,0)         = noahmpio->EMISS(i,j);
            h_albsfcdir_vis_arr(i,j,0) = noahmpio->ALBSFCDIRXY(i,1,j);
            h_albsfcdir_nir_arr(i,j,0) = noahmpio->ALBSFCDIRXY(i,2,j);
            h_albsfcdif_vis_arr(i,j,0) = noahmpio->ALBSFCDIFXY(i,1,j);
            h_albsfcdif_nir_arr(i,j,0) = noahmpio->ALBSFCDIFXY(i,2,j);
        });

        // Copy results back to output arrays
        auto const& tmp_hfx_arr           = tmp_hfx.array();
        auto const& tmp_lh_arr            = tmp_lh.array();
        auto const& tmp_tau_ew_arr        = tmp_tau_ew.array();
        auto const& tmp_tau_ns_arr        = tmp_tau_ns.array();
        auto const& tmp_tsk_arr           = tmp_tsk.array();
        auto const& tmp_emiss_arr         = tmp_emiss.array();
        auto const& tmp_albsfcdir_vis_arr = tmp_albsfcdir_vis.array();
        auto const& tmp_albsfcdir_nir_arr = tmp_albsfcdir_nir.array();
        auto const& tmp_albsfcdif_vis_arr = tmp_albsfcdif_vis.array();
        auto const& tmp_albsfcdif_nir_arr = tmp_albsfcdif_nir.array();


        // Copy forcing data from Noahmp to ERF
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
        {
            // SurfaceLayer fluxes at CC
            t_flux_arr(i,j,0)    = tmp_hfx_arr(i,j,0)/(CONS(i,j,0,Rho_comp)*Cp_d);
            q_flux_arr(i,j,0)    = tmp_lh_arr(i,j,0)/(CONS(i,j,0,Rho_comp)*L_v);

            // NOTE: The following fluxes are nodal in xz/yz.
            //       The 2D MFs have 1 ghost cell so we can average these
            //       when using them in the surface layer class.
            tau13_arr(i,j,0)  = tmp_tau_ew_arr(i,j,0)/CONS(i,j,0,Rho_comp);
            tau23_arr(i,j,0)  = tmp_tau_ns_arr(i,j,0)/CONS(i,j,0,Rho_comp);

            // RRTMGP variables
            TSK(i,j,0)           = tmp_tsk_arr(i,j,0);
            EMISS(i,j,0)         = tmp_emiss_arr(i,j,0);
            ALBSFCDIR_VIS(i,j,0) = tmp_albsfcdir_vis_arr(i,j,0);
            ALBSFCDIR_NIR(i,j,0) = tmp_albsfcdir_nir_arr(i,j,0);
            ALBSFCDIF_VIS(i,j,0) = tmp_albsfcdif_vis_arr(i,j,0);
            ALBSFCDIF_NIR(i,j,0) = tmp_albsfcdif_nir_arr(i,j,0);
        });
    }

    // Fill the ghost cells
    for (auto ivar = 0; ivar < LsmFlux_NOAHMP::NumVars; ++ivar) {
        lsm_fab_flux[ivar]->FillBoundary(m_geom.periodicity());
    }
    Print () << "Noah-MP driver completed" << std::endl;
};
