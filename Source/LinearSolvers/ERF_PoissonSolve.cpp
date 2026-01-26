#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

void
solve_with_mlmg    (int lev,
                    Vector<amrex::MultiFab>& rhs, Vector<MultiFab>& p,
                    Vector<amrex::Array<MultiFab,AMREX_SPACEDIM>>& fluxes,
                    const Geometry& geom,
                    const amrex::Vector<amrex::IntVect>& ref_ratio,
                    Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
                    int mg_verbose, Real reltol, Real abstol);
void
solve_with_EB_mlmg (int lev,
                    Vector<amrex::MultiFab>& rhs, Vector<MultiFab>& p,
                    Vector<amrex::Array<MultiFab,AMREX_SPACEDIM>>& fluxes,
                    EBFArrayBoxFactory const& ebfact,
                    eb_aux_ const& ebfact_u,
                    eb_aux_ const& ebfact_v,
                    eb_aux_ const& ebfact_w,
                    const Geometry& geom,
                    const amrex::Vector<amrex::IntVect>& ref_ratio,
                    Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
                    int mg_verbose, Real reltol, Real abstol);

/**
 * Project the single-level velocity field to enforce the anelastic constraint
 * Note that the level may or may not be level 0.
 */
void ERF::project_initial_velocity (int lev, Real time, Real l_dt)
{
    BL_PROFILE("ERF::project_initial_velocity()");
    // Impose FillBoundary on density since we use it in the conversion of velocity to momentum
    vars_new[lev][Vars::cons].FillBoundary(geom[lev].periodicity());

    const MultiFab* c_vfrac = nullptr;
    if (solverChoice.terrain_type == TerrainType::EB) {
        c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
    }

    VelocityToMomentum(vars_new[lev][Vars::xvel], IntVect{0},
                       vars_new[lev][Vars::yvel], IntVect{0},
                       vars_new[lev][Vars::zvel], IntVect{0},
                       vars_new[lev][Vars::cons],
                       rU_new[lev], rV_new[lev], rW_new[lev],
                       Geom(lev).Domain(), domain_bcs_type, c_vfrac);

    Vector<MultiFab> tmp_mom;

    tmp_mom.push_back(MultiFab(vars_new[lev][Vars::cons],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rU_new[lev],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rV_new[lev],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rW_new[lev],make_alias,0,1));

    // If at lev > 0 we must first fill the velocities at the c/f interface -- this must
    //    be done *after* the projection at lev-1
    if (lev > 0) {
        int levc = lev-1;

        const MultiFab* c_vfrac_crse = nullptr;
        if (solverChoice.terrain_type == TerrainType::EB) {
            c_vfrac_crse = &((get_eb(levc).get_const_factory())->getVolFrac());
        }

        MultiFab& S_new_crse = vars_new[levc][Vars::cons];
        MultiFab& U_new_crse = vars_new[levc][Vars::xvel];
        MultiFab& V_new_crse = vars_new[levc][Vars::yvel];
        MultiFab& W_new_crse = vars_new[levc][Vars::zvel];

        VelocityToMomentum(U_new_crse, IntVect{0}, V_new_crse, IntVect{0}, W_new_crse, IntVect{0}, S_new_crse,
                           rU_new[levc], rV_new[levc], rW_new[levc],
                           Geom(levc).Domain(), domain_bcs_type, c_vfrac_crse);

        rU_new[levc].FillBoundary(geom[levc].periodicity());
        FPr_u[levc].RegisterCoarseData({&rU_new[levc], &rU_new[levc]}, {time, time+l_dt});

        rV_new[levc].FillBoundary(geom[levc].periodicity());
        FPr_v[levc].RegisterCoarseData({&rV_new[levc], &rV_new[levc]}, {time, time+l_dt});

        rW_new[levc].FillBoundary(geom[levc].periodicity());
        FPr_w[levc].RegisterCoarseData({&rW_new[levc], &rW_new[levc]}, {time, time+l_dt});
    }

    Real l_time = 0.0;
    project_momenta(lev, l_time, l_dt, tmp_mom);

    MomentumToVelocity(vars_new[lev][Vars::xvel],
                       vars_new[lev][Vars::yvel],
                       vars_new[lev][Vars::zvel],
                       vars_new[lev][Vars::cons],
                       rU_new[lev], rV_new[lev], rW_new[lev],
                       Geom(lev).Domain(), domain_bcs_type, c_vfrac);
 }

/**
 * Project the single-level momenta to enforce the anelastic constraint
 * Note that the level may or may not be level 0.
 */
void ERF::project_momenta (int lev, Real l_time, Real l_dt, Vector<MultiFab>& mom_mf)
{
    BL_PROFILE("ERF::project_momenta()");
    //
    // If at lev > 0 we must first fill the momenta at the c/f interface with interpolated coarse values
    //
    if (lev > 0) {
        PhysBCFunctNoOp null_bc;
        FPr_u[lev-1].FillSet(mom_mf[IntVars::xmom], l_time, null_bc, domain_bcs_type);
        FPr_v[lev-1].FillSet(mom_mf[IntVars::ymom], l_time, null_bc, domain_bcs_type);
        FPr_w[lev-1].FillSet(mom_mf[IntVars::zmom], l_time, null_bc, domain_bcs_type);
    }

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(mom_mf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(mom_mf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

    Box domain = geom[lev].Domain();

    MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp, 1);

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;

    if (solverChoice.terrain_type == TerrainType::EB)
    {
        rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0, MFInfo(), EBFactory(lev));
        phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1, MFInfo(), EBFactory(lev));
    } else {
        rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
        phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1);
    }

    MultiFab rhs_lev(rhs[0], make_alias, 0, 1);
    MultiFab phi_lev(phi[0], make_alias, 0, 1);

    auto dx    = geom[lev].CellSizeArray();
    auto dxInv = geom[lev].InvCellSizeArray();

    // Inflow on an x-face -- note only the normal velocity is used in the projection
    if (domain_bc_type[0] == "Inflow" || domain_bc_type[3] == "Inflow") {
        (*physbcs_u[lev])(vars_new[lev][Vars::xvel],vars_new[lev][Vars::xvel],vars_new[lev][Vars::yvel],
                        IntVect{1,0,0},t_new[lev],BCVars::xvel_bc,false);
    }

    // Inflow on a  y-face -- note only the normal velocity is used in the projection
    if (domain_bc_type[1] == "Inflow" || domain_bc_type[4] == "Inflow") {
        (*physbcs_v[lev])(vars_new[lev][Vars::yvel],vars_new[lev][Vars::xvel],vars_new[lev][Vars::yvel],
                          IntVect{0,1,0},t_new[lev],BCVars::yvel_bc,false);
    }

    if (domain_bc_type[0] == "Inflow" || domain_bc_type[3] == "Inflow" ||
        domain_bc_type[1] == "Inflow" || domain_bc_type[4] == "Inflow") {

        const MultiFab* c_vfrac = nullptr;
        if (solverChoice.terrain_type == TerrainType::EB) {
            c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
        }

        VelocityToMomentum(vars_new[lev][Vars::xvel], IntVect{0},
                            vars_new[lev][Vars::yvel], IntVect{0},
                            vars_new[lev][Vars::zvel], IntVect{0},
                            vars_new[lev][Vars::cons],
                            mom_mf[IntVars::xmom],
                            mom_mf[IntVars::ymom],
                            mom_mf[IntVars::zmom],
                            Geom(lev).Domain(),
                            domain_bcs_type, c_vfrac);
    }

    // If !fixed_density, we must convert (rho u) which came in
    // to (rho0 u) which is what we will project
    if (!solverChoice.fixed_density[lev]) {
        ConvertForProjection(mom_mf[Vars::cons], r_hse,
                             mom_mf[IntVars::xmom],
                             mom_mf[IntVars::ymom],
                             mom_mf[IntVars::zmom],
                             Geom(lev).Domain(),
                             domain_bcs_type);
    }

    //
    // ****************************************************************************
    // Now convert the rho0w MultiFab to hold Omega rather than rhow
    // ****************************************************************************
    //
    if (solverChoice.mesh_type == MeshType::VariableDz)
    {
        for ( MFIter mfi(rhs_lev,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Array4<Real const>& rho0u_arr = mom_mf[IntVars::xmom].const_array(mfi);
            const Array4<Real const>& rho0v_arr = mom_mf[IntVars::ymom].const_array(mfi);
            const Array4<Real      >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);

            const Array4<Real const>&     z_nd = z_phys_nd[lev]->const_array(mfi);
            const Array4<Real const>&     mf_u =  mapfac[lev][MapFacType::u_x]->const_array(mfi);
            const Array4<Real const>&     mf_v =  mapfac[lev][MapFacType::v_y]->const_array(mfi);

            //
            // Define Omega from (rho0 W) but store it in the same array
            //
            Box tbz = mfi.nodaltilebox(2);
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (k == 0) {
                    rho0w_arr(i,j,k) = Real(0.0);
                } else {
                    Real rho0w = rho0w_arr(i,j,k);
                    rho0w_arr(i,j,k) = OmegaFromW(i,j,k,rho0w,
                                                  rho0u_arr,rho0v_arr,
                                                  mf_u,mf_v,z_nd,dxInv);
                }
            });
        } // mfi
    }

    // ****************************************************************************
    // Allocate fluxes
    // ****************************************************************************
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;
    fluxes.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (solverChoice.terrain_type == TerrainType::EB) {
            fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0, MFInfo(), EBFactory(lev));
        } else {
            fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
        }
    }

    // ****************************************************************************
    // Initialize phi to 0
    // (It is essential that we do this in order to fill the corners; these are never
    //  used but the Saxpy requires the values to be initialized.)
    // ****************************************************************************
    phi_lev.setVal(0.0);

    // ****************************************************************************
    // Break into subdomains
    // ****************************************************************************

    std::map<int,int> index_map;

    BoxArray ba(grids[lev]);

    Vector<MultiFab> rhs_sub; rhs_sub.resize(1);
    Vector<MultiFab> phi_sub; phi_sub.resize(1);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> fluxes_sub; fluxes_sub.resize(1);

    MultiFab ax_sub, ay_sub, az_sub, dJ_sub, znd_sub;
    MultiFab mfmx_sub, mfmy_sub;

    Array<MultiFab,AMREX_SPACEDIM> rho0_u_sub;
    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;

    // If we are going to solve with MLMG then we do not need to break this into subdomains
    bool will_solve_with_mlmg = false;
    if (solverChoice.mesh_type == MeshType::ConstantDz) {
        will_solve_with_mlmg = true;
#ifdef ERF_USE_FFT
        if (use_fft) {
            bool all_boxes_ok = true;
            for (int isub = 0; isub < subdomains[lev].size(); ++isub) {
                Box my_region(subdomains[lev][isub].minimalBox());
                bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
                if (!boxes_make_rectangle) {
                    all_boxes_ok = false;
                }
            } // isub
            if (all_boxes_ok) {
                will_solve_with_mlmg = false;
            }
        } // use_fft
#else
        if (use_fft) {
            amrex::Warning("You set use_fft=true but didn't build with USE_FFT = TRUE; defaulting to MLMG");
        }
#endif
    } // No terrain or grid stretching

    for (int isub = 0; isub < subdomains[lev].size(); ++isub)
    {
        BoxList bl_sub;
        Vector<int> dm_sub;

        for (int j = 0; j < ba.size(); j++)
        {
            if (subdomains[lev][isub].intersects(ba[j]))
            {
                //
                // Note that bl_sub.size() is effectively a counter which is
                // incremented above
                //
                // if (ParallelDescriptor::MyProc() == j) {
                // }
                index_map[bl_sub.size()] = j;

                bl_sub.push_back(grids[lev][j]);
                dm_sub.push_back(dmap[lev][j]);
            } // intersects
        } // loop over ba (j)

        BoxArray ba_sub(bl_sub);

        BoxList bl2d_sub = ba_sub.boxList();
        for (auto& b : bl2d_sub) {
            b.setRange(2,0);
        }
        BoxArray ba2d_sub(std::move(bl2d_sub));

        // Define MultiFabs that hold only the data in this particular subdomain
        if (solverChoice.terrain_type == TerrainType::EB) {
            if (ba_sub != ba) {
                amrex::Print() << "EB Solves with multiple regions is not yet supported" << std::endl;
            }
            rhs_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, rhs_lev.nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            phi_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, phi_lev.nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));

            mfmx_sub.define(ba2d_sub, DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::m_x]->nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            mfmy_sub.define(ba2d_sub, DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::m_y]->nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
              dJ_sub.define(ba_sub, DistributionMapping(dm_sub), 1, detJ_cc[lev]->nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                fluxes_sub[0][idim].define(convert(ba_sub, IntVect::TheDimensionVector(idim)), DistributionMapping(dm_sub), 1,
                                                   IntVect::TheZeroVector(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            }
            rho0_u_sub[0].define(convert(ba_sub, IntVect::TheDimensionVector(0)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::xmom].nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            rho0_u_sub[1].define(convert(ba_sub, IntVect::TheDimensionVector(1)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::ymom].nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            rho0_u_sub[2].define(convert(ba_sub, IntVect::TheDimensionVector(2)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::zmom].nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
        } else {
            rhs_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, rhs_lev.nGrowVect(), MFInfo{}.SetAlloc(false));
            phi_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, phi_lev.nGrowVect(), MFInfo{}.SetAlloc(false));

            mfmx_sub.define(ba2d_sub, DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::m_x]->nGrowVect(), MFInfo{}.SetAlloc(false));
            mfmy_sub.define(ba2d_sub, DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::m_y]->nGrowVect(), MFInfo{}.SetAlloc(false));
              dJ_sub.define(ba_sub, DistributionMapping(dm_sub), 1, detJ_cc[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                fluxes_sub[0][idim].define(convert(ba_sub, IntVect::TheDimensionVector(idim)), DistributionMapping(dm_sub), 1,
                                                   IntVect::TheZeroVector(), MFInfo{}.SetAlloc(false));
            }
            rho0_u_sub[0].define(convert(ba_sub, IntVect::TheDimensionVector(0)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::xmom].nGrowVect(), MFInfo{}.SetAlloc(false));
            rho0_u_sub[1].define(convert(ba_sub, IntVect::TheDimensionVector(1)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::ymom].nGrowVect(), MFInfo{}.SetAlloc(false));
            rho0_u_sub[2].define(convert(ba_sub, IntVect::TheDimensionVector(2)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::zmom].nGrowVect(), MFInfo{}.SetAlloc(false));
        }

        // Link the new MultiFabs to the FABs in the original MultiFabs (no copy required)
        for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi)
        {
            int orig_index = index_map[mfi.index()];
            rhs_sub[0].setFab(mfi, FArrayBox(rhs_lev[orig_index], amrex::make_alias, 0, 1));
            phi_sub[0].setFab(mfi, FArrayBox(phi_lev[orig_index], amrex::make_alias, 0, 1));

            mfmx_sub.setFab(mfi, FArrayBox((*mapfac[lev][MapFacType::m_x])[orig_index], amrex::make_alias, 0, 1));
            mfmy_sub.setFab(mfi, FArrayBox((*mapfac[lev][MapFacType::m_y])[orig_index], amrex::make_alias, 0, 1));

            fluxes_sub[0][0].setFab(mfi,FArrayBox(fluxes[0][0][orig_index], amrex::make_alias, 0, 1));
            fluxes_sub[0][1].setFab(mfi,FArrayBox(fluxes[0][1][orig_index], amrex::make_alias, 0, 1));
            fluxes_sub[0][2].setFab(mfi,FArrayBox(fluxes[0][2][orig_index], amrex::make_alias, 0, 1));

            rho0_u_sub[0].setFab(mfi,FArrayBox(mom_mf[IntVars::xmom][orig_index], amrex::make_alias, 0, 1));
            rho0_u_sub[1].setFab(mfi,FArrayBox(mom_mf[IntVars::ymom][orig_index], amrex::make_alias, 0, 1));
            rho0_u_sub[2].setFab(mfi,FArrayBox(mom_mf[IntVars::zmom][orig_index], amrex::make_alias, 0, 1));
        }

        rho0_u_const[0] = &rho0_u_sub[0];
        rho0_u_const[1] = &rho0_u_sub[1];
        rho0_u_const[2] = &rho0_u_sub[2];

        if (solverChoice.mesh_type != MeshType::ConstantDz) {
            ax_sub.define(convert(ba_sub,IntVect(1,0,0)), DistributionMapping(dm_sub), 1,
                          ax[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            ay_sub.define(convert(ba_sub,IntVect(0,1,0)), DistributionMapping(dm_sub), 1,
                          ay[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            az_sub.define(convert(ba_sub,IntVect(0,0,1)), DistributionMapping(dm_sub), 1,
                          az[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            znd_sub.define(convert(ba_sub,IntVect(1,1,1)), DistributionMapping(dm_sub), 1,
                           z_phys_nd[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));

            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi) {
                int orig_index = index_map[mfi.index()];
                ax_sub.setFab(mfi, FArrayBox((*ax[lev])[orig_index], amrex::make_alias, 0, 1));
                ay_sub.setFab(mfi, FArrayBox((*ay[lev])[orig_index], amrex::make_alias, 0, 1));
                az_sub.setFab(mfi, FArrayBox((*az[lev])[orig_index], amrex::make_alias, 0, 1));
                znd_sub.setFab(mfi, FArrayBox((*z_phys_nd[lev])[orig_index], amrex::make_alias, 0, 1));
                dJ_sub.setFab(mfi, FArrayBox((*detJ_cc[lev])[orig_index], amrex::make_alias, 0, 1));
            }
        }

        if (solverChoice.terrain_type == TerrainType::EB) {
            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi) {
                int orig_index = index_map[mfi.index()];
                dJ_sub.setFab(mfi, FArrayBox((*detJ_cc[lev])[orig_index], amrex::make_alias, 0, 1));
            }
        }

        // ****************************************************************************
        // Compute divergence which will form RHS
        // Note that we replace "rho0w" with the contravariant momentum, Omega
        // ****************************************************************************

        compute_divergence(lev, rhs_sub[0], rho0_u_const, geom_tmp[0]);

        Real rhsnorm;

        // Max norm over the entire MultiFab
        rhsnorm = rhs_sub[0].norm0();

        if (mg_verbose > 0) {
            bool local = false;
            Real sum = volWgtSumMF(lev,rhs_sub[0],0,dJ_sub,mfmx_sub,mfmy_sub,false,local);
            Print() << "Max/L2 norm of divergence before solve in subdomain " << isub << " at level " << lev << " : " << rhsnorm << " " <<
                        rhs_sub[0].norm2() << " and volume-weighted sum " << sum << std::endl;
        }

        if (lev == 0 && solverChoice.use_real_bcs)
        {
            // We always use VariableDz if use_real_bcs is true
            AMREX_ALWAYS_ASSERT(solverChoice.mesh_type == MeshType::VariableDz);

            // Note that we always impose the projections one level at a time so this will always be a vector of length 1
            Array<MultiFab*, AMREX_SPACEDIM> rho0_u_vec =
               {&mom_mf[IntVars::xmom], &mom_mf[IntVars::ymom], &mom_mf[IntVars::zmom]};
            Array<MultiFab*, AMREX_SPACEDIM> area_vec = {ax[lev].get(), ay[lev].get(), az[lev].get()};
            //
            // Modify ax,ay,ax to include the map factors as used in the divergence calculation
            // We do this here so that it is seen in the call to enforceInOutSolvability
            //
            for (MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box xbx = mfi.nodaltilebox(0);
                Box ybx = mfi.nodaltilebox(1);
                Box zbx = mfi.nodaltilebox(2);
                const Array4<Real      >& ax_ar = ax[lev]->array(mfi);
                const Array4<Real      >& ay_ar = ay[lev]->array(mfi);
                const Array4<Real      >& az_ar = az[lev]->array(mfi);
                const Array4<Real const>& mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);
                const Array4<Real const>& mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                ParallelFor(xbx,ybx,zbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ax_ar(i,j,k) /= mf_uy(i,j,0);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ay_ar(i,j,k) /= mf_vx(i,j,0);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    az_ar(i,j,k) /= (mf_mx(i,j,0)*mf_my(i,j,0));
                });
            } // mfi

            if (mg_verbose > 0) {
                Print() << "Calling enforceInOutSolvability" << std::endl;
            }
            enforceInOutSolvability(lev, rho0_u_vec, area_vec, geom[lev]);

            //
            // Return ax,ay,ax to their original definition
            //
            for (MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box xbx = mfi.nodaltilebox(0);
                Box ybx = mfi.nodaltilebox(1);
                Box zbx = mfi.nodaltilebox(2);
                const Array4<Real      >& ax_ar = ax[lev]->array(mfi);
                const Array4<Real      >& ay_ar = ay[lev]->array(mfi);
                const Array4<Real      >& az_ar = az[lev]->array(mfi);
                const Array4<Real const>& mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);
                const Array4<Real const>& mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                ParallelFor(xbx,ybx,zbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ax_ar(i,j,k) *= mf_uy(i,j,0);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ay_ar(i,j,k) *= mf_vx(i,j,0);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    az_ar(i,j,k) *= (mf_mx(i,j,0)*mf_my(i,j,0));
                });
            } // mfi

            compute_divergence(lev, rhs_lev, rho0_u_const, geom_tmp[0]);

            // Re-define max norm over the entire MultiFab
            rhsnorm = rhs_lev.norm0();

            if (mg_verbose > 0)
            {
                bool local = false;
                Real sum = volWgtSumMF(lev,rhs_sub[0],0,dJ_sub,mfmx_sub,mfmy_sub,false,local);
                Print() << "Max/L2 norm of divergence before solve at level " << lev << " : " << rhsnorm << " " <<
                            rhs_lev.norm2() << " and volume-weighted sum " << sum << std::endl;
            }
        } // lev 0 && use_real_bcs

        // *******************************************************************************************
        // Enforce solvability if the problem is singular (i.e all sides Neumann or periodic)
        // Note that solves at lev > 0 are always singular because we impose Neumann bc's on all sides
        // *******************************************************************************************
        bool is_singular = true;
        if (lev == 0) {
            if ( (domain_bc_type[0] == "Outflow" || domain_bc_type[0] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
            if ( (domain_bc_type[1] == "Outflow" || domain_bc_type[1] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
            if ( (domain_bc_type[3] == "Outflow" || domain_bc_type[3] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
            if ( (domain_bc_type[4] == "Outflow" || domain_bc_type[4] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
            if ( (domain_bc_type[5] == "Outflow" || domain_bc_type[5] == "Open")                               ) is_singular = false;
        } else {
            Box my_region(subdomains[lev][isub].minimalBox());
            if ( (domain_bc_type[5] == "Outflow" || domain_bc_type[5] == "Open") && (my_region.bigEnd(2) == domain.bigEnd(2)) ) is_singular = false;
        }

        if (is_singular)
        {
            bool local = false;
            Real sum = volWgtSumMF(lev,rhs_sub[0],0,dJ_sub,mfmx_sub,mfmy_sub,false,local);

            Real vol;
            if (solverChoice.mesh_type == MeshType::ConstantDz) {
                vol = rhs_sub[0].boxArray().numPts();
            } else {
                vol = dJ_sub.sum();
            }

            sum /= (vol * dx[0] * dx[1] * dx[2]);

            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi)
            {
                rhs_sub[0][mfi.index()].template minus<RunOn::Device>(sum);
            }
            if (mg_verbose > 0) {
                amrex::Print() << " Subtracting " << sum << " from rhs in subdomain " << isub << std::endl;

                sum = volWgtSumMF(lev,rhs_sub[0],0,dJ_sub,mfmx_sub,mfmy_sub,false,local);
                Print() << "Sum after subtraction " << sum << " in subdomain " << isub << std::endl;
            }

        } // if is_singular

        rhsnorm = rhs_sub[0].norm0();

        // ****************************************************************************
        // No need to build the solver if RHS == 0
        // ****************************************************************************
        if (rhsnorm <= solverChoice.poisson_abstol) return;

        Real start_step = static_cast<Real>(ParallelDescriptor::second());

        if (mg_verbose > 0) {
            amrex::Print() << " Solving in subdomain " << isub << " of " << subdomains[lev].size() << " bins at level " << lev << std::endl;
        }

        if (solverChoice.mesh_type == MeshType::VariableDz) {
            //
            // Modify ax,ay,ax to include the map factors as used in the divergence calculation
            // We do this here to set the coefficients used in the stencil -- the extra factor
            // of the mapfac comes from the gradient
            //
            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi)
            {
                Box xbx = mfi.nodaltilebox(0);
                Box ybx = mfi.nodaltilebox(1);
                Box zbx = mfi.nodaltilebox(2);
                const Array4<Real      >& ax_ar = ax_sub.array(mfi);
                const Array4<Real      >& ay_ar = ay_sub.array(mfi);
                const Array4<Real      >& az_ar = az_sub.array(mfi);
                const Array4<Real const>& mf_ux = mapfac[lev][MapFacType::u_x]->const_array(mfi);
                const Array4<Real const>& mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);
                const Array4<Real const>& mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
                const Array4<Real const>& mf_vy = mapfac[lev][MapFacType::v_y]->const_array(mfi);
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                ParallelFor(xbx,ybx,zbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ax_ar(i,j,k) *= (mf_ux(i,j,0) / mf_uy(i,j,0));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ay_ar(i,j,k) *= (mf_vy(i,j,0) / mf_vx(i,j,0));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    az_ar(i,j,k) /= (mf_mx(i,j,0)*mf_my(i,j,0));
                });
            } // mfi
        }

        if (solverChoice.terrain_type != TerrainType::EB) {

#ifdef ERF_USE_FFT
        Box my_region(subdomains[lev][isub].minimalBox());
#endif

        // ****************************************************************************
        // No terrain or grid stretching
        // ****************************************************************************
        if (solverChoice.mesh_type == MeshType::ConstantDz) {
            if (will_solve_with_mlmg) {
                solve_with_mlmg(lev, rhs_sub, phi_sub, fluxes_sub, geom[lev], ref_ratio, domain_bc_type,
                                mg_verbose, solverChoice.poisson_reltol, solverChoice.poisson_abstol);
            } else {
#ifdef ERF_USE_FFT
                solve_with_fft(lev, my_region, rhs_sub[0], phi_sub[0], fluxes_sub[0]);
#endif
            }
        } // No terrain or grid stretching
        // ****************************************************************************
        // Grid stretching (flat terrain)
        // ****************************************************************************
        else if (solverChoice.mesh_type == MeshType::StretchedDz) {
#ifndef ERF_USE_FFT
            amrex::Abort("Rebuild with USE_FFT = TRUE so you can use the FFT solver");
#else
            bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
            if (!boxes_make_rectangle) {
                amrex::Abort("FFT won't work unless the union of boxes is rectangular");
            } else {
                if (!use_fft) {
                    amrex::Warning("Using FFT even though you didn't set use_fft to true; it's the best choice");
                }
                solve_with_fft(lev, my_region, rhs_sub[0], phi_sub[0], fluxes_sub[0]);
            }
#endif
        } // grid stretching

        // ****************************************************************************
        // General terrain
        // ****************************************************************************
        else if (solverChoice.mesh_type == MeshType::VariableDz) {
#ifdef ERF_USE_FFT
            bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
            if (!boxes_make_rectangle) {
                amrex::Abort("FFT preconditioner for GMRES won't work unless the union of boxes is rectangular");
            } else {
                solve_with_gmres(lev, my_region, rhs_sub[0], phi_sub[0], fluxes_sub[0], ax_sub, ay_sub, az_sub, dJ_sub, znd_sub);
            }
#else
            amrex::Abort("Rebuild with USE_FFT = TRUE so you can use the FFT preconditioner for GMRES");
#endif

            //
            // Restore ax,ay,ax to their original definitions
            //
            for (MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box xbx = mfi.nodaltilebox(0);
                Box ybx = mfi.nodaltilebox(1);
                Box zbx = mfi.nodaltilebox(2);
                const Array4<Real      >& ax_ar = ax_sub.array(mfi);
                const Array4<Real      >& ay_ar = ay_sub.array(mfi);
                const Array4<Real      >& az_ar = az_sub.array(mfi);
                const Array4<Real const>& mf_ux = mapfac[lev][MapFacType::u_x]->const_array(mfi);
                const Array4<Real const>& mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);
                const Array4<Real const>& mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
                const Array4<Real const>& mf_vy = mapfac[lev][MapFacType::v_y]->const_array(mfi);
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                ParallelFor(xbx,ybx,zbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ax_ar(i,j,k) *= (mf_uy(i,j,0) / mf_ux(i,j,0));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ay_ar(i,j,k) *= (mf_vx(i,j,0) / mf_vy(i,j,0));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    az_ar(i,j,k) *= (mf_mx(i,j,0)*mf_my(i,j,0));
                });
            } // mfi

        } // MeshType::VariableDz

        // ****************************************************************************
        // Print time in solve
        // ****************************************************************************
        Real end_step = static_cast<Real>(ParallelDescriptor::second());
        if (mg_verbose > 0) {
            amrex::Print() << "Time in solve " << end_step - start_step << std::endl;
        }

        } // not EB
    } // loop over subdomains (i)

    // ****************************************************************************
    // When using multigrid we can solve for all of the level at once, even if there
    //      are disjoint regions
    // ****************************************************************************
    if (solverChoice.terrain_type == TerrainType::EB) {
        Real start_step_eb = static_cast<Real>(ParallelDescriptor::second());
        solve_with_EB_mlmg(lev, rhs_sub, phi_sub, fluxes_sub,
                           *(get_eb(lev).get_const_factory()),
                           *(get_eb(lev).get_u_const_factory()),
                           *(get_eb(lev).get_v_const_factory()),
                           *(get_eb(lev).get_w_const_factory()),
                           geom[lev], ref_ratio, domain_bc_type,
                           mg_verbose, solverChoice.poisson_reltol, solverChoice.poisson_abstol);
        Real end_step_eb = static_cast<Real>(ParallelDescriptor::second());
        if (mg_verbose > 0) {
            amrex::Print() << "Time in solve " << end_step_eb - start_step_eb << std::endl;
        }
    }

    // ****************************************************************************
    // Subtract dt grad(phi) from the momenta (rho0u, rho0v, Omega)
    // ****************************************************************************
    MultiFab::Add(mom_mf[IntVars::xmom],fluxes[0][0],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::ymom],fluxes[0][1],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::zmom],fluxes[0][2],0,0,1,0);

    // ****************************************************************************
    // Define gradp from fluxes -- note that fluxes is dt * change in Gp
    //   (weighted by map factor!)
    // ****************************************************************************
    MultiFab::Saxpy(gradp[lev][GpVars::gpx],-1.0/l_dt,fluxes[0][0],0,0,1,0);
    MultiFab::Saxpy(gradp[lev][GpVars::gpy],-1.0/l_dt,fluxes[0][1],0,0,1,0);
    MultiFab::Saxpy(gradp[lev][GpVars::gpz],-1.0/l_dt,fluxes[0][2],0,0,1,0);

    gradp[lev][GpVars::gpx].FillBoundary(geom_tmp[0].periodicity());
    gradp[lev][GpVars::gpy].FillBoundary(geom_tmp[0].periodicity());
    gradp[lev][GpVars::gpz].FillBoundary(geom_tmp[0].periodicity());

    //
    // This call is only to verify the divergence after the solve
    // It is important we do this before computing the rho0w_arr from Omega back to rho0w
    //
    // ****************************************************************************
    // THIS IS SIMPLY VERIFYING THE DIVERGENCE AFTER THE SOLVE
    // ****************************************************************************
    //
    if (mg_verbose > 0)
    {
        rho0_u_const[0] = &mom_mf[IntVars::xmom];
        rho0_u_const[1] = &mom_mf[IntVars::ymom];
        rho0_u_const[2] = &mom_mf[IntVars::zmom];

        compute_divergence(lev, rhs_lev, rho0_u_const, geom_tmp[0]);

        bool local = false;
        Real sum = volWgtSumMF(lev,rhs_lev,0,*detJ_cc[lev],*mapfac[lev][MapFacType::m_x],*mapfac[lev][MapFacType::m_y],false,local);

        if (mg_verbose > 0) {
            Print() << "Max/L2 norm of divergence after  solve at level " << lev << " : " << rhs_lev.norm0() << " " <<
                        rhs_lev.norm2() << " and volume-weighted sum " << sum << std::endl;
        }

#if 0
         // FOR DEBUGGING ONLY
         for ( MFIter mfi(rhs_lev,TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Array4<Real const>& rhs_arr = rhs_lev.const_array(mfi);
            Box bx = mfi.validbox();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (std::abs(rhs_arr(i,j,k)) > 1.e-10) {
                    amrex::AllPrint() << "RHS after solve at " <<
                                          IntVect(i,j,k) << " " << rhs_arr(i,j,k) << std::endl;
                }
            });
         } // mfi
#endif

    } // mg_verbose

    //
    // ****************************************************************************
    // Now convert the rho0w MultiFab back to holding (rho0w) rather than Omega
    // ****************************************************************************
    //
    if (solverChoice.mesh_type == MeshType::VariableDz)
    {
        for (MFIter mfi(mom_mf[Vars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
             Box tbz = mfi.nodaltilebox(2);
             const Array4<Real      >& rho0u_arr = mom_mf[IntVars::xmom].array(mfi);
             const Array4<Real      >& rho0v_arr = mom_mf[IntVars::ymom].array(mfi);
             const Array4<Real      >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);
             const Array4<Real const>&      z_nd = z_phys_nd[lev]->const_array(mfi);
             const Array4<Real const>&      mf_u =  mapfac[lev][MapFacType::u_x]->const_array(mfi);
             const Array4<Real const>&      mf_v =  mapfac[lev][MapFacType::v_y]->const_array(mfi);
             ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                 Real omega = rho0w_arr(i,j,k);
                 rho0w_arr(i,j,k) = WFromOmega(i,j,k,omega,
                                               rho0u_arr,rho0v_arr,
                                               mf_u,mf_v,z_nd,dxInv);
             });
        } // mfi
    }

    // If !fixed_density, we must convert (rho0 u) back
    // to (rho0 u) which is what we will pass back out
    if (!solverChoice.fixed_density[lev]) {
        ConvertForProjection(r_hse, mom_mf[Vars::cons],
                             mom_mf[IntVars::xmom],
                             mom_mf[IntVars::ymom],
                             mom_mf[IntVars::zmom],
                             Geom(lev).Domain(),
                             domain_bcs_type);
    }

    // ****************************************************************************
    // Update pressure variable with phi -- note that phi is dt * change in pressure
    // ****************************************************************************
    MultiFab::Saxpy(pp_inc[lev], 1.0/l_dt, phi_lev,0,0,1,1);
}
