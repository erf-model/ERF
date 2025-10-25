#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Project the single-level velocity field to enforce the anelastic constraint
 * Note that the level may or may not be level 0.
 */
void ERF::project_velocity (int lev, Real l_dt)
{
    // Impose FillBoundary on density since we use it in the conversion of velocity to momentum
    vars_new[lev][Vars::cons].FillBoundary(geom[lev].periodicity());

    BL_PROFILE("ERF::project_velocity()");
    VelocityToMomentum(vars_new[lev][Vars::xvel], IntVect{0},
                       vars_new[lev][Vars::yvel], IntVect{0},
                       vars_new[lev][Vars::zvel], IntVect{0},
                       vars_new[lev][Vars::cons],
                       rU_new[lev], rV_new[lev], rW_new[lev],
                       Geom(lev).Domain(), domain_bcs_type);

    Vector<MultiFab> tmp_mom;

    tmp_mom.push_back(MultiFab(vars_new[lev][Vars::cons],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rU_new[lev],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rV_new[lev],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rW_new[lev],make_alias,0,1));

    project_momenta(lev, l_dt, tmp_mom);

    MomentumToVelocity(vars_new[lev][Vars::xvel],
                       vars_new[lev][Vars::yvel],
                       vars_new[lev][Vars::zvel],
                       vars_new[lev][Vars::cons],
                       rU_new[lev], rV_new[lev], rW_new[lev],
                       Geom(lev).Domain(), domain_bcs_type);
 }

/**
 * Project the single-level momenta to enforce the anelastic constraint
 * Note that the level may or may not be level 0.
 */
void ERF::project_momenta (int lev, Real l_dt, Vector<MultiFab>& mom_mf)
{
    BL_PROFILE("ERF::project_momenta()");

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(mom_mf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(mom_mf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

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
            VelocityToMomentum(vars_new[lev][Vars::xvel], IntVect{0},
                               vars_new[lev][Vars::yvel], IntVect{0},
                               vars_new[lev][Vars::zvel], IntVect{0},
                               vars_new[lev][Vars::cons],
                               mom_mf[IntVars::xmom],
                               mom_mf[IntVars::ymom],
                               mom_mf[IntVars::zmom],
                               Geom(lev).Domain(),
                               domain_bcs_type);
    }

    // If !fixed_density, we must convert (rho u) which came in
    // to (rho0 u) which is what we will project
    if (!solverChoice.fixed_density) {
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
    // Compute divergence which will form RHS
    // Note that we replace "rho0w" with the contravariant momentum, Omega
    // ****************************************************************************
    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;
    rho0_u_const[0] = &mom_mf[IntVars::xmom];
    rho0_u_const[1] = &mom_mf[IntVars::ymom];
    rho0_u_const[2] = &mom_mf[IntVars::zmom];

    compute_divergence(lev, rhs_lev, rho0_u_const, geom_tmp[0]);

    Real rhsnorm, sum;

    // Max norm over the entire MultiFab
    rhsnorm = rhs_lev.norm0();

    sum = volWgtSumMF(lev,rhs_lev,0,false);

    if (mg_verbose > 0) {
        Print() << "Max/L2 norm of divergence before solve at level " << lev << " : " << rhsnorm << " " <<
                    rhs_lev.norm2() << " and volume-weighted sum " << sum << std::endl;
    }

    if (lev == 0 && solverChoice.use_real_bcs) {
        // Note that we always impose the projections one level at a time so this will always be a vector of length 1
        Array<MultiFab*, AMREX_SPACEDIM> rho0_u_vec =
           {&mom_mf[IntVars::xmom], &mom_mf[IntVars::ymom], &mom_mf[IntVars::zmom]};
        Array<MultiFab*, AMREX_SPACEDIM> area_vec = {ax[lev].get(), ay[lev].get(), az[lev].get()};
#if 0
        //
        // Modify ax,ay,ax to include the map factors as used in the divergence calculation
        // We do this here so that it is seen in the call to enforceInOutSolvability
        //
        if (solverChoice.mesh_type == MeshType::VariableDz) {
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
        } // variable dz
#endif

        if (mg_verbose > 0) {
            Print() << "Calling enforceInOutSolvability" << std::endl;
        }
        enforceInOutSolvability(lev, rho0_u_vec, area_vec, geom[lev]);

#if 0
        //
        // Return ax,ay,ax to their original definition
        //
        if (solverChoice.mesh_type == MeshType::VariableDz) {
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
        } // variabledz
#endif
        compute_divergence(lev, rhs_lev, rho0_u_const, geom_tmp[0]);

        // Max norm over the entire MultiFab
        rhsnorm = rhs_lev.norm0();

        sum = volWgtSumMF(lev,rhs_lev,0,false);

        if (mg_verbose > 0) {
            Print() << "Max/L2 norm of divergence before solve at level " << lev << " : " << rhsnorm << " " <<
                        rhs_lev.norm2() << " and volume-weighted sum " << sum << std::endl;
        }
    } // lev 0 && use_real_bcs

    // ****************************************************************************
    // Enforce solvability if the problem is singular (i.e all sides Neumann or periodic)
    // ****************************************************************************
    bool is_singular = true;
    if ( (domain_bc_type[0] == "Outflow" || domain_bc_type[0] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
    if ( (domain_bc_type[1] == "Outflow" || domain_bc_type[1] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
    if ( (domain_bc_type[3] == "Outflow" || domain_bc_type[3] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
    if ( (domain_bc_type[4] == "Outflow" || domain_bc_type[4] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
    if ( (domain_bc_type[5] == "Outflow" || domain_bc_type[5] == "Open")                               ) is_singular = false;

    if (is_singular) {
        if (lev > 0)
        {
            Vector<Real> sum_sub; sum_sub.resize(subdomains[lev].size(),Real(0.));

            for (MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.validbox();
                for (int i = 0; i < subdomains[lev].size(); ++i) {
                    if (subdomains[lev][i].intersects(bx)) {
                        sum_sub[i] += rhs_lev[mfi.index()].template sum<RunOn::Device>(0);
                    }
                }
            }
            ParallelDescriptor::ReduceRealSum(sum_sub.data(), sum_sub.size());

            for (int i = 0; i < subdomains[lev].size(); ++i) {
                sum_sub[i] /= static_cast<Real>(subdomains[lev][i].numPts());
            }

            for ( MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.validbox();
                for (int i = 0; i < subdomains[lev].size(); ++i) {
                    if (subdomains[lev][i].intersects(bx)) {
                        rhs_lev[mfi.index()].template minus<RunOn::Device>(sum_sub[i]);
                        if (mg_verbose > 1) {
                            amrex::Print() << " Subtracting " << sum_sub[i] << " in " << rhs_lev[mfi.index()].box() << std::endl;
                        }
                    }
                }
            }
        } else {

            sum = volWgtSumMF(lev,rhs_lev,0,false);

            Real vol = detJ_cc[lev]->sum() / (dxInv[0] * dxInv[1] * dxInv[2]);

            Print() << "Vol wgt sum " << sum << std::endl;
            Print() << "Vol         " << vol << std::endl;

            sum /= vol;

            for ( MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.validbox();
                for (int i = 0; i < subdomains[lev].size(); ++i) {
                    rhs_lev[mfi.index()].template minus<RunOn::Device>(sum);
                    amrex::Print() << " Subtracting " << sum << " in " << rhs_lev[mfi.index()].box() << std::endl;
                }
            }

            sum = volWgtSumMF(lev,rhs_lev,0,false);

            Print() << "SUM AFTER SUBTRACTION " << sum << std::endl;
        }
    } // if is_singular

    // ****************************************************************************
    //
    // No need to build the solver if RHS == 0
    //
    if (rhsnorm <= solverChoice.poisson_abstol) return;
    // ****************************************************************************

    // ****************************************************************************
    // Initialize phi to 0
    // (It is essential that we do this in order to fill the corners; these are never
    //  used but the Saxpy requires the values to be initialized.)
    // ****************************************************************************
    phi[0].setVal(0.0);

    Real start_step = static_cast<Real>(ParallelDescriptor::second());

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
    // Choose the solver and solve
    // ****************************************************************************

    std::map<int,int> index_map;

    BoxArray ba(grids[lev]);

    Vector<MultiFab> rhs_sub; rhs_sub.resize(1);
    Vector<MultiFab> phi_sub; phi_sub.resize(1);
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes_sub; fluxes_sub.resize(1);

    MultiFab ax_sub, ay_sub, az_sub, dJ_sub, znd_sub;

    for (int isub = 0; isub < subdomains[lev].size(); ++isub)
    {
        if (mg_verbose > 0) {
            amrex::Print() << " Solving in subdomain " << isub << " of " << subdomains[lev].size() << " bins at level " << lev << std::endl;
        }

        BoxList bl_sub;
        Vector<int> dm_sub;

        for (int j = 0; j < ba.size(); j++)
        {
            if (subdomains[lev][isub].intersects(ba[j]))
            {
                // amrex::Print() <<" INTERSECTS I " << isub << " " << j << " " << grids[lev][j] << std::endl;
                //
                // Note that bl_sub.size() is effectively a counter which is
                // incremented above
                //
                // if (ParallelDescriptor::MyProc() == j) {
                // }
                index_map[bl_sub.size()] = j;

                // amrex::Print() <<" PUSHING BACK " << j << " " << index_map[bl_sub.size()] << std::endl;
                bl_sub.push_back(grids[lev][j]);
                dm_sub.push_back(dmap[lev][j]);
            } // intersects

        } // loop over ba (j)

        BoxArray ba_sub(bl_sub);

        // Define MultiFabs that hold only the data in this particular subdomain
        rhs_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, rhs[0].nGrowVect(), MFInfo{}.SetAlloc(false));
        phi_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, phi[0].nGrowVect(), MFInfo{}.SetAlloc(false));

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            fluxes_sub[0][idim].define(convert(ba_sub, IntVect::TheDimensionVector(idim)), DistributionMapping(dm_sub), 1,
                                               IntVect::TheZeroVector(), MFInfo{}.SetAlloc(false));
        }

        // Link the new MultiFabs to the FABs in the original MultiFabs (no copy required)
        for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi) {
            int orig_index = index_map[mfi.index()];
            // amrex::Print() << " INDEX        " << orig_index << " TO " << mfi.index() << std::endl;
            rhs_sub[0].setFab(mfi, FArrayBox(rhs[0][orig_index], amrex::make_alias, 0, 1));
            phi_sub[0].setFab(mfi, FArrayBox(phi[0][orig_index], amrex::make_alias, 0, 1));
            fluxes_sub[0][0].setFab(mfi,FArrayBox(fluxes[0][0][orig_index], amrex::make_alias, 0, 1));
            fluxes_sub[0][1].setFab(mfi,FArrayBox(fluxes[0][1][orig_index], amrex::make_alias, 0, 1));
            fluxes_sub[0][2].setFab(mfi,FArrayBox(fluxes[0][2][orig_index], amrex::make_alias, 0, 1));
        }

        if (solverChoice.mesh_type == MeshType::VariableDz) {
            ax_sub.define(convert(ba_sub,IntVect(1,0,0)), DistributionMapping(dm_sub), 1,
                          ax[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            ay_sub.define(convert(ba_sub,IntVect(0,1,0)), DistributionMapping(dm_sub), 1,
                          ay[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            az_sub.define(convert(ba_sub,IntVect(0,0,1)), DistributionMapping(dm_sub), 1,
                          az[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            znd_sub.define(convert(ba_sub,IntVect(1,1,1)), DistributionMapping(dm_sub), 1,
                           z_phys_nd[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
             dJ_sub.define(ba_sub, DistributionMapping(dm_sub), 1,
                           detJ_cc[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));

            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi) {
                int orig_index = index_map[mfi.index()];
                ax_sub.setFab(mfi, FArrayBox((*ax[lev])[orig_index], amrex::make_alias, 0, 1));
                ay_sub.setFab(mfi, FArrayBox((*ay[lev])[orig_index], amrex::make_alias, 0, 1));
                az_sub.setFab(mfi, FArrayBox((*az[lev])[orig_index], amrex::make_alias, 0, 1));
                znd_sub.setFab(mfi, FArrayBox((*z_phys_nd[lev])[orig_index], amrex::make_alias, 0, 1));
                 dJ_sub.setFab(mfi, FArrayBox((*detJ_cc[lev])[orig_index], amrex::make_alias, 0, 1));
            }

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
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi); const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
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

        if (lev > 0) {
           amrex::Print() << "RHSSUB BA " << rhs_sub[0].boxArray() << std::endl;
        }

        // ****************************************************************************
        // EB
        // ****************************************************************************
        if (solverChoice.terrain_type == TerrainType::EB) {
            solve_with_EB_mlmg(lev, rhs_sub, phi_sub, fluxes_sub);
        } else {

        // ****************************************************************************
        // No terrain or grid stretching
        // ****************************************************************************
        if (solverChoice.mesh_type == MeshType::ConstantDz) {
#ifdef ERF_USE_FFT
            if (use_fft) {
                Box my_region(subdomains[lev][isub].minimalBox());
                bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
                if (boxes_make_rectangle) {
                    solve_with_fft(lev, my_region, rhs_sub[0], phi_sub[0], fluxes_sub[0]);
                } else {
                    amrex::Warning("FFT won't work unless the union of boxes is rectangular: defaulting to MLMG");
                    solve_with_mlmg(lev, rhs_sub, phi_sub, fluxes_sub);
                }
        } else {
            solve_with_mlmg(lev, rhs, phi, fluxes);
        }
#else
        if (use_fft) {
            amrex::Warning("You set use_fft=true but didn't build with USE_FFT = TRUE; defaulting to MLMG");
        }
        solve_with_mlmg(lev, rhs_sub, phi_sub, fluxes_sub);
#endif
    } // No terrain or grid stretching

    // ****************************************************************************
    // Grid stretching (flat terrain)
    // ****************************************************************************
    else if (solverChoice.mesh_type == MeshType::StretchedDz) {
#ifndef ERF_USE_FFT
        amrex::Abort("Rebuild with USE_FFT = TRUE so you can use the FFT solver");
#else
        Box my_region(subdomains[lev][isub].minimalBox());
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
            Box my_region(subdomains[lev][isub].minimalBox());
            bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
            if (!boxes_make_rectangle) {
                amrex::Abort("FFT preconditioner for GMRES won't work unless the union of boxes is rectangular");
            } else {
                solve_with_gmres(lev, my_region, rhs_sub[0], phi_sub[0], fluxes_sub[0], ax_sub, ay_sub, dJ_sub, znd_sub);
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

        } // not EB

    } // loop over subdomains (i)

    // ****************************************************************************
    // Print time in solve
    // ****************************************************************************
    Real end_step = static_cast<Real>(ParallelDescriptor::second());
    if (mg_verbose > 0) {
        amrex::Print() << "Time in solve " << end_step - start_step << std::endl;
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
        compute_divergence(lev, rhs_lev, rho0_u_const, geom_tmp[0]);

        sum = volWgtSumMF(lev,rhs_lev,0,false);

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
                    amrex::AllPrint() << "RHS AFTER SOLVE AT " <<
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
    if (!solverChoice.fixed_density) {
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
    MultiFab::Saxpy(pp_inc[lev], 1.0/l_dt, phi[0],0,0,1,1);
}
