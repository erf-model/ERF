#include <EOS.H>
#include <ERF.H>
#include "AMReX_Interp_3D_C.H"
#include "AMReX_PlotFileUtil.H"
#include "TerrainMetrics.H"
#include "ERF_Constants.H"

using namespace amrex;

template<typename V, typename T>
bool containerHasElement (const V& iterable, const T& query) {
    return std::find(iterable.begin(), iterable.end(), query) != iterable.end();
}

void
ERF::setPlotVariables (const std::string& pp_plot_var_names, Vector<std::string>& plot_var_names)
{
    ParmParse pp(pp_prefix);

    if (pp.contains(pp_plot_var_names.c_str()))
    {
        std::string nm;

        int nPltVars = pp.countval(pp_plot_var_names.c_str());

        for (int i = 0; i < nPltVars; i++)
        {
            pp.get(pp_plot_var_names.c_str(), nm, i);

            // Add the named variable to our list of plot variables
            // if it is not already in the list
            if (!containerHasElement(plot_var_names, nm)) {
                plot_var_names.push_back(nm);
            }
        }
    } else {
        //
        // The default is to add none of the variables to the list
        //
        plot_var_names.clear();
    }

    // Get state variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_plot_names;

    int n_qstate   = micro.Get_Qstate_Size();
    int ncomp_cons = NVAR_max - (NMOIST_max - n_qstate);
    for (int i = 0; i < ncomp_cons; ++i) {
        if ( containerHasElement(plot_var_names, cons_names[i]) ) {
            tmp_plot_names.push_back(cons_names[i]);
        }
    }
    // check for velocity since it's not in cons_names
    // if we are asked for any velocity component, we will need them all
    if (containerHasElement(plot_var_names, "x_velocity") ||
        containerHasElement(plot_var_names, "y_velocity") ||
        containerHasElement(plot_var_names, "z_velocity")) {
        tmp_plot_names.push_back("x_velocity");
        tmp_plot_names.push_back("y_velocity");
        tmp_plot_names.push_back("z_velocity");
    }
    for (int i = 0; i < derived_names.size(); ++i) {
        if ( containerHasElement(plot_var_names, derived_names[i]) ) {
            if (solverChoice.use_terrain || (derived_names[i] != "z_phys" && derived_names[i] != "detJ") ) {
               tmp_plot_names.push_back(derived_names[i]);
            }
        }
    }
#ifdef ERF_USE_PARTICLES
    const auto& particles_namelist( particleData.getNamesUnalloc() );
    for (auto it = particles_namelist.cbegin(); it != particles_namelist.cend(); ++it) {
        std::string tmp( (*it)+"_count" );
        if (containerHasElement(plot_var_names, tmp) ) {
            tmp_plot_names.push_back(tmp);
        }
    }
#endif

    plot_var_names = tmp_plot_names;
}

void
ERF::appendPlotVariables (const std::string& pp_plot_var_names, Vector<std::string>& a_plot_var_names)
{
    ParmParse pp(pp_prefix);

    Vector<std::string> plot_var_names(0);
    if (pp.contains(pp_plot_var_names.c_str())) {
        std::string nm;
        int nPltVars = pp.countval(pp_plot_var_names.c_str());
        for (int i = 0; i < nPltVars; i++) {
            pp.get(pp_plot_var_names.c_str(), nm, i);
            // Add the named variable to our list of plot variables
            // if it is not already in the list
            if (!containerHasElement(plot_var_names, nm)) {
                plot_var_names.push_back(nm);
            }
        }
    }

    Vector<std::string> tmp_plot_names(0);
#ifdef ERF_USE_PARTICLES
    Vector<std::string> particle_mesh_plot_names;
    particleData.GetMeshPlotVarNames( particle_mesh_plot_names );
    for (int i = 0; i < particle_mesh_plot_names.size(); i++) {
        std::string tmp(particle_mesh_plot_names[i]);
        if (containerHasElement(plot_var_names, tmp) ) {
            tmp_plot_names.push_back(tmp);
        }
    }
#endif

    for (int i = 0; i < tmp_plot_names.size(); i++) {
        a_plot_var_names.push_back( tmp_plot_names[i] );
    }

    // Finally, check to see if we found all the requested variables
    for (const auto& plot_name : plot_var_names) {
        if (!containerHasElement(a_plot_var_names, plot_name)) {
             if (amrex::ParallelDescriptor::IOProcessor()) {
                 Warning("\nWARNING: Requested to plot variable '" + plot_name + "' but it is not available");
             }
        }
    }
}

// set plotfile variable names
Vector<std::string>
ERF::PlotFileVarNames (Vector<std::string> plot_var_names )
{
    Vector<std::string> names;

    names.insert(names.end(), plot_var_names.begin(), plot_var_names.end());

    return names;

}

// Write plotfile to disk
void
ERF::WritePlotFile (int which, Vector<std::string> plot_var_names)
{
    const Vector<std::string> varnames = PlotFileVarNames(plot_var_names);
    const int ncomp_mf = varnames.size();

    int ncomp_cons = vars_new[0][Vars::cons].nComp();

    if (ncomp_mf == 0) return;

    // We Fillpatch here because some of the derived quantities require derivatives
    //     which require ghost cells to be filled.  We do not need to call FillPatcher
    //     because we don't need to set interior fine points.
    for (int lev = 0; lev <= finest_level; ++lev) {
        bool fillset = false;
        FillPatch(lev, t_new[lev], {&vars_new[lev][Vars::cons], &vars_new[lev][Vars::xvel],
                                    &vars_new[lev][Vars::yvel], &vars_new[lev][Vars::zvel]},
                                   {&vars_new[lev][Vars::cons], &rU_new[lev],
                                    &rV_new[lev], &rW_new[lev]}, fillset);
    }

    // Get qmoist pointers if using moisture
    bool use_moisture = (solverChoice.moisture_type != MoistureType::None);
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int mvar(0); mvar<qmoist[lev].size(); ++mvar) {
            qmoist[lev][mvar] = micro.Get_Qmoist_Ptr(lev,mvar);
        }
    }

    Vector<MultiFab> mf(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp_mf, 0);
    }

    Vector<MultiFab> mf_nd(finest_level+1);
    if (solverChoice.use_terrain) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            BoxArray nodal_grids(grids[lev]); nodal_grids.surroundingNodes();
            mf_nd[lev].define(nodal_grids, dmap[lev], 3, 0);
            mf_nd[lev].setVal(0.);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        int mf_comp = 0;

        // First, copy any of the conserved state variables into the output plotfile
        AMREX_ALWAYS_ASSERT(cons_names.size() >= ncomp_cons);
        for (int i = 0; i < ncomp_cons; ++i) {
            if (containerHasElement(plot_var_names, cons_names[i])) {
                MultiFab::Copy(mf[lev],vars_new[lev][Vars::cons],i,mf_comp,1,0);
                mf_comp++;
            }
        }

        // Next, check for velocities and if desired, output them -- note we output none or all, not just some
        if (containerHasElement(plot_var_names, "x_velocity") ||
            containerHasElement(plot_var_names, "y_velocity") ||
            containerHasElement(plot_var_names, "z_velocity")) {

            average_face_to_cellcenter(mf[lev],mf_comp,
                Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],&vars_new[lev][Vars::yvel],&vars_new[lev][Vars::zvel]});
            mf_comp += AMREX_SPACEDIM;
        }

        // Finally, check for any derived quantities and compute them, inserting
        // them into our output multifab
        auto calculate_derived = [&](const std::string& der_name,
                                     decltype(derived::erf_dernull)& der_function)
        {
            if (containerHasElement(plot_var_names, der_name)) {
                MultiFab dmf(mf[lev], make_alias, mf_comp, 1);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    auto& dfab = dmf[mfi];
                    auto& sfab = vars_new[lev][Vars::cons][mfi];
                    der_function(bx, dfab, 0, 1, sfab, Geom(lev), t_new[0], nullptr, lev);
                }

                mf_comp++;
            }
        };

        // Note: All derived variables must be computed in order of "derived_names" defined in ERF.H
        calculate_derived("soundspeed",  derived::erf_dersoundspeed);
        calculate_derived("temp",        derived::erf_dertemp);
        calculate_derived("theta",       derived::erf_dertheta);
        calculate_derived("KE",          derived::erf_derKE);
        calculate_derived("QKE",         derived::erf_derQKE);
        calculate_derived("scalar",      derived::erf_derscalar);

        MultiFab r_hse(base_state[lev], make_alias, 0, 1); // r_0 is first  component
        MultiFab p_hse(base_state[lev], make_alias, 1, 1); // p_0 is second component
        if (containerHasElement(plot_var_names, "pres_hse"))
        {
            // p_0 is second component of base_state
            MultiFab::Copy(mf[lev],p_hse,0,mf_comp,1,0);
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names, "dens_hse"))
        {
            // r_0 is first component of base_state
            MultiFab::Copy(mf[lev],r_hse,0,mf_comp,1,0);
            mf_comp += 1;
        }

        if (containerHasElement(plot_var_names, "pressure"))
        {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real      >& derdat = mf[lev].array(mfi);
                const Array4<Real const>&  S_arr = vars_new[lev][Vars::cons].const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    Real qv_for_p = (use_moisture) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0;
                    const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                    derdat(i, j, k, mf_comp) = getPgivenRTh(rhotheta,qv_for_p);
                });
            }
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names, "pert_pres"))
        {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<Real const>& p0_arr = p_hse.const_array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    Real qv_for_p = (use_moisture) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0;
                    const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                    derdat(i, j, k, mf_comp) = getPgivenRTh(rhotheta,qv_for_p) - p0_arr(i,j,k);
                });
            }
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names, "pert_dens"))
        {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat  = mf[lev].array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                const Array4<Real const>& r0_arr = r_hse.const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i, j, k, mf_comp) = S_arr(i,j,k,Rho_comp) - r0_arr(i,j,k);
                });
            }
            mf_comp ++;
        }

        if (containerHasElement(plot_var_names, "eq_pot_temp"))
        {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat  = mf[lev].array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    Real qv = (use_moisture) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0.0;
                    Real qc = (use_moisture) ? S_arr(i,j,k,RhoQ2_comp)/S_arr(i,j,k,Rho_comp) : 0.0;
                    Real T = getTgivenRandRTh(S_arr(i,j,k,Rho_comp), S_arr(i,j,k,RhoTheta_comp), qv);
                    Real pressure = getPgivenRTh(S_arr(i,j,k,RhoTheta_comp), qv);
                    Real fac = Cp_d + Cp_l*(qv + qc);
                    Real pv = erf_esatw(T)*100.0;

                    derdat(i, j, k, mf_comp) = T*std::pow((pressure - pv)/p_0, -R_d/fac)*std::exp(L_v*qv/(fac*T)) ;
                });
            }
            mf_comp ++;
        }

        int klo = geom[lev].Domain().smallEnd(2);
        int khi = geom[lev].Domain().bigEnd(2);

        if (containerHasElement(plot_var_names, "dpdx"))
        {
            auto dxInv = geom[lev].InvCellSizeArray();
            MultiFab pres(vars_new[lev][Vars::cons].boxArray(), vars_new[lev][Vars::cons].DistributionMap(), 1, 1);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // First define pressure on grown box
                const Box& gbx = mfi.growntilebox(1);
                const Array4<Real> & p_arr  = pres.array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    p_arr(i,j,k) = getPgivenRTh(S_arr(i,j,k,RhoTheta_comp));
                });
            }
            pres.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Now compute pressure gradient on valid box
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<Real> & p_arr  = pres.array(mfi);

                if (solverChoice.use_terrain) {
                    const Array4<Real const>& z_nd = z_phys_nd[lev]->const_array(mfi);

                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                        // Pgrad at lower I face
                        Real met_h_xi_lo   = Compute_h_xi_AtIface  (i, j, k, dxInv, z_nd);
                        Real met_h_zeta_lo = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);
                        Real gp_xi_lo = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
                        Real gp_zeta_on_iface_lo;
                        if(k == klo) {
                            gp_zeta_on_iface_lo = 0.5 * dxInv[2] * (
                                p_arr(i-1,j,k+1) + p_arr(i,j,k+1)
                                - p_arr(i-1,j,k  ) - p_arr(i,j,k  ) );
                        } else if (k == khi) {
                            gp_zeta_on_iface_lo = 0.5 * dxInv[2] * (
                                p_arr(i-1,j,k  ) + p_arr(i,j,k  )
                                - p_arr(i-1,j,k-1) - p_arr(i,j,k-1) );
                        } else {
                            gp_zeta_on_iface_lo = 0.25 * dxInv[2] * (
                                p_arr(i-1,j,k+1) + p_arr(i,j,k+1)
                                - p_arr(i-1,j,k-1) - p_arr(i,j,k-1) );
                        }
                        amrex::Real gpx_lo = gp_xi_lo - (met_h_xi_lo/ met_h_zeta_lo) * gp_zeta_on_iface_lo;

                        // Pgrad at higher I face
                        Real met_h_xi_hi   = Compute_h_xi_AtIface  (i+1, j, k, dxInv, z_nd);
                        Real met_h_zeta_hi = Compute_h_zeta_AtIface(i+1, j, k, dxInv, z_nd);
                        Real gp_xi_hi = dxInv[0] * (p_arr(i+1,j,k) - p_arr(i,j,k));
                        Real gp_zeta_on_iface_hi;
                        if(k == klo) {
                            gp_zeta_on_iface_hi = 0.5 * dxInv[2] * (
                                p_arr(i+1,j,k+1) + p_arr(i,j,k+1)
                                - p_arr(i+1,j,k  ) - p_arr(i,j,k  ) );
                        } else if (k == khi) {
                            gp_zeta_on_iface_hi = 0.5 * dxInv[2] * (
                                p_arr(i+1,j,k  ) + p_arr(i,j,k  )
                                - p_arr(i+1,j,k-1) - p_arr(i,j,k-1) );
                        } else {
                            gp_zeta_on_iface_hi = 0.25 * dxInv[2] * (
                                p_arr(i+1,j,k+1) + p_arr(i,j,k+1)
                                - p_arr(i+1,j,k-1) - p_arr(i,j,k-1) );
                        }
                        amrex::Real gpx_hi = gp_xi_hi - (met_h_xi_hi/ met_h_zeta_hi) * gp_zeta_on_iface_hi;

                        // Average P grad to CC
                        derdat(i ,j ,k, mf_comp) = 0.5 * (gpx_lo + gpx_hi);
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        derdat(i ,j ,k, mf_comp) = 0.5 * (p_arr(i+1,j,k) - p_arr(i-1,j,k)) * dxInv[0];
                    });
                }
            } // mfi
            mf_comp ++;
        } // dpdx

        if (containerHasElement(plot_var_names, "dpdy"))
        {
            auto dxInv = geom[lev].InvCellSizeArray();

            MultiFab pres(vars_new[lev][Vars::cons].boxArray(), vars_new[lev][Vars::cons].DistributionMap(), 1, 1);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // First define pressure on grown box
                const Box& gbx = mfi.growntilebox(1);
                const Array4<Real> & p_arr  = pres.array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    p_arr(i,j,k) = getPgivenRTh(S_arr(i,j,k,RhoTheta_comp));
                });
            }
            pres.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Now compute pressure gradient on valid box
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<Real> & p_arr  = pres.array(mfi);

                if (solverChoice.use_terrain) {
                    const Array4<Real const>& z_nd = z_phys_nd[lev]->const_array(mfi);

                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                        Real met_h_eta_lo  = Compute_h_eta_AtJface (i, j, k, dxInv, z_nd);
                        Real met_h_zeta_lo = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);
                        Real gp_eta_lo = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
                        Real gp_zeta_on_jface_lo;
                        if (k == klo) {
                            gp_zeta_on_jface_lo = 0.5 * dxInv[2] * (
                                p_arr(i,j,k+1) + p_arr(i,j-1,k+1)
                                - p_arr(i,j,k  ) - p_arr(i,j-1,k  ) );
                        } else if (k == khi) {
                            gp_zeta_on_jface_lo = 0.5 * dxInv[2] * (
                                p_arr(i,j,k  ) + p_arr(i,j-1,k  )
                                - p_arr(i,j,k-1) - p_arr(i,j-1,k-1) );
                        } else {
                            gp_zeta_on_jface_lo = 0.25 * dxInv[2] * (
                                p_arr(i,j,k+1) + p_arr(i,j-1,k+1)
                                - p_arr(i,j,k-1) - p_arr(i,j-1,k-1) );
                        }
                        amrex::Real gpy_lo = gp_eta_lo - (met_h_eta_lo / met_h_zeta_lo) * gp_zeta_on_jface_lo;

                        Real met_h_eta_hi  = Compute_h_eta_AtJface (i, j+1, k, dxInv, z_nd);
                        Real met_h_zeta_hi = Compute_h_zeta_AtJface(i, j+1, k, dxInv, z_nd);
                        Real gp_eta_hi = dxInv[1] * (p_arr(i,j+1,k) - p_arr(i,j,k));
                        Real gp_zeta_on_jface_hi;
                        if (k == klo) {
                            gp_zeta_on_jface_hi = 0.5 * dxInv[2] * (
                                p_arr(i,j+1,k+1) + p_arr(i,j,k+1)
                                - p_arr(i,j+1,k  ) - p_arr(i,j,k  ) );
                        } else if (k == khi) {
                            gp_zeta_on_jface_hi = 0.5 * dxInv[2] * (
                                p_arr(i,j+1,k  ) + p_arr(i,j,k  )
                                - p_arr(i,j+1,k-1) - p_arr(i,j,k-1) );
                        } else {
                            gp_zeta_on_jface_hi = 0.25 * dxInv[2] * (
                                p_arr(i,j+1,k+1) + p_arr(i,j,k+1)
                                - p_arr(i,j+1,k-1) - p_arr(i,j,k-1) );
                        }
                        amrex::Real gpy_hi = gp_eta_hi - (met_h_eta_hi / met_h_zeta_hi) * gp_zeta_on_jface_hi;

                        derdat(i ,j ,k, mf_comp) = 0.5 * (gpy_lo + gpy_hi);
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        derdat(i ,j ,k, mf_comp) = 0.5 * (p_arr(i,j+1,k) - p_arr(i,j-1,k)) * dxInv[1];
                    });
                }
            } // mf
            mf_comp ++;
        } // dpdy

        if (containerHasElement(plot_var_names, "pres_hse_x"))
        {
            auto dxInv = geom[lev].InvCellSizeArray();
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real      >&  derdat = mf[lev].array(mfi);
                const Array4<Real const>&   p_arr = p_hse.const_array(mfi);

                //USE_TERRAIN POSSIBLE ISSUE HERE
                const Array4<Real const>& z_nd  = z_phys_nd[lev]->const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    Real met_h_xi_lo   = Compute_h_xi_AtIface  (i, j, k, dxInv, z_nd);
                    Real met_h_zeta_lo = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);
                    Real gp_xi_lo = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
                    Real gp_zeta_on_iface_lo;
                    if (k == klo) {
                      gp_zeta_on_iface_lo = 0.5 * dxInv[2] * (
                                                              p_arr(i-1,j,k+1) + p_arr(i,j,k+1)
                                                            - p_arr(i-1,j,k  ) - p_arr(i,j,k  ) );
                    } else if (k == khi) {
                      gp_zeta_on_iface_lo = 0.5 * dxInv[2] * (
                                                              p_arr(i-1,j,k  ) + p_arr(i,j,k  )
                                                            - p_arr(i-1,j,k-1) - p_arr(i,j,k-1) );
                    } else {
                      gp_zeta_on_iface_lo = 0.25 * dxInv[2] * (
                                                               p_arr(i-1,j,k+1) + p_arr(i,j,k+1)
                                                             - p_arr(i-1,j,k-1) - p_arr(i,j,k-1) );
                    }
                    amrex::Real gpx_lo = gp_xi_lo - (met_h_xi_lo/ met_h_zeta_lo) * gp_zeta_on_iface_lo;

                    Real met_h_xi_hi   = Compute_h_xi_AtIface  (i+1, j, k, dxInv, z_nd);
                    Real met_h_zeta_hi = Compute_h_zeta_AtIface(i+1, j, k, dxInv, z_nd);
                    Real gp_xi_hi = dxInv[0] * (p_arr(i+1,j,k) - p_arr(i,j,k));
                    Real gp_zeta_on_iface_hi;
                    if (k == klo) {
                      gp_zeta_on_iface_hi = 0.5 * dxInv[2] * (
                                                              p_arr(i+1,j,k+1) + p_arr(i,j,k+1)
                                                            - p_arr(i+1,j,k  ) - p_arr(i,j,k  ) );
                    } else if (k == khi) {
                      gp_zeta_on_iface_hi = 0.5 * dxInv[2] * (
                                                              p_arr(i+1,j,k  ) + p_arr(i,j,k  )
                                                            - p_arr(i+1,j,k-1) - p_arr(i,j,k-1) );
                    } else {
                      gp_zeta_on_iface_hi = 0.25 * dxInv[2] * (
                                                               p_arr(i+1,j,k+1) + p_arr(i,j,k+1)
                                                             - p_arr(i+1,j,k-1) - p_arr(i,j,k-1) );
                    }
                    amrex::Real gpx_hi = gp_xi_hi - (met_h_xi_hi/ met_h_zeta_hi) * gp_zeta_on_iface_hi;

                    derdat(i ,j ,k, mf_comp) = 0.5 * (gpx_lo + gpx_hi);
                });
            }
            mf_comp += 1;
        } // pres_hse_x

        if (containerHasElement(plot_var_names, "pres_hse_y"))
        {
            auto dxInv = geom[lev].InvCellSizeArray();
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real      >& derdat = mf[lev].array(mfi);
                const Array4<Real const>&   p_arr = p_hse.const_array(mfi);
                const Array4<Real const>& z_nd    = z_phys_nd[lev]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    Real met_h_eta_lo  = Compute_h_eta_AtJface (i, j, k, dxInv, z_nd);
                    Real met_h_zeta_lo = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);
                    Real gp_eta_lo = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
                    Real gp_zeta_on_jface_lo;
                    if (k == klo) {
                      gp_zeta_on_jface_lo = 0.5 * dxInv[2] * (
                                                              p_arr(i,j,k+1) + p_arr(i,j-1,k+1)
                                                            - p_arr(i,j,k  ) - p_arr(i,j-1,k  ) );
                    } else if (k == khi) {
                      gp_zeta_on_jface_lo = 0.5 * dxInv[2] * (
                                                              p_arr(i,j,k  ) + p_arr(i,j-1,k  )
                                                            - p_arr(i,j,k-1) - p_arr(i,j-1,k-1) );
                    } else {
                      gp_zeta_on_jface_lo = 0.25 * dxInv[2] * (
                                                               p_arr(i,j,k+1) + p_arr(i,j-1,k+1)
                                                             - p_arr(i,j,k-1) - p_arr(i,j-1,k-1) );
                    }
                    amrex::Real gpy_lo = gp_eta_lo - (met_h_eta_lo / met_h_zeta_lo) * gp_zeta_on_jface_lo;

                    Real met_h_eta_hi  = Compute_h_eta_AtJface (i, j+1, k, dxInv, z_nd);
                    Real met_h_zeta_hi = Compute_h_zeta_AtJface(i, j+1, k, dxInv, z_nd);
                    Real gp_eta_hi = dxInv[1] * (p_arr(i,j+1,k) - p_arr(i,j,k));
                    Real gp_zeta_on_jface_hi;
                    if (k == klo) {
                      gp_zeta_on_jface_hi = 0.5 * dxInv[2] * (
                                                              p_arr(i,j+1,k+1) + p_arr(i,j,k+1)
                                                            - p_arr(i,j+1,k  ) - p_arr(i,j,k  ) );
                    } else if (k == khi) {
                      gp_zeta_on_jface_hi = 0.5 * dxInv[2] * (
                                                              p_arr(i,j+1,k  ) + p_arr(i,j,k  )
                                                            - p_arr(i,j+1,k-1) - p_arr(i,j,k-1) );
                    } else {
                      gp_zeta_on_jface_hi = 0.25 * dxInv[2] * (
                                                               p_arr(i,j+1,k+1) + p_arr(i,j,k+1)
                                                             - p_arr(i,j+1,k-1) - p_arr(i,j,k-1) );
                    }
                    amrex::Real gpy_hi = gp_eta_hi - (met_h_eta_hi / met_h_zeta_hi) * gp_zeta_on_jface_hi;

                    derdat(i ,j ,k, mf_comp) = 0.5 * (gpy_lo + gpy_hi);
                });
            }
            mf_comp += 1;
        } // pres_hse_y

        if (solverChoice.use_terrain) {
            if (containerHasElement(plot_var_names, "z_phys"))
            {
                MultiFab::Copy(mf[lev],*z_phys_cc[lev],0,mf_comp,1,0);
                mf_comp ++;
            }

            if (containerHasElement(plot_var_names, "detJ"))
            {
                MultiFab::Copy(mf[lev],*detJ_cc[lev],0,mf_comp,1,0);
                mf_comp ++;
            }
        } // use_terrain

        if (containerHasElement(plot_var_names, "mapfac")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<Real>& mf_m   = mapfac_m[lev]->array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   derdat(i ,j ,k, mf_comp) = mf_m(i,j,0);
                });
            }
            mf_comp ++;
        }

        if (containerHasElement(plot_var_names, "Kmv")) {
            MultiFab::Copy(mf[lev],*eddyDiffs_lev[lev],EddyDiff::Mom_v,mf_comp,1,0);
            mf_comp ++;
        }
        if (containerHasElement(plot_var_names, "Kmh")) {
            MultiFab::Copy(mf[lev],*eddyDiffs_lev[lev],EddyDiff::Mom_h,mf_comp,1,0);
            mf_comp ++;
        }
        if (containerHasElement(plot_var_names, "Khv")) {
            MultiFab::Copy(mf[lev],*eddyDiffs_lev[lev],EddyDiff::Theta_v,mf_comp,1,0);
            mf_comp ++;
        }
        if (containerHasElement(plot_var_names, "Khh")) {
            MultiFab::Copy(mf[lev],*eddyDiffs_lev[lev],EddyDiff::Theta_h,mf_comp,1,0);
            mf_comp ++;
        }
        if (containerHasElement(plot_var_names, "Lpbl")) {
            MultiFab::Copy(mf[lev],*eddyDiffs_lev[lev],EddyDiff::PBL_lengthscale,mf_comp,1,0);
            mf_comp ++;
        }

        // TODO: The size of the q variables can vary with different
        //       moisture models. Therefore, certain components may
        //       reside at different indices. For example, Kessler is
        //       warm but precipitating. This puts qp at index 3.
        //       However, SAM is cold and precipitating so qp is index 4.
        //       Need to built an external enum struct or a better pathway.

        // NOTE: Protect against accessing non-existent data
        if (use_moisture) {
            int q_size = qmoist[lev].size();

            if (containerHasElement(plot_var_names, "qt") && (q_size >= 1))
            {
                MultiFab qt_mf(*(qmoist[lev][0]), make_alias, 0, 1);
                MultiFab::Copy(mf[lev],qt_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }

            if (containerHasElement(plot_var_names, "qv") && (q_size >= 2))
            {
                MultiFab qv_mf(*(qmoist[lev][1]), make_alias, 0, 1);
                MultiFab::Copy(mf[lev],qv_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }

            if (containerHasElement(plot_var_names, "qc") && (q_size >= 3))
            {
                MultiFab qc_mf(*(qmoist[lev][2]), make_alias, 0, 1);
                MultiFab::Copy(mf[lev],qc_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }

            if (containerHasElement(plot_var_names, "qi") && (q_size > 4))
            {
                MultiFab qi_mf(*(qmoist[lev][3]), make_alias, 0, 1);
                MultiFab::Copy(mf[lev],qi_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }

            if (containerHasElement(plot_var_names, "qp") &&
               ((q_size >= 5) || (q_size==4)))
            {
                int q_ind;
                if (q_size >= 5) { // Cold moisture physics (has ice)
                    q_ind = 4;
                } else {           // Warm moisture physics (no  ice)
                    q_ind = 3;
                }
                MultiFab qp_mf(*(qmoist[lev][q_ind]), make_alias, 0, 1);
                MultiFab::Copy(mf[lev],qp_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }

            if (containerHasElement(plot_var_names, "qrain") && (q_size >= 6))
            {
                MultiFab qr_mf(*(qmoist[lev][5]), make_alias, 0, 1);
                MultiFab::Copy(mf[lev],qr_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }

            if (containerHasElement(plot_var_names, "qsnow") && (q_size >= 7))
            {
                MultiFab qs_mf(*(qmoist[lev][6]), make_alias, 0, 1);
                MultiFab::Copy(mf[lev],qs_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }

            if (containerHasElement(plot_var_names, "qgraup") && (q_size >= 8))
            {
                MultiFab qg_mf(*(qmoist[lev][7]), make_alias, 0, 1);
                MultiFab::Copy(mf[lev],qg_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }
        }

#ifdef ERF_USE_PARTICLES
        const auto& particles_namelist( particleData.getNames() );
        for (ParticlesNamesVector::size_type i = 0; i < particles_namelist.size(); i++) {
            if (containerHasElement(plot_var_names, std::string(particles_namelist[i]+"_count"))) {
                MultiFab temp_dat(mf[lev].boxArray(), mf[lev].DistributionMap(), 1, 0);
                temp_dat.setVal(0);
                particleData[particles_namelist[i]]->Increment(temp_dat, lev);
                MultiFab::Copy(mf[lev], temp_dat, 0, mf_comp, 1, 0);
                mf_comp += 1;
            }
        }

        Vector<std::string> particle_mesh_plot_names(0);
        particleData.GetMeshPlotVarNames( particle_mesh_plot_names );
        for (int i = 0; i < particle_mesh_plot_names.size(); i++) {
            std::string plot_var_name(particle_mesh_plot_names[i]);
            if (containerHasElement(plot_var_names, plot_var_name) ) {
                MultiFab temp_dat(mf[lev].boxArray(), mf[lev].DistributionMap(), 1, 1);
                temp_dat.setVal(0);
                particleData.GetMeshPlotVar(plot_var_name, temp_dat, lev);
                MultiFab::Copy(mf[lev], temp_dat, 0, mf_comp, 1, 0);
                mf_comp += 1;
            }
        }
#endif

#ifdef ERF_COMPUTE_ERROR
        // Next, check for error in velocities and if desired, output them -- note we output none or all, not just some
        if (containerHasElement(plot_var_names, "xvel_err") ||
            containerHasElement(plot_var_names, "yvel_err") ||
            containerHasElement(plot_var_names, "zvel_err"))
        {
            //
            // Moving terrain ANALYTICAL
            //
            Real H           = geom[lev].ProbHi()[2];
            Real Ampl        = 0.16;
            Real wavelength  = 100.;
            Real kp          = 2. * PI / wavelength;
            Real g           = CONST_GRAV;
            Real omega       = std::sqrt(g * kp);
            Real omega_t     = omega * t_new[lev];

            const auto dx = geom[lev].CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Box xbx(bx); xbx.surroundingNodes(0);
                const Array4<Real> xvel_arr = vars_new[lev][Vars::xvel].array(mfi);
                const Array4<Real> zvel_arr = vars_new[lev][Vars::zvel].array(mfi);

                const Array4<Real const>& z_nd = z_phys_nd[lev]->const_array(mfi);

                ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real x = i * dx[0];
                    Real z = 0.25 * (z_nd(i,j,k) + z_nd(i,j+1,k) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1));

                    Real z_base = Ampl * std::sin(kp * x - omega_t);
                    z -= z_base;

                    Real fac = std::cosh( kp * (z - H) ) / std::sinh(kp * H);

                    xvel_arr(i,j,k) -= -Ampl * omega * fac * std::sin(kp * x - omega_t);
                });

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real x   = (i + 0.5) * dx[0];
                    Real z   = 0.25 * ( z_nd(i,j,k) + z_nd(i+1,j,k) + z_nd(i,j+1,k) + z_nd(i+1,j+1,k));

                    Real z_base = Ampl * std::sin(kp * x - omega_t);
                    z -= z_base;

                    Real fac = std::sinh( kp * (z - H) ) / std::sinh(kp * H);

                    zvel_arr(i,j,k) -= Ampl * omega * fac * std::cos(kp * x - omega_t);
                });
            }

            MultiFab temp_mf(mf[lev].boxArray(), mf[lev].DistributionMap(), AMREX_SPACEDIM, 0);
            average_face_to_cellcenter(temp_mf,0,
                Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],&vars_new[lev][Vars::yvel],&vars_new[lev][Vars::zvel]});

            if (containerHasElement(plot_var_names, "xvel_err")) {
                MultiFab::Copy(mf[lev],temp_mf,0,mf_comp,1,0);
                mf_comp += 1;
            }
            if (containerHasElement(plot_var_names, "yvel_err")) {
                MultiFab::Copy(mf[lev],temp_mf,1,mf_comp,1,0);
                mf_comp += 1;
            }
            if (containerHasElement(plot_var_names, "zvel_err")) {
                MultiFab::Copy(mf[lev],temp_mf,2,mf_comp,1,0);
                mf_comp += 1;
            }

            // Now restore the velocities to what they were
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Box xbx(bx); xbx.surroundingNodes(0);

                const Array4<Real> xvel_arr = vars_new[lev][Vars::xvel].array(mfi);
                const Array4<Real> zvel_arr = vars_new[lev][Vars::zvel].array(mfi);

                const Array4<Real const>& z_nd = z_phys_nd[lev]->const_array(mfi);

                ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real x = i * dx[0];
                    Real z = 0.25 * (z_nd(i,j,k) + z_nd(i,j+1,k) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1));
                    Real z_base = Ampl * std::sin(kp * x - omega_t);

                    z -= z_base;

                    Real fac = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
                    xvel_arr(i,j,k) += -Ampl * omega * fac * std::sin(kp * x - omega_t);
                });
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real x   = (i + 0.5) * dx[0];
                    Real z   = 0.25 * ( z_nd(i,j,k) + z_nd(i+1,j,k) + z_nd(i,j+1,k) + z_nd(i+1,j+1,k));
                    Real z_base = Ampl * std::sin(kp * x - omega_t);

                    z -= z_base;
                    Real fac = std::sinh( kp * (z - H) ) / std::sinh(kp * H);

                    zvel_arr(i,j,k) += Ampl * omega * fac * std::cos(kp * x - omega_t);
                });
            }
        } // end xvel_err, yvel_err, zvel_err

        if (containerHasElement(plot_var_names, "pp_err"))
        {
            // Moving terrain ANALYTICAL
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<Real const>& p0_arr = p_hse.const_array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);

                const auto dx = geom[lev].CellSizeArray();
                const Array4<Real const>& z_nd = z_phys_nd[lev]->const_array(mfi);
                const Array4<Real const>& r0_arr = r_hse.const_array(mfi);

                Real H           = geom[lev].ProbHi()[2];
                Real Ampl        = 0.16;
                Real wavelength  = 100.;
                Real kp          = 2. * PI / wavelength;
                Real g           = CONST_GRAV;
                Real omega       = std::sqrt(g * kp);
                Real omega_t     = omega * t_new[lev];

                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                    derdat(i, j, k, mf_comp) = getPgivenRTh(rhotheta) - p0_arr(i,j,k);

                    Real rho_hse     = r0_arr(i,j,k);

                    Real x   = (i + 0.5) * dx[0];
                    Real z   = 0.125 * ( z_nd(i,j,k  ) + z_nd(i+1,j,k  ) + z_nd(i,j+1,k  ) + z_nd(i+1,j+1,k  )
                                        +z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1) );
                    Real z_base = Ampl * std::sin(kp * x - omega_t);

                    z -= z_base;
                    Real fac = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
                    Real pprime_exact = -(Ampl * omega * omega / kp) * fac *
                                              std::sin(kp * x - omega_t) * r0_arr(i,j,k);

                    derdat(i,j,k,mf_comp) -= pprime_exact;
                });
            }
            mf_comp += 1;
        }
#endif
    }

    // Fill terrain distortion MF
    if (solverChoice.use_terrain) {
        for (int lev(0); lev <= finest_level; ++lev) {
            MultiFab::Copy(mf_nd[lev],*z_phys_nd[lev],0,2,1,0);
            Real dz = Geom()[lev].CellSizeArray()[2];
            for (MFIter mfi(mf_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                Array4<      Real> mf_arr = mf_nd[lev].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    mf_arr(i,j,k,2) -= k * dz;
                });
            }
        }
    }

    std::string plotfilename;
    if (which == 1)
       plotfilename = Concatenate(plot_file_1, istep[0], 5);
    else if (which == 2)
       plotfilename = Concatenate(plot_file_2, istep[0], 5);

    // LSM writes it's own data
    if (which==1 && plot_lsm) {
        lsm.Plot_Lsm_Data(t_new[0], istep, refRatio());
    }

    if (finest_level == 0)
    {
        if (plotfile_type == "amrex") {
            amrex::Print() << "Writing native plotfile " << plotfilename << "\n";
            if (solverChoice.use_terrain) {
                WriteMultiLevelPlotfileWithTerrain(plotfilename, finest_level+1,
                                                   GetVecOfConstPtrs(mf),
                                                   GetVecOfConstPtrs(mf_nd),
                                                   varnames,
                                                   t_new[0], istep);
            } else {
                WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                        GetVecOfConstPtrs(mf),
                                        varnames,
                                        Geom(), t_new[0], istep, refRatio());
            }
            writeJobInfo(plotfilename);

#ifdef ERF_USE_PARTICLES
            particleData.Checkpoint(plotfilename);
#endif
#ifdef ERF_USE_HDF5
        } else if (plotfile_type == "hdf5" || plotfile_type == "HDF5") {
            amrex::Print() << "Writing plotfile " << plotfilename+"d01.h5" << "\n";
            WriteMultiLevelPlotfileHDF5(plotfilename, finest_level+1,
                                        GetVecOfConstPtrs(mf),
                                        varnames,
                                        Geom(), t_new[0], istep, refRatio());
#endif
#ifdef ERF_USE_NETCDF
        } else if (plotfile_type == "netcdf" || plotfile_type == "NetCDF") {
             int lev   = 0;
             int l_which = 0;
             writeNCPlotFile(lev, l_which, plotfilename, GetVecOfConstPtrs(mf), varnames, istep, t_new[0]);
#endif
        } else {
            amrex::Print() << "User specified plot_filetype = " << plotfile_type << std::endl;
            amrex::Abort("Dont know this plot_filetype");
        }

    } else { // multilevel

        Vector<IntVect>   r2(finest_level);
        Vector<Geometry>  g2(finest_level+1);
        Vector<MultiFab> mf2(finest_level+1);

        mf2[0].define(grids[0], dmap[0], ncomp_mf, 0);

        // Copy level 0 as is
        MultiFab::Copy(mf2[0],mf[0],0,0,mf[0].nComp(),0);

        // Define a new multi-level array of Geometry's so that we pass the new "domain" at lev > 0
        Array<int,AMREX_SPACEDIM> periodicity =
                     {Geom()[0].isPeriodic(0),Geom()[0].isPeriodic(1),Geom()[0].isPeriodic(2)};
        g2[0].define(Geom()[0].Domain(),&(Geom()[0].ProbDomain()),0,periodicity.data());

        if (plotfile_type == "amrex") {
            r2[0] = IntVect(1,1,ref_ratio[0][0]);
            for (int lev = 1; lev <= finest_level; ++lev) {
                if (lev > 1) {
                    r2[lev-1][0] = 1;
                    r2[lev-1][1] = 1;
                    r2[lev-1][2] = r2[lev-2][2] * ref_ratio[lev-1][0];
                }

                mf2[lev].define(refine(grids[lev],r2[lev-1]), dmap[lev], ncomp_mf, 0);

                // Set the new problem domain
                Box d2(Geom()[lev].Domain());
                d2.refine(r2[lev-1]);

                g2[lev].define(d2,&(Geom()[lev].ProbDomain()),0,periodicity.data());
            }

            // Do piecewise interpolation of mf into mf2
            for (int lev = 1; lev <= finest_level; ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf2[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    pcinterp_interp(bx,mf2[lev].array(mfi), 0, mf[lev].nComp(), mf[lev].const_array(mfi),0,r2[lev-1]);
                }
            }

            // Define an effective ref_ratio which is isotropic to be passed into WriteMultiLevelPlotfile
            Vector<IntVect> rr(finest_level);
            for (int lev = 0; lev < finest_level; ++lev) {
                rr[lev] = IntVect(ref_ratio[lev][0],ref_ratio[lev][1],ref_ratio[lev][0]);
            }

            amrex::Print() << "Writing plotfile " << plotfilename << "\n";
            if (solverChoice.use_terrain) {
                WriteMultiLevelPlotfileWithTerrain(plotfilename, finest_level+1,
                                                   GetVecOfConstPtrs(mf),
                                                   GetVecOfConstPtrs(mf_nd),
                                                   varnames,
                                                   t_new[0], istep);
            } else {
                WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                        GetVecOfConstPtrs(mf2), varnames,
                                        g2, t_new[0], istep, rr);
            }

            writeJobInfo(plotfilename);

#ifdef ERF_USE_PARTICLES
            particleData.Checkpoint(plotfilename);
#endif
#ifdef ERF_USE_NETCDF
        } else if (plotfile_type == "netcdf" || plotfile_type == "NetCDF") {
             for (int lev = 0; lev <= finest_level; ++lev) {
                 for (int which_box = 0; which_box < num_boxes_at_level[lev]; which_box++) {
                     writeNCPlotFile(lev, which_box, plotfilename, GetVecOfConstPtrs(mf), varnames, istep, t_new[0]);
                 }
             }
#endif
        }
    } // end multi-level
}

void
ERF::WriteMultiLevelPlotfileWithTerrain (const std::string& plotfilename, int nlevels,
                                         const Vector<const MultiFab*>& mf,
                                         const Vector<const MultiFab*>& mf_nd,
                                         const Vector<std::string>& varnames,
                                         Real time,
                                         const Vector<int>& level_steps,
                                         const std::string &versionName,
                                         const std::string &levelPrefix,
                                         const std::string &mfPrefix,
                                         const Vector<std::string>& extra_dirs) const
{
    BL_PROFILE("WriteMultiLevelPlotfileWithTerrain()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::MyProc() == ParallelDescriptor::NProcs()-1) {
        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        auto f = [=]() {
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
            std::string HeaderFileName(plotfilename + "/Header");
            std::ofstream HeaderFile;
            HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) FileOpenFailed(HeaderFileName);
            WriteGenericPlotfileHeaderWithTerrain(HeaderFile, nlevels, boxArrays, varnames,
                                                  time, level_steps, versionName,
                                                  levelPrefix, mfPrefix);
        };

        if (AsyncOut::UseAsyncOut()) {
            AsyncOut::Submit(std::move(f));
        } else {
            f();
        }
    }

    std::string mf_nodal_prefix = "Nu_nd";
    for (int level = 0; level <= finest_level; ++level)
    {
        if (AsyncOut::UseAsyncOut()) {
            VisMF::AsyncWrite(*mf[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix),
                              true);
            VisMF::AsyncWrite(*mf_nd[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix),
                              true);
        } else {
            const MultiFab* data;
            std::unique_ptr<MultiFab> mf_tmp;
            if (mf[level]->nGrowVect() != 0) {
                mf_tmp = std::make_unique<MultiFab>(mf[level]->boxArray(),
                                                    mf[level]->DistributionMap(),
                                                    mf[level]->nComp(), 0, MFInfo(),
                                                    mf[level]->Factory());
                MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
                data = mf_tmp.get();
            } else {
                data = mf[level];
            }
            VisMF::Write(*data        , MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
            VisMF::Write(*mf_nd[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix));
        }
    }
}

void
ERF::WriteGenericPlotfileHeaderWithTerrain (std::ostream &HeaderFile,
                                            int nlevels,
                                            const Vector<BoxArray> &bArray,
                                            const Vector<std::string> &varnames,
                                            Real time,
                                            const Vector<int> &level_steps,
                                            const std::string &versionName,
                                            const std::string &levelPrefix,
                                            const std::string &mfPrefix) const
{
        BL_ASSERT(nlevels <= bArray.size());
        BL_ASSERT(nlevels <= ref_ratio.size()+1);
        BL_ASSERT(nlevels <= level_steps.size());

        HeaderFile.precision(17);

        // ---- this is the generic plot file type name
        HeaderFile << versionName << '\n';

        HeaderFile << varnames.size() << '\n';

        for (int ivar = 0; ivar < varnames.size(); ++ivar) {
            HeaderFile << varnames[ivar] << "\n";
        }
        HeaderFile << AMREX_SPACEDIM << '\n';
        HeaderFile << time << '\n';
        HeaderFile << finest_level << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbLo(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbHi(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < finest_level; ++i) {
            HeaderFile << ref_ratio[i][0] << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << geom[i].Domain() << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << level_steps[i] << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
            }
            HeaderFile << '\n';
        }
        HeaderFile << (int) geom[0].Coord() << '\n';
        HeaderFile << "0\n";

        for (int level = 0; level <= finest_level; ++level) {
            HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
            HeaderFile << level_steps[level] << '\n';

            const IntVect& domain_lo = geom[level].Domain().smallEnd();
            for (int i = 0; i < bArray[level].size(); ++i)
            {
                // Need to shift because the RealBox ctor we call takes the
                // physical location of index (0,0,0).  This does not affect
                // the usual cases where the domain index starts with 0.
                const Box& b = shift(bArray[level][i], -domain_lo);
                RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
                }
            }

            HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
        }
        HeaderFile << "1" << "\n";
        HeaderFile << "3" << "\n";
        HeaderFile << "amrexvec_nu_x" << "\n";
        HeaderFile << "amrexvec_nu_y" << "\n";
        HeaderFile << "amrexvec_nu_z" << "\n";
        std::string mf_nodal_prefix = "Nu_nd";
        for (int level = 0; level <= finest_level; ++level) {
            HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_nodal_prefix) << '\n';
        }
}
