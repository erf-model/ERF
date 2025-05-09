#include <ERF_EOS.H>
#include <ERF.H>
#include "AMReX_Interp_3D_C.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_Constants.H"
#include "ERF_SrcHeaders.H"
#include "ERF_Container.H"

using namespace amrex;

PhysBCFunctNoOp null_bc_for_fill;

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

    for (int i = 0; i < cons_names.size(); ++i) {
        if ( containerHasElement(plot_var_names, cons_names[i]) ) {
            if (solverChoice.moisture_type == MoistureType::None) {
                if (cons_names[i] != "rhoQ1" && cons_names[i] != "rhoQ2" && cons_names[i] != "rhoQ3" &&
                    cons_names[i] != "rhoQ4" && cons_names[i] != "rhoQ5" && cons_names[i] != "rhoQ6")
                {
                    tmp_plot_names.push_back(cons_names[i]);
                }
            } else if (solverChoice.moisture_type == MoistureType::Kessler) { // allow rhoQ1, rhoQ2, rhoQ3
                if (cons_names[i] != "rhoQ4" && cons_names[i] != "rhoQ5" && cons_names[i] != "rhoQ6")
                {
                    tmp_plot_names.push_back(cons_names[i]);
                }
            } else if ( (solverChoice.moisture_type == MoistureType::SatAdj) ||
                        (solverChoice.moisture_type == MoistureType::SAM_NoPrecip_NoIce) ||
                        (solverChoice.moisture_type == MoistureType::Kessler_NoRain) ) { // allow rhoQ1, rhoQ2
                if (cons_names[i] != "rhoQ3" && cons_names[i] != "rhoQ4" &&
                    cons_names[i] != "rhoQ5" && cons_names[i] != "rhoQ6")
                {
                    tmp_plot_names.push_back(cons_names[i]);
                }
            } else if ( (solverChoice.moisture_type == MoistureType::Morrison_NoIce) ||
                        (solverChoice.moisture_type == MoistureType::SAM_NoIce     ) ) { // allow rhoQ1, rhoQ2, rhoQ4
                if (cons_names[i] != "rhoQ3" && cons_names[i] != "rhoQ5" && cons_names[i] != "rhoQ6")
                {
                    tmp_plot_names.push_back(cons_names[i]);
                }
            } else
            {
                // For moisture_type SAM and Morrison we have all six variables
                tmp_plot_names.push_back(cons_names[i]);
            }
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

    //
    // If the model we are running doesn't have the variable listed in the inputs file,
    //     just ignore it rather than aborting
    //
    for (int i = 0; i < derived_names.size(); ++i) {
        if ( containerHasElement(plot_var_names, derived_names[i]) ) {
            bool ok_to_add = ( (solverChoice.terrain_type == TerrainType::ImmersedForcing) ||
                               (derived_names[i] != "terrain_IB_mask") );
            ok_to_add     &= ( (SolverChoice::terrain_type == TerrainType::StaticFittedMesh) ||
                               (SolverChoice::terrain_type == TerrainType::MovingFittedMesh) ||
                               (derived_names[i] != "detJ") );
            ok_to_add     &= ( (SolverChoice::terrain_type == TerrainType::StaticFittedMesh) ||
                               (SolverChoice::terrain_type == TerrainType::MovingFittedMesh) ||
                               (derived_names[i] != "z_phys") );
#ifndef ERF_USE_WINDFARM
            ok_to_add     &= (derived_names[i] != "SMark0" && derived_names[i] != "SMark1");
#endif
            if (ok_to_add)
            {
                if (solverChoice.moisture_type == MoistureType::None) { // no moist quantities allowed
                    if (derived_names[i] != "qv" && derived_names[i] != "qc"    && derived_names[i] != "qrain"  &&
                        derived_names[i] != "qi" && derived_names[i] != "qsnow" && derived_names[i] != "qgraup" &&
                        derived_names[i] != "rain_accum" && derived_names[i] != "snow_accum" && derived_names[i] != "graup_accum")
                    {
                        tmp_plot_names.push_back(derived_names[i]);
                    }
                } else if ( (solverChoice.moisture_type == MoistureType::Kessler       ) ||
                            (solverChoice.moisture_type == MoistureType::Morrison_NoIce) ||
                            (solverChoice.moisture_type == MoistureType::SAM_NoIce     ) ) { // allow qv, qc, qrain
                    if (derived_names[i] != "qi" && derived_names[i] != "qsnow" && derived_names[i] != "qgraup" &&
                        derived_names[i] != "snow_accum" && derived_names[i] != "graup_accum")
                    {
                        tmp_plot_names.push_back(derived_names[i]);
                    }
                } else if ( (solverChoice.moisture_type == MoistureType::SatAdj) ||
                            (solverChoice.moisture_type == MoistureType::SAM_NoPrecip_NoIce) ||
                            (solverChoice.moisture_type == MoistureType::Kessler_NoRain) ) { // allow qv, qc
                    if (derived_names[i] != "qrain"  &&
                        derived_names[i] != "qi" && derived_names[i] != "qsnow" && derived_names[i] != "qgraup" &&
                        derived_names[i] != "rain_accum" && derived_names[i] != "snow_accum" && derived_names[i] != "graup_accum")
                    {
                        tmp_plot_names.push_back(derived_names[i]);
                    }
                } else
                {
                    // For moisture_type SAM and Morrison we have all moist quantities
                    tmp_plot_names.push_back(derived_names[i]);
                }
            } // use_terrain?
        } // hasElement
    }

#ifdef ERF_USE_WINDFARM
    for (int i = 0; i < derived_names.size(); ++i) {
        if ( containerHasElement(plot_var_names, derived_names[i]) ) {
            if(solverChoice.windfarm_type == WindFarmType::Fitch or solverChoice.windfarm_type == WindFarmType::EWP) {
                if(derived_names[i] == "num_turb" or derived_names[i] == "SMark0") {
                    tmp_plot_names.push_back(derived_names[i]);
                }
            }
            if( solverChoice.windfarm_type == WindFarmType::SimpleAD or
                solverChoice.windfarm_type == WindFarmType::GeneralAD ) {
                if(derived_names[i] == "num_turb" or derived_names[i] == "SMark0" or derived_names[i] == "SMark1") {
                    tmp_plot_names.push_back(derived_names[i]);
                }
            }
        }
    }
#endif

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

    {
        Vector<std::string> microphysics_plot_names;
        micro->GetPlotVarNames(microphysics_plot_names);
        if (microphysics_plot_names.size() > 0) {
            static bool first_call = true;
            if (first_call) {
                Print() << getEnumNameString(solverChoice.moisture_type)
                        << ": the following additional variables are available to plot:\n";
                for (int i = 0; i < microphysics_plot_names.size(); i++) {
                    Print() << "    " << microphysics_plot_names[i] << "\n";
                }
                first_call = false;
            }
            for (auto& plot_name : microphysics_plot_names) {
                if (containerHasElement(plot_var_names, plot_name)) {
                    tmp_plot_names.push_back(plot_name);
                }
            }
        }
    }

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
ERF::WritePlotFile (int which, PlotFileType plotfile_type, Vector<std::string> plot_var_names)
{
    const Vector<std::string> varnames = PlotFileVarNames(plot_var_names);
    const int ncomp_mf = varnames.size();

    int ncomp_cons = vars_new[0][Vars::cons].nComp();

    if (ncomp_mf == 0) return;

    // We Fillpatch here because some of the derived quantities require derivatives
    //     which require ghost cells to be filled.  We do not need to call FillPatcher
    //     because we don't need to set interior fine points.
    // NOTE: the momenta here are only used as scratch space, the momenta themselves are not fillpatched

    // Level 0 FilLPatch
    FillPatch(0, t_new[0], {&vars_new[0][Vars::cons], &vars_new[0][Vars::xvel],
                            &vars_new[0][Vars::yvel], &vars_new[0][Vars::zvel]});

    for (int lev = 1; lev <= finest_level; ++lev) {
        bool fillset = false;
        FillPatch(lev, t_new[lev], {&vars_new[lev][Vars::cons], &vars_new[lev][Vars::xvel],
                                    &vars_new[lev][Vars::yvel], &vars_new[lev][Vars::zvel]},
                                   {&vars_new[lev][Vars::cons],
                                    &rU_new[lev], &rV_new[lev], &rW_new[lev]},
                                    base_state[lev], base_state[lev], fillset);
    }

    // Get qmoist pointers if using moisture
    bool use_moisture = (solverChoice.moisture_type != MoistureType::None);
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int mvar(0); mvar<qmoist[lev].size(); ++mvar) {
            qmoist[lev][mvar] = micro->Get_Qmoist_Ptr(lev,mvar);
        }
    }

    // Vector of MultiFabs for cell-centered data
    Vector<MultiFab> mf(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp_mf, 0);
    }

    // Vector of MultiFabs for nodal data
    Vector<MultiFab> mf_nd(finest_level+1);
    if ( SolverChoice::mesh_type != MeshType::ConstantDz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            BoxArray nodal_grids(grids[lev]); nodal_grids.surroundingNodes();
            mf_nd[lev].define(nodal_grids, dmap[lev], 3, 0);
            mf_nd[lev].setVal(0.);
        }
    }

    // Vector of MultiFabs for face-centered velocity
    Vector<MultiFab> mf_u(finest_level+1);
    Vector<MultiFab> mf_v(finest_level+1);
    Vector<MultiFab> mf_w(finest_level+1);
    if (m_plot_face_vels) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            BoxArray grid_stag_u(grids[lev]); grid_stag_u.surroundingNodes(0);
            BoxArray grid_stag_v(grids[lev]); grid_stag_v.surroundingNodes(1);
            BoxArray grid_stag_w(grids[lev]); grid_stag_w.surroundingNodes(2);
            mf_u[lev].define(grid_stag_u, dmap[lev], 1, 0);
            mf_v[lev].define(grid_stag_v, dmap[lev], 1, 0);
            mf_w[lev].define(grid_stag_w, dmap[lev], 1, 0);
            MultiFab::Copy(mf_u[lev],vars_new[lev][Vars::xvel],0,0,1,0);
            MultiFab::Copy(mf_v[lev],vars_new[lev][Vars::yvel],0,0,1,0);
            MultiFab::Copy(mf_w[lev],vars_new[lev][Vars::zvel],0,0,1,0);
        }
    }

    // Array of MultiFabs for cell-centered velocity
    Vector<MultiFab> mf_cc_vel(finest_level+1);

    if (containerHasElement(plot_var_names, "x_velocity" ) ||
        containerHasElement(plot_var_names, "y_velocity" ) ||
        containerHasElement(plot_var_names, "z_velocity" ) ||
        containerHasElement(plot_var_names, "magvel"     ) ||
        containerHasElement(plot_var_names, "vorticity_x") ||
        containerHasElement(plot_var_names, "vorticity_y") ||
        containerHasElement(plot_var_names, "vorticity_z") ) {

        for (int lev = 0; lev <= finest_level; ++lev) {
            mf_cc_vel[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(1,1,1));
            mf_cc_vel[lev].setVal(-1.e20);
            average_face_to_cellcenter(mf_cc_vel[lev],0,
                                       Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],
                                                                &vars_new[lev][Vars::yvel],
                                                                &vars_new[lev][Vars::zvel]});
        } // lev
    } // if (vel or vort)

    // We need ghost cells if computing vorticity
    if ( containerHasElement(plot_var_names, "vorticity_x")||
         containerHasElement(plot_var_names, "vorticity_y") ||
         containerHasElement(plot_var_names, "vorticity_z") )
    {
        amrex::Interpolater* mapper = &cell_cons_interp;
        for (int lev = 1; lev <= finest_level; ++lev)
        {
            Vector<MultiFab*> fmf = {&(mf_cc_vel[lev]), &(mf_cc_vel[lev])};
            Vector<Real> ftime    = {t_new[lev], t_new[lev]};
            Vector<MultiFab*> cmf = {&mf_cc_vel[lev-1], &mf_cc_vel[lev-1]};
            Vector<Real> ctime    = {t_new[lev], t_new[lev]};

            FillBdyCCVels(mf_cc_vel,lev-1);

            // Call FillPatch which ASSUMES that all ghost cells at lev-1 have already been filled
            FillPatchTwoLevels(mf_cc_vel[lev], mf_cc_vel[lev].nGrowVect(), IntVect(0,0,0),
                               t_new[lev], cmf, ctime, fmf, ftime,
                               0, 0, mf_cc_vel[lev].nComp(), geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BaseBCVars::rho0_bc_comp);
        } // lev
        FillBdyCCVels(mf_cc_vel);
    } // if (vort)


    for (int lev = 0; lev <= finest_level; ++lev)
    {
        int mf_comp = 0;

        BoxArray ba(vars_new[lev][Vars::cons].boxArray());
        DistributionMapping dm = vars_new[lev][Vars::cons].DistributionMap();

        // First, copy any of the conserved state variables into the output plotfile
        for (int i = 0; i < cons_names.size(); ++i) {
            if (containerHasElement(plot_var_names, cons_names[i])) {
                MultiFab::Copy(mf[lev],vars_new[lev][Vars::cons],i,mf_comp,1,0);
                mf_comp++;
            }
        }

        // Next, check for velocities
        if (containerHasElement(plot_var_names, "x_velocity")) {
            MultiFab::Copy(mf[lev], mf_cc_vel[lev], 0, mf_comp, 1, 0);
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names, "y_velocity")) {
            MultiFab::Copy(mf[lev], mf_cc_vel[lev], 1, mf_comp, 1, 0);
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names, "z_velocity")) {
            MultiFab::Copy(mf[lev], mf_cc_vel[lev], 2, mf_comp, 1, 0);
            mf_comp += 1;
        }

        // Finally, check for any derived quantities and compute them, inserting
        // them into our output multifab
        auto calculate_derived = [&](const std::string& der_name,
                                     MultiFab& src_mf,
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
                    auto& sfab = src_mf[mfi];
                    der_function(bx, dfab, 0, 1, sfab, Geom(lev), t_new[0], nullptr, lev);
                }

                mf_comp++;
            }
        };

        bool ismoist = (solverChoice.moisture_type != MoistureType::None);

        // Note: All derived variables must be computed in order of "derived_names" defined in ERF.H
        calculate_derived("soundspeed",  vars_new[lev][Vars::cons], derived::erf_dersoundspeed);
        if (ismoist) {
            calculate_derived("temp",        vars_new[lev][Vars::cons], derived::erf_dermoisttemp);
        } else {
            calculate_derived("temp",        vars_new[lev][Vars::cons], derived::erf_dertemp);
        }
        calculate_derived("theta",       vars_new[lev][Vars::cons], derived::erf_dertheta);
        calculate_derived("KE",          vars_new[lev][Vars::cons], derived::erf_derKE);
        calculate_derived("scalar",      vars_new[lev][Vars::cons], derived::erf_derscalar);
        calculate_derived("vorticity_x", mf_cc_vel[lev]           , derived::erf_dervortx);
        calculate_derived("vorticity_y", mf_cc_vel[lev]           , derived::erf_dervorty);
        calculate_derived("vorticity_z", mf_cc_vel[lev]           , derived::erf_dervortz);
        calculate_derived("magvel"     , mf_cc_vel[lev]           , derived::erf_dermagvel);

        MultiFab  r_hse(base_state[lev], make_alias, BaseState::r0_comp , 1);
        MultiFab  p_hse(base_state[lev], make_alias, BaseState::p0_comp , 1);
        MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);

        MultiFab pressure;

        if (solverChoice.anelastic[lev] == 0) {
            if (containerHasElement(plot_var_names, "pressure")  ||
                containerHasElement(plot_var_names, "pert_pres") ||
                containerHasElement(plot_var_names, "dpdx")      ||
                containerHasElement(plot_var_names, "dpdy")      ||
                containerHasElement(plot_var_names, "dpdz"))
            {
                int ng = (containerHasElement(plot_var_names, "dpdx") || containerHasElement(plot_var_names, "dpdy") ||
                          containerHasElement(plot_var_names, "dpdz")) ? 1 : 0;

                // Allocate space for pressure
                pressure.define(ba,dm,1,ng);

                if (ng > 0) {
                    // Default to p_hse as a way of filling ghost cells at domain boundaries
                    MultiFab::Copy(pressure,p_hse,0,0,1,1);
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    #endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& gbx = mfi.growntilebox(IntVect(1,1,0));

                    const Array4<Real      >& p_arr = pressure.array(mfi);
                    const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                    const int ncomp = vars_new[lev][Vars::cons].nComp();

                    ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        Real qv_for_p = (use_moisture && (ncomp > RhoQ1_comp)) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0;
                        const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                        p_arr(i, j, k) = getPgivenRTh(rhotheta,qv_for_p);
                    });
               } // mfi
               pressure.FillBoundary(geom[lev].periodicity());
            } // compute compressible pressure
        } // not anelastic

        if (containerHasElement(plot_var_names, "divU"))
        {
            MultiFab dmf(mf[lev], make_alias, mf_comp, 1);
            Array<MultiFab const*, AMREX_SPACEDIM> u;
            u[0] = &(vars_new[lev][Vars::xvel]);
            u[1] = &(vars_new[lev][Vars::yvel]);
            u[2] = &(vars_new[lev][Vars::zvel]);
            compute_divergence (lev, dmf, u, geom[lev]);
            mf_comp += 1;
        }

        if (containerHasElement(plot_var_names, "pres_hse"))
        {
            MultiFab::Copy(mf[lev],p_hse,0,mf_comp,1,0);
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names, "dens_hse"))
        {
            MultiFab::Copy(mf[lev],r_hse,0,mf_comp,1,0);
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names, "theta_hse"))
        {
            MultiFab::Copy(mf[lev],th_hse,0,mf_comp,1,0);
            mf_comp += 1;
        }

        if (containerHasElement(plot_var_names, "pressure"))
        {
            if (solverChoice.anelastic[lev] == 1) {
                MultiFab::Copy(mf[lev], p_hse, 0, mf_comp, 1, 0);
            } else {
                MultiFab::Copy(mf[lev], pressure, 0, mf_comp, 1, 0);
            }

            mf_comp += 1;
        } // pressure

        if (containerHasElement(plot_var_names, "pert_pres"))
        {
            if (solverChoice.anelastic[lev] == 1) {
                MultiFab::Copy(mf[lev], pp_inc[lev], 0, mf_comp, 1, 0);
            } else {
                MultiFab::Copy(mf[lev], pressure, 0, mf_comp, 1, 0);
                MultiFab::Subtract(mf[lev],p_hse,0,mf_comp,1,IntVect{0});
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
                const Array4<Real const>& p_arr = pressure.const_array(mfi);
                const int ncomp = vars_new[lev][Vars::cons].nComp();
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    Real qv = (use_moisture && (ncomp > RhoQ1_comp)) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0.0;
                    Real qc = (use_moisture && (ncomp > RhoQ2_comp)) ? S_arr(i,j,k,RhoQ2_comp)/S_arr(i,j,k,Rho_comp) : 0.0;
                    Real T = getTgivenRandRTh(S_arr(i,j,k,Rho_comp), S_arr(i,j,k,RhoTheta_comp), qv);
                    Real fac = Cp_d + Cp_l*(qv + qc);
                    Real pv = erf_esatw(T)*100.0;

                    derdat(i, j, k, mf_comp) = T*std::pow((p_arr(i,j,k) - pv)/p_0, -R_d/fac)*std::exp(L_v*qv/(fac*T)) ;
                });
            }
            mf_comp ++;
        }

        if (containerHasElement(plot_var_names, "terrain_IB_mask"))
        {
            MultiFab* terrain_blank = terrain_blanking[lev].get();
            MultiFab::Copy(mf[lev],*terrain_blank,0,mf_comp,1,0);
            mf_comp ++;
        }

#ifdef ERF_USE_WINDFARM
        if ( containerHasElement(plot_var_names, "num_turb") and
             (solverChoice.windfarm_type == WindFarmType::Fitch or solverChoice.windfarm_type == WindFarmType::EWP or
              solverChoice.windfarm_type == WindFarmType::SimpleAD or solverChoice.windfarm_type == WindFarmType::GeneralAD) )
        {
            MultiFab::Copy(mf[lev],Nturb[lev],0,mf_comp,1,0);
            mf_comp ++;
        }

        if ( containerHasElement(plot_var_names, "SMark0") and
             (solverChoice.windfarm_type == WindFarmType::Fitch or solverChoice.windfarm_type == WindFarmType::EWP or
              solverChoice.windfarm_type == WindFarmType::SimpleAD or solverChoice.windfarm_type == WindFarmType::GeneralAD) )
        {
            MultiFab::Copy(mf[lev],SMark[lev],0,mf_comp,1,0);
            mf_comp ++;
        }

        if (containerHasElement(plot_var_names, "SMark1") and
           (solverChoice.windfarm_type == WindFarmType::SimpleAD or solverChoice.windfarm_type == WindFarmType::GeneralAD))
        {
            MultiFab::Copy(mf[lev],SMark[lev],1,mf_comp,1,0);
            mf_comp ++;
        }

#endif

        // **********************************************************************************************
        // Allocate space if we are computing any pressure gradients
        // **********************************************************************************************

        Vector<MultiFab> gradp_temp;  gradp_temp.resize(AMREX_SPACEDIM);
        if (containerHasElement(plot_var_names, "dpdx")       ||
            containerHasElement(plot_var_names, "dpdy")       ||
            containerHasElement(plot_var_names, "dpdz")       ||
            containerHasElement(plot_var_names, "pres_hse_x") ||
            containerHasElement(plot_var_names, "pres_hse_y"))
        {
            gradp_temp[GpVars::gpx].define(convert(ba, IntVect(1,0,0)), dm, 1, 1); gradp_temp[GpVars::gpx].setVal(0.);
            gradp_temp[GpVars::gpy].define(convert(ba, IntVect(0,1,0)), dm, 1, 1); gradp_temp[GpVars::gpy].setVal(0.);
            gradp_temp[GpVars::gpz].define(convert(ba, IntVect(0,0,1)), dm, 1, 1); gradp_temp[GpVars::gpz].setVal(0.);
        }

        // **********************************************************************************************
        // These are based on computing gradient of full pressure
        // **********************************************************************************************

        if ( (containerHasElement(plot_var_names, "dpdx")) ||
             (containerHasElement(plot_var_names, "dpdy")) ||
             (containerHasElement(plot_var_names, "dpdz")) ) {
            compute_gradp(pressure, geom[lev], z_phys_nd[lev], z_phys_cc[lev], gradp_temp, solverChoice);
        }

        if (containerHasElement(plot_var_names, "dpdx"))
        {
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real      >&  derdat  = mf[lev].array(mfi);
                const Array4<Real const>&  gpx_arr = gradp_temp[GpVars::gpx].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i ,j ,k, mf_comp) = 0.5 * (gpx_arr(i+1,j,k) + gpx_arr(i,j,k));
                });
            }
            mf_comp ++;
        } // dpdx
        if (containerHasElement(plot_var_names, "dpdy"))
        {
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real      >&  derdat  = mf[lev].array(mfi);
                const Array4<Real const>&  gpy_arr = gradp_temp[GpVars::gpy].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i ,j ,k, mf_comp) = 0.5 * (gpy_arr(i,j+1,k) + gpy_arr(i,j,k));
                });
            }
            mf_comp ++;
        } // dpdy
        if (containerHasElement(plot_var_names, "dpdz"))
        {
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real      >&  derdat  = mf[lev].array(mfi);
                const Array4<Real const>&  gpz_arr = gradp_temp[GpVars::gpz].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i ,j ,k, mf_comp) = 0.5 * (gpz_arr(i,j,k+1) + gpz_arr(i,j,k));
                });
            }
            mf_comp ++;
        } // dpdz

        // **********************************************************************************************
        // These are based on computing gradient of basestate pressure
        // **********************************************************************************************

        if ( (containerHasElement(plot_var_names, "pres_hse_x")) ||
             (containerHasElement(plot_var_names, "pres_hse_y")) ) {
            compute_gradp(p_hse, geom[lev], z_phys_nd[lev], z_phys_cc[lev], gradp_temp, solverChoice);
        }

        if (containerHasElement(plot_var_names, "pres_hse_x"))
        {
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real      >&  derdat  = mf[lev].array(mfi);
                const Array4<Real const>&  gpx_arr = gradp_temp[0].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i ,j ,k, mf_comp) = 0.5 * (gpx_arr(i+1,j,k) + gpx_arr(i,j,k));
                });
            }
            mf_comp += 1;
        } // pres_hse_x

        if (containerHasElement(plot_var_names, "pres_hse_y"))
        {
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real      >&  derdat  = mf[lev].array(mfi);
                const Array4<Real const>&  gpy_arr = gradp_temp[1].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i ,j ,k, mf_comp) = 0.5 * (gpy_arr(i,j+1,k) + gpy_arr(i,j,k));
                });
            }
            mf_comp += 1;
        } // pres_hse_y

        // **********************************************************************************************
        // Metric terms
        // **********************************************************************************************

        if (SolverChoice::mesh_type != MeshType::ConstantDz) {
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

#ifdef ERF_USE_NETCDF
        if (solverChoice.use_real_bcs) {
            if (containerHasElement(plot_var_names, "lat_m")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat = mf[lev].array(mfi);
                    const Array4<Real>& data   = lat_m[lev]->array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = data(i,j,0);
                    });
                }
                mf_comp ++;
            } // lat_m

            if (containerHasElement(plot_var_names, "lon_m")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat = mf[lev].array(mfi);
                    const Array4<Real>& data   = lon_m[lev]->array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = data(i,j,0);
                    });
                }
                mf_comp ++;
            } // lon_m
        } // use_real_bcs
#endif


        if (solverChoice.time_avg_vel) {
            if (containerHasElement(plot_var_names, "u_t_avg")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat = mf[lev].array(mfi);
                    const Array4<Real>& data   = vel_t_avg[lev]->array(mfi);
                    const Real norm = t_avg_cnt[lev];
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        derdat(i ,j ,k, mf_comp) = data(i,j,k,0) / norm;
                    });
                }
                mf_comp ++;
            }

            if (containerHasElement(plot_var_names, "v_t_avg")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat = mf[lev].array(mfi);
                    const Array4<Real>& data   = vel_t_avg[lev]->array(mfi);
                    const Real norm = t_avg_cnt[lev];
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        derdat(i ,j ,k, mf_comp) = data(i,j,k,1) / norm;
                    });
                }
                mf_comp ++;
            }

            if (containerHasElement(plot_var_names, "w_t_avg")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat = mf[lev].array(mfi);
                    const Array4<Real>& data   = vel_t_avg[lev]->array(mfi);
                    const Real norm = t_avg_cnt[lev];
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        derdat(i ,j ,k, mf_comp) = data(i,j,k,2) / norm;
                    });
                }
                mf_comp ++;
            }

            if (containerHasElement(plot_var_names, "umag_t_avg")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat = mf[lev].array(mfi);
                    const Array4<Real>& data   = vel_t_avg[lev]->array(mfi);
                    const Real norm = t_avg_cnt[lev];
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        derdat(i ,j ,k, mf_comp) = data(i,j,k,3) / norm;
                    });
                }
                mf_comp ++;
            }
        }

        if (containerHasElement(plot_var_names, "nut")) {
            MultiFab dmf(mf[lev], make_alias, mf_comp, 1);
            MultiFab cmf(vars_new[lev][Vars::cons], make_alias, 0, 1); // to provide rho only
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                auto       prim = dmf[mfi].array();
                auto const cons = cmf[mfi].const_array();
                auto const diff = (*eddyDiffs_lev[lev])[mfi].const_array();
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    const Real rho = cons(i, j, k, Rho_comp);
                    const Real Kmv = diff(i, j, k, EddyDiff::Mom_v);
                    prim(i,j,k) = Kmv / rho;
                });
            }

            mf_comp++;
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
        if (containerHasElement(plot_var_names, "Lturb")) {
            MultiFab::Copy(mf[lev],*eddyDiffs_lev[lev],EddyDiff::Turb_lengthscale,mf_comp,1,0);
            mf_comp ++;
        }
        if (containerHasElement(plot_var_names, "walldist")) {
            MultiFab::Copy(mf[lev],*walldist[lev],0,mf_comp,1,0);
            mf_comp ++;
        }
        if (containerHasElement(plot_var_names, "diss")) {
            MultiFab::Copy(mf[lev],*SFS_diss_lev[lev],0,mf_comp,1,0);
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
            int n_qstate_moist   = micro->Get_Qstate_Moist_Size();

            // Moist density
            if(containerHasElement(plot_var_names, "moist_density"))
            {
                int n_start = RhoQ1_comp; // qv
                int n_end   = RhoQ2_comp; // qc
                if (n_qstate_moist > 3) n_end = RhoQ3_comp; // qi
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], Rho_comp, mf_comp, 1, 0);
                for (int n_comp(n_start); n_comp <= n_end; ++n_comp) {
                    MultiFab::Add(mf[lev], vars_new[lev][Vars::cons], n_comp, mf_comp, 1, 0);
                }
                mf_comp += 1;
            }

            // Non-precipitating components
            //--------------------------------------------------------------------------
            if(containerHasElement(plot_var_names, "qt"))
            {
                int n_start = RhoQ1_comp; // qv
                int n_end   = RhoQ2_comp; // qc
                if (n_qstate_moist > 3) n_end = RhoQ3_comp; // qi
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], n_start, mf_comp, 1, 0);
                for (int n_comp(n_start+1); n_comp <= n_end; ++n_comp) {
                    MultiFab::Add(mf[lev], vars_new[lev][Vars::cons], n_comp, mf_comp, 1, 0);
                }
                MultiFab::Divide(mf[lev], vars_new[lev][Vars::cons], Rho_comp  , mf_comp, 1, 0);
                mf_comp += 1;
            }

            if(containerHasElement(plot_var_names, "qv") && (n_qstate_moist >= 1))
            {
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], RhoQ1_comp, mf_comp, 1, 0);
                MultiFab::Divide(mf[lev], vars_new[lev][Vars::cons], Rho_comp  , mf_comp, 1, 0);
                mf_comp += 1;
            }

            if(containerHasElement(plot_var_names, "qc") && (n_qstate_moist >= 2))
            {
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], RhoQ2_comp, mf_comp, 1, 0);
                MultiFab::Divide(mf[lev], vars_new[lev][Vars::cons], Rho_comp  , mf_comp, 1, 0);
                mf_comp += 1;
            }

            if(containerHasElement(plot_var_names, "qi") && (n_qstate_moist >= 4))
            {
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], RhoQ3_comp, mf_comp, 1, 0);
                MultiFab::Divide(mf[lev], vars_new[lev][Vars::cons], Rho_comp  , mf_comp, 1, 0);
                mf_comp += 1;
            }

            // Precipitating components
            //--------------------------------------------------------------------------
            if(containerHasElement(plot_var_names, "qp") && (n_qstate_moist >= 3))
            {
                int n_start = (n_qstate_moist > 3) ? RhoQ4_comp : RhoQ3_comp;
                int n_end   = ncomp_cons - 1;
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], n_start, mf_comp, 1, 0);
                for (int n_comp(n_start+1); n_comp <= n_end; ++n_comp) {
                    MultiFab::Add(  mf[lev], vars_new[lev][Vars::cons], n_comp, mf_comp, 1, 0);
                }
                MultiFab::Divide(mf[lev], vars_new[lev][Vars::cons], Rho_comp  , mf_comp, 1, 0);
                mf_comp += 1;
            }

            if(containerHasElement(plot_var_names, "qrain") && (n_qstate_moist >= 3))
            {
                int n_start = (n_qstate_moist > 3) ? RhoQ4_comp : RhoQ3_comp;
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], n_start , mf_comp, 1, 0);
                MultiFab::Divide(mf[lev], vars_new[lev][Vars::cons], Rho_comp, mf_comp, 1, 0);
                mf_comp += 1;
            }

            if(containerHasElement(plot_var_names, "qsnow") && (n_qstate_moist >= 5))
            {
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], RhoQ5_comp, mf_comp, 1, 0);
                MultiFab::Divide(mf[lev], vars_new[lev][Vars::cons],   Rho_comp, mf_comp, 1, 0);
                mf_comp += 1;
            }

            if(containerHasElement(plot_var_names, "qgraup") && (n_qstate_moist >= 6))
            {
                MultiFab::Copy(  mf[lev], vars_new[lev][Vars::cons], RhoQ6_comp, mf_comp, 1, 0);
                MultiFab::Divide(mf[lev], vars_new[lev][Vars::cons],   Rho_comp, mf_comp, 1, 0);
                mf_comp += 1;
            }

            if (containerHasElement(plot_var_names, "qsat"))
            {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat  = mf[lev].array(mfi);
                    const Array4<Real const>& p_arr = pressure.array(mfi);
                    const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        Real       qv = S_arr(i,j,k,RhoQ1_comp) / S_arr(i,j,k,Rho_comp);
                        Real       T  = getTgivenRandRTh(S_arr(i,j,k,Rho_comp), S_arr(i,j,k,RhoTheta_comp), qv);
                        Real pressure = p_arr(i,j,k) * Real(0.01);
                        erf_qsatw(T, pressure, derdat(i,j,k,mf_comp));
                    });
                }
                mf_comp ++;
            }

            if ( (solverChoice.moisture_type == MoistureType::Kessler) ||
                 (solverChoice.moisture_type == MoistureType::Morrison_NoIce) ||
                 (solverChoice.moisture_type == MoistureType::SAM_NoIce) )
            {
                int offset = (solverChoice.moisture_type == MoistureType::Morrison_NoIce) ? 5 : 0;
                if (containerHasElement(plot_var_names, "rain_accum"))
                {
                    MultiFab::Copy(mf[lev],*(qmoist[lev][offset]),0,mf_comp,1,0);
                    mf_comp += 1;
                }
            }
            else if ( (solverChoice.moisture_type == MoistureType::SAM) ||
                      (solverChoice.moisture_type == MoistureType::Morrison) )
            {
                int offset = (solverChoice.moisture_type == MoistureType::Morrison) ? 5 : 0;
                if (containerHasElement(plot_var_names, "rain_accum"))
                {
                    MultiFab::Copy(mf[lev],*(qmoist[lev][offset]),0,mf_comp,1,0);
                    mf_comp += 1;
                }
                if (containerHasElement(plot_var_names, "snow_accum"))
                {
                    MultiFab::Copy(mf[lev],*(qmoist[lev][offset+1]),0,mf_comp,1,0);
                    mf_comp += 1;
                }
                if (containerHasElement(plot_var_names, "graup_accum"))
                {
                    MultiFab::Copy(mf[lev],*(qmoist[lev][offset+2]),0,mf_comp,1,0);
                    mf_comp += 1;
                }
            }
        } // use_moisture

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

        {
            Vector<std::string> microphysics_plot_names;
            micro->GetPlotVarNames(microphysics_plot_names);
            for (auto& plot_name : microphysics_plot_names) {
                if (containerHasElement(plot_var_names, plot_name)) {
                    MultiFab temp_dat(mf[lev].boxArray(), mf[lev].DistributionMap(), 1, 1);
                    temp_dat.setVal(0);
                    micro->GetPlotVar(plot_name, temp_dat, lev);
                    MultiFab::Copy(mf[lev], temp_dat, 0, mf_comp, 1, 0);
                    mf_comp += 1;
                }
            }
        }

        if (containerHasElement(plot_var_names, "volfrac")) {
            if ( solverChoice.terrain_type == TerrainType::EB ||
                 solverChoice.terrain_type == TerrainType::ImmersedForcing)
            {
                MultiFab::Copy(mf[lev], EBFactory(lev).getVolFrac(), 0, mf_comp, 1, 0);
            } else {
                mf[lev].setVal(1.0, mf_comp, 1, 0);
            }
            mf_comp += 1;
        }

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
                const Array4<Real const>&  p_arr = pressure.const_array(mfi);
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
                    derdat(i, j, k, mf_comp) = p_arr(i,j,k) - p0_arr(i,j,k);

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

#ifdef ERF_USE_RRTMGP
    if (containerHasElement(plot_var_names, "qsrc_sw")) {
        MultiFab::Copy(mf[lev], *(qheating_rates[lev]), 0, mf_comp, 1, 0);
        mf_comp += 1;
    }
    if (containerHasElement(plot_var_names, "qsrc_lw")) {
        MultiFab::Copy(mf[lev], *(qheating_rates[lev]), 1, mf_comp, 1, 0);
        mf_comp += 1;
    }
#endif
    }

    if (solverChoice.terrain_type == TerrainType::EB)
    {
        for (int lev = 0; lev <= finest_level; ++lev) {
            EB_set_covered(mf[lev], 0.0);
        }
    }

    // Fill terrain distortion MF (nu_nd)
    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
        for (int lev(0); lev <= finest_level; ++lev) {
            MultiFab::Copy(mf_nd[lev],*z_phys_nd[lev],0,2,1,0);
            Real dz = Geom()[lev].CellSizeArray()[2];
            for (MFIter mfi(mf_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                Array4<Real> mf_arr = mf_nd[lev].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    mf_arr(i,j,k,2) -= k * dz;
                });
            }
        }
    }

    std::string plotfilename;
    std::string plotfilenameU;
    std::string plotfilenameV;
    std::string plotfilenameW;
    if (which == 1) {
       plotfilename = Concatenate(plot_file_1, istep[0], 5);
       plotfilenameU = Concatenate(plot_file_1+"U", istep[0], 5);
       plotfilenameV = Concatenate(plot_file_1+"V", istep[0], 5);
       plotfilenameW = Concatenate(plot_file_1+"W", istep[0], 5);
    } else if (which == 2) {
       plotfilename = Concatenate(plot_file_2, istep[0], 5);
       plotfilenameU = Concatenate(plot_file_2+"U", istep[0], 5);
       plotfilenameV = Concatenate(plot_file_2+"V", istep[0], 5);
       plotfilenameW = Concatenate(plot_file_2+"W", istep[0], 5);
    }

    // LSM writes it's own data
    if (which==1 && plot_lsm) {
        lsm.Plot_Lsm_Data(t_new[0], istep, refRatio());
    }

#ifdef ERF_USE_RRTMGP
    /*
    // write additional RRTMGP data
    // TODO: currently single level only
    if (which==1 && plot_rad) {
        rad[0]->writePlotfile(plot_file_1, t_new[0], istep[0]);
    }
    */
#endif

    // Single level
    if (finest_level == 0)
    {
        if (plotfile_type == PlotFileType::Amrex)
        {
            Print() << "Writing native plotfile " << plotfilename << "\n";
            if (SolverChoice::mesh_type != MeshType::ConstantDz) {
                WriteMultiLevelPlotfileWithTerrain(plotfilename, finest_level+1,
                                                   GetVecOfConstPtrs(mf),
                                                   GetVecOfConstPtrs(mf_nd),
                                                   varnames,
                                                   Geom(), t_new[0], istep, refRatio());
            } else {
                WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                        GetVecOfConstPtrs(mf),
                                        varnames,
                                        Geom(), t_new[0], istep, refRatio());
            }
            writeJobInfo(plotfilename);

            if (m_plot_face_vels) {
                Print() << "Writing face velocities" << std::endl;
                WriteMultiLevelPlotfile(plotfilenameU, finest_level+1,
                                        GetVecOfConstPtrs(mf_u),
                                        {"x_velocity_stag"},
                                        Geom(), t_new[0], istep, refRatio());
                WriteMultiLevelPlotfile(plotfilenameV, finest_level+1,
                                        GetVecOfConstPtrs(mf_v),
                                        {"y_velocity_stag"},
                                        Geom(), t_new[0], istep, refRatio());
                WriteMultiLevelPlotfile(plotfilenameW, finest_level+1,
                                        GetVecOfConstPtrs(mf_w),
                                        {"z_velocity_stag"},
                                        Geom(), t_new[0], istep, refRatio());
            }

#ifdef ERF_USE_PARTICLES
            particleData.writePlotFile(plotfilename);
#endif
#ifdef ERF_USE_NETCDF
        } else if (plotfile_type == PlotFileType::Netcdf) {
             int lev   = 0;
             int l_which = 0;
             writeNCPlotFile(lev, l_which, plotfilename, GetVecOfConstPtrs(mf), varnames, istep, t_new[0]);
#endif
        } else {
            // Here we assume the plotfile_type is PlotFileType::None
            Print() << "Writing no plotfile since plotfile_type is none" << std::endl;
        }

    } else { // Multilevel

        if (plotfile_type == PlotFileType::Amrex) {

            int lev0 = 0;
            int desired_ratio = std::max(std::max(ref_ratio[lev0][0],ref_ratio[lev0][1]),ref_ratio[lev0][2]);
            bool any_ratio_one = ( ( (ref_ratio[lev0][0] == 1) || (ref_ratio[lev0][1] == 1) ) ||
                                     (ref_ratio[lev0][2] == 1) );
            for (int lev = 1; lev < finest_level; lev++) {
                any_ratio_one = any_ratio_one ||
                                     ( ( (ref_ratio[lev][0] == 1) || (ref_ratio[lev][1] == 1) ) ||
                                         (ref_ratio[lev][2] == 1) );
            }

            if (any_ratio_one && m_expand_plotvars_to_unif_rr)
            {
                Vector<IntVect>   r2(finest_level);
                Vector<Geometry>  g2(finest_level+1);
                Vector<MultiFab> mf2(finest_level+1);

                mf2[0].define(grids[0], dmap[0], ncomp_mf, 0);

                // Copy level 0 as is
                MultiFab::Copy(mf2[0],mf[0],0,0,mf[0].nComp(),0);

                // Define a new multi-level array of Geometry's so that we pass the new "domain" at lev > 0
                Array<int,AMREX_SPACEDIM> periodicity =
                             {Geom()[lev0].isPeriodic(0),Geom()[lev0].isPeriodic(1),Geom()[lev0].isPeriodic(2)};
                g2[lev0].define(Geom()[lev0].Domain(),&(Geom()[lev0].ProbDomain()),0,periodicity.data());

                r2[0] = IntVect(desired_ratio/ref_ratio[lev0][0],
                                desired_ratio/ref_ratio[lev0][1],
                                desired_ratio/ref_ratio[lev0][2]);

                for (int lev = 1; lev <= finest_level; ++lev) {
                    if (lev > 1) {
                        r2[lev-1][0] = r2[lev-2][0] * desired_ratio / ref_ratio[lev-1][0];
                        r2[lev-1][1] = r2[lev-2][1] * desired_ratio / ref_ratio[lev-1][1];
                        r2[lev-1][2] = r2[lev-2][2] * desired_ratio / ref_ratio[lev-1][2];
                    }

                    mf2[lev].define(refine(grids[lev],r2[lev-1]), dmap[lev], ncomp_mf, 0);

                    // Set the new problem domain
                    Box d2(Geom()[lev].Domain());
                    d2.refine(r2[lev-1]);

                    g2[lev].define(d2,&(Geom()[lev].ProbDomain()),0,periodicity.data());
                }

                //
                // We need to make a temporary that is the size of ncomp_mf
                // in order to not get an out of bounds error
                // even though the values will not be used
                //
                Vector<BCRec> temp_domain_bcs_type;
                temp_domain_bcs_type.resize(ncomp_mf);

                //
                // Do piecewise constant interpolation of mf into mf2
                //
                for (int lev = 1; lev <= finest_level; ++lev) {
                    Interpolater* mapper_c = &pc_interp;
                    InterpFromCoarseLevel(mf2[lev], t_new[lev], mf[lev],
                                          0, 0, ncomp_mf,
                                          geom[lev], g2[lev],
                                          null_bc_for_fill, 0, null_bc_for_fill, 0,
                                          r2[lev-1], mapper_c, temp_domain_bcs_type, 0);
                }

                // Define an effective ref_ratio which is isotropic to be passed into WriteMultiLevelPlotfile
                Vector<IntVect> rr(finest_level);
                for (int lev = 0; lev < finest_level; ++lev) {
                    rr[lev] = IntVect(desired_ratio);
                }

               Print() << "Writing plotfile " << plotfilename << "\n";
               if (SolverChoice::mesh_type != MeshType::ConstantDz) {
                   WriteMultiLevelPlotfileWithTerrain(plotfilename, finest_level+1,
                                                      GetVecOfConstPtrs(mf2),
                                                      GetVecOfConstPtrs(mf_nd),
                                                      varnames,
                                                      g2, t_new[0], istep, rr);
               } else {
                   WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                           GetVecOfConstPtrs(mf2), varnames,
                                           g2, t_new[0], istep, rr);
               }

            } else {
               if (SolverChoice::mesh_type != MeshType::ConstantDz) {
                    WriteMultiLevelPlotfileWithTerrain(plotfilename, finest_level+1,
                                                       GetVecOfConstPtrs(mf),
                                                       GetVecOfConstPtrs(mf_nd),
                                                       varnames,
                                                       geom, t_new[0], istep, ref_ratio);
                } else {
                    WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                            GetVecOfConstPtrs(mf), varnames,
                                            geom, t_new[0], istep, ref_ratio);
                }
                if (m_plot_face_vels) {
                    Print() << "Writing face velocities" << std::endl;
                    WriteMultiLevelPlotfile(plotfilenameU, finest_level+1,
                                            GetVecOfConstPtrs(mf_u),
                                            {"x_velocity_stag"},
                                            geom, t_new[0], istep, ref_ratio);
                    WriteMultiLevelPlotfile(plotfilenameV, finest_level+1,
                                            GetVecOfConstPtrs(mf_v),
                                            {"y_velocity_stag"},
                                            geom, t_new[0], istep, ref_ratio);
                    WriteMultiLevelPlotfile(plotfilenameW, finest_level+1,
                                            GetVecOfConstPtrs(mf_w),
                                            {"z_velocity_stag"},
                                            geom, t_new[0], istep, ref_ratio);
                }
            } // ref_ratio test

            writeJobInfo(plotfilename);

#ifdef ERF_USE_PARTICLES
            particleData.writePlotFile(plotfilename);
#endif

#ifdef ERF_USE_NETCDF
        } else if (plotfile_type == PlotFileType::Netcdf) {
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
                                         const Vector<Geometry>& my_geom,
                                         Real time,
                                         const Vector<int>& level_steps,
                                         const Vector<IntVect>& rr,
                                         const std::string &versionName,
                                         const std::string &levelPrefix,
                                         const std::string &mfPrefix,
                                         const Vector<std::string>& extra_dirs) const
{
    BL_PROFILE("WriteMultiLevelPlotfileWithTerrain()");

    AMREX_ALWAYS_ASSERT(nlevels <= mf.size());
    AMREX_ALWAYS_ASSERT(nlevels <= rr.size()+1);
    AMREX_ALWAYS_ASSERT(nlevels <= level_steps.size());
    AMREX_ALWAYS_ASSERT(mf[0]->nComp() == varnames.size());

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
                                                  my_geom, time, level_steps, rr, versionName,
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
                                            const Vector<Geometry>& my_geom,
                                            Real my_time,
                                            const Vector<int>& level_steps,
                                            const Vector<IntVect>& my_ref_ratio,
                                            const std::string &versionName,
                                            const std::string &levelPrefix,
                                            const std::string &mfPrefix) const
{
    AMREX_ALWAYS_ASSERT(nlevels <= bArray.size());
    AMREX_ALWAYS_ASSERT(nlevels <= my_ref_ratio.size()+1);
    AMREX_ALWAYS_ASSERT(nlevels <= level_steps.size());

    HeaderFile.precision(17);

    // ---- this is the generic plot file type name
    HeaderFile << versionName << '\n';

    HeaderFile << varnames.size() << '\n';

    for (int ivar = 0; ivar < varnames.size(); ++ivar) {
        HeaderFile << varnames[ivar] << "\n";
    }
    HeaderFile << AMREX_SPACEDIM << '\n';
    HeaderFile << my_time << '\n';
    HeaderFile << finest_level << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        HeaderFile << my_geom[0].ProbLo(i) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        HeaderFile << my_geom[0].ProbHi(i) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < finest_level; ++i) {
        HeaderFile << my_ref_ratio[i][0] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        HeaderFile << my_geom[i].Domain() << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        HeaderFile << level_steps[i] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        for (int k = 0; k < AMREX_SPACEDIM; ++k) {
            HeaderFile << my_geom[i].CellSize()[k] << ' ';
        }
        HeaderFile << '\n';
    }
    HeaderFile << (int) my_geom[0].Coord() << '\n';
    HeaderFile << "0\n";

    for (int level = 0; level <= finest_level; ++level) {
        HeaderFile << level << ' ' << bArray[level].size() << ' ' << my_time << '\n';
        HeaderFile << level_steps[level] << '\n';

        const IntVect& domain_lo = my_geom[level].Domain().smallEnd();
        for (int i = 0; i < bArray[level].size(); ++i)
        {
            // Need to shift because the RealBox ctor we call takes the
            // physical location of index (0,0,0).  This does not affect
            // the usual cases where the domain index starts with 0.
            const Box& b = shift(bArray[level][i], -domain_lo);
            RealBox loc = RealBox(b, my_geom[level].CellSize(), my_geom[level].ProbLo());
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
