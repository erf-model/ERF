#include <ERF_EOS.H>
#include <ERF.H>
#include <ERF_EpochTime.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
ERF::setSubVolVariables (const std::string& pp_subvol_var_names,
                         Vector<std::string>&  subvol_var_names)
{
    ParmParse pp(pp_prefix);

    std::string nm;

    int nSubVolVars = pp.countval(pp_subvol_var_names.c_str());

    // We pre-populate the list with velocities, but allow these to be over-written
    //    by user input
    if (nSubVolVars == 0)
    {
        subvol_var_names.push_back("x_velocity");
        subvol_var_names.push_back("y_velocity");
        subvol_var_names.push_back("z_velocity");

    } else {
        for (int i = 0; i < nSubVolVars; i++)
        {
            pp.get(pp_subvol_var_names.c_str(), nm, i);

            // Add the named variable to our list of subvol variables
            // if it is not already in the list
            if (!containerHasElement(subvol_var_names, nm)) {
                subvol_var_names.push_back(nm);
            }
        }
    }

    // Get state variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_plot_names;

    for (int i = 0; i < cons_names.size(); ++i) {
        if ( containerHasElement(subvol_var_names, cons_names[i]) ) {
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

    // Check for velocity since it's not in cons_names
    if (containerHasElement(subvol_var_names, "x_velocity")) {
        tmp_plot_names.push_back("x_velocity");
    }
    if (containerHasElement(subvol_var_names, "y_velocity")) {
        tmp_plot_names.push_back("y_velocity");
    }
    if (containerHasElement(subvol_var_names, "z_velocity")) {
        tmp_plot_names.push_back("z_velocity");
    }

    //
    // If the model we are running doesn't have the variable listed in the inputs file,
    //     just ignore it rather than aborting
    //
    for (int i = 0; i < derived_subvol_names.size(); ++i) {
        if ( containerHasElement(subvol_var_names, derived_names[i]) ) {
            bool ok_to_add = ( (solverChoice.terrain_type == TerrainType::ImmersedForcing) ||
                               (derived_names[i] != "terrain_IB_mask") );
            ok_to_add     &= ( (SolverChoice::terrain_type == TerrainType::StaticFittedMesh) ||
                               (SolverChoice::terrain_type == TerrainType::MovingFittedMesh) ||
                               (derived_names[i] != "detJ") );
            ok_to_add     &= ( (SolverChoice::terrain_type == TerrainType::StaticFittedMesh) ||
                               (SolverChoice::terrain_type == TerrainType::MovingFittedMesh) ||
                               (derived_names[i] != "z_phys") );
            if (ok_to_add)
            {
                if (solverChoice.moisture_type == MoistureType::None) { // no moist quantities allowed
                    if (derived_names[i] != "qv" && derived_names[i] != "qc"    && derived_names[i] != "qrain"  &&
                        derived_names[i] != "qi" && derived_names[i] != "qsnow" && derived_names[i] != "qgraup" &&
                        derived_names[i] != "qt" && derived_names[i] != "qn"    && derived_names[i] != "qp" &&
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
                        derived_names[i] != "qp" &&
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

    subvol_var_names = tmp_plot_names;
}

// Write plotfile to disk
void
ERF::WriteSubvolume (int isub,Vector<std::string> subvol_var_names)
{
    ParmParse pp("erf.subvol");

    Vector<Real> origin;
    Vector< int> ncell;
    Vector<Real> delta;

    // **************************************************************
    // Read in the origin, number of cells in each dir, and resolution
    // **************************************************************

    int lev_for_sub = 0;
    int offset = isub * AMREX_SPACEDIM;

    pp.getarr("origin",origin,offset,AMREX_SPACEDIM);
    pp.getarr("nxnynz", ncell,offset,AMREX_SPACEDIM);
    pp.getarr("dxdydz", delta,offset,AMREX_SPACEDIM);

    bool found = false;
    for (int i = 0; i <= finest_level; i++) {
        if (!found) {
            if (almostEqual(delta[offset+0],geom[i].CellSize(0)) &&
                almostEqual(delta[offset+1],geom[i].CellSize(1)) &&
                almostEqual(delta[offset+2],geom[i].CellSize(2)) ) {

                amrex::Print() << "WriteSubvolume:Resolution specified matches that of level " << i << std::endl;
                found = true;
                lev_for_sub = i;
            }
        }
    }

    if (!found) {
        amrex::Abort("Resolution specified for subvol does not match the resolution of any of the levels.");
    }


    // **************************************************************
    // Now that we know which level we're at, we can figure out which (i,j,k) the origin corresponds to
    // Note we use 1.0001 as a fudge factor since the division of two reals --> integer will do a floor
    // **************************************************************
    int i0 = static_cast<int>((origin[offset+0] - geom[lev_for_sub].ProbLo(0)) * 1.0001 / delta[offset+0]);
    int j0 = static_cast<int>((origin[offset+1] - geom[lev_for_sub].ProbLo(1)) * 1.0001 / delta[offset+1]);
    int k0 = static_cast<int>((origin[offset+2] - geom[lev_for_sub].ProbLo(2)) * 1.0001 / delta[offset+2]);

    found = false;
    if (almostEqual(geom[lev_for_sub].ProbLo(0)+i0*delta[offset+0],origin[offset+0]) &&
        almostEqual(geom[lev_for_sub].ProbLo(1)+j0*delta[offset+1],origin[offset+1]) &&
        almostEqual(geom[lev_for_sub].ProbLo(2)+k0*delta[offset+2],origin[offset+2]) )
    {
        amrex::Print() << "WriteSubvolume:Specified origin is the lower left corner of cell " << IntVect(i0,j0,k0) << std::endl;
        found = true;
    }

    if (!found) {
        amrex::Abort("Origin specified does not correspond to a node at this level.");
    }

    Box domain(geom[lev_for_sub].Domain());

    Box bx(IntVect(i0,j0,k0),IntVect(i0+ncell[offset+0]-1,j0+ncell[offset+1]-1,k0+ncell[offset+2]-1));
    amrex::Print() << "WriteSubvolume:Box requested is " << bx << std::endl;

    if (!domain.contains(bx))
    {
        amrex::Abort("WriteSubvolume:Box requested is larger than the existing domain");
    }

    Vector<int> cs(AMREX_SPACEDIM);
    int count = pp.countval("chunk_size");
    if (count > 0) {
        pp.queryarr("chunk_size",cs,0,AMREX_SPACEDIM);
    } else {
        cs[0] = max_grid_size[0][0];
        cs[1] = max_grid_size[0][1];
        cs[2] = max_grid_size[0][2];
    }
    IntVect chunk_size(cs[0],cs[1],cs[2]);

    BoxArray ba(bx);
    ba.maxSize(chunk_size);

    amrex::Print() << "WriteSubvolume:BoxArray is " << ba << std::endl;

    Vector<std::string> varnames;
    varnames.insert(varnames.end(), subvol_var_names.begin(), subvol_var_names.end());

    int ncomp_mf = subvol_var_names.size();

    DistributionMapping dm(ba);

    MultiFab mf(ba, dm, ncomp_mf, 0);

    int mf_comp = 0;

    // *****************************************************************************************

    // First, copy any of the conserved state variables into the output plotfile
    for (int i = 0; i < cons_names.size(); ++i) {
        if (containerHasElement(subvol_var_names, cons_names[i])) {
            mf.ParallelCopy(vars_new[lev_for_sub][Vars::cons],i,mf_comp,1,1,0);
            mf_comp++;
        }
    }

    // *****************************************************************************************

    if (containerHasElement(subvol_var_names, "x_velocity") ||
        containerHasElement(subvol_var_names, "y_velocity") ||
        containerHasElement(subvol_var_names, "z_velocity"))
    {
        MultiFab mf_cc_vel(grids[lev_for_sub], dmap[lev_for_sub], AMREX_SPACEDIM, 0);
        average_face_to_cellcenter(mf_cc_vel,0,
                                   Array<const MultiFab*,3>{&vars_new[lev_for_sub][Vars::xvel],
                                                            &vars_new[lev_for_sub][Vars::yvel],
                                                            &vars_new[lev_for_sub][Vars::zvel]});
        if (containerHasElement(subvol_var_names, "x_velocity")) {
            mf.ParallelCopy(mf_cc_vel,0,mf_comp,1,0,0);
            mf_comp++;
        }
        if (containerHasElement(subvol_var_names, "y_velocity")) {
            mf.ParallelCopy(mf_cc_vel,1,mf_comp,1,0,0);
            mf_comp++;
        }
        if (containerHasElement(subvol_var_names, "z_velocity")) {
            mf.ParallelCopy(mf_cc_vel,2,mf_comp,1,0,0);
            mf_comp++;
        }
    }

    // *****************************************************************************************

    // Finally, check for any derived quantities and compute them, inserting
    // them into our output multifab
    auto calculate_derived = [&](const std::string& der_name,
                                 MultiFab& src_mf,
                                 decltype(derived::erf_dernull)& der_function)
    {
        if (containerHasElement(subvol_var_names, der_name)) {
            MultiFab dmf(src_mf.boxArray(), src_mf.DistributionMap(), 1, 0);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.tilebox();
                auto& dfab = dmf[mfi];
                auto& sfab = src_mf[mfi];
                der_function(tbx, dfab, 0, 1, sfab, Geom(lev_for_sub), t_new[0], nullptr, lev_for_sub);
            }
            mf.ParallelCopy(dmf,0,mf_comp,1,0,0);
            mf_comp++;
        }
    };

    // *****************************************************************************************
    // NOTE: All derived variables computed below **MUST MATCH THE ORDER** of "derived_names"
    //       defined in ERF.H
    // *****************************************************************************************

    calculate_derived("soundspeed",  vars_new[lev_for_sub][Vars::cons], derived::erf_dersoundspeed);
    if (solverChoice.moisture_type != MoistureType::None) {
        calculate_derived("temp",    vars_new[lev_for_sub][Vars::cons], derived::erf_dermoisttemp);
    } else {
        calculate_derived("temp",    vars_new[lev_for_sub][Vars::cons], derived::erf_dertemp);
    }
    calculate_derived("theta",       vars_new[lev_for_sub][Vars::cons], derived::erf_dertheta);
    calculate_derived("KE",          vars_new[lev_for_sub][Vars::cons], derived::erf_derKE);
    calculate_derived("scalar",      vars_new[lev_for_sub][Vars::cons], derived::erf_derscalar);

    // *****************************************************************************************

    Real time = t_new[lev_for_sub];

    std::string sf = subvol_file + "_" + std::to_string(isub);
    std::string subvol_filename;

    if (use_real_time_in_pltname) {
        const std::string dt_format = "%Y-%m-%d_%H:%M:%S"; // ISO 8601 standard
        subvol_filename = sf + getTimestamp(start_time+time, dt_format);
    } else {
       subvol_filename = Concatenate(sf + "_", istep[0], file_name_digits);
    }

    amrex::Print() <<"Writing subvolume into " << subvol_filename << std::endl;
    WriteSingleLevelPlotfile(subvol_filename,mf,varnames,geom[lev_for_sub],time,istep[0]);

}
