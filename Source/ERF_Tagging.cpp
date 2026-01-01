#include <ERF.H>
#include <ERF_Derive.H>

using namespace amrex;

/**
 * Function to tag cells for refinement -- this overrides the pure virtual function in AmrCore
 *
 * @param[in] levc level of refinement at which we tag cells (0 is coarsest level)
 * @param[out] tags array of tagged cells
 * @param[in] time current time
*/

#ifdef ERF_USE_NETCDF
Box read_subdomain_from_wrfinput (int lev, const std::string& fname, int& ratio);
Real read_start_time_from_wrfinput (int lev, const std::string& fname);
#endif

void
tag_on_distance_from_eye(const Geometry& cgeom, TagBoxArray* tags,
                         const Real eye_x, const Real eye_y, const Real rad_tag);

void
ERF::ErrorEst (int levc, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

#ifdef ERF_USE_NETCDF
    if (solverChoice.init_type == InitType::WRFInput) {
        int ratio;
        Box subdomain;

        if (!nc_init_file[levc+1].empty())
        {
            Real levc_start_time = read_start_time_from_wrfinput(levc  , nc_init_file[levc  ][0]);
            amrex::Print() << " WRFInput       time at level " << levc << " is " << levc_start_time << std::endl;

            for (int isub = 0; isub < nc_init_file[levc+1].size(); isub++) {
                if (!have_read_nc_init_file[levc+1][isub])
                {
                    Real levf_start_time = read_start_time_from_wrfinput(levc+1, nc_init_file[levc+1][isub]);
                    amrex::Print() << " WRFInput start_time at level " << levc+1 << " is " << levf_start_time << std::endl;

                    // We assume there is only one subdomain at levc; otherwise we don't know
                    //     which one is the parent of the fine region we are trying to create
                    AMREX_ALWAYS_ASSERT(subdomains[levc].size() == 1);

                    if ( (ref_ratio[levc][2]) != 1) {
                        amrex::Abort("The ref_ratio specified in the inputs file must have 1 in the z direction; please use ref_ratio_vect rather than ref_ratio");
                    }

                    if ( levf_start_time <= (levc_start_time + t_new[levc]) ) {
                        amrex::Print() << " WRFInput file to read: " << nc_init_file[levc+1][isub] << std::endl;
                        subdomain = read_subdomain_from_wrfinput(levc, nc_init_file[levc+1][isub], ratio);
                        amrex::Print() << " WRFInput subdomain " << isub << " at level " << levc+1 << " is " << subdomain << std::endl;

                        if ( (ratio != ref_ratio[levc][0]) || (ratio != ref_ratio[levc][1]) ) {
                            amrex::Print() << "File " << nc_init_file[levc+1][0] << " has refinement ratio = " << ratio << std::endl;
                            amrex::Print() << "The inputs file has refinement ratio = " << ref_ratio[levc] << std::endl;
                            amrex::Abort("These must be the same -- please edit your inputs file and try again.");
                        }

                        subdomain.coarsen(IntVect(ratio,ratio,1));

                        Box coarser_level(subdomains[levc][isub].minimalBox());
                        subdomain.shift(coarser_level.smallEnd());

                        if (verbose > 0) {
                            amrex::Print() << " Crse subdomain to be tagged is" << subdomain << std::endl;
                        }

                        Box new_fine(subdomain); new_fine.refine(IntVect(ratio,ratio,1));
                        num_boxes_at_level[levc+1] = 1;
                        boxes_at_level[levc+1].push_back(new_fine);

                        for (MFIter mfi(tags); mfi.isValid(); ++mfi) {
                            auto tag_arr = tags.array(mfi);  // Get device-accessible array

                            Box bx = mfi.validbox(); bx &= subdomain;

                            if (!bx.isEmpty()) {
                                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                    tag_arr(i,j,k) = TagBox::SET;
                                });
                            }
                        }
                    } // time is right
                } else {
                    // Re-tag this region
                    for (MFIter mfi(tags); mfi.isValid(); ++mfi)
                    {
                        auto tag_arr = tags.array(mfi);  // Get device-accessible array

                        Box existing_bx_coarsened(boxes_at_level[levc+1][isub]);
                        existing_bx_coarsened.coarsen(ref_ratio[levc]);

                        Box bx = mfi.validbox(); bx &= existing_bx_coarsened;

                        if (!bx.isEmpty()) {
                            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                tag_arr(i,j,k) = TagBox::SET;
                            });
                        }
                    }
                } // has file been read?
            } // isub
            return;
        } // file not empty
    }
#endif

    //
    // Make sure the ghost cells of the level we are tagging at are filled
    //    in case we take differences that require them
    // NOTE: We are Fillpatching only the cell-centered variables here
    //
    MultiFab& S_new = vars_new[levc][Vars::cons];
    MultiFab& U_new = vars_new[levc][Vars::xvel];
    MultiFab& V_new = vars_new[levc][Vars::yvel];
    MultiFab& W_new = vars_new[levc][Vars::zvel];
    //
    if (levc == 0) {
        FillPatchCrseLevel(levc, time, {&S_new, &U_new, &V_new, &W_new});
    } else {
        FillPatchFineLevel(levc, time, {&S_new, &U_new, &V_new, &W_new},
                           {&S_new, &rU_new[levc], &rV_new[levc], &rW_new[levc]},
                           base_state[levc], base_state[levc],
                           false, true);
    }

    for (int j=0; j < ref_tags.size(); ++j)
    {
        //
        // This mf must have ghost cells because we may take differences between adjacent values
        //
        std::unique_ptr<MultiFab> mf = std::make_unique<MultiFab>(grids[levc], dmap[levc], 1, 1);
        mf->setVal(0.0);

        // This allows dynamic refinement based on the value of the density
        if (ref_tags[j].Field() == "density")
        {
            MultiFab::Copy(*mf,vars_new[levc][Vars::cons],Rho_comp,0,1,1);

        // This allows dynamic refinement based on the value of qv
        } else if ( ref_tags[j].Field() == "qv" ) {
            MultiFab::Copy(  *mf, vars_new[levc][Vars::cons], RhoQ1_comp, 0, 1, 1);
            MultiFab::Divide(*mf, vars_new[levc][Vars::cons],   Rho_comp, 0, 1, 1);


        // This allows dynamic refinement based on the value of qc
        } else if (ref_tags[j].Field() == "qc" ) {
            MultiFab::Copy(  *mf, vars_new[levc][Vars::cons], RhoQ2_comp, 0, 1, 1);
            MultiFab::Divide(*mf, vars_new[levc][Vars::cons],   Rho_comp, 0, 1, 1);

        // This allows dynamic refinement based on the value of the z-component of vorticity
        } else if (ref_tags[j].Field() == "vorticity" ) {
            Vector<MultiFab> mf_cc_vel(1);
            mf_cc_vel[0].define(grids[levc], dmap[levc], AMREX_SPACEDIM, IntVect(1,1,1));
            average_face_to_cellcenter(mf_cc_vel[0],0,Array<const MultiFab*,3>{&U_new, &V_new, &W_new});

            // Impose bc's at domain boundaries at all levels
            FillBdyCCVels(mf_cc_vel,levc);

            mf->setVal(0.);

            for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                auto& dfab = (*mf)[mfi];
                auto& sfab = mf_cc_vel[0][mfi];
                derived::erf_dervortz(bx, dfab, 0, 1, sfab, Geom(levc), time, nullptr, levc);
            }

        // This allows dynamic refinement based on the value of the scalar/theta
        } else if ( (ref_tags[j].Field() == "scalar"  ) ||
                    (ref_tags[j].Field() == "theta"   ) )
        {
            for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox();
                auto& dfab = (*mf)[mfi];
                auto& sfab = vars_new[levc][Vars::cons][mfi];
                if (ref_tags[j].Field() == "scalar") {
                    derived::erf_derscalar(bx, dfab, 0, 1, sfab, Geom(levc), time, nullptr, levc);
                } else if (ref_tags[j].Field() == "theta") {
                    derived::erf_dertheta(bx, dfab, 0, 1, sfab, Geom(levc), time, nullptr, levc);
                }
            } // mfi
        // This allows dynamic refinement based on the value of the density
        } else if ( (SolverChoice::terrain_type == TerrainType::ImmersedForcing) &&
                    (ref_tags[j].Field() == "terrain_blanking") )
        {
            MultiFab::Copy(*mf,*terrain_blanking[levc],0,0,1,1);
        }
        else if (ref_tags[j].Field() == "velmag")
        {
            ParmParse pp(pp_prefix);
            Vector<std::string> refinement_indicators;
            pp.queryarr("refinement_indicators",refinement_indicators,0,pp.countval("refinement_indicators"));
            Real velmag_threshold;
            bool is_hurricane_tracker = false;
            for (int i=0; i<refinement_indicators.size(); ++i)
            {
                if (refinement_indicators[i]=="hurricane_tracker") {
                    is_hurricane_tracker = true;
                    std::string ref_prefix = pp_prefix + "." + refinement_indicators[i];
                    ParmParse ppr(ref_prefix);
                    ppr.get("value_greater", velmag_threshold);
                    break;
                }
            }

            Vector<MultiFab> mf_cc_vel(1);
            mf_cc_vel[0].define(grids[levc], dmap[levc], AMREX_SPACEDIM, IntVect(0,0,0));
            average_face_to_cellcenter(mf_cc_vel[0],0,Array<const MultiFab*,3>{&U_new, &V_new, &W_new});

            if (is_hurricane_tracker) {
                HurricaneTracker(levc, time, mf_cc_vel[0], velmag_threshold, &tags);
            } else {
                for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    auto& dfab = (*mf)[mfi];
                    auto& sfab = mf_cc_vel[0][mfi];
                    derived::erf_dermagvel(bx, dfab, 0, 1, sfab, Geom(levc), time, nullptr, levc);
                }
            }

#ifdef ERF_USE_PARTICLES
        } else {
            //
            // This allows dynamic refinement based on the number of particles per cell
            //
            // Note that we must count all the particles in levels both at and above the current,
            //      since otherwise, e.g., if the particles are all at level 1, counting particles at
            //      level 0 will not trigger refinement when regridding so level 1 will disappear,
            //      then come back at the next regridding
            //
            const auto& particles_namelist( particleData.getNames() );
            mf->setVal(0.0);
            for (ParticlesNamesVector::size_type i = 0; i < particles_namelist.size(); i++)
            {
                std::string tmp_string(particles_namelist[i]+"_count");
                IntVect rr = IntVect::TheUnitVector();
                if (ref_tags[j].Field() == tmp_string) {
                    for (int lev = levc; lev <= finest_level; lev++)
                    {
                        MultiFab temp_dat(grids[lev], dmap[lev], 1, 0); temp_dat.setVal(0);
                        particleData[particles_namelist[i]]->IncrementWithTotal(temp_dat, lev);

                        MultiFab temp_dat_crse(grids[levc], dmap[levc], 1, 0); temp_dat_crse.setVal(0);

                        if (lev == levc) {
                            MultiFab::Copy(*mf, temp_dat, 0, 0, 1, 0);
                        } else {
                            for (int d = 0; d < AMREX_SPACEDIM; d++) {
                                rr[d] *= ref_ratio[levc][d];
                            }
                            average_down(temp_dat, temp_dat_crse, 0, 1, rr);
                            MultiFab::Add(*mf, temp_dat_crse, 0, 0, 1, 0);
                        }
                    }
                }
            }
#endif
        }

        ref_tags[j](tags,mf.get(),clearval,tagval,time,levc,geom[levc]);
    } // loop over j

    // ********************************************************************************************
    // Refinement based on 2d distance from the "eye" which is defined here as the (x,y) location of
    //    the integrated qv
    // ********************************************************************************************
    ParmParse pp(pp_prefix);
    Vector<std::string> refinement_indicators;
    pp.queryarr("refinement_indicators",refinement_indicators,0,pp.countval("refinement_indicators"));
    for (int i=0; i<refinement_indicators.size(); ++i)
    {
        if ( (refinement_indicators[i]=="storm_tracker") && (solverChoice.moisture_type != MoistureType::None) )
        {
            std::string ref_prefix = pp_prefix + "." + refinement_indicators[i];
            ParmParse ppr(ref_prefix);

            Real ref_start_time = -1.0;
            ppr.query("start_time",ref_start_time);

            if (time >= ref_start_time) {

                Real max_radius = -1.0;
                ppr.get("max_radius", max_radius);

                // Create the volume-weighted sum of (rho qv) in each column
                MultiFab mf_qv_int(ba2d[levc], dmap[levc], 1, 0); mf_qv_int.setVal(0.);

                // Define the 2D MultiFab holding the column-integrated (rho qv)
                volWgtColumnSum(levc, S_new, RhoQ1_comp, mf_qv_int, *detJ_cc[levc]);

                // Find the max value in the domain
                IntVect eye = mf_qv_int.maxIndex(0);

                const auto dx      = geom[levc].CellSizeArray();
                const auto prob_lo = geom[levc].ProbLoArray();

                Real eye_x = prob_lo[0] + (eye[0] + 0.5) * dx[0];
                Real eye_y = prob_lo[1] + (eye[1] + 0.5) * dx[1];

                tag_on_distance_from_eye(geom[levc], &tags, eye_x, eye_y, max_radius);
            }
        }
    }
}

/**
 * Function to define the refinement criteria based on user input
*/

void
ERF::refinement_criteria_setup ()
{
    if (max_level > 0)
    {
        ParmParse pp(pp_prefix);
        Vector<std::string> refinement_indicators;
        pp.queryarr("refinement_indicators",refinement_indicators,0,pp.countval("refinement_indicators"));

        for (int i=0; i<refinement_indicators.size(); ++i)
        {
            std::string ref_prefix = pp_prefix + "." + refinement_indicators[i];

            ParmParse ppr(ref_prefix);
            RealBox realbox;
            int lev_for_box;

            int num_real_lo = ppr.countval("in_box_lo");
            int num_indx_lo = ppr.countval("in_box_lo_indices");
            int num_real_hi = ppr.countval("in_box_hi");
            int num_indx_hi = ppr.countval("in_box_hi_indices");

            AMREX_ALWAYS_ASSERT(num_real_lo == num_real_hi);
            AMREX_ALWAYS_ASSERT(num_indx_lo == num_indx_hi);

            if ( !((num_real_lo >= AMREX_SPACEDIM-1 && num_indx_lo == 0) ||
                   (num_indx_lo >= AMREX_SPACEDIM-1 && num_real_lo == 0) ||
                   (num_indx_lo ==              0   && num_real_lo == 0)) )
            {
                amrex::Abort("Must only specify box for refinement using real OR index space");
            }

            if (num_real_lo > 0) {
                std::vector<Real> rbox_lo(3), rbox_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box <= max_level)
                {
                    if (n_error_buf[0] != IntVect::TheZeroVector()) {
                        amrex::Abort("Don't use n_error_buf > 0 when setting the box explicitly");
                    }

                    const Real* plo = geom[lev_for_box].ProbLo();
                    const Real* phi = geom[lev_for_box].ProbHi();

                    ppr.getarr("in_box_lo",rbox_lo,0,num_real_lo);
                    ppr.getarr("in_box_hi",rbox_hi,0,num_real_hi);

                    if (rbox_lo[0] < plo[0]) rbox_lo[0] = plo[0];
                    if (rbox_lo[1] < plo[1]) rbox_lo[1] = plo[1];
                    if (rbox_hi[0] > phi[0]) rbox_hi[0] = phi[0];
                    if (rbox_hi[1] > phi[1]) rbox_hi[1] = phi[1];
                    if (num_real_lo < AMREX_SPACEDIM) {
                        rbox_lo[2] = plo[2];
                        rbox_hi[2] = phi[2];
                    }

                    realbox = RealBox(&(rbox_lo[0]),&(rbox_hi[0]));

                    Print() << "Realbox read in and intersected laterally with domain is " << realbox << std::endl;

                    num_boxes_at_level[lev_for_box] += 1;

                    int ilo, jlo, klo;
                    int ihi, jhi, khi;
                    const auto* dx  = geom[lev_for_box].CellSize();
                    ilo = static_cast<int>((rbox_lo[0] - plo[0])/dx[0]);
                    jlo = static_cast<int>((rbox_lo[1] - plo[1])/dx[1]);
                    ihi = static_cast<int>((rbox_hi[0] - plo[0])/dx[0]-1);
                    jhi = static_cast<int>((rbox_hi[1] - plo[1])/dx[1]-1);
                    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
                        // Search for k indices corresponding to nominal grid
                        // AGL heights
                        const Box& domain = geom[lev_for_box].Domain();
                        klo = domain.smallEnd(2) - 1;
                        khi = domain.smallEnd(2) - 1;

                        if (rbox_lo[2] <= zlevels_stag[lev_for_box][domain.smallEnd(2)])
                        {
                            klo = domain.smallEnd(2);
                        }
                        else
                        {
                            for (int k=domain.smallEnd(2); k<=domain.bigEnd(2)+1; ++k) {
                                if (zlevels_stag[lev_for_box][k] > rbox_lo[2]) {
                                    klo = k-1;
                                    break;
                                }
                            }
                        }
                        AMREX_ASSERT(klo >= domain.smallEnd(2));

                        if (rbox_hi[2] >= zlevels_stag[lev_for_box][domain.bigEnd(2)+1])
                        {
                            khi = domain.bigEnd(2);
                        }
                        else
                        {
                            for (int k=klo+1; k<=domain.bigEnd(2)+1; ++k) {
                                if (zlevels_stag[lev_for_box][k] > rbox_hi[2]) {
                                    khi = k-1;
                                    break;
                                }
                            }
                        }
                        AMREX_ASSERT((khi <= domain.bigEnd(2)) && (khi > klo));

                        // Need to update realbox because tagging is based on
                        // the initial _un_deformed grid
                        realbox = RealBox(plo[0]+ ilo   *dx[0], plo[1]+ jlo   *dx[1], plo[2]+ klo   *dx[2],
                                          plo[0]+(ihi+1)*dx[0], plo[1]+(jhi+1)*dx[1], plo[2]+(khi+1)*dx[2]);
                    } else {
                        klo = static_cast<int>((rbox_lo[2] - plo[2])/dx[2]);
                        khi = static_cast<int>((rbox_hi[2] - plo[2])/dx[2]-1);
                    }

                    Box bx(IntVect(ilo,jlo,klo),IntVect(ihi,jhi,khi));
                    if ( (ilo%ref_ratio[lev_for_box-1][0] != 0) || ((ihi+1)%ref_ratio[lev_for_box-1][0] != 0) ||
                         (jlo%ref_ratio[lev_for_box-1][1] != 0) || ((jhi+1)%ref_ratio[lev_for_box-1][1] != 0) ||
                         (klo%ref_ratio[lev_for_box-1][2] != 0) || ((khi+1)%ref_ratio[lev_for_box-1][2] != 0) )
                    {
                        amrex::Print() << "Box : " << bx << std::endl;
                        amrex::Print() << "RealBox : " << realbox << std::endl;
                        amrex::Print() << "ilo, ihi+1, jlo, jhi+1, klo, khi+1 by ref_ratio : "
                                       << ilo%ref_ratio[lev_for_box-1][0] << " " << (ihi+1)%ref_ratio[lev_for_box-1][0] << " "
                                       << jlo%ref_ratio[lev_for_box-1][1] << " " << (jhi+1)%ref_ratio[lev_for_box-1][1] << " "
                                       << klo%ref_ratio[lev_for_box-1][2] << " " << (khi+1)%ref_ratio[lev_for_box-1][2] << std::endl;
                        amrex::Error("Fine box is not legit with this ref_ratio");
                    }
                    boxes_at_level[lev_for_box].push_back(bx);
                    Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev

                if (solverChoice.init_type == InitType::WRFInput) {
                    if (num_boxes_at_level[lev_for_box] != num_files_at_level[lev_for_box]) {
                        amrex::Error("Number of boxes doesn't match number of input files");

                    }
                }

            } else if (num_indx_lo > 0) {

                std::vector<int> box_lo(3), box_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box <= max_level)
                {
                    if (n_error_buf[0] != IntVect::TheZeroVector()) {
                        amrex::Abort("Don't use n_error_buf > 0 when setting the box explicitly");
                    }

                    ppr.getarr("in_box_lo_indices",box_lo,0,AMREX_SPACEDIM);
                    ppr.getarr("in_box_hi_indices",box_hi,0,AMREX_SPACEDIM);

                    Box bx(IntVect(box_lo[0],box_lo[1],box_lo[2]),IntVect(box_hi[0],box_hi[1],box_hi[2]));
                    amrex::Print() << "BOX " << bx << std::endl;

                    const auto* dx  = geom[lev_for_box].CellSize();
                    const Real* plo = geom[lev_for_box].ProbLo();
                    realbox = RealBox(plo[0]+ box_lo[0]   *dx[0], plo[1]+ box_lo[1]   *dx[1], plo[2]+ box_lo[2]   *dx[2],
                                      plo[0]+(box_hi[0]+1)*dx[0], plo[1]+(box_hi[1]+1)*dx[1], plo[2]+(box_hi[2]+1)*dx[2]);

                    Print() << "Reading " << bx << " at level " << lev_for_box << std::endl;
                    num_boxes_at_level[lev_for_box] += 1;

                    if ( (box_lo[0]%ref_ratio[lev_for_box-1][0] != 0) || ((box_hi[0]+1)%ref_ratio[lev_for_box-1][0] != 0) ||
                         (box_lo[1]%ref_ratio[lev_for_box-1][1] != 0) || ((box_hi[1]+1)%ref_ratio[lev_for_box-1][1] != 0) ||
                         (box_lo[2]%ref_ratio[lev_for_box-1][2] != 0) || ((box_hi[2]+1)%ref_ratio[lev_for_box-1][2] != 0) )
                         amrex::Error("Fine box is not legit with this ref_ratio");
                    boxes_at_level[lev_for_box].push_back(bx);
                    Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev

                if (solverChoice.init_type == InitType::WRFInput) {
                    if (num_boxes_at_level[lev_for_box] != num_files_at_level[lev_for_box]) {
                        amrex::Error("Number of boxes doesn't match number of input files");

                    }
                }
            }

            AMRErrorTagInfo info;

            if (realbox.ok()) {
                info.SetRealBox(realbox);
            }
            if (ppr.countval("start_time") > 0) {
                Real ref_min_time; ppr.get("start_time",ref_min_time);
                info.SetMinTime(ref_min_time);
            }
            if (ppr.countval("end_time") > 0) {
                Real ref_max_time; ppr.get("end_time",ref_max_time);
                info.SetMaxTime(ref_max_time);
            }
            if (ppr.countval("max_level") > 0) {
                int ref_max_level; ppr.get("max_level",ref_max_level);
                info.SetMaxLevel(ref_max_level);
            }

            if (ppr.countval("value_greater")) {
                int num_val = ppr.countval("value_greater");
                Vector<Real> value(num_val);
                ppr.getarr("value_greater",value,0,num_val);
                std::string field; ppr.get("field_name",field);
                ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::GREATER,field,info));
            }
            else if (ppr.countval("value_less")) {
                int num_val = ppr.countval("value_less");
                Vector<Real> value(num_val);
                ppr.getarr("value_less",value,0,num_val);
                std::string field; ppr.get("field_name",field);
                ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,field,info));
            }
            else if (ppr.countval("adjacent_difference_greater")) {
                int num_val = ppr.countval("adjacent_difference_greater");
                Vector<Real> value(num_val);
                ppr.getarr("adjacent_difference_greater",value,0,num_val);
                std::string field; ppr.get("field_name",field);
                ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,field,info));
            }
            else if (realbox.ok())
            {
                ref_tags.push_back(AMRErrorTag(info));
            } else if (refinement_indicators[i] != "storm_tracker") {
                Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
            }
        } // loop over criteria
    } // if max_level > 0
}

bool
ERF::FindInitialEye(int levc,
                    const MultiFab& mf_cc_vel,
                    const Real velmag_threshold,
                    Real& eye_x, Real& eye_y)
{
    const auto dx = geom[levc].CellSizeArray();
    const auto prob_lo = geom[levc].ProbLoArray();

    Gpu::DeviceVector<Real> d_coords(2, 0.0);
    Gpu::DeviceVector<int>  d_found(1,0);

    Real* d_coords_ptr = d_coords.data();
    int*   d_found_ptr = d_found.data();

    for (MFIter mfi(mf_cc_vel); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const Array4<const Real>& vel_arr = mf_cc_vel.const_array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real magnitude = std::sqrt(vel_arr(i,j,k,0) * vel_arr(i,j,k,0) +
                                       vel_arr(i,j,k,1) * vel_arr(i,j,k,1) +
                                       vel_arr(i,j,k,2) * vel_arr(i,j,k,2));

            magnitude *= 3.6;

            Real z = prob_lo[2] + (k + 0.5) * dx[2];

            // Check if magnitude exceeds threshold
            if (z < 2000. && magnitude > velmag_threshold) {
                // Use atomic operations to set found flag and store coordinates
                Gpu::Atomic::Add(&d_found_ptr[0], 1); // Mark as found

                Real x = prob_lo[0] + (i + 0.5) * dx[0];
                Real y = prob_lo[1] + (j + 0.5) * dx[1];

                // Store coordinates
                Gpu::Atomic::Add(&d_coords_ptr[0],x); // Store x index
                Gpu::Atomic::Add(&d_coords_ptr[1],y); // Store x index
            }
        });
    }

    // Synchronize to ensure all threads complete their execution
    amrex::Gpu::streamSynchronize(); // Wait for all GPU threads to finish

    Vector<int> h_found(1,0);
    Gpu::copy(Gpu::deviceToHost, d_found.begin(), d_found.end(), h_found.begin());
    ParallelAllReduce::Sum(h_found.data(), h_found.size(), ParallelContext::CommunicatorAll());

    // Broadcast coordinates if found
    if (h_found[0] > 0) {
        Vector<Real> h_coords(2,-1e10);
        Gpu::copy(Gpu::deviceToHost, d_coords.begin(), d_coords.end(), h_coords.begin());

        ParallelAllReduce::Sum(h_coords.data(), h_coords.size(), ParallelContext::CommunicatorAll());

        eye_x = h_coords[0]/h_found[0];
        eye_y = h_coords[1]/h_found[0];

    } else {
        // Random large negative numbers so we don't trigger refinement in this case
        eye_x = -1.e20;
        eye_y = -1.e20;
    }

    return (h_found[0] > 0);
}

void
tag_on_distance_from_eye(const Geometry& cgeom, TagBoxArray* tags,
                         const Real eye_x, const Real eye_y, const Real rad_tag)
{
    const auto dx      = cgeom.CellSizeArray();
    const auto prob_lo = cgeom.ProbLoArray();

    for (MFIter mfi(*tags); mfi.isValid(); ++mfi) {
        TagBox& tag = (*tags)[mfi];
        auto tag_arr = tag.array();  // Get device-accessible array

        const Box& tile_box = mfi.tilebox(); // The box for this tile

        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Compute cell center coordinates
            Real x = prob_lo[0] + (i + 0.5) * dx[0];
            Real y = prob_lo[1] + (j + 0.5) * dx[1];

            Real dist = std::sqrt((x - eye_x)*(x - eye_x) + (y - eye_y)*(y - eye_y));

            if (dist < rad_tag) {
                tag_arr(i,j,k) = TagBox::SET;
            } else {
                tag_arr(i,j,k) = TagBox::CLEAR;
            }
        });
    }
}

void
ERF::HurricaneTracker(int levc,
                      Real time,
                      const MultiFab& mf_cc_vel,
                      const Real velmag_threshold,
                      TagBoxArray* tags)
{
    bool is_found;

    Real eye_x, eye_y;

    if (time==0.0) {
        is_found = FindInitialEye(levc, mf_cc_vel, velmag_threshold, eye_x, eye_y);
    } else {
        is_found = true;
        const auto& last = hurricane_eye_track_xy.back();
        eye_x = last[0];
        eye_y = last[1];
    }

    if (is_found) {
        Real rad_tag = 4.e5 * std::pow(2, max_level-1-levc);
        tag_on_distance_from_eye(geom[levc], tags, eye_x, eye_y, rad_tag);
    }
}
