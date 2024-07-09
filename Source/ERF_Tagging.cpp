#include <ERF.H>
#include <Derive.H>

using namespace amrex;

/**
 * Function to tag cells for refinement -- this overrides the pure virtual function in AmrCore
 *
 * @param[in] levc level of refinement at which we tag cells (0 is coarsest level)
 * @param[out] tags array of tagged cells
 * @param[in] time current time
*/

void
ERF::ErrorEst (int levc, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    for (int j=0; j < ref_tags.size(); ++j)
    {
        std::unique_ptr<MultiFab> mf = std::make_unique<MultiFab>(grids[levc], dmap[levc], 1, 0);

        // This allows dynamic refinement based on the value of the density
        if (ref_tags[j].Field() == "density")
        {
            MultiFab::Copy(*mf,vars_new[levc][Vars::cons],Rho_comp,0,1,0);

        // This allows dynamic refinement based on the value of qv
        } else if ( ref_tags[j].Field() == "qv" ) {
            MultiFab qv(vars_new[levc][Vars::cons],make_alias,0,RhoQ1_comp+1);
            MultiFab::Copy(  *mf, qv, RhoQ1_comp, 0, 1, 0);
            MultiFab::Divide(*mf, qv, Rho_comp  , 0, 1, 0);


        // This allows dynamic refinement based on the value of qc
        } else if (ref_tags[j].Field() == "qc" ) {
            MultiFab qc(vars_new[levc][Vars::cons],make_alias,0,RhoQ2_comp+1);
            MultiFab::Copy(  *mf, qc, RhoQ2_comp, 0, 1, 0);
            MultiFab::Divide(*mf, qc, Rho_comp  , 0, 1, 0);


        // This allows dynamic refinement based on the value of the scalar/pressure/theta
        } else if ( (ref_tags[j].Field() == "scalar"  ) ||
                    (ref_tags[j].Field() == "pressure") ||
                    (ref_tags[j].Field() == "theta"   ) )
        {
            for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                auto& dfab = (*mf)[mfi];
                auto& sfab = vars_new[levc][Vars::cons][mfi];
                if (ref_tags[j].Field() == "scalar") {
                    derived::erf_derscalar(bx, dfab, 0, 1, sfab, Geom(levc), time, nullptr, levc);
                } else if (ref_tags[j].Field() == "theta") {
                    derived::erf_dertheta(bx, dfab, 0, 1, sfab, Geom(levc), time, nullptr, levc);
                }
            } // mfi
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

            if ( !((num_real_lo == AMREX_SPACEDIM && num_indx_lo == 0) ||
                   (num_indx_lo == AMREX_SPACEDIM && num_real_lo == 0) ||
                   (num_indx_lo ==              0 && num_real_lo == 0)) )
            {
                amrex::Abort("Must only specify box for refinement using real OR index space");
            }

            if (num_real_lo > 0) {
                std::vector<Real> box_lo(3), box_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box <= max_level)
                {
                    if (n_error_buf[0] != IntVect::TheZeroVector()) {
                        amrex::Abort("Don't use n_error_buf > 0 when setting the box explicitly");
                    }

                    ppr.getarr("in_box_lo",box_lo,0,AMREX_SPACEDIM);
                    ppr.getarr("in_box_hi",box_hi,0,AMREX_SPACEDIM);
                    realbox = RealBox(&(box_lo[0]),&(box_hi[0]));

                    Print() << "Reading " << realbox << " at level " << lev_for_box << std::endl;
                    num_boxes_at_level[lev_for_box] += 1;

                    const auto* dx  = geom[lev_for_box].CellSize();
                    const Real* plo = geom[lev_for_box].ProbLo();
                    int ilo = static_cast<int>((box_lo[0] - plo[0])/dx[0]);
                    int jlo = static_cast<int>((box_lo[1] - plo[1])/dx[1]);
                    int klo = static_cast<int>((box_lo[2] - plo[2])/dx[2]);
                    int ihi = static_cast<int>((box_hi[0] - plo[0])/dx[0]-1);
                    int jhi = static_cast<int>((box_hi[1] - plo[1])/dx[1]-1);
                    int khi = static_cast<int>((box_hi[2] - plo[2])/dx[2]-1);
                    Box bx(IntVect(ilo,jlo,klo),IntVect(ihi,jhi,khi));
                    if ( (ilo%ref_ratio[lev_for_box-1][0] != 0) || ((ihi+1)%ref_ratio[lev_for_box-1][0] != 0) ||
                         (jlo%ref_ratio[lev_for_box-1][1] != 0) || ((jhi+1)%ref_ratio[lev_for_box-1][1] != 0) ||
                         (klo%ref_ratio[lev_for_box-1][2] != 0) || ((khi+1)%ref_ratio[lev_for_box-1][2] != 0) )
                         amrex::Error("Fine box is not legit with this ref_ratio");
                    boxes_at_level[lev_for_box].push_back(bx);
                    Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev
                if (init_type == "real" || init_type == "metgrid") {
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
                    realbox = RealBox(plo[0]+ box_lo[0]   *dx[0],plo[1] +box_lo[1]   *dx[1],plo[2] +box_lo[2]   *dx[2],
                                      plo[0]+(box_hi[0]+1)*dx[0],plo[1]+(box_hi[1]+1)*dx[1],plo[2]+(box_hi[2]+1)*dx[2]);

                    Print() << "Reading " << bx << " at level " << lev_for_box << std::endl;
                    num_boxes_at_level[lev_for_box] += 1;

                    if ( (box_lo[0]%ref_ratio[lev_for_box-1][0] != 0) || ((box_hi[0]+1)%ref_ratio[lev_for_box-1][0] != 0) ||
                         (box_lo[1]%ref_ratio[lev_for_box-1][1] != 0) || ((box_hi[1]+1)%ref_ratio[lev_for_box-1][1] != 0) ||
                         (box_lo[2]%ref_ratio[lev_for_box-1][2] != 0) || ((box_hi[2]+1)%ref_ratio[lev_for_box-1][2] != 0) )
                         amrex::Error("Fine box is not legit with this ref_ratio");
                    boxes_at_level[lev_for_box].push_back(bx);
                    Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev
                if (init_type == "real" || init_type == "metgrid") {
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
            } else {
                Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
            }
        } // loop over criteria
    } // if max_level > 0
}
