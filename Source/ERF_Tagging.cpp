#include <ERF.H>
#include <Derive.H>

using namespace amrex;

/**
 * Function to tag cells for refinement -- this overrides the pure virtual function in AmrCore
 *
 * @param[in] level level of refinement (0 is coarsest leve)
 * @param[out] tags array of tagged cells
 * @param[in] time current time
*/

void
ERF::ErrorEst (int level, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;
    for (int j=0; j < ref_tags.size(); ++j)
    {
        std::unique_ptr<MultiFab> mf;
        bool tag_found = false;

        // This allows dynamic refinement based on the value of the density
        if (ref_tags[j].Field() == "density")
        {
            mf = std::make_unique<MultiFab>(grids[level], dmap[level], 1, 0);
            MultiFab::Copy(*mf,vars_new[level][Vars::cons],Rho_comp,0,1,0);
            tag_found = true;

        // This allows dynamic refinement based on the value of the scalar/pressure/theta
        } else if ( (ref_tags[j].Field() == "scalar"  ) ||
                    (ref_tags[j].Field() == "pressure") ||
                    (ref_tags[j].Field() == "theta"   ) )
        {
            mf = std::make_unique<MultiFab>(grids[level], dmap[level], 1, 0);
            for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                auto& dfab = (*mf)[mfi];
                auto& sfab = vars_new[level][Vars::cons][mfi];
                if (ref_tags[j].Field() == "scalar") {
                    derived::erf_derscalar(bx, dfab, 0, 1, sfab, Geom(level), time, nullptr, level);
                } else if (ref_tags[j].Field() == "theta") {
                    derived::erf_dertheta(bx, dfab, 0, 1, sfab, Geom(level), time, nullptr, level);
                }
            } // mfi
            tag_found = true;
#ifdef ERF_USE_PARTICLES
        } else {
            // This allows dynamic refinement based on the number of particles per cell
            const auto& particles_namelist( particleData.getNames() );
            for (ParticlesNamesVector::size_type i = 0; i < particles_namelist.size(); i++) {
                std::string tmp_string(particles_namelist[i]+"_count");
                if (ref_tags[j].Field() == tmp_string) {
                    MultiFab temp_dat(grids[level], dmap[level], 1, 0); temp_dat.setVal(0);
                    particleData[particles_namelist[i]]->Increment(temp_dat, level);
                    mf = std::make_unique<MultiFab>(grids[level], dmap[level], 1, 0);
                    MultiFab::Copy(*mf, temp_dat, 0, 0, 1, 0);
                    tag_found = true;
                }
            }
#endif
        }
        if (tag_found) {
            ref_tags[j](tags,mf.get(),clearval,tagval,time,level,geom[level]);
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
            if (ppr.countval("in_box_lo")) {
                std::vector<Real> box_lo(3), box_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box <= max_level)
                {
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
                         (jlo%ref_ratio[lev_for_box-1][1] != 0) || ((jhi+1)%ref_ratio[lev_for_box-1][1] != 0) )
                         amrex::Error("Fine box is not legit with this ref_ratio");
                    boxes_at_level[lev_for_box].push_back(bx);
                    Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev
                if (init_type == "real" || init_type == "metgrid") {
                    if (num_boxes_at_level[lev_for_box] != num_files_at_level[lev_for_box]) {
                        amrex::Error("Number of boxes doesnt match number of input files");

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

