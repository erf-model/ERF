#include <ERF.H>

using namespace amrex;

//
// Tag cells for refinement -- this overrides the pure virtual function in AmrCore
//
void
ERF::ErrorEst (int level, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;
    for (int j=0; j < ref_tags.size(); ++j)
    {
        std::unique_ptr<MultiFab> mf;

        // This will work for static refinement
        ref_tags[j](tags,mf.get(),clearval,tagval,time,level,geom[level]);
  }
}

void
ERF::refinement_criteria_setup()
{
    if (max_level > 0)
    {
        std::string amr_prefix = "amr";
        ParmParse ppamr(amr_prefix);
        Vector<std::string> refinement_indicators;
        ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));

        for (int i=0; i<refinement_indicators.size(); ++i)
        {
            std::string ref_prefix = amr_prefix + "." + refinement_indicators[i];

            ParmParse ppr(ref_prefix);
            RealBox realbox;
            int lev_for_box;
            if (ppr.countval("in_box_lo")) {
                std::vector<Real> box_lo(3), box_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box <= max_level)
                {
                    ppr.getarr("in_box_lo",box_lo,0,2);
                    ppr.getarr("in_box_hi",box_hi,0,2);
                    box_lo[2] = geom[0].ProbLo(2);
                    box_hi[2] = geom[0].ProbHi(2);
                    realbox = RealBox(&(box_lo[0]),&(box_hi[0]));

                    amrex::Print() << "Reading " << realbox << " at level lev " << lev_for_box << std::endl;
                    num_boxes_at_level[lev_for_box] += 1;

                    auto dx = geom[lev_for_box].CellSize();
                    int ilo = box_lo[0]/dx[0]; int jlo = box_lo[1]/dx[1]; int klo = box_lo[2]/dx[2];
                    int ihi = box_hi[0]/dx[0]; int jhi = box_hi[1]/dx[1]; int khi = box_hi[2]/dx[2];
                    Box bx(IntVect(ilo,jlo,klo),IntVect(ihi,jhi,khi));
                    boxes_at_level[lev_for_box].push_back(bx);
                    amrex::Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev
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

