#include <algorithm>

#include <ERF.H>
#include <ERF_Derive.H>
#include <ERF_TerrainMetrics.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF
Box read_subdomain_from_wrfinput (int lev, const std::string& fname, int& ratio);
double read_start_time_from_wrfinput (int lev, const std::string& fname);
Box read_subdomain_from_metgrid (int lev, const std::string& fname, int& ratio, int& klo, int& khi);
#endif

/**
 * @brief Tag cells based on 2D distance from the storm eye.
 *
 * @param[in] cgeom Geometry defining the domain.
 * @param[out] tags Array of tagged cells for refinement.
 * @param[in] eye_x X-coordinate of the eye center.
 * @param[in] eye_y Y-coordinate of the eye center.
 * @param[in] rad_tag Radius within which cells are tagged.
 */
void
tag_on_distance_from_eye (const Geometry& cgeom, TagBoxArray* tags,
                          const Real eye_x, const Real eye_y, const Real rad_tag);

namespace {

/**
 * Whether a field name, already lowercased, is one of the built-in refinement fields
 */
bool
is_builtin_tag_field (const std::string& f)
{
    return f == "density"   || f == "theta"    || f == "scalar"           ||
           f == "qv"        || f == "qc"       || f == "qi"               ||
           f == "qr"        || f == "qs"       || f == "qg"               || f == "qt" ||
           f == "vorticity" || f == "helicity" || f == "max_reflectivity" ||
           f == "velmag"    || f == "u"        || f == "v"                || f == "w"  ||
           f == "omega"     || f == "terrain_blanking";
}

/**
 * Whether a field name is one of the moisture mixing ratios qv, qc, qi, qr, qs, qg
 */
bool
is_moist_species_field (const std::string& f)
{
    return f == "qv" || f == "qc" || f == "qi" ||
           f == "qr" || f == "qs" || f == "qg";
}

/**
 * RhoQ component of a moisture mixing ratio field, or -1 if the moisture model does not carry it
 */
int
moist_comp_for_field (const MoistureComponentIndices& mi, const std::string& f)
{
    if      (f == "qv") { return mi.qv; }
    else if (f == "qc") { return mi.qc; }
    else if (f == "qi") { return mi.qi; }
    else if (f == "qr") { return mi.qr; }
    else if (f == "qs") { return mi.qs; }
    else if (f == "qg") { return mi.qg; }
    return -1;
}

} // namespace

/**
 * Function to tag cells for refinement -- this overrides the pure virtual function in AmrCore
 *
 * @param[in ] levc level of refinement at which we tag cells (0 is coarsest level)
 * @param[out] tags array of tagged cells
 * @param[in ] time current time
 * @param[in ] ngrow number of ghost cells (not used here)
*/
void
ERF::ErrorEst (int levc, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const char tagval = TagBox::SET;

#ifdef ERF_USE_NETCDF
    if ((solverChoice.init_type == InitType::WRFInput) || (solverChoice.init_type == InitType::Metgrid)) {
        int ratio;
        Box subdomain;

        // Box isub at level levc+1 belongs to file isub. Boxes may already have been defined in the
        // refinement_criteria_setup routine; a box that exists is used as long as it is contained
        // inside the subdomain its file gives, and a file without a box gets the file's subdomain.
        //
        // The met_em files at a level are successive times of a single nest, so only the first
        // one defines a box there.
        const int num_subdomains = (solverChoice.init_type == InitType::Metgrid) ?
                                   amrex::min(1, static_cast<int>(nc_init_file[levc+1].size())) :
                                   static_cast<int>(nc_init_file[levc+1].size());

        if (!nc_init_file[levc+1].empty())
        {
            double levc_start_time = read_start_time_from_wrfinput(levc  , nc_init_file[levc  ][0]);
            if (solverChoice.init_type == InitType::WRFInput) {
                amrex::Print() << " WRFInput       time at level " << levc << " is " << levc_start_time << std::endl;
            } else if (solverChoice.init_type == InitType::Metgrid) {
                amrex::Print() << " met_em         time at level " << levc << " is " << levc_start_time << std::endl;
            }

            for (int isub = 0; isub < num_subdomains; isub++) {
                if (!have_read_nc_init_file[levc+1][isub])
                {
                    double levf_start_time = read_start_time_from_wrfinput(levc+1, nc_init_file[levc+1][isub]);
                    if (solverChoice.init_type == InitType::WRFInput) {
                        amrex::Print() << " WRFInput start_time at level " << levc+1 << " is " << levf_start_time << std::endl;
                    } else if (solverChoice.init_type == InitType::Metgrid) {
                        amrex::Print() << " met_em   start time at level " << levc+1 << " is " << levf_start_time << std::endl;
                    }

                    // We assume there is only one subdomain at levc; otherwise we don't know
                    //     which one is the parent of the fine region we are trying to create
                    AMREX_ALWAYS_ASSERT(subdomains[levc].size() == 1);

                    if ( levf_start_time <= (levc_start_time + t_new[levc]) ) {
                        if (solverChoice.init_type == InitType::WRFInput) {
                            amrex::Print() << " WRFInput file to read: " << nc_init_file[levc+1][isub] << std::endl;
                            // WRF starts are relative to the full parent file.  Account
                            // for every ancestor before locating this nest in ERF.
                            subdomain = wrfinput_file_box(levc+1, isub, &ratio);
                            amrex::Print() << " WRFInput subdomain " << isub << " at level " << levc+1 << " is " << subdomain << std::endl;
                        } else if (solverChoice.init_type == InitType::Metgrid) {
                            amrex::Print() << "met_em file to read: " << nc_init_file[levc+1][0] << std::endl;
                            const Box& domain = geom[levc].Domain();
                            int klo = domain.smallEnd(2);
                            int khi = domain.bigEnd(2);
                            subdomain = read_subdomain_from_metgrid(levc, nc_init_file[levc+1][0], ratio, klo, khi);
                            amrex::Print() << " met_em subdomain at level " << levc+1 << " is " << subdomain << std::endl;
                        }

                        if ( (ratio != ref_ratio[levc][0]) || (ratio != ref_ratio[levc][1]) ) {
                            amrex::Print() << "File " << nc_init_file[levc+1][isub] << " has refinement ratio = " << ratio << std::endl;
                            amrex::Print() << "The inputs file has refinement ratio = " << ref_ratio[levc] << std::endl;
                            amrex::Abort("These must be the same in the horizontal -- please edit your inputs file and try again.");
                        }

                        if (solverChoice.init_type == InitType::WRFInput) {
                            //
                            // A WRF nest is refined in the horizontal only: wrfinput_d02 carries
                            // exactly the eta levels of wrfinput_d01.  The z extent just read is
                            // therefore already in the PARENT's vertical index space, so only the
                            // horizontal directions are coarsened here; any vertical refinement is
                            // applied when the box is refined below.  Clip that extent to the
                            // parent's domain in case the file carries more layers than this run
                            // was given cells in z.
                            //
                            subdomain.setSmall(2, amrex::max(subdomain.smallEnd(2),
                                                             geom[levc].Domain().smallEnd(2)));
                            subdomain.setBig  (2, amrex::min(subdomain.bigEnd(2),
                                                             geom[levc].Domain().bigEnd(2)));
                            subdomain.coarsen(IntVect(ref_ratio[levc][0],ref_ratio[levc][1],1));
                        } else {
                            subdomain.coarsen(ref_ratio[levc]);
                        }

                        // Metgrid boxes are relative to the parent ERF region.  WRF
                        // boxes have already been placed in global coordinates above.
                        if (solverChoice.init_type == InitType::Metgrid) {
                            Box coarser_level(subdomains[levc][0].minimalBox());
                            subdomain.shift(coarser_level.smallEnd());
                        }

                        if (verbose > 0) {
                            amrex::Print() << " Crse version of subdomain available for tagging is" << subdomain << std::endl;
                        }

                        Box new_fine(subdomain);
                        if (solverChoice.init_type == InitType::WRFInput) {
                            // The horizontal ratio comes from the file (and was just checked
                            // against the inputs file); the vertical one is whatever the inputs
                            // file asked for, since the file itself has nothing to say about it.
                            new_fine.refine(IntVect(ratio,ratio,ref_ratio[levc][2]));
                        } else if (solverChoice.init_type == InitType::Metgrid) {
                            new_fine.refine(ref_ratio[levc]);
                        }
                        const int num_boxes = static_cast<int>(boxes_at_level[levc+1].size());
                        if (isub < num_boxes) {
                            if (!new_fine.contains(boxes_at_level[levc+1][isub])) {
                                amrex::Print() << "\n";
                                amrex::Print() << "Box available in wrfinputs file             " << new_fine << std::endl;
                                amrex::Print() << "Box requested for refinement in inputs file " << boxes_at_level[levc+1][isub] << std::endl;
                                amrex::Abort("Specified boxes must be contained within boxes specified in wrfinput at this level");
                            }
                        } else if (isub == num_boxes) {
                            num_boxes_at_level[levc+1] += 1;
                            boxes_at_level[levc+1].push_back(new_fine);
                        } else {
                            // An earlier file has no box yet, so this file's box cannot be stored at index isub
                            amrex::Print() << "File " << nc_init_file[levc+1][isub] << " starts before "
                                           << nc_init_file[levc+1][num_boxes] << std::endl;
                            amrex::Abort("The nc_init_file files at a level must be listed in order of start time");
                        }

                        Box coarsened_bx(boxes_at_level[levc+1][isub]); coarsened_bx.coarsen(ref_ratio[levc]);

                        // A nest above level 1 must fit inside its parent's ERF patch, which may
                        //    cover only part of the parent's file, with the proper-nesting buffer
                        //    (except at the domain boundary).  Otherwise tagging would clip it
                        //    quietly while the initial data were still read for the whole box.
                        if (solverChoice.init_type == InitType::WRFInput && levc > 0) {
                            const Box needed = amrex::grow(coarsened_bx, nProper()) & geom[levc].Domain();
                            if (!grids[levc].contains(needed)) {
                                amrex::Print() << "\n";
                                amrex::Print() << "Box to refine at level " << levc+1 << "             "
                                               << boxes_at_level[levc+1][isub] << std::endl;
                                amrex::Print() << "The same at level " << levc << ", with a buffer of "
                                               << nProper() << " cells   " << needed << std::endl;
                                amrex::Print() << "Grids of level " << levc << " cover            "
                                               << grids[levc].minimalBox() << std::endl;
                                amrex::Abort("The level-" + std::to_string(levc+1) + " region of " +
                                             nc_init_file[levc+1][isub] + " does not fit inside the level-" +
                                             std::to_string(levc) + " grids with the proper-nesting buffer;"
                                             " enlarge the level-" + std::to_string(levc) + " region, or crop"
                                             " level " + std::to_string(levc+1) + " with a refinement indicator");
                            }
                        }

                        for (MFIter mfi(tags); mfi.isValid(); ++mfi)
                        {
                            auto tag_arr = tags.array(mfi);  // Get device-accessible array

                            Box bx = mfi.validbox() & coarsened_bx;

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
                           false, true);
    }

    const auto dx      = geom[levc].CellSizeArray();
    const auto prob_lo = geom[levc].ProbLoArray();

    for (auto& crit : ref_criteria)
    {
        // Move the box if the inputs file asks for it
        if (crit.real_box.ok()) {
            std::string ref_prefix = pp_prefix + "." + crit.name;
            if (crit.box_z != RefinementCriterion::BoxZ::Grid) {
                Real offset_x, offset_y;
                refinement_box_offset(ref_prefix, time, offset_x, offset_y);
                crit.real_box = crit.base_box;
                crit.real_box.setLo(0, crit.base_box.lo(0) + offset_x);
                crit.real_box.setHi(0, crit.base_box.hi(0) + offset_x);
                crit.real_box.setLo(1, crit.base_box.lo(1) + offset_y);
                crit.real_box.setHi(1, crit.base_box.hi(1) + offset_y);
            } else {
                int lev_for_box;
                update_box_for_refinement(ref_prefix, lev_for_box, crit.real_box, time);
            }
        }

        if (!crit.is_active(levc, time)) { continue; }

        // Cells that satisfy the conditions -- all of them, or any one, as the criterion asks;
        //    with no conditions every cell does, and only the box matters
        iMultiFab mask(grids[levc], dmap[levc], 1, 0);
        mask.setVal((crit.require_all || crit.conditions.empty()) ? 1 : 0);

        for (const auto& cond : crit.conditions)
        {
            //
            // This mf must have ghost cells because we may take differences between adjacent values
            //
            MultiFab mf(grids[levc], dmap[levc], 1, 1);
            mf.setVal(0.0);

            compute_tag_field(levc, cond.field, time, mf, tags);

            const auto test        = cond.test;
            const Real threshold   = cond.threshold(levc);
            const bool require_all = crit.require_all;

            auto const& mask_ma = mask.arrays();
            auto const& dat_ma  = mf.const_arrays();
            ParallelFor(mask, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
            {
                auto const& dat = dat_ma[bi];
                bool pass;
                if (test == TagCondition::Test::Greater) {
                    pass = (dat(i,j,k) >= threshold);
                } else if (test == TagCondition::Test::Less) {
                    pass = (dat(i,j,k) <= threshold);
                } else {
                    // Largest difference from any of the six face neighbors
                    const Real d = amrex::max(std::abs(dat(i+1,j,k) - dat(i,j,k)), std::abs(dat(i,j,k) - dat(i-1,j,k)),
                                              std::abs(dat(i,j+1,k) - dat(i,j,k)), std::abs(dat(i,j,k) - dat(i,j-1,k)),
                                              std::abs(dat(i,j,k+1) - dat(i,j,k)), std::abs(dat(i,j,k) - dat(i,j,k-1)));
                    pass = (d >= threshold);
                }
                int& m = mask_ma[bi](i,j,k);
                m = require_all ? (m && pass) : (m || pass);
            });
        }

        // Tag the cells that pass and whose centers lie inside the box. In z the center is at
        //    its height above ground or its physical height when the box restricts conditions,
        //    and otherwise (a static nest, possibly given by indices) at its height on the
        //    undeformed grid
        const bool has_box   = crit.real_box.ok();
        const bool z_grid    = (crit.box_z == RefinementCriterion::BoxZ::Grid);
        const bool z_aground = (crit.box_z == RefinementCriterion::BoxZ::AboveGround);
        const Real xlo = crit.real_box.lo(0), xhi = crit.real_box.hi(0);
        const Real ylo = crit.real_box.lo(1), yhi = crit.real_box.hi(1);
        const Real zlo = crit.real_box.lo(2), zhi = crit.real_box.hi(2);

        auto const& tag_ma  = tags.arrays();
        auto const& mask_ma = mask.const_arrays();
        auto const& zcc_ma  = z_phys_cc[levc]->const_arrays();
        auto const& zsfc_ma = z_sfc_nd[levc]->const_arrays();
        ParallelFor(mask, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
        {
            if (!mask_ma[bi](i,j,k)) { return; }
            if (has_box) {
                const Real x = prob_lo[0] + (Real(i) + myhalf) * dx[0];
                const Real y = prob_lo[1] + (Real(j) + myhalf) * dx[1];
                const Real z = z_grid    ? prob_lo[2] + (Real(k) + myhalf) * dx[2]
                             : z_aground ? zcc_ma[bi](i,j,k) - Compute_Z_AtWFace(i,j,0,zsfc_ma[bi])
                                         : zcc_ma[bi](i,j,k);
                if (x < xlo || x > xhi || y < ylo || y > yhi || z < zlo || z > zhi) { return; }
            }
            tag_ma[bi](i,j,k) = tagval;
        });
        Gpu::streamSynchronize();
    } // loop over criteria

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

            double ref_start_time = -1.0;
            ppr.queryAdd("start_time",ref_start_time);

            if (time >= ref_start_time) {

                Real max_radius = -one;
                ppr.get("max_radius", max_radius);

                // Create the volume-weighted sum of (rho qv) in each column
                MultiFab mf_qv_int(ba2d[levc], dmap[levc], 1, 0); mf_qv_int.setVal(0.);

                // Define the 2D MultiFab holding the column-integrated (rho qv)
                volWgtColumnSum(levc, S_new, RhoQ1_comp, mf_qv_int, *detJ_cc[levc]);

                // Find the max value in the domain
                IntVect eye = mf_qv_int.maxIndex(0);

                Real eye_x = prob_lo[0] + (eye[0] + myhalf) * dx[0];
                Real eye_y = prob_lo[1] + (eye[1] + myhalf) * dx[1];

                tag_on_distance_from_eye(geom[levc], &tags, eye_x, eye_y, max_radius);
            }
        }
    }
}

/**
 * Fill mf with the refinement field named field at level lev
 *
 * @param[in ] levc  level at which we tag cells
 * @param[in ] field name of the field, as given by erf.<indicator>.field_name
 * @param[in ] time  current time
 * @param[out] mf    one-component MultiFab with one ghost cell, zeroed by the caller
 * @param[out] tags  tags set directly by the hurricane tracker, which "velmag" may invoke
 */
void
ERF::compute_tag_field (int levc, const std::string& field, Real time, MultiFab& mf, TagBoxArray& tags)
{
    MultiFab& U_new = vars_new[levc][Vars::xvel];
    MultiFab& V_new = vars_new[levc][Vars::yvel];
    MultiFab& W_new = vars_new[levc][Vars::zvel];

    // Built-in field names are matched case-insensitively; the particle
    // fields below keep the name as given since species names may carry case
    const std::string field_lc = amrex::toLower(field);

    // This allows dynamic refinement based on the value of the density
    if (field_lc == "density")
    {
        MultiFab::Copy(mf,vars_new[levc][Vars::cons],Rho_comp,0,1,1);

    // Refinement based on a moisture mixing ratio (qv, qc, qi, qr, qs, qg).
    // The map from name to RhoQ component depends on the moisture model and
    // is held in solverChoice.moisture_indices.
    } else if (is_moist_species_field(field_lc))
    {
        const int qcomp = moist_comp_for_field(solverChoice.moisture_indices, field_lc);
        AMREX_ALWAYS_ASSERT(qcomp >= 0);
        MultiFab::Copy(  mf, vars_new[levc][Vars::cons], qcomp,    0, 1, 1);
        MultiFab::Divide(mf, vars_new[levc][Vars::cons], Rho_comp, 0, 1, 1);

    // qt = total condensed water mixing ratio (qc + qi + qr + qs + qg).
    // Excludes qv by convention here.
    } else if (field_lc == "qt") {
        const auto& mi = solverChoice.moisture_indices;
        const int idx_qc = mi.qc, idx_qi = mi.qi, idx_qr = mi.qr,
                  idx_qs = mi.qs, idx_qg = mi.qg;
        AMREX_ALWAYS_ASSERT(idx_qc >= 0 || idx_qi >= 0 || idx_qr >= 0 ||
                            idx_qs >= 0 || idx_qg >= 0);
        mf.setVal(0.0);
        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox();
            auto qt_arr   = mf.array(mfi);
            auto cons_arr = vars_new[levc][Vars::cons].const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real rho_inv = Real(1.0) / cons_arr(i,j,k,Rho_comp);
                Real q = Real(0.0);
                if (idx_qc >= 0) { q += cons_arr(i,j,k,idx_qc) * rho_inv; }
                if (idx_qi >= 0) { q += cons_arr(i,j,k,idx_qi) * rho_inv; }
                if (idx_qr >= 0) { q += cons_arr(i,j,k,idx_qr) * rho_inv; }
                if (idx_qs >= 0) { q += cons_arr(i,j,k,idx_qs) * rho_inv; }
                if (idx_qg >= 0) { q += cons_arr(i,j,k,idx_qg) * rho_inv; }
                qt_arr(i,j,k) = q;
            });
        }

    // This allows dynamic refinement based on the value of the z-component of vorticity
    } else if (field_lc == "vorticity" ) {
        MultiFab mf_cc_vel(grids[levc], dmap[levc], AMREX_SPACEDIM, IntVect(1,1,1));
        average_face_to_cellcenter(mf_cc_vel,0,Array<const MultiFab*,3>{&U_new, &V_new, &W_new}, 1);

        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto& dfab = mf[mfi];
            auto& sfab = mf_cc_vel[mfi];
            auto& zfab = (*z_phys_cc[levc])[mfi];
            derived::erf_dervortz(bx, dfab, 0, 1, sfab, zfab, Geom(levc), time, nullptr, levc);
        }

    // This allows dynamic refinement based on the value of the scalar/theta
    } else if ( (field_lc == "scalar"  ) ||
                (field_lc == "theta"   ) )
    {
        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto& dfab = mf[mfi];
            auto& sfab = vars_new[levc][Vars::cons][mfi];
            auto& zfab = (*z_phys_cc[levc])[mfi];
            if (field_lc == "scalar") {
                derived::erf_derscalar(bx, dfab, 0, 1, sfab, zfab, Geom(levc), time, nullptr, levc);
            } else if (field_lc == "theta") {
                derived::erf_dertheta(bx, dfab, 0, 1, sfab, zfab, Geom(levc), time, nullptr, levc);
            }
        } // mfi
    // This allows dynamic refinement based on the value of updraft helicity
    } else if (field_lc == "helicity")
    {
        MultiFab mf_cc_vel(grids[levc], dmap[levc], AMREX_SPACEDIM, IntVect(1,1,1));
        average_face_to_cellcenter(mf_cc_vel,0,Array<const MultiFab*,3>{&U_new, &V_new, &W_new}, 1);

        // NOTE: helicity is a vertical integral so we must not tile in z
        for (MFIter mfi(mf, TileNoZ()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto& dfab = mf[mfi];
            auto& sfab = mf_cc_vel[mfi];
            auto& zfab = (*z_phys_cc[levc])[mfi];

            derived::erf_derhelicity(bx, dfab, 0, 1, sfab, zfab, Geom(levc), time, nullptr, levc);
        }
    } else if (field_lc == "max_reflectivity")
    {
        // NOTE: this tests the data layout the reflectivity kernel assumes rather
        //       than enumerating moisture models, so it picks up any scheme that
        //       carries rain/snow/graupel in RhoQ4/RhoQ5/RhoQ6
        if (solverChoice.moisture_indices.has_reflectivity_species()) {
            // NOTE: this takes a max over the column so we must not tile in z
            for (MFIter mfi(mf, TileNoZ()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                auto& dfab = mf[mfi];
                auto& sfab = vars_new[levc][Vars::cons][mfi];
                auto& zfab = (*z_phys_cc[levc])[mfi];

                derived::erf_dermaxreflectivity(bx, dfab, 0, 1, sfab, zfab, Geom(levc), time, nullptr, levc);
            }
        } else {
            Abort("Max reflectivity requires a microphysics model carrying rain, snow and graupel.");
        }
    // This allows dynamic refinement based on the terrain blanking
    } else if ( (SolverChoice::terrain_type == TerrainType::ImmersedForcing) &&
                (field_lc == "terrain_blanking") )
    {
        MultiFab::Copy(mf,*terrain_blanking[levc],0,0,1,1);
    }
    else if (field_lc == "velmag")
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
            for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                auto& dfab = mf[mfi];
                auto& sfab = mf_cc_vel[0][mfi];
                auto& zfab = (*z_phys_cc[levc])[mfi];
                derived::erf_dermagvel(bx, dfab, 0, 1, sfab, zfab, Geom(levc), time, nullptr, levc);
            }
        }
    }
    // This allows dynamic refinement based on a cell-centered velocity component;
    //    omega is the contravariant vertical velocity, which differs from w only with terrain
    else if (field_lc == "u" ||
             field_lc == "v" ||
             field_lc == "w" ||
             field_lc == "omega")
    {
        // Omega equals w unless the coordinate surfaces follow the terrain
        const bool fitted_mesh = (SolverChoice::mesh_type != MeshType::ConstantDz) &&
                                 (SolverChoice::terrain_type == TerrainType::StaticFittedMesh ||
                                  SolverChoice::terrain_type == TerrainType::MovingFittedMesh);

        if (field_lc != "omega" || !fitted_mesh) {
            // Average only the component we need to the cell centers, ghost cells included
            const int dir = (field_lc == "u") ? 0 : (field_lc == "v") ? 1 : 2;
            const MultiFab& vel = (dir == 0) ? U_new : (dir == 1) ? V_new : W_new;
            const int ioff = (dir == 0), joff = (dir == 1), koff = (dir == 2);
            for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(1);
                const Array4<      Real>& d_arr = mf.array(mfi);
                const Array4<const Real>& v_arr = vel.const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    d_arr(i,j,k) = myhalf * (v_arr(i,j,k) + v_arr(i+ioff,j+joff,k+koff));
                });
            }
        } else {
            // Omega on the z-faces of the box grown by one cell. As in the dycore, it is zero
            //    on (and below) the bottom face and equal to w on (and above) the top face,
            //    which validate_flat_top makes flat.
            // NOTE: with MovingFittedMesh the mesh velocity z_t is not subtracted, so this is
            //       the velocity normal to the coordinate surfaces as if they were at rest
            const int klo = geom[levc].Domain().smallEnd(2);
            const int khi = geom[levc].Domain().bigEnd(2);
            const auto dxInv = geom[levc].InvCellSizeArray();

            MultiFab omega(W_new.boxArray(), dmap[levc], 1, IntVect(1,1,1));
            for (MFIter mfi(omega, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(1);
                const Array4<      Real>& om_arr = omega.array(mfi);
                const Array4<const Real>& u_arr  = U_new.const_array(mfi);
                const Array4<const Real>& v_arr  = V_new.const_array(mfi);
                const Array4<const Real>& w_arr  = W_new.const_array(mfi);
                const Array4<const Real>& z_nd   = z_phys_nd[levc]->const_array(mfi);
                const Array4<const Real>& mf_u   = mapfac[levc][MapFacType::u_x]->const_array(mfi);
                const Array4<const Real>& mf_v   = mapfac[levc][MapFacType::v_y]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    om_arr(i,j,k) = (k <= klo  ) ? zero :
                                    (k >= khi+1) ? w_arr(i,j,k) :
                        OmegaFromW(i,j,k,w_arr(i,j,k),u_arr,v_arr,mf_u,mf_v,z_nd,dxInv);
                });
            }

            for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(1);
                const Array4<      Real>& d_arr  = mf.array(mfi);
                const Array4<const Real>& om_arr = omega.const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    d_arr(i,j,k) = myhalf * (om_arr(i,j,k) + om_arr(i,j,k+1));
                });
            }
        }

#ifdef ERF_USE_PARTICLES
    } else {
        //
        // Particle-derived refinement.  Two forms of `field_name` are
        // supported:
        //   <species>_count          : particle count per cell
        //   <species>_<mesh_var>     : Eulerian mesh variable provided
        //                              by the species' `meshPlotVarNames()`
        //                              (e.g. `super_droplets_moisture_mass_density`,
        //                              `tracer_particles_mass_density`).
        //
        // In both cases the field is deposited at every level in
        // [levc, finest_level] and averaged down level-by-level so that
        // a signal present only on a finer level still triggers
        // refinement at the coarser level being tagged.  Without this,
        // particles localised at level 1 would not register at level 0
        // and the fine grid would disappear at the next regrid.
        //
        const auto& particles_namelist( particleData.getNames() );
        mf.setVal(0.0);
        for (ParticlesNamesVector::size_type i = 0; i < particles_namelist.size(); i++)
        {
            auto* pc = particleData[particles_namelist[i]];
            const std::string& sp_name = particles_namelist[i];

            const std::string count_str = sp_name + "_count";
            const std::string prefix    = sp_name + "_";
            std::string mesh_var;
            if (field != count_str
                && field.size() > prefix.size()
                && field.compare(0, prefix.size(), prefix) == 0)
            {
                const std::string suffix = field.substr(prefix.size());
                for (const auto& v : pc->meshPlotVarNames()) {
                    if (v == suffix) { mesh_var = v; break; }
                }
            }
            if (field != count_str && mesh_var.empty()) { continue; }

            pc->resizeData();
            const int pc_nlevs = static_cast<int>(pc->GetParticles().size());

            // Deposit at each level into per-level MultiFabs.
            Vector<MultiFab> per_lev(finest_level+1);
            for (int lev = levc; lev <= finest_level; lev++) {
                per_lev[lev].define(grids[lev], dmap[lev], 1, 0);
                per_lev[lev].setVal(0);
                if (field == count_str) {
                    if (lev < pc_nlevs) {
                        pc->IncrementWithTotal(per_lev[lev], lev);
                    }
                } else {
                    pc->computeMeshVar(mesh_var, per_lev[lev],
                                       *z_phys_nd[lev], lev);
                }
            }

            // Average down level-by-level from finest to levc.  This
            // avoids multi-level coarsening (e.g. L2->L0 with ratio
            // (4,1,4)) which can fail when fine-level boxes are not
            // aligned to the composite refinement ratio.
            for (int lev = finest_level; lev > levc; lev--) {
                MultiFab temp_crse(grids[lev-1], dmap[lev-1], 1, 0);
                temp_crse.setVal(0);
                average_down(per_lev[lev], temp_crse,
                             0, 1, ref_ratio[lev-1]);
                MultiFab::Add(per_lev[lev-1], temp_crse, 0, 0, 1, 0);
            }

            MultiFab::Copy(mf, per_lev[levc], 0, 0, 1, 0);
        }
#endif
    }
}

/**
 * Function to define the refinement criteria based on user input
 *
 * An indicator erf.<name> tests either one condition, given by field_name and one of
 * value_greater, value_less or adjacent_difference_greater, or several, listed by name in
 * erf.<name>.conditions and each given by those same keys under erf.<name>.<condition>.
 * Several conditions must all hold (erf.<name>.combine = and, the default) or any one of
 * them (combine = or). An in_box_lo/in_box_hi box then restricts the tags to that region,
 * with z compared against the height of the cell centers above ground or, with
 * erf.<name>.box_z_coordinate = physical, against z_phys_cc.
*/

void
ERF::refinement_criteria_setup ()
{
    if (max_level > 0)
    {
        ParmParse pp(pp_prefix);
        Vector<std::string> refinement_indicators;
        pp.queryarr("refinement_indicators",refinement_indicators,0,pp.countval("refinement_indicators"));

        // With hurricane_tracker listed, a velmag field tags through HurricaneTracker directly
        //    rather than through the conditions, combine mode and box below
        const bool has_hurricane_tracker =
            std::find(refinement_indicators.begin(), refinement_indicators.end(), "hurricane_tracker")
            != refinement_indicators.end();

        // Check each field_name here so that an unknown name, or a field the configuration
        //    cannot provide, aborts at setup rather than tagging against a field that stays
        //    zero. Built-in names, listed in lowercase in is_builtin_tag_field, are matched
        //    case-insensitively; particle fields are matched as given.
        auto check_field = [&] (const std::string& field, const std::string& indicator)
        {
            const std::string field_lc = amrex::toLower(field);
            if (is_builtin_tag_field(field_lc)) {
                const auto& mi = solverChoice.moisture_indices;
                bool available = true;
                if (is_moist_species_field(field_lc)) {
                    available = (moist_comp_for_field(mi, field_lc) >= 0);
                } else if (field_lc == "qt") {
                    available = (mi.qc >= 0 || mi.qi >= 0 || mi.qr >= 0 ||
                                 mi.qs >= 0 || mi.qg >= 0);
                }
                if (!available) {
                    amrex::Abort("Refinement field_name '" + field +
                                 "' is not available for the configured moisture model");
                }
                if (field_lc == "terrain_blanking" &&
                    solverChoice.terrain_type != TerrainType::ImmersedForcing) {
                    amrex::Abort("Refinement field_name '" + field +
                                 "' requires erf.terrain_type = ImmersedForcing");
                }
            } else {
                bool is_particle_field = false;
#ifdef ERF_USE_PARTICLES
                // Species containers are not built yet, so only the species name is checked
                //    here; the mesh variable is matched in compute_tag_field
                auto matches_species = [&field] (const std::string& sp_name) {
                    const std::string prefix = sp_name + "_";
                    return field.size() > prefix.size() &&
                           field.compare(0, prefix.size(), prefix) == 0;
                };
                for (const auto& sp_name : particleData.getNames()) {
                    if (matches_species(sp_name)) { is_particle_field = true; }
                }
                for (const auto& sp_name : particleData.getNamesUnalloc()) {
                    if (matches_species(sp_name)) { is_particle_field = true; }
                }
#endif
                if (!is_particle_field) {
                    amrex::Abort("Unrecognized refinement field_name '" + field +
                                 "' for indicator " + indicator);
                }
            }
        };

        // Read the condition under prefix, if it gives a test; with more than one test the
        //    first in the order below is used
        auto read_condition = [&] (const std::string& prefix, const std::string& indicator,
                                   TagCondition& cond) -> bool
        {
            ParmParse ppc(prefix);
            const std::pair<std::string, TagCondition::Test> tests[] = {
                {"value_greater",               TagCondition::Test::Greater},
                {"value_less",                  TagCondition::Test::Less},
                {"adjacent_difference_greater", TagCondition::Test::AdjacentDifferenceGreater}};

            bool found = false;
            for (const auto& [key, test] : tests) {
                const int num_val = ppc.countval(key.c_str());
                if (num_val > 0) {
                    cond.test = test;
                    ppc.getarr(key.c_str(), cond.values, 0, num_val);
                    found = true;
                    break;
                }
            }

            if (ppc.countval("field_name") > 0) {
                ppc.get("field_name", cond.field);
                check_field(cond.field, indicator);
            } else if (found) {
                amrex::Abort("Refinement condition " + prefix + " gives a test but no field_name");
            }
            return found;
        };

        for (int i=0; i<refinement_indicators.size(); ++i)
        {
            std::string ref_prefix = pp_prefix + "." + refinement_indicators[i];

            ParmParse ppr(ref_prefix);

            RefinementCriterion crit;
            crit.name = refinement_indicators[i];

            if (ppr.countval("start_time") > 0) {
                double ref_min_time; ppr.get("start_time",ref_min_time);
                crit.min_time = static_cast<Real>(ref_min_time);
            }

            if (ppr.countval("end_time") > 0) {
                double ref_max_time; ppr.get("end_time",ref_max_time);
                crit.max_time = static_cast<Real>(ref_max_time);
            }

            if (ppr.countval("max_level") > 0) {
                ppr.get("max_level",crit.max_level);
            }

            if (ppr.countval("combine") > 0) {
                std::string combine; ppr.get("combine", combine);
                combine = amrex::toLower(combine);
                if (combine != "and" && combine != "or") {
                    amrex::Abort("erf." + crit.name + ".combine must be 'and' or 'or'");
                }
                crit.require_all = (combine == "and");
            }

            if (ppr.countval("conditions") > 0) {
                Vector<std::string> condition_names;
                ppr.getarr("conditions", condition_names, 0, ppr.countval("conditions"));
                for (const auto& cname : condition_names) {
                    TagCondition cond;
                    if (!read_condition(ref_prefix + "." + cname, crit.name, cond)) {
                        amrex::Abort("Refinement condition " + ref_prefix + "." + cname +
                                     " must give field_name and one of value_greater, value_less"
                                     " or adjacent_difference_greater");
                    }
                    crit.conditions.push_back(cond);
                }
                TagCondition unused;
                if (ppr.countval("field_name") > 0 || read_condition(ref_prefix, crit.name, unused)) {
                    amrex::Abort("Refinement indicator " + crit.name + " gives conditions, so its"
                                 " tests must be given under erf." + crit.name + ".<condition>");
                }
            } else {
                TagCondition cond;
                if (read_condition(ref_prefix, crit.name, cond)) {
                    crit.conditions.push_back(cond);
                } else if (ppr.countval("field_name") > 0) {
                    amrex::Abort("Refinement indicator " + crit.name + " gives field_name but none of"
                                 " value_greater, value_less or adjacent_difference_greater");
                }
            }

            if (has_hurricane_tracker) {
                const bool velmag = std::any_of(crit.conditions.begin(), crit.conditions.end(),
                    [] (const TagCondition& c) { return amrex::toLower(c.field) == "velmag"; });
                if ( (velmag || crit.name == "hurricane_tracker") &&
                     (ppr.countval("conditions") > 0 || ppr.countval("in_box_lo") > 0) ) {
                    amrex::Abort("Refinement indicator " + crit.name + ": with hurricane_tracker listed,"
                                 " an indicator testing velmag (and hurricane_tracker itself) must give a"
                                 " single field_name and value_greater, with no conditions and no box");
                }
            }

            if (ppr.countval("box_z_coordinate") > 0 &&
                (crit.conditions.empty() || ppr.countval("in_box_lo") == 0)) {
                amrex::Abort("erf." + crit.name + ".box_z_coordinate applies only to an in_box_lo/in_box_hi"
                             " box on an indicator that tests a field");
            }

            if (!crit.conditions.empty())
            {
                // The box only restricts where the conditions may tag, so it is read here rather
                //    than registered as a static refinement box by read_box_for_refinement
                if (ppr.countval("in_box_lo_indices")      > 0 || ppr.countval("in_box_hi_indices")      > 0 ||
                    ppr.countval("in_box_lo_indices_crse") > 0 || ppr.countval("in_box_hi_indices_crse") > 0) {
                    amrex::Abort("Refinement indicator " + crit.name + " tests a field, so its region"
                                 " must be given by in_box_lo and in_box_hi");
                }
                const int num_lo = ppr.countval("in_box_lo");
                const int num_hi = ppr.countval("in_box_hi");
                if (num_lo != num_hi || (num_lo > 0 && num_lo < AMREX_SPACEDIM-1)) {
                    amrex::Abort("erf." + crit.name + ".in_box_lo and in_box_hi must each give"
                                 " 2 (x,y) or 3 (x,y,z) values");
                }
                if (num_lo > 0) {
                    Vector<Real> rbox_lo, rbox_hi;
                    ppr.getarr("in_box_lo", rbox_lo, 0, num_lo);
                    ppr.getarr("in_box_hi", rbox_hi, 0, num_hi);
                    if (num_lo < AMREX_SPACEDIM) {
                        rbox_lo.push_back(std::numeric_limits<Real>::lowest());
                        rbox_hi.push_back(std::numeric_limits<Real>::max());
                    }
                    crit.base_box = RealBox(rbox_lo.data(), rbox_hi.data());
                    if (!crit.base_box.ok()) {
                        amrex::Abort("erf." + crit.name + ".in_box_hi must exceed in_box_lo in every direction");
                    }
                    crit.real_box = crit.base_box;

                    crit.box_z = RefinementCriterion::BoxZ::AboveGround;
                    if (ppr.countval("box_z_coordinate") > 0) {
                        std::string zcoord; ppr.get("box_z_coordinate", zcoord);
                        zcoord = amrex::toLower(zcoord);
                        if (zcoord == "physical") {
                            crit.box_z = RefinementCriterion::BoxZ::Physical;
                        } else if (zcoord != "above_ground") {
                            amrex::Abort("erf." + crit.name + ".box_z_coordinate must be"
                                         " 'above_ground' or 'physical'");
                        }
                    }

                    // The ground is taken from the bottom of z_phys_nd, which is the flat bottom of the
                    //    domain rather than the terrain when the terrain is immersed
                    if (crit.box_z == RefinementCriterion::BoxZ::AboveGround && num_lo == AMREX_SPACEDIM &&
                        (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                         solverChoice.terrain_type == TerrainType::EB)) {
                        amrex::Abort("erf." + crit.name + ".box_z_coordinate = above_ground does not follow"
                                     " ImmersedForcing or EB terrain; set it to physical");
                    }
                }
                ref_criteria.push_back(crit);
            }
            else
            {
                RealBox real_box;
                int lev_for_box;

                read_box_for_refinement(ref_prefix, lev_for_box, real_box);

                if (real_box.ok())
                {
                    crit.real_box = real_box;
                    ref_criteria.push_back(crit);
                }
                else if ( (lev_for_box > 0) && (refinement_indicators[i] != "storm_tracker") )
                {
                    Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
                }
            }
        } // loop over criteria

        // With WRFInput or Metgrid, a level read from nc_init_file_<lev> is the region of its nest,
        //    given by the file or by a box without a field test, and ErrorEst tags nothing else
        //    there, so an indicator that tests a field does not refine that level
        if (solverChoice.init_type == InitType::WRFInput ||
            solverChoice.init_type == InitType::Metgrid) {
            for (const auto& crit : ref_criteria) {
                if (crit.conditions.empty()) { continue; }
                for (int lev = 1; lev <= std::min(crit.max_level, max_level); ++lev) {
                    if (num_files_at_level[lev] > 0) {
                        amrex::Print() << "WARNING: refinement indicator " << crit.name << " tests a field,"
                                       << " so it does not refine level " << lev << ", which is the nest"
                                       << " read from nc_init_file_" << lev << std::endl;
                    }
                }
            }
        }
    } // if max_level > 0
}
