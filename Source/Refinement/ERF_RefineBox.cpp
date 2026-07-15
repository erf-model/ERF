#include <ERF.H>

using namespace amrex;

void
ERF::read_box_for_refinement (std::string& ref_prefix, int& lev_for_box, RealBox& real_box)
{
    ParmParse ppr(ref_prefix);

    int num_real_lo      = ppr.countval("in_box_lo");
    int num_real_hi      = ppr.countval("in_box_hi");

    int num_indx_lo      = ppr.countval("in_box_lo_indices");
    int num_indx_hi      = ppr.countval("in_box_hi_indices");

    int num_indx_lo_crse = ppr.countval("in_box_lo_indices_crse");
    int num_indx_hi_crse = ppr.countval("in_box_hi_indices_crse");

    AMREX_ALWAYS_ASSERT( (num_real_lo      == num_real_hi)      && (num_real_lo      == 0 || num_real_lo      >= 2) );
    AMREX_ALWAYS_ASSERT( (num_indx_lo      == num_indx_hi)      && (num_indx_lo      == 0 || num_indx_lo      >= 2) );
    AMREX_ALWAYS_ASSERT( (num_indx_lo_crse == num_indx_hi_crse) && (num_indx_lo_crse == 0 || num_indx_lo_crse >= 2) );

    // Problem low and high (in real not index space) are the same at all levels
    const Real* plo = geom[0].ProbLo();
    const Real* phi = geom[0].ProbHi();
    if ( !((num_real_lo >= AMREX_SPACEDIM-1 && num_indx_lo == 0 && num_indx_lo_crse == 0) ||
           (num_indx_lo >= AMREX_SPACEDIM-1 && num_real_lo == 0 && num_indx_lo_crse == 0) ||
           (num_indx_lo ==              0   && num_real_lo == 0 && num_indx_lo_crse == 0) ||
           (num_indx_lo_crse >= AMREX_SPACEDIM-1 && num_real_lo == 0 && num_indx_lo == 0)
        ) )
    {
        amrex::Abort("Must only specify box for refinement using real OR index space with fine/coarse grid indices");
    }

    // Clear stale data
    lev_for_box = max_level;
    if (num_real_lo > 0) {
        ppr.query("max_level",lev_for_box);
    } else if (num_indx_lo > 0 || num_indx_lo_crse > 0) {
        ppr.get("max_level",lev_for_box);
    }
    num_boxes_at_level[lev_for_box] = 0;
    boxes_at_level[lev_for_box].clear();

    if (num_real_lo > 0) {

        std::vector<Real> rbox_lo(3), rbox_hi(3);
        if (lev_for_box > 0 && lev_for_box <= max_level)
        {
            if (n_error_buf[0] != IntVect::TheZeroVector()) {
                amrex::Abort("Don't use n_error_buf > 0 when setting the box explicitly");
            }

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

            const Box& domain = geom[lev_for_box].Domain();

            real_box = RealBox(&(rbox_lo[0]),&(rbox_hi[0]));

            Print() << "Realbox read in and intersected laterally with domain is " << real_box << std::endl;

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

                // Need to update real_box because tagging is based on
                // the initial _un_deformed grid
                real_box = RealBox(plo[0]+ ilo   *dx[0], plo[1]+ jlo   *dx[1], plo[2]+ klo   *dx[2],
                                   plo[0]+(ihi+1)*dx[0], plo[1]+(jhi+1)*dx[1], plo[2]+(khi+1)*dx[2]);
            } else {
                klo = static_cast<int>((rbox_lo[2] - plo[2])/dx[2]);
                khi = static_cast<int>((rbox_hi[2] - plo[2])/dx[2]-1);
            }

            // Snap box indices to ref_ratio alignment (round lo down, hi up)
            {
                const auto& rr = ref_ratio[lev_for_box-1];
                auto snap_lo = [](int idx, int r) { return idx - (idx % r + r) % r; };
                auto snap_hi = [](int idx_p1, int r) { // idx_p1 = ihi+1
                    int rem = idx_p1 % r;
                    return (rem == 0) ? idx_p1 - 1 : idx_p1 + (r - rem) - 1;
                };
                int ilo_old = ilo, jlo_old = jlo, klo_old = klo;
                int ihi_old = ihi, jhi_old = jhi, khi_old = khi;
                ilo = snap_lo(ilo, rr[0]);
                jlo = snap_lo(jlo, rr[1]);
                klo = snap_lo(klo, rr[2]);
                ihi = snap_hi(ihi+1, rr[0]);
                jhi = snap_hi(jhi+1, rr[1]);
                khi = snap_hi(khi+1, rr[2]);
                if (ilo != ilo_old || ihi != ihi_old ||
                    jlo != jlo_old || jhi != jhi_old ||
                    klo != klo_old || khi != khi_old) {
                    amrex::Print() << "Refinement box indices snapped to ref_ratio alignment:\n"
                                   << "  ilo: " << ilo_old << " -> " << ilo
                                   << "  ihi: " << ihi_old << " -> " << ihi
                                   << "  jlo: " << jlo_old << " -> " << jlo
                                   << "  jhi: " << jhi_old << " -> " << jhi
                                   << "  klo: " << klo_old << " -> " << klo
                                   << "  khi: " << khi_old << " -> " << khi << "\n";
                }
            }

            Box bx(IntVect(ilo,jlo,klo),IntVect(ihi,jhi,khi));

            bool using_pbl = (solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYJ      ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYNN25   ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYNNEDMF ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::YSU      ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MRF);

            if ( using_pbl && ( (rbox_lo[2] > plo[2]) || (rbox_hi[2] < phi[2]) ) ) {
                amrex::Print() << "PBL models need refinement boxes that go from the bottom to the top of the domain for calculation of PBLH" << std::endl;
                amrex::Print() << "Please set in_box_lo to geometry.prob_lo in z and in_box_hi to geometry.prob_hi in z and try again" << std::endl;
                amrex::Abort();
            }

            boxes_at_level[lev_for_box].push_back(bx);
            Print() << "Saving in 'boxes at level' as " << bx << std::endl;
        } // (lev_for_box > 0 && lev_for_box <= max_level)

        if (solverChoice.init_type == InitType::WRFInput) {
            amrex::Print() << "Size of num_boxes " << num_boxes_at_level.size() << std::endl;
            amrex::Print() << "Size of num_files " << num_files_at_level.size() << std::endl;
            amrex::Print() << "Consider lev_for_box = " << lev_for_box << std::endl;
            amrex::Print() << "Number of boxes at level  " << num_boxes_at_level[lev_for_box] << std::endl;
            amrex::Print() << "Number of available files " << num_files_at_level[lev_for_box] << std::endl;
            if (num_boxes_at_level[lev_for_box] != num_files_at_level[lev_for_box]) {
                amrex::Print() << "Will need to rely on refinement criteria from inputs file" << std::endl;
            }
        }

    } else if (num_indx_lo > 0) {

        std::vector<int> box_lo(3), box_hi(3);
        if (lev_for_box > 0 && lev_for_box <= max_level)
        {
            if (n_error_buf[0] != IntVect::TheZeroVector()) {
                amrex::Abort("Don't use n_error_buf > 0 when setting the box explicitly");
            }

            ppr.getarr("in_box_lo_indices",box_lo,0,num_indx_lo);
            ppr.getarr("in_box_hi_indices",box_hi,0,num_indx_hi);

            if (num_indx_lo < AMREX_SPACEDIM) {
                box_lo[2] = geom[lev_for_box].Domain().smallEnd(2);
                box_hi[2] = geom[lev_for_box].Domain().bigEnd(2);
            }

            Box bx(IntVect(box_lo[0],box_lo[1],box_lo[2]),IntVect(box_hi[0],box_hi[1],box_hi[2]));
            const Box& domain = geom[lev_for_box].Domain();

            if (!domain.contains(bx)) {
                amrex::Print() << "\n";
                amrex::Print() << "Box specified       is " << bx << std::endl;
                amrex::Print() << "But domain at level is " << domain << std::endl;
                amrex::Error("Specified box doesn't fit in the domain");
            }

            const auto* dx  = geom[lev_for_box].CellSize();
            real_box = RealBox(plo[0]+ box_lo[0]   *dx[0], plo[1]+ box_lo[1]   *dx[1], plo[2]+ box_lo[2]   *dx[2],
                               plo[0]+(box_hi[0]+1)*dx[0], plo[1]+(box_hi[1]+1)*dx[1], plo[2]+(box_hi[2]+1)*dx[2]);

            Print() << "Reading " << bx << " at level " << lev_for_box << std::endl;
            num_boxes_at_level[lev_for_box] += 1;

            // Snap box indices to ref_ratio alignment (round lo down, hi up)
            {
                const auto& rr = ref_ratio[lev_for_box-1];
                auto snap_lo_fn = [](int idx, int r) { return idx - (idx % r + r) % r; };
                auto snap_hi_fn = [](int idx_p1, int r) {
                    int rem = idx_p1 % r;
                    return (rem == 0) ? idx_p1 - 1 : idx_p1 + (r - rem) - 1;
                };
                int lo_old[3] = {box_lo[0], box_lo[1], box_lo[2]};
                int hi_old[3] = {box_hi[0], box_hi[1], box_hi[2]};
                box_lo[0] = snap_lo_fn(box_lo[0], rr[0]);
                box_lo[1] = snap_lo_fn(box_lo[1], rr[1]);
                box_lo[2] = snap_lo_fn(box_lo[2], rr[2]);
                box_hi[0] = snap_hi_fn(box_hi[0]+1, rr[0]);
                box_hi[1] = snap_hi_fn(box_hi[1]+1, rr[1]);
                box_hi[2] = snap_hi_fn(box_hi[2]+1, rr[2]);
                if (box_lo[0] != lo_old[0] || box_hi[0] != hi_old[0] ||
                    box_lo[1] != lo_old[1] || box_hi[1] != hi_old[1] ||
                    box_lo[2] != lo_old[2] || box_hi[2] != hi_old[2]) {
                    amrex::Print() << "Refinement box indices snapped to ref_ratio alignment:\n"
                                   << "  ilo: " << lo_old[0] << " -> " << box_lo[0]
                                   << "  ihi: " << hi_old[0] << " -> " << box_hi[0]
                                   << "  jlo: " << lo_old[1] << " -> " << box_lo[1]
                                   << "  jhi: " << hi_old[1] << " -> " << box_hi[1]
                                   << "  klo: " << lo_old[2] << " -> " << box_lo[2]
                                   << "  khi: " << hi_old[2] << " -> " << box_hi[2] << "\n";
                }
                bx = Box(IntVect(box_lo[0],box_lo[1],box_lo[2]),
                         IntVect(box_hi[0],box_hi[1],box_hi[2]));
            }

            bool using_pbl = (solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYJ      ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYNN25   ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYNNEDMF ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::YSU      ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MRF);

            if ( using_pbl && ( (box_lo[2] > 0) || (box_hi[2] < domain.bigEnd(2)) ) ) {
                amrex::Print() << "PBL models need refinement boxes that go from the bottom to the top of the domain for calculation of PBLH" << std::endl;
                amrex::Print() << "Please set in_box_lo_indices to 0 in z and in_box_hi_indices to amr.n_cell-1 in z and try again" << std::endl;
                amrex::Abort();
            }

            boxes_at_level[lev_for_box].push_back(bx);
            Print() << "Saving in 'boxes at level' as " << bx << std::endl;
        } // lev

        if (solverChoice.init_type == InitType::WRFInput) {
            if ( (num_files_at_level[lev_for_box] > 0) &&
                 (num_boxes_at_level[lev_for_box] != num_files_at_level[lev_for_box]) ) {
                amrex::Error("Number of boxes doesn't match number of input files");

            }
        }

    } else if (num_indx_lo_crse > 0) {

        std::vector<int> box_lo(3), box_hi(3);
        if (lev_for_box > 0 && lev_for_box <= max_level)
        {
            if (n_error_buf[0] != IntVect::TheZeroVector()) {
                amrex::Abort("Don't use n_error_buf > 0 when setting the box explicitly");
            }

            ppr.getarr("in_box_lo_indices_crse",box_lo,0,num_indx_lo_crse);
            ppr.getarr("in_box_hi_indices_crse",box_hi,0,num_indx_hi_crse);

            if (num_indx_lo_crse < AMREX_SPACEDIM) {
                box_lo[2] = geom[lev_for_box-1].Domain().smallEnd(2);
                box_hi[2] = geom[lev_for_box-1].Domain().bigEnd(2);
            }

            Box bx(IntVect(box_lo[0],box_lo[1],box_lo[2]),IntVect(box_hi[0],box_hi[1],box_hi[2]));

            if (!geom[lev_for_box-1].Domain().contains(bx)) {
                amrex::Print() << "\n";
                amrex::Print() << "(Coarse) Box specified       is " << bx << std::endl;
                amrex::Print() << "But (coarse) domain at level is " << geom[lev_for_box-1].Domain() << std::endl;
                amrex::Error("Specified box doesn't fit in the domain");
            }

            bx.refine(ref_ratio[lev_for_box-1]);

            const auto* dx  = geom[lev_for_box-1].CellSize();

            real_box = RealBox(plo[0]+ box_lo[0]   *dx[0], plo[1]+ box_lo[1]   *dx[1], plo[2]+ box_lo[2]   *dx[2],
                               plo[0]+(box_hi[0]+1)*dx[0], plo[1]+(box_hi[1]+1)*dx[1], plo[2]+(box_hi[2]+1)*dx[2]);

            Print() << "Reading " << bx << " at level " << lev_for_box << std::endl;
            num_boxes_at_level[lev_for_box] += 1;
            bool using_pbl = (solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYJ      ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYNN25   ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MYNNEDMF ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::YSU      ||
                              solverChoice.turbChoice[lev_for_box].pbl_type == PBLType::MRF);

            const Box& domain = geom[lev_for_box].Domain();
            if ( using_pbl && ( (box_lo[2] > 0) || (box_hi[2] < domain.bigEnd(2)) ) ) {
                amrex::Print() << "PBL models need refinement boxes that go from the bottom to the top of the domain for calculation of PBLH" << std::endl;
                amrex::Print() << "Please set in_box_lo_indices_crse to 0 in z and in_box_hi_indices_crse  to amr.n_cell-1 in z and try again" << std::endl;
                amrex::Abort();
            }

            boxes_at_level[lev_for_box].push_back(bx);
            Print() << "Saving in 'boxes at level' as " << bx << std::endl;
        } // lev

        if (solverChoice.init_type == InitType::WRFInput) {
            if ( (num_files_at_level[lev_for_box] > 0) &&
                 (num_boxes_at_level[lev_for_box] != num_files_at_level[lev_for_box]) ) {
                amrex::Error("Number of boxes doesn't match number of input files");

            }
        }
    }
}

void
ERF::update_box_for_refinement (std::string& ref_prefix, int& lev_for_box, RealBox& real_box, const Real time)
{
    ParmParse ppr(ref_prefix);

    Vector<Real> move_start_time, move_stop_time;
    int ni = ppr.queryarr("move_start_time", move_start_time);
    int nj = ppr.queryarr("move_stop_time" ,  move_stop_time);
    if (ni != nj) {
        amrex::Print() << "Must be same number of start times as stop times for moving grids" << std::endl;
        amrex::Abort();
    }
    for (int i = 0; i < ni; i++) {
        if (move_stop_time[i] <= move_start_time[i]) {
            amrex::Print() << "start time for interval " << i << " is " << move_start_time[i] << std::endl;
            amrex::Print() << "stop  time for interval " << i << " is " << move_stop_time[i]  << std::endl;
            amrex::Abort("moving grid: stop time must be greater than start time");
        }
    }
    for (int i = 1; i < ni; i++) {
        if (move_start_time[i] < move_stop_time[i-1]) {
            amrex::Print() << "start time for interval " << i   << " is " << move_start_time[i] << std::endl;
            amrex::Print() << "stop  time for interval " << i-1 << " is " << move_stop_time[i-1]  << std::endl;
            amrex::Abort("moving grid: stop time must be less than start time of the next interval");
        }
    }

    Vector<Real> move_speed_x, move_speed_y;
    int ni2 = ppr.queryarr("move_speed_x", move_speed_x);
    int nj2 = ppr.queryarr("move_speed_y", move_speed_y);
    if (ni2 != nj2 ) {
        amrex::Print() << "Must be same number of speeds in x- and y-directions" << std::endl;
        amrex::Abort();
    }
    if (ni != ni2 ) {
        amrex::Print() << "Must be same number of speeds as time intervals" << std::endl;
        amrex::Abort();
    }

    Real offset_x = zero;
    Real offset_y = zero;

    for (int i = 0; i < ni; i++) {
        if (time > move_start_time[i]) {
            offset_x += move_speed_x[i] * (std::min(time,move_stop_time[i]) - move_start_time[i]);
            offset_y += move_speed_y[i] * (std::min(time,move_stop_time[i]) - move_start_time[i]);
        }
    }

    RealBox orig_real_box;
    read_box_for_refinement (ref_prefix, lev_for_box, orig_real_box);

    Real xlo = orig_real_box.lo(0) + offset_x; Real ylo = orig_real_box.lo(1) + offset_y;
    Real xhi = orig_real_box.hi(0) + offset_x; Real yhi = orig_real_box.hi(1) + offset_y;
    Real zlo = orig_real_box.lo(2); Real zhi = orig_real_box.hi(2);

    real_box.setLo(RealVect(xlo,ylo,zlo));
    real_box.setHi(RealVect(xhi,yhi,zhi));
}
