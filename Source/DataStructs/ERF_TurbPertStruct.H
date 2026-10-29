#ifndef ERF_TURB_PERT_STRUCT_H_
#define ERF_TURB_PERT_STRUCT_H_

#include <ERF_DataStruct.H>
#include "ERF_Constants.H"
#include <AMReX_MultiFabUtil.H>
#include <ERF_TileNoZ.H>
#include <ERF_HashRNG.H>
#include <cstdint>
/**
 * Container holding quantities related to turbulent perturbation parameters
 */

/*  The general rule of thumb is to create a perturbation box size of 1/8th of the boundary layer height.
    The boundary layer height can't be the height of the domain. The length and width of the box should
    be twice the height of the box. If meandering flow is present, the width of the box should be take
    the angle of the inflow into consideration.
*/

AMREX_ENUM(PerturbationType,
    Source, Direct, CPM, CPM_W, None
);

// One lateral edge of a rectangular nest. Faces are xlo, ylo, xhi, yhi.
struct AutomaticPerturbationFace {
    amrex::Box region;
    int face = 0;
    amrex::Real u = zero;
    amrex::Real v = zero;
    amrex::Real Ug = zero; // CPM wind scale sampled above the PBL
    amrex::Real height = zero; // perturbation ceiling, meters above local terrain

    bool inflow () const {
        const amrex::Real normal = (face % 2 == 0) ? u : v;
        return height > zero && ((face < 2) ? normal > zero : normal < zero);
    }
};

/**
 * @brief Runtime state and operations for turbulent perturbation forcing.
 */
struct TurbulentPerturbation {

  public:

    /**
     * @brief Destroy turbulent perturbation state.
     */
    ~TurbulentPerturbation () {}

    /**
     * @brief Store the integer perturbation-type selector for one AMR level.
     * @param lev AMR level index.
     * @param pert_type Perturbation type selected for this level.
     * @param max_level Maximum AMR level configured for the run.
     */
    void init_tpi_type (const int lev,
                        const PerturbationType& pert_type,
                        const int max_level)
    {
        if (pt_type.size() < max_level + 1) {
            pt_type.resize(max_level + 1, -1);
        }

        if (pert_type == PerturbationType::Source) {
            pt_type[lev] = 0;
        } else if (pert_type == PerturbationType::Direct) {
            pt_type[lev] = 1;
        } else if (pert_type == PerturbationType::CPM) {
            pt_type[lev] = 2;
        } else if (pert_type == PerturbationType::CPM_W) {
            pt_type[lev] = 3;
        } else {
            pt_type[lev] = -1;
        }
    }

    /**
     * @brief Initialize turbulent perturbation regions and per-box storage.
     * @param lev AMR level index.
     * @param subdomains_lev Rectangular level subdomains used to define perturbation regions.
     * @param dx Cell spacing for this level.
     * @param ba State BoxArray used for perturbation-cell storage.
     * @param dm DistributionMapping for perturbation-cell storage.
     * @param ngrow_state Number of ghost cells for perturbation-cell storage.
     * @param pp_prefix ParmParse prefix for the ERF input namespace.
     * @param refRatio Refinement ratios between AMR levels.
     * @param max_level Maximum AMR level configured for the run.
     * @param automatic Select runtime inflow/height masks driven by the parent.
     */
    void init_tpi (const int lev,
                   const amrex::Vector<amrex::BoxArray>& subdomains_lev,
                   const amrex::GpuArray<amrex::Real,3> dx,
                   const amrex::BoxArray& ba,
                   const amrex::DistributionMapping& dm,
                   const int ngrow_state,
                   std::string pp_prefix,
                   const amrex::Vector<amrex::IntVect> refRatio,
                   const int max_level,
                   const bool automatic = false)

    {
        ref_ratio = refRatio;
        auto_mode = automatic;

        amrex::ParmParse pp(pp_prefix);

        uniform_mode = true;
        pp.queryAdd("perturbation_uniform_refresh", uniform_mode);
        uniform_height = amrex::Real(20);
        pp.queryAdd("perturbation_refresh_height", uniform_height);
        uniform_wind_period = amrex::Real(600);
        pp.queryAdd("perturbation_refresh_wind_period", uniform_wind_period);
        if (uniform_refresh_enabled(lev) &&
            (!(uniform_height > zero) || !std::isfinite(uniform_height) ||
             !(uniform_wind_period > zero) || !std::isfinite(uniform_wind_period))) {
            amrex::Abort("Uniform CPM refresh requires positive, finite "
                         "perturbation_refresh_height and perturbation_refresh_wind_period");
        }

        // Reading inputs, and placing assertion for the perturbation inflow to work
        pp.getarr("perturbation_box_dims",tpi_boxDim);
        if (auto_mode) {
            tpi_direction = {1, 1, 0, 1, 1, 0};
        } else {
            pp.getarr("perturbation_direction",tpi_direction);
        }
        if (tpi_boxDim.size() != 3) {
            amrex::Abort("perturbation_box_dims must have 3 entries (Nx Ny Nz)");
        }
        if (tpi_direction.size() != 6) {
            amrex::Abort("perturbation_direction must have 6 entries (west south bottom east north top)");
        }
        pp.get("perturbation_layers",tpi_layers);
        pp.get("perturbation_offset",tpi_offset);

        tpi_nonDim = zero;
        pp.queryAdd("perturbation_nondimensional",tpi_nonDim);

        tpi_Tinf = amrex::Real(300.);
        pp.queryAdd("perturbation_T_infinity",tpi_Tinf);

        tpi_Ti = zero;
        pp.queryAdd("perturbation_T_intensity",tpi_Ti);

        input_Ug = zero;
        pp.queryAdd("perturbation_Ug",input_Ug);

        input_w_amp = zero;
        pp.queryAdd("perturbation_w_amp",input_w_amp);

        // Seed of the perturbation's random numbers, read the way erf.deterministic_ic_pert reads
        // it. The numbers are hashed from the global cell index (or the perturbation-box index) and
        // the box's refresh count, so they do not depend on the domain decomposition (#4317).
        {
            int fix_random_seed = 0;
            long random_seed = -1;
            amrex::ParmParse pp_erf("erf");
            pp_erf.query("fix_random_seed", fix_random_seed);
            pp_erf.query("random_seed", random_seed);
            tpi_seed = fix_random_seed ? std::uint64_t(1024)
                     : ((random_seed >= 0) ? static_cast<std::uint64_t>(random_seed) : std::uint64_t(0));
        }

        // Negative sentinels rather than the queryAdd return values: k = 0 is a legal
        // bound, so the natural default cannot double as "unset", and the return value
        // only reports whether the key existed before the call.  init_tpi runs once per
        // level, so testing the return would make every level above 0 see "specified"
        // with the value level 0 planted, skipping the valid_box_lo/hi defaults below.
        perturbation_klo = -1;
        pp.queryAdd("perturbation_klo", perturbation_klo);
        const bool have_klo = (perturbation_klo >= 0);

        perturbation_khi = -1;
        pp.queryAdd("perturbation_khi", perturbation_khi);
        const bool have_khi = (perturbation_khi >= 0);

        // Check variables message
        if (tpi_offset < 0) { amrex::Abort("Please provide a valid inflow cell offset value for perturbation region (ie. 0-5)"); }
        if (tpi_layers < 0) { amrex::Abort("Please provide a valid perturbation layer value (ie. 3-5)"); }
        if ((auto_mode || uniform_refresh_enabled(lev)) && tpi_layers == 0) {
            amrex::Abort("Automatic perturbations and uniform CPM refresh require perturbation_layers > 0");
        }
        if (tpi_nonDim < zero) { amrex::Abort("Please provide a valid nondimensional number (ie. Ri = amrex::Real(0.042))"); }
        for (int i = 0; i < tpi_boxDim.size(); i++) {
            if (tpi_boxDim[i] < 3) { amrex::Abort("Please provide valid dimensions for perturbation boxes."); }
        }
        if (input_Ug   < zero) { amrex::Abort("Please provide a valid geostrophic wind speed (ie. Ug = amrex::Real(10.0) m/s)"); }
        if (tpi_Tinf   < zero) { amrex::Abort("Please provide a valid ambient temperature value (ie. T_0 = T_infty)"); }
        if (tpi_Ti     < zero) { amrex::Abort("Please provide a valid temperature intensity value (ie. 0-one)"); }

        // Create a temporary box list to accumulate all the perturbation regions after box modification
        amrex::BoxList tmp_bl;

        // boxSize for individual boxes
        amrex::IntVect boxSize(tpi_boxDim[0],tpi_boxDim[1],tpi_boxDim[2]);

        if (tpi_direction[2] || tpi_direction[5]) { amrex::Abort("Currently not supporting z-direction flow perturbation"); }

        for (int isub = 0; isub < subdomains_lev.size(); ++isub) {
            const amrex::BoxArray& subdomain = subdomains_lev[isub];
            amrex::Box subdomain_box(subdomain.minimalBox());

            if (subdomain_box.numPts() != subdomain.numPts()) {
                amrex::Abort("Turbulent perturbations require rectangular subdomains. "
                             "Level " + std::to_string(lev) +
                             ", subdomain " + std::to_string(isub) +
                             " is not a rectangular region fully covered by grids.");
            }

            const amrex::IntVect& valid_box_lo = subdomain_box.smallEnd();
            const amrex::IntVect& valid_box_hi = subdomain_box.bigEnd();

            if (auto_mode) {
                // Allocate all possible edges once. A mask chooses inflow and height
                // at runtime, including corners when either adjoining edge is inflow.
                for (int face = 0; face < 4; ++face) {
                    const amrex::Box region = automatic_face_box(subdomain_box, face);
                    if (region.ok()) { tmp_bl.push_back(region); }
                }
                continue;
            }

            // default perturbation region to be the entire z extent if the user does not specify bounds
            if (!have_klo) { perturbation_klo = valid_box_lo[2]; }
            if (!have_khi) { perturbation_khi = valid_box_hi[2]; }

            // Creating perturbation regions and initializing with generic size.
            amrex::Box lo_x_bx(amrex::IntVect(0), amrex::IntVect(1), amrex::IntVect(0));
            amrex::Box hi_x_bx(amrex::IntVect(0), amrex::IntVect(1), amrex::IntVect(0));
            amrex::Box lo_y_bx(amrex::IntVect(0), amrex::IntVect(1), amrex::IntVect(0));
            amrex::Box hi_y_bx(amrex::IntVect(0), amrex::IntVect(1), amrex::IntVect(0));

            // Starting logic to set the size of the perturbation region(s)
            //amrex::PrintToFile("BoxPerturbationOutput") << "Setting perturbation region in:";
            // ***** X-direction perturbation *****
            if (tpi_direction[0]) { // West
                lo_x_bx.setSmall(amrex::IntVect(valid_box_lo[0]+tpi_offset, valid_box_lo[1]+tpi_direction[1]*tpi_offset, perturbation_klo));
                lo_x_bx.setBig  (amrex::IntVect(valid_box_lo[0]+(tpi_layers*tpi_boxDim[0]-1)+tpi_offset, valid_box_hi[1]-(tpi_direction[4]*tpi_offset), perturbation_khi));
                amrex::PrintToFile("BoxPerturbationOutput") << " West face";
            }

            if (tpi_direction[3]) { // East
                hi_x_bx.setSmall(amrex::IntVect(valid_box_hi[0]-((tpi_layers*tpi_boxDim[0]-1)+tpi_offset), valid_box_lo[1]+tpi_direction[1]*tpi_offset, perturbation_klo));
                hi_x_bx.setBig(amrex::IntVect(valid_box_hi[0]-tpi_offset, valid_box_hi[1]-(tpi_direction[4]*tpi_offset), perturbation_khi));
                amrex::PrintToFile("BoxPerturbationOutput") << " East face";
            }

            // ***** Y-direction Perturbation *****
            if (tpi_direction[1]) { // North
                lo_y_bx.setSmall(amrex::IntVect(valid_box_lo[0]+tpi_direction[0]*tpi_offset, valid_box_lo[1]+tpi_offset, perturbation_klo));
                lo_y_bx.setBig(amrex::IntVect(valid_box_hi[0]-tpi_direction[3]*tpi_offset, valid_box_lo[1]+((tpi_layers*tpi_boxDim[1])-1)+tpi_offset, perturbation_khi));
                amrex::PrintToFile("BoxPerturbationOutput") << " North face";
            }

            if (tpi_direction[4]) { // South
                hi_y_bx.setSmall(amrex::IntVect(valid_box_lo[0]+tpi_direction[0]*tpi_offset, valid_box_hi[1]-((tpi_layers*tpi_boxDim[1]-1)+tpi_offset), perturbation_klo));
                hi_y_bx.setBig(amrex::IntVect(valid_box_hi[0]-tpi_direction[3]*tpi_offset, valid_box_hi[1]-tpi_offset, perturbation_khi));
                amrex::PrintToFile("BoxPerturbationOutput") << " South face";
            }

            // Performing box union for intersecting perturbation regions to avoid overlapping sections (double counting at corners)
            if (tpi_direction[0] && tpi_direction[1]) {  // Reshaping South smallEnd
                amrex::Box lo_x_lo_y_u = lo_x_bx & lo_y_bx;
                lo_y_bx.setSmall(amrex::IntVect(lo_x_lo_y_u.bigEnd(0)+1, lo_x_lo_y_u.smallEnd(1), lo_x_lo_y_u.smallEnd(2)));
            }

            if (tpi_direction[3] && tpi_direction[1]) { // Reshaping South bigEnd
                 amrex::Box hi_x_lo_y_u = hi_x_bx & lo_y_bx;
                 lo_y_bx.setBig(amrex::IntVect(hi_x_lo_y_u.smallEnd(0)-1, hi_x_lo_y_u.bigEnd(1), hi_x_lo_y_u.bigEnd(2)));
            }

            if (tpi_direction[0] && tpi_direction[4]) {  // Reshaping North smallEnd
                amrex::Box lo_x_hi_y_u = lo_x_bx & hi_y_bx;
                hi_y_bx.setSmall(amrex::IntVect(lo_x_hi_y_u.bigEnd(0)+1, lo_x_hi_y_u.smallEnd(1), lo_x_hi_y_u.smallEnd(2)));
            }

            if (tpi_direction[3] && tpi_direction[4]) { // Reshaping North bigEnd
                 amrex::Box hi_x_hi_y_u = hi_x_bx & hi_y_bx;
                 hi_y_bx.setBig(amrex::IntVect(hi_x_hi_y_u.smallEnd(0)-1, hi_x_hi_y_u.bigEnd(1), hi_x_hi_y_u.bigEnd(2)));
            }

            // Creating structure box array for conserved quantity
            if (tpi_direction[0]) { tmp_bl.push_back(lo_x_bx); }
            if (tpi_direction[1]) { tmp_bl.push_back(lo_y_bx); }
            if (tpi_direction[3]) { tmp_bl.push_back(hi_x_bx); }
            if (tpi_direction[4]) { tmp_bl.push_back(hi_y_bx); }
        }

        //amrex::PrintToFile("BoxPerturbationOutput") << "\nBoxList: " << tmp_bl << "\n";
        amrex::BoxArray tmp_ba(tmp_bl);
        if (auto_mode) { tmp_ba.removeOverlap(); }
        tmp_ba.maxSize(boxSize);
        if (auto_mode && tmp_ba.empty()) {
            amrex::Abort("Automatic perturbation strips are empty; reduce perturbation_offset");
        }

        const int num_levels = max_level + 1;
        if (pb_ba.size() < num_levels) {
            pb_ba.resize(num_levels);
            pb_mag.resize(num_levels);
            pb_Ug.resize(num_levels);
            pb_dir.resize(num_levels);
            pb_netZero.resize(num_levels);
            pb_interval.resize(num_levels);
            pb_local_etime.resize(num_levels);
            pb_amp.resize(num_levels);
            pb_nupdate.resize(num_levels);
            pb_cell.resize(num_levels);
            tpi_Lpb.resize(num_levels);
            tpi_Wpb.resize(num_levels);
            tpi_Hpb.resize(num_levels);
            tpi_lref.resize(num_levels);
            auto_mask.resize(num_levels);
            auto_next_update.resize(num_levels, -1.0);
            auto_ready.resize(num_levels, 0);
            uniform_clock.resize(num_levels);
        }
        uniform_clock[lev] = UniformClock{};

        pb_ba[lev] = tmp_ba;

        // Initializing mean magnitude and direction vectors
        pb_mag[lev].resize(pb_ba[lev].size(), zero);
        pb_Ug[lev].resize(pb_ba[lev].size(), zero);
        pb_dir[lev].resize(pb_ba[lev].size(), zero);
        pb_netZero[lev].resize(pb_ba[lev].size(), zero);

        // Set size of vector and initialize
        pb_interval[lev].resize(pb_ba[lev].size(), -one);
        pb_local_etime[lev].resize(pb_ba[lev].size(), zero);
        pb_amp[lev].resize(pb_ba[lev].size(), zero);
        pb_nupdate[lev].resize(pb_ba[lev].size(), 0);

        // Creating data array for perturbation amplitude storage
        if (pt_type[lev] == 3) { // CPM_W converts to the k-face ba.
            pb_cell[lev].define(convert(ba, amrex::IntVect(0,0,1)), dm, 1, ngrow_state);
        } else {
            pb_cell[lev].define(ba, dm, 1, ngrow_state); // this is the only place ba is used. Maybe we can print here to determine what's valid...
        }
        pb_cell[lev].setVal(0.);
        if (auto_mode) {
            auto_mask[lev].define(pb_cell[lev].boxArray(), dm, 1, 0);
            auto_mask[lev].setVal(zero);
            auto_next_update[lev] = -1.0;
            auto_ready[lev] = 0;
            // Regridding can leave vectors at the same size; resize alone would retain
            // clocks/amplitudes belonging to a different set of boxes.
            pb_interval[lev].assign(pb_ba[lev].size(), -one);
            pb_local_etime[lev].assign(pb_ba[lev].size(), zero);
            pb_nupdate[lev].assign(pb_ba[lev].size(), 0);
        }

        // Computing perturbation reference length
        tpi_Lpb[lev] = tpi_boxDim[0]*dx[0];
        tpi_Wpb[lev] = tpi_boxDim[1]*dx[1];
        tpi_Hpb[lev] = tpi_boxDim[2]*dx[2];
        tpi_lref[lev] = std::sqrt(tpi_Lpb[lev]*tpi_Lpb[lev] + tpi_Wpb[lev]*tpi_Wpb[lev]);

        tpi_pert_adjust = zero;
        tpi_net_buoyant = zero;

        /*
        // Function check point message
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_box_dims: "
                                                    << tpi_boxDim[0] << " "
                                                    << tpi_boxDim[1] << " "
                                                    << tpi_boxDim[2] << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_direction: "
                                                    << tpi_direction[0] << " "
                                                    << tpi_direction[1] << " "
                                                    << tpi_direction[2] << "\n\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_layers: " << tpi_layers << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_offset: " << tpi_offset << "\n\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_nondimensional: " << tpi_nonDim << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_T_infinity: " << tpi_Tinf << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_T_intensity: " << tpi_Ti << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "Reference length per box = " << tpi_lref[lev] << "\n\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "Turbulent perturbation BoxArray:\n" << pb_ba[lev] << "\n";
        */
    }


    // Geometry of a candidate edge strip, clipped to the nest interior.
    amrex::Box automatic_face_box (const amrex::Box& subdomain, int face) const
    {
        amrex::Box region = subdomain;
        region.grow(0, -tpi_offset);
        region.grow(1, -tpi_offset);
        const int dir = face % 2;
        const int width = tpi_layers * tpi_boxDim[dir];
        if (face < 2) {
            region.setBig(dir, std::min(region.bigEnd(dir), region.smallEnd(dir)+width-1));
        } else {
            region.setSmall(dir, std::max(region.smallEnd(dir), region.bigEnd(dir)-width+1));
        }
        return region;
    }

    bool automatic_update_due (int lev, double time) const
    {
        return !auto_ready[lev] || time >= auto_next_update[lev];
    }

    bool uniform_refresh_enabled (int lev) const
    {
        return uniform_mode && (pt_type[lev] == 2 || pt_type[lev] == 3);
    }

    bool uniform_wind_update_due (int lev, double time) const
    {
        return time >= uniform_clock[lev].next_wind_update;
    }

    amrex::Real uniform_sampling_height () const { return uniform_height; }
    amrex::Real uniform_sampling_period () const { return uniform_wind_period; }

    bool uniform_sampling_face (int face) const
    {
        return tpi_direction[face < 2 ? face : face+1] != 0;
    }

    // speed is max(|mean u|, |mean v|) over the sampled edges. This is the
    // existing CPM speed*cos(angle) correction, with one clock for the level.
    void set_uniform_advection (int lev, double time, amrex::Real speed, double wind_period)
    {
        AMREX_ALWAYS_ASSERT(uniform_refresh_enabled(lev));
        AMREX_ALWAYS_ASSERT(speed >= zero && std::isfinite(speed));
        auto& clock = uniform_clock[lev];
        clock.interval = speed > zero
            ? double(tpi_layers) * double(std::min(tpi_Lpb[lev], tpi_Wpb[lev])) / double(speed)
            : 0.0;
        clock.next_refresh = time;
        clock.next_wind_update = time + wind_period;
        clock.apply = false;
    }

    double uniform_refresh_interval (int lev) const { return uniform_clock[lev].interval; }

    // Called once per level step. No box scans or field clears between events;
    // apply_tpi also checks this step's flag so stored increments cannot repeat.
    bool calc_tpi_uniform_update (int lev, double time)
    {
        AMREX_ALWAYS_ASSERT(uniform_refresh_enabled(lev));
        auto& clock = uniform_clock[lev];
        clock.apply = false;
        if ((auto_mode && !auto_ready[lev]) || clock.interval <= 0.0 ||
            time < clock.next_refresh || pb_ba[lev].empty()) { return false; }

        clock.next_refresh = time + clock.interval;
        clock.apply = true;
        pb_cell[lev].setVal(zero);
        const auto m_ixtype = pb_cell[lev].ixType();
        for (int b = 0; b < pb_ba[lev].size(); ++b) {
            pb_interval[lev][b] = clock.interval / double(tpi_layers);
            pb_local_etime[lev][b] = zero;
            calc_tpi_amp(lev, b, pb_interval[lev][b]);
            pseudoRandomPert(b, lev, m_ixtype);
            ++pb_nupdate[lev][b];
        }
        return true;
    }

    // Refresh the parent controls without changing the perturbation boxes or
    // their random refresh counts. z_ground is a nodal bottom slab per grid box.
    void configure_automatic (int lev, double time, double period,
                              const amrex::Vector<AutomaticPerturbationFace>& faces,
                              const amrex::Geometry& geom,
                              const amrex::MultiFab* z_cc,
                              const amrex::MultiFab* z_nd,
                              const amrex::MultiFab* z_ground)
    {
        auto_mask[lev].setVal(zero);
        pb_cell[lev].setVal(zero);
        pb_mag[lev].assign(pb_ba[lev].size(), zero);
        pb_Ug[lev].assign(pb_ba[lev].size(), zero);
        pb_dir[lev].assign(pb_ba[lev].size(), zero);
        const bool on_w_faces = (pt_type[lev] == 3);
        const bool terrain = (z_cc != nullptr);
        AMREX_ALWAYS_ASSERT(!terrain || (z_nd && z_ground));
        const int klo = geom.Domain().smallEnd(2);
        const auto dx = geom.CellSizeArray();
        const auto problo = geom.ProbLoArray();
        for (const auto& face : faces) {
            if (!face.inflow()) { continue; }
            const amrex::Real speed = std::hypot(face.u, face.v);
            for (const auto& isect : pb_ba[lev].intersections(face.region)) {
                const int b = isect.first;
                pb_Ug[lev][b] = std::max(pb_Ug[lev][b], face.Ug);
                if (speed > pb_mag[lev][b]) {
                    pb_mag[lev][b] = speed;
                    pb_dir[lev][b] = std::atan2(std::abs(face.u), std::abs(face.v));
                }
            }
            const amrex::Real height = face.height;
            const amrex::Box region = amrex::convert(face.region, auto_mask[lev].ixType());
            for (amrex::MFIter mfi(auto_mask[lev]); mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.validbox() & region;
                if (!bx.ok()) { continue; }
                const auto mask = auto_mask[lev].array(mfi);
                const auto cc = terrain ? z_cc->const_array(mfi) : amrex::Array4<const amrex::Real>{};
                const auto nd = terrain ? z_nd->const_array(mfi) : amrex::Array4<const amrex::Real>{};
                const auto ground = terrain ? z_ground->const_array(mfi) : amrex::Array4<const amrex::Real>{};
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    amrex::Real z = problo[2] + (k-klo + (on_w_faces ? zero : myhalf))*dx[2];
                    amrex::Real zs = problo[2];
                    if (terrain) {
                        zs = fourth*(ground(i,j,klo) + ground(i+1,j,klo) +
                                     ground(i,j+1,klo) + ground(i+1,j+1,klo));
                        z = on_w_faces ? fourth*(nd(i,j,k) + nd(i+1,j,k) +
                                                 nd(i,j+1,k) + nd(i+1,j+1,k))
                                       : cc(i,j,k);
                    }
                    if (z > zs && z-zs <= height) { mask(i,j,k) = one; }
                });
            }
        }
        // Reconsider even boxes that previously had calm winds. They must be able
        // to start again after a later parent update changes the inflow.
        pb_interval[lev].assign(pb_ba[lev].size(), -one);
        pb_local_etime[lev].assign(pb_ba[lev].size(), zero);
        auto_ready[lev] = 1;
        auto_next_update[lev] = time + period;
    }

    /**
     * @brief Update perturbation amplitudes and intervals when each box is ready.
     * @param lev AMR level index.
     * @param dt Time step for elapsed-time accumulation.
     * @param mf_xvel X-velocity field used to compute box mean velocities.
     * @param mf_yvel Y-velocity field used to compute box mean velocities.
     * @param mf_cons Conserved state used for iteration layout.
     */
    void calc_tpi_update (const int lev,
                          const double dt,
                          amrex::MultiFab& mf_xvel,
                          amrex::MultiFab& mf_yvel,
                          amrex::MultiFab& mf_cons)
    {
        auto m_ixtype = mf_cons.boxArray().ixType(); // safety step
        if (pt_type[lev] == 3) { m_ixtype = amrex::IndexType(amrex::IntVect(0,0,1)); }

        const int nboxes = pb_ba[lev].size();
        if (nboxes == 0 || (auto_mode && !auto_ready[lev])) { return; }

        // Which boxes refresh this step. pb_interval and pb_local_etime hold the same values on
        // every rank, so every rank makes the same choice.
        amrex::Vector<int> update_box(nboxes, 0);
        for (int boxIdx = 0; boxIdx < nboxes; boxIdx++) {
            if (pt_type[lev] == 2 || pt_type[lev] == 3) {
                // for CPM, perturbation boxes/cells refresh after boxes have advected box width/length * num_layers (advective time scale)
                update_box[boxIdx] = ((pb_local_etime[lev][boxIdx] >= pb_interval[lev][boxIdx]*tpi_layers) && pb_local_etime[lev][boxIdx] != -one);
            } else {
                update_box[boxIdx] = ( pb_local_etime[lev][boxIdx] >= pb_interval[lev][boxIdx] && pb_local_etime[lev][boxIdx] != -one );
            }
        }

        // Compute mean velocity of each refreshing perturbation box
        if (!auto_mode) { calc_tpi_meanMag(lev, update_box, mf_cons, mf_xvel, mf_yvel); }

        for (int boxIdx = 0; boxIdx < nboxes; boxIdx++) {
            if ( update_box[boxIdx] ) {

                if (pb_mag[lev][boxIdx] !=zero) {
                    double interval = 0.0;
                    if (pt_type[lev] == 2 || pt_type[lev] == 3) {
                        //  Wind direction correction for angled wind
                        amrex::Real wind_direction = pb_dir[lev][boxIdx];
                        if (wind_direction > PI / 4) { wind_direction = PI / 2 - wind_direction; }
                        // CPM only cares about the side length of a box, min call maintains flexibility.
                        interval = one / std::cos(wind_direction) * std::min(tpi_Lpb[lev], tpi_Wpb[lev]) / pb_mag[lev][boxIdx];
                    } else {
                        interval = tpi_lref[lev] / pb_mag[lev][boxIdx];
                    }
                    // 10% variation
                    const amrex::Real r = erf_hash_rng::hash_uniform(boxIdx, pb_nupdate[lev][boxIdx], 0,
                                                                     hash_comp_interval, lev, tpi_seed);
                    pb_interval[lev][boxIdx] = (amrex::Real(0.9) + amrex::Real(0.2)*r) * interval;

                    // Reset local elapsed time
                    pb_local_etime[lev][boxIdx] = zero;
                } else {
                    // No mean wind through this box, so no advective time scale: it is not refreshed again
                    pb_local_etime[lev][boxIdx] = -one;
                }

                // Trigger amplitude calculation per perturbation box
                calc_tpi_amp(lev, boxIdx, pb_interval[lev][boxIdx]);

                // Trigger random amplitude storage per cell within perturbation box
                pseudoRandomPert(boxIdx, lev, m_ixtype);

                pb_nupdate[lev][boxIdx]++;

            } else {
                // set perturbation amplitudes to 0 for CPM. A little inefficient but leverages as much as existing code as possible.
                if (pt_type[lev] == 2 || pt_type[lev] == 3) { zero_amp(boxIdx, lev, m_ixtype); }

                // Increase by timestep of level 0
                if (pb_local_etime[lev][boxIdx] != -one) { pb_local_etime[lev][boxIdx] += dt; }
            } // if
        } // for

        if (pt_type[lev] < 2) { // box perturbation method only
            netZeroBuoyant(lev);
        }
    }

    /**
     * @brief Apply stored turbulent perturbations to a source or state array.
     * @param lev AMR level index.
     * @param vbx Valid box over which perturbations may be applied.
     * @param comp Component to modify, or -1 for vertical velocity perturbations.
     * @param m_ixtype Index type of the destination array.
     * @param src_arr Destination array receiving perturbation increments.
     * @param pert_cell Stored perturbation amplitudes.
     * @param scale Factor applied to the stored amplitudes (the time step for the direct
     *              perturbation, whose stored amplitudes are rates).
     */
    void apply_tpi (const int& lev,
                    const amrex::Box& vbx,                      // box union from upper level
                    const int& comp,                            // Component to modify
                    const amrex::IndexType& m_ixtype,           // IntVect type of src_arr
                    const amrex::Array4<amrex::Real>& src_arr,  // Array to apply perturbation
                    const amrex::Array4<amrex::Real const>& pert_cell,
                    const amrex::Real scale = one)
    {
        if (uniform_refresh_enabled(lev) && !uniform_clock[lev].apply) { return; }
        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            amrex::Box pbx = amrex::convert(pb_ba[lev][boxIdx], m_ixtype);
            if (pt_type[lev] == 3) { pbx.setBig(2, pbx.bigEnd(2) - 1); } // prevent double counting after converting to k-faces
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                if (comp == -1) { // vertical velocity perturbations
                    ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        src_arr(i,j,k) += scale * pert_cell(i,j,k);
                    });
                } else {
                    ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        src_arr(i,j,k,comp) += scale * pert_cell(i,j,k);

                        // For box region debug only
                        #ifdef INDEX_PERTURB
                        src_arr(i,j,k,comp) = (amrex::Real) (boxIdx + 5.);
                        #endif
                    });
                }
            }
        }
    }

    /**
     * @brief Compute the perturbation amplitude for one perturbation box.
     * @param lev AMR level index.
     * @param boxIdx Perturbation-box index on the level.
     * @param interval Current perturbation update interval [s].
     */
    void calc_tpi_amp (const int& lev,
                       const int& boxIdx,
                       const double& interval)
    {
        pb_amp[lev][boxIdx] = zero; // Safety step
        if (pt_type[lev] == 2) { // CPM
            constexpr amrex::Real cp = amrex::Real(1004); // specific heat of air [J/(kg K)]
            constexpr amrex::Real Ec = amrex::Real(0.2); // Eckert number
            const amrex::Real Ug = auto_mode ? pb_Ug[lev][boxIdx] : input_Ug;
            pb_amp[lev][boxIdx] = (Ug * Ug) / (Ec * cp);
        } else if (pt_type[lev] == 3) {
            pb_amp[lev][boxIdx] = input_w_amp;
        } else { // box perturbation
            amrex::Real Um = pb_mag[lev][boxIdx];

            // No mean wind: nothing to trip, and with tpi_Ti > 0 the amplitude below would be
            //    0/0. The box is not refreshed again, so a zero amplitude also clears the
            //    perturbation it held
            if (Um == zero) { return; }

            amrex::Real beta = one/tpi_Tinf; // Thermal expansion coefficient

            // Pseudo Random temperature the ignores scale when mechanically tripping turbulence
            amrex::Real g = CONST_GRAV;
            // get total refinement ratio on each level
            int total_ref_ratio = 1;
            for (int level = lev; level >= 1; level--) {
                total_ref_ratio *= ref_ratio[level-1][2];
            }
            // calculation needs to be scale-aware since the formulation relies on the physical size of the box
            if (tpi_Ti > zero) g = (tpi_nonDim * Um * Um) / (tpi_Ti * tpi_Hpb[lev]) * 1 / total_ref_ratio;

            // Ma and Senocak (2023) Eq. 8, solving for delta phi
            pb_amp[lev][boxIdx] = (tpi_nonDim * Um * Um) / (g * beta * tpi_Hpb[lev]) * 1 / total_ref_ratio;

            // Divide by the interval to get a forcing rate [K/s], Ma & Senocak (2023) Eq. 7.
            //    The source perturbation adds this rate as a source term; the direct
            //    perturbation adds it times the time step to the state at every step. Either
            //    way the total over one refresh interval is delta phi, even if dt changes
            pb_amp[lev][boxIdx] /= static_cast<amrex::Real>(interval);
        }
    }

    /**
     * @brief Assign pseudo-random perturbations to cells in one perturbation box.
     * @param boxIdx Perturbation-box index on the level.
     * @param lev AMR level index.
     * @param m_ixtype Index type used to align perturbation storage.
     */
    void pseudoRandomPert (const int& boxIdx,
                           const int& lev,
                           const amrex::IndexType& m_ixtype)
    {
        for (amrex::MFIter mfi(pb_cell[lev],TileNoZ()); mfi.isValid(); ++mfi) {
            amrex::Box vbx = mfi.tilebox();
            amrex::Box pbx = amrex::convert(pb_ba[lev][boxIdx], m_ixtype);
            if (pt_type[lev] == 3) { pbx.setBig(2, pbx.bigEnd(2) - 1); } // prevent double counting after converting to k-faces
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                amrex::Real amp_copy = pb_amp[lev][boxIdx];
                const amrex::Array4<amrex::Real>& pert_cell = pb_cell[lev].array(mfi);
                const bool use_mask = auto_mode;
                const auto mask = use_mask ? auto_mask[lev].const_array(mfi)
                                           : amrex::Array4<const amrex::Real>{};

                if (pt_type[lev] == 2  || pt_type[lev] == 3) { // CPM: one value per perturbation box
                    const amrex::Real rand_number_const = erf_hash_rng::hash_symmetric(
                        boxIdx, pb_nupdate[lev][boxIdx], 0, hash_comp_cpm, lev, tpi_seed);
                    ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        pert_cell(i,j,k) = (use_mask && mask(i,j,k) == zero) ? zero : rand_number_const * amp_copy;
                    });
                } else { // one value per cell, keyed by its global index and the box's refresh count
                    const std::uint64_t cell_seed = erf_hash_rng::splitmix64(
                        tpi_seed ^ erf_hash_rng::splitmix64(static_cast<std::uint64_t>(pb_nupdate[lev][boxIdx])));
                    ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        pert_cell(i,j,k) = (use_mask && mask(i,j,k) == zero) ? zero
                            : erf_hash_rng::hash_symmetric(i, j, k, hash_comp_cell, lev, cell_seed) * amp_copy;
                    });
                }
            }
        }
    }

    /**
     * @brief Reset CPM perturbation amplitudes in one perturbation box to zero.
     * @param boxIdx Perturbation-box index on the level.
     * @param lev AMR level index.
     * @param m_ixtype Index type used to align perturbation storage.
     */
    void zero_amp (const int& boxIdx,
                   const int& lev,
                   const amrex::IndexType& m_ixtype)
    {

        for (amrex::MFIter mfi(pb_cell[lev],TileNoZ()); mfi.isValid(); ++mfi) {
            amrex::Box vbx = mfi.tilebox();
            amrex::Box pbx = amrex::convert(pb_ba[lev][boxIdx], m_ixtype);
            if (pt_type[lev] == 3) { pbx.setBig(2, pbx.bigEnd(2) - 1); } // prevent double counting after converting to k-faces
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                const amrex::Array4<amrex::Real>& pert_cell = pb_cell[lev].array(mfi);
                ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    pert_cell(i,j,k) = zero;
                });
            }
        }
    }

    /**
     * @brief Remove the mean of the box-perturbation field, so the net buoyant forcing is zero.
     *
     * Each perturbation box's mean is summed over every grid box it overlaps and over all
     * ranks, so the correction does not depend on the domain decomposition.
     * @param lev AMR level index.
     */
    void netZeroBuoyant (const int& lev)
    {
        const int nboxes = pb_ba[lev].size();
        if (nboxes == 0) { return; }

        // Sum and cell count of each perturbation box, then one reduction over the ranks
        amrex::Vector<amrex::Real> sums(2*nboxes, zero);
        for (amrex::MFIter mfi(pb_cell[lev]); mfi.isValid(); ++mfi) {
            const amrex::Box& vbx = mfi.validbox();
            const amrex::Array4<const amrex::Real>& pert_cell = pb_cell[lev].const_array(mfi);
            for (const auto& isect : pb_ba[lev].intersections(vbx)) {
                sums[2*isect.first  ] += box_sum(pert_cell, isect.second);
                sums[2*isect.first+1] += auto_mode
                    ? box_sum(auto_mask[lev].const_array(mfi), isect.second)
                    : static_cast<amrex::Real>(isect.second.numPts());
            }
        }
        amrex::ParallelAllReduce::Sum(sums.data(), static_cast<int>(sums.size()),
                                      amrex::ParallelContext::CommunicatorSub());

        // Normalizing the adjustment based on how many boxes there are; each box's
        // value is already normalized by the number of cells within it
        tpi_net_buoyant = zero;
        for (int boxIdx = 0; boxIdx < nboxes; boxIdx++) {
            pb_netZero[lev][boxIdx] = (sums[2*boxIdx+1] > zero) ? sums[2*boxIdx] / sums[2*boxIdx+1] : zero;
            tpi_net_buoyant += pb_netZero[lev][boxIdx];
        }
        tpi_pert_adjust = tpi_net_buoyant / static_cast<amrex::Real>(nboxes);
        if (auto_mode) {
            amrex::Real sum = zero, count = zero;
            for (int b = 0; b < nboxes; ++b) { sum += sums[2*b]; count += sums[2*b+1]; }
            tpi_pert_adjust = count > zero ? sum / count : zero;
        }

        const amrex::Real adjust = tpi_pert_adjust;
        for (amrex::MFIter mfi(pb_cell[lev]); mfi.isValid(); ++mfi) {
            const amrex::Box& vbx = mfi.validbox();
            const amrex::Array4<amrex::Real>& pert_cell = pb_cell[lev].array(mfi);
            const bool use_mask = auto_mode;
            const auto mask = use_mask ? auto_mask[lev].const_array(mfi)
                                       : amrex::Array4<const amrex::Real>{};
            for (const auto& isect : pb_ba[lev].intersections(vbx)) {
                ParallelFor(isect.second, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (!use_mask || mask(i,j,k) > zero) { pert_cell(i,j,k) -= adjust; }
                });
            }
        }
    }

    /**
     * @brief Compute the mean horizontal velocity magnitude and direction of the refreshing perturbation boxes.
     *
     * The velocities are summed over every grid box a perturbation box overlaps, each face once,
     * and over all ranks, so every rank gets the same mean whatever the domain decomposition.
     * @param lev AMR level index.
     * @param update_box Nonzero for the perturbation boxes that refresh this step.
     * @param mf_cons Conserved state used for iteration layout.
     * @param mf_xvel X-velocity field.
     * @param mf_yvel Y-velocity field.
     */
    void calc_tpi_meanMag (const int& lev,
                           const amrex::Vector<int>& update_box,
                           amrex::MultiFab& mf_cons,
                           amrex::MultiFab& mf_xvel,
                           amrex::MultiFab& mf_yvel)
    {
        const amrex::BoxArray& m_pb_ba = pb_ba[lev];
        const int nboxes = m_pb_ba.size();

        // Per box: sum of u, number of u faces, sum of v, number of v faces
        amrex::Vector<amrex::Real> sums(4*nboxes, zero);
        for (amrex::MFIter mfi(mf_cons); mfi.isValid(); ++mfi) {
            const amrex::Box& vbx = mfi.validbox(); // cell centered
            const amrex::Array4<const amrex::Real>& xvel_arry = mf_xvel.const_array(mfi);
            const amrex::Array4<const amrex::Real>& yvel_arry = mf_yvel.const_array(mfi);
            for (const auto& isect : m_pb_ba.intersections(vbx)) {
                const int boxIdx = isect.first;
                if (!update_box[boxIdx]) { continue; }
                const amrex::Box& pbx = m_pb_ba[boxIdx];

                // The faces of the overlap's cells; the high face only where the overlap ends
                // at the perturbation box's high side, so a face shared by two grid boxes counts once
                amrex::Box ubx_u = amrex::surroundingNodes(isect.second, 0);
                if (isect.second.bigEnd(0) < pbx.bigEnd(0)) { ubx_u.growHi(0, -1); }
                amrex::Box ubx_v = amrex::surroundingNodes(isect.second, 1);
                if (isect.second.bigEnd(1) < pbx.bigEnd(1)) { ubx_v.growHi(1, -1); }

                sums[4*boxIdx  ] += box_sum(xvel_arry, ubx_u);
                sums[4*boxIdx+1] += static_cast<amrex::Real>(ubx_u.numPts());
                sums[4*boxIdx+2] += box_sum(yvel_arry, ubx_v);
                sums[4*boxIdx+3] += static_cast<amrex::Real>(ubx_v.numPts());
            }
        }
        amrex::ParallelAllReduce::Sum(sums.data(), static_cast<int>(sums.size()),
                                      amrex::ParallelContext::CommunicatorSub());

        for (int boxIdx = 0; boxIdx < nboxes; boxIdx++) {
            if (!update_box[boxIdx]) { continue; }
            const amrex::Real u_avg = (sums[4*boxIdx+1] > zero) ? sums[4*boxIdx  ] / sums[4*boxIdx+1] : zero;
            const amrex::Real v_avg = (sums[4*boxIdx+3] > zero) ? sums[4*boxIdx+2] / sums[4*boxIdx+3] : zero;
            pb_mag[lev][boxIdx] = std::sqrt(u_avg*u_avg + v_avg*v_avg);
            pb_dir[lev][boxIdx] = std::atan(std::abs(u_avg) / std::abs(v_avg+std::numeric_limits<amrex::Real>::epsilon()));
        }
    }

    /**
     * @brief Sum an array over a box.
     * @param arr Array to sum.
     * @param bx Box to sum over.
     * @return Sum of arr over bx.
     */
    static amrex::Real box_sum (const amrex::Array4<const amrex::Real>& arr, const amrex::Box& bx)
    {
        amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            return {arr(i,j,k)};
        });
        return amrex::get<0>(reduce_data.value(reduce_op));
    }

    /**
     * @brief Write perturbation debug information when debug output is enabled.
     */
    void debug (double /*time*/)
    {
        /*
        amrex::PrintToFile("BoxPerturbationOutput") << "#################### PB output at time = "
            << time << " ####################\n";
        amrex::PrintToFile("BoxPerturbationOutput") << " Using type: " <<  pt_type << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << " Net: " << tpi_net_buoyant << " Adjust : " <<  tpi_pert_adjust << "\n";
        for (int i = 0; i < pb_mag.size(); i++) {
            amrex::PrintToFile("BoxPerturbationOutput") << "[" << i
                << "] pb_Umag=" << pb_mag[i]
                << " | pb_interval=" << pb_interval[i]
                << " (" << pb_local_etime[i]
                << ") | pb_amp=" << pb_amp[i] << "\n";
        }
        amrex::PrintToFile("BoxPerturbationOutput") << "\n";
        */
    }

    amrex::Vector<int> pt_type;

    // Public data members
    amrex::Vector<amrex::BoxArray>                 pb_ba;  // PB box array
    amrex::Vector<amrex::Vector<amrex::Real>>      pb_mag; // BP mean magnitude [m/s]
    amrex::Vector<amrex::Vector<amrex::Real>>      pb_dir; // BP mean direction [deg]

    // Perturbation amplitude cell storage
    // This is after random assignment of equation (10) in Ma and Senocak 2023
    amrex::Vector<amrex::MultiFab>  pb_cell;

  private:

    struct UniformClock {
        double interval = 0.0;
        double next_refresh = -1.0;
        double next_wind_update = -1.0;
        bool apply = false;
    };
    bool uniform_mode{true};
    amrex::Real uniform_height{amrex::Real(20)};
    amrex::Real uniform_wind_period{amrex::Real(600)};
    amrex::Vector<UniformClock> uniform_clock;

    bool auto_mode{false};
    amrex::Vector<amrex::Vector<amrex::Real>> pb_Ug; // CPM amplitude scale; independent of refresh winds
    amrex::Vector<amrex::MultiFab> auto_mask;
    amrex::Vector<double> auto_next_update;
    // Integer flags avoid std::vector<bool>'s proxy references in amrex::Vector.
    amrex::Vector<int> auto_ready;

    // Private data members
    int tpi_layers;                   // Number of layers of perturbation boxes
    int tpi_offset;                   // Cells to offset the start of the perturbation region

    amrex::Vector<int> tpi_boxDim;    // Dimensions of each perturbation box
    amrex::Vector<int> tpi_direction; // Direction of the perturbation

    // Richardson Formulation
    amrex::Real tpi_nonDim;           // Richardson number
    amrex::Real tpi_Ti;               // Temperature intensity value
    amrex::Real tpi_Tinf;             // Reference temperature [K]

    // Physical dimensions
    amrex::Vector<amrex::Real> tpi_Hpb;              // PB height [m]
    amrex::Vector<amrex::Real> tpi_Lpb;              // PB length [m]
    amrex::Vector<amrex::Real> tpi_Wpb;              // PB width  [m]
    amrex::Vector<amrex::Real> tpi_lref;             // PB reference length [m]

    amrex::Real tpi_net_buoyant;      // Perturbation net-zero calculation storage
    amrex::Real tpi_pert_adjust;      // Perturbation adjust for net-zero per cell adjustment

    amrex::Vector<amrex::IntVect> ref_ratio;  // ref_ratio for multilevel run
    amrex::Real input_Ug;                     // input geostrophic wind speed to scale perturbations when using CPM
    amrex::Real input_w_amp;                  // input vertical velocity amplitude for w
    int         perturbation_klo;             // upper bound in Nz for perturbations
    int         perturbation_khi;             // lower bound in Nz for perturbations

    // Perturbation data arrays
    amrex::Vector<amrex::Vector<double>> pb_interval;         // PB update time [s]
    amrex::Vector<amrex::Vector<double>> pb_local_etime;      // PB local elapsed time [s]
    amrex::Vector<amrex::Vector<amrex::Real>> pb_amp;         // PB perturbation amplitude Ri:[K]
    amrex::Vector<amrex::Vector<amrex::Real>> pb_netZero;     // PB array used for net zero sum calculation
    amrex::Vector<amrex::Vector<int>>         pb_nupdate;     // PB number of refreshes so far

    // Random numbers
    std::uint64_t tpi_seed = 0;                // seed of the hashed random numbers
    static constexpr int hash_comp_cell     = 0; // hash salts of the three uses
    static constexpr int hash_comp_cpm      = 1;
    static constexpr int hash_comp_interval = 2;

};

/**
 * @var TurbulentPerturbation::pt_type
 * @brief Integer perturbation type for each AMR level.
 * @var TurbulentPerturbation::pb_ba
 * @brief Perturbation-box BoxArray for each AMR level.
 * @var TurbulentPerturbation::pb_mag
 * @brief Mean velocity magnitude for each perturbation box [m/s].
 * @var TurbulentPerturbation::pb_dir
 * @brief Mean velocity direction for each perturbation box.
 * @var TurbulentPerturbation::pb_cell
 * @brief Per-cell perturbation amplitude storage.
 * @var TurbulentPerturbation::tpi_layers
 * @brief Number of layers of perturbation boxes.
 * @var TurbulentPerturbation::tpi_offset
 * @brief Cell offset for the start of the perturbation region.
 * @var TurbulentPerturbation::tpi_boxDim
 * @brief Dimensions of each perturbation box.
 * @var TurbulentPerturbation::tpi_direction
 * @brief Boundary directions where perturbations are applied.
 * @var TurbulentPerturbation::tpi_nonDim
 * @brief Nondimensional number used by the perturbation formulation.
 * @var TurbulentPerturbation::tpi_Ti
 * @brief Temperature intensity used by the perturbation formulation.
 * @var TurbulentPerturbation::tpi_Tinf
 * @brief Reference temperature used by the perturbation formulation [K].
 * @var TurbulentPerturbation::tpi_Hpb
 * @brief Perturbation-box height for each AMR level [m].
 * @var TurbulentPerturbation::tpi_Lpb
 * @brief Perturbation-box length for each AMR level [m].
 * @var TurbulentPerturbation::tpi_Wpb
 * @brief Perturbation-box width for each AMR level [m].
 * @var TurbulentPerturbation::tpi_lref
 * @brief Perturbation-box reference length for each AMR level [m].
 * @var TurbulentPerturbation::tpi_net_buoyant
 * @brief Accumulated net buoyant perturbation used for correction.
 * @var TurbulentPerturbation::tpi_pert_adjust
 * @brief Per-cell perturbation adjustment used for net-zero buoyancy.
 * @var TurbulentPerturbation::ref_ratio
 * @brief Refinement ratios used by multilevel perturbation scaling.
 * @var TurbulentPerturbation::input_Ug
 * @brief Input geostrophic wind speed used to scale CPM perturbations.
 * @var TurbulentPerturbation::input_w_amp
 * @brief Input vertical velocity perturbation amplitude.
 * @var TurbulentPerturbation::perturbation_klo
 * @brief Lower vertical index bound for perturbations.
 * @var TurbulentPerturbation::perturbation_khi
 * @brief Upper vertical index bound for perturbations.
 * @var TurbulentPerturbation::pb_interval
 * @brief Perturbation update interval for each box [s].
 * @var TurbulentPerturbation::pb_local_etime
 * @brief Local elapsed time for each perturbation box [s].
 * @var TurbulentPerturbation::pb_amp
 * @brief Perturbation amplitude for each perturbation box.
 * @var TurbulentPerturbation::pb_netZero
 * @brief Per-box storage used for net-zero buoyancy calculation.
 * @var TurbulentPerturbation::pb_nupdate
 * @brief Number of times each perturbation box has been refreshed; keys its random numbers.
 * @var TurbulentPerturbation::tpi_seed
 * @brief Seed of the hashed random numbers (1024 with erf.fix_random_seed, else erf.random_seed, else 0).
 * @var TurbulentPerturbation::hash_comp_cell
 * @brief Hash salt of the per-cell perturbation values.
 * @var TurbulentPerturbation::hash_comp_cpm
 * @brief Hash salt of the per-box CPM values.
 * @var TurbulentPerturbation::hash_comp_interval
 * @brief Hash salt of the refresh-interval variation.
 */
#endif
