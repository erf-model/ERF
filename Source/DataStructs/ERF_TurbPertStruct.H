#ifndef ERF_TURB_PERT_STRUCT_H_
#define ERF_TURB_PERT_STRUCT_H_

#include <ERF_DataStruct.H>
#include <AMReX_MultiFabUtil.H>
#include <ERF_TileNoZ.H>
#include <time.h>
/**
 * Container holding quantities related to turbulent perturbation parameters
 */

/*  The general rule of thumb is to create a perturbation box size of 1/8th of the boundary layer height.
    The boundary layer height can't be the height of the domain. The length and width of the box should
    be twice the height of the box. If meandering flow is present, the width of the box should be take
    the angle of the inflow into consideration.
*/

AMREX_ENUM(PerturbationType,
    Source, Direct, CPM, None
);

struct TurbulentPerturbation {

  public:

    ~TurbulentPerturbation () {}

    void init_tpi_type (const PerturbationType& pert_type)
    {
        if (pert_type == PerturbationType::Source) {
            pt_type = 0;
        } else if (pert_type == PerturbationType::Direct) {
            pt_type = 1;
        } else if (pert_type == PerturbationType::CPM) {
            pt_type = 2;
        }
        AMREX_ALWAYS_ASSERT(pt_type >= 0);
    }

    // Initializing Perturbation Region
    // Currently only support perturbation in the x and y direction
    // This portion of the code is only called once for initialization,
    // therefore efficiency is not required. Dont spend too much time
    // here other than to correctly setup the perturbation box region
    void init_tpi (const int lev,
                   const amrex::IntVect& valid_box_lo, // start points of the valid box
                   const amrex::IntVect& valid_box_hi, // end points of the valid box
                   const amrex::GpuArray<amrex::Real,3> dx,
                   const amrex::BoxArray& ba,
                   const amrex::DistributionMapping& dm,
                   const int ngrow_state,
                   std::string pp_prefix,
                   const amrex::Vector<amrex::IntVect> refRatio,
                   const int max_level)

    {
        // Initialization for some 0 dependent terms
        tpi_Ti = 0.;
        tpi_Tinf = 300.;

        ref_ratio = refRatio;

        amrex::ParmParse pp(pp_prefix);

        // Reading inputs, and placing assertion for the perturbation inflow to work
        pp.queryarr("perturbation_box_dims",tpi_boxDim);
        pp.queryarr("perturbation_direction",tpi_direction);
        pp.query("perturbation_layers",tpi_layers);
        pp.query("perturbation_offset",tpi_offset);

        pp.query("perturbation_nondimensional",tpi_nonDim);
        pp.query("perturbation_T_infinity",tpi_Tinf);
        pp.query("perturbation_T_intensity",tpi_Ti);
        pp.query("perturbation_Ug",input_Ug);

        // Check variables message
        if (tpi_layers < 0) { amrex::Abort("Please provide a valid perturbation layer value (ie. 3-5)"); }
        if (tpi_offset < 0) { amrex::Abort("Please provide a valid inflow cell offset value for perturbation region (ie. 0-5)"); }
        for (int i = 0; i < tpi_boxDim.size(); i++) {
            if (tpi_boxDim[i] == 0) { amrex::Abort("Please provide valid dimensions for perturbation boxes."); }
        }
        if (tpi_nonDim < 0.) { amrex::Abort("Please provide a valid nondimensional number (ie. Ri = 0.042)"); }
        if (input_Ug < 0.) { amrex::Abort("Please provide a valid geostrophic wind speed (ie. Ug = 10.0 m/s)"); }
        if (tpi_Tinf < 0.) { amrex::Abort("Please provide a valid ambient temperature value (ie. T_0 = T_infty)"); }
        if (tpi_Ti < 0.) { amrex::Abort("Please provide a valid temperature intensity value (ie. 0-1.0)"); }

        // Creating perturbation regions and initializing with generic size. Temporary size for now
        amrex::Box lo_x_bx(amrex::IntVect(0), amrex::IntVect(1), amrex::IntVect(0));
        amrex::Box hi_x_bx(amrex::IntVect(0), amrex::IntVect(1), amrex::IntVect(0));
        amrex::Box lo_y_bx(amrex::IntVect(0), amrex::IntVect(1), amrex::IntVect(0));
        amrex::Box hi_y_bx(amrex::IntVect(0), amrex::IntVect(1), amrex::IntVect(0));

        // Create a temporary box list to accumulate all the perturbation regions after box modification
        amrex::BoxList tmp_bl;

        // boxSize for individual boxes
        amrex::IntVect boxSize(tpi_boxDim[0],tpi_boxDim[1],tpi_boxDim[2]);

        // Starting logic to set the size of the perturbation region(s)
        //amrex::PrintToFile("BoxPerturbationOutput") << "Setting perturbation region in:";
        // ***** X-direction perturbation *****
        if (tpi_direction[0]) { // West
            lo_x_bx.setSmall(amrex::IntVect(valid_box_lo[0]+tpi_offset, valid_box_lo[1]+tpi_direction[1]*tpi_offset, valid_box_lo[2]));
            lo_x_bx.setBig(amrex::IntVect(valid_box_lo[0]+(tpi_layers*tpi_boxDim[0]-1)+tpi_offset, valid_box_hi[1]-(tpi_direction[4]*tpi_offset), valid_box_hi[2]));
            //amrex::PrintToFile("BoxPerturbationOutput") << " West face";
        }

        if (tpi_direction[3]) { // East
            hi_x_bx.setSmall(amrex::IntVect(valid_box_hi[0]-((tpi_layers*tpi_boxDim[0]-1)+tpi_offset), valid_box_lo[1]+tpi_direction[1]*tpi_offset, valid_box_lo[2]));
            hi_x_bx.setBig(amrex::IntVect(valid_box_hi[0]-tpi_offset, valid_box_hi[1]-(tpi_direction[4]*tpi_offset), valid_box_hi[2]));
            //amrex::PrintToFile("BoxPerturbationOutput") << " East face";
        }

        // ***** Y-direction Perturbation *****
        if (tpi_direction[1]) { // North
            lo_y_bx.setSmall(amrex::IntVect(valid_box_lo[0]+tpi_direction[0]*tpi_offset, valid_box_lo[1]+tpi_offset, valid_box_lo[2]));
            lo_y_bx.setBig(amrex::IntVect(valid_box_hi[0]-tpi_direction[3]*tpi_offset, valid_box_lo[1]+((tpi_layers*tpi_boxDim[1])-1)+tpi_offset, valid_box_hi[2]));
            //amrex::PrintToFile("BoxPerturbationOutput") << " North face";
        }

        if (tpi_direction[4]) { // South
            hi_y_bx.setSmall(amrex::IntVect(valid_box_lo[0]+tpi_direction[0]*tpi_offset, valid_box_hi[1]-((tpi_layers*tpi_boxDim[1]-1)+tpi_offset), valid_box_lo[2]));
            hi_y_bx.setBig(amrex::IntVect(valid_box_hi[0]-tpi_direction[3]*tpi_offset, valid_box_hi[1]-tpi_offset, valid_box_hi[2]));
            //amrex::PrintToFile("BoxPerturbationOutput") << " South face";
        }

        if (tpi_direction[2] || tpi_direction[5]) { amrex::Abort("Currently not supporting z-direction flow perturbation"); }

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
        //amrex::PrintToFile("BoxPerturbationOutput") << "\nBoxList: " << tmp_bl << "\n";
        amrex::BoxArray tmp_ba(tmp_bl);
        tmp_ba.maxSize(boxSize);
        pb_ba.push_back(tmp_ba);

        if (lev == 0) {
            pb_mag.resize(max_level+1);
            pb_dir.resize(max_level+1);
            pb_netZero.resize(max_level+1);
            pb_interval.resize(max_level+1);
            pb_local_etime.resize(max_level+1);
            pb_amp.resize(max_level+1);
            pb_cell.resize(max_level+1);
            tpi_Lpb.resize(max_level+1);
            tpi_Wpb.resize(max_level+1);
            tpi_Hpb.resize(max_level+1);
            tpi_lref.resize(max_level+1);
        }

        // Initializing mean magnitude and direction vectors
        pb_mag[lev].resize(pb_ba[lev].size(), 0.);
        pb_dir[lev].resize(pb_ba[lev].size(), 0.);
        pb_netZero[lev].resize(pb_ba[lev].size(), 0.);

        // Set size of vector and initialize
        pb_interval[lev].resize(pb_ba[lev].size(), -1.0);
        pb_local_etime[lev].resize(pb_ba[lev].size(), 0.0);
        pb_amp[lev].resize(pb_ba[lev].size(), 0.0);

        // Creating data array for perturbation amplitude storage
        pb_cell[lev].define(ba, dm, 1, ngrow_state); // this is the only place ba is used. Maybe we can print here to determine what's valid...
        pb_cell[lev].setVal(0.);

        // Computing perturbation reference length
        tpi_Lpb[lev] = tpi_boxDim[0]*dx[0];
        tpi_Wpb[lev] = tpi_boxDim[1]*dx[1];
        tpi_Hpb[lev] = tpi_boxDim[2]*dx[2];
        tpi_lref[lev] = sqrt(tpi_Lpb[lev]*tpi_Lpb[lev] + tpi_Wpb[lev]*tpi_Wpb[lev]);

        tpi_pert_adjust = 0.;
        tpi_net_buoyant = 0.;

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


    // Perturbation update frequency check.
    // This function will trigger calc_tpi_meanMag() and calc_tpi_amp(), when
    // the local elapsed time is greater than the perturbation frequency.
    // At each timestep a total buoyant force sanity check is made to ensure
    // that the sum of the buoyant force introduced into the system is net-zero
    void calc_tpi_update (const int lev,
                          const amrex::Real dt,
                          amrex::MultiFab& mf_xvel,
                          amrex::MultiFab& mf_yvel,
                          amrex::MultiFab& mf_cons)
    {
        // Resettubg the net buoyant force value
        tpi_net_buoyant = 0.;

        // Setting random number generator for update interval
        srand( (unsigned) time(NULL) );

        auto m_ixtype = mf_cons.boxArray().ixType(); // safety step

        // Seed the random generator at 1024UL for regression testing
        int fix_random_seed = 0;
        amrex::ParmParse pp("erf");
        pp.query("fix_random_seed", fix_random_seed);
        if (fix_random_seed) {
            // We need this one for the ParalleForRNG used in calc_tpi
            amrex::InitRandom(1024UL);

            // We need this one for the RandomReal below
            srand(1024UL);
        }

        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {

            bool update_box = true; // initialize flag
            // Check if the local elapsed time is greater than the update interval and don't go into boxes the rank doesn't own
            if (pt_type == 2) {
                // for CPM, perturbation boxes/cells refresh after boxes have advected box width/length * num_layers (advective time scale)
                update_box = ((pb_local_etime[lev][boxIdx] >= pb_interval[lev][boxIdx]*tpi_layers) && pb_local_etime[lev][boxIdx] != -1.0);
            } else {
                update_box = ( pb_local_etime[lev][boxIdx] >= pb_interval[lev][boxIdx] && pb_local_etime[lev][boxIdx] != -1.0 );
            }
            if ( update_box ) {

                // Compute mean velocity of each perturbation box
                calc_tpi_meanMag_perBox(boxIdx, lev, mf_cons, mf_xvel, mf_yvel);

                // Only the rank owning the box will be able to access the storage location
                // Done for parallelism to avoid Inf being stored in array
                if (pb_mag[lev][boxIdx] !=0.) {
                    amrex::Real interval = 0.0;
                    if (pt_type == 2) {
                        //  Wind direction correction for angled wind
                        amrex::Real wind_direction = pb_dir[lev][boxIdx];
                        if (wind_direction > PI / 4) { wind_direction = PI / 2 - wind_direction; }
                        // CPM only cares about the side length of a box, min call maintains flexibility.
                        interval = 1.0 / std::cos(wind_direction) * std::min(tpi_Lpb[lev], tpi_Wpb[lev]) / pb_mag[lev][boxIdx];
                    } else {
                        interval = tpi_lref[lev] / pb_mag[lev][boxIdx];
                    }
                    pb_interval[lev][boxIdx] = RandomReal(0.9*interval,1.1*interval); // 10% variation

                    // Reset local elapsed time
                    pb_local_etime[lev][boxIdx] = 0.;
                } else {
                    // this box is not on this rank, we shouldn't enter it again.
                    // Technically, all boxes are looped through on the very first step of a simulation and this is when this is set
                    pb_local_etime[lev][boxIdx] = -1.0;
                }

                // Trigger amplitude calculation per perturbation box
                calc_tpi_amp(lev, boxIdx, pb_interval[lev][boxIdx]);

                // Trigger random amplitude storage per cell within perturbation box
                pseudoRandomPert(boxIdx, lev, m_ixtype);

            } else {
                // set perturbation amplitudes to 0 for CPM. A little inefficient but leverages as much as existing code as possible.
                if (pt_type == 2) { zero_amp(boxIdx, lev, m_ixtype); }

                // Increase by timestep of level 0 (but only if the box is owned by this rank)
                if (pb_local_etime[lev][boxIdx] != -1.0) { pb_local_etime[lev][boxIdx] += dt; }
            } // if

            if (pt_type < 2) { // box perturbation method only
                // Per iteration operation of net-zero buoyant force check
                if (pb_mag[lev][boxIdx] !=0.) netZeroBuoyantAdd(boxIdx, lev);
                tpi_net_buoyant += pb_netZero[lev][boxIdx];
            }
        } // for

        if (pt_type < 2) { // box perturbation method only
            // Normalizing the adjustment based on how many boxes there are
            // the values within the array is already normalized by the number
            // of cells within each box
            tpi_pert_adjust = tpi_net_buoyant / (amrex::Real) pb_ba[lev].size();

            // Per iteration operation of net-zero buoyant force adjustment
            for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
                if (pb_mag[lev][boxIdx] !=0.) netZeroBuoyantAdjust(boxIdx, lev);
            }
        }
    }

    // Applying perturbation amplitude onto source term (Umphrey and Senocak 2016)
    // This function does per cell application within the box union. Random perturbation
    // is assigned in calc_tpi_update.
    void apply_tpi (const int& lev,
                    const amrex::Box& vbx,                      // box union from upper level
                    const int& comp,                            // Component to modify
                    const amrex::IndexType& m_ixtype,           // IntVect type of src_arr
                    const amrex::Array4<amrex::Real>& src_arr,  // Array to apply perturbation
                    const amrex::Array4<amrex::Real const>& pert_cell)
    {
        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            amrex::Box pbx = amrex::convert(pb_ba[lev][boxIdx], m_ixtype);
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    src_arr(i,j,k,comp) += pert_cell(i,j,k);

                    // For box region debug only
                    #ifdef INDEX_PERTURB
                    src_arr(i,j,k,comp) = (amrex::Real) (boxIdx + 5.);
                    #endif
                });
            }
        }
    }

    // Perturbation amplitude calculation
    void calc_tpi_amp (const int& lev,
                       const int& boxIdx,
                       const amrex::Real& interval)
    {
        pb_amp[lev][boxIdx] = 0.; // Safety step
        if (pt_type == 2) { // CPM
            amrex::Real cp = 1004; // specific heat of air [J/(kg K)]
            amrex::Real Ec = 0.2; // Eckert number
            pb_amp[lev][boxIdx] = (input_Ug * input_Ug) / (Ec * cp);
        } else { // box perturbation
            amrex::Real Um = pb_mag[lev][boxIdx];
            amrex::Real beta = 1./tpi_Tinf; // Thermal expansion coefficient

            // Pseudo Random temperature the ignores scale when mechanically tripping turbulence
            amrex::Real g = CONST_GRAV;
            // get total refinement ratio on each level
            int total_ref_ratio = 1;
            for (int level = lev; level >= 1; level--) {
                total_ref_ratio *= ref_ratio[level-1][2];
            }
            // calculation needs to be scale-aware since the formulation relies on the physical size of the box
            if (tpi_Ti > 0.) g = (tpi_nonDim * Um * Um) / (tpi_Ti * tpi_Hpb[lev]) * 1 / total_ref_ratio;

            // Ma and Senocak (2023) Eq. 8, solving for delta phi
            pb_amp[lev][boxIdx] = (tpi_nonDim * Um * Um) / (g * beta * tpi_Hpb[lev]) * 1 / total_ref_ratio;

            if (pt_type == 0) {
                // Performing this step converts the perturbation proportionality into
                // the forcing term
                // Ma & Senocak (2023) Eq. 7
                pb_amp[lev][boxIdx] /= interval;
            }
        }
    }

    // Assigns pseudo-random (ie. white noise) perturbation to a storage cell, this
    // value is then held constant for the duration of the update interval and assigned onto
    // the source term
    void pseudoRandomPert (const int& boxIdx,
                           const int& lev,
                           const amrex::IndexType& m_ixtype)
    {
        for (amrex::MFIter mfi(pb_cell[lev],TileNoZ()); mfi.isValid(); ++mfi) {
            amrex::Box vbx = mfi.validbox();
            amrex::Box pbx = amrex::convert(pb_ba[lev][boxIdx], m_ixtype);
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                amrex::Real amp_copy = pb_amp[lev][boxIdx];
                const amrex::Array4<amrex::Real>& pert_cell = pb_cell[lev].array(mfi);

                if (pt_type == 2) { // CPM
                    amrex::Real rand_number_const = RandomReal(-1.0, 1.0);
                    ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        pert_cell(i,j,k) = rand_number_const * amp_copy;
                    });
                } else {
                    ParallelForRNG(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                        amrex::Real rand_double = amrex::Random(engine);
                        pert_cell(i,j,k) = (rand_double*2.0 - 1.0) * amp_copy;
                    });
                }

            }
        }
    }

    // convert pert_cell value to be 0 so that perturbation amplitudes do not compound for CPM
    void zero_amp (const int& boxIdx,
                   const int& lev,
                   const amrex::IndexType& m_ixtype)
    {

        for (amrex::MFIter mfi(pb_cell[lev],TileNoZ()); mfi.isValid(); ++mfi) {
            amrex::Box vbx = mfi.validbox();
            amrex::Box pbx = amrex::convert(pb_ba[lev][boxIdx], m_ixtype);
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                const amrex::Array4<amrex::Real>& pert_cell = pb_cell[lev].array(mfi);
                ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    pert_cell(i,j,k) = 0.0;
                });
            }
        }
    }

    // Checks for net-zero buoyant force introduction into the system
    void netZeroBuoyantAdd (const int& boxIdx,
                            const int& lev)
    {
        // Creating local copy of PB box array and magnitude
        const amrex::BoxArray m_pb_ba = pb_ba[lev];
        amrex::Real* m_pb_netZero = pb_netZero[lev].data();

        // Create device array for summation
        amrex::Vector<amrex::Real> avg_h(1,0.);
        amrex::Gpu::DeviceVector<amrex::Real> avg_d(1,0.);
        amrex::Real* avg = avg_d.data();

        // Iterates through the cells of each box and sum the white noise perturbation
        for (amrex::MFIter mfi(pb_cell[lev], TileNoZ()) ; mfi.isValid(); ++mfi) {
            const amrex::Box& vbx = mfi.validbox();
            amrex::Box pbx = amrex::convert(m_pb_ba[boxIdx], vbx.ixType());
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                const amrex::Array4<const amrex::Real>& pert_cell = pb_cell[lev].const_array(mfi);
                amrex::Real norm = 1.0 / static_cast<amrex::Real>(ubx.numPts());
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&avg[0], pert_cell(i,j,k)*norm, handler);
                });
                amrex::Gpu::copy(amrex::Gpu::deviceToHost, avg_d.begin(), avg_d.end(), avg_h.begin());

                // Assigning onto storage array
                m_pb_netZero[boxIdx] = avg_h[0];
            }
        }
    }

    // If it's not a net-zero buoyant force, then adjust all cell by a normalized value
    // to achieve this logic
    void netZeroBuoyantAdjust (const int& boxIdx,
                               const int& lev)
    {
        // Creating local copy of PB box array and magnitude
        const amrex::BoxArray m_pb_ba = pb_ba[lev];
        for (amrex::MFIter mfi(pb_cell[lev], TileNoZ()) ; mfi.isValid(); ++mfi) {
            const amrex::Box& vbx = mfi.validbox();
            amrex::Box pbx = amrex::convert(m_pb_ba[boxIdx], vbx.ixType());
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
               const amrex::Real adjust = tpi_pert_adjust;
               const amrex::Array4<amrex::Real>& pert_cell = pb_cell[lev].array(mfi);
               ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   pert_cell(i,j,k) -= adjust;
               });
            }
        }
    }

// TODO: Test the difference between these two for Source term perturbation
#define USE_VOLUME_AVERAGE
    // Perturbation box mean velocity magnitude calculation
    // This is pulled into the structure to also utilize during runtime
    void calc_tpi_meanMag_perBox (const int& boxIdx,
                                  const int& lev,
                                  amrex::MultiFab& mf_cons,
                                  amrex::MultiFab& mf_xvel,
                                  amrex::MultiFab& mf_yvel)

    {
        // Creating local copy of PB box array and magnitude
        const amrex::BoxArray m_pb_ba = pb_ba[lev];
        amrex::Real* m_pb_mag = pb_mag[lev].data();
        amrex::Real* m_pb_dir = pb_dir[lev].data();
        m_pb_mag[boxIdx] = 0.; // Safety step
        m_pb_dir[boxIdx] = 0.; // Safety step

        // Storage of averages per PB
        // Index: 0=u (vol/slab_lo), 1=v (vol/slab_lo)
        //        2=u (slab_hi),     3=v (slab_hi)
        int n_avg = 4;
        amrex::Vector<amrex::Real> avg_h(n_avg,0.);
        amrex::Gpu::DeviceVector<amrex::Real> avg_d(n_avg,0.);
        amrex::Real* avg = avg_d.data();

        // Averaging u & v components in single MFIter
        for (amrex::MFIter mfi(mf_cons, TileNoZ()); mfi.isValid(); ++mfi) {

            // CC valid box (inherited from mf_cons)
            const amrex::Box& vbx = mfi.validbox();

            // Box logic for u velocity
            auto ixtype_u    = mf_xvel.boxArray().ixType();
            amrex::Box vbx_u = amrex::convert(vbx,ixtype_u);
            amrex::Box pbx_u = amrex::convert(m_pb_ba[boxIdx], ixtype_u);
            amrex::Box ubx_u = pbx_u & vbx_u;

            // Box logic for v velocity
            auto ixtype_v    = mf_yvel.boxArray().ixType();
            amrex::Box vbx_v = amrex::convert(vbx,ixtype_v);
            amrex::Box pbx_v = amrex::convert(m_pb_ba[boxIdx], ixtype_v);
            amrex::Box ubx_v = pbx_v & vbx_v;

            // Operation over box union (U)
            if (ubx_u.ok()) {
                const amrex::Array4<const amrex::Real>& xvel_arry = mf_xvel.const_array(mfi);

                #ifdef USE_VOLUME_AVERAGE
                amrex::Real norm = 1.0 / static_cast<amrex::Real>(ubx_u.numPts());
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubx_u, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&avg[0], xvel_arry(i,j,k)*norm, handler);
                });
                #endif // USE_VOLUME_AVERAGE

                #ifdef USE_SLAB_AVERAGE
                amrex::Box ubxSlab_lo = makeSlab(ubx_u,2,ubx_u.smallEnd(2));
                amrex::Box ubxSlab_hi = makeSlab(ubx_u,2,ubx_u.bigEnd(2));
                amrex::Real norm_lo = 1.0 / static_cast<amrex::Real>(ubxSlab_lo.numPts());
                amrex::Real norm_hi = 1.0 / static_cast<amrex::Real>(ubxSlab_hi.numPts());

                // Average u in the low slab
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubxSlab_lo, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&avg[0], xvel_arry(i,j,k)*norm_lo, handler);
                });

                // Average u in the high slab
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubxSlab_hi, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&avg[2], xvel_arry(i,j,k)*norm_hi, handler);
                });
                #endif // USE_SLAB_AVERAGE
            } // if

            // Operation over box union (V)
            if (ubx_v.ok()) {
                const amrex::Array4<const amrex::Real>& yvel_arry = mf_yvel.const_array(mfi);

                #ifdef USE_VOLUME_AVERAGE
                amrex::Real norm = 1.0 / static_cast<amrex::Real>(ubx_v.numPts());
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubx_v, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&avg[1], yvel_arry(i,j,k)*norm, handler);
                });
                #endif // USE_VOLUME_AVERAGE

                #ifdef USE_SLAB_AVERAGE
                amrex::Box ubxSlab_lo = makeSlab(ubx_v,2,ubx_v.smallEnd(2));
                amrex::Box ubxSlab_hi = makeSlab(ubx_v,2,ubx_v.bigEnd(2));
                amrex::Real norm_lo = 1.0 / static_cast<amrex::Real>(ubxSlab_lo.numPts());
                amrex::Real norm_hi = 1.0 / static_cast<amrex::Real>(ubxSlab_hi.numPts());

                // Average v in the low slab
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubxSlab_lo, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&avg[1], yvel_arry(i,j,k)*norm_lo, handler);
                });

                // Average v in the high slab
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubxSlab_hi, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&avg[3], yvel_arry(i,j,k)*norm_hi, handler);
                });
                #endif // USE_SLAB_AVERAGE
            } // if
        } // MFIter

        // Copy from device back to host
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, avg_d.begin(), avg_d.end(), avg_h.begin());

        // Computing the average magnitude within PB
        #ifdef USE_VOLUME_AVERAGE
        m_pb_mag[boxIdx] = sqrt(avg_h[0]*avg_h[0] + avg_h[1]*avg_h[1]);
        m_pb_dir[boxIdx] = std::atan(std::abs(avg_h[0]) / std::abs(avg_h[1]+std::numeric_limits<amrex::Real>::epsilon()));
        #endif

        #ifdef USE_SLAB_AVERAGE
        m_pb_mag[boxIdx] = 0.5*( sqrt(avg_h[0]*avg_h[0] + avg_h[1]*avg_h[1])
                               + sqrt(avg_h[2]*avg_h[2] + avg_h[3]*avg_h[3]));
        m_pb_dir[boxIdx] = std::atan(std::abs(0.5*(avg_h[0]+avg_h[2])) / std::abs(0.5*(avg_h[1]+avg_h[3])+std::numeric_limits<amrex::Real>::epsilon()));
        #endif
    }

    // Output debug message into a file
    void debug (amrex::Real /*time*/)
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

    int pt_type = -1;

    // Public data members
    amrex::Vector<amrex::BoxArray>                 pb_ba;  // PB box array
    amrex::Vector<amrex::Vector<amrex::Real>>      pb_mag; // BP mean magnitude [m/s]
    amrex::Vector<amrex::Vector<amrex::Real>>      pb_dir; // BP mean direction [deg]

    // Perturbation amplitude cell storage
    // This is after random assignment of equation (10) in Ma and Senocak 2023
    amrex::Vector<amrex::MultiFab>  pb_cell;

  private:

    // Private data members
    int tpi_layers;                   // Number of layers of perturbation boxes
    int tpi_offset;                   // Cells to offset the start of the perturbation region

    amrex::Vector<int> tpi_boxDim;    // Dimensions of each perturbation box
    amrex::Vector<int> tpi_direction; // Direction of the peturbation

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

    // Perturbation data arrays
    amrex::Vector<amrex::Vector<amrex::Real>> pb_interval;    // PB update time [s]
    amrex::Vector<amrex::Vector<amrex::Real>> pb_local_etime; // PB local elapsed time [s]
    amrex::Vector<amrex::Vector<amrex::Real>> pb_amp;         // PB perturbation amplitude Ri:[K]
    amrex::Vector<amrex::Vector<amrex::Real>> pb_netZero;     // PB array used for net zero sum calculation

    // Random number generation between range (used for interval calculation)
    amrex::Real RandomReal (const amrex::Real min, const amrex::Real max)
    {
        amrex::Real r = (amrex::Real) rand() / (amrex::Real) RAND_MAX;
        return min + r * (max - min);
    }

};
#endif