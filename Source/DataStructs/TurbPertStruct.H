#ifndef _TURB_PERT_STRUCT_H_
#define _TURB_PERT_STRUCT_H_

#include <DataStruct.H>
#include <AMReX_MultiFabUtil.H>
#include <TileNoZ.H>
#include <time.h>
/**
 * Container holding quantities related to turbulent perturbation parameters
 */

/*  The general rule of thumb is to create a perturbation box size of 1/8th of the boundary layer height.
    The boundary layer height can't be the height of the domain. The length and width of the box should
    be twice the height of the box. If meandering flow is present, the width of the box should be take
    the angle of the inflow into consideration.
*/

struct TurbulentPerturbation {

  public:

    ~TurbulentPerturbation () {}

    // Initializing Perturbation Region
    // Currently only support perturbation in the x and y direction
    // This portion of the code is only called once for initialization,
    // therefore efficiency is not required. Dont spend too much time
    // here other than to correctly setup the perturbation box region
    void init_tpi (const int lev,
                   const amrex::IntVect& nx,
                   const amrex::GpuArray<amrex::Real,3> dx,
                   const amrex::BoxArray& ba,
                   const amrex::DistributionMapping& dm,
                   const int ngrow_state)

    {
        // Initialization for some 0 dependent terms
        tpi_Ti = 0.;
        tpi_Tinf = 300.;

        amrex::ParmParse pp(pp_prefix);

        // Reading inputs, and placing assertion for the perturbation inflow to work
        pp.queryarr("perturbation_box_dims",tpi_boxDim);
        pp.queryarr("perturbation_direction",tpi_direction);
        pp.query("perturbation_layers",tpi_layers);
        pp.query("perturbation_offset",tpi_offset);

        pp.query("perturbation_nondimensional",tpi_nonDim);
        pp.query("perturbation_T_infinity",tpi_Tinf);
        pp.query("perturbation_T_intensity",tpi_Ti);

        // Check variables message
        if (tpi_layers < 0) { amrex::Abort("Please provide a valid perturbation layer value (ie. 3-5)"); }
        if (tpi_offset < 0) { amrex::Abort("Please provide a valid inflow cell offset value for perturbation region (ie. 0-5)"); }
        for (int i = 0; i < tpi_boxDim.size(); i++) {
            if (tpi_boxDim[i] == 0) { amrex::Abort("Please provide valid dimensions for perturbation boxes."); }
        }
        if (tpi_nonDim < 0.) { amrex::Abort("Please provide a valid nondimensional number (ie. Ri = 0.042)"); }
        if (tpi_Tinf < 0.) { amrex::Abort("Please provide a valid ambient temperature value (ie. T_0 = T_infty)"); }
        if (tpi_Ti < 0.) { amrex::Abort("Please provide a valid temperature intensity value (ie. 0-1.0)"); }

        // Creating perturbation regions and initializing with generic size. Temporary size for now
        amrex::Box lo_x_bx(amrex::IntVect(0,0,0), amrex::IntVect(1,1,1), amrex::IntVect(0,0,0));
        amrex::Box hi_x_bx(amrex::IntVect(0,0,0), amrex::IntVect(1,1,1), amrex::IntVect(0,0,0));
        amrex::Box lo_y_bx(amrex::IntVect(0,0,0), amrex::IntVect(1,1,1), amrex::IntVect(0,0,0));
        amrex::Box hi_y_bx(amrex::IntVect(0,0,0), amrex::IntVect(1,1,1), amrex::IntVect(0,0,0));

        // Create a temporary box list to accumulate all the perturbation regions after box modification
        amrex::BoxList tmp_bl;

        // boxSize for individual boxes
        amrex::IntVect boxSize(tpi_boxDim[0],tpi_boxDim[1],tpi_boxDim[2]);

        // Starting logic to set the size of the perturbation region(s)
        amrex::PrintToFile("BoxPerturbationOutput") << "Setting perturbation region in:";
        // ***** X-direction perturbation *****
        if (tpi_direction[0]) { // West
            lo_x_bx.setSmall(amrex::IntVect(tpi_offset, tpi_direction[1]*tpi_offset, 0));
            lo_x_bx.setBig(amrex::IntVect((tpi_layers*tpi_boxDim[0]-1)+tpi_offset, nx[1]-(tpi_direction[4]*tpi_offset), nx[2]));
            amrex::PrintToFile("BoxPerturbationOutput") << " West face";
        }

        if (tpi_direction[3]) { // East
            hi_x_bx.setSmall(amrex::IntVect(nx[0]-((tpi_layers*tpi_boxDim[0]-1)+tpi_offset), tpi_direction[1]*tpi_offset, 0));
            hi_x_bx.setBig(amrex::IntVect(nx[0]-tpi_offset, nx[1]-(tpi_direction[4]*tpi_offset), nx[2]));
            amrex::PrintToFile("BoxPerturbationOutput") << " East face";
        }

        // ***** Y-direction Perturbation *****
        if (tpi_direction[1]) { // North
            lo_y_bx.setSmall(amrex::IntVect(tpi_direction[0]*tpi_offset, tpi_offset, 0));
            lo_y_bx.setBig(amrex::IntVect(nx[0]-tpi_direction[3]*tpi_offset, ((tpi_layers*tpi_boxDim[1])-1)+tpi_offset, nx[2]));
            amrex::PrintToFile("BoxPerturbationOutput") << " North face";
        }

        if (tpi_direction[4]) { // South
            hi_y_bx.setSmall(amrex::IntVect(tpi_direction[0]*tpi_offset, nx[1]-((tpi_layers*tpi_boxDim[1]-1)+tpi_offset), 0));
            hi_y_bx.setBig(amrex::IntVect(nx[0]-tpi_direction[3]*tpi_offset, nx[1]-tpi_offset, nx[2]));
            amrex::PrintToFile("BoxPerturbationOutput") << " South face";
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
        amrex::PrintToFile("BoxPerturbationOutput") << "\nBoxList: " << tmp_bl << "\n";
        amrex::BoxArray tmp_ba(tmp_bl);
        tmp_ba.maxSize(boxSize);
        pb_ba.push_back(tmp_ba);

        // Initializing mean magnitude vector
        pb_mag.resize(pb_ba[lev].size(), 0.);
        pb_netZero.resize(pb_ba[lev].size(), 0.);

        // Set size of vector and initialize
        pb_interval.resize(pb_ba[lev].size(), -1.0);
        pb_local_etime.resize(pb_ba[lev].size(), 0.0);
        pb_amp.resize(pb_ba[lev].size(), 0.0);

        // Creating data array for perturbation amplitude storage
        pb_cell.define(ba, dm, 1, ngrow_state);
        pb_cell.setVal(0.);

        // Computing perturbation reference length
        tpi_Lpb = tpi_boxDim[0]*dx[0];
        tpi_Wpb = tpi_boxDim[1]*dx[1];
        tpi_Hpb = tpi_boxDim[2]*dx[2];
        tpi_lref = sqrt(tpi_Lpb*tpi_Lpb + tpi_Wpb*tpi_Wpb);

        tpi_pert_adjust = 0.;
        tpi_net_buoyant = 0.;

        // Function check point message
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_box_dims: " << tpi_boxDim[0] << " "
                                                                                << tpi_boxDim[1] << " "
                                                                                << tpi_boxDim[2] << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_direction: " << tpi_direction[0] << " "
                                                                                 << tpi_direction[1] << " "
                                                                                 << tpi_direction[2] << "\n\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_layers: " << tpi_layers << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_offset: " << tpi_offset << "\n\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_nondimensional: " << tpi_nonDim << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_T_infinity: " << tpi_Tinf << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_T_intensity: " << tpi_Ti << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "Reference length per box = " << tpi_lref << "\n\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "Turbulent perturbation BoxArray:\n" << pb_ba[lev] << "\n";
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

        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {

            // Check if the local elapsed time is greater than the update interval
            if ( pb_interval[boxIdx] <= pb_local_etime[boxIdx] ) {

                // Compute mean velocity of each perturbation box
                calc_tpi_meanMag_perBox(boxIdx, lev, mf_xvel, mf_yvel);

                // Only the rank owning the box will be able to access the storage location
                // Done for parallelism to avoid Inf being stored in array
                if (pb_mag[boxIdx] !=0.) {
                    amrex::Real interval = tpi_lref / pb_mag[boxIdx];
                    pb_interval[boxIdx] = RandomReal(0.9*interval,1.1*interval); // 10% variation
                }

                // Trigger amplitude calculation per perturbation box
                calc_tpi_amp(boxIdx, pb_interval[boxIdx]);

                // Trigger random amplitude storage per cell within perturbation box
                pseudoRandomPert(boxIdx, lev, m_ixtype);

                // Reset local elapsed time
                pb_local_etime[boxIdx] = 0.;
            } else {
                // Increase by timestep of level 0
                pb_local_etime[boxIdx] += dt;
            } // if

            // Per iteration operation of net-zero buoyant force check
            if (pb_mag[boxIdx] !=0.) netZeroBuoyantAdd(boxIdx, lev);
            tpi_net_buoyant += pb_netZero[boxIdx];
        } // for

        // Normalizing the adjustment based on how many boxes there are
        // the values within the array is already normalized by the number
        // of cells within each box
        tpi_pert_adjust = tpi_net_buoyant / (amrex::Real) pb_ba[lev].size();

        // Per iteration operation of net-zero buoyant force adjustment
        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            if (pb_mag[boxIdx] !=0.) netZeroBuoyantAdjust(boxIdx, lev);
        }
    }

//#define INDEX_PERTURB

    // Applying perturbation amplitude onto source term (Umphrey and Senocak 2016)
    // This function does per cell application within the box union. Random perturbation
    // is assigned in calc_tpi_update.
    void apply_tpi (const int lev,
                    const amrex::Box& vbx,                     // box union from upper level
                    const int comp,                            // Component to modify
                    const amrex::IndexType m_ixtype,           // IntVect type of src_arr
                    const amrex::Array4<amrex::Real> &src_arr, // Array to apply perturbation
                    const amrex::Array4<amrex::Real const> pert_cell)
    {
        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            amrex::Box pbx = convert(pb_ba[lev][boxIdx], m_ixtype);
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    src_arr(i,j,k,comp) += pert_cell(i,j,k,0);

                    // For box region debug only
                    #ifdef INDEX_PERTURB
                    src_arr(i,j,k,comp) = (amrex::Real) (boxIdx + 5.);
                    #endif
                });
            }
        }
    }

    // Quality of life function definitions
    amrex::Real* get_pb_mag() { return pb_mag.data(); }
    amrex::Real* get_pb_netZero() { return pb_netZero.data(); }

    // Output debug message into a file
    void debug (amrex::Real time)
    {
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
    }

    std::string pp_prefix {"erf"};

    int pt_type; // TODO: May need to pass in solverChoice to replace this

    // Public data members
    amrex::Vector<amrex::BoxArray> pb_ba;  // PB box array
    amrex::Vector<amrex::Real>     pb_mag; // BP mean magnitude [m/s]

    // Perturbation amplitude cell storage
    // This is after random assignment of equation (10) in Ma and Senocak 2023
    amrex::MultiFab pb_cell;

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
    amrex::Real tpi_Hpb;              // PB height [m]
    amrex::Real tpi_Lpb;              // PB length [m]
    amrex::Real tpi_Wpb;              // PB width  [m]
    amrex::Real tpi_lref;             // PB reference length [m]

    amrex::Real tpi_net_buoyant;      // Perturbation net-zero calculation storage
    amrex::Real tpi_pert_adjust;      // Perturbation adjust for net-zero per cell adjustment

    // Perturbation data arrays
    amrex::Vector<amrex::Real> pb_interval;    // PB update time [s]
    amrex::Vector<amrex::Real> pb_local_etime; // PB local elapsed time [s]
    amrex::Vector<amrex::Real> pb_amp;         // PB perturbation amplitude Ri:[K]
    amrex::Vector<amrex::Real> pb_netZero;     // PB array used for net zero sum calculation

    // Perturbation amplitude calculation
    void calc_tpi_amp (const int boxIdx, const amrex::Real interval)
    {
        amrex::Real Um = pb_mag[boxIdx];
        pb_amp[boxIdx] = 0.; // Safety step

        amrex::Real beta = 1./tpi_Tinf; // Thermal expansion coefficient

        // Pseudo Random temperature the ignores scale when mechanically tripping turbulence
        amrex::Real g = CONST_GRAV;
        if (tpi_Ti > 0.) g = (tpi_nonDim * Um * Um) / (tpi_Ti * tpi_Hpb);

        // Ma and Senocak (2023) Eq. 8, solving for delta phi
        pb_amp[boxIdx] = (tpi_nonDim * Um * Um) / (g * beta * tpi_Hpb);

        if (pt_type == 0) {
            // Performing this step converts the perturbation proportionality into
            // the forcing term
            // Ma & Senocak (2023) Eq. 7
            pb_amp[boxIdx] /= interval;
        }
    }

    // Assigns pseudo-random (ie. white noise) perturbation to a storage cell, this
    // value is then held constant for the duration of the update interval and assigned onto
    // the source term
    void pseudoRandomPert (const int boxIdx, const int lev, const amrex::IndexType m_ixtype)
    {
        for (amrex::MFIter mfi(pb_cell,TileNoZ()); mfi.isValid(); ++mfi) {
            amrex::Box vbx = mfi.validbox();
            amrex::Box pbx = convert(pb_ba[lev][boxIdx], m_ixtype);
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                amrex::Real amp_copy = pb_amp[boxIdx];
                const amrex::Array4<amrex::Real> &pert_cell = pb_cell.array(mfi);
                ParallelForRNG(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                    amrex::Real rand_double = amrex::Random(engine);
                    pert_cell(i,j,k,0) = (rand_double*2.0 - 1.0) * amp_copy;
                });
            }
        }
    }

    // Checks for net-zero buoyant force introduction into the system
    void netZeroBuoyantAdd (const int boxIdx, const int lev)
    {
        // Creating local copy of PB box array and magnitude
        const amrex::BoxArray m_pb_ba = pb_ba[lev];
        amrex::Real* m_pb_netZero = get_pb_netZero();

        // Create device array for summation
        amrex::Gpu::DeviceVector<amrex::Real> tmp_d(1,0.);
        amrex::Real* tmp = tmp_d.data();

        // Iterates through the cells of each box and sum the white noise perturbation
        for (amrex::MFIter mfi(pb_cell, TileNoZ()) ; mfi.isValid(); ++mfi) {
            const amrex::Box &vbx = mfi.validbox();
            amrex::Box pbx = convert(m_pb_ba[boxIdx], vbx.ixType());
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                const amrex::Array4<const amrex::Real> &pert_cell = pb_cell.array(mfi);
                int npts = ubx.numPts(); // Used in normalization
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[0], pert_cell(i,j,k), handler);
                });
                tmp[0] /= (amrex::Real) npts;

                // Assigning onto storage array
                m_pb_netZero[boxIdx] = tmp[0];
            }
        }
    }

    // If it's not a net-zero buoyant force, then adjust all cell by a normalized value
    // to achieve this logic
    void netZeroBuoyantAdjust (const int boxIdx, const int lev)
    {
        // Creating local copy of PB box array and magnitude
        const amrex::BoxArray m_pb_ba = pb_ba[lev];
        for (amrex::MFIter mfi(pb_cell, TileNoZ()) ; mfi.isValid(); ++mfi) {
            const amrex::Box &vbx = mfi.validbox();
            amrex::Box pbx = convert(m_pb_ba[boxIdx], vbx.ixType());
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
               const amrex::Real adjust = tpi_pert_adjust;
               const amrex::Array4<amrex::Real> &pert_cell = pb_cell.array(mfi);
               ParallelFor(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   pert_cell(i,j,k,0) -= adjust;
               });
            }
        }
    }

// TODO: Test the difference between these two for Source term perturbation
//#define USE_SLAB_AVERAGE
#define USE_VOLUME_AVERAGE

    // Perturbation box mean velocity magnitude calculation
    // This is pulled into the structure to also utilize during runtime
    void calc_tpi_meanMag_perBox (const int boxIdx,
                                  const int lev,
                                  amrex::MultiFab& mf_xvel,
                                  amrex::MultiFab& mf_yvel)

    {
        // Creating local copy of PB box array and magnitude
        const amrex::BoxArray m_pb_ba = pb_ba[lev];
        amrex::Real* m_pb_mag = get_pb_mag();

        // Storage of averages per PB
        // Index: 0=u (vol/slab_lo), 1=v (vol/slab_lo)
        //        2=u (slab_hi),     3=v (slab_hi)
        amrex::Gpu::DeviceVector<amrex::Real> tmp_d(4,0.);
        amrex::Real* tmp = tmp_d.data();

        // Averaging u components
        for (amrex::MFIter mfi(mf_xvel, TileNoZ()) ; mfi.isValid(); ++mfi) {
            const amrex::Box &vbx = mfi.validbox();
            amrex::Box pbx = convert(m_pb_ba[boxIdx], vbx.ixType());
            amrex::Box ubx = pbx & vbx;

            // Operation over box union
            if (ubx.ok()) {
                const amrex::Array4<const amrex::Real> &xvel_arry = mf_xvel.array(mfi);

                #ifdef USE_VOLUME_AVERAGE
                int npts = ubx.numPts();
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[0], xvel_arry(i,j,k), handler);
                });
                tmp[0] /= (amrex::Real) npts;
                #endif // USE_VOLUME_AVERAGE

                #ifdef USE_SLAB_AVERAGE
                amrex::Box ubxSlab_lo = makeSlab(ubx,2,ubx.smallEnd(2));
                amrex::Box ubxSlab_hi = makeSlab(ubx,2,ubx.bigEnd(2));
                int npts_lo = ubxSlab_lo.numPts();
                int npts_hi = ubxSlab_hi.numPts();

                // Average u in the low slab
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubxSlab_lo, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[0], xvel_arry(i,j,k), handler);
                });
                tmp[0] /= (amrex::Real) npts_lo;

                // Average u in the high slab
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubxSlab_hi, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[2], xvel_arry(i,j,k), handler);
                });
                tmp[2] /= (amrex::Real) npts_hi;
                #endif // USE_SLAB_AVERAGE
            } // if
        } // MFIter

        // Averaging v components
        for (amrex::MFIter mfi(mf_yvel, TileNoZ()); mfi.isValid(); ++mfi) {
            amrex::Box const& vbx = mfi.validbox();
            amrex::Box pbx = convert(m_pb_ba[boxIdx], vbx.ixType());
            amrex::Box ubx = pbx & vbx;

            // Operation over box union
            if (ubx.ok()) {
                const amrex::Array4<const amrex::Real> &yvel_arry = mf_yvel.array(mfi);

                #ifdef USE_VOLUME_AVERAGE
                int npts = ubx.numPts();
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[1], yvel_arry(i,j,k), handler);
                });
                tmp[1] /= (amrex::Real) npts;
                #endif // USE_VOLUME_AVERAGE

                #ifdef USE_SLAB_AVERAGE
                amrex::Box ubxSlab_lo = makeSlab(ubx,2,ubx.smallEnd(2));
                amrex::Box ubxSlab_hi = makeSlab(ubx,2,ubx.bigEnd(2));
                int npts_lo = ubxSlab_lo.numPts();
                int npts_hi = ubxSlab_hi.numPts();

                // Average v in the low slab
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubxSlab_lo, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[1], yvel_arry(i,j,k), handler);
                });
                tmp[1] /= (amrex::Real) npts_lo;

                // Average v in the high slab
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubxSlab_hi, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[3], yvel_arry(i,j,k), handler);
                });
                tmp[3] /= (amrex::Real) npts_hi;
                #endif // USE_SLAB_AVERAGE
            } // if
        } // MFIter

        // Computing the average magnitude within PB
        #ifdef USE_SLAB_AVERAGE
        m_pb_mag[boxIdx] = 0.5*(sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]) + sqrt(tmp[2]*tmp[2] + tmp[3]*tmp[3]));
        #endif

        #ifdef USE_VOLUME_AVERAGE
        m_pb_mag[boxIdx] = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
        #endif
    }

    // Random number generation between range (used for interval calculation)
    amrex::Real RandomReal (const amrex::Real min, const amrex::Real max)
    {
        amrex::Real r = (amrex::Real) rand() / (amrex::Real) RAND_MAX;
        return min + r * (max - min);
    }

};
#endif
