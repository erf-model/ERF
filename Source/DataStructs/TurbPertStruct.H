#ifndef _TURB_PERT_STRUCT_H_
#define _TURB_PERT_STRUCT_H_

#include <AMReX_MultiFabUtil.H>
#include <TileNoZ.H>
#include <time.h>

/**
 * Container holding quantities related to turbulent perturbation parameters
 */

struct TurbulentPerturbation {

  public:

    ~TurbulentPerturbation () {}

    // Initializing Perturbation Region
    // Currently only support perturbation in the x and y direction
    // Note: Box will always be initialized at cell centered
    void init_tpi (const int lev,
                   const amrex::IntVect& nx,
                   amrex::GpuArray<amrex::Real,3> const dx)

    {
        amrex::ParmParse pp(pp_prefix);

        // Reading inputs, and placing assertion for the perturbation inflow to work
        pp.queryarr("perturbation_box_dims",tpi_boxDim);
        pp.queryarr("perturbation_direction",tpi_direction);
        pp.query("perturbation_layers",tpi_layers);
        pp.query("perturbation_offset",tpi_offset);

        pp.query("perturbation_nondimensional",tpi_nonDim);
        pp.query("perturbation_T_infinity",tpi_Tinf);
        pp.query("perturbation_T_intensity",tpi_Ti);
        pp.query("perturbation_rho_0",tpi_rho);
        pp.query("perturbation_cp",tpi_cp);

        // Check variables message
        if (tpi_layers < 0) { amrex::Abort("Please provide a valid perturbation layer value (ie. 3-5)"); }
        if (tpi_offset < 0) { amrex::Abort("Please provide a valid inflow cell offset value for perturbation region (ie. 0-5)"); }
        for (int i = 0; i < tpi_boxDim.size(); i++) {
            if (tpi_boxDim[i] == 0) { amrex::Abort("Please provide valid dimensions for perturbation boxes."); }
        }
        if (tpi_nonDim < 0.) { amrex::Abort("Please provide a valid nondimensional number (ie. Ri = 0.042)"); }
        if (tpi_Tinf < 0.) { amrex::Abort("Please provide a valid ambient temperature value (ie. T_0)"); }
        if (tpi_Ti <= 0.) { amrex::Abort("Please provide a valid temperature intensity value (ie. 0-1.0)"); }
        if (tpi_rho < 0.) { amrex::Abort("Please provide a valid density (ie. rho_0)"); }
        if (tpi_cp < 0.) { amrex::Abort("Please provide a valid specific head capacity"); }

        // Creating perturbation region
        amrex::IntVect regionLo;
        amrex::IntVect regionHi;
        amrex::IntVect boxSize(tpi_boxDim[0],tpi_boxDim[1],tpi_boxDim[2]); // boxSize for individual boxes

        // X-direction perturbation
        if (tpi_direction[0]) {
            amrex::IntVect tmpLo(0+tpi_offset, 0, 0);
            amrex::IntVect tmpHi(((tpi_layers*tpi_boxDim[0])-1)+tpi_offset, nx[1], nx[2]);
            regionLo = tmpLo;
            regionHi = tmpHi;
        }

        // Y-direction Perturbation
        if (tpi_direction[1]) {
            amrex::IntVect tmpLo(0, 0+tpi_offset, 0);
            amrex::IntVect tmpHi(nx[0], ((tpi_layers*tpi_boxDim[1])-1)+tpi_offset, nx[2]);
            regionLo = tmpLo;
            regionHi = tmpHi;
        }

        if(tpi_direction[2]) { amrex::Abort("Currently not supporting z-direction flow perturbation"); }

        // Creating structure box array for conserved quantity
        amrex::Box tmp_bx(regionLo, regionHi, amrex::IntVect(0,0,0));
        amrex::BoxArray tmp_ba(tmp_bx);
        tmp_ba.maxSize(boxSize);
        pb_ba.push_back(tmp_ba);

        // Initializing mean magnitude vector
        pb_mag.resize(pb_ba[lev].size(), 0.);

        // Set size of vector and initialize
        pb_interval.resize(pb_ba[lev].size(), -1.0);
        pb_local_etime.resize(pb_ba[lev].size(), 0.0);
        pb_amp.resize(pb_ba[lev].size(), 0.0);

        // Computing perturbation reference length
        tpi_Lpb = tpi_boxDim[0]*dx[0];
        tpi_Wpb = tpi_boxDim[1]*dx[1];
        tpi_Hpb = tpi_boxDim[2]*dx[2];
        tpi_lref = sqrt(tpi_Lpb*tpi_Lpb + tpi_Wpb*tpi_Wpb);

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
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_rho_0: " << tpi_rho << "\n";
        amrex::PrintToFile("BoxPerturbationOutput") << "perturbation_cp: " << tpi_cp << "\n\n";

        amrex::PrintToFile("BoxPerturbationOutput") << "Reference length per box = " << tpi_lref << "\n\n";

        amrex::PrintToFile("BoxPerturbationOutput") << "Turbulent perturbation BoxArray:\n" << pb_ba[lev] << "\n";
    }

    // Perturbation update frequency check.
    // This function will trigger calc_tpi_meanMag() and calc_tpi_amp(), when
    // the local elapsed time is greater than the perturbation frequency
    void calc_tpi_update (const int lev,
                          const amrex::Real dt,
                          amrex::MultiFab& mf_xvel,
                          amrex::MultiFab& mf_yvel)
    {
        // Setting random number generator for update interval
        srand( (unsigned) time(NULL) );

        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            if ( pb_interval[boxIdx] <= pb_local_etime[boxIdx] ) {
                // Compute mean velocity magnitude within perturbation box
                calc_tpi_meanMag_perBox(boxIdx, lev, mf_xvel, mf_yvel);
                //calc_tpi_meanMag_firstCell(boxIdx, lev, mf_xvel, mf_yvel);

                if (pb_mag[boxIdx] !=0.) {
                    amrex::Real interval = tpi_lref / pb_mag[boxIdx];
                    pb_interval[boxIdx] = RandomReal(0.9*interval,1.1*interval); // 10% variation
                }

                // Trigger amplitude calculation per perturbation box
                calc_tpi_amp(boxIdx, pb_interval[boxIdx]);

                // Reset local elapsed time
                pb_local_etime[boxIdx] = 0.;
            } else {
                // Increase by timestep
                pb_local_etime[boxIdx] += dt;
            } // if
        } // for
    }

//#define ECKERT_FORMULATION
#define BULK_RICHARDSON_FORMULATION

    // Perturbation amplitude calculation
    void calc_tpi_amp (const int boxIdx, const amrex::Real interval)
    {
        amrex::Real Um = pb_mag[boxIdx];
        pb_amp[boxIdx] = 0.; // Safety step

        #ifdef ECKERT_FORMULATION
        // Muñoz-Esparza et al. (2015) Eq. 9
        pb_amp[boxIdx] = Um * Um / (tpi_rho * tpi_cp * tpi_nonDim);
        #endif

        #ifdef BULK_RICHARDSON_FORMULATION
        amrex::Real beta = 1./tpi_Tinf;

        // Pseudo Random temperature the ignores scale when mechanically tripping turbulence
        amrex::Real g = (tpi_nonDim * Um * Um) / (tpi_Ti * tpi_Hpb);
        //amrex::Real g = CONST_GRAV;

        // Ma and Senocak (2023) Eq. 8, solving for delta phi
        if (g != 0) pb_amp[boxIdx] = (tpi_nonDim * Um * Um) / (g * beta * tpi_Hpb);
        //amrex::Print()<<"tpi_nonDim: "<<tpi_nonDim<<" Um: "<<Um<<" g: "<<g<<" beta: "<<beta<<" tpi_Hpb: "<<tpi_Hpb<<" theta["<<boxIdx<<"]: "<<pb_amp[boxIdx]
        //              <<" interval: "<<interval<<" amp: "<<pb_amp[boxIdx]/interval<< "\n";
        #endif

        // Ma & Senocak (2023) Eq. 7
        pb_amp[boxIdx] /= interval;

    }

//#define INDEX_PERTURB

    // Applying perturbation amplitude onto source term (Umphrey and Senocak 2016)
    // Random perturbation is applied as a white noise on the temperature source term
    // it becomes colored through the filtering of temperature transport equation from
    // the eddy diffusivity.
    void apply_tpi (const int lev,
                    const amrex::Box& vbx,                     // box union from upper level
                    const int comp,                            // Component to modify
                    const amrex::IndexType m_ixtype,           // IntVect type of src_arr
                    const amrex::Array4<amrex::Real>& src_arr) // Array to apply perturbation
    {
        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            amrex::Box pbx = convert(pb_ba[lev][boxIdx], m_ixtype);
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                amrex::Real amp_copy = pb_amp[boxIdx];
                ParallelForRNG(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                    // Random number generation between 0 and 1
                    amrex::Real rand_double = amrex::Random(engine);

                    // Rescale between -1 and 1 and multiple with amplitude
                    src_arr(i,j,k,comp) += (rand_double*2.0 - 1.0) * amp_copy;

                    // For box region debug only
                    #ifdef INDEX_PERTURB
                    src_arr(i,j,k,comp) = (amrex::Real) (boxIdx + 5.);
                    #endif
                });
            }
        }
    }

// TODO: Test the difference between these two
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

    // This emulates calc_tpi_meanMag_perBox() but is heavily simplified for
    // Eckert number formulation
    void calc_tpi_meanMag_firstCell (const int boxIdx,
                                     const int lev,
                                     amrex::MultiFab& mf_xvel,
                                     amrex::MultiFab& mf_yvel)
    {
        // Creating local copy of PB box array and magnitude
        const amrex::BoxArray m_pb_ba = pb_ba[lev];
        amrex::Real* m_pb_mag = get_pb_mag();

        // Storage of averages per PB
        // Index: 0=u (vol/slab_lo), 1=v (vol/slab_lo)
        amrex::Gpu::DeviceVector<amrex::Real> tmp_d(2,0.);
        amrex::Real* tmp = tmp_d.data();

        // Averaging u components
        for (amrex::MFIter mfi(mf_xvel, TileNoZ()) ; mfi.isValid(); ++mfi) {
            const amrex::Box &vbx = mfi.validbox();
            // This is hardcoded so that it's always on the 1st cell level (May need to generalize for EB)
            amrex::Box pbx(amrex::IntVect(m_pb_ba[boxIdx].smallEnd(0),m_pb_ba[boxIdx].smallEnd(1),0),
                           amrex::IntVect(m_pb_ba[boxIdx].bigEnd(0)  ,m_pb_ba[boxIdx].bigEnd(1)  ,1),
                           vbx.ixType());
            amrex::Box ubx = pbx & vbx;

            // Operation over box union
            if (ubx.ok()) {
                const amrex::Array4<const amrex::Real> &xvel_arry = mf_xvel.array(mfi);
                int npts = ubx.numPts();
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[0], xvel_arry(i,j,k), handler);
                });
                tmp[0] /= (amrex::Real) npts;
            } // if
        } // MFIter

        // Averaging v components
        for (amrex::MFIter mfi(mf_yvel, TileNoZ()) ; mfi.isValid(); ++mfi) {
            const amrex::Box &vbx = mfi.validbox();
            // This is hardcoded so that it's always on the 1st cell level (May need to generalize for EB)
            amrex::Box pbx(amrex::IntVect(m_pb_ba[boxIdx].smallEnd(0),m_pb_ba[boxIdx].smallEnd(1),0),
                           amrex::IntVect(m_pb_ba[boxIdx].bigEnd(0)  ,m_pb_ba[boxIdx].bigEnd(1)  ,1),
                           vbx.ixType());
            amrex::Box ubx = pbx & vbx;

            // Operation over box union
            if (ubx.ok()) {
                const amrex::Array4<const amrex::Real> &yvel_arry = mf_yvel.array(mfi);
                int npts = ubx.numPts();
                ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), ubx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept {
                    amrex::Gpu::deviceReduceSum(&tmp[1], yvel_arry(i,j,k), handler);
                });
                tmp[1] /= (amrex::Real) npts;
            } // if
        } // MFIter

        // Computing the average magnitude within 1st z level slab
        m_pb_mag[boxIdx] = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
    }

    // Quality of life function definitions
    amrex::Real* get_pb_mag() { return pb_mag.data(); }

    // Output debug message into a file
    void debug (amrex::Real time)
    {
        amrex::PrintToFile("BoxPerturbationOutput") << "#################### PB output at time = "
            << time << " ####################\n";

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

    // Public data members
    amrex::Vector<amrex::BoxArray> pb_ba;  // PB box array
    amrex::Vector<amrex::Real>     pb_mag; // BP mean magnitude [m/s]

  private:

    // Random number generation between range (used for interval calculation)
    amrex::Real RandomReal (const amrex::Real min, const amrex::Real max)
    {
        amrex::Real r = (amrex::Real) rand() / (amrex::Real) RAND_MAX;
        return min + r * (max - min);
    }

    // Private data members
    int tpi_layers;                   // Number of layers of perturbation boxes
    int tpi_offset;                   // Cells to offset the start of the perturbation region

    amrex::Vector<int> tpi_boxDim;    // Dimensions of each perturbation box
    amrex::Vector<int> tpi_direction; // Direction of the peturbation

    amrex::Real tpi_nonDim;           // Nondimensional Value

    // Richardson Formulation
    amrex::Real tpi_Ti;               // Temperature intensity value
    amrex::Real tpi_Tinf;             // Reference temperature [K]

    // Eckert Formulation
    amrex::Real tpi_rho;              // Reference density
    amrex::Real tpi_cp;               // Specific heat capaticty

    // Physical dimensions
    amrex::Real tpi_Hpb;              // PB height [m]
    amrex::Real tpi_Lpb;              // PB length [m]
    amrex::Real tpi_Wpb;              // PB width  [m]
    amrex::Real tpi_lref;             // PB reference length [m]

    // Perturbation data arrays
    amrex::Vector<amrex::Real> pb_interval;        // PB update time [s]
    amrex::Vector<amrex::Real> pb_local_etime; // PB local elapsed time [s]
    amrex::Vector<amrex::Real> pb_amp;         // PB perturbation amplitude Ri:[K] Ec:[K/s]

};
#endif
