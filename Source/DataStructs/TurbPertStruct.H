#ifndef _TURB_PERT_STRUCT_H_
#define _TURB_PERT_STRUCT_H_

#include <AMReX_MultiFabUtil.H>
#include <TileNoZ.H>

/**
 * Container holding quantities related to turbulent perturbation parameters
 */

struct TurbulentPerturbation {

  public:

    ~TurbulentPerturbation () {
        amrex::Print() << "~TurbulentPerturbation destructor\n";
    }

//#define DIAGONAL_REF
#define LENGTH_REF

/*    void init_tpi_const(const amrex::Real T_0)
    {
       // Bulk Richardson number
       tpi_nonDim = 0.042;

       // Ambient Temperature
       tpi_Tinfty = T_0;
    }*/

    // Initializing Perturbation Region
    // Note: Box will always be initialized at cell centered
    void init_tpi (const int lev,
                   const amrex::IntVect& nx,
                   amrex::GpuArray<amrex::Real,3> const dx)

    {
        // TODO: Use user Param to define this (Currently HARDCODED)
        // Perturbation box region setup
        int pertBox_offset = 0;
        amrex::IntVect regionLo(0+pertBox_offset, 0, 0);
        amrex::IntVect regionHi(7+pertBox_offset, nx[1], nx[2]);
        amrex::IntVect boxSize(8,8,8);
        //amrex::IntVect boxSize(4,4,4); // HACK

        // TODO: This will need to be changed later on
        tpi_Lpb = 8*dx[0];
        tpi_Wpb = 8*dx[1];
        tpi_Hpb = 8*dx[2];

        #ifdef DIAGONAL_REF
        tpi_lref = sqrt(boxSize(0)*boxSize(0)*dx[0]*dx[0] + boxSize(1)*boxSize(1)*dx[1]*dx[1]);
        #elif defined(LENGTH_REF)
        tpi_lref = 8*dx[0];
        #endif

        // Creating structure box array for conserved quantity
        amrex::Box tmp_bx(regionLo, regionHi, amrex::IntVect(0,0,0));
        amrex::BoxArray tmp_ba(tmp_bx);
        tmp_ba.maxSize(boxSize);
        pb_ba.push_back(tmp_ba);

        // Initializing mean magnitude vector
        pb_mag.resize(pb_ba[lev].size(), 0.);

        // Set size of vector and initialize
        pb_freq.resize(pb_ba[lev].size(), -1.0);
        pb_local_etime.resize(pb_ba[lev].size(), 0.0);
        pb_amp.resize(pb_ba[lev].size(), 0.0);

        // Function check point message
        amrex::Print() << "Turbulent perturbation region initialized with ba " << pb_ba[lev].size() << " boxes\n";

        // Function method check point message
        amrex::Print() << "Using ";
        #ifdef DIAGONAL_REF
        amrex::Print() << "box diagonal ";
        #elif defined(LENGTH_REF)
        amrex::Print() << "box length ";
        #endif
        amrex::Print() << "distance as reference for turbulent inlet perturbation: tpi_lref = " << tpi_lref << "\n";
    }

    // Perturbation update frequency check.
    // This function will trigger calc_tpi_meanMag() and calc_tpi_amp(), when
    // the local elasped time is greater than the perturbation frequency
    void calc_tpi_update (const int lev,
                          const amrex::Real dt,
                          amrex::MultiFab& mf_xvel,
                          amrex::MultiFab& mf_yvel,
                          amrex::MultiFab& mf_cons)
    {
        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            if ( pb_freq[boxIdx] <= pb_local_etime[boxIdx] ) {
                // Compute mean velocity magnitude within perturbation box
                calc_tpi_meanMag(boxIdx, lev, mf_xvel, mf_yvel, mf_cons);

                // Ma and Senocak (2015) Eq. 15
                if (pb_mag[boxIdx] !=0.) {
                    pb_freq[boxIdx] = tpi_lref / pb_mag[boxIdx];
                }

                // Trigger amplitude calculation per perturbation box
                calc_tpi_amp(boxIdx);

                // Reset local elapsed time
                pb_local_etime[boxIdx] = 0.;
            } else {
                // Increase by timestep
                pb_local_etime[boxIdx] += dt;
            } // if
        } // for
    }

// TODO: Test the difference between these two
//#define USE_SLAB_AVERAGE
#define USE_VOLUME_AVERAGE

    // Perturbation box mean velocity magnitude calculation
    // This is pulled into the strucutre to also utilize during runtime
    void calc_tpi_meanMag (const int boxIdx,
                           const int lev,
                           amrex::MultiFab& mf_xvel,
                           amrex::MultiFab& mf_yvel,
                           amrex::MultiFab& mf_cons)

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
                #endif

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
                #endif
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
                #endif

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
                #endif
            } // if
        } // MFIter

        // Computing the average magnitude within PB
        for (amrex::MFIter mfi(mf_cons, TileNoZ()); mfi.isValid(); ++mfi) {
            amrex::Box const& vbx = mfi.validbox();
            amrex::Box pbx = convert(m_pb_ba[boxIdx], vbx.ixType());
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                #ifdef USE_SLAB_AVERAGE
                m_pb_mag[boxIdx] = 0.5*(sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]) + sqrt(tmp[2]*tmp[2] + tmp[3]*tmp[3]));
                #endif

                #ifdef USE_VOLUME_AVERAGE
                m_pb_mag[boxIdx] = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
                #endif
            } // if
        } // MFIter

    }

//#define ECKERT_FORMULATION
#define BULK_RICHARDSON_FORMULATION

    // Perturbation amplitude calcluation
    void calc_tpi_amp (const int boxIdx)
    {
        amrex::Real g = CONST_GRAV;
        amrex::Real beta = 1./tpi_Tinfty;
        amrex::Real Um = pb_mag[boxIdx];

        // Reset the perturbation amplitude
        pb_amp[boxIdx] = 0.;

        #ifdef BULK_RICHARDSON_FORMULATION
        // Ma and Senocak (2023) Eq. 17
        pb_amp[boxIdx] = tpi_nonDim * Um * Um * Um / (g * beta * tpi_lref * tpi_Hpb);
        #endif

    }

//#define INDEX_PERTURB

    // Applying perturbation amplitude onto source term (Umphrey and Senocak 2016)
    // Random perturbation is applied as a white noise on the temperature source term
    // it becomes colored through the filtering of temperature transport equation from
    // the eddy diffusivity.
    void apply_tpi (const int lev,
                    const amrex::Box& vbx,                    // box union from upper level
                    const int comp,                           // Component to modify
                    const amrex::IndexType m_ixtype,          // IntVect type of src_arr
                    amrex::Array4<amrex::Real> const& src_arr)// Source Array to apply perturbation
    {
        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            amrex::Box pbx = convert(pb_ba[lev][boxIdx], m_ixtype);
            amrex::Box ubx = pbx & vbx;
            if (ubx.ok()) {
                amrex::Real amp_copy = pb_amp[boxIdx]; // This exposes a copy so the GPU can capture by value
                ParallelForRNG(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                    // Random number generation between 0 and 1
                    amrex::Real rand_double = amrex::Random(engine);

                    // Rescale between -1 and 1 and multiple with amplitude
                    src_arr(i,j,k,comp) = (rand_double*2.0 - 1.0) * amp_copy;

                    // For box region debug only
                    #ifdef INDEX_PERTURB
                    src_arr(i,j,k,comp) = (amrex::Real) (boxIdx + 1.);
                    #endif
                });
            }
        }
    }


    // Quality of life function definitions
    amrex::Real* get_pb_mag() { return pb_mag.data(); }

    void debug ()
    {
        for (int i = 0; i < pb_mag.size(); i++) {
            amrex::Print() << "[" << i<< "] pb_mag=" << pb_mag[i]
                           << " |  pb_freq=" << pb_freq[i]
                           << " |  pb_local_etime=" << pb_local_etime[i] <<"\n";
        }
    }

    // Public constant constants
    amrex::Real tpi_Tinfty; // [K]
    amrex::Real tpi_nonDim;

    // Public data members
    amrex::Vector<amrex::BoxArray> pb_ba; // PB box array
    amrex::Vector<amrex::Real> pb_mag;    // BP mean magnitude [m/s]

  private:

    // Private constant constants
    amrex::Real tpi_Hpb;  // PB height [m]
    amrex::Real tpi_Lpb;  // PB length [m]
    amrex::Real tpi_Wpb;  // PB width  [m]
    amrex::Real tpi_lref; // PB reference length [m]

    // Storage Arrays
    amrex::Vector<amrex::Real> pb_freq;        // PB frequency [1/s]
    amrex::Vector<amrex::Real> pb_local_etime; // PB local elapsed time [s]
    amrex::Vector<amrex::Real> pb_amp;         // PB perturbation amplitude Ri:[K] Ec:[K/s]

};
#endif
