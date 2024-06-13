#ifndef _TURB_PERT_STRUCT_H_
#define _TURB_PERT_STRUCT_H_

/**
 * Container holding quantities related to turbulent perturbation parameters
 */

#define DEBUG_PERTBOX_MSG

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

        // TODO: This will need to be changed later on
        tpi_Lpb = 8*dx[0];
        tpi_Wpb = 8*dx[1];
        tpi_Hpb = 8*dx[2];

        #ifdef DIAGONAL_REF
        tpi_lref = sqrt(8*8*dx[0]*dx[0] + 8*8*dx[1]*dx[1]);
        #elif defined(LENGTH_REF)
        tpi_lref = 8*dx[0];
        #endif

        // Creating structure box array for conserved quantity
        amrex::Box tmp_bx(regionLo, regionHi, amrex::IntVect(0,0,0));
        amrex::BoxArray tmp_ba(tmp_bx);
        tmp_ba.maxSize(boxSize);
        pb_ba.push_back(tmp_ba);

        // Initializing mean U vector
        pb_um.resize(pb_ba[lev].size(), 0.);
        pb_vm.resize(pb_ba[lev].size(), 0.);

        // Set size of vector and initialize
        pb_freq.resize(pb_ba[lev].size(), -1.0);
        pb_local_etime.resize(pb_ba[lev].size(), 0.0);
        pb_amp.resize(pb_ba[lev].size(), 0.0);

        // Function finalizing output message
        amrex::Print() << "Turbulent perturbation region initialized with ba " << pb_ba[lev].size() << " boxes\n";

        amrex::Print() << "Using ";
        #ifdef DIAGONAL_REF
        amrex::Print() << "box diagonal ";
        #elif defined(LENGTH_REF)
        amrex::Print() << "box length ";
        #endif
        amrex::Print() << "distance as reference for turbulent inlet perturbation: tpi_lref = " << tpi_lref << "\n";
    }

    // perturbation box update freqnecy is calculate and perturbation amplitude
    // function is triggered through local elapsed time
    void calc_tpi_update (const int lev,
                          const amrex::Real dt)
    {
        for (int boxIdx = 0; boxIdx < pb_ba[lev].size(); boxIdx++) {
            if ( (pb_freq[boxIdx] <= pb_local_etime[boxIdx]) && 
                ((pb_um[boxIdx] != 0.) || (pb_vm[boxIdx] != 0.)) ) { // TODO: Second part of condition may need to change
                // Ma and Senocak (2015) Eq. 15
                pb_freq[boxIdx] = tpi_lref / sqrt(pb_um[boxIdx]*pb_um[boxIdx] + pb_vm[boxIdx]*pb_vm[boxIdx]);

                // Trigger amplitude calculation per perturbation box
                calc_tpi_amp(boxIdx);

                // Reset local elapsed time
                pb_local_etime[boxIdx] = 0.;
            } else {
                // Increase by timestep
                pb_local_etime[boxIdx] += dt;
            } // if

            #ifdef DEBUG_PERTBOX_MSG
            amrex::Print() << "pb_freq[" << boxIdx << "] = " << pb_freq[boxIdx] << "\n";
            #endif
        } // for
    }

//#define EKERT_FORMULATION
#define BULK_RICHARDSON_FORMULATION 
//#define GRAD_RICHARDSON_FORMULATION 

    // Perturbation amplitude calcluation
    // TODO: This is where all of the different ways of computing perturbation amplitude will be defined
    void calc_tpi_amp (const int boxIdx)
    {
        beta = 1./tpi_Tinfty;

        Um = pb_um[boxIdx];
        Vm = pb_vm[boxIdx];

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
                ParallelForRNG(ubx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                    // Random number generation between 0 and 1
                    amrex::Real rand_double = amrex::Random(engine);

                    // Rescale between -1 and 1 and multiple with amplitude
                    src_arr(i,j,k,comp) = (rand_double*2.0 - 1.0) * pb_amp[boxIdx];

                    // For box region debug only
                    #ifdef INDEX_PERTURB
                    src_arr(i,j,k,comp) = (amrex::Real) (boxIdx +1.);
                    #endif
                });
            }
        }
    }

    amrex::Real* get_pb_um() { return pb_um.data(); }

    void debug()
    {
        for (int i=0; i < pb_um.size(); i++){
            amrex::Print()<< "pb_um[" << i << "] = " << pb_um[i] << "\n";
        }
    }
 
    // Public data members 
    amrex::Vector<amrex::BoxArray> pb_ba;

    // Perturbation box mean velocities
    amrex::Vector<amrex::Real> pb_um;
    amrex::Vector<amrex::Real> pb_vm;

    // Public constant constants
    amrex::Real tpi_Tinfty;
    amrex::Real tpi_nonDim;

  private:

    // Private constant constants
    amrex::Real g = CONST_GRAV;
    amrex::Real beta;

    amrex::Real tpi_Hpb;
    amrex::Real tpi_Lpb;
    amrex::Real tpi_Wpb;

    amrex::Real tpi_lref;

    amrex::Real Um;
    amrex::Real Vm;

    // Perturbation update time
    amrex::Vector<amrex::Real> pb_freq;
    amrex::Vector<amrex::Real> pb_local_etime;
    amrex::Vector<amrex::Real> pb_amp;

};
#endif
