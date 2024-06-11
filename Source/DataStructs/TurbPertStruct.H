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

    // Initializing Perturbation Region
    void init_tpi (const int lev,                             // Level of current data
                                  const amrex::IntVect& nx,
                                  amrex::GpuArray<amrex::Real,3> const dx,
                                  const amrex::BoxArray& data_ba,            // Box array from AMReX domain decomposition
                                  const amrex::DistributionMapping& data_dm) // Distribution map from AMReX domain decomposition
    {
        // TODO: Use user Param to define this (Currently HARDCODED)
        // Perturbation box region setup
        int pertBox_offset = 0;
        amrex::IntVect regionLo(0+pertBox_offset, 0, 0);
        amrex::IntVect regionHi(7+pertBox_offset, nx[1], nx[2]);
        amrex::IntVect boxSize(8,8,8);

        #ifdef DIAGONAL_REF
        tpi_lref = sqrt(8*8*dx[0]*dx[0] + 8*8*dx[1]*dx[1]);
        #elif defined(LENGTH_REF)
        tpi_lref = 8*dx[0];
        #endif

        // TODO: Generalize for xbx, ybx, and zbx using amrex::convert()
        // Creating structure box array for conserved quantity
        amrex::Box cbx(regionLo, regionHi, amrex::IntVect(0,0,0)); // Cell center
        amrex::BoxArray cba(cbx);
        cba.maxSize(boxSize);

        // N^2 search to build local distribution map
        const amrex::Vector<int> d_pmap = data_dm.ProcessorMap();
        amrex::Vector<int> p_pmap(cba.size()); // Initialization of p_pmap with npbx elements

        for (int pbx(0); pbx < cba.size(); ++pbx) {
            bool found = false;

            for (int dbx(0); dbx < (int) data_ba.size(); ++dbx) {
                amrex::Box m_dbx = data_ba[dbx];

                if (m_dbx.contains(cba[pbx])) {
                    found = true;
                    p_pmap[pbx] = d_pmap[dbx];
                    continue;
                } // found
            } // data box
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(found, "Perturbation box not contained in a valid box!");
        } // pert box

        amrex::DistributionMapping dm(p_pmap);

        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "[TurbPert::init_tpi]\nmf_tpi_freq.DistributionMap():\n" << dm << "\n";
        #endif

        // TODO: The following 3 mf can probably be combined into 1
        // Creating multifab for update time
        mf_tpi_freq.define(cba, dm, 1, 0);
        mf_tpi_etime.define(cba, dm, 1, 0);
        mf_tpi_amp.define(cba, dm, 1, 0);

        // Initializing mf data
        mf_tpi_freq.setVal(1.);
        mf_tpi_etime.setVal(0.);
        mf_tpi_amp.setVal(0.);

        // Function finalizing output message
        amrex::Print() << "Turbulent perturbation region initialized with ba " << mf_tpi_freq.boxArray().size() << " boxes\n";

        amrex::Print() << "Using ";
        #ifdef DIAGONAL_REF
        amrex::Print() << "box diagonal ";
        #elif defined(LENGTH_REF)
        amrex::Print() << "box length ";
        #endif
        amrex::Print() << "distance as reference for turbulent inlet perturbation: tpi_lref = " << tpi_lref << "\n";
    }

//#define INDEX_PERTURB // Used for box debugging
#define RANDOM_PERTURB

    // TODO: Write in function to call specific box to calculate new amplitude if update time is larger than the elapsed time
    // Calculate perturbation update time
    void calc_tpi_updateTime (const int lev,
                                   const amrex::Real dt_lev0,
                                   const amrex::Box& bx,
                                   amrex::Array4<amrex::Real> const& xvel,
                                   amrex::Array4<amrex::Real> const& yvel)
    {
       // Reading in Multifab boxArray
       amrex::BoxArray ba = mf_tpi_freq.boxArray();

       for (int boxIdx = 0; boxIdx < ba.size(); boxIdx++) {
           //bool check = false;

           // Checking if perturbation boxes are within current rank box
           if (bx.contains(ba[boxIdx])) {
               // Reading in current boxIdx for perturbation boxes
               amrex::Array4<amrex::Real> const& pb_freq  = mf_tpi_freq.array(boxIdx);
               amrex::Array4<amrex::Real> const& pb_etime = mf_tpi_etime.array(boxIdx);

               amrex::Array4<amrex::Real> const& pb_amp   = mf_tpi_amp.array(boxIdx); // XXX HACK

               ParallelFor(ba[boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   #ifdef RANDOM_PERTURB
                   if (pb_freq(i,j,k,0) >= pb_etime(i,j,k,0)) { 
                       // Eq. 10 from Muñoz-Esparza et al. (2015)
                       pb_freq(i,j,k,0) = tpi_lref / sqrt(xvel(i,j,k,0)*xvel(i,j,k,0) + yvel(i,j,k,0*yvel(i,j,k,0)));

                       pb_amp(i,j,k,0) = pb_freq(i,j,k,0); // XXX HACK

                       // Set elapsed time back to 0.
                       pb_etime(i,j,k,0) = 0.;

                       // Bool to trigger new amplitdue calculation
                       //check = true;
                   }
                   #endif

                   #ifdef INDEX_PERTURB
                   pb_freq(i,j,k,0) = (amrex::Real) boxIdx+1.;
                   #endif

                   // Incrementing time value for each box
                   pb_etime(i,j,k,0) += dt_lev0;
               });

               // Calculate new amplitude
               //if (check) calc_tpi_amp(lev, bx, ba, 0);
           } // if bx.contain()
       } // for
    }

    //TODO: Inside of this function I am operating from ba being already within the geom box
    // therfore another would be redudant. We can directly access data from u and v velocity
    // to compute (1) average of velocity between top and bottom perturbaiton box planes, or
    // (2) compute the box volume average of the velocity components.
    // Will need to find out how to compute this plane average first.
    void calc_tpi_amp (const int lev,         // Current level
                       const amrex::Box& bx,  // TileBox bx
                       const amrex::BoxArray& ba, 
                       const int comp)        // Component to modify
    {

    }

    // TODO: Access velocity data to calculate perturbation based on simulation values
    // Calculate perturbation amplitude
    void apply_tpi (const int lev,                            // Current level
                    const amrex::Box& bx,                     // TileBox bx
                    const int comp,                           // Component to modify
                    amrex::Array4<amrex::Real> const& src_arr)// Source Array to apply perturbation
    {
        amrex::BoxArray ba = mf_tpi_amp.boxArray();
        for (int boxIdx = 0; boxIdx < ba.size(); boxIdx++) {

            // Populating array with values
            if (bx.contains(ba[boxIdx])) {
                amrex::Array4<amrex::Real> const& pb_amp = mf_tpi_amp.array(boxIdx);

                #ifdef DEBUG_PERTBOX_MSG
                amrex::Print() << "  bx: " << bx << " -- contains perturbation box #" << boxIdx << ": " << ba[boxIdx] << "\n";
                #endif

                // Applying random perturbation to perturbation amplitude
                ParallelForRNG(ba[boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                    #ifdef RANDOM_PERTURB
                    amrex::Real rand_double = amrex::Random(engine);  // Between 0.0 and 1.0
                    src_arr(i,j,k,comp) = (rand_double*2.0 - 1.0) * pb_amp(i,j,k,0);
                    #endif

                    #ifdef INDEX_PERTURB
                    src_arr(i, j, k, comp) = pb_amp(i,j,k,0);
                    #endif
                });
            } // if
        } // boxIdx
    }


  private:

    // Turbulent perturbation inleft reference length
    amrex::Real tpi_lref;

    // Perturbation box mf
    amrex::MultiFab mf_tpi_freq;  // Update Time
    amrex::MultiFab mf_tpi_etime; // Elapsed Time
    amrex::MultiFab mf_tpi_amp;   // Amplitdue

};
#endif
