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

//#define MANUAL_DISTRIBUTION_MAP_ASSIGNMENT

    // Initializing Perturbation Region
    void init_PerturbationRegion (const int lev,                             // Level of current data
                                  const amrex::IntVect& nx,
                                  const amrex::BoxArray& data_ba,            // Box array from AMReX domain decomposition
                                  const amrex::DistributionMapping& data_dm) // Distribution map from AMReX domain decomposition
    {
        // TODO: Use user Param to define this (Currently HARDCODED)
        // Perturbation box region setup
        int pertBox_offset = 0;
        amrex::IntVect regionLo(0+pertBox_offset, 0, 0);
        amrex::IntVect regionHi(7+pertBox_offset, nx[1], nx[2]);
        amrex::IntVect boxSize(8,8,8);

        // TODO: Generalize for xbx, ybx, and zbx using amrex::convert()
        // Creating structure box array for conserved quantity
        amrex::Box cbx(regionLo, regionHi, amrex::IntVect(0,0,0)); // Cell center
        amrex::BoxArray cba(cbx);
        cba.maxSize(boxSize);

        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "[TurbPert::init_PerturbationRegion] C subdividing into smaller boxes: "
                       << cba.size() << "\n" << cba;
        #endif

        #ifdef MANUAL_DISTRIBUTION_MAP_ASSIGNMENT
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
        #else
        // TODO: Once mf is accessing data (ie. u, and v velocitices), come back and test if this fails... for science
        amrex::DistributionMapping dm;
        dm.define(cba,2);
        #endif

        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "[TurbPert::init_PerturbationRegion]\nmf_perturb_updateTime.DistributionMap():\n" << dm << "\n";
        #endif

        // Creating multifab for update time
        mf_perturb_updateTime.define(cba, dm, 1, 0); // Storing only update time
        mf_perturb_updateTime.setVal(-1.);
    }

    // TODO: Write in function to call specific box to calculate new amplitude if update time is larger than the elapsed time
    // Calculate perturbation update time
    void calc_TurbPert_updateTime (const int lev,
                                   const amrex::Box& bx,
                                   amrex::Array4<amrex::Real> const& xvel,
                                   amrex::Array4<amrex::Real> const& yvel)
    {
       // Reading in Multifab boxArray
       amrex::BoxArray ba = mf_perturb_updateTime.boxArray();

       for (int boxIdx = 0; boxIdx < ba.size(); boxIdx++) {
           // Reading in current boxIdx for perturbation boxes
           amrex::Array4<amrex::Real> const& pb_t = mf_perturb_updateTime.array(boxIdx);

           // Checking if perturbation boxes are within current rank box
           if (bx.contains(ba[boxIdx])) {
               ParallelFor(ba[boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   if (pb_t(i,j,k,0) < 0.) { 
                       //pb_t(i,j,k,0) = sqrt(xvel(i,j,k,0)*xvel(i,j,k,0) + yvel(i,j,k,0*yvel(i,j,k,0)));
                       pb_t(i,j,k,0) = (amrex::Real) boxIdx+1.;
                   }

                   #ifdef DEBUG_PERTBOX_MSG
                   /*
                   amrex::Print() << "bx = " << bx << " ba[" << boxIdx << "] = " << ba[boxIdx]
                                  << " ----- bx.contain(ba[boxIdx]) = " << bx.contains(ba[boxIdx])
                                  << " xvel = " << xvel(i,j,k,0) << " yvel = " << yvel(i,j,k,0)
                                  << " updateTime = " << pb_t(i,j,k,0) << "\n";
                   */
                   #endif
               });

           } // if

       } // for

    }

#define INDEX_PERTURB
//#define RANDOM_PERTURB

    // TODO: Access velocity data to calculate perturbation based on simulation values
    // Calculate perturbation amplitude
    void calc_TurbPert_amplitude (const int lev,                            // Current level
                                  const amrex::Box& bx,                     // TileBox bx
                                  const int comp,                           // Component to modify
                                  amrex::Array4<amrex::Real> const& src_arr)// Source Array to apply perturbation
    {
        amrex::BoxArray ba = mf_perturb_updateTime.boxArray();
        for (int boxIdx = 0; boxIdx < ba.size(); boxIdx++) {
            amrex::Array4<amrex::Real> const& pb_t = mf_perturb_updateTime.array(boxIdx);

            // Populating array with values
            if (bx.contains(ba[boxIdx])) {
                #ifdef DEBUG_PERTBOX_MSG
                amrex::Print() << "  bx: " << bx << " -- contains perturbation box #" << boxIdx << ": " << ba[boxIdx] << "\n";
                #endif

                #ifdef INDEX_PERTURB
                // This exists for debugging purpose
                ParallelFor(ba[boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    src_arr(i, j, k, comp) = pb_t(i,j,k,0);
                });
                #endif

                #ifdef RANDOM_PERTURB
                ParallelForRNG(ba[boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                    amrex::Real rand_double = amrex::Random(engine);  // Between 0.0 and 1.0
                    src_arr(i,j,k,comp) = (rand_double*2.0 - 1.0) * pb_t(i,j,k,0);
                });
                #endif

            } // if

        } // boxIdx

    }

    // Perturbation box mf
    amrex::MultiFab mf_perturb_updateTime;

  private:

};
#endif
