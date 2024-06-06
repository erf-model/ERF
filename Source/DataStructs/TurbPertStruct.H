#ifndef _TURB_PERT_STRUCT_H_
#define _TURB_PERT_STRUCT_H_

/**
 * Container holding quantities related to turbulent perturbation parameters
 */

#include <AMReX_MultiFab.H>

#define DEBUG_PERTBOX_MSG

#define MANUAL_DISTRIBUTION_MAP_ASSIGNMENT

struct TurbulentPerturbation {

  public:

    ~TurbulentPerturbation () { 
        amrex::Print() << "~TurbulentPerturbation destructor\n";
    }

    // Initializing Perturbation Region
    void init_PerturbationRegion (const int lev,                             // Level of current data
                                  const amrex::BoxArray& data_ba,            // Box array from AMReX domain decomposition
                                  const amrex::DistributionMapping& data_dm) // Distribution map from AMReX domain decomposition
    {
        // TODO: Need to use user Param to define this (Currently HARDCODED)
        // Perturbation box region setup
        int pertBox_offset = 0;
        amrex::IntVect regionLo(0+pertBox_offset, 0, 0);
        amrex::IntVect regionHi(7+pertBox_offset,31,31);

        // Defining perturbation box size
        amrex::IntVect boxSize(8,8,8);

        // Creating structure box array for conserved quantity
        amrex::Box tmp_cbx(regionLo, regionHi, amrex::IntVect(0,0,0)); // Cell center
        amrex::BoxArray tmp_cba(tmp_cbx);
        turbPert_cba.push_back(tmp_cba);
        turbPert_cba[lev].maxSize(boxSize);
        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "[TurbPert::init_PerturbationRegion] C subdividing into smaller boxes: "
                       << turbPert_cba[lev].size() << "\n"
                       << "turbPert_cba[" << lev << "] constains: " << turbPert_cba[lev];
        #endif

        // Creating structure box array for z velocity
        amrex::BoxArray tmp_zba = amrex::convert(tmp_cba, amrex::IntVect(0,0,1)); // Convert for Z
        turbPert_zba.push_back(tmp_zba);
        turbPert_zba[lev].maxSize(boxSize);
        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "[TurbPert::init_PerturbationRegion] Z subdividing into smaller boxes: "
                       << turbPert_zba[lev].size() << "\n"
                       << "turbPert_zba[" << lev << "] constains: " << turbPert_zba[lev];
        #endif

        // Initialize perturbation update time and amplitude array
        pb_updateTime.resize(turbPert_cba[lev].size());
        pb_amplitude.resize(turbPert_cba[lev].size());

        amrex:: Print() << "Turbulence perturbation region initialized\n";
    }

#define LENGTH_SCALING
#define WIDTH_SCALING
#define DIAGONAL_SCALING

    // TODO: Write in function to call specific box to calculate new amplitude if update time is larger than the elapsed time
    // Calculate perturbation update time
    void calc_TurbPert_updateTime (const int lev, amrex::Vector<amrex::Real>& pb_t)
    {
       for (int boxIdx = 0; boxIdx < (int) pb_t.size(); boxIdx++) {
           pb_t[boxIdx] = (amrex::Real) (boxIdx+1.0);
       }
    }

#define INDEX_PERTURB
//#define RANDOM_PERTURB

    // TODO: Access velocity data to calculate perturbation based on simulation values
    // Calculate perturbation amplitude
    void calc_TurbPert_amplitude (const int lev,                            // Current level
                                  const amrex::Box& bx,                     // TileBox bx
                                  const amrex::Vector<amrex::BoxArray>& ba, // Perturbation boxArray
                                  amrex::Vector<amrex::Real>& amp,
                                  const int comp,                           // Component to modify
                                  const auto &src_arr)                      // Source Array to apply perturbation
    {
        for (int boxIdx = 0; boxIdx < (int) ba[lev].size(); boxIdx++) {
            // Populating array with values
            amp[boxIdx] = (amrex::Real) (boxIdx+1.0);
            if (bx.contains(ba[lev][boxIdx])) {
                #ifdef DEBUG_PERTBOX_MSG
                amrex::Print() << "  bx: " << bx << " -- contains perturbation box #" << boxIdx << ": " << ba[lev][boxIdx] << "\n";
                #endif

                #ifdef INDEX_PERTURB
                // This exists for debugging purpose
                ParallelFor(ba[lev][boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    //amrex::Real tmp_Src = (i + j*bx.size()[0] + k*bx.size()[0]*bx.size()[1]); // continuous count
                    amrex::Real tmp_Src = amp[boxIdx]; // distinguish each box
                    src_arr(i, j, k, comp) = tmp_Src;
                });
                #endif

                #ifdef RANDOM_PERTURB
                ParallelForRNG(ba[lev][boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                    amrex::Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
                    amrex::Real tmp_Src = amp[boxIdx] * (rand_double*2.0 - 1.0); // Automatically rescales between -1.0 and 1.0
                    src_arr(i,j,k,comp) = tmp_Src;
                });
                #endif
            } // if
        } // boxIdx
    }

    // Turbulent perturbation box array
    amrex::Vector<amrex::BoxArray> turbPert_cba;
    amrex::Vector<amrex::BoxArray> turbPert_zba;

    // Perturbation box quantities
    amrex::Vector<amrex::Real> pb_updateTime;
    amrex::Vector<amrex::Real> pb_amplitude;

  private:

    // Turbulent perturbation box array
    //amrex::Vector<amrex::BoxArray> turbPert_cba;
    //amrex::Vector<amrex::BoxArray> turbPert_zba;
};
#endif
