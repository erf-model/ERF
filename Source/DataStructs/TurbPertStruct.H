#ifndef _TURB_PERT_STRUCT_H_
#define _TURB_PERT_STRUCT_H_

/**
 * Container holding quantities related to turbulent perturbation parameters
 */

#include <AMReX_MultiFab.H>

#define DEBUG_PERTBOX_MSG

#define MANUAL_DISTRIBUTION_MAP_ASSIGNMENT

//#define INDEX_PERTURB
#define RANDOM_PERTURB

struct TurbulentPerturbation {
  public:

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
        amrex::Print() << "  [TurbPert::init_PerturbationRegion] C subdividing into smaller boxes: "
                       << turbPert_cba[lev].size() << "\n"
                       << "  turbPert_cba[" << lev << "] constains: " << turbPert_cba[lev];
        #endif

        // TODO: Once mf is accessing data (ie. u, and v velocitices), come back and test if this fails... for science
        #ifdef MANUAL_DISTRIBUTION_MAP_ASSIGNMENT
        // N^2 search to build local distribution map
        const amrex::Vector<int> d_pmap = data_dm.ProcessorMap();
        amrex::Vector<int> p_pmap(turbPert_cba[lev].size()); // Initialization of p_pmap with npbx elements

        for (int pbx(0); pbx < turbPert_cba[lev].size(); ++pbx) {
            bool found = false;
            for (int dbx(0); dbx < (int) data_ba.size(); ++dbx) {
                amrex::Box m_dbx = data_ba[dbx];
                if (m_dbx.contains(turbPert_cba[lev][pbx])) {
                    found = true;
                    p_pmap[pbx] = d_pmap[dbx];
                    continue;
                } // found
            } // data box
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(found, "Perturbation box not contained in a valid box!");
        } // pert box

        amrex::DistributionMapping dm(p_pmap);
        #else
        amrex::DistributionMapping dm;
        dm.define(turbPert_cba[lev],2);
        #endif

        // Creating multiFab for conserved quantity
        // XXX: TESTING NOW
        int ncomp = PertBox::NumTypes;
        int ngrow = 0;
        amrex::MultiFab mf_turbPert_c(turbPert_cba[lev], dm, ncomp, ngrow); // XXX This works, but mf_turbPert_c is only seen inside this function
        //mf_turbPert_c.define(turbPert_cba[lev], dm, ncomp, ngrow); // XXX This gives a "pure virtual method called" at runtime

        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "--------------------------DEBUGGING--------------------------\n"
                       << mf_turbPert_c.boxArray() << "\n" << mf_turbPert_c.DistributionMap() << "\n";
        #endif

        // Creating structure box array for z velocity
        amrex::BoxArray tmp_zba = amrex::convert(tmp_cba, amrex::IntVect(0,0,1)); // Convert for Z
        turbPert_zba.push_back(tmp_zba);
        turbPert_zba[lev].maxSize(boxSize);
        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "  [TurbPert::init_PerturbationRegion] Z subdividing into smaller boxes: "
                       << turbPert_zba[lev].size() << "\n"
                       << "  turbPert_zba[" << lev << "] constains: " << turbPert_zba[lev];
        #endif

        amrex:: Print() << "Turbulence perturbation region initialized\n";
    }

    // TODO: Write in function to call specific box to calculate new amplitude if update time is larger than the elapsed time
    // Calculate perturbation update time
    void calc_TurbPert_updateTime (const int lev)
    {

    }

    // TODO: Access velocity data to calculate perturbation based on simulation values
    // Calculate perturbation amplitude
    void calc_TurbPert_amplitude (const int lev,
                                  const amrex::Box& bx,
                                  const amrex::Vector<amrex::BoxArray>& ba,
                                  const int comp,
                                  const auto &src_arr)
    {
        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "  [TurbPert::calc_TurbPert_amplitude] Inside of structure function\n";
        #endif

        for (int boxIdx = 0; boxIdx < (int) ba[lev].size(); boxIdx++) {
            if (bx.contains(ba[lev][boxIdx])) {
                #ifdef DEBUG_PERTBOX_MSG
                amrex::Print() << "  bx: " << bx << " -- contains perturbation box #" << boxIdx << ": " << ba[lev][boxIdx] << "\n";
                #endif

                #ifdef INDEX_PERTURB
                // This exists for debugging purpose
                ParallelFor(ba[lev][boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    //amrex::Real tmp_Src = (i + j*bx.size()[0] + k*bx.size()[0]*bx.size()[1]); // continuous count
                    amrex::Real tmp_Src = (amrex::Real) (boxIdx+1.0); // distinguish each box
                    src_arr(i, j, k, comp) = tmp_Src;
                });
                #endif

                #ifdef RANDOM_PERTURB
                ParallelForRNG(ba[lev][boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
                    amrex::Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
                    amrex::Real tmp_Src = (amrex::Real) (boxIdx+1.0) * (rand_double*2.0 - 1.0); // Automatically rescales between -1.0 and 1.0
                    src_arr(i,j,k,comp) = tmp_Src;
                });
                #endif
            } // if
        } // boxIdx
    }

    // Turbulent perturbation box array
    amrex::Vector<amrex::BoxArray> turbPert_cba;
    amrex::Vector<amrex::BoxArray> turbPert_zba;

    // Creating MultiFab for Data
    amrex::MultiFab mf_turbPert_c;

  private:

    // Turbulent perturbation box array
    //amrex::Vector<amrex::BoxArray> turbPert_cba;
    //amrex::Vector<amrex::BoxArray> turbPert_zba;
};
#endif
