#ifndef _TURB_PERT_STRUCT_H_
#define _TURB_PERT_STRUCT_H_

/**
 * Container holding quantities related to turbulent perturbation parameters
 */

//#define DEBUG_PERTBOX_MSG

struct TurbulentPerturbation {
  public:

    // Initializing Perturbation Region
    void init_PerturbationRegion (int lev)
    {
        int pertBox_offset = 0;
        amrex::Box turbPert_bx(amrex::IntVect(0+pertBox_offset,0,0), amrex::IntVect(7+pertBox_offset,31,31));
        amrex::BoxArray tmp_ba(turbPert_bx);
        turbPert_ba.push_back(tmp_ba);
        turbPert_ba[lev].maxSize(amrex::IntVect(8,8,8));
        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "  [TurbPert::init_PerturbationRegion] Subdividing into smaller boxes: "  << turbPert_ba[lev].size() << "\n"
                << "  turbPert_ba[" << lev << "] constains: " << turbPert_ba[lev];
        #endif
        amrex:: Print() << "Turbulence perturbation region initialized\n";
    }

    // Calculate perturbation update time
    void calc_TurbPert_updateTime (int lev)
    {

    }

    // Calculate perturbation amplitude
    void calc_TurbPert_amplitude (int lev, const amrex::Box bx, const auto &cons_pert_arr)
    {
        #ifdef DEBUG_PERTBOX_MSG
        amrex::Print() << "  [TurbPert::calc_TurbPert_amplitude] Inside of structure function\n";
        #endif
        for (int boxIdx = 0; boxIdx < turbPert_ba[lev].size(); boxIdx++) {
            if (bx.contains(turbPert_ba[lev][boxIdx])) {
                #ifdef DEBUG_PERTBOX_MSG
                amrex::Print() << "bx: " << bx << " -- contains perturbation box #" << boxIdx << ": " << turbPert_ba[lev][boxIdx] << "\n";
                #endif
                ParallelFor(turbPert_ba[lev][boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    //Real tmp_Src = (i + j*bx.size()[0] + k*bx.size()[0]*bx.size()[1]); // continuous count
                    amrex::Real tmp_Src = (amrex::Real) (boxIdx+1.0); // distinguish each box

                    // Adding temperature source onto RHS
                    cons_pert_arr(i, j, k, RhoTheta_comp) = tmp_Src;
                });
            }
        } // boxIdx
    }

  private:

    // Turbulent perturbation box array
    amrex::Vector<amrex::BoxArray> turbPert_ba;
};
#endif
