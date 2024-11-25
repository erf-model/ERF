#include <ERF.H>
#include <ERF_ParFunctions.H>

using namespace amrex;

void ERF::advance_microphysics (int lev,
                                MultiFab& cons,
                                const Real& dt_advance,
                                const int& iteration,
                                const Real& time )
{
    if (solverChoice.moisture_type != MoistureType::None) {

        // Test the PS interface.
        // These type of ops need to be done on GPU
        // but outside the micro class. Note, we can't
        // operate through the class pointers on GPU.
        int size(0);
        int* ind_map = micro->NonPrecip_Index_Map(size);
        Print() << "CHECK: " << size << ' '
                << (ind_map==nullptr) << "\n";
        for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const auto& tbx  = mfi.tilebox();
            const auto& data = cons.const_array(mfi);
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real partial_sum = comp_partial_sum (i, j, k,
                                                     size, ind_map, data);
                if (i==0 && j==0) {
                    printf("Partial Micro Sum: %i %e\n",k,partial_sum);
                }
            });
        }


        micro->Update_Micro_Vars_Lev(lev, cons);
        micro->Advance(lev, dt_advance, iteration, time, solverChoice, vars_new, z_phys_nd);
        micro->Update_State_Vars_Lev(lev, cons);
    }
}
