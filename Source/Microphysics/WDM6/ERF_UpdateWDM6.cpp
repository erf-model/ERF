#include "ERF_WDM6.H"
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <iomanip>

using namespace amrex;

void
WDM6::Copy_Micro_to_State(MultiFab& cons)
{
    // DEBUG: Track calls to Copy_Micro_to_State
    static int copy_to_state_count = 0;
    copy_to_state_count++;

    // Conservative update of all fields
    for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto& box3d = mfi.tilebox();
        auto states = cons.array(mfi);

        const auto& rho   = mic_fab_vars[MicVar_WDM6::rho]->array(mfi);
        const auto& theta = mic_fab_vars[MicVar_WDM6::theta]->array(mfi);
        const auto& qv    = mic_fab_vars[MicVar_WDM6::qv]->array(mfi);
        const auto& qc    = mic_fab_vars[MicVar_WDM6::qc]->array(mfi);
        const auto& qi    = mic_fab_vars[MicVar_WDM6::qi]->array(mfi);
        const auto& qr    = mic_fab_vars[MicVar_WDM6::qr]->array(mfi);
        const auto& qs    = mic_fab_vars[MicVar_WDM6::qs]->array(mfi);
        const auto& qg    = mic_fab_vars[MicVar_WDM6::qg]->array(mfi);
        const auto& nn    = mic_fab_vars[MicVar_WDM6::nn]->array(mfi);
        const auto& nc    = mic_fab_vars[MicVar_WDM6::nc]->array(mfi);
        const auto& nr    = mic_fab_vars[MicVar_WDM6::nr]->array(mfi);

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            states(i,j,k,RhoTheta_comp) = rho(i,j,k) * theta(i,j,k);

            states(i,j,k,RhoQ1_comp) = rho(i,j,k) * amrex::max(Real(0), qv(i,j,k));
            states(i,j,k,RhoQ2_comp) = rho(i,j,k) * amrex::max(Real(0), qc(i,j,k));
            states(i,j,k,RhoQ3_comp) = rho(i,j,k) * amrex::max(Real(0), qi(i,j,k));
            states(i,j,k,RhoQ4_comp) = rho(i,j,k) * amrex::max(Real(0), qr(i,j,k));
            states(i,j,k,RhoQ5_comp) = rho(i,j,k) * amrex::max(Real(0), qs(i,j,k));
            states(i,j,k,RhoQ6_comp) = rho(i,j,k) * amrex::max(Real(0), qg(i,j,k));

            // Number concentrations
            states(i,j,k,RhoQ7_comp) = rho(i,j,k) * amrex::max(Real(0), nc(i,j,k));
            states(i,j,k,RhoQ8_comp) = rho(i,j,k) * amrex::max(Real(0), nn(i,j,k));
            states(i,j,k,RhoQ9_comp) = rho(i,j,k) * amrex::max(Real(0), nr(i,j,k));
        });
    }

    // After first Copy_Micro_to_State, nn has been written to state (RhoQ8).
    // Now it's safe to clear the flag so subsequent Copy_State_to_Micro calls
    // will read nn from state like any other variable.
    if (m_nn_initialized) {
        amrex::Print() << "Copy_Micro_to_State call #" << copy_to_state_count
                      << ": nn written to state, clearing m_nn_initialized flag\n";
        m_nn_initialized = false;
    }

    // Diagnostic: Print comprehensive hydrometeor statistics every timestep
    // Global domain statistics (across all MPI ranks)
    static int update_call_count = 0;
    update_call_count++;

    {
        // Simple approach: use MultiFab min/max/sum methods
        Real min_qv = mic_fab_vars[MicVar_WDM6::qv]->min(0);
        Real max_qv = mic_fab_vars[MicVar_WDM6::qv]->max(0);
        Real sum_qv = mic_fab_vars[MicVar_WDM6::qv]->sum(0);

        Real min_qc = mic_fab_vars[MicVar_WDM6::qc]->min(0);
        Real max_qc = mic_fab_vars[MicVar_WDM6::qc]->max(0);
        Real sum_qc = mic_fab_vars[MicVar_WDM6::qc]->sum(0);

        Real min_qi = mic_fab_vars[MicVar_WDM6::qi]->min(0);
        Real max_qi = mic_fab_vars[MicVar_WDM6::qi]->max(0);
        Real sum_qi = mic_fab_vars[MicVar_WDM6::qi]->sum(0);

        Real min_qr = mic_fab_vars[MicVar_WDM6::qr]->min(0);
        Real max_qr = mic_fab_vars[MicVar_WDM6::qr]->max(0);
        Real sum_qr = mic_fab_vars[MicVar_WDM6::qr]->sum(0);

        Real min_qs = mic_fab_vars[MicVar_WDM6::qs]->min(0);
        Real max_qs = mic_fab_vars[MicVar_WDM6::qs]->max(0);
        Real sum_qs = mic_fab_vars[MicVar_WDM6::qs]->sum(0);

        Real min_qg = mic_fab_vars[MicVar_WDM6::qg]->min(0);
        Real max_qg = mic_fab_vars[MicVar_WDM6::qg]->max(0);
        Real sum_qg = mic_fab_vars[MicVar_WDM6::qg]->sum(0);

        Real min_nn = mic_fab_vars[MicVar_WDM6::nn]->min(0);
        Real max_nn = mic_fab_vars[MicVar_WDM6::nn]->max(0);
        Real sum_nn = mic_fab_vars[MicVar_WDM6::nn]->sum(0);

        Real min_nc = mic_fab_vars[MicVar_WDM6::nc]->min(0);
        Real max_nc = mic_fab_vars[MicVar_WDM6::nc]->max(0);
        Real sum_nc = mic_fab_vars[MicVar_WDM6::nc]->sum(0);

        Real min_nr = mic_fab_vars[MicVar_WDM6::nr]->min(0);
        Real max_nr = mic_fab_vars[MicVar_WDM6::nr]->max(0);
        Real sum_nr = mic_fab_vars[MicVar_WDM6::nr]->sum(0);

        // Compute Qp = qr + qs + qg by creating temporary MultiFab
        MultiFab qp_mf(mic_fab_vars[MicVar_WDM6::qr]->boxArray(),
                       mic_fab_vars[MicVar_WDM6::qr]->DistributionMap(), 1, 0);
        MultiFab::LinComb(qp_mf, 1.0, *mic_fab_vars[MicVar_WDM6::qr], 0,
                                  1.0, *mic_fab_vars[MicVar_WDM6::qs], 0, 0, 1, 0);
        MultiFab::Add(qp_mf, *mic_fab_vars[MicVar_WDM6::qg], 0, 0, 1, 0);
        Real min_qp = qp_mf.min(0);
        Real max_qp = qp_mf.max(0);
        Real sum_qp = qp_mf.sum(0);

        // Get total cell count
        Long total_cells = mic_fab_vars[MicVar_WDM6::qv]->boxArray().numPts();
        ParallelDescriptor::ReduceLongSum(total_cells);

        // MPI reductions
        ParallelDescriptor::ReduceRealMin(min_qv);
        ParallelDescriptor::ReduceRealMax(max_qv);
        ParallelDescriptor::ReduceRealSum(sum_qv);
        ParallelDescriptor::ReduceRealMin(min_qc);
        ParallelDescriptor::ReduceRealMax(max_qc);
        ParallelDescriptor::ReduceRealSum(sum_qc);
        ParallelDescriptor::ReduceRealMin(min_qi);
        ParallelDescriptor::ReduceRealMax(max_qi);
        ParallelDescriptor::ReduceRealSum(sum_qi);
        ParallelDescriptor::ReduceRealMin(min_qr);
        ParallelDescriptor::ReduceRealMax(max_qr);
        ParallelDescriptor::ReduceRealSum(sum_qr);
        ParallelDescriptor::ReduceRealMin(min_qs);
        ParallelDescriptor::ReduceRealMax(max_qs);
        ParallelDescriptor::ReduceRealSum(sum_qs);
        ParallelDescriptor::ReduceRealMin(min_qg);
        ParallelDescriptor::ReduceRealMax(max_qg);
        ParallelDescriptor::ReduceRealSum(sum_qg);
        ParallelDescriptor::ReduceRealMin(min_qp);
        ParallelDescriptor::ReduceRealMax(max_qp);
        ParallelDescriptor::ReduceRealSum(sum_qp);
        ParallelDescriptor::ReduceRealMin(min_nn);
        ParallelDescriptor::ReduceRealMax(max_nn);
        ParallelDescriptor::ReduceRealSum(sum_nn);
        ParallelDescriptor::ReduceRealMin(min_nc);
        ParallelDescriptor::ReduceRealMax(max_nc);
        ParallelDescriptor::ReduceRealSum(sum_nc);
        ParallelDescriptor::ReduceRealMin(min_nr);
        ParallelDescriptor::ReduceRealMax(max_nr);
        ParallelDescriptor::ReduceRealSum(sum_nr);

        // Compute means
        Real mean_qv = sum_qv / total_cells;
        Real mean_qc = sum_qc / total_cells;
        Real mean_qi = sum_qi / total_cells;
        Real mean_qr = sum_qr / total_cells;
        Real mean_qs = sum_qs / total_cells;
        Real mean_qg = sum_qg / total_cells;
        Real mean_qp = sum_qp / total_cells;
        Real mean_nn = sum_nn / total_cells;
        Real mean_nc = sum_nc / total_cells;
        Real mean_nr = sum_nr / total_cells;

        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "\n=== ERF WDM6 Statistics: Timestep " << update_call_count << " ===\n";
            amrex::Print() << std::scientific << std::setprecision(4);
            amrex::Print() << "qv (kg/kg):  min=" << min_qv*1e3
                          << "  max=" << max_qv*1e3
                          << "  mean=" << mean_qv*1e3 << " g/kg\n";
            amrex::Print() << "qc (kg/kg):  min=" << min_qc*1e3
                          << "  max=" << max_qc*1e3
                          << "  mean=" << mean_qc*1e3 << " g/kg\n";
            amrex::Print() << "qi (kg/kg):  min=" << min_qi*1e3
                          << "  max=" << max_qi*1e3
                          << "  mean=" << mean_qi*1e3 << " g/kg\n";
            amrex::Print() << "qr (kg/kg):  min=" << min_qr*1e3
                          << "  max=" << max_qr*1e3
                          << "  mean=" << mean_qr*1e3 << " g/kg\n";
            amrex::Print() << "qs (kg/kg):  min=" << min_qs*1e3
                          << "  max=" << max_qs*1e3
                          << "  mean=" << mean_qs*1e3 << " g/kg\n";
            amrex::Print() << "qg (kg/kg):  min=" << min_qg*1e3
                          << "  max=" << max_qg*1e3
                          << "  mean=" << mean_qg*1e3 << " g/kg\n";
            amrex::Print() << "Qp (kg/kg):  min=" << min_qp*1e3
                          << "  max=" << max_qp*1e3
                          << "  mean=" << mean_qp*1e3 << " g/kg  (qr+qs+qg)\n";
            amrex::Print() << "nn (#/kg):   min=" << min_nn
                          << "  max=" << max_nn
                          << "  mean=" << mean_nn << "\n";
            amrex::Print() << "nc (#/kg):   min=" << min_nc
                          << "  max=" << max_nc
                          << "  mean=" << mean_nc << "\n";
            amrex::Print() << "nr (#/kg):   min=" << min_nr
                          << "  max=" << max_nr
                          << "  mean=" << mean_nr << "\n";
            amrex::Print() << "========================================\n\n";
        }
    }

    cons.FillBoundary(m_geom.periodicity());
}
