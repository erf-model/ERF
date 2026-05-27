#include <ERF.H>
#include <ERF_Constants.H>

using namespace amrex;

void ERF::advance_microphysics (int lev,
                                MultiFab& cons,
                                const Real& dt_advance,
                                const int& iteration,
                                const Real& time )
{
    if (solverChoice.moisture_type != MoistureType::None) {
        micro->Set_RealWidth(lev, real_width);
        cons.FillBoundary(geom[lev].periodicity());

        // Get number of moisture species for the microphysics scheme.
        int n_qstate = micro->Get_Qstate_Moist_Size();

        // Zero condensed moisture in boundary cells when using real BCs.
        // Only zero if we have hydrometeors (n_qstate > 1 means QV + at least one hydrometeor)
        if (real_width > 0 && lev == 0 && n_qstate > 1) {
            amrex::Print() << "DEBUG: Zeroing boundary hydrometeors at time=" << time
                          << " n_qstate=" << n_qstate << " real_width=" << real_width << "\n";
            const auto& domain = geom[lev].Domain();
            int i_lo = domain.smallEnd(0);
            int i_hi = domain.bigEnd(0);
            int j_lo = domain.smallEnd(1);
            int j_hi = domain.bigEnd(1);

            for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                auto cons_arr = cons.array(mfi);
                Box tbx = mfi.tilebox();

                // Zero boundary cells on domain faces
                if (tbx.smallEnd(0) == i_lo) {
                    Box bx_xlo(tbx); bx_xlo.setBig(0, i_lo + real_width - 1);
                    amrex::Print() << "DEBUG: Zeroing xlo boundary, box=" << bx_xlo << "\n";
                    ParallelFor(bx_xlo, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Zero only hydrometeors used by the current microphysics scheme
                        if (n_qstate >= 2) cons_arr(i,j,k,RhoQ2_comp) = zero;
                        if (n_qstate >= 3) cons_arr(i,j,k,RhoQ3_comp) = zero;
                        if (n_qstate >= 4) cons_arr(i,j,k,RhoQ4_comp) = zero;
                        if (n_qstate >= 5) cons_arr(i,j,k,RhoQ5_comp) = zero;
                        if (n_qstate >= 6) cons_arr(i,j,k,RhoQ6_comp) = zero;
                    });
                }
                if (tbx.bigEnd(0) == i_hi) {
                    Box bx_xhi(tbx); bx_xhi.setSmall(0, i_hi - real_width + 1);
                    amrex::Print() << "DEBUG: Zeroing xhi boundary, box=" << bx_xhi << "\n";
                    ParallelFor(bx_xhi, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        if (n_qstate >= 2) cons_arr(i,j,k,RhoQ2_comp) = zero;
                        if (n_qstate >= 3) cons_arr(i,j,k,RhoQ3_comp) = zero;
                        if (n_qstate >= 4) cons_arr(i,j,k,RhoQ4_comp) = zero;
                        if (n_qstate >= 5) cons_arr(i,j,k,RhoQ5_comp) = zero;
                        if (n_qstate >= 6) cons_arr(i,j,k,RhoQ6_comp) = zero;
                    });
                }
                if (tbx.smallEnd(1) == j_lo) {
                    Box bx_ylo(tbx); bx_ylo.setBig(1, j_lo + real_width - 1);
                    amrex::Print() << "DEBUG: Zeroing ylo boundary, box=" << bx_ylo << "\n";
                    ParallelFor(bx_ylo, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        if (n_qstate >= 2) cons_arr(i,j,k,RhoQ2_comp) = zero;
                        if (n_qstate >= 3) cons_arr(i,j,k,RhoQ3_comp) = zero;
                        if (n_qstate >= 4) cons_arr(i,j,k,RhoQ4_comp) = zero;
                        if (n_qstate >= 5) cons_arr(i,j,k,RhoQ5_comp) = zero;
                        if (n_qstate >= 6) cons_arr(i,j,k,RhoQ6_comp) = zero;
                    });
                }
                if (tbx.bigEnd(1) == j_hi) {
                    Box bx_yhi(tbx); bx_yhi.setSmall(1, j_hi - real_width + 1);
                    amrex::Print() << "DEBUG: Zeroing yhi boundary, box=" << bx_yhi << "\n";
                    ParallelFor(bx_yhi, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        if (n_qstate >= 2) cons_arr(i,j,k,RhoQ2_comp) = zero;
                        if (n_qstate >= 3) cons_arr(i,j,k,RhoQ3_comp) = zero;
                        if (n_qstate >= 4) cons_arr(i,j,k,RhoQ4_comp) = zero;
                        if (n_qstate >= 5) cons_arr(i,j,k,RhoQ5_comp) = zero;
                        if (n_qstate >= 6) cons_arr(i,j,k,RhoQ6_comp) = zero;
                    });
                }
            }
        }

        // Check boundary values before microphysics
        if (real_width > 0 && lev == 0 && n_qstate > 1) {
            const auto& domain = geom[lev].Domain();
            int i_lo = domain.smallEnd(0);
            Real max_q2 = cons.max(RhoQ2_comp);
            Real max_q3 = cons.max(RhoQ3_comp);
            amrex::Print() << "DEBUG: Before micro->Advance at time=" << time
                          << " max(RhoQ2)=" << max_q2 << " max(RhoQ3)=" << max_q3 << "\n";
            // Sample boundary cell
            for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
                auto cons_arr = cons.array(mfi);
                Box tbx = mfi.tilebox();
                if (tbx.smallEnd(0) == i_lo) {
                    amrex::Print() << "DEBUG: Sample xlo boundary at i=" << i_lo
                                  << " RhoQ2=" << cons_arr(i_lo,tbx.smallEnd(1),tbx.smallEnd(2),RhoQ2_comp)
                                  << " RhoQ3=" << cons_arr(i_lo,tbx.smallEnd(1),tbx.smallEnd(2),RhoQ3_comp) << "\n";
                    break;
                }
            }
        }

        micro->Update_Micro_Vars_Lev(lev, cons);
        micro->Advance(lev, dt_advance, iteration, time, solverChoice, vars_new, z_phys_nd, phys_bc_type);
        micro->Update_State_Vars_Lev(lev, cons, *z_phys_nd[lev]);

        // Check boundary values after microphysics
        if (real_width > 0 && lev == 0 && n_qstate > 1) {
            Real max_q2 = cons.max(RhoQ2_comp);
            Real max_q3 = cons.max(RhoQ3_comp);
            amrex::Print() << "DEBUG: After micro->Advance at time=" << time
                          << " max(RhoQ2)=" << max_q2 << " max(RhoQ3)=" << max_q3 << "\n";
        }
    }
}
