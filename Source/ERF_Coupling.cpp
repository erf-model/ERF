#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_MicrophysicsUtils.H>

#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

amrex::Real
ERF::EvolveOneStep (amrex::Real /*time*/, amrex::Real /*dt_request*/)
{
    amrex::Real cur_time = t_new[0];
    const int step = istep[0];

    if (start_time + cur_time >= stop_time) {
        return amrex::Real(0.0);
    }

    ComputeDt(step);

    int iteration = 1;
    timeStep(0, cur_time, iteration);
    cur_time += dt[0];

    post_timestep(step, cur_time, dt[0]);

    // ****************************************************************************************
    // Write plotfiles at intermediate times
    // ****************************************************************************************

    if (writeNow(cur_time, step+1, m_plot3d_int_1, m_plot3d_per_1, dt[0], last_plot3d_file_time_1)) {
        last_plot3d_file_step_1 = step+1;
        Write3DPlotFile(1,plotfile3d_type_1,plot3d_var_names_1);
        for (int lev = 0; lev <= finest_level; ++lev) {lsm.Plot(lev, step+1);}
        if (m_plot3d_per_1 > amrex::Real(0.0)) {last_plot3d_file_time_1 += m_plot3d_per_1;}
    }
    if (writeNow(cur_time, step+1, m_plot3d_int_2, m_plot3d_per_2, dt[0], last_plot3d_file_time_2)) {
        last_plot3d_file_step_2 = step+1;
        Write3DPlotFile(2,plotfile3d_type_2,plot3d_var_names_2);
        for (int lev = 0; lev <= finest_level; ++lev) {lsm.Plot(lev, step+1);}
        if (m_plot3d_per_2 > amrex::Real(0.0)) {last_plot3d_file_time_2 += m_plot3d_per_2;}
    }
    if (writeNow(cur_time, step+1, m_plot2d_int_1, m_plot2d_per_1, dt[0], last_plot2d_file_time_1)) {
        last_plot2d_file_step_1 = step+1;
        Write2DPlotFile(1,plotfile2d_type_1,plot2d_var_names_1);
        if (m_plot2d_per_1 > amrex::Real(0.0)) {last_plot2d_file_time_1 += m_plot2d_per_1;}
    }
    if (writeNow(cur_time, step+1, m_plot2d_int_2, m_plot2d_per_2, dt[0], last_plot2d_file_time_2)) {
        last_plot2d_file_step_2 = step+1;
        Write2DPlotFile(2,plotfile2d_type_2,plot2d_var_names_2);
        if (m_plot2d_per_2 > amrex::Real(0.0)) {last_plot2d_file_time_2 += m_plot2d_per_2;}
    }
    for (int i = 0; i < m_subvol_int.size(); i++) {
        if (writeNow(cur_time, step+1, m_subvol_int[i], m_subvol_per[i], dt[0], last_subvol_time[i])) {
            last_subvol_step[i] = step+1;
            WriteSubvolume(i,subvol3d_var_names);
            if (m_subvol_per[i] > amrex::Real(0.0)) {last_subvol_time[i] += m_subvol_per[i];}
        }
    }
    if (writeNow(cur_time, step+1, m_check_int, m_check_per, dt[0], last_check_file_time)) {
        last_check_file_step = step+1;
        WriteCheckpointFile();
        if (m_check_per > amrex::Real(0.0)) {last_check_file_time += m_check_per;}
    }

    WriteAtFinalTime();

    return dt[0];
}

/*
  Coupling reference context (implementation-side):

  1) Legacy state-passing contract:
     Warner et al. (2010), COAWST, Fig. 5 / Block B.
     Atmosphere -> ocean fields are passed as states
     (Uwind, Vwind, Patm, RH, Tair, cloud, rain, SWrad, LWrad), while
     ocean -> atmosphere provides SST.

  2) Future direct flux-passing roadmap:
     COAWST's ATM2OCN_FLUXES pathway (documented in COAWST manuals/workshops
     and exercised in Zambon et al., 2014) motivates a direct
     stress/heat/moisture flux mode.
     COAWST code anchors:
       - Master/mct_roms_wrf.h
       - ROMS/Nonlinear/atm2ocn_flux.F
       - ROMS/Nonlinear/bulk_flux.F
     This file currently implements the legacy state-passing test path only.
*/

void
ERF::PackAtmosphericStates (amrex::Vector<amrex::MultiFab*>& states,
                            amrex::Real /*time*/)
{
    using namespace amrex;

    // Contract slot indices (mirrors ERFRemoraCouplingContract.H; repeated here
    // to avoid a driver→submodule header dependency).
    constexpr int iUwind = 0, iVwind = 1, iPatm = 2, iRH = 3, iTair = 4;
    constexpr int iCloud = 5, iRain  = 6, iSWrad = 7, iLWrad = 8;

    const int lev   = 0;

    auto& cons = vars_new[lev][Vars::cons];
    auto& xvel = vars_new[lev][Vars::xvel]; // XFace
    auto& yvel = vars_new[lev][Vars::yvel]; // YFace

    const bool has_moisture  = (solverChoice.moisture_type != MoistureType::None);
    const bool has_radiation = (!rad_fluxes.empty() && rad_fluxes[lev] != nullptr);

    amrex::ignore_unused(has_moisture, has_radiation);

    const auto& ba = cons.boxArray();
    const auto& dm = cons.DistributionMap();
    const auto& ba2d_lev = ba2d[lev];

    // --- Uwind + Vwind: use AMReX's average_face_to_cellcenter which correctly
    // handles tile boundaries via growntilebox(1) internally. ---
    if ((iUwind < static_cast<int>(states.size()) && states[iUwind] != nullptr) ||
        (iVwind < static_cast<int>(states.size()) && states[iVwind] != nullptr)) {

        auto& zvel = vars_new[lev][Vars::zvel];
        MultiFab cc_vel(ba, dm, AMREX_SPACEDIM, 0);
        amrex::average_face_to_cellcenter(cc_vel, 0,
            Array<const MultiFab*, AMREX_SPACEDIM>{&xvel, &yvel, &zvel});

        // Collapse to 2D slab
        MultiFab uv_slab(ba2d_lev, dm, 2, 0);  // comp0=u, comp1=v
        uv_slab.ParallelCopy(cc_vel, 0, 0, 2);

        if (iUwind < static_cast<int>(states.size()) && states[iUwind] != nullptr) {
            MultiFab u_alias(uv_slab, amrex::make_alias, 0, 1); // alias u component
            IntVect ratio = ba2d_lev.minimalBox().length()
                          / states[iUwind]->boxArray().minimalBox().length();
            amrex::average_down(u_alias, *states[iUwind], 0, 1, ratio);
        }

        if (iVwind < static_cast<int>(states.size()) && states[iVwind] != nullptr) {
            MultiFab v_alias(uv_slab, amrex::make_alias, 1, 1); // alias v component
            IntVect ratio = ba2d_lev.minimalBox().length()
                          / states[iVwind]->boxArray().minimalBox().length();
            amrex::average_down(v_alias, *states[iVwind], 0, 1, ratio);
        }
    }

    // --- Patm: getPgivenRTh(RhoTheta, qv) at k=0 ---
    if (iPatm < static_cast<int>(states.size()) && states[iPatm] != nullptr) {
        MultiFab tmp(ba2d_lev, dm, 1, 0);
        for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            auto const& c = cons.const_array(mfi);
            auto         t = tmp.array(mfi);
            if (has_moisture) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    const Real qv = c(i,j,k,RhoQ1_comp) / c(i,j,k,Rho_comp);
                    t(i,j,k) = getPgivenRTh(c(i,j,k,RhoTheta_comp), qv);
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = getPgivenRTh(c(i,j,k,RhoTheta_comp));
                });
            }
        }
        IntVect ratio = ba2d_lev.minimalBox().length() / states[iPatm]->boxArray().minimalBox().length();
        amrex::average_down(tmp, *states[iPatm], 0, 1, ratio);
    }

    // --- Tair: getTgivenRandRTh(rho, RhoTheta, qv) at k=0 [K] ---
    if (iTair < static_cast<int>(states.size()) && states[iTair] != nullptr) {
        MultiFab tmp(ba2d_lev, dm, 1, 0);
        for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            auto const& c = cons.const_array(mfi);
            auto         t = tmp.array(mfi);
            if (has_moisture) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    const Real qv = c(i,j,k,RhoQ1_comp) / c(i,j,k,Rho_comp);
                    t(i,j,k) = getTgivenRandRTh(c(i,j,k,Rho_comp), c(i,j,k,RhoTheta_comp), qv);
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = getTgivenRandRTh(c(i,j,k,Rho_comp), c(i,j,k,RhoTheta_comp));
                });
            }
        }
        IntVect ratio = ba2d_lev.minimalBox().length() / states[iTair]->boxArray().minimalBox().length();
        amrex::average_down(tmp, *states[iTair], 0, 1, ratio);
    }

    // --- Humidity lane: export relative humidity [0-1] for REMORA bulk fluxes ---
    if (has_moisture) {
#if 0
        if (iRH < static_cast<int>(states.size()) && states[iRH] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                auto const& c = cons.const_array(mfi);
                auto         t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    const Real qv = c(i,j,k,RhoQ1_comp) / c(i,j,k,Rho_comp);
                    const Real p_pa = getPgivenRTh(c(i,j,k,RhoTheta_comp), qv);
                    const Real temp = getTgivenRandRTh(c(i,j,k,Rho_comp), c(i,j,k,RhoTheta_comp), qv);
                    Real qsat = Real(0.0);
                    erf_qsatw(temp, p_pa * Real(0.01), qsat);
                    t(i,j,k) = amrex::max(Real(0.0), amrex::min(Real(1.0),
                                                                 qv / amrex::max(qsat, Real(1.0e-12))));
                });
            }
            IntVect ratio = ba2d_lev.minimalBox().length() / states[iRH]->boxArray().minimalBox().length();
            amrex::average_down(tmp, *states[iRH], 0, 1, ratio);
        }
#else
        amrex::ignore_unused(iRH);
        // Avoiding unexercised moist assumptions: retain the driver-prefilled
        // humidity lane until live ERF->REMORA Qair semantics are made explicit.
#endif
        if (iCloud < static_cast<int>(states.size()) && states[iCloud] != nullptr) {
            const int qc_idx = solverChoice.moisture_indices.qc;
            const int qi_idx = solverChoice.moisture_indices.qi;
            if (qc_idx != -1 || qi_idx != -1) {
                MultiFab tmp(ba2d_lev, dm, 1, 0);
                const Real cf = amrex::max(Real(0.0), amrex::min(Real(1.0), cloud_fraction(0.0)));
                amrex::ignore_unused(cf); // keep diagnostic computation active for consistency with ERF scalar stats
                for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box bx = mfi.tilebox();
                    auto const& c = cons.const_array(mfi);
                    auto t = tmp.array(mfi);
                    const int klo = ba.minimalBox().smallEnd(2);
                    const int khi = ba.minimalBox().bigEnd(2);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        int cloudy = 0;
                        for (int kk = klo; kk <= khi; ++kk) {
                            const Real rho = c(i,j,kk,Rho_comp);
                            const Real qc = (qc_idx != -1) ? c(i,j,kk,qc_idx) / rho : Real(0.0);
                            const Real qi = (qi_idx != -1) ? c(i,j,kk,qi_idx) / rho : Real(0.0);
                            if (qc + qi > Real(0.0)) { cloudy = 1; break; }
                        }
                        t(i,j,k) = static_cast<Real>(cloudy);
                    });
                }
                IntVect ratio = ba2d_lev.minimalBox().length() / states[iCloud]->boxArray().minimalBox().length();
                amrex::average_down(tmp, *states[iCloud], 0, 1, ratio);
            }
        }
#if 0
        if (iRain < static_cast<int>(states.size()) && states[iRain] != nullptr) {
            int qr_idx = solverChoice.moisture_indices.qr;
            if (qr_idx != -1) {
                MultiFab tmp(ba2d_lev, dm, 1, 0);
                for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box bx = mfi.tilebox();
                    auto const& c = cons.const_array(mfi);
                    auto         t = tmp.array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        t(i,j,k) = c(i,j,k,qr_idx) / c(i,j,k,Rho_comp);
                    });
                }
                IntVect ratio = ba2d_lev.minimalBox().length() / states[iRain]->boxArray().minimalBox().length();
                amrex::average_down(tmp, *states[iRain], 0, 1, ratio);
            }
        }
#else
        amrex::ignore_unused(iRain);
        // Avoiding unexercised moist assumptions: retain the driver-prefilled
        // rain lane until live ERF->REMORA rain semantics are made explicit.
#endif
    }
    // No moisture: leave RH/Cloud/Rain slabs at their driver-pre-filled values.

    // --- Radiation: sw_flux_dn (comp=1) and lw_flux_dn (comp=3) from rad_fluxes ---
    // When absent, leave slabs at their driver-pre-filled values.
    if (has_radiation) {
        if (iSWrad < static_cast<int>(states.size()) && states[iSWrad] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            tmp.ParallelCopy(*rad_fluxes[lev], 1, 0, 1);
            IntVect ratio = ba2d_lev.minimalBox().length() / states[iSWrad]->boxArray().minimalBox().length();
            amrex::average_down(tmp, *states[iSWrad], 0, 1, ratio);
        }
        if (iLWrad < static_cast<int>(states.size()) && states[iLWrad] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            tmp.ParallelCopy(*rad_fluxes[lev], 3, 0, 1);
            IntVect ratio = ba2d_lev.minimalBox().length() / states[iLWrad]->boxArray().minimalBox().length();
            amrex::average_down(tmp, *states[iLWrad], 0, 1, ratio);
        }
    }
}

void
ERF::ApplyOceanSurfaceState (const amrex::Vector<amrex::MultiFab*>& state,
                             amrex::Real time)
{
    if (solverChoice.lsm_type != LandSurfaceType::OceanSurf) {
        return;
    }

    if (!state.empty() && state[0] != nullptr && lsm.Get_Data_Ptr(0, 0) != nullptr) {
        auto* dst = lsm.Get_Data_Ptr(0, 0);
        const auto lsm_geom = lsm.Get_Lsm_Geom(0);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            lsm_geom.isPeriodic(0) == Geom(0).isPeriodic(0) &&
            lsm_geom.isPeriodic(1) == Geom(0).isPeriodic(1) &&
            lsm_geom.isPeriodic(2) == Geom(0).isPeriodic(2),
            "OceanSurf t_surf geometry lost ERF periodic flags.");
        const int dst_k = dst->boxArray().minimalBox().smallEnd(2);
        const int src_k = state[0]->boxArray().minimalBox().bigEnd(2);
        for (amrex::MFIter mfi(*dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx = amrex::makeSlab(mfi.validbox(), 2, dst_k);
            auto dst_arr = dst->array(mfi);
            auto src_arr = state[0]->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
                dst_arr(i,j,dst_k) = src_arr(i,j,src_k);
            });
        }
        dst->FillBoundary(lsm_geom.periodicity());
        amrex::Gpu::streamSynchronize();

        const amrex::Real src_min = state[0]->min(0);
        const amrex::Real src_max = state[0]->max(0);
        const amrex::Real dst_min = dst->min(0);
        const amrex::Real dst_max = dst->max(0);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "OceanSurf apply at t=" << time
                           << " s from SST slab: source min/max = "
                           << src_min << " / " << src_max
                           << " K, cache min/max = "
                           << dst_min << " / " << dst_max
                           << " K" << std::endl;
        }

        if (dst_min < amrex::Real(260.0) || dst_max > amrex::Real(320.0)) {
            amrex::Warning("OceanSurf t_surf is outside the expected [260, 320] K range");
        }
    }
}
