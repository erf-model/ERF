#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_MicrophysicsUtils.H>
#include <Diagnostics/ERF_SurfaceFluxDiagnostics.H>

#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

namespace
{
void
PrintFluxLaneStats (const char* label, const amrex::MultiFab& mf, int comp = 0)
{
    amrex::Print() << "ERF flux-pack " << label
                   << ": min=" << mf.min(comp)
                   << " max=" << mf.max(comp) << "\n";
}

void
AverageDownThenRemap (const amrex::MultiFab& src,
                      amrex::MultiFab& dst)
{
    using namespace amrex;

    AMREX_ALWAYS_ASSERT(src.boxArray().ixType() == dst.boxArray().ixType());

    const Box src_cells = enclosedCells(src.boxArray().minimalBox());
    const Box dst_cells = enclosedCells(dst.boxArray().minimalBox());
    const IntVect src_len = src_cells.length();
    const IntVect dst_len = dst_cells.length();

    const bool same_layout =
        (src.boxArray() == dst.boxArray()) &&
        (src.DistributionMap() == dst.DistributionMap());
    if (same_layout) {
        dst.ParallelCopy(src, 0, 0, dst.nComp());
        return;
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_len[2] == 1 && dst_len[2] == 1,
        "AverageDownThenRemap expects one-cell-thick source and destination slabs.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_cells.smallEnd(0) == dst_cells.smallEnd(0) &&
        src_cells.smallEnd(1) == dst_cells.smallEnd(1),
        "AverageDownThenRemap requires aligned source/destination slab origins.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_len[0] >= dst_len[0] && src_len[1] >= dst_len[1] &&
        src_len[0] % dst_len[0] == 0 && src_len[1] % dst_len[1] == 0,
        "AverageDownThenRemap requires source/destination slab extents to be integer-ratio compatible.");

    const IntVect ratio(src_len[0] / dst_len[0], src_len[1] / dst_len[1], 1);

    amrex::BoxArray coarsened_src_ba = src.boxArray();
    for (int i = 0; i < coarsened_src_ba.size(); ++i) {
        const amrex::Box original = coarsened_src_ba[i];
        amrex::Box coarsened = original;
        coarsened.coarsen(ratio);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            coarsened.refine(ratio) == original,
            "AverageDownThenRemap source boxes are not evenly coarsenable by the source/destination ratio.");
    }
    coarsened_src_ba.coarsen(ratio);

    const bool direct_average_down_safe = (coarsened_src_ba == dst.boxArray());

    if (direct_average_down_safe) {
        if (src.boxArray().ixType().cellCentered()) {
            amrex::average_down(src, dst, 0, dst.nComp(), ratio);
        } else {
            amrex::average_down_faces(src, dst, ratio, 0);
        }
    } else {
        MultiFab dst_avg(coarsened_src_ba, src.DistributionMap(), dst.nComp(), 0);
        if (src.boxArray().ixType().cellCentered()) {
            amrex::average_down(src, dst_avg, 0, dst.nComp(), ratio);
        } else {
            amrex::average_down_faces(src, dst_avg, ratio, 0);
        }

        dst.ParallelCopy(dst_avg, 0, 0, dst.nComp());
    }
}
}

double
ERF::EvolveOneStep (double /*time*/, double /*dt_request*/)
{
    double cur_time = t_new[0];
    const int step = istep[0];

    if (start_time + cur_time >= stop_time) {
        return 0.0;
    }

    ComputeDt(step);

    int iteration = 1;
    timeStep(0, cur_time, iteration);
    cur_time += dt[0];

    post_timestep(step, cur_time, dt[0]);

    // ****************************************************************************************
    // Write plotfiles at intermediate times
    // ****************************************************************************************
    WriteAtIntermediateTime(step, cur_time);

    return dt[0];
}

void
ERF::ConfigureDriverAtmosToOceanCoupling (bool use_coupling_driver,
                                          bool use_two_way_coupling,
                                          bool use_state_contract)
{
    m_driver_has_atm2ocn_coupling = use_coupling_driver;
    m_driver_uses_two_way_coupling = use_two_way_coupling;
    m_driver_atm2ocn_uses_state_contract = use_state_contract;
}

void
ERF::SetDriverAtmosToOceanStateContract (bool use_state_contract)
{
    m_driver_atm2ocn_uses_state_contract = use_state_contract;
}

void
ERF::GetOceanToAtmosSurfaceLayout (amrex::BoxArray& ba,
                                   amrex::DistributionMapping& dm)
{
    auto* sst_ptr = lsm.Get_Data_Ptr(0, 0);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        sst_ptr != nullptr,
        "ERF::GetOceanToAtmosSurfaceLayout requires OceanSurf level-0 surface storage after InitData.");
    ba = sst_ptr->boxArray();
    dm = sst_ptr->DistributionMap();
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
     This file uses one pack entrypoint for both contract modes. The driver
     passes either the 9-lane state view or the 8-lane flux view; this routine
     branches on the active view shape so the coupling boundary keeps one
     stable solver-facing API.
*/

void
ERF::PackAtmosphericStates (amrex::Vector<amrex::MultiFab*>& states,
                            double /*time*/)
{
    using namespace amrex;

    // Contract slot indices (mirrors ERFRemoraCouplingContract.H; repeated here
    // to avoid a driver→submodule header dependency).
    constexpr int iUwind = 0, iVwind = 1, iPatm = 2, iRH = 3, iTair = 4;
    constexpr int iCloud = 5, iRain  = 6, iSWrad = 7, iLWrad = 8;
    constexpr int nFluxLanes = 8;

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
    const auto& xba2d_lev = amrex::convert(ba2d_lev, IntVect(1,0,0));
    const auto& yba2d_lev = amrex::convert(ba2d_lev, IntVect(0,1,0));
    const int klo = ba.minimalBox().smallEnd(2);

    const bool flux_mode = (states.size() == nFluxLanes);
    if (flux_mode) {
        constexpr int iTauX = 0, iTauY = 1, iSHflux = 2, iLHflux = 3;
        constexpr int iFluxSWrad = 4, iFluxLWrad = 5, iFluxRain = 6, iFluxEvap = 7;
        if (verbose) {
            amrex::Print() << "ERF flux-pack: states.size()=" << states.size()
                           << " has_radiation=" << has_radiation
                           << " tau13_ptr=" << (Tau[lev][TauType::tau13] != nullptr)
                           << " tau23_ptr=" << (Tau[lev][TauType::tau23] != nullptr)
                           << " sfs_hfx3_ptr="
                           << (!SFS_hfx3_lev.empty() && SFS_hfx3_lev[lev] != nullptr)
                           << " sfs_q1fx3_ptr="
                           << (!SFS_q1fx3_lev.empty() && SFS_q1fx3_lev[lev] != nullptr)
                           << " rad_fluxes_ptr="
                           << (!rad_fluxes.empty() && rad_fluxes[lev] != nullptr)
                           << "\n";
            if (Tau[lev][TauType::tau13] != nullptr) {
                PrintFluxLaneStats("tau13_src", *Tau[lev][TauType::tau13]);
            }
            if (Tau[lev][TauType::tau23] != nullptr) {
                PrintFluxLaneStats("tau23_src", *Tau[lev][TauType::tau23]);
            }
            if (!SFS_hfx3_lev.empty() && SFS_hfx3_lev[lev] != nullptr) {
                PrintFluxLaneStats("SFS_hfx3_src", *SFS_hfx3_lev[lev]);
            }
            if (!SFS_q1fx3_lev.empty() && SFS_q1fx3_lev[lev] != nullptr) {
                PrintFluxLaneStats("SFS_q1fx3_src", *SFS_q1fx3_lev[lev]);
            }
            if (has_radiation) {
                PrintFluxLaneStats("rad_fluxes_swup_src", *rad_fluxes[lev], 0);
                PrintFluxLaneStats("rad_fluxes_swdn_src", *rad_fluxes[lev], 1);
                PrintFluxLaneStats("rad_fluxes_lwup_src", *rad_fluxes[lev], 2);
                PrintFluxLaneStats("rad_fluxes_lwdn_src", *rad_fluxes[lev], 3);
            }
        }

        //
        // NOTE: Tau tau13/tau23 already hold CONSERVATIVE stress in [N m-2] --
        //       compute_u_flux in ERF_MOSTStress.H returns rho*<u'w'>, and
        //       ComputeStress_*_N/T scale the strain by rho_bar*mu_eff or by
        //       mu_turb (= rho*K).  They must therefore be exported as-is,
        //       without another factor of density.  This also matches the
        //       SH/LH lanes below, which apply only Cp_d / L_v.
        //
        // NOTE: Tau is only allocated when diffusion is active, so the pointers
        //       have to be checked before they are dereferenced.
        //
        if (iTauX < static_cast<int>(states.size()) && states[iTauX] != nullptr) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Tau[lev][TauType::tau13] != nullptr,
                "Flux-mode coupling of tau_x requires Tau; enable diffusion or a surface_layer bc");
            MultiFab tmp(xba2d_lev, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                auto const& tau13 = Tau[lev][TauType::tau13]->const_array(mfi);
                auto t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = tau13(i,j,klo);
                });
            }
            AverageDownThenRemap(tmp, *states[iTauX]);
            if (verbose) { PrintFluxLaneStats("tau_x_lane", *states[iTauX]); }
        }

        if (iTauY < static_cast<int>(states.size()) && states[iTauY] != nullptr) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Tau[lev][TauType::tau23] != nullptr,
                "Flux-mode coupling of tau_y requires Tau; enable diffusion or a surface_layer bc");
            MultiFab tmp(yba2d_lev, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                auto const& tau23 = Tau[lev][TauType::tau23]->const_array(mfi);
                auto t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = tau23(i,j,klo);
                });
            }
            AverageDownThenRemap(tmp, *states[iTauY]);
            if (verbose) { PrintFluxLaneStats("tau_y_lane", *states[iTauY]); }
        }

        if (iSHflux < static_cast<int>(states.size()) && states[iSHflux] != nullptr &&
            !SFS_hfx3_lev.empty() && SFS_hfx3_lev[lev] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                auto const& hfx = SFS_hfx3_lev[lev]->const_array(mfi);
                auto t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = surface_flux_diagnostics::sensible_heat_flux_wm2_from_rhotheta_flux(
                        hfx(i,j,klo));
                });
            }
            AverageDownThenRemap(tmp, *states[iSHflux]);
            if (verbose) { PrintFluxLaneStats("SHflux_lane", *states[iSHflux]); }
        }

        if (iLHflux < static_cast<int>(states.size()) && states[iLHflux] != nullptr &&
            !SFS_q1fx3_lev.empty() && SFS_q1fx3_lev[lev] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                auto const& qfx = SFS_q1fx3_lev[lev]->const_array(mfi);
                auto t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = surface_flux_diagnostics::latent_heat_flux_wm2_from_rhoqv_flux(
                        qfx(i,j,klo));
                });
            }
            AverageDownThenRemap(tmp, *states[iLHflux]);
            if (verbose) { PrintFluxLaneStats("LHflux_lane", *states[iLHflux]); }
        }

        if (has_radiation) {
            if (iFluxSWrad < static_cast<int>(states.size()) && states[iFluxSWrad] != nullptr) {
                MultiFab tmp(ba2d_lev, dm, 1, 0);
                tmp.ParallelCopy(*rad_fluxes[lev], 1, 0, 1);
                AverageDownThenRemap(tmp, *states[iFluxSWrad]);
                if (verbose) { PrintFluxLaneStats("SWrad_lane", *states[iFluxSWrad]); }
            }
            if (iFluxLWrad < static_cast<int>(states.size()) && states[iFluxLWrad] != nullptr) {
                MultiFab tmp(ba2d_lev, dm, 1, 0);
                for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box bx = mfi.tilebox();
                    auto const& rad_flux = rad_fluxes[lev]->const_array(mfi);
                    auto t = tmp.array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        t(i,j,k) = rad_flux(i,j,k,3) - rad_flux(i,j,k,2);
                    });
                }
                AverageDownThenRemap(tmp, *states[iFluxLWrad]);
                if (verbose) { PrintFluxLaneStats("LWrad_lane", *states[iFluxLWrad]); }
            }
        }

        if (iFluxEvap < static_cast<int>(states.size()) && states[iFluxEvap] != nullptr &&
            !SFS_q1fx3_lev.empty() && SFS_q1fx3_lev[lev] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                auto const& qfx = SFS_q1fx3_lev[lev]->const_array(mfi);
                auto t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = qfx(i,j,klo);
                });
            }
            AverageDownThenRemap(tmp, *states[iFluxEvap]);
            if (verbose) { PrintFluxLaneStats("evap_lane", *states[iFluxEvap]); }
        }

        amrex::ignore_unused(iFluxRain);
        return;
    }

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
            AverageDownThenRemap(u_alias, *states[iUwind]);
        }

        if (iVwind < static_cast<int>(states.size()) && states[iVwind] != nullptr) {
            MultiFab v_alias(uv_slab, amrex::make_alias, 1, 1); // alias v component
            AverageDownThenRemap(v_alias, *states[iVwind]);
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
        AverageDownThenRemap(tmp, *states[iPatm]);
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
        AverageDownThenRemap(tmp, *states[iTair]);
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
                AverageDownThenRemap(tmp, *states[iCloud]);
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
            AverageDownThenRemap(tmp, *states[iSWrad]);
        }
        if (iLWrad < static_cast<int>(states.size()) && states[iLWrad] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            tmp.ParallelCopy(*rad_fluxes[lev], 3, 0, 1);
            AverageDownThenRemap(tmp, *states[iLWrad]);
        }
    }
}

void
ERF::ApplyOceanSurfaceState (const amrex::Vector<amrex::MultiFab*>& state,
                             double time)
{
    if (solverChoice.lsm_type != LandSurfaceType::OceanSurf) {
        return;
    }

    if (!state.empty() && state[0] != nullptr && lsm.Get_Data_Ptr(0, 0) != nullptr) {
        auto* dst = lsm.Get_Data_Ptr(0, 0);
        amrex::MultiFab src_remapped(dst->boxArray(), dst->DistributionMap(), 1, 0);
        src_remapped.ParallelCopy(*state[0], 0, 0, 1);
        const auto lsm_geom = lsm.Get_Lsm_Geom(0);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            lsm_geom.isPeriodic(0) == Geom(0).isPeriodic(0) &&
            lsm_geom.isPeriodic(1) == Geom(0).isPeriodic(1) &&
            lsm_geom.isPeriodic(2) == Geom(0).isPeriodic(2),
            "OceanSurf t_surf geometry lost ERF periodic flags.");
        dst->ParallelCopy(src_remapped, 0, 0, 1);
        dst->FillBoundary(lsm_geom.periodicity());

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
