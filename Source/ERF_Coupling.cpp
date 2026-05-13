#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>

#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

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
    const int k_atm = 0; // lowest ERF cell = atmospheric surface layer

    auto& cons = vars_new[lev][Vars::cons];
    auto& xvel = vars_new[lev][Vars::xvel]; // XFace
    auto& yvel = vars_new[lev][Vars::yvel]; // YFace

    const bool has_moisture  = (solverChoice.moisture_type != MoistureType::None);
    const bool has_radiation = (!rad_fluxes.empty() && rad_fluxes[lev] != nullptr);

    const auto& ba = cons.boxArray();
    const auto& dm = cons.DistributionMap();

    // --- Uwind: average xvel XFace → cell center at k=0 ---
    if (iUwind < static_cast<int>(states.size()) && states[iUwind] != nullptr) {
        MultiFab tmp(ba, dm, 1, 0);
        for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = makeSlab(mfi.validbox(), 2, k_atm);
            auto const& u = xvel.const_array(mfi);
            auto         w = tmp.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                w(i,j,k) = 0.5_rt * (u(i,j,k) + u(i+1,j,k));
            });
        }
        states[iUwind]->ParallelCopy(tmp, 0, 0, 1);
    }

    // --- Vwind: average yvel YFace → cell center at k=0 ---
    if (iVwind < static_cast<int>(states.size()) && states[iVwind] != nullptr) {
        MultiFab tmp(ba, dm, 1, 0);
        for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = makeSlab(mfi.validbox(), 2, k_atm);
            auto const& v = yvel.const_array(mfi);
            auto         w = tmp.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                w(i,j,k) = 0.5_rt * (v(i,j,k) + v(i,j+1,k));
            });
        }
        states[iVwind]->ParallelCopy(tmp, 0, 0, 1);
    }

    // --- Patm: getPgivenRTh(RhoTheta, qv) at k=0 ---
    if (iPatm < static_cast<int>(states.size()) && states[iPatm] != nullptr) {
        MultiFab tmp(ba, dm, 1, 0);
        for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = makeSlab(mfi.validbox(), 2, k_atm);
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
        states[iPatm]->ParallelCopy(tmp, 0, 0, 1);
    }

    // --- Tair: getTgivenRandRTh(rho, RhoTheta, qv) at k=0 [K] ---
    if (iTair < static_cast<int>(states.size()) && states[iTair] != nullptr) {
        MultiFab tmp(ba, dm, 1, 0);
        for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = makeSlab(mfi.validbox(), 2, k_atm);
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
        states[iTair]->ParallelCopy(tmp, 0, 0, 1);
    }

    // --- Moisture fields: from cons when available, else REMORA inputs-file constants ---
    if (has_moisture) {
        if (iRH < static_cast<int>(states.size()) && states[iRH] != nullptr) {
            MultiFab tmp(ba, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = makeSlab(mfi.validbox(), 2, k_atm);
                auto const& c = cons.const_array(mfi);
                auto         t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = c(i,j,k,RhoQ1_comp) / c(i,j,k,Rho_comp);
                });
            }
            states[iRH]->ParallelCopy(tmp, 0, 0, 1);
        }
        if (iCloud < static_cast<int>(states.size()) && states[iCloud] != nullptr) {
            MultiFab tmp(ba, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = makeSlab(mfi.validbox(), 2, k_atm);
                auto const& c = cons.const_array(mfi);
                auto         t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = (c(i,j,k,RhoQ2_comp) + c(i,j,k,RhoQ3_comp)) / c(i,j,k,Rho_comp);
                });
            }
            states[iCloud]->ParallelCopy(tmp, 0, 0, 1);
        }
        if (iRain < static_cast<int>(states.size()) && states[iRain] != nullptr) {
            MultiFab tmp(ba, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = makeSlab(mfi.validbox(), 2, k_atm);
                auto const& c = cons.const_array(mfi);
                auto         t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = c(i,j,k,RhoQ4_comp) / c(i,j,k,Rho_comp);
                });
            }
            states[iRain]->ParallelCopy(tmp, 0, 0, 1);
        }
    } else {
        // No prognostic moisture in ERF: fill slabs with REMORA's inputs-file
        // bulk-flux constants so the coupled path sees physically consistent values.
        amrex::Real air_humidity_fallback = 0.0, cloud_fallback = 0.0, rain_fallback = 0.0;
        amrex::ParmParse("remora").query("air_humidity", air_humidity_fallback);
        amrex::ParmParse("remora").query("cloud",        cloud_fallback);
        amrex::ParmParse("remora").query("rain",         rain_fallback);
        if (iRH    < static_cast<int>(states.size()) && states[iRH]    != nullptr)
            states[iRH]->setVal(air_humidity_fallback);
        if (iCloud < static_cast<int>(states.size()) && states[iCloud] != nullptr)
            states[iCloud]->setVal(cloud_fallback);
        if (iRain  < static_cast<int>(states.size()) && states[iRain]  != nullptr)
            states[iRain]->setVal(rain_fallback);
    }

    // --- Radiation: sw_flux_dn (comp=1) and lw_flux_dn (comp=3) from rad_fluxes ---
    // When radiation is absent, zero-fill so REMORA retains its own initialized defaults.
    if (has_radiation) {
        if (iSWrad < static_cast<int>(states.size()) && states[iSWrad] != nullptr)
            states[iSWrad]->ParallelCopy(*rad_fluxes[lev], 1, 0, 1);
        if (iLWrad < static_cast<int>(states.size()) && states[iLWrad] != nullptr)
            states[iLWrad]->ParallelCopy(*rad_fluxes[lev], 3, 0, 1);
    } else {
        if (iSWrad < static_cast<int>(states.size()) && states[iSWrad] != nullptr)
            states[iSWrad]->setVal(0.0);
        if (iLWrad < static_cast<int>(states.size()) && states[iLWrad] != nullptr)
            states[iLWrad]->setVal(0.0);
    }
}

void
ERF::ApplyOceanSurfaceState (const amrex::Vector<amrex::MultiFab*>& state,
                             amrex::Real time)
{
    amrex::ignore_unused(time);

    // Example (legacy state-passing): state[0] carries SST [K].
    // We intentionally copy only one horizontal slab, derived from the
    // interface face convention:
    // - ocean/atmos interface face index on source: k_face = src.bigEnd(2) + 1
    // - source cell below that face (ocean top cell): src_k = k_face - 1
    // - destination LSM surface-facing slab uses its bottom-most k index
    //   (for current level-0 matched-grid tests this is the interface slab).
    // This encodes the physical alignment note from coupled discussions:
    // ERF k=0 atmospheric cell lies above REMORA k=nz-1 ocean cell.
    // For initial matched-grid tests, require identical horizontal BoxArray/DM.
    if (solverChoice.lsm_type == LandSurfaceType::None) {
        return;
    }

    if (!state.empty() && state[0] != nullptr && lsm.Get_Data_Ptr(0, 0) != nullptr) {
        auto* dst = lsm.Get_Data_Ptr(0, 0);
        const int dst_k = dst->boxArray().minimalBox().smallEnd(2);
        const int src_k  = state[0]->boxArray().minimalBox().bigEnd(2);
        for (amrex::MFIter mfi(*dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx = amrex::makeSlab(mfi.validbox(), 2, dst_k);
            auto dst_arr = dst->array(mfi);
            auto src_arr = state[0]->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
                dst_arr(i,j,dst_k) = src_arr(i,j,src_k);
            });
        }
    }
}
