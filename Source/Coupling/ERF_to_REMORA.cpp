#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_MicrophysicsUtils.H>
#include <Diagnostics/ERF_SurfaceFluxDiagnostics.H>

#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <limits>

namespace
{
void
PrintFluxLaneStats (const char* label, const amrex::MultiFab& mf, int comp = 0)
{
    amrex::Print() << "ERF flux-pack " << label
                   << ": min=" << mf.min(comp)
                   << " max=" << mf.max(comp) << "\n";
}

// -------------------------------------------------------------------------
// NEW: Conservative Sparse Matrix Remap Engine
// -------------------------------------------------------------------------

// Source index region each destination box's stencils actually reference.
//
// The entries in index_mf are *source* indices. For non-conformal grids they
// routinely fall far outside the destination box's own index range: 637 ERF cells
// against 563 REMORA cells in x puts the high-corner stencil roughly 70 cells
// past the destination box, so staging the source on `dst.boxArray()` plus a
// fixed ghost halo leaves those reads outside the valid region entirely - zero at
// best, out of bounds at worst. Deriving the staged region from the indices in
// hand is what makes this correct for an arbitrary grid pair.
//
// The result becomes a BoxArray, which must be globally consistent, hence the
// reduction. Recomputed per call rather than cached: the scan is one pass over the
// index arrays, small next to a solver step, and a cache would need invalidating
// on every regrid and weight rebuild.
amrex::BoxArray
StagedSourceBoxArray (const amrex::MultiFab& src,
                      const amrex::MultiFab& dst,
                      const amrex::iMultiFab& index_mf,
                      int max_stencil_size)
{
    using namespace amrex;

    const Box src_domain = src.boxArray().minimalBox();
    const int nboxes = static_cast<int>(dst.boxArray().size());
    constexpr int int_big = std::numeric_limits<int>::max();

    Vector<int> lo(2 * nboxes,  int_big);
    Vector<int> hi(2 * nboxes, -int_big);

    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        const int b = mfi.index();
        const Box bx = mfi.validbox();
        auto const& idx = index_mf.const_array(mfi);

        ReduceOps<ReduceOpMin, ReduceOpMin, ReduceOpMax, ReduceOpMax> reduce_op;
        ReduceData<int, int, int, int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int i_min = int_big, j_min = int_big, i_max = -int_big, j_max = -int_big;
            for (int m = 0; m < max_stencil_size; ++m) {
                const int si = idx(i, j, k, m * 3);
                const int sj = idx(i, j, k, m * 3 + 1);
                if (si < 0 || sj < 0) { continue; }   // -1 marks an unused slot
                i_min = amrex::min(i_min, si); i_max = amrex::max(i_max, si);
                j_min = amrex::min(j_min, sj); j_max = amrex::max(j_max, sj);
            }
            return {i_min, j_min, i_max, j_max};
        });

        auto const& hv = reduce_data.value(reduce_op);
        lo[2*b  ] = amrex::get<0>(hv); lo[2*b+1] = amrex::get<1>(hv);
        hi[2*b  ] = amrex::get<2>(hv); hi[2*b+1] = amrex::get<3>(hv);
    }

    ParallelDescriptor::ReduceIntMin(lo.dataPtr(), static_cast<int>(lo.size()));
    ParallelDescriptor::ReduceIntMax(hi.dataPtr(), static_cast<int>(hi.size()));

    // The staged boxes must carry the source's index type from the start: the tau
    // lanes hand us a face-staggered source, and intersecting a cell-centered box
    // with a face-staggered one trips Box::operator&='s sameType assertion.
    const IndexType src_ixtype = src.boxArray().ixType();

    BoxList bl(src_ixtype);
    for (int b = 0; b < nboxes; ++b) {
        Box need;
        if (lo[2*b] > hi[2*b] || lo[2*b+1] > hi[2*b+1]) {
            // No stencil anywhere in this destination box. A BoxArray cannot hold
            // an empty box, so stage a single cell; nothing reads it.
            need = Box(src_domain.smallEnd(), src_domain.smallEnd(), src_ixtype);
        } else {
            need = Box(IntVect(lo[2*b], lo[2*b+1], src_domain.smallEnd(2)),
                       IntVect(hi[2*b], hi[2*b+1], src_domain.bigEnd(2)),
                       src_ixtype);
            // The intersection below silently drops any index outside the source
            // domain, but the stencil kernel dereferences it regardless - so a
            // weight generator that emitted an out-of-domain index would turn into
            // an out-of-bounds read instead of a diagnosable failure. Assert first.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                src_domain.contains(need),
                "StagedSourceBoxArray: a remap stencil references a source index outside the source domain.");
            need &= src_domain;
        }
        bl.push_back(need);
    }
    return BoxArray(std::move(bl));
}

// fallback_val fills destination cells that no source cell overlaps.
//
// Zero is the wrong default for a physical field: the receiving model cannot tell
// "no data here" from "the atmosphere says zero", and for an intensive quantity
// zero is not merely inaccurate, it drives the receiver's flux formulas outside
// their valid domain. A lane of 0 Pa / 0 K arrives in REMORA as PairM = 0 mb and
// (after the K->C conversion) TairK = 0.01 K, which poisons the bulk-flux chain -
// candidates include sqrt(vap_p) at REMORA_bulk_flux.cpp:145 and
// pow(Bf*blk_Zabl, 1/3) at :314 - and any NaN produced there survives every
// subsequent land-mask multiplication. Uncovered cells are
// unavoidable wherever the two land masks disagree along a coastline or the ocean
// grid reaches past the atmosphere's footprint, so every lane carrying an
// intensive field needs a physically admissible value here.
void
ApplyConservativeRemap (const amrex::MultiFab& src,
                        amrex::MultiFab& dst,
                        const amrex::MultiFab& weight_mf,
                        const amrex::iMultiFab& index_mf,
                        int max_stencil_size,
                        const amrex::MultiFab* dst_mask = nullptr,
                        amrex::Real fallback_val = amrex::Real(0.0))
{
    using namespace amrex;

    // 1. Data Routing: Get the ERF source data onto the REMORA destination layout,
    // staging exactly the source region the stencils reference so every indexed
    // read below lands in valid data.
    MultiFab src_on_dst(StagedSourceBoxArray(src, dst, index_mf, max_stencil_size),
                        dst.DistributionMap(), src.nComp(), 0);
    src_on_dst.setVal(0.0);
    src_on_dst.ParallelCopy(src);

    // Ensure the destination is zeroed out before accumulating the weighted sum
    dst.setVal(0.0);

    // 2. Stencil Application: Execute the local sparse dot product
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();

        auto const& w_arr   = weight_mf.const_array(mfi);
        auto const& idx_arr = index_mf.const_array(mfi);
        auto const& src_arr = src_on_dst.const_array(mfi);
        auto        dst_arr = dst.array(mfi);
        const bool has_mask = (dst_mask != nullptr);
        auto const& mask_arr = has_mask ? dst_mask->const_array(mfi) : Array4<const Real>{};

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // No stencil entry at all means no source cell overlapped this
            // destination cell. Hand back the fallback rather than an accumulated
            // zero, and do not apply the mask to it: a masked-out cell still gets
            // evaluated by the receiving model's flux formulas before the mask is
            // applied, so it too must hold an admissible value.
            if (idx_arr(i, j, k, 0) < 0) {
                dst_arr(i, j, k) = fallback_val;
                return;
            }

            Real sum = 0.0;
            for (int m = 0; m < max_stencil_size; ++m) {
                Real w = w_arr(i, j, k, m);
                if (w > 0.0) {
                    int src_i = idx_arr(i, j, k, m * 3);
                    int src_j = idx_arr(i, j, k, m * 3 + 1);
                    int src_k = idx_arr(i, j, k, m * 3 + 2);

                    sum += w * src_arr(src_i, src_j, src_k);
                }
            }
            // Blend toward the fallback rather than multiplying by the mask.
            //
            // Multiplying drives masked cells to exactly zero, which is right for
            // a flux lane (fallback_val = 0, so this reduces to sum *= mask and is
            // bit-identical) but destructive for an intensive state lane: a land
            // cell would hold 0 Pa and 0 K, and the receiving model evaluates its
            // bulk formulas over land before applying its own mask, so those
            // zeros become NaN and the NaN outlives the masking. Land cells still
            // receive no coupled information either way; they just stay
            // physically admissible.
            if (has_mask) {
                const Real mask = mask_arr(i, j, k);
                sum = mask * sum + (Real(1.0) - mask) * fallback_val;
            }
            dst_arr(i, j, k) = sum;
        });
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

void
ERF::GetOceanToAtmosCornerCoordinates (const amrex::MultiFab*& x_corner,
                                       const amrex::MultiFab*& y_corner) const
{
    if (lon_m.empty() || lat_m.empty() ||
        lon_m[0] == nullptr || lat_m[0] == nullptr) {
        x_corner = nullptr;
        y_corner = nullptr;
        return;
    }
    x_corner = lon_m[0].get();
    y_corner = lat_m[0].get();
}

void
ERF::GetLandMask (const amrex::iMultiFab*& lmask) const
{
    if (lmask_lev.empty() || lmask_lev[0].empty() || lmask_lev[0][0] == nullptr) {
        lmask = nullptr;
        return;
    }
    lmask = lmask_lev[0][0].get();
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
                            double /*time*/,
                            const amrex::Vector<const amrex::MultiFab*>& weight_mf_by_family,
                            const amrex::Vector<const amrex::iMultiFab*>& index_mf_by_family,
                            int max_stencil_size,
                            const amrex::MultiFab* dst_mskr,
                            const amrex::MultiFab* dst_msku,
                            const amrex::MultiFab* dst_mskv)
{
    using namespace amrex;

    // WeightFamily slot order, mirrored from ERFRemoraMultiBlockContainer.H:
    // {SCALAR_CENTERED, TAUX_FACE, TAUY_FACE, WIND_FACE_TO_CENTER}.
    constexpr int kScalarFamily = 0, kTauXFamily = 1, kTauYFamily = 2, kWindFamily = 3;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        weight_mf_by_family.size() == 4 && index_mf_by_family.size() == 4,
        "PackAtmosphericStates requires 4 WeightFamily-keyed weight/index pairs.");
    for (int f = 0; f < 4; ++f) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            weight_mf_by_family[f] != nullptr && index_mf_by_family[f] != nullptr,
            "PackAtmosphericStates requires valid driver remap stencil components for every WeightFamily slot.");
        // The stencil loops read max_stencil_size weights and 3*max_stencil_size
        // index components. Nothing else checks that the driver allocated that
        // many, so a driver/submodule version skew would read past the component
        // count and silently return garbage rather than failing.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            weight_mf_by_family[f]->nComp() >= max_stencil_size &&
            index_mf_by_family[f]->nComp() >= 3 * max_stencil_size,
            "PackAtmosphericStates: driver remap stencil arrays have too few components for max_stencil_size.");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            weight_mf_by_family[f]->boxArray() == index_mf_by_family[f]->boxArray() &&
            weight_mf_by_family[f]->DistributionMap() == index_mf_by_family[f]->DistributionMap(),
            "PackAtmosphericStates: weight and index arrays for a WeightFamily slot disagree on layout.");
    }
    const amrex::MultiFab* weight_mf = weight_mf_by_family[kScalarFamily];
    const amrex::iMultiFab* index_mf = index_mf_by_family[kScalarFamily];

    // Contract slot indices (mirrors ERFRemoraCouplingContract.H; repeated here
    // to avoid a driver→submodule header dependency).
    constexpr int iUwind = 0, iVwind = 1, iPatm = 2, iRH = 3, iTair = 4;

    // Values for destination cells no source cell overlaps (see
    // ApplyConservativeRemap). These are standard-atmosphere placeholders, not
    // physics: they exist so the receiving model's flux formulas stay in their
    // valid domain over cells the atmosphere grid never covered. Wind is the one
    // field where zero is genuinely meaningful (calm), so it needs no special
    // value. Units are the driver contract's: Pa and K, converted downstream.
    constexpr Real fallback_pressure_pa = Real(101325.0);   // sea-level standard
    constexpr Real fallback_temp_k      = Real(288.15);     // 15 C, standard
    constexpr Real fallback_wind_ms     = Real(0.0);        // calm
    constexpr int iCloud = 5, iRain  = 6, iSWrad = 7, iLWrad = 8;
    constexpr int nFluxLanes = 8;

    const int lev   = 0;

    auto& cons = vars_new[lev][Vars::cons];
    auto& xvel = vars_new[lev][Vars::xvel]; // XFace
    auto& yvel = vars_new[lev][Vars::yvel]; // YFace

    const bool has_moisture  = (solverChoice.moisture_type != MoistureType::None);
    // Level 0 only, matching HasRadiation(); lev is 0 throughout this routine.
    const bool has_radiation = HasRadiation();

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

        // TAU LANE STAGGER: unresolved, see the parked block below.
        //
        // ERF tau13/tau23 live on x-/y-faces and REMORA vec_sustr/vec_svstr live on
        // u-/v-faces, so the physically meaningful transfer is face -> face: a
        // stagger-matched transfer with no averaging on a conformal grid. That is
        // what the live code below does, and what the driver's WeightFamily comment
        // describes.
        //
        // The defect is on the weight side, not here: the driver builds its remap
        // stencils from source *cell* polygons and stores ERF *cell* indices for
        // every family (ERFRemoraMultiBlockContainer.cpp, the ExtractCellQuadFromCoords
        // / CellBoundsFromGeom calls in BuildAtmosToOceanWeights). So the index this
        // face-staggered source is read with was computed for a cell, displacing tau
        // by ~dx/2 and leaving the outermost face column unreachable.
        //
        // The parked block below "fixed" that by averaging tau to cell centres so the
        // source would match the cell-based weights. That is the wrong direction: it
        // smooths the stress twice (face -> cell, then cell -> face via area weights)
        // and throws away the face alignment. Keeping it only as a record of the
        // attempt. The correct fix is to build ERF *face* control-volume source
        // polygons for the TAUX/TAUY families in the driver, so the weights match the
        // face source that is already being passed here.
#if 0
        if (iTauX < static_cast<int>(states.size()) && states[iTauX] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                auto const& tau13 = Tau[lev][TauType::tau13]->const_array(mfi);
                auto const& c = cons.const_array(mfi);
                auto t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = c(i,j,klo,Rho_comp)
                             * Real(0.5) * (tau13(i,j,klo) + tau13(i+1,j,klo));
                });
            }
            ApplyConservativeRemap(tmp, *states[iTauX], *weight_mf_by_family[kTauXFamily],
                                   *index_mf_by_family[kTauXFamily], max_stencil_size, dst_msku);
        }
        if (iTauY < static_cast<int>(states.size()) && states[iTauY] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                auto const& tau23 = Tau[lev][TauType::tau23]->const_array(mfi);
                auto const& c = cons.const_array(mfi);
                auto t = tmp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    t(i,j,k) = c(i,j,klo,Rho_comp)
                             * Real(0.5) * (tau23(i,j,klo) + tau23(i,j+1,klo));
                });
            }
            ApplyConservativeRemap(tmp, *states[iTauY], *weight_mf_by_family[kTauYFamily],
                                   *index_mf_by_family[kTauYFamily], max_stencil_size, dst_mskv);
        }
#endif

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
            ApplyConservativeRemap(tmp, *states[iTauX], *weight_mf_by_family[kTauXFamily],
                                   *index_mf_by_family[kTauXFamily], max_stencil_size, dst_msku);
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
            ApplyConservativeRemap(tmp, *states[iTauY], *weight_mf_by_family[kTauYFamily],
                                   *index_mf_by_family[kTauYFamily], max_stencil_size, dst_mskv);
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
            ApplyConservativeRemap(tmp, *states[iSHflux], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
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
            ApplyConservativeRemap(tmp, *states[iLHflux], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
            if (verbose) { PrintFluxLaneStats("LHflux_lane", *states[iLHflux]); }
        }

        if (has_radiation) {
            if (iFluxSWrad < static_cast<int>(states.size()) && states[iFluxSWrad] != nullptr) {
                MultiFab tmp(ba2d_lev, dm, 1, 0);
                tmp.ParallelCopy(*rad_fluxes[lev], 1, 0, 1);
                ApplyConservativeRemap(tmp, *states[iFluxSWrad], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
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
                ApplyConservativeRemap(tmp, *states[iFluxLWrad], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
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
            ApplyConservativeRemap(tmp, *states[iFluxEvap], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
            if (verbose) { PrintFluxLaneStats("evap_lane", *states[iFluxEvap]); }
        }

#if 0
        // Extract Rain Flux using surface qr and a proxy terminal velocity
        if (has_moisture && iFluxRain < static_cast<int>(states.size()) && states[iFluxRain] != nullptr) {
            int qr_idx = solverChoice.moisture_indices.qr;
            if (qr_idx != -1) {
                MultiFab tmp(ba2d_lev, dm, 1, 0);
                for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box bx = mfi.tilebox();
                    auto const& c = cons.const_array(mfi);
                    auto t = tmp.array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        // Rain flux proxy: (rho * qr) * V_terminal (approx 5.0 m/s for rain) -> kg m^-2 s^-1
                        t(i,j,k) = c(i,j,klo,qr_idx) * Real(5.0);
                    });
                }
                ApplyConservativeRemap(tmp, *states[iFluxRain], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
                if (verbose) { PrintFluxLaneStats("rain_lane", *states[iFluxRain]); }
            }
        }
#endif
        amrex::ignore_unused(iFluxRain);
        return;
    }

    // --- Uwind + Vwind ---
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
            ApplyConservativeRemap(u_alias, *states[iUwind], *weight_mf_by_family[kWindFamily],
                                   *index_mf_by_family[kWindFamily], max_stencil_size,
                                   nullptr, fallback_wind_ms);
        }

        if (iVwind < static_cast<int>(states.size()) && states[iVwind] != nullptr) {
            MultiFab v_alias(uv_slab, amrex::make_alias, 1, 1); // alias v component
            ApplyConservativeRemap(v_alias, *states[iVwind], *weight_mf_by_family[kWindFamily],
                                   *index_mf_by_family[kWindFamily], max_stencil_size,
                                   nullptr, fallback_wind_ms);
        }
    }

    // --- Patm ---
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
        ApplyConservativeRemap(tmp, *states[iPatm], *weight_mf, *index_mf, max_stencil_size, dst_mskr,
                               fallback_pressure_pa);
    }

    // --- Tair ---
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
        ApplyConservativeRemap(tmp, *states[iTair], *weight_mf, *index_mf, max_stencil_size, dst_mskr,
                               fallback_temp_k);
    }

    // --- Humidity/Cloud/Rain ---
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
                ApplyConservativeRemap(tmp, *states[iCloud], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
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

    // --- Radiation ---
    if (has_radiation) {
        if (iSWrad < static_cast<int>(states.size()) && states[iSWrad] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            tmp.ParallelCopy(*rad_fluxes[lev], 1, 0, 1);
            ApplyConservativeRemap(tmp, *states[iSWrad], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
        }
        if (iLWrad < static_cast<int>(states.size()) && states[iLWrad] != nullptr) {
            MultiFab tmp(ba2d_lev, dm, 1, 0);
            tmp.ParallelCopy(*rad_fluxes[lev], 3, 0, 1);
            ApplyConservativeRemap(tmp, *states[iLWrad], *weight_mf, *index_mf, max_stencil_size, dst_mskr);
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
