/*
 * Per-step ERF <-> Noah-MP state exchange: Advance_With_State drives the box loop
 * via interp_from_lev0/time_to_fire/stage_forcing/read_results. Mechanical
 * host<->NoahmpIO copies come from the X-macros in ERF_NOAHMP_Fields.H.
 */

#include <iostream>
#include <cmath>
#include <limits>

#include <AMReX_Print.H>

#include <ERF_NOAHMP.H>
#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include "ERF_NOAHMP_ResultPolicy.H"

using namespace amrex;

// ---------------------------------------------------------------------------
//  Fine-level interpolation (no per-level NetCDF land file)
// ---------------------------------------------------------------------------
void
NOAHMP::interp_from_lev0 (const int& lev,
                          MultiFab& cons_in,
                          const int& nstep)
{
    amrex::ignore_unused(cons_in);

    Print () << "Noah-MP interpolation at level " << lev << " started at time step: " << nstep+1 << std::endl;
    m_updated = true;
    for (int ivar(0); ivar<m_lsm_data_size; ++ivar) {
        InterpFromCoarseLevel(*lsm_fab_data[ivar], lsm_fab_data[ivar]->nGrowVect(),
                              IntVect(0,0,0), // do NOT fill ghost cells outside the domain
                              *lsm_lev0_data[ivar], 0, 0, 1,
                              m_geom0, m_geom,
                              m_refRatio, &cell_cons_interp,
                              m_domain_bcs_type, BCVars::cons_bc);
    }

    // SurfaceLayer filled the noah flux where lsm_undefined existed, so the
    // coarse-grid flux is complete and safe to interpolate.
    for (int ivar(0); ivar<LsmFlux_NOAHMP::NumVars; ++ivar) {
        InterpFromCoarseLevel(*lsm_fab_flux[ivar], lsm_fab_flux[ivar]->nGrowVect(),
                              IntVect(0,0,0), // do NOT fill ghost cells outside the domain
                              *lsm_lev0_flux[ivar], 0, 0, 1,
                              m_geom0, m_geom,
                              m_refRatio, &cell_cons_interp,
                              m_domain_bcs_type, BCVars::cons_bc);
    }
    Print () << "Noah-MP interpolation at level " << lev << " completed" << std::endl;
}

// ---------------------------------------------------------------------------
//  Subcycling gate. Uses the broadcast class members so every rank (including
//  land-free ones) decides identically and enters FillBoundary collectively.
// ---------------------------------------------------------------------------
bool
NOAHMP::time_to_fire (const Real& elapsed_time)
{
    Real NOAH_time = static_cast<Real>(m_itimestep-1) * m_dtbl;
    if (elapsed_time < NOAH_time) {
        m_updated = false;
        return false;
    }

    m_updated = true;
    m_itimestep += 1;   // advance once per firing, in lockstep on every rank
    return true;
}

// ---------------------------------------------------------------------------
//  Steps 1-3: ERF forcing -> pinned input buffer -> NoahmpIO arrays
// ---------------------------------------------------------------------------
void
NOAHMP::stage_forcing (const MFIter& mfi,
                       const erf_noahmp::NoahmpBlockViews& blk,
                       const Box& bx,
                       const int& klo,
                       const int& lev,
                       const bool is_moist,
                       MultiFab& cons_in,
                       MultiFab& xvel_in,
                       MultiFab& yvel_in,
                       const erf_noahmp::PrecipSlots& precip,
                       Vector<erf_noahmp::ClampedPrecipCell>& clamped_cells,
                       Vector<erf_noahmp::InvariantPrecipCell>& invariant_cells)
{
    NoahmpIO_type* noahmpio = blk.io;

    const Array4<const Real>& U_PHY  = xvel_in.const_array(mfi);
    const Array4<const Real>& V_PHY  = yvel_in.const_array(mfi);
    const Array4<const Real>& CONS   = cons_in.const_array(mfi);

    // Forcing pulled from the coupling data fields
    const Array4<const Real>& SWDOWN = lsm_fab_data[LsmData_NOAHMP::sw_flux_dn]->const_array(mfi);
    const Array4<const Real>& GLW    = lsm_fab_data[LsmData_NOAHMP::lw_flux_dn]->const_array(mfi);
    const Array4<const Real>& COSZEN = lsm_fab_data[LsmData_NOAHMP::cos_zenith_angle]->const_array(mfi);

    // Pinned (host-accessible, GPU-usable) staging buffer
    Array4<Real> noah_input_arr = blk.input;

    // Per-slot precip views (current + previous snapshot) and their native->kg/m^2
    // factors. Read only inside the device kernel; absent slots guarded by slot_present.
    const bool hp = precip.have_any;
    GpuArray<Array4<const Real>, NoahmpPrecipSlot::NumSlots> accum_now, accum_prv;
    GpuArray<Real, NoahmpPrecipSlot::NumSlots> accum_fac;
    GpuArray<int,  NoahmpPrecipSlot::NumSlots> slot_present;
    for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
        slot_present[s] = (hp && precip.has[s]) ? 1 : 0;
        accum_fac[s]    = precip.factor[s];
        if (slot_present[s]) {
            accum_now[s] = precip.accum[s]->const_array(mfi);
            accum_prv[s] = m_precip_accum_prev[lev][s]->const_array(mfi);
        }
    }
    const int kklo = klo;

    // (1) Stage ERF forcing into the pinned buffer (device).
    ParallelFor(bx, [=,zero_d=zero] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real qv = (is_moist) ? CONS(i,j,k,RhoQ1_comp)/CONS(i,j,k,Rho_comp) : zero_d;
        noah_input_arr(i,j,0,NoahmpInputComp::u_phy)   = myhalf*(U_PHY(i,j,k)+U_PHY(i+1,j,k));
        noah_input_arr(i,j,0,NoahmpInputComp::v_phy)   = myhalf*(V_PHY(i,j,k)+V_PHY(i  ,j+1,k));
        noah_input_arr(i,j,0,NoahmpInputComp::t_phy)   = getTgivenRandRTh(CONS(i,j,k,Rho_comp),CONS(i,j,k,RhoTheta_comp),qv);
        noah_input_arr(i,j,0,NoahmpInputComp::qv_curr) = qv;
        noah_input_arr(i,j,0,NoahmpInputComp::p8w)     = getPgivenRTh(CONS(i,j,k,RhoTheta_comp),qv);
        noah_input_arr(i,j,0,NoahmpInputComp::swdown)  = SWDOWN(i,j,0);
        noah_input_arr(i,j,0,NoahmpInputComp::glw)     = GLW(i,j,0);
        noah_input_arr(i,j,0,NoahmpInputComp::coszen)  = COSZEN(i,j,0);

        // RAW water-equivalent interval precip [mm]: per slot d = max(0, (now-prev)
        // * native_to_kg_m2). Host applies the guard and derives SR / MP_RAINNC.
        Real drain = zero_d, dsnow = zero_d, dgraup = zero_d;
        if (hp) {
            Real dd[NoahmpPrecipSlot::NumSlots];
            for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
                dd[s] = slot_present[s]
                    ? amrex::max(zero_d, (accum_now[s](i,j,kklo) - accum_prv[s](i,j,kklo)) * accum_fac[s])
                    : zero_d;
            }
            // Hail is folded into MP_GRAUP (no ERF scheme fills the hail slot
            // today); to feed it separately, stage dd[hail] into MP_HAIL instead.
            dsnow  = dd[NoahmpPrecipSlot::snow];
            dgraup = dd[NoahmpPrecipSlot::graupel] + dd[NoahmpPrecipSlot::hail];
            const Real dfroz = dsnow + dgraup;
            // Total: scheme's `total` slot, else rain + frozen.
            drain = slot_present[NoahmpPrecipSlot::total]
                  ? dd[NoahmpPrecipSlot::total]
                  : (dd[NoahmpPrecipSlot::rain] + dfroz);
        }
        noah_input_arr(i,j,0,NoahmpInputComp::rainbl)   = drain;
        noah_input_arr(i,j,0,NoahmpInputComp::mp_snow)  = dsnow;
        noah_input_arr(i,j,0,NoahmpInputComp::mp_graup) = dgraup;
    });

    // (2) Wait for the kernel before host access.
    Gpu::streamSynchronize();

    // (3) Copy pinned staged data to NoahmpIO on the host (3D members take the k/j
    // transpose (i,1,j), 2D take (i,j)); precip is derived below.
    LoopOnCpu(bx, [&,zero_d=zero,one_d=one,lsm_max_precip_interval_d=lsm_max_precip_interval]
                      (int i, int j, int ) noexcept
    {
#define NOAHMP_STAGE_IN_3D(comp,member) noahmpio->member(i,1,j) = noah_input_arr(i,j,0,NoahmpInputComp::comp);
#define NOAHMP_STAGE_IN_2D(comp,member) noahmpio->member(i,j)   = noah_input_arr(i,j,0,NoahmpInputComp::comp);
        NOAHMP_INPUT_3D_FIELDS(NOAHMP_STAGE_IN_3D)
        NOAHMP_INPUT_2D_FIELDS(NOAHMP_STAGE_IN_2D)
#undef NOAHMP_STAGE_IN_3D
#undef NOAHMP_STAGE_IN_2D
        // Guard the RAW total, then derive SR / MP_RAINNC. Clamped cells are
        // reported after the box loop (never silent).
        Real dsnow_h  = noah_input_arr(i,j,0,NoahmpInputComp::mp_snow);
        Real dgraup_h = noah_input_arr(i,j,0,NoahmpInputComp::mp_graup);
        Real drain_h  = noah_input_arr(i,j,0,NoahmpInputComp::rainbl);
        if (drain_h > lsm_max_precip_interval_d) {
            clamped_cells.push_back(erf_noahmp::ClampedPrecipCell{i, j, drain_h});
            drain_h = lsm_max_precip_interval_d;
        }
        // Frozen/total invariant: MP_SNOW+MP_GRAUP are subsets of MP_RAINNC; if
        // their sum exceeds it, rescale the frozen components down and record the cell.
        Real dfroz_h = dsnow_h + dgraup_h;
        const Real inv_tol = Real(1.0e-6) * (one_d + drain_h);
        if (dfroz_h > drain_h + inv_tol) {
            invariant_cells.push_back(erf_noahmp::InvariantPrecipCell{i, j, dfroz_h, drain_h});
            if (dfroz_h > zero_d) {
                const Real scale = drain_h / dfroz_h;   // <= 1
                dsnow_h  *= scale;
                dgraup_h *= scale;
                dfroz_h   = dsnow_h + dgraup_h;         // == drain_h (to round-off)
            }
        }
        noahmpio->RAINBL(i,j)    = drain_h;                                       // [mm]
        noahmpio->SR(i,j)        = (drain_h > zero_d) ? amrex::min(one_d, dfroz_h/drain_h) : zero_d; // [-]
        noahmpio->MP_RAINNC(i,j) = drain_h;                                       // [mm]
        noahmpio->MP_SNOW(i,j)   = dsnow_h;                                       // [mm]
        noahmpio->MP_GRAUP(i,j)  = dgraup_h;                                      // [mm] (includes folded hail)
        // 0 = "no hail" (a real value, not a sentinel; lsm_undefined would corrupt opt_snf=4)
        noahmpio->MP_HAIL(i,j)   = zero_d;                                        // [mm]
    });
}

// ---------------------------------------------------------------------------
//  Steps 5-6: NoahmpIO results -> pinned output buffer -> ERF coupling fields
// ---------------------------------------------------------------------------
void
NOAHMP::read_results (const MFIter& mfi,
                      const erf_noahmp::NoahmpBlockViews& blk,
                      const Box& bx,
                      const Box& gbx,
                      MultiFab& cons_in)
{
    NoahmpIO_type* noahmpio = blk.io;

    // For limiting when populating ghost cells
    int i_lo = bx.smallEnd(0); int i_hi = bx.bigEnd(0);
    int j_lo = bx.smallEnd(1); int j_hi = bx.bigEnd(1);

    const Array4<const Real>& CONS = cons_in.const_array(mfi);

    // Mechanical result destinations (RRTMGP coupling + return-term diagnostics),
    // bound from NOAHMP_RESULT_FIELDS to stay in sync with the scatter below.
#define NOAHMP_BIND_RESULT(alias,lsm,out) Array4<Real> alias = lsm_fab_data[LsmData_NOAHMP::lsm]->array(mfi);
    NOAHMP_RESULT_FIELDS(NOAHMP_BIND_RESULT)
#undef NOAHMP_BIND_RESULT
    // Sentinels (this core does not compute them; always lsm_undefined below).
    Array4<Real> SMSTAV_o    = lsm_fab_data[LsmData_NOAHMP::smstav]->array(mfi);
    Array4<Real> SMSTOT_o    = lsm_fab_data[LsmData_NOAHMP::smstot]->array(mfi);

    // Per-layer soil output fabs (3 groups x m_nsoil) in a device-visible vector so
    // the scatter loops over any NSOIL. Layout: g in {0:smois,1:sh2o,2:tslb}, k in [0,nsoil).
    const int nsoil = m_nsoil;
    const int n_soil_fld = m_num_soil_groups*nsoil;
    Gpu::DeviceVector<Array4<Real>> soil_arr_d(n_soil_fld);
    {
        Vector<Array4<Real>> soil_arr_h(n_soil_fld);
        for (int g(0); g < 3; ++g) {
            for (int k(0); k < nsoil; ++k) {
                soil_arr_h[g*nsoil+k] = lsm_fab_data[soil_data_idx(g,k)]->array(mfi);
            }
        }
        Gpu::copyAsync(Gpu::hostToDevice, soil_arr_h.begin(), soil_arr_h.end(), soil_arr_d.begin());
        Gpu::streamSynchronize();
    }
    Array4<Real>* soil_arr = soil_arr_d.data();
    const int soil_out_base = NoahmpOutputComp::NumComps; // soil outputs start here

    // SurfaceLayer flux destinations
    Array4<Real> q_flux_arr    = lsm_fab_flux[LsmFlux_NOAHMP::q_flux]->array(mfi);
    Array4<Real> t_flux_arr    = lsm_fab_flux[LsmFlux_NOAHMP::t_flux]->array(mfi);
    Array4<Real> tau13_arr     = lsm_fab_flux[LsmFlux_NOAHMP::tau13]->array(mfi);
    Array4<Real> tau23_arr     = lsm_fab_flux[LsmFlux_NOAHMP::tau23]->array(mfi);

    Array4<Real> noah_output_arr = blk.output;

    // (5) Copy NoahmpIO results into the pinned output buffer; banded albedos and
    // the soil profile are read explicitly. SMSTAV/SMSTOT are emitted lsm_undefined.
    LoopOnCpu(bx, [&] (int i, int j, int ) noexcept
    {
#define NOAHMP_READ_OUT_2D(comp,member) noah_output_arr(i,j,0,NoahmpOutputComp::comp) = noahmpio->member(i,j);
        NOAHMP_OUTPUT_2D_FIELDS(NOAHMP_READ_OUT_2D)
#undef NOAHMP_READ_OUT_2D
        // Banded radiation reads (VIS=1, NIR=2)
        noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdir_vis) = noahmpio->ALBSFCDIRXY(i,1,j);
        noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdir_nir) = noahmpio->ALBSFCDIRXY(i,2,j);
        noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdif_vis) = noahmpio->ALBSFCDIFXY(i,1,j);
        noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdif_nir) = noahmpio->ALBSFCDIFXY(i,2,j);
        // Soil profile, one 2D component per layer (Fortran layer index is 1-based).
        for (int k(0); k < nsoil; ++k) {
            noah_output_arr(i,j,0,soil_out_base + 0*nsoil + k) = noahmpio->SMOIS(i,k+1,j);
            noah_output_arr(i,j,0,soil_out_base + 1*nsoil + k) = noahmpio->SH2O (i,k+1,j);
            noah_output_arr(i,j,0,soil_out_base + 2*nsoil + k) = noahmpio->TSLB (i,k+1,j);
        }
    });

    // (6) Copy Noahmp results to ERF (device; from pinned buffer): per-field math
    // (flux ÷rho, -9999 fill guard, sentinels, soil loop), deliberately not table-driven.
    ParallelFor(gbx, [=,Cp_d_d=Cp_d,L_v_d=L_v,lsm_undefined_d=lsm_undefined]
                         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Limit indices to the valid box. FillBoundary will pick these up below.
        int ii = std::min(std::max(i,i_lo),i_hi);
        int jj = std::min(std::max(j,j_lo),j_hi);

        // Noah-MP returns -9999 for cells it does not process (sea-ice/open-water);
        // HFX is Noah-MP's processed-cell gate; it is not a statement that
        // every physical output is valid. tau13/tau23 are nodal in xz/yz;
        // the 2D MFs carry 1 ghost cell for averaging.
        const Real hfx_lsm = noah_output_arr(ii,jj,0,NoahmpOutputComp::hfx);
        const noahmp_result_policy::CellPolicy cell_policy(hfx_lsm);
        // Noah-MP returns HFX,LH [W/m2] and TAU [N/m2] (rho included); divide by
        // rho for the kinematic MOST form ERF stores. Exner factor on t_flux is
        // omitted to match WRF (<~1% error near surface, ~6-7% at p~800 hPa).
        const Real rho_l = CONS(ii,jj,k,Rho_comp);
        const auto fluxes = cell_policy.select_fluxes(
            noah_output_arr(ii,jj,0,NoahmpOutputComp::lh),
            noah_output_arr(ii,jj,0,NoahmpOutputComp::tau_ew),
            noah_output_arr(ii,jj,0,NoahmpOutputComp::tau_ns),
            rho_l, Cp_d_d, L_v_d, lsm_undefined_d);
        t_flux_arr(i,j,k) = fluxes.temperature;
        q_flux_arr(i,j,k) = fluxes.mixing_ratio;
        tau13_arr(i,j,k)  = fluxes.tau_ew;
        tau23_arr(i,j,k)  = fluxes.tau_ns;

        // Processed cells retain field-specific validity checks. Unprocessed
        // cells invalidate the complete provider result set.
        if (cell_policy.processed) {
#define NOAHMP_COPY_RESULT(alias,lsm,out)                                      \
        {                                                                       \
            const Real value = noah_output_arr(ii,jj,0,NoahmpOutputComp::out);  \
            alias(i,j,0) = cell_policy.select(                                  \
                value, lsm_undefined_d, NoahmpOutputComp::out);                \
        }
        NOAHMP_RESULT_FIELDS(NOAHMP_COPY_RESULT)
#undef NOAHMP_COPY_RESULT
            // Not computed by this core -> sentinel, not a misleading 0.
            SMSTAV_o(i,j,0) = lsm_undefined_d;
            SMSTOT_o(i,j,0) = lsm_undefined_d;
            for (int s(0); s < n_soil_fld; ++s) {
                const Real value = noah_output_arr(ii,jj,0,soil_out_base + s);
                soil_arr[s](i,j,0) = cell_policy.select(value, lsm_undefined_d, -1);
            }
        } else {
#define NOAHMP_UNDEF_RESULT(alias,lsm,out) alias(i,j,0) = lsm_undefined_d;
            NOAHMP_RESULT_FIELDS(NOAHMP_UNDEF_RESULT)
#undef NOAHMP_UNDEF_RESULT
            SMSTAV_o(i,j,0) = lsm_undefined_d;
            SMSTOT_o(i,j,0) = lsm_undefined_d;
            for (int s(0); s < n_soil_fld; ++s) {
                soil_arr[s](i,j,0) = lsm_undefined_d;
            }
        }
    });
    // The scatter kernel reads soil_arr_d; sync before it is destroyed.
    Gpu::streamSynchronize();
}

// ---------------------------------------------------------------------------
//  The orchestrator
// ---------------------------------------------------------------------------
void
NOAHMP::Advance_With_State (const int& lev,
                            MultiFab& cons_in,
                            MultiFab& xvel_in,
                            MultiFab& yvel_in,
                            MultiFab* /*hfx3_out*/,
                            MultiFab* /*qfx3_out*/,
                            const SurfacePrecipAccumulationSources& precip_sources,
                            const Real& elapsed_time,
                            const Real& dt,
                            const int& nstep,
                            const bool updated_lev0)
{
    amrex::ignore_unused(dt);

    if (!m_has_nc_file) {
        // Skip interpolation if lev 0 was just updated: Noah runs post-step, so
        // the fluxes & data still hold lsm_undefined values.
        if (!updated_lev0) {
            interp_from_lev0(lev, cons_in, nstep);
        } else {
            m_updated = false;
        }
    } else {
        // Gate on the subcycling schedule (identical on every rank).
        if (!time_to_fire(elapsed_time)) { return; }

        Box domain = m_geom.Domain();

        Print () << "Noah-MP driver at level " << lev << " started at time step: " << nstep+1 << std::endl;

        bool is_moist = (cons_in.nComp() > RhoQ1_comp);

        int klo = domain.smallEnd(2);

        // Collect the typed precip slots, then seed the per-level previous-accumulation
        // snapshots the delta kernel differences against (ERF_NOAHMP_Precip.cpp).
        const erf_noahmp::PrecipSlots precip = collect_precip_sources(precip_sources);
        prepare_precip_snapshots(lev, precip);

        // Cells the guards touched this call; reported after the box loop.
        Vector<erf_noahmp::ClampedPrecipCell>   clamped_cells;
        Vector<erf_noahmp::InvariantPrecipCell> invariant_cells;

        // Loop over blocks: ERF -> Noahmp, drive the land model, Noahmp -> ERF.
        int idb = 0;
        for (MFIter mfi(cons_in); mfi.isValid(); ++mfi, ++idb) {

            Box bx  = mfi.tilebox();
            Box gbx = mfi.tilebox(IntVect(0,0,0),IntVect(1,1,0));

            // Only tiles at the lower z boundary
            if (bx.smallEnd(2) != klo) { continue; }

            bx.makeSlab(2,klo);
            gbx.makeSlab(2,klo);

            // The NoahmpIO block and its pinned staging buffers for this box,
            // all co-indexed by idb, bundled into one handle.
            const erf_noahmp::NoahmpBlockViews blk {
                &noahmpio_vect[idb],
                noahmp_input_tmp[idb]->array(),
                noahmp_output_tmp[idb]->array()
            };

            // (1-3) ERF forcing -> pinned input -> NoahmpIO arrays.
            stage_forcing(mfi, blk, bx, klo, lev, is_moist,
                          cons_in, xvel_in, yvel_in,
                          precip, clamped_cells, invariant_cells);

            // (4) Drive Noah-MP. Mirror the authoritative counter into the block first.
            blk.io->itimestep = m_itimestep;
            blk.io->DriverMain();

            // (5-6) NoahmpIO results -> pinned output -> ERF coupling fields.
            read_results(mfi, blk, bx, gbx, cons_in);
        }

        // Advance the snapshots now that this call's RAINBL/SR have been consumed.
        advance_precip_snapshots(lev, precip);

        // Report any cells the guards clamped / rescaled (never silent).
        report_precip_diagnostics(lev, clamped_cells, invariant_cells);

        for (auto ivar = 0; ivar < LsmFlux_NOAHMP::NumVars; ++ivar) {
            lsm_fab_flux[ivar]->FillBoundary(m_geom.periodicity());
        }
        Print () << "Noah-MP driver at level " << lev << " completed" << std::endl;
    }
}
