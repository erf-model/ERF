/*
 * The per-step ERF <-> Noah-MP state exchange: NOAHMP::Advance_With_State and
 * its pipeline helpers.
 *
 *   Advance_With_State  -- the orchestrator (firing gate, box loop, ghost fill)
 *   interp_from_lev0    -- fine-level fill from level 0 (no per-level land file)
 *   time_to_fire        -- the Noah-MP subcycling gate (dev/spec-noahmp-api.md §5)
 *   stage_forcing       -- gpu-spec steps 1-3: ERF -> pinned input -> NoahmpIO
 *   read_results        -- gpu-spec steps 5-6: NoahmpIO -> pinned output -> ERF
 *
 * The mechanical, 1:1 host<->NoahmpIO copies are generated from the X-macro
 * field registry in ERF_NOAHMP_Fields.H so the enum and the copy loops cannot
 * drift. Computed forcing (winds/EOS), precip, banded albedos, the flux ÷rho +
 * -9999 fill guard, and the soil-layer loops are NOT table-driven and stay
 * explicit -- see dev/spec-noahmp-gpu.md and dev/spec-noahmp-reorg.md §4.
 */

#include <iostream>
#include <limits>

#include <AMReX_Print.H>

#include <ERF_NOAHMP.H>
#include <ERF_Constants.H>
#include <ERF_EOS.H>

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

    // NOTE: Do not try to interpolate if lev 0 was just updated. Since Noah is
    //       called post-step the fluxes & data will contain lsm_undefined values.
    Print () << "Noah-MP interpolation at level " << lev << " started at time step: " << nstep+1 << std::endl;
    m_updated = true;
    for (int ivar(0); ivar<m_lsm_data_size; ++ivar) {
        // Interpolate from lev 0 to obtain the lsm data
        InterpFromCoarseLevel(*lsm_fab_data[ivar], lsm_fab_data[ivar]->nGrowVect(),
                              IntVect(0,0,0), // do NOT fill ghost cells outside the domain
                              *lsm_lev0_data[ivar], 0, 0, 1,
                              m_geom0, m_geom,
                              m_refRatio, &cell_cons_interp,
                              m_domain_bcs_type, BCVars::cons_bc);
    } // ivar

    // NOTE: The surface layer class wrote into the noah flux data structures
    //       where lsm_undefined values existed. This makes the noah flux
    //       complete on the coarse grid and we can safely interpolate.
    for (int ivar(0); ivar<LsmFlux_NOAHMP::NumVars; ++ivar) {
        // Interpolate from lev 0 to obtain the lsm fluxes
        InterpFromCoarseLevel(*lsm_fab_flux[ivar], lsm_fab_flux[ivar]->nGrowVect(),
                              IntVect(0,0,0), // do NOT fill ghost cells outside the domain
                              *lsm_lev0_flux[ivar], 0, 0, 1,
                              m_geom0, m_geom,
                              m_refRatio, &cell_cons_interp,
                              m_domain_bcs_type, BCVars::cons_bc);
    } // ivar
    Print () << "Noah-MP interpolation at level " << lev << " completed" << std::endl;
}

// ---------------------------------------------------------------------------
//  Subcycling gate. Identical on every rank (uses the broadcast class members)
//  so the FillBoundary at the end of Advance is entered collectively by all.
// ---------------------------------------------------------------------------
bool
NOAHMP::time_to_fire (const Real& elapsed_time)
{
    // Verify we need to take another LSM step. Use the class-level counter/dtbl
    // (valid on land-free ranks) so every rank decides identically.
    Real NOAH_time = static_cast<Real>(m_itimestep-1) * m_dtbl;
    if (elapsed_time < NOAH_time) {
        m_updated = false;
        return false;
    }

    // We are updating
    m_updated = true;

    // Advance the counter once per firing, in lockstep on every rank.
    m_itimestep += 1;
    return true;
}

// ---------------------------------------------------------------------------
//  Steps 1-3: ERF forcing -> pinned input buffer -> NoahmpIO arrays
// ---------------------------------------------------------------------------
void
NOAHMP::stage_forcing (const MFIter& mfi,
                       const int& idb,
                       const Box& bx,
                       const int& klo,
                       const int& lev,
                       const bool is_moist,
                       MultiFab& cons_in,
                       MultiFab& xvel_in,
                       MultiFab& yvel_in,
                       const noahmp_detail::PrecipSlots& precip,
                       Vector<noahmp_detail::ClampedPrecipCell>& clamped_cells,
                       Vector<noahmp_detail::InvariantPrecipCell>& invariant_cells)
{
    NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

    const Array4<const Real>& U_PHY  = xvel_in.const_array(mfi);
    const Array4<const Real>& V_PHY  = yvel_in.const_array(mfi);
    const Array4<const Real>& CONS   = cons_in.const_array(mfi);

    // Into NOAH-MP (forcing pulled from the coupling data fields)
    const Array4<const Real>& SWDOWN = lsm_fab_data[LsmData_NOAHMP::sw_flux_dn]->const_array(mfi);
    const Array4<const Real>& GLW    = lsm_fab_data[LsmData_NOAHMP::lw_flux_dn]->const_array(mfi);
    const Array4<const Real>& COSZEN = lsm_fab_data[LsmData_NOAHMP::cos_zenith_angle]->const_array(mfi);

    // Use The_Pinned_Arena() for host-accessible memory that can be used with GPU
    Array4<Real> noah_input_arr = noahmp_input_tmp[idb]->array();

    // Per-slot precip accumulation views (current + previous-call snapshot) and
    // their native->kg/m^2 factors, for building the water-equivalent interval
    // precip. These MultiFabs live in the (device) default arena, so they are read
    // ONLY inside the device ParallelFor below and staged into the pinned
    // noah_input_arr -- never dereferenced on the host (that segfaults on GPU).
    // have_any is false for schemes with no precip (RAINBL stays 0). GpuArrays are
    // captured by value into the device kernel; absent slots keep null Array4s that
    // are guarded by slot_present before any dereference.
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

    // (1) Copy forcing data from ERF to Noahmp (device; stage into pinned buffer).
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real qv = (is_moist) ? CONS(i,j,k,RhoQ1_comp)/CONS(i,j,k,Rho_comp) : zero;
        noah_input_arr(i,j,0,NoahmpInputComp::u_phy)   = myhalf*(U_PHY(i,j,k)+U_PHY(i+1,j,k));
        noah_input_arr(i,j,0,NoahmpInputComp::v_phy)   = myhalf*(V_PHY(i,j,k)+V_PHY(i  ,j+1,k));
        noah_input_arr(i,j,0,NoahmpInputComp::t_phy)   = getTgivenRandRTh(CONS(i,j,k,Rho_comp),CONS(i,j,k,RhoTheta_comp),qv);
        noah_input_arr(i,j,0,NoahmpInputComp::qv_curr) = qv;
        noah_input_arr(i,j,0,NoahmpInputComp::p8w)     = getPgivenRTh(CONS(i,j,k,RhoTheta_comp),qv);
        noah_input_arr(i,j,0,NoahmpInputComp::swdown)  = SWDOWN(i,j,0);
        noah_input_arr(i,j,0,NoahmpInputComp::glw)     = GLW(i,j,0);
        noah_input_arr(i,j,0,NoahmpInputComp::coszen)  = COSZEN(i,j,0);

        // RAW water-equivalent interval precip [mm]: per slot d = (now - prev) *
        // native_to_kg_m2 (kg/m^2 == mm water). Host applies the guard and derives
        // SR / MP_RAINNC so clamped cells can be reported.
        Real drain = zero, dsnow = zero, dgraup = zero;
        if (hp) {
            Real dd[NoahmpPrecipSlot::NumSlots];
            for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
                dd[s] = slot_present[s]
                      ? amrex::max(zero, (accum_now[s](i,j,kklo) - accum_prv[s](i,j,kklo)) * accum_fac[s])
                      : zero;
            }
            // Hail folded into MP_GRAUP (mass kept in the frozen total); no ERF scheme
            // fills the hail slot today. MP_HAIL is coupled but ERF sets it 0 below --
            // to feed hail separately, stage dd[hail] into MP_HAIL instead.
            dsnow  = dd[NoahmpPrecipSlot::snow];
            dgraup = dd[NoahmpPrecipSlot::graupel] + dd[NoahmpPrecipSlot::hail];
            const Real dfroz = dsnow + dgraup;
            // Total: scheme's `total` slot, else rain + frozen (SAM/Kessler have no total).
            drain = slot_present[NoahmpPrecipSlot::total]
                  ? dd[NoahmpPrecipSlot::total]
                  : (dd[NoahmpPrecipSlot::rain] + dfroz);
        }
        noah_input_arr(i,j,0,NoahmpInputComp::rainbl)   = drain;
        noah_input_arr(i,j,0,NoahmpInputComp::mp_snow)  = dsnow;
        noah_input_arr(i,j,0,NoahmpInputComp::mp_graup) = dgraup;
    });

    // (2) Synchronize to ensure GPU kernel is complete before host access
    Gpu::streamSynchronize();

    // (3) Now on the host, copy the pinned staged data to NoahmpIO arrays. The
    // mechanical 1:1 copies are generated from the registry (3-D members take the
    // k/j transpose (i,1,j); 2-D members take (i,j)); precip is derived below.
    LoopOnCpu(bx, [&] (int i, int j, int ) noexcept
    {
#define NOAHMP_X(name,member) noahmpio->member(i,1,j) = noah_input_arr(i,j,0,NoahmpInputComp::name);
        NOAHMP_INPUT_3D(NOAHMP_X)
#undef NOAHMP_X
#define NOAHMP_X(name,member) noahmpio->member(i,j)   = noah_input_arr(i,j,0,NoahmpInputComp::name);
        NOAHMP_INPUT_2D(NOAHMP_X)
#undef NOAHMP_X
        // Guard the RAW total against non-physical values, then derive SR / MP_RAINNC.
        // Clamped cells are recorded and reported after the box loop (never silent).
        Real dsnow_h  = noah_input_arr(i,j,0,NoahmpInputComp::mp_snow);
        Real dgraup_h = noah_input_arr(i,j,0,NoahmpInputComp::mp_graup);
        Real drain_h  = noah_input_arr(i,j,0,NoahmpInputComp::rainbl);
        if (drain_h > lsm_max_precip_interval) {
            clamped_cells.push_back(noahmp_detail::ClampedPrecipCell{i, j, drain_h});
            drain_h = lsm_max_precip_interval;
        }
        // Enforce the frozen/total invariant AFTER the clamp: MP_SNOW and MP_GRAUP
        // are subsets of MP_RAINNC, so their sum must not exceed it (violated by an
        // inconsistent accumulator, or by the clamp lowering drain_h below the frozen
        // sum). Rescale the frozen components to the total to keep the breakdown
        // self-consistent, and record the cell so a bad source stays diagnosable.
        Real dfroz_h = dsnow_h + dgraup_h;
        const Real inv_tol = Real(1.0e-6) * (Real(1.0) + drain_h);
        if (dfroz_h > drain_h + inv_tol) {
            invariant_cells.push_back(noahmp_detail::InvariantPrecipCell{i, j, dfroz_h, drain_h});
            if (dfroz_h > zero) {
                const Real scale = drain_h / dfroz_h;   // <= 1
                dsnow_h  *= scale;
                dgraup_h *= scale;
                dfroz_h   = dsnow_h + dgraup_h;         // == drain_h (to round-off)
            }
        }
        noahmpio->RAINBL(i,j)    = drain_h;                                       // [mm]
        noahmpio->SR(i,j)        = (drain_h > zero) ? amrex::min(Real(1.0), dfroz_h/drain_h) : zero; // [-]
        noahmpio->MP_RAINNC(i,j) = drain_h;                                       // [mm]
        noahmpio->MP_SNOW(i,j)   = dsnow_h;                                       // [mm]
        noahmpio->MP_GRAUP(i,j)  = dgraup_h;                                      // [mm] (includes folded hail)
        // MP_HAIL is a consumed forcing input; 0 = "no hail" (correct value, not a
        // placeholder). Not lsm_undefined -- a sentinel would corrupt opt_snf=4.
        noahmpio->MP_HAIL(i,j)   = zero;                                          // [mm] ERF-owned; hail folded into MP_GRAUP
    });
}

// ---------------------------------------------------------------------------
//  Steps 5-6: NoahmpIO results -> pinned output buffer -> ERF coupling fields
// ---------------------------------------------------------------------------
void
NOAHMP::read_results (const MFIter& mfi,
                      const int& idb,
                      const Box& bx,
                      const Box& gbx,
                      MultiFab& cons_in)
{
    NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

    // For limiting when populating ghost cells
    int i_lo = bx.smallEnd(0); int i_hi = bx.bigEnd(0);
    int j_lo = bx.smallEnd(1); int j_hi = bx.bigEnd(1);

    const Array4<const Real>& CONS = cons_in.const_array(mfi);

    // Out of NOAH-MP (RRTMGP coupling destinations)
    Array4<Real> TSK           = lsm_fab_data[LsmData_NOAHMP::t_sfc]->array(mfi);
    Array4<Real> EMISS         = lsm_fab_data[LsmData_NOAHMP::sfc_emis]->array(mfi);
    Array4<Real> ALBSFCDIR_VIS = lsm_fab_data[LsmData_NOAHMP::sfc_alb_dir_vis]->array(mfi);
    Array4<Real> ALBSFCDIR_NIR = lsm_fab_data[LsmData_NOAHMP::sfc_alb_dir_nir]->array(mfi);
    Array4<Real> ALBSFCDIF_VIS = lsm_fab_data[LsmData_NOAHMP::sfc_alb_dif_vis]->array(mfi);
    Array4<Real> ALBSFCDIF_NIR = lsm_fab_data[LsmData_NOAHMP::sfc_alb_dif_nir]->array(mfi);
    // Land return-term diagnostic destinations (output only)
    Array4<Real> GRDFLX_o    = lsm_fab_data[LsmData_NOAHMP::grdflx]->array(mfi);
    Array4<Real> FIRA_o      = lsm_fab_data[LsmData_NOAHMP::fira]->array(mfi);
    Array4<Real> SAV_o       = lsm_fab_data[LsmData_NOAHMP::sav]->array(mfi);
    Array4<Real> SAG_o       = lsm_fab_data[LsmData_NOAHMP::sag]->array(mfi);
    Array4<Real> ALBEDO_o    = lsm_fab_data[LsmData_NOAHMP::albedo]->array(mfi);
    Array4<Real> SFCRUNOFF_o = lsm_fab_data[LsmData_NOAHMP::sfcrunoff]->array(mfi);
    Array4<Real> UDRUNOFF_o  = lsm_fab_data[LsmData_NOAHMP::udrunoff]->array(mfi);
    Array4<Real> SMSTAV_o    = lsm_fab_data[LsmData_NOAHMP::smstav]->array(mfi);
    Array4<Real> SMSTOT_o    = lsm_fab_data[LsmData_NOAHMP::smstot]->array(mfi);

    // Per-layer soil output fabs (3 groups x m_nsoil). Gather their Array4s
    // into a device-visible vector so the scatter kernel below can loop over
    // an arbitrary NSOIL. Layout matches soil_data_idx / soil_out_idx:
    // group g in {0:smois,1:sh2o,2:tslb}, layer k in [0,nsoil).
    const int nsoil = m_nsoil;
    const int n_soil_fld = 3*nsoil;
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

    Array4<Real> noah_output_arr = noahmp_output_tmp[idb]->array();

    // (5) Copy results from NoahmpIO back to the pinned output buffer. The
    // mechanical 2-D reads are generated from the registry; banded albedos and
    // the soil profile are read explicitly. SMSTAV/SMSTOT are NOT staged here:
    // this Noah-MP core does not compute them, so they are emitted as
    // lsm_undefined downstream rather than as a misleading physical 0.
    LoopOnCpu(bx, [&] (int i, int j, int ) noexcept
    {
#define NOAHMP_X(name,member) noah_output_arr(i,j,0,NoahmpOutputComp::name) = noahmpio->member(i,j);
        NOAHMP_OUTPUT_2D(NOAHMP_X)
#undef NOAHMP_X
        // Banded radiation reads (VIS=1, NIR=2)
        noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdir_vis) = noahmpio->ALBSFCDIRXY(i,1,j);
        noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdir_nir) = noahmpio->ALBSFCDIRXY(i,2,j);
        noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdif_vis) = noahmpio->ALBSFCDIFXY(i,1,j);
        noah_output_arr(i,j,0,NoahmpOutputComp::albsfcdif_nir) = noahmpio->ALBSFCDIFXY(i,2,j);
        // 3D soil-profile state, one 2D component per soil layer, for any NSOIL.
        // Fortran layer index is 1-based; output components are laid out as
        // [NumComps .. +nsoil) smois, [+nsoil .. +2*nsoil) sh2o, [+2*nsoil ..) tslb.
        for (int k(0); k < nsoil; ++k) {
            noah_output_arr(i,j,0,soil_out_base + 0*nsoil + k) = noahmpio->SMOIS(i,k+1,j);
            noah_output_arr(i,j,0,soil_out_base + 1*nsoil + k) = noahmpio->SH2O (i,k+1,j);
            noah_output_arr(i,j,0,soil_out_base + 2*nsoil + k) = noahmpio->TSLB (i,k+1,j);
        }
    });

    // (6) Copy forcing data from Noahmp to ERF (device; from pinned buffer). This
    // path carries per-field math (flux ÷rho, the -9999 fill guard, sentinels,
    // soil loop) and is deliberately NOT table-driven.
    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Limit indices to the valid box. FillBoundary will pick these up below.
        int ii = std::min(std::max(i,i_lo),i_hi);
        int jj = std::min(std::max(j,j_lo),j_hi);

        // SurfaceLayer fluxes at CC.
        // Noah-MP returns the -9999 fill value for cells it does NOT process
        // (sea-ice / open-water points, which still have LANDMASK=1). Applying
        // that as a flux gives -9999/(rho*Cp) ~ -7.6 K*m/s and crashes the
        // lowest cell to ~200 K. Detect the fill and instead write the
        // lsm_undefined sentinel; the surface layer then falls back to
        // the MOST flux for those cells (see ERF_SurfaceLayer.cpp).
        // NOTE: tau13/tau23 are nodal in xz/yz; the 2D MFs have 1 ghost cell
        //       so the surface layer can average them.
        Real hfx_lsm = noah_output_arr(ii,jj,0,NoahmpOutputComp::hfx);
        if (hfx_lsm > Real(-9990.0)) {
            // Convert Noah-MP surface fluxes to the KINEMATIC form ERF's
            // surface layer stores for its MOST parameter updates (u*, t*, L).
            // The LSM value is interchanged with the MOST value <x'w'> in the
            // same array; ERF's surface-layer/PBL convention is the KINEMATIC
            // MOST form (no rho):
            //   t_flux = <theta' w'>   [K m/s]   (ERF_MOSTStress.H, surf_temp_flux)
            //   q_flux = <qv' w'>      [kg/kg m/s]
            //   tau13/23 = u*^2 comps  [m2/s2]   (SurfaceLayer.cpp: u*=sqrt(tau))
            // Noah-MP returns HFX,LH [W/m2] and TAU_EW/NS [N/m2] (rho INCLUDED),
            // so we divide by rho. rho is dry-air density (~1.14 kg/m3).
            //
            // NOTE on Exner: strictly, ERF's energy variable is POTENTIAL
            // temperature, so the exact conversion of the Noah-MP sensible-heat
            // flux (built from ACTUAL temperatures) to a theta-flux would carry an
            // Exner factor:  theta-flux = HFX/(rho*Cp) * (p0/p)^(Rd/Cp).
            // We deliberately OMIT it here to match WRF, which applies HFX/(rho*Cp)
            // directly into its potential-temperature PBL equation (WRF
            // module_pbl_driver.F: b_t = hfx/rho/CP, no Exner). This is an
            // approximation relative to ERF's exact potential-temperature flux, not
            // an identity: the error scales with (p0/p)^(Rd/Cp)-1, which is small
            // near sea level (p ~= p0 -> Exner ~= 1.00-1.01, <~1%) but grows over
            // low-pressure / elevated terrain (e.g. p ~ 800 hPa -> ~6-7%). Reinstate
            // by multiplying t_flux by std::pow(p_0/getPgivenRTh(
            // CONS(ii,jj,k,RhoTheta_comp),
            // (is_moist)?CONS(ii,jj,k,RhoQ1_comp)/CONS(ii,jj,k,Rho_comp):zero),
            // R_d/Cp_d) if that bias matters for the case.
            Real rho_l  = CONS(ii,jj,k,Rho_comp);
            t_flux_arr(i,j,k) = hfx_lsm/(rho_l*Cp_d);
            q_flux_arr(i,j,k) = noah_output_arr(ii,jj,0,NoahmpOutputComp::lh)/(rho_l*L_v);
            tau13_arr(i,j,k)  = noah_output_arr(ii,jj,0,NoahmpOutputComp::tau_ew)/rho_l;
            tau23_arr(i,j,k)  = noah_output_arr(ii,jj,0,NoahmpOutputComp::tau_ns)/rho_l;
        } else {
            t_flux_arr(i,j,k) = lsm_undefined;
            q_flux_arr(i,j,k) = lsm_undefined;
            tau13_arr(i,j,k)  = lsm_undefined;
            tau23_arr(i,j,k)  = lsm_undefined;
        }

        // RRTMGP variables
        if (hfx_lsm > Real(-9990.0)) {
            TSK(i,j,0)           = noah_output_arr(ii,jj,0,NoahmpOutputComp::tsk);
            EMISS(i,j,0)         = noah_output_arr(ii,jj,0,NoahmpOutputComp::emiss);
            ALBSFCDIR_VIS(i,j,0) = noah_output_arr(ii,jj,0,NoahmpOutputComp::albsfcdir_vis);
            ALBSFCDIR_NIR(i,j,0) = noah_output_arr(ii,jj,0,NoahmpOutputComp::albsfcdir_nir);
            ALBSFCDIF_VIS(i,j,0) = noah_output_arr(ii,jj,0,NoahmpOutputComp::albsfcdif_vis);
            ALBSFCDIF_NIR(i,j,0) = noah_output_arr(ii,jj,0,NoahmpOutputComp::albsfcdif_nir);
            // Land return-term diagnostics (output only; not fed to dynamics)
            GRDFLX_o(i,j,0)      = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_grdflx);
            FIRA_o(i,j,0)        = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_fira);
            SAV_o(i,j,0)         = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_sav);
            SAG_o(i,j,0)         = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_sag);
            ALBEDO_o(i,j,0)      = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_albedo);
            SFCRUNOFF_o(i,j,0)   = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_sfcrunoff);
            UDRUNOFF_o(i,j,0)    = noah_output_arr(ii,jj,0,NoahmpOutputComp::o_udrunoff);
            // SMSTAV/SMSTOT are not computed by this Noah-MP core, so emit the
            // lsm_undefined sentinel rather than the uncomputed core value: a
            // physical-looking 0 would misrepresent them as valid diagnostics.
            SMSTAV_o(i,j,0)      = lsm_undefined;
            SMSTOT_o(i,j,0)      = lsm_undefined;
            // Per-layer soil profile (any NSOIL): copy each output component
            // into its dedicated 2D fab (soil_arr laid out group-major).
            for (int s(0); s < n_soil_fld; ++s) {
                soil_arr[s](i,j,0) = noah_output_arr(ii,jj,0,soil_out_base + s);
            }
        } else {
            TSK(i,j,0)           = lsm_undefined;
            EMISS(i,j,0)         = lsm_undefined;
            ALBSFCDIR_VIS(i,j,0) = lsm_undefined;
            ALBSFCDIR_NIR(i,j,0) = lsm_undefined;
            ALBSFCDIF_VIS(i,j,0) = lsm_undefined;
            ALBSFCDIF_NIR(i,j,0) = lsm_undefined;
            GRDFLX_o(i,j,0)      = lsm_undefined;
            FIRA_o(i,j,0)        = lsm_undefined;
            SAV_o(i,j,0)         = lsm_undefined;
            SAG_o(i,j,0)         = lsm_undefined;
            ALBEDO_o(i,j,0)      = lsm_undefined;
            SFCRUNOFF_o(i,j,0)   = lsm_undefined;
            UDRUNOFF_o(i,j,0)    = lsm_undefined;
            SMSTAV_o(i,j,0)      = lsm_undefined;
            SMSTOT_o(i,j,0)      = lsm_undefined;
            // Per-layer soil profile: undefined sentinel on non-land cells.
            for (int s(0); s < n_soil_fld; ++s) {
                soil_arr[s](i,j,0) = lsm_undefined;
            }
        }
    });
    // The scatter kernel reads soil_arr_d (a per-iteration DeviceVector);
    // synchronize before it is destroyed at the end of this method.
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
        // NOTE: Do not try to interpolate if lev 0 was just updated. Since Noah is
        //       called post-step the fluxes & data will contain lsm_undefined values.
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

        // -------------------------------------------------------------------------
        // Precipitation forcing into Noah-MP (RAINBL / MP_RAINNC / MP_SNOW / MP_GRAUP /
        // SR). Works for ANY microphysics scheme through the typed source interface;
        // see ERF_NOAHMP_Precip.cpp. Collect the slots, then lazily allocate/seed the
        // per-level previous-accumulation snapshots the delta kernel differences against.
        const noahmp_detail::PrecipSlots precip = collect_precip_sources(precip_sources);
        prepare_precip_snapshots(lev, precip);

        // Cells the guards touched this call (per rank); reported after the box loop.
        Vector<noahmp_detail::ClampedPrecipCell>   clamped_cells;
        Vector<noahmp_detail::InvariantPrecipCell> invariant_cells;

        // Loop over blocks to copy forcing data to Noahmp, drive the land model,
        // and copy data back to ERF Multifabs.
        int idb = 0;
        for (MFIter mfi(cons_in); mfi.isValid(); ++mfi, ++idb) {

            Box bx  = mfi.tilebox();
            Box gbx = mfi.tilebox(IntVect(0,0,0),IntVect(1,1,0));

            // Check if tile is at the lower boundary in lower z direction
            if (bx.smallEnd(2) != klo) { continue; }

            bx.makeSlab(2,klo);
            gbx.makeSlab(2,klo);

            NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

            // (1-3) ERF forcing -> pinned input -> NoahmpIO arrays.
            stage_forcing(mfi, idb, bx, klo, lev, is_moist,
                          cons_in, xvel_in, yvel_in,
                          precip, clamped_cells, invariant_cells);

            // (4) Call the noahmpio driver code. Mirror the authoritative counter
            // into each block just before firing.
            noahmpio->itimestep = m_itimestep;
            noahmpio->DriverMain();

            // The soil-layer count was fixed in Init (erf.lsm_nsoil, asserted to match
            // the namelist NSOIL there), so the per-layer reads are always in range.

            // (5-6) NoahmpIO results -> pinned output -> ERF coupling fields.
            read_results(mfi, idb, bx, gbx, cons_in);
        }

        // Advance the precip snapshots to the current accumulation now that this
        // land call's RAINBL/SR have been consumed.
        advance_precip_snapshots(lev, precip);

        // Report any cells the guards clamped / rescaled (per rank; never silent).
        report_precip_diagnostics(lev, clamped_cells, invariant_cells);

        // Fill the ghost cells
        for (auto ivar = 0; ivar < LsmFlux_NOAHMP::NumVars; ++ivar) {
            lsm_fab_flux[ivar]->FillBoundary(m_geom.periodicity());
        }
        Print () << "Noah-MP driver at level " << lev << " completed" << std::endl;
    } // lev == 0
}
