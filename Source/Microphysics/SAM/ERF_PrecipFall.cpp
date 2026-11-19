#include "ERF_NumericalConstants.H"
#include "ERF_MicrophysicsConstants.H"
#include "ERF_SAM.H"
#include "ERF_SAMUtils.H"
#include "ERF_TileNoZ.H"
#include <cmath>

using namespace amrex;

/**
 * Precipitation fluxes P_{r/s/g} (A19)
 *
 * Code modified from SAMXX, the C++ version of the SAM code.
 *
 * @param[in] hydro_type Type selection for the precipitation advection hydrodynamics scheme (0-3)
 */
void
SAM::PrecipFall (const SolverChoice& sc)
{
    if (sam_is_no_precip(sc.moisture_type)) return;

    // Conservative component-sedimentation contract: rain, snow, and graupel
    // sediment with separate component fluxes. Each face flux is limited by
    // the donor cell's available component mass, including detJ when present.
    // The same limited bottom fluxes update surface accumulation and the
    // in-column qpr/qps/qpg flux-divergence update. This closes component and
    // total precipitation column budgets up to roundoff. Sedimentation does
    // not update theta.

    Real rho_0 = Real(1.29);

    Real gamr3 = erf_gammafff(Real(4.0)+b_rain);
    Real gams3 = erf_gammafff(Real(4.0)+b_snow);
    Real gamg3 = erf_gammafff(Real(4.0)+b_grau);

    Real vrain = (a_rain*gamr3/Real(6.0))*std::pow((PI*rhor*nzeror),-crain);
    Real vsnow = (a_snow*gams3/Real(6.0))*std::pow((PI*rhos*nzeros),-csnow);
    Real vgrau = (a_grau*gamg3/Real(6.0))*std::pow((PI*rhog*nzerog),-cgrau);

    Real dtn  = dt;
    // NOTE: coef carries only the reference vertical spacing. The physical
    //       spacing is supplied by the inverse Jacobian in the tendency below,
    //       since detJ = dz_phys/CellSize(2). Dividing by dz_phys here as well
    //       would apply the vertical metric twice.
    Real coef = dtn * m_geom.InvCellSize(2);

    auto domain = m_geom.Domain();
    int k_lo = domain.smallEnd(2);
    int k_hi = domain.bigEnd(2);

    auto qpr   = mic_fab_vars[MicVar::qpr];
    auto qps   = mic_fab_vars[MicVar::qps];
    auto qpg   = mic_fab_vars[MicVar::qpg];
    auto qp    = mic_fab_vars[MicVar::qp];
    auto rho   = mic_fab_vars[MicVar::rho];
    auto tabs  = mic_fab_vars[MicVar::tabs];
    auto rain_accum = mic_fab_vars[MicVar::rain_accum];
    auto snow_accum = mic_fab_vars[MicVar::snow_accum];
    auto graupel_accum = mic_fab_vars[MicVar::graupel_accum];
    auto precip_accum = mic_fab_vars[MicVar::precip_accum];

    auto ba    = tabs->boxArray();
    auto dm    = tabs->DistributionMap();
    auto ngrow = tabs->nGrowVect();

    MultiFab fz;
    fz.define(convert(ba, IntVect(0,0,1)), dm, 3, ngrow);

    // The face states read tabs, rho and the precipitating species at k-1, a ghost cell
    // whenever a box boundary is interior in z. Cloud and Precip have just updated them in
    // valid cells only, so refresh the ghosts before the fall speeds that size the substeps
    // are computed; the substep loop below refreshes the species again as they fall. tabs
    // does not change here, so once is enough for it.
    tabs->FillBoundary(m_geom.periodicity());
    rho->FillBoundary(m_geom.periodicity());
    qpr->FillBoundary(m_geom.periodicity());
    qps->FillBoundary(m_geom.periodicity());
    qpg->FillBoundary(m_geom.periodicity());

    // Size the substeps from the largest fall Courant number over the faces (#5168): the
    // fastest of rain, snow and graupel at each face, from the same face state as the fluxes
    // below, over the smaller of the heights of the cell that face drains and the cell it
    // drains into (sam_fall_sizing_detJ says why the cell below counts too).
    int n_substep;
    {
        auto const& ma_rho  = rho->const_arrays();
        auto const& ma_tabs = tabs->const_arrays();
        auto const& ma_qpr  = qpr->const_arrays();
        auto const& ma_qps  = qps->const_arrays();
        auto const& ma_qpg  = qpg->const_arrays();
        // On an EB mesh detJ is the volume fraction, not a cell height, so the sizing uses the
        // reference dz there (J = 1).  The donor cap and the tendency below still use detJ, the
        // volume fraction, so a small cut cell can be donor-capped rather than substepped.
        // Off EB, detJ_cc has one ghost cell filled (make_J), so the receiver cell k-1 at the
        // bottom of a box that is interior in z can be read here.
        const auto ma_dJ = (m_detJ_cc && sc.terrain_type != TerrainType::EB)
                         ? m_detJ_cc->const_arrays() : MultiArray4<const Real>{};
        const Real dz_ref = m_geom.CellSize(2);
        GpuTuple<Real> max = ParReduce(TypeList<ReduceOpMax>{},
                                       TypeList<Real>{},
                                       fz, IntVect(0),
                             [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                             -> GpuTuple<Real>
                             {
                                 const SAMPrecipComponentFaceState face_state =
                                     sam_precip_component_face_state(ma_rho[box_no], ma_tabs[box_no],
                                                                     ma_qpr[box_no], ma_qps[box_no],
                                                                     ma_qpg[box_no],
                                                                     i, j, k, k_lo, k_hi);
                                 const int donor_k = sam_precip_face_donor_k(k, k_lo, k_hi);
                                 const int receiver_k = sam_precip_face_receiver_k(k, k_lo, k_hi);
                                 const Real detJ_size = (ma_dJ)
                                     ? sam_fall_sizing_detJ(ma_dJ[box_no](i,j,donor_k),
                                                            ma_dJ[box_no](i,j,receiver_k))
                                     : one;
                                 return { sam_fall_courant(sam_precip_component_max_fall_speed(face_state,
                                                                                               vrain, vsnow, vgrau,
                                                                                               rho_0),
                                                           dtn, dz_ref, detJ_size) };
                             });
        Real max_courant = get<0>(max);
        // ParReduce over a MultiFab gives the local rank maximum here. PrecipFall
        // needs one global substep count across ranks, including ranks that own no
        // active tiles for this reduction, so keep the explicit MPI max reduction.
        ParallelDescriptor::ReduceRealMax(max_courant);
        sam_check_sedimentation_substeps("SAM::PrecipFall", max_courant);
        n_substep = sam_sedimentation_substep_count(max_courant);
    }
    AMREX_ALWAYS_ASSERT(n_substep >= 1);
    coef /= Real(n_substep);
    dtn  /= Real(n_substep);

    // Substep the vertical advection
    for (int nsub(0); nsub<n_substep; ++nsub)
    {
        // Refresh ghosts so both sides of a shared z face use the same donor values
        qpr->FillBoundary(m_geom.periodicity());
        qps->FillBoundary(m_geom.periodicity());
        qpg->FillBoundary(m_geom.periodicity());
        mic_fab_vars[MicVar::rho]->FillBoundary(m_geom.periodicity());

        for (MFIter mfi(*qp, TileNoZ()); mfi.isValid(); ++mfi) {
            auto qpr_array    = qpr->array(mfi);
            auto qps_array    = qps->array(mfi);
            auto qpg_array    = qpg->array(mfi);
            auto qp_array     = qp->array(mfi);
            auto rho_array    = rho->const_array(mfi);
            auto tabs_array   = tabs->const_array(mfi);
            auto fz_array     = fz.array(mfi);

            auto rain_accum_array = rain_accum->array(mfi);
            auto snow_accum_array = snow_accum->array(mfi);
            auto graupel_accum_array = graupel_accum->array(mfi);
            auto precip_accum_array = precip_accum->array(mfi);
            auto sed_array = m_sed.down->array(mfi);
            const auto inflow_array = (m_sed_inflow) ? m_sed_inflow->const_array(mfi)
                                                     : Array4<const Real>{};
            const int inflow_flag = m_sed.inflow_flag_comp();
            const Real inv_dt = one / dt;

            const auto dJ_array = (m_detJ_cc) ? m_detJ_cc->const_array(mfi) : Array4<const Real>{};

            const auto& tbx = mfi.tilebox();
            const auto& tbz = mfi.tilebox(IntVect(0,0,1),IntVect(0));

            // Update vertical flux every substep
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const SAMPrecipComponentFaceState face_state =
                    sam_precip_component_face_state(rho_array, tabs_array,
                                                    qpr_array, qps_array, qpg_array,
                                                    i, j, k, k_lo, k_hi);
                const SAMPrecipFluxComponents raw_fluxes =
                    sam_precip_component_fluxes_from_face_state(face_state,
                                                                vrain, vsnow, vgrau);
                const SAMPrecipFluxComponents corrected_fluxes =
                    sam_precip_flux_components_density_corrected(raw_fluxes,
                                                                 rho_0,
                                                                 face_state.rho_avg);
                const int donor_k = sam_precip_face_donor_k(k, k_lo, k_hi);
                const Real detJ_donor = (dJ_array) ? dJ_array(i,j,donor_k) : one;

                // NOTE: Fz is the sedimentation flux from the advective operator.
                //       In the terrain-following coordinate system, the z-deriv in
                //       the divergence uses the normal velocity (Omega). However,
                //       there are no u/v components to the sedimentation velocity.
                //       Therefore, we simply end up with a division by detJ when
                //       evaluating the source term: dJinv * (flux_hi - flux_lo) * dzinv.
                fz_array(i,j,k,0) = sam_limit_precip_component_flux(corrected_fluxes.rain,
                                                                    rho_array(i,j,donor_k),
                                                                    qpr_array(i,j,donor_k),
                                                                    detJ_donor,
                                                                    coef);
                fz_array(i,j,k,1) = sam_limit_precip_component_flux(corrected_fluxes.snow,
                                                                    rho_array(i,j,donor_k),
                                                                    qps_array(i,j,donor_k),
                                                                    detJ_donor,
                                                                    coef);
                fz_array(i,j,k,2) = sam_limit_precip_component_flux(corrected_fluxes.graupel,
                                                                    rho_array(i,j,donor_k),
                                                                    qpg_array(i,j,donor_k),
                                                                    detJ_donor,
                                                                    coef);

                // At the top of a fine patch the precipitation coming in is the coarse
                // level's flux through that face, not one built from interpolated ghost cells
                if (inflow_array && inflow_array(i,j,k,inflow_flag) > Real(0.5)) {
                    fz_array(i,j,k,0) = inflow_array(i,j,k,0) * inv_dt;
                    fz_array(i,j,k,1) = inflow_array(i,j,k,1) * inv_dt;
                    fz_array(i,j,k,2) = inflow_array(i,j,k,2) * inv_dt;
                }

                sed_array(i,j,k,0) += fz_array(i,j,k,0) * dtn;
                sed_array(i,j,k,1) += fz_array(i,j,k,1) * dtn;
                sed_array(i,j,k,2) += fz_array(i,j,k,2) * dtn;

                if (k == k_lo) {
                    const SAMSurfaceAccumulation surface_accum =
                        sam_surface_accumulation_from_component_fluxes(
                            {fz_array(i,j,k,0), fz_array(i,j,k,1), fz_array(i,j,k,2)},
                            dtn);
                    rain_accum_array(i,j,k)  = rain_accum_array(i,j,k)  + surface_accum.rain;
                    snow_accum_array(i,j,k)  = snow_accum_array(i,j,k)  + surface_accum.snow;
                    graupel_accum_array(i,j,k) = graupel_accum_array(i,j,k) + surface_accum.graupel;
                    // The total is in mm of liquid water, whatever the species' own densities
                    precip_accum_array(i,j,k) += (fz_array(i,j,k,0) + fz_array(i,j,k,1) + fz_array(i,j,k,2))
                                               * dtn / rhor * Real(1000.0);
                }
            });

            // Update precip every substep
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Jacobian determinant
                Real dJinv = (dJ_array) ? one/dJ_array(i,j,k) : one;

                //==================================================
                // Precipitating sedimentation (A19)
                //==================================================
                const Real dqpr = sam_sedimentation_tendency(fz_array(i,j,k+1,0), fz_array(i,j,k,0),
                                                             rho_array(i,j,k), dJinv, coef);
                const Real dqps = sam_sedimentation_tendency(fz_array(i,j,k+1,1), fz_array(i,j,k,1),
                                                             rho_array(i,j,k), dJinv, coef);
                const Real dqpg = sam_sedimentation_tendency(fz_array(i,j,k+1,2), fz_array(i,j,k,2),
                                                             rho_array(i,j,k), dJinv, coef);

                qpr_array(i,j,k) = std::max(Real(0.0), qpr_array(i,j,k) + dqpr);
                qps_array(i,j,k) = std::max(Real(0.0), qps_array(i,j,k) + dqps);
                qpg_array(i,j,k) = std::max(Real(0.0), qpg_array(i,j,k) + dqpg);
                 qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                 // NOTE: Sedimentation does not affect the potential temperature,
                 //       but it does affect the liquid/ice static energy.
                 //       No source to Theta occurs here.
            });
        } // mfi
    } // nsub
}
