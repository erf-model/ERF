#include "ERF_Constants.H"
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

    Real rho_0 = Real(1.29);

    Real gamr3 = erf_gammafff(Real(4.0)+b_rain);
    Real gams3 = erf_gammafff(Real(4.0)+b_snow);
    Real gamg3 = erf_gammafff(Real(4.0)+b_grau);

    Real vrain = (a_rain*gamr3/Real(6.0))*std::pow((PI*rhor*nzeror),-crain);
    Real vsnow = (a_snow*gams3/Real(6.0))*std::pow((PI*rhos*nzeros),-csnow);
    Real vgrau = (a_grau*gamg3/Real(6.0))*std::pow((PI*rhog*nzerog),-cgrau);

    Real dtn  = dt;
    Real coef = dtn/m_dzmin;

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
    auto graup_accum = mic_fab_vars[MicVar::graup_accum];

    auto ba    = tabs->boxArray();
    auto dm    = tabs->DistributionMap();
    auto ngrow = tabs->nGrowVect();

    MultiFab fz;
    fz.define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow);

    int SAM_moisture_type = 1;
    if (sc.moisture_type == MoistureType::SAM_NoIce) {
        SAM_moisture_type = 2;
    }

    //  Precompute the vertical fluxes for CFL constraint
    for (MFIter mfi(fz, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto qp_array   = qp->const_array(mfi);
        auto rho_array  = rho->const_array(mfi);
        auto tabs_array = tabs->const_array(mfi);
        auto fz_array   = fz.array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const SAMPrecipFaceState face_state =
                sam_precip_face_state(SAM_moisture_type,
                                      rho_array, tabs_array, qp_array,
                                      i, j, k, k_lo, k_hi);
            const Real Pprecip = sam_precip_flux_from_face_state(face_state,
                                                                 vrain, vsnow, vgrau);

            // NOTE: Fz is the sedimentation flux from the advective operator.
            //       In the terrain-following coordinate system, the z-deriv in
            //       the divergence uses the normal velocity (Omega). However,
            //       there are no u/v components to the sedimentation velocity.
            //       Therefore, we simply end up with a division by detJ when
            //       evaluating the source term: dJinv * (flux_hi - flux_lo) * dzinv.
            fz_array(i,j,k) = sam_precip_flux_density_corrected(Pprecip,
                                                                rho_0,
                                                                face_state.rho_avg);
        });
    }

    // Compute number of substeps from maximum terminal velocity
    Real wt_max;
    int n_substep;
    auto const& ma_fz_arr = fz.const_arrays();
    GpuTuple<Real> max = ParReduce(TypeList<ReduceOpMax>{},
                                   TypeList<Real>{},
                                   fz, IntVect(0),
                         [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                         -> GpuTuple<Real>
                         {
                             return { ma_fz_arr[box_no](i,j,k) };
                         });
    wt_max = get<0>(max) + std::numeric_limits<Real>::epsilon();
    n_substep = sam_substep_count_from_reduced_flux(wt_max, dtn, m_dzmin);
    AMREX_ALWAYS_ASSERT(n_substep >= 1);
    coef /= Real(n_substep);
    dtn  /= Real(n_substep);

    // Substep the vertical advection
    for (int nsub(0); nsub<n_substep; ++nsub) {
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
            auto graup_accum_array = graup_accum->array(mfi);

            const auto dJ_array = (m_detJ_cc) ? m_detJ_cc->const_array(mfi) : Array4<const Real>{};

            const auto& tbx = mfi.tilebox();
            const auto& tbz = mfi.tilebox(IntVect(0,0,1),IntVect(0));

            // Update vertical flux every substep
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const SAMPrecipFaceState face_state =
                    sam_precip_face_state(SAM_moisture_type,
                                          rho_array, tabs_array, qp_array,
                                          i, j, k, k_lo, k_hi);
                const Real Pprecip = sam_precip_flux_from_face_state(face_state,
                                                                     vrain, vsnow, vgrau);

                // NOTE: Fz is the sedimentation flux from the advective operator.
                //       In the terrain-following coordinate system, the z-deriv in
                //       the divergence uses the normal velocity (Omega). However,
                //       there are no u/v components to the sedimentation velocity.
                //       Therefore, we simply end up with a division by detJ when
                //       evaluating the source term: dJinv * (flux_hi - flux_lo) * dzinv.
                fz_array(i,j,k) = sam_precip_flux_density_corrected(Pprecip,
                                                                    rho_0,
                                                                    face_state.rho_avg);

                if(k==k_lo){
                    const SAMSurfaceAccumulation surface_accum =
                        sam_surface_accumulation(face_state.rho_avg, face_state.qp_avg,
                                                 face_state.omp, face_state.omg,
                                                 vrain, vsnow, vgrau, dtn);
                    rain_accum_array(i,j,k)  = rain_accum_array(i,j,k)  + surface_accum.rain;
                    snow_accum_array(i,j,k)  = snow_accum_array(i,j,k)  + surface_accum.snow;
                    graup_accum_array(i,j,k) = graup_accum_array(i,j,k) + surface_accum.graupel;
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
                Real dqp = sam_sedimentation_tendency(fz_array(i,j,k+1), fz_array(i,j,k),
                                                      rho_array(i,j,k), dJinv, coef);
                const Real omp = sam_precip_rain_fraction(SAM_moisture_type,
                                                          tabs_array(i,j,k));
                const Real omg = sam_graupel_fraction(SAM_moisture_type,
                                                      tabs_array(i,j,k));

                qpr_array(i,j,k) = std::max(Real(0), qpr_array(i,j,k) + dqp*omp);
                qps_array(i,j,k) = std::max(Real(0), qps_array(i,j,k) + dqp*(one-omp)*(one-omg));
                qpg_array(i,j,k) = std::max(Real(0), qpg_array(i,j,k) + dqp*(one-omp)*omg);
                 qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                 // NOTE: Sedimentation does not affect the potential temperature,
                 //       but it does affect the liquid/ice static energy.
                 //       No source to Theta occurs here.
            });
        } // mfi
    } // nsub
}

