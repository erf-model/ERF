#include "ERF_Constants.H"
#include "SAM.H"
#include "TileNoZ.H"

using namespace amrex;

/**
 * Precipitation fluxes P_{r/s/g} (A19)
 *
 * Code modified from SAMXX, the C++ version of the SAM code.
 *
 * @param[in] hydro_type Type selection for the precipitation advection hydrodynamics scheme (0-3)
 */
void SAM::PrecipFall (int hydro_type)
{
    Real eps = std::numeric_limits<Real>::epsilon();
    bool constexpr nonos = true;

    Real gamr3 = erf_gammafff(4.0+b_rain);
    Real gams3 = erf_gammafff(4.0+b_snow);
    Real gamg3 = erf_gammafff(4.0+b_grau);

    Real vrain = (a_rain*gamr3/6.0)*pow((PI*rhor*nzeror),-crain);
    Real vsnow = (a_snow*gams3/6.0)*pow((PI*rhos*nzeros),-csnow);
    Real vgrau = (a_grau*gamg3/6.0)*pow((PI*rhog*nzerog),-cgrau);

    Real dt_advance = dt;
    int nz = nlev;

    auto qpr   = mic_fab_vars[MicVar::qpr];
    auto qps   = mic_fab_vars[MicVar::qps];
    auto qpg   = mic_fab_vars[MicVar::qpg];
    auto qp    = mic_fab_vars[MicVar::qp];
    auto omega = mic_fab_vars[MicVar::omega];
    auto rho   = mic_fab_vars[MicVar::rho];
    auto tabs  = mic_fab_vars[MicVar::tabs];
    auto theta = mic_fab_vars[MicVar::theta];

    auto ba    = tabs->boxArray();
    auto dm    = tabs->DistributionMap();
    auto ngrow = tabs->nGrowVect();

    MultiFab mx;
    MultiFab mn;
    MultiFab lfac;
    MultiFab www;
    MultiFab fz;
    MultiFab wp;
    MultiFab tmp_qp;
    MultiFab prec_cfl_fab;

    mx.define(ba,dm, 1, ngrow);
    mn.define(ba,dm, 1, ngrow);
    lfac.define(ba, dm, 1, ngrow);
    www.define(ba, dm, 1, ngrow);
    fz.define(ba, dm, 1, ngrow);
    wp.define(ba, dm, 1, ngrow);
    tmp_qp.define(ba, dm, 1, ngrow);
    prec_cfl_fab.define(ba, dm, 1, ngrow);

    TableData<Real, 1> irho;
    TableData<Real, 1> iwmax;
    TableData<Real, 1> rhofac;

    irho.resize({zlo},{zhi});
    iwmax.resize({zlo},{zhi});
    rhofac.resize({zlo},{zhi});

    auto irho_t   = irho.table();
    auto iwmax_t  = iwmax.table();
    auto rhofac_t = rhofac.table();
    auto rho1d_t  = rho1d.table();

    auto dz = m_geom.CellSize(2);

    // Precompute omega variable
    for ( MFIter mfi(*(mic_fab_vars[MicVar::omega]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto omega_array = mic_fab_vars[MicVar::omega]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar::tabs]->array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            omega_array(i,j,k) = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
        });
    }

    ParallelFor(nz, [=] AMREX_GPU_DEVICE (int k) noexcept
    {
        rhofac_t(k)  = std::sqrt(1.29/rho1d_t(k));
        irho_t(k)    = 1.0/rho1d_t(k);
        Real wmax    = dz/dt_advance;   // Velocity equivalent to a cfl of 1.0.
        iwmax_t(k)   = 1.0/wmax;
    });

    //  Add sedimentation of precipitation field to the vert. vel.
    for ( MFIter mfi(lfac, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto lfac_array     = lfac.array(mfi);
        auto omega_array    = omega->array(mfi);
        auto qp_array       = qp->array(mfi);
        auto rho_array      = rho->array(mfi);
        auto tabs_array     = tabs->array(mfi);
        auto wp_array       = wp.array(mfi);
        auto www_array      = www.array(mfi);
        auto fz_array       = fz.array(mfi);
        auto prec_cfl_array = prec_cfl_fab.array(mfi);

        const auto& box3d = mfi.tilebox();

        const Real fac_cond = m_fac_cond;
        const Real fac_sub  = m_fac_sub;
        const Real fac_fus  = m_fac_fus;

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (hydro_type == 0) {
                lfac_array(i,j,k) = fac_cond;
            }
            else if (hydro_type == 1) {
                lfac_array(i,j,k) = fac_sub;
            }
            else if (hydro_type == 2) {
                lfac_array(i,j,k) = fac_cond + (1.0-omega_array(i,j,k))*fac_fus;
            }
            else if (hydro_type == 3) {
                lfac_array(i,j,k) = 0.0;
            }
            Real Veff = term_vel_qp(qp_array(i,j,k), vrain, vsnow, vgrau,
                                    rho_array(i,j,k), tabs_array(i,j,k));
            wp_array(i,j,k) = -Veff * std::sqrt(1.29/rho_array(i,j,k));
            prec_cfl_array(i,j,k) = wp_array(i,j,k) * iwmax_t(k);
            wp_array(i,j,k) *= rho_array(i,j,k) * dt_advance/dz;
            if (k == 0) {
                fz_array(i,j,nz-1)   = 0.0;
                www_array(i,j,nz-1)  = 0.0;
                lfac_array(i,j,nz-1) = 0.0;
            }
        });
    }

    auto const& cfl_arrays = prec_cfl_fab.const_arrays();
    Real prec_cfl = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                              prec_cfl_fab, IntVect(0),
                              [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                              ->GpuTuple<Real>
    {
        return { cfl_arrays[box_no](i,j,k) };
    });

    // If maximum CFL due to precipitation velocity is greater than 0.9,
    // take more than one advection step to maintain stability.
    int nprec;
    if (prec_cfl > 0.9) {
        nprec = static_cast<int>(std::ceil(prec_cfl/0.9));
        for (MFIter mfi(wp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto wp_array = wp.array(mfi);
            const auto& box3d = mfi.tilebox();

            ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int k, int j, int i)
            {
                // wp already includes factor of dt, so reduce it by a
                // factor equal to the number of precipitation steps.
                wp_array(i,j,k) = wp_array(i,j,k)/Real(nprec);
            });
        }
    } else {
        nprec = 1;
    }

#ifdef ERF_FIXED_SUBCYCLE
    nprec = 4;
#endif

    for(int iprec = 1; iprec<=nprec; iprec++) {
        for ( MFIter mfi(tmp_qp, TileNoZ()); mfi.isValid(); ++mfi) {
            auto qpr_array    = qpr->array(mfi);
            auto qps_array    = qps->array(mfi);
            auto qpg_array    = qpg->array(mfi);
            auto qp_array     = qp->array(mfi);
            auto rho_array    = rho->array(mfi);
            auto tabs_array   = tabs->array(mfi);
            auto theta_array  = theta->array(mfi);
            auto tmp_qp_array = tmp_qp.array(mfi);
            auto mx_array     = mx.array(mfi);
            auto mn_array     = mn.array(mfi);
            auto fz_array     = fz.array(mfi);
            auto wp_array     = wp.array(mfi);
            auto lfac_array   = lfac.array(mfi);
            auto www_array    = www.array(mfi);

            const auto& box3d = mfi.tilebox();

            ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                tmp_qp_array(i,j,k) = qp_array(i,j,k); // Temporary array for qp in this column
            });

            ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (nonos) {
                    int kc=min(nz-1,k+1);
                    int kb=max(0,k-1);
                    mx_array(i,j,k) = max(tmp_qp_array(i,j,kb), max(tmp_qp_array(i,j,kc), tmp_qp_array(i,j,k)));
                    mn_array(i,j,k) = min(tmp_qp_array(i,j,kb), min(tmp_qp_array(i,j,kc), tmp_qp_array(i,j,k)));
                }
                // Define upwind precipitation flux
                fz_array(i,j,k) = tmp_qp_array(i,j,k)*wp_array(i,j,k);
            });

            ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int kc = min(k+1, nz-1);

                Real dqp = (fz_array(i,j,kc)-fz_array(i,j,k)) / rho_array(i,j,k);
                tmp_qp_array(i,j,k) = tmp_qp_array(i,j,k) + (fz_array(i,j,kc)-fz_array(i,j,k)) / rho_array(i,j,k); //Update temporary qp
            });

            ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Also, compute anti-diffusive correction to previous
                // (upwind) approximation to the flux
                int kb=max(0,k-1);
                // The precipitation velocity is a cell-centered quantity,
                // since it is computed from the cell-centered
                // precipitation mass fraction.  Therefore, a reformulated
                // anti-diffusive flux is used here which accounts for
                // this and results in reduced numerical diffusion.
                www_array(i,j,k) = 0.5*(1.0+wp_array(i,j,k)/rho_array(i,j,k))*(tmp_qp_array(i,j,kb)*wp_array(i,j,kb) -
                                                                        tmp_qp_array(i,j,k)*wp_array(i,j,k)); // works for wp(k)<0
            });

            if (nonos) {
                ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    int kc=min(nz-1,k+1);
                    int kb=max(0,k-1);
                    mx_array(i,j,k) = max(tmp_qp_array(i,j,kb),max(tmp_qp_array(i,j,kc), max(tmp_qp_array(i,j,k), mx_array(i,j,k))));
                    mn_array(i,j,k) = min(tmp_qp_array(i,j,kb),min(tmp_qp_array(i,j,kc), min(tmp_qp_array(i,j,k), mn_array(i,j,k))));
                    kc = min(nz-1,k+1);
                    mx_array(i,j,k) = rho_array(i,j,k)*(mx_array(i,j,k)-tmp_qp_array(i,j,k))/(pn(www_array(i,j,kc)) +
                                                                                        pp(www_array(i,j,k))+eps);
                    mn_array(i,j,k) = rho_array(i,j,k)*(tmp_qp_array(i,j,k)-mn_array(i,j,k))/(pp(www_array(i,j,kc)) +
                                                                                        pn(www_array(i,j,k))+eps);
                });

                ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    int kb=max(0,k-1);
                    // Add limited flux correction to fz(k).
                    fz_array(i,j,k) = fz_array(i,j,k)
                                    + pp(www_array(i,j,k))*std::min(1.0,std::min(mx_array(i,j,k), mn_array(i,j,kb)))
                                    - pn(www_array(i,j,k))*std::min(1.0,std::min(mx_array(i,j,kb),mn_array(i,j,k))); // Anti-diffusive flux
                });
            }

            // Update precipitation mass fraction and liquid-ice static
            // energy using precipitation fluxes computed in this column.
            ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int kc=min(k+1, nz-1);
                // Update precipitation mass fraction.
                // Note that fz is the total flux, including both the
                // upwind flux and the anti-diffusive correction.
                //      Real flagstat = 1.0;
                Real dqp = ( fz_array(i,j,kc) - fz_array(i,j,k) ) / rho_array(i,j,k);
                Real omp = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
                Real omg = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tgrmin)*a_gr));

                qpr_array(i,j,k) = std::max(0.0, qpr_array(i,j,k) + dqp*omp);
                qps_array(i,j,k) = std::max(0.0, qps_array(i,j,k) + dqp*(1.0-omp)*(1.0-omg));
                qpg_array(i,j,k) = std::max(0.0, qpg_array(i,j,k) + dqp*(1.0-omp)*omg);
                 qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                /*
                // NOTE: Sedimentation does not affect the potential temperature,
                //       but it does affect the liquid/ice static  energy

                Real lat_heat = -( lfac_array(i,j,kc)*fz_array(i,j,kc) - lfac_array(i,j,k)*fz_array(i,j,k) ) / rho_array(i,j,k);
                amrex::Gpu::Atomic::Add(&theta_array(i,j,k), -lat_heat);
                */
            });

            if (iprec < nprec) {
                // Re-compute precipitation velocity using new value of qp.
                ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real tmp = term_vel_qp(qp_array(i,j,k),
                                           vrain, vsnow, vgrau,
                                           rho_array(i,j,k), tabs_array(i,j,k));
                    wp_array(i,j,k) = std::sqrt(1.29/rho_array(i,j,k))*tmp;
                    // Decrease precipitation velocity by factor of nprec
                    wp_array(i,j,k) = -wp_array(i,j,k)*rho_array(i,j,k)*dt_advance/dz/nprec;
                    // Note: Don't bother checking CFL condition at each
                    // substep since it's unlikely that the CFL will
                    // increase very much between substeps when using
                    // monotonic advection schemes.
                    if (k == 0) {
                          fz_array(i,j,nz-1) = 0.0;
                         www_array(i,j,nz-1) = 0.0;
                        lfac_array(i,j,nz-1) = 0.0;
                    }
                });
            }
        } // mfi
    } // iprec loop
}

