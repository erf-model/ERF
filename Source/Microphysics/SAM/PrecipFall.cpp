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
void SAM::PrecipFall ()
{
    Real eps = std::numeric_limits<Real>::epsilon();
    Real rho_0 = 1.29;

    Real gamr3 = erf_gammafff(4.0+b_rain);
    Real gams3 = erf_gammafff(4.0+b_snow);
    Real gamg3 = erf_gammafff(4.0+b_grau);

    Real vrain = (a_rain*gamr3/6.0)*pow((PI*rhor*nzeror),-crain);
    Real vsnow = (a_snow*gams3/6.0)*pow((PI*rhos*nzeros),-csnow);
    Real vgrau = (a_grau*gamg3/6.0)*pow((PI*rhog*nzerog),-cgrau);

    auto dz   = m_geom.CellSize(2);
    Real dtn  = dt;
    Real coef = dtn/dz;

    auto domain = m_geom.Domain();
    int k_lo = domain.smallEnd(2);
    int k_hi = domain.bigEnd(2);

    auto qpr   = mic_fab_vars[MicVar::qpr];
    auto qps   = mic_fab_vars[MicVar::qps];
    auto qpg   = mic_fab_vars[MicVar::qpg];
    auto qp    = mic_fab_vars[MicVar::qp];
    auto rho   = mic_fab_vars[MicVar::rho];
    auto tabs  = mic_fab_vars[MicVar::tabs];
    auto theta = mic_fab_vars[MicVar::theta];

    auto ba    = tabs->boxArray();
    auto dm    = tabs->DistributionMap();
    auto ngrow = tabs->nGrowVect();

    MultiFab fz;
    fz.define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow);

    //  Add sedimentation of precipitation field to the vert. vel.
    for (MFIter mfi(fz, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto qp_array   = qp->array(mfi);
        auto rho_array  = rho->array(mfi);
        auto tabs_array = tabs->array(mfi);
        auto fz_array   = fz.array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real rho_avg, tab_avg, qp_avg;
            if (k==k_lo) {
                rho_avg =  rho_array(i,j,k);
                tab_avg = tabs_array(i,j,k);
                 qp_avg =   qp_array(i,j,k);
            } else if (k==k_hi+1) {
                rho_avg =  rho_array(i,j,k-1);
                tab_avg = tabs_array(i,j,k-1);
                 qp_avg =   qp_array(i,j,k-1);
            } else {
                rho_avg = 0.5*( rho_array(i,j,k-1) +  rho_array(i,j,k));
                tab_avg = 0.5*(tabs_array(i,j,k-1) + tabs_array(i,j,k));
                 qp_avg = 0.5*(  qp_array(i,j,k-1) +   qp_array(i,j,k));
            }

            Real Pprecip = 0.0;
            if(qp_avg > qp_threshold) {
                Real omp = std::max(0.0,std::min(1.0,(tab_avg-tprmin)*a_pr));
                Real omg = std::max(0.0,std::min(1.0,(tab_avg-tgrmin)*a_gr));
                Real qrr = omp*qp_avg;
                Real qss = (1.0-omp)*(1.0-omg)*qp_avg;
                Real qgg = (1.0-omp)*(omg)*qp_avg;
                Pprecip = omp*vrain*std::pow(rho_avg*qrr,1.0+crain)
                        + (1.0-omp)*( (1.0-omg)*vsnow*std::pow(rho_avg*qss,1.0+csnow)
                                    +      omg *vgrau*std::pow(rho_avg*qgg,1.0+cgrau) );
            }
            fz_array(i,j,k) = Pprecip * std::sqrt(rho_0/rho_avg);
        });
    }

    for (MFIter mfi(*qp, TileNoZ()); mfi.isValid(); ++mfi) {
        auto qpr_array    = qpr->array(mfi);
        auto qps_array    = qps->array(mfi);
        auto qpg_array    = qpg->array(mfi);
        auto qp_array     = qp->array(mfi);
        auto rho_array    = rho->array(mfi);
        auto tabs_array   = tabs->array(mfi);
        auto fz_array     = fz.array(mfi);

        const auto& box3d = mfi.tilebox();

        // Update precipitation mass fraction and liquid-ice static
        // energy using precipitation fluxes computed in this column.
        ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //==================================================
            // Precipitating sedimentation (A19)
            //==================================================
            Real dqp = (1.0/rho_array(i,j,k)) * ( fz_array(i,j,k+1) - fz_array(i,j,k) ) * coef;
            Real omp = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
            Real omg = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tgrmin)*a_gr));

            qpr_array(i,j,k) = std::max(0.0, qpr_array(i,j,k) + dqp*omp);
            qps_array(i,j,k) = std::max(0.0, qps_array(i,j,k) + dqp*(1.0-omp)*(1.0-omg));
            qpg_array(i,j,k) = std::max(0.0, qpg_array(i,j,k) + dqp*(1.0-omp)*omg);
             qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

            // NOTE: Sedimentation does not affect the potential temperature,
            //       but it does affect the liquid/ice static energy.
            //       No source to Theta occurs here.
        });
    } // mfi
}

