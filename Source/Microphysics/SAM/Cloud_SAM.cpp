
#include "Microphysics.H"
#include "IndexDefines.H"
#include "TileNoZ.H"

using namespace amrex;

/**
 * Compute Cloud-related Microphysics quantities.
 */
void SAM::Cloud () {

    constexpr Real an   = 1.0/(tbgmax-tbgmin);
    constexpr Real bn   = tbgmin*an;
    constexpr Real ap   = 1.0/(tprmax-tprmin);
    constexpr Real bp   = tprmin*ap;

    auto pres1d_t = pres1d.table();

    auto qt    = mic_fab_vars[MicVar::qt];
    auto qp    = mic_fab_vars[MicVar::qp];
    auto qn    = mic_fab_vars[MicVar::qn];
    auto rho   = mic_fab_vars[MicVar::rho];
    auto tabs  = mic_fab_vars[MicVar::tabs];

    Real fac_cond = m_fac_cond;
    Real fac_sub  = m_fac_sub;
    Real fac_fus  = m_fac_fus;

    for ( MFIter mfi(*tabs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto qt_array    = qt->array(mfi);
        auto qp_array    = qp->array(mfi);
        auto qn_array    = qn->array(mfi);
        auto tabs_array  = tabs->array(mfi);

        const auto& box3d = mfi.tilebox() & m_gtoe[mfi.index()];

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            qt_array(i,j,k) = std::max(0.0,qt_array(i,j,k));
            // Initial guess for temperature assuming no cloud water/ice:
            Real tabs1 = tabs_array(i,j,k);

            Real qsatt;
            Real om;
            Real qsatt1;
            Real qsatt2;

            // Warm cloud:
            if(tabs1 > tbgmax) {
                tabs1 = tabs_array(i,j,k)+fac_cond*qp_array(i,j,k);
                erf_qsatw(tabs1, pres1d_t(k), qsatt);
            }
            // Ice cloud:
            else if(tabs1 <= tbgmin) {
                tabs1 = tabs_array(i,j,k)+fac_sub*qp_array(i,j,k);
                erf_qsati(tabs1, pres1d_t(k), qsatt);
            }
            // Mixed-phase cloud:
            else {
                om = an*tabs1-bn;
                erf_qsatw(tabs1, pres1d_t(k), qsatt1);
                erf_qsati(tabs1, pres1d_t(k), qsatt2);
                qsatt = om*qsatt1 + (1.-om)*qsatt2;
            }

            int niter;
            Real dtabs, lstarn, dlstarn, omp, lstarp, dlstarp, fff, dfff, dqsat;
            //  Test if condensation is possible:
            if(qt_array(i,j,k) > qsatt) {
                niter = 0;
                dtabs = 1;
                do {
                    if(tabs1 >= tbgmax) {
                        om=1.0;
                        lstarn  = fac_cond;
                        dlstarn = 0.0;
                        erf_qsatw(tabs1, pres1d_t(k), qsatt);
                        erf_dtqsatw(tabs1, pres1d_t(k), dqsat);
                    }
                    else if(tabs1 <= tbgmin) {
                        om      = 0.0;
                        lstarn  = fac_sub;
                        dlstarn = 0.0;
                        erf_qsati(tabs1, pres1d_t(k), qsatt);
                        erf_dtqsati(tabs1, pres1d_t(k), dqsat);
                    }
                    else {
                        om=an*tabs1-bn;
                        lstarn  = fac_cond+(1.0-om)*fac_fus;
                        dlstarn = an*fac_fus;
                        erf_qsatw(tabs1, pres1d_t(k), qsatt1);
                        erf_qsati(tabs1, pres1d_t(k), qsatt2);

                        qsatt = om*qsatt1+(1.-om)*qsatt2;
                        erf_dtqsatw(tabs1, pres1d_t(k), qsatt1);
                        erf_dtqsati(tabs1, pres1d_t(k), qsatt2);
                        dqsat = om*qsatt1+(1.-om)*qsatt2;
                    }

                    if(tabs1 >= tprmax) {
                        omp = 1.0;
                        lstarp  = fac_cond;
                        dlstarp = 0.0;
                    }
                    else if(tabs1 <= tprmin) {
                        omp     = 0.0;
                        lstarp  = fac_sub;
                        dlstarp = 0.0;
                    }
                    else {
                        omp=ap*tabs1-bp;
                        lstarp  = fac_cond+(1.0-omp)*fac_fus;
                        dlstarp = ap*fac_fus;
                    }
                    fff   = tabs_array(i,j,k)-tabs1+lstarn*(qt_array(i,j,k)-qsatt)+lstarp*qp_array(i,j,k);
                    dfff  = dlstarn*(qt_array(i,j,k)-qsatt)+dlstarp*qp_array(i,j,k)-lstarn*dqsat-1.0;
                    dtabs = -fff/dfff;
                    niter = niter+1;
                    tabs1 = tabs1+dtabs;
                } while(std::abs(dtabs) > 0.01 && niter < 10);
                qsatt = qsatt + dqsat*dtabs;
                qn_array(i,j,k) = std::max(0.0, qt_array(i,j,k)-qsatt);
            }
            else {
                qn_array(i,j,k) = 0.0;
            }
            tabs_array(i,j,k) = tabs1;
            qp_array(i,j,k)   = std::max(0.0, qp_array(i,j,k)); // just in case
        });
    }
}
