#ifndef ERF_MACROPHYS_H
#define ERF_MACROPHYS_H

#include "ERF_EOS.H"
#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"
#include "ERF_MicrophysicsUtils.H"

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
static void
NewtonIterSat (const amrex::Real& fac_cond,
               amrex::Real& tabs,
               const amrex::Real& pres,
               amrex::Real& qv)
{
    // Solution tolerance
    amrex::Real tol = amrex::Real(1.0e-6);

    // Saturation moisture fractions
    amrex::Real qsat, dqsat;

    // Newton iteration vars
    int  niter;
    amrex::Real fff, dfff;
    amrex::Real dtabs;

    // Initial guess for temperature
    amrex::Real tabs_old = tabs;

    niter = 0;
    dtabs = 1;

    //==================================================
    // Newton iteration to qv=qsat (cloud phase only)
    //==================================================
    do {
        // Saturation moisture fractions
        erf_qsatw(tabs, pres, qsat);
        erf_dtqsatw(tabs, pres, dqsat);

        // Function for root finding:
        // 0 = -T_new + T_old + L_v/C_p * (qv - qsat)
        fff   = -tabs + tabs_old + fac_cond*(qv - qsat);

        // Derivative of function (T_new iterated on)
        dfff  = -one - fac_cond*dqsat;

        // Update the temperature
        dtabs = -fff/dfff;
        tabs += dtabs;

        // Update iteration
        niter = niter+1;
    } while(std::abs(dtabs) > tol && niter < 20);

    // Update from last iteration (dq = dq/dt * dt)
    qv  = qsat + dqsat*dtabs;
}

void
GetMacroPhysicsSrc (int& real_width,
                    amrex::Real& dt,
                    amrex::Geometry& geom,
                    amrex::MultiFab* macro_src,
                    amrex::MultiFab& cons_data)
{
    amrex::Real idt        = amrex::Real(1.) / dt;
    amrex::Real d_rdOcp    = R_d / Cp_d;
    amrex::Real d_fac_cond = L_v / Cp_d;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(cons_data,TileNoZ()); mfi.isValid(); ++mfi)
    {
        amrex::Box tbx = mfi.tilebox();

        const amrex::Array4<amrex::Real>& mac_src   = macro_src->array(mfi);
        const amrex::Array4<amrex::Real>& cell_data = cons_data.array(mfi);

        // No macrophysics in the nudging region
        const amrex::Box& domain = geom.Domain();
        int i_lo = domain.smallEnd(0);
        int i_hi = domain.bigEnd(0);
        int j_lo = domain.smallEnd(1);
        int j_hi = domain.bigEnd(1);
        if (tbx.smallEnd(0) == i_lo) { tbx.growLo(0,-real_width); }
        if (tbx.bigEnd(0)   == i_hi) { tbx.growHi(0,-real_width); }
        if (tbx.smallEnd(1) == j_lo) { tbx.growLo(1,-real_width); }
        if (tbx.bigEnd(1)   == j_hi) { tbx.growHi(1,-real_width); }

        // Theta tendency and update qv/qc in place
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Conserved CC vars at start of step
            amrex::Real rho_s = cell_data(i,j,k,Rho_comp);
            amrex::Real rt_s  = cell_data(i,j,k,RhoTheta_comp);
            amrex::Real th_s  = cell_data(i,j,k,RhoTheta_comp) / rho_s;
            amrex::Real qv_s  = cell_data(i,j,k,RhoQ1_comp   ) / rho_s;
            amrex::Real qc_s  = cell_data(i,j,k,RhoQ2_comp   ) / rho_s;

            // To be updated in place
            amrex::Real qs;
            amrex::Real qv_f   = qv_s;
            amrex::Real tabs_f = getTgivenRandRTh(rho_s, rt_s, qv_s);

            // Invariants
            amrex::Real qt   = qv_s + qc_s;
            amrex::Real pres = getPgivenRTh(rt_s, qv_s);

            // Valid Newton iteration
            erf_qsatw(tabs_f, pres*amrex::Real(0.01), qs);
            if (qt > qs) {
                NewtonIterSat(d_fac_cond, tabs_f, pres*amrex::Real(0.01), qv_f);
            }
            // Put qc into qv and double check Newton iter criteria
            else {
                qv_f   += qc_s;
                tabs_f -= d_fac_cond * qc_s;
                erf_qsatw(tabs_f, pres*amrex::Real(0.01), qs);
                if (qv_f > qs) {
                    NewtonIterSat(d_fac_cond, tabs_f, pres*amrex::Real(0.01), qv_f);
                }
            }

            // Populate tendency and data
            amrex::Real th_f = getThgivenTandP(tabs_f, pres, d_rdOcp);
            mac_src(i,j,k) = rho_s * idt * ( th_f - th_s );
            cell_data(i,j,k,RhoQ1_comp) = rho_s * qv_f;
            cell_data(i,j,k,RhoQ2_comp) = rho_s * std::max(amrex::Real(0.0),qt - qv_f);
        });
    } // mfi
}

#endif
