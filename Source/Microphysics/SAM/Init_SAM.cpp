#include <AMReX_GpuContainers.H>
#include "SAM.H"
#include "IndexDefines.H"
#include "PlaneAverage.H"
#include "EOS.H"
#include "TileNoZ.H"

using namespace amrex;

/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 * @param[in] qc_in Cloud variables input
 * @param[in,out] qv_in Vapor variables input
 * @param[in] qi_in Ice variables input
 * @param[in] grids The boxes on which we will evolve the solution
 * @param[in] geom Geometry associated with these MultiFabs and grids
 * @param[in] dt_advance Timestep for the advance
 */
void SAM::Init (const MultiFab& cons_in,
                const BoxArray& grids,
                const Geometry& geom,
                const Real& dt_advance)
{
    dt = dt_advance;
    m_geom = geom;
    m_gtoe = grids;

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar::qt, MicVar::qv , MicVar::qcl, MicVar::qci,
                 MicVar::qp, MicVar::qpr, MicVar::qps, MicVar::qpg};

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.);
    }

    // Set class data members
    for ( MFIter mfi(cons_in, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        const auto& lo = amrex::lbound(box3d);
        const auto& hi = amrex::ubound(box3d);

        nlev = box3d.length(2);
        zlo  = lo.z;
        zhi  = hi.z;

        // parameters
        accrrc.resize({zlo},  {zhi});
        accrsi.resize({zlo},  {zhi});
        accrsc.resize({zlo},  {zhi});
        coefice.resize({zlo}, {zhi});
        evaps1.resize({zlo},  {zhi});
        evaps2.resize({zlo},  {zhi});
        accrgi.resize({zlo},  {zhi});
        accrgc.resize({zlo},  {zhi});
        evapg1.resize({zlo},  {zhi});
        evapg2.resize({zlo},  {zhi});
        evapr1.resize({zlo},  {zhi});
        evapr2.resize({zlo},  {zhi});

        // data (input)
        rho1d.resize({zlo}, {zhi});
        pres1d.resize({zlo}, {zhi});
        tabs1d.resize({zlo}, {zhi});
        gamaz.resize({zlo}, {zhi});
        zmid.resize({zlo}, {zhi});
    }
}


/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 */
void SAM::Copy_State_to_Micro (const MultiFab& cons_in)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto states_array = cons_in.array(mfi);

        // Non-precipitating
        auto qv_array    = mic_fab_vars[MicVar::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar::qcl]->array(mfi);
        auto qi_array    = mic_fab_vars[MicVar::qci]->array(mfi);
        auto qn_array    = mic_fab_vars[MicVar::qn]->array(mfi);
        auto qt_array    = mic_fab_vars[MicVar::qt]->array(mfi);

        // Precipitating
        auto qpr_array   = mic_fab_vars[MicVar::qpr]->array(mfi);
        auto qps_array   = mic_fab_vars[MicVar::qps]->array(mfi);
        auto qpg_array   = mic_fab_vars[MicVar::qpg]->array(mfi);
        auto qp_array    = mic_fab_vars[MicVar::qp]->array(mfi);

        auto rho_array   = mic_fab_vars[MicVar::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar::pres]->array(mfi);

        // Get pressure, theta, temperature, density, and qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            theta_array(i,j,k) = states_array(i,j,k,RhoTheta_comp)/states_array(i,j,k,Rho_comp);

            qv_array(i,j,k)    = states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp);
            qc_array(i,j,k)    = states_array(i,j,k,RhoQ2_comp)/states_array(i,j,k,Rho_comp);
            qi_array(i,j,k)    = states_array(i,j,k,RhoQ3_comp)/states_array(i,j,k,Rho_comp);
            qn_array(i,j,k)    = qc_array(i,j,k) + qi_array(i,j,k);
            qt_array(i,j,k)    = qv_array(i,j,k) + qn_array(i,j,k);

            qpr_array(i,j,k)   = states_array(i,j,k,RhoQ4_comp)/states_array(i,j,k,Rho_comp);
            qps_array(i,j,k)   = states_array(i,j,k,RhoQ5_comp)/states_array(i,j,k,Rho_comp);
            qpg_array(i,j,k)   = states_array(i,j,k,RhoQ6_comp)/states_array(i,j,k,Rho_comp);
             qp_array(i,j,k)   = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

            tabs_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),
                                                  states_array(i,j,k,RhoTheta_comp),
                                                  qv_array(i,j,k));
            pres_array(i,j,k)  = getPgivenRTh(states_array(i,j,k,RhoTheta_comp), qv_array(i,j,k))/100.;
        });
    }
}


void SAM::Compute_Coefficients ()
{
    auto dz   = m_geom.CellSize(2);
    auto lowz = m_geom.ProbLo(2);

    auto accrrc_t  = accrrc.table();
    auto accrsi_t  = accrsi.table();
    auto accrsc_t  = accrsc.table();
    auto coefice_t = coefice.table();
    auto evaps1_t  = evaps1.table();
    auto evaps2_t  = evaps2.table();
    auto accrgi_t  = accrgi.table();
    auto accrgc_t  = accrgc.table();
    auto evapg1_t  = evapg1.table();
    auto evapg2_t  = evapg2.table();
    auto evapr1_t  = evapr1.table();
    auto evapr2_t  = evapr2.table();

    auto rho1d_t  = rho1d.table();
    auto pres1d_t = pres1d.table();
    auto tabs1d_t = tabs1d.table();

    auto gamaz_t  = gamaz.table();
    auto zmid_t   = zmid.table();

    Real gam3  = erf_gammafff(3.0             );
    Real gamr1 = erf_gammafff(3.0+b_rain      );
    Real gamr2 = erf_gammafff((5.0+b_rain)/2.0);
    Real gams1 = erf_gammafff(3.0+b_snow      );
    Real gams2 = erf_gammafff((5.0+b_snow)/2.0);
    Real gamg1 = erf_gammafff(3.0+b_grau      );
    Real gamg2 = erf_gammafff((5.0+b_grau)/2.0);

    // calculate the plane average variables
    PlaneAverage rho_ave(mic_fab_vars[MicVar::rho].get(), m_geom, m_axis);
    PlaneAverage theta_ave(mic_fab_vars[MicVar::theta].get(), m_geom, m_axis);
    PlaneAverage qv_ave(mic_fab_vars[MicVar::qv].get(), m_geom, m_axis);
    rho_ave.compute_averages(ZDir(), rho_ave.field());
    theta_ave.compute_averages(ZDir(), theta_ave.field());
    qv_ave.compute_averages(ZDir(), qv_ave.field());

    // get host variable rho, and rhotheta
    int ncell = rho_ave.ncell_line();

    Gpu::HostVector<Real> rho_h(ncell), theta_h(ncell), qv_h(ncell);
    rho_ave.line_average(0, rho_h);
    theta_ave.line_average(0, theta_h);
    qv_ave.line_average(0, qv_h);

    // copy data to device
    Gpu::DeviceVector<Real> rho_d(ncell), theta_d(ncell), qv_d(ncell);
    Gpu::copyAsync(Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
    Gpu::copyAsync(Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());
    Gpu::copyAsync(Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    Gpu::streamSynchronize();

    Real* rho_dptr   = rho_d.data();
    Real* theta_dptr = theta_d.data();
    Real* qv_dptr    = qv_d.data();

    Real gOcp = m_gOcp;

    ParallelFor(nlev, [=] AMREX_GPU_DEVICE (int k) noexcept
    {
        Real RhoTheta = rho_dptr[k]*theta_dptr[k];
        Real pressure = getPgivenRTh(RhoTheta, qv_dptr[k]);
        rho1d_t(k)    = rho_dptr[k];
        pres1d_t(k)   = pressure/100.;
        tabs1d_t(k)   = getTgivenRandRTh(rho_dptr[k], RhoTheta, qv_dptr[k]);
        zmid_t(k)     = lowz + (k+0.5)*dz;
        gamaz_t(k)    = gOcp*zmid_t(k);
    });

    if(round(gam3) != 2) {
        std::cout << "cannot compute gamma-function in Microphysics::Init" << std::endl;
        std::exit(-1);
    }

    // Populate all the coefficients
    ParallelFor(nlev, [=] AMREX_GPU_DEVICE (int k) noexcept
    {
        Real pratio = sqrt(1.29 / rho1d_t(k));
        Real rrr1   = 393.0/(tabs1d_t(k)+120.0)*std::pow((tabs1d_t(k)/273.0),1.5);
        Real rrr2   = std::pow((tabs1d_t(k)/273.0),1.94)*(1000.0/pres1d_t(k));
        Real estw   = 100.0*erf_esatw(tabs1d_t(k));
        Real esti   = 100.0*erf_esati(tabs1d_t(k));

        // accretion by snow:
        Real coef1   = 0.25 * PI * nzeros * a_snow * gams1 * pratio/pow((PI * rhos * nzeros/rho1d_t(k) ) , ((3.0+b_snow)/4.0));
        Real coef2   = exp(0.025*(tabs1d_t(k) - 273.15));
        accrsi_t(k)  =  coef1 * coef2 * esicoef;
        accrsc_t(k)  =  coef1 * esccoef;
        coefice_t(k) =  coef2;

        // evaporation of snow:
        coef1 = (lsub/(tabs1d_t(k)*rv)-1.0)*lsub/(therco*rrr1*tabs1d_t(k));
        coef2 = rv*tabs1d_t(k)/(diffelq*rrr2*esti);
        evaps1_t(k) = 0.65*4.0*nzeros/sqrt(PI*rhos*nzeros)/(coef1+coef2)/sqrt(rho1d_t(k));
        evaps2_t(k) = 0.49*4.0*nzeros*gams2*sqrt(a_snow/(muelq*rrr1))/pow((PI*rhos*nzeros) , ((5.0+b_snow)/8.0)) /
                      (coef1+coef2) * pow(rho1d_t(k) , ((1.0+b_snow)/8.0))*sqrt(pratio);

        // accretion by graupel:
        coef1 = 0.25*PI*nzerog*a_grau*gamg1*pratio/pow((PI*rhog*nzerog/rho1d_t(k)) , ((3.0+b_grau)/4.0));
        coef2 = exp(0.025*(tabs1d_t(k) - 273.15));
        accrgi_t(k) = coef1 * coef2 * egicoef;
        accrgc_t(k) = coef1 * egccoef;

        // evaporation of graupel:
        coef1 = (lsub/(tabs1d_t(k)*rv)-1.0)*lsub/(therco*rrr1*tabs1d_t(k));
        coef2 = rv*tabs1d_t(k)/(diffelq*rrr2*esti);
        evapg1_t(k) = 0.65*4.0*nzerog/sqrt(PI*rhog*nzerog)/(coef1+coef2)/sqrt(rho1d_t(k));
        evapg2_t(k) = 0.49*4.0*nzerog*gamg2*sqrt(a_grau/(muelq*rrr1))/pow((PI * rhog * nzerog) , ((5.0+b_grau)/8.0)) /
                     (coef1+coef2) * pow(rho1d_t(k) , ((1.0+b_grau)/8.0))*sqrt(pratio);

        // accretion by rain:
        accrrc_t(k) = 0.25 * PI * nzeror * a_rain * gamr1 * pratio/pow((PI * rhor * nzeror / rho1d_t(k)) , ((3+b_rain)/4.))* erccoef;

        // evaporation of rain:
        coef1 = (lcond/(tabs1d_t(k)*rv)-1.0)*lcond/(therco*rrr1*tabs1d_t(k));
        coef2 = rv*tabs1d_t(k)/(diffelq * rrr2 * estw);
        evapr1_t(k) = 0.78 * 2.0 * PI * nzeror /
                      sqrt(PI * rhor * nzeror) / (coef1+coef2) / sqrt(rho1d_t(k));
        evapr2_t(k) = 0.31 * 2.0 * PI  * nzeror * gamr2 * 0.89 * sqrt(a_rain/(muelq*rrr1))/
                      pow((PI * rhor * nzeror) , ((5.0+b_rain)/8.0)) /
                      (coef1+coef2) * pow(rho1d_t(k) , ((1.0+b_rain)/8.0))*sqrt(pratio);
    });
}
