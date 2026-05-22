#include <AMReX_GpuContainers.H>
#include "ERF_Morrison.H"
#include "ERF_IndexDefines.H"
#include "ERF_PlaneAverage.H"
#include "ERF_EOS.H"
#include "ERF_TileNoZ.H"

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
void
Morrison::Init (const MultiFab& cons_in,
                const BoxArray& /*grids*/,
                const Geometry& geom,
                const Real& dt_advance,
                std::unique_ptr<MultiFab>& z_phys_nd,
                std::unique_ptr<MultiFab>& detJ_cc)
{
    [[maybe_unused]] amrex::Real dt     = dt_advance;
    m_geom = geom;

    m_z_phys_nd = z_phys_nd.get();
    m_detJ_cc   = detJ_cc.get();

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar_Morr::nc, MicVar_Morr::nr, MicVar_Morr::ni, MicVar_Morr::ns, MicVar_Morr::ng,
                 MicVar_Morr::rain_accum, MicVar_Morr::snow_accum, MicVar_Morr::graup_accum};

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar_Morr::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.);
    }

    // NOTE: For multi-level not all ranks will own a box.
    //       Furthermore, the plane average allocates space
    //       for the entire domain. We make this consistent.
    nlev = m_geom.Domain().length(2);
    zlo  = m_geom.Domain().smallEnd(2);
    zhi  = m_geom.Domain().bigEnd(2);

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

#ifdef ERF_USE_MORR_FORT
    int morr_rimed_ice = 0; // This is used to set something called "ihail"
    amrex::ParmParse pp("erf");
    MoistureType moisture_type;
    pp.query_enum_case_insensitive("moisture_model",moisture_type);
    int morr_noice = (moisture_type == MoistureType::Morrison_NoIce);
    Print()<<"Setting No Ice flag in fortran to "<<morr_noice<<std::endl;
    morr_two_moment_init_c(morr_rimed_ice, morr_noice);
#endif
}


/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 */
void
Morrison::Copy_State_to_Micro (const MultiFab& cons_in)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.growntilebox();

        auto states_array = cons_in.array(mfi);

        // Non-precipitating
        auto qv_array    = mic_fab_vars[MicVar_Morr::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_Morr::qcl]->array(mfi);
        auto qi_array    = mic_fab_vars[MicVar_Morr::qci]->array(mfi);
        auto qn_array    = mic_fab_vars[MicVar_Morr::qn]->array(mfi);
        auto qt_array    = mic_fab_vars[MicVar_Morr::qt]->array(mfi);

        // Precipitating
        auto qpr_array   = mic_fab_vars[MicVar_Morr::qpr]->array(mfi);
        auto qps_array   = mic_fab_vars[MicVar_Morr::qps]->array(mfi);
        auto qpg_array   = mic_fab_vars[MicVar_Morr::qpg]->array(mfi);
        auto qp_array    = mic_fab_vars[MicVar_Morr::qp]->array(mfi);

        auto rho_array   = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_Morr::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_Morr::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_Morr::pres]->array(mfi);

        // Get pressure, theta, temperature, density, and qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            theta_array(i,j,k) = states_array(i,j,k,RhoTheta_comp)/states_array(i,j,k,Rho_comp);

            qv_array(i,j,k)    = std::max(Real(0),states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp));
            qc_array(i,j,k)    = std::max(Real(0),states_array(i,j,k,RhoQ2_comp)/states_array(i,j,k,Rho_comp));
            qi_array(i,j,k)    = std::max(Real(0),states_array(i,j,k,RhoQ3_comp)/states_array(i,j,k,Rho_comp));
            qn_array(i,j,k)    = qc_array(i,j,k) + qi_array(i,j,k);
            qt_array(i,j,k)    = qv_array(i,j,k) + qn_array(i,j,k);

            qpr_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ4_comp)/states_array(i,j,k,Rho_comp));
            qps_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ5_comp)/states_array(i,j,k,Rho_comp));
            qpg_array(i,j,k)   = std::max(Real(0),states_array(i,j,k,RhoQ6_comp)/states_array(i,j,k,Rho_comp));
             qp_array(i,j,k)   = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

            tabs_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),
                                                  states_array(i,j,k,RhoTheta_comp),
                                                  qv_array(i,j,k));

            // NOTE: the Morrison Fortran version uses Pa not hPa so we don't divideby 100!
            pres_array(i,j,k)  = getPgivenRTh(states_array(i,j,k,RhoTheta_comp), qv_array(i,j,k)); //  * Real(0.01);
        });
    }
}


void Morrison::Compute_Coefficients ()
{
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

    Real gam3  = erf_gammafff(three             );
    Real gamr1 = erf_gammafff(three+b_rain      );
    Real gamr2 = erf_gammafff((Real(5.0)+b_rain)/two);
    Real gams1 = erf_gammafff(three+b_snow      );
    Real gams2 = erf_gammafff((Real(5.0)+b_snow)/two);
    Real gamg1 = erf_gammafff(three+b_grau      );
    Real gamg2 = erf_gammafff((Real(5.0)+b_grau)/two);

    // calculate the plane average variables
    PlaneAverage rho_ave(mic_fab_vars[MicVar_Morr::rho].get(), m_geom, m_axis);
    PlaneAverage theta_ave(mic_fab_vars[MicVar_Morr::theta].get(), m_geom, m_axis);
    PlaneAverage qv_ave(mic_fab_vars[MicVar_Morr::qv].get(), m_geom, m_axis);
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

    ParallelFor(nlev, [=] AMREX_GPU_DEVICE (int k) noexcept
    {
        Real RhoTheta = rho_dptr[k]*theta_dptr[k];
        Real pressure = getPgivenRTh(RhoTheta, qv_dptr[k]);
        rho1d_t(k)    = rho_dptr[k];
        pres1d_t(k)   = pressure*Real(0.01);
        // NOTE: Limit the temperature to the melting point of ice to avoid a divide by
        //       0 condition when computing the cold evaporation coefficients. This should
        //       not affect results since evaporation requires snow/graupel to be present
        //       and thus T<Real(273.16)
        tabs1d_t(k)   = std::min(getTgivenRandRTh(rho_dptr[k], RhoTheta, qv_dptr[k]),Real(273.16));
    });

    if(round(gam3) != 2) {
        std::cout << "cannot compute gamma-function in Microphysics::Init" << std::endl;
        std::exit(-1);
    }

    // Populate all the coefficients
    ParallelFor(nlev, [=] AMREX_GPU_DEVICE (int k) noexcept
    {
        Real Prefactor;
        Real pratio = std::sqrt(Real(1.29) / rho1d_t(k));
        //Real rrr1   = Real(393.0)/(tabs1d_t(k)+Real(120.0))*std::pow((tabs1d_t(k)/Real(273.0)),Real(1.5));
        //Real rrr2   = std::pow((tabs1d_t(k)/Real(273.0)),Real(1.94))*(Real(1000.0)/pres1d_t(k));
        Real estw   = Real(100.0)*erf_esatw(tabs1d_t(k));
        Real esti   = Real(100.0)*erf_esati(tabs1d_t(k));

        // accretion by snow:
        Real coef1   = fourth * PI * nzeros * a_snow * gams1 * pratio/std::pow((PI * rhos * nzeros/rho1d_t(k) ) , ((three+b_snow)/Real(4.0)));
        Real coef2   = std::exp(Real(0.025)*(tabs1d_t(k) - Real(273.15)));
        accrsi_t(k)  =  coef1 * coef2 * esicoef;
        accrsc_t(k)  =  coef1 * esccoef;
        coefice_t(k) =  coef2;

        // evaporation of snow:
        coef1 = (lsub/(tabs1d_t(k)*R_v)-one)*lsub/(therco*tabs1d_t(k));
        coef2 = R_v * R_d / (diffelq * esti);
        Prefactor = two * PI * nzeros / (rho1d_t(k) * (coef1 + coef2));
        Prefactor *= (two/PI); // Shape factor snow
        evaps1_t(k) = Prefactor * Real(0.65) * std::sqrt(rho1d_t(k) / (PI * rhos * nzeros));
        evaps2_t(k) = Prefactor * Real(0.44) * std::sqrt(a_snow * rho1d_t(k) / muelq) * gams2
                    * std::sqrt(pratio) * std::pow(rho1d_t(k) / (PI * rhos * nzeros) , ((Real(5.0)+b_snow)/Real(8.0)));

        // accretion by graupel:
        coef1 = fourth*PI*nzerog*a_grau*gamg1*pratio/std::pow((PI*rhog*nzerog/rho1d_t(k)) , ((three+b_grau)/Real(4.0)));
        coef2 = std::exp(Real(0.025)*(tabs1d_t(k) - Real(273.15)));
        accrgi_t(k) = coef1 * coef2 * egicoef;
        accrgc_t(k) = coef1 * egccoef;

        // evaporation of graupel:
        coef1 = (lsub/(tabs1d_t(k)*R_v)-one)*lsub/(therco*tabs1d_t(k));
        coef2 = R_v * R_d / (diffelq * esti);
        Prefactor = two * PI * nzerog / (rho1d_t(k) * (coef1 + coef2)); // Shape factor for graupel is 1
        evapg1_t(k) = Prefactor * Real(0.78) * std::sqrt(rho1d_t(k) / (PI * rhog * nzerog));
        evapg2_t(k) = Prefactor * Real(0.31) * std::sqrt(a_grau * rho1d_t(k) / muelq) * gamg2
                    * std::sqrt(pratio) * std::pow(rho1d_t(k) / (PI * rhog * nzerog) , ((Real(5.0)+b_grau)/Real(8.0)));

        // accretion by rain:
        accrrc_t(k) = fourth * PI * nzeror * a_rain * gamr1 * pratio/std::pow((PI * rhor * nzeror / rho1d_t(k)) , ((three+b_rain)/Real(4.)))* erccoef;

        // evaporation of rain:
        coef1 = (lcond/(tabs1d_t(k)*R_v)-one)*lcond/(therco*tabs1d_t(k));
        coef2 = R_v * R_d / (diffelq * estw);
        Prefactor = two * PI * nzeror / (rho1d_t(k) * (coef1 + coef2)); // Shape factor for rain is 1
        evapr1_t(k) = Prefactor * Real(0.78) * std::sqrt(rho1d_t(k) / (PI * rhor * nzeror));
        evapr2_t(k) = Prefactor * Real(0.31) * std::sqrt(a_rain * rho1d_t(k) / muelq) * gamr2
                    * std::sqrt(pratio) * std::pow(rho1d_t(k) / (PI * rhor * nzeror) , ((Real(5.0)+b_rain)/Real(8.0)));
    });
}
