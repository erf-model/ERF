#include <AMReX_GpuContainers.H>
#include "ERF_SAM.H"
#include "ERF_IndexDefines.H"
#include "ERF_PlaneAverage.H"
#include "ERF_EOS.H"
#include "ERF_TileNoZ.H"

using namespace amrex;

namespace {

amrex::Real copy_table_value_to_host (const amrex::TableData<amrex::Real, 1>& src,
                                      const int lo,
                                      const int hi,
                                      const int k)
{
    amrex::TableData<amrex::Real, 1> host({lo}, {hi}, amrex::The_Pinned_Arena());
    host.copy(src);
    amrex::Gpu::streamSynchronize();
    return host.const_table()(k);
}

} // namespace

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
SAM::Init (const MultiFab& cons_in,
           const BoxArray& /*grids*/,
           const Geometry& geom,
           const Real& dt_advance,
           std::unique_ptr<MultiFab>& z_phys_nd,
           std::unique_ptr<MultiFab>& detJ_cc)
{
    dt = dt_advance;
    m_geom = geom;

    m_z_phys_nd = z_phys_nd.get();
    m_detJ_cc   = detJ_cc.get();

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar::rain_accum, MicVar::snow_accum, MicVar::graup_accum};

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar::NumVars; ++ivar) {
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
}


/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 */
void
SAM::Copy_State_to_Micro (const MultiFab& cons_in)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.growntilebox();

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
        ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const SAMPrimitiveCell primitive =
                sam_cons_to_primitive(states_array(i,j,k,Rho_comp),
                                      states_array(i,j,k,RhoTheta_comp),
                                      states_array(i,j,k,RhoQ1_comp),
                                      states_array(i,j,k,RhoQ2_comp),
                                      states_array(i,j,k,RhoQ3_comp),
                                      states_array(i,j,k,RhoQ4_comp),
                                      states_array(i,j,k,RhoQ5_comp),
                                      states_array(i,j,k,RhoQ6_comp));
            rho_array(i,j,k) = primitive.rho;
            theta_array(i,j,k) = primitive.theta;
            qv_array(i,j,k) = primitive.qv;
            qc_array(i,j,k) = primitive.qcl;
            qi_array(i,j,k) = primitive.qci;
            qn_array(i,j,k) = primitive.qn;
            qt_array(i,j,k) = primitive.qt;
            qpr_array(i,j,k) = primitive.qpr;
            qps_array(i,j,k) = primitive.qps;
            qpg_array(i,j,k) = primitive.qpg;
            qp_array(i,j,k) = primitive.qp;
            tabs_array(i,j,k) = primitive.tabs;
            pres_array(i,j,k) = primitive.pres_mbar;
        });
    }
}


void SAM::Compute_Coefficients ()
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

    ParallelFor(nlev, [=] AMREX_GPU_DEVICE (int k) noexcept
    {
        Real RhoTheta = rho_dptr[k]*theta_dptr[k];
        Real pressure = getPgivenRTh(RhoTheta, qv_dptr[k]);
        rho1d_t(k)    = rho_dptr[k];
        pres1d_t(k)   = sam_pa_to_mbar(pressure);
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
        const SAMCoefficientRow row =
            sam_compute_coefficient_row(rho1d_t(k), tabs1d_t(k),
                                        gamr1, gamr2, gams1, gams2, gamg1, gamg2);
        accrrc_t(k) = row.accrrc;
        accrsi_t(k) = row.accrsi;
        accrsc_t(k) = row.accrsc;
        coefice_t(k) = row.coefice;
        evaps1_t(k) = row.evaps1;
        evaps2_t(k) = row.evaps2;
        accrgi_t(k) = row.accrgi;
        accrgc_t(k) = row.accrgc;
        evapg1_t(k) = row.evapg1;
        evapg2_t(k) = row.evapg2;
        evapr1_t(k) = row.evapr1;
        evapr2_t(k) = row.evapr2;
    });
}

SAMCoefficientRow
SAM::CoefficientRowAt (const int k) const
{
    AMREX_ALWAYS_ASSERT(k >= zlo && k <= zhi);

    SAMCoefficientRow row{};
    row.accrrc = copy_table_value_to_host(accrrc, zlo, zhi, k);
    row.accrsi = copy_table_value_to_host(accrsi, zlo, zhi, k);
    row.accrsc = copy_table_value_to_host(accrsc, zlo, zhi, k);
    row.coefice = copy_table_value_to_host(coefice, zlo, zhi, k);
    row.evaps1 = copy_table_value_to_host(evaps1, zlo, zhi, k);
    row.evaps2 = copy_table_value_to_host(evaps2, zlo, zhi, k);
    row.accrgi = copy_table_value_to_host(accrgi, zlo, zhi, k);
    row.accrgc = copy_table_value_to_host(accrgc, zlo, zhi, k);
    row.evapg1 = copy_table_value_to_host(evapg1, zlo, zhi, k);
    row.evapg2 = copy_table_value_to_host(evapg2, zlo, zhi, k);
    row.evapr1 = copy_table_value_to_host(evapr1, zlo, zhi, k);
    row.evapr2 = copy_table_value_to_host(evapr2, zlo, zhi, k);
    return row;
}
