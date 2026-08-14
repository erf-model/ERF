#include "ERF_SuperDropletsMoist.H"
#include "ERF_EOS.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! \brief Compute saturation ratio */
static void saturation_ratio (MultiFab& a_sr,
                              const MultiFab& a_qv,
                              const MultiFab& a_temperature,
                              const MultiFab& a_pressure,
                              const MaterialProperties& a_species_mat )
{
    a_species_mat.computeSaturationVapFrac(a_sr, a_temperature, a_pressure);
    for (MFIter mfi(a_sr); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        bx.grow( a_sr.nGrowVect() );
        const Array4<Real>& sr_arr = a_sr.array(mfi);
        const Array4<Real const>& qv_arr = a_qv.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        { sr_arr(i,j,k,0) = (sr_arr(i,j,k,0) > Real(0) ? qv_arr(i,j,k,0) / sr_arr(i,j,k,0) : Real(0)); });
    }
    a_sr.FillBoundary();
}

/*! Compute phase change for a timestep */
void SuperDropletsMoist::phaseChange ( const Real& a_dt,
                                       const Vector<std::unique_ptr<MultiFab>>& a_z,
                                       const int a_lev )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange()");
    auto dt_s = a_dt / static_cast<Real>(m_num_substeps_phase_change);

    // Mask of cells not covered by a finer level.  Built once and shared by
    // every sub-phase-change call so that cell-centred q/T updates skip
    // coarse cells whose values will be overwritten by averaging-down.
    iMultiFab fine_mask = buildFineMask(a_lev);

    for (int substep = 0; substep < m_num_substeps_phase_change; substep++) {

        // freezing/melting (solid <--> liquid) -  water
        if (m_with_ice) { phaseChange_SL_w(dt_s, a_z, a_lev, fine_mask); }

        // evaporation and condensation (vapour <--> liquid) - water
        phaseChange_LV_w(dt_s, a_z, a_lev, fine_mask);

        // deposition and sublimation (vapour <--> solid) - water (ice)
        if (m_with_ice) { phaseChange_SV_i(dt_s, a_z, a_lev, fine_mask); }

        // evaporation and condensation (vapour <--> liquid) - other species
        for (int is = m_istart_sp; is < m_num_species; is++) {
            phaseChange_LV_s(is, dt_s, a_z, a_lev, fine_mask);
        }

        // Update pressure and temperature
        const auto& gvec = (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]).nGrowVect();
        for (MFIter mfi(*m_mic_fab_vars[a_lev][MicVar_SD::temperature], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            Box bx = mfi.tilebox();
            bx.grow(gvec);

            const Array4<Real>& t_arr = m_mic_fab_vars[a_lev][MicVar_SD::temperature]->array(mfi);
            const Array4<Real>& p_arr = m_mic_fab_vars[a_lev][MicVar_SD::pressure]->array(mfi);

            const Array4<Real const>& rho_arr = m_mic_fab_vars[a_lev][MicVar_SD::rho]->const_array(mfi);
            const Array4<Real const>& theta_arr = m_mic_fab_vars[a_lev][MicVar_SD::theta]->const_array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[a_lev][MicVar_SD::q_v]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                t_arr(i,j,k,0) = getTgivenRandRTh( rho_arr(i,j,k,0),
                                                   rho_arr(i,j,k,0)*theta_arr(i,j,k,0),
                                                   qv_arr(i,j,k,0) );
                p_arr(i,j,k,0) = getPgivenRTh( rho_arr(i,j,k,0)*theta_arr(i,j,k,0),
                                               qv_arr(i,j,k,0) );
            });
        }

        // Update saturation ratio for each species
        for (int is = 0; is < m_num_species; is++) {
            auto& species = m_species[is];
            auto& species_mat = m_super_droplets->getSpeciesMaterial(species);

            MultiFab *qv_ptr, *sr_ptr;
            if (is == m_idx_w) {
                AMREX_ALWAYS_ASSERT(species==Species::Name::H2O);
                qv_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_v].get();
                sr_ptr = m_mic_fab_vars[a_lev][MicVar_SD::rh_w].get();
            } else if (is == m_idx_i) {
                AMREX_ALWAYS_ASSERT(species==Species::Name::ice);
                qv_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_v].get();
                sr_ptr = m_mic_fab_vars[a_lev][MicVar_SD::rh_i].get();
            } else {
                qv_ptr = m_mic_fab_vars[a_lev][s_qv_idx(is,m_istart_sp)].get();
                sr_ptr = m_mic_fab_vars[a_lev][s_sr_idx(is,m_istart_sp)].get();
            }

            saturation_ratio( (*sr_ptr),
                              (*qv_ptr),
                              (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                              (*m_mic_fab_vars[a_lev][MicVar_SD::pressure]),
                              species_mat );

        }

    }
}

/*! Condensation/evaporation (liquid <--> vapour) for water: shrink/grow the
 * super-droplet particles for a timestep, depending on the ambient flow conditions
 * (saturation ratio,  saturation pressure, and temperature). Update the Eulerian vapour
 * and condensate mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange_LV_w ( const Real& a_dt,
                                            const Vector<std::unique_ptr<MultiFab>>& a_z,
                                            const int a_lev,
                                            const iMultiFab& fine_mask )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange_LV_w()");

    auto& species = m_species[m_idx_w];
    auto& species_mat = m_super_droplets->getSpeciesMaterial(species);
    bool is_water = species_mat.m_is_water;

    AMREX_ALWAYS_ASSERT(species == Species::Name::H2O);
    AMREX_ALWAYS_ASSERT(is_water);

    auto* qv_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_v].get();
    auto* qc_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_c].get();
    auto* qr_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_r].get();
    auto* qt_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_t].get();
    auto* sr_ptr = m_mic_fab_vars[a_lev][MicVar_SD::rh_w].get();

    // Compute saturation pressure
    MultiFab mf_sat_pressure(   m_mic_fab_vars[a_lev][MicVar_SD::pressure]->boxArray(),
                                m_mic_fab_vars[a_lev][MicVar_SD::pressure]->DistributionMap(),
                                1,
                                m_mic_fab_vars[a_lev][MicVar_SD::pressure]->nGrowVect() );
    species_mat.computeSaturationPressure( mf_sat_pressure,
                                           (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]) );
    mf_sat_pressure.FillBoundary();

    // Compute saturation ratio
    saturation_ratio( (*sr_ptr),
                      (*qv_ptr),
                      (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                      (*m_mic_fab_vars[a_lev][MicVar_SD::pressure]),
                      species_mat );

    // Recompute qc/qr from the current particles so qt below is built from the
    // particle-carried condensate, not the dycore-advected condensate left in
    // the micro fabs by Copy_State_to_Micro (which leaks the Eulerian-vs-
    // particle advection mismatch into qv as spurious supersaturation).
    computeQcQrWater(*a_z[a_lev]);

    // Compute total liquid and vapour water content qt
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->array(mfi);
        auto qv_arr = qv_ptr->const_array(mfi);
        auto qc_arr = qc_ptr->const_array(mfi);
        auto qr_arr = qr_ptr->const_array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { qt_arr(i,j,k) = qv_arr(i,j,k) + qc_arr(i,j,k) + qr_arr(i,j,k); });
    }

    // Compute super-droplets mass change
    m_super_droplets->MassChange_LV ( a_lev,
                                      a_dt,
                                      m_species[m_idx_w],
                                      (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                                      (*m_mic_fab_vars[a_lev][MicVar_SD::pressure]),
                                      mf_sat_pressure,
                                      (*sr_ptr),
                                      a_z,
                                      is_water );

    // Compute new condensate mixing ratio
    computeQcQrWater(*a_z[a_lev]);

    // Update vapour mixing ratio
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->const_array(mfi);
        auto qv_arr = qv_ptr->array(mfi);
        auto qc_arr = qc_ptr->array(mfi);
        auto qr_arr = qr_ptr->array(mfi);

        auto theta_arr = m_mic_fab_vars[a_lev][MicVar_SD::theta]->array(mfi);
        auto T_arr = m_mic_fab_vars[a_lev][MicVar_SD::temperature]->const_array(mfi);
        auto dqc_arr = m_mic_fab_vars[a_lev][MicVar_SD::dqcdt]->array(mfi);

        auto mask_arr = fine_mask.const_array(mfi);

        auto fac_cond = species_mat.m_lat_vap / m_Cp;

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mask_arr(i,j,k) == 0) { return; }
            auto old_qv = qv_arr(i,j,k);
            auto qw = qc_arr(i,j,k) + qr_arr(i,j,k);
            if (qw > qt_arr(i,j,k)) {
                qv_arr(i,j,k) = 0.0;
                if (qr_arr(i,j,k) > qt_arr(i,j,k)) {
                    qc_arr(i,j,k) = 0.0;
                    qr_arr(i,j,k) = qt_arr(i,j,k);
                } else {
                    qc_arr(i,j,k) = qt_arr(i,j,k) - qr_arr(i,j,k);
                }
            } else {
                qv_arr(i,j,k) = qt_arr(i,j,k) - qw;
            }

            AMREX_ALWAYS_ASSERT(qr_arr(i,j,k) >= 0.0);
            AMREX_ALWAYS_ASSERT(qv_arr(i,j,k) >= 0.0);
            AMREX_ALWAYS_ASSERT(qc_arr(i,j,k) >= 0.0);

            dqc_arr(i,j,k) = - (qv_arr(i,j,k) - old_qv) / a_dt;

            auto theta_over_T = theta_arr(i,j,k)/T_arr(i,j,k);
            theta_arr(i,j,k) += theta_over_T * fac_cond * (old_qv-qv_arr(i,j,k));
        });
    }

}

/*! Condensation/evaporation (liquid <--> vapour) for other species: shrink/grow the
 * super-droplet particles for a timestep, depending on the ambient flow conditions
 * (saturation ratio,  saturation pressure, and temperature). Update the Eulerian vapour
 * and condensate mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange_LV_s ( const int a_idx,
                                            const Real& a_dt,
                                            const Vector<std::unique_ptr<MultiFab>>& a_z,
                                            const int a_lev,
                                            const iMultiFab& fine_mask )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange_LV_s()");

    auto& species = m_species[a_idx];
    auto& species_mat = m_super_droplets->getSpeciesMaterial(species);
    bool is_water = species_mat.m_is_water; // may be true for dummy water species used for debugging

    AMREX_ALWAYS_ASSERT((species!=Species::Name::ice) || (species!=Species::Name::H2O));

    auto* qv_ptr = m_mic_fab_vars[a_lev][s_qv_idx(a_idx,m_istart_sp)].get();
    auto* qc_ptr = m_mic_fab_vars[a_lev][s_qc_idx(a_idx,m_istart_sp)].get();
    auto* qt_ptr = m_mic_fab_vars[a_lev][s_qt_idx(a_idx,m_istart_sp)].get();
    auto* sr_ptr = m_mic_fab_vars[a_lev][s_sr_idx(a_idx,m_istart_sp)].get();

    // Compute saturation pressure
    MultiFab mf_sat_pressure(   m_mic_fab_vars[a_lev][MicVar_SD::pressure]->boxArray(),
                                m_mic_fab_vars[a_lev][MicVar_SD::pressure]->DistributionMap(),
                                1,
                                m_mic_fab_vars[a_lev][MicVar_SD::pressure]->nGrowVect() );
    species_mat.computeSaturationPressure( mf_sat_pressure,
                                           (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]) );
    mf_sat_pressure.FillBoundary();

    // Compute saturation ratio
    saturation_ratio( (*sr_ptr),
                      (*qv_ptr),
                      (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                      (*m_mic_fab_vars[a_lev][MicVar_SD::pressure]),
                      species_mat );

    // Recompute the species condensate from the current particles before qt.
    computeQcSpecies(a_idx, *a_z[a_lev]);

    // Compute total species content qt
    computeQtSpecies(a_idx);

    // Compute super-droplets mass change
    m_super_droplets->MassChange_LV ( a_lev,
                                      a_dt,
                                      m_species[a_idx],
                                      (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                                      (*m_mic_fab_vars[a_lev][MicVar_SD::pressure]),
                                      mf_sat_pressure,
                                      (*sr_ptr),
                                      a_z,
                                      is_water );

    // Compute new condensate mixing ratio
    computeQcSpecies(a_idx, *a_z[a_lev]);

    // Update vapour mixing ratio
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->const_array(mfi);
        auto qv_arr = qv_ptr->array(mfi);
        auto qc_arr = qc_ptr->array(mfi);

        auto theta_arr = m_mic_fab_vars[a_lev][MicVar_SD::theta]->array(mfi);
        auto T_arr = m_mic_fab_vars[a_lev][MicVar_SD::temperature]->const_array(mfi);
        auto L_by_Cp = species_mat.m_lat_vap / m_Cp;

        auto mask_arr = fine_mask.const_array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mask_arr(i,j,k) == 0) { return; }
            auto old_qv = qv_arr(i,j,k);
            if (qc_arr(i,j,k) > qt_arr(i,j,k)) {
                qv_arr(i,j,k) = 0.0;
                qc_arr(i,j,k) = qt_arr(i,j,k);
            } else {
                qv_arr(i,j,k) = qt_arr(i,j,k) - qc_arr(i,j,k);
            }
            AMREX_ALWAYS_ASSERT(qv_arr(i,j,k) >= 0.0);
            AMREX_ALWAYS_ASSERT(qc_arr(i,j,k) >= 0.0);

            auto theta_over_T = theta_arr(i,j,k)/T_arr(i,j,k);
            theta_arr(i,j,k) += theta_over_T * L_by_Cp * (old_qv-qv_arr(i,j,k));
        });

    }
}

/*! Freezing/melting (solid <--> liquid) for water: transfer mass between solid
 * and liquid phases for the super-droplet particles for a timestep, depending on
 * the ambient flow conditions (saturation ratio and temperature). Update the Eulerian
 * ice and cloud/rain mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange_SL_w ( const Real& a_dt,
                                            const Vector<std::unique_ptr<MultiFab>>& a_z,
                                            const int a_lev,
                                            const iMultiFab& fine_mask )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange_SL_w()");

    auto& species = m_species[m_idx_w];
    auto& species_mat = m_super_droplets->getSpeciesMaterial(species);
    bool is_water = species_mat.m_is_water;

    AMREX_ALWAYS_ASSERT(species == Species::Name::H2O);
    AMREX_ALWAYS_ASSERT(is_water);

    auto* qv_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_v].get();
    auto* qc_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_c].get();
    auto* qr_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_r].get();
    auto* qi_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_i].get();
    auto* qs_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_s].get();
    auto* qg_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_g].get();
    auto* qt_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_t].get();
    auto* sr_ptr = m_mic_fab_vars[a_lev][MicVar_SD::rh_w].get();

    // Compute saturation ratio
    saturation_ratio( (*sr_ptr),
                      (*qv_ptr),
                      (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                      (*m_mic_fab_vars[a_lev][MicVar_SD::pressure]),
                      species_mat );

    // Recompute the condensate from the current particles before forming qt,
    // so it reflects the particle-carried water/ice and not the dycore-advected
    // values left in the micro fabs by Copy_State_to_Micro.
    computeQcQrWater(*a_z[a_lev]);
    computeQiQgQsWater(*a_z[a_lev]);

    // Compute total liquid and ice/snow/graupel water content in qt
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->array(mfi);
        auto qc_arr = qc_ptr->const_array(mfi);
        auto qr_arr = qr_ptr->const_array(mfi);
        auto qi_arr = qi_ptr->const_array(mfi);
        auto qs_arr = qs_ptr->const_array(mfi);
        auto qg_arr = qg_ptr->const_array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            qt_arr(i,j,k) =   qc_arr(i,j,k) + qr_arr(i,j,k)
                            + qi_arr(i,j,k) + qs_arr(i,j,k) + qg_arr(i,j,k);
        });
    }

    // Moist air density for the melt-rate ventilation factor
    MultiFab mf_moist_density(  m_mic_fab_vars[a_lev][MicVar_SD::rho]->boxArray(),
                                m_mic_fab_vars[a_lev][MicVar_SD::rho]->DistributionMap(),
                                1,
                                m_mic_fab_vars[a_lev][MicVar_SD::rho]->nGrowVect() );
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());
        auto qv_arr = qv_ptr->const_array(mfi);
        auto rho_arr = m_mic_fab_vars[a_lev][MicVar_SD::rho]->const_array(mfi);
        auto rhom_arr = mf_moist_density.array(mfi);
        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     { rhom_arr(i,j,k) = (1.0+qv_arr(i,j,k)) * rho_arr(i,j,k); } );
    }
    mf_moist_density.FillBoundary();

    // Compute super-droplets mass change
    m_super_droplets->MassChange_SL ( a_lev,
                                      a_dt,
                                      (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                                      (*sr_ptr),
                                      mf_moist_density,
                                      a_z );

    // Compute new cloud/rain mixing ratios
    computeQcQrWater(*a_z[a_lev]);

    // Compute latent heat
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->const_array(mfi);
        auto qc_arr = qc_ptr->const_array(mfi);
        auto qr_arr = qr_ptr->const_array(mfi);
        auto qi_arr = qi_ptr->const_array(mfi);
        auto qs_arr = qs_ptr->const_array(mfi);
        auto qg_arr = qg_ptr->const_array(mfi);

        auto theta_arr = m_mic_fab_vars[a_lev][MicVar_SD::theta]->array(mfi);
        auto T_arr = m_mic_fab_vars[a_lev][MicVar_SD::temperature]->const_array(mfi);

        auto mask_arr = fine_mask.const_array(mfi);

        auto fac_cond = species_mat.m_lat_fus / m_Cp;

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mask_arr(i,j,k) == 0) { return; }
            auto q_liquid = qc_arr(i,j,k) + qr_arr(i,j,k);

            auto q_solid  = qi_arr(i,j,k) + qs_arr(i,j,k) + qg_arr(i,j,k);
            auto q_liquid_old = qt_arr(i,j,k) - q_solid;

            auto dq_liquid = q_liquid_old - q_liquid;

            auto theta_over_T = theta_arr(i,j,k)/T_arr(i,j,k);
            theta_arr(i,j,k) += theta_over_T * fac_cond * dq_liquid;
        });
    }

    // Now compute new ice/snow/graupel mixing ratio
    computeQiQgQsWater(*a_z[a_lev]);
}

/*! Deposition/sublimation (solid <--> vapour) for water (ice): shrink/grow the
 * super-droplet particles for a timestep, depending on the ambient flow conditions
 * (saturation ratio,  saturation pressure, and temperature). Update the Eulerian vapour
 * and ice/graupel mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange_SV_i ( const Real& a_dt,
                                            const Vector<std::unique_ptr<MultiFab>>& a_z,
                                            const int a_lev,
                                            const iMultiFab& fine_mask )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange_SV_i()");

    auto& species = m_species[m_idx_i];
    auto& species_mat = m_super_droplets->getSpeciesMaterial(species);
    AMREX_ALWAYS_ASSERT(species == Species::Name::ice);

    auto* qv_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_v].get();
    auto* qi_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_i].get();
    auto* qg_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_g].get();
    auto* qs_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_s].get();
    auto* qt_ptr = m_mic_fab_vars[a_lev][MicVar_SD::q_t].get();
    auto* sr_ptr = m_mic_fab_vars[a_lev][MicVar_SD::rh_i].get();

    // Compute moist density
    MultiFab mf_moist_density(  m_mic_fab_vars[a_lev][MicVar_SD::rho]->boxArray(),
                                m_mic_fab_vars[a_lev][MicVar_SD::rho]->DistributionMap(),
                                1,
                                m_mic_fab_vars[a_lev][MicVar_SD::rho]->nGrowVect() );
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qv_arr = qv_ptr->const_array(mfi);
        auto rho_arr = m_mic_fab_vars[a_lev][MicVar_SD::rho]->const_array(mfi);
        auto rhom_arr = mf_moist_density.array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     { rhom_arr(i,j,k) = (1.0+qv_arr(i,j,k)) * rho_arr(i,j,k); } );
    }
    mf_moist_density.FillBoundary();

    // Compute saturation pressure
    MultiFab mf_sat_pressure(   m_mic_fab_vars[a_lev][MicVar_SD::pressure]->boxArray(),
                                m_mic_fab_vars[a_lev][MicVar_SD::pressure]->DistributionMap(),
                                1,
                                m_mic_fab_vars[a_lev][MicVar_SD::pressure]->nGrowVect() );
    species_mat.computeSaturationPressure( mf_sat_pressure,
                                           (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]) );
    mf_sat_pressure.FillBoundary();

    // Compute ratio of saturation pressure for water to that of ice
    MultiFab mf_sat_pressure_wi(  m_mic_fab_vars[a_lev][MicVar_SD::pressure]->boxArray(),
                                  m_mic_fab_vars[a_lev][MicVar_SD::pressure]->DistributionMap(),
                                  1,
                                  m_mic_fab_vars[a_lev][MicVar_SD::pressure]->nGrowVect() );
    {
        auto& species_w_mat = m_super_droplets->getSpeciesMaterial(m_species[m_idx_w]);
        species_w_mat.computeSaturationPressure( mf_sat_pressure_wi,
                                               (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]) );
    }
    for (MFIter mfi(mf_sat_pressure_wi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        bx.grow( mf_sat_pressure_wi.nGrowVect() );
        auto esat_ratio_wi_arr = mf_sat_pressure_wi.array(mfi);
        auto esat_i_arr = mf_sat_pressure.array(mfi);
        ParallelFor( bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                     { esat_ratio_wi_arr(i,j,k) /= esat_i_arr(i,j,k); } );
    }
    mf_sat_pressure_wi.FillBoundary();

    // Compute saturation ratio
    saturation_ratio( (*sr_ptr),
                      (*qv_ptr),
                      (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                      (*m_mic_fab_vars[a_lev][MicVar_SD::pressure]),
                      species_mat );

    // Recompute ice/snow/graupel from the current particles so qt below is
    // built from the particle-carried condensate, not the dycore-advected
    // values left in the micro fabs by Copy_State_to_Micro.
    computeQiQgQsWater(*a_z[a_lev]);

    // Compute total liquid and vapour water content qt
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->array(mfi);
        auto qv_arr = qv_ptr->const_array(mfi);
        auto qi_arr = qi_ptr->const_array(mfi);
        auto qg_arr = qg_ptr->const_array(mfi);
        auto qs_arr = qs_ptr->const_array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { qt_arr(i,j,k) = qv_arr(i,j,k) + qi_arr(i,j,k) + qg_arr(i,j,k) + qs_arr(i,j,k); });
    }

    // Compute super-droplets mass change
    m_super_droplets->MassChange_SV ( a_lev,
                                      a_dt,
                                      (*m_mic_fab_vars[a_lev][MicVar_SD::rho]),
                                      (*m_mic_fab_vars[a_lev][MicVar_SD::temperature]),
                                      (*m_mic_fab_vars[a_lev][MicVar_SD::pressure]),
                                      mf_moist_density,
                                      mf_sat_pressure,
                                      mf_sat_pressure_wi,
                                      (*sr_ptr),
                                      a_z );

    // Compute new ice/graupel/snow mixing ratios
    computeQiQgQsWater(*a_z[a_lev]);

    // Update vapour mixing ratio
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->const_array(mfi);
        auto qv_arr = qv_ptr->array(mfi);
        auto qi_arr = qi_ptr->array(mfi);
        auto qg_arr = qg_ptr->array(mfi);
        auto qs_arr = qs_ptr->array(mfi);

        auto theta_arr = m_mic_fab_vars[a_lev][MicVar_SD::theta]->array(mfi);
        auto T_arr = m_mic_fab_vars[a_lev][MicVar_SD::temperature]->const_array(mfi);

        auto mask_arr = fine_mask.const_array(mfi);

        auto fac_cond = species_mat.m_lat_vap / m_Cp;

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mask_arr(i,j,k) == 0) { return; }
            auto old_qv = qv_arr(i,j,k);
            auto qw = qi_arr(i,j,k) + qg_arr(i,j,k) + qs_arr(i,j,k);
            if (qw > qt_arr(i,j,k)) {
                qv_arr(i,j,k) = 0.0;
                if (qg_arr(i,j,k) + qs_arr(i,j,k) > qt_arr(i,j,k)) {
                    qi_arr(i,j,k) = 0.0;
                    if (qs_arr(i,j,k) > qt_arr(i,j,k)) {
                        qg_arr(i,j,k) = 0.0;
                        qs_arr(i,j,k) = qt_arr(i,j,k);
                    } else {
                        qg_arr(i,j,k) = qt_arr(i,j,k) - qs_arr(i,j,k);
                    }
                } else {
                    qi_arr(i,j,k) = qt_arr(i,j,k) - (qg_arr(i,j,k) + qs_arr(i,j,k));
                }
            } else {
                qv_arr(i,j,k) = qt_arr(i,j,k) - qw;
            }

            AMREX_ALWAYS_ASSERT(qv_arr(i,j,k) >= 0.0);
            AMREX_ALWAYS_ASSERT(qi_arr(i,j,k) >= 0.0);
            AMREX_ALWAYS_ASSERT(qg_arr(i,j,k) >= 0.0);
            AMREX_ALWAYS_ASSERT(qs_arr(i,j,k) >= 0.0);

            auto theta_over_T = theta_arr(i,j,k)/T_arr(i,j,k);
            theta_arr(i,j,k) += theta_over_T * fac_cond * (old_qv-qv_arr(i,j,k));
        });
    }

}
#endif
