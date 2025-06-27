#include "ERF_SuperDropletsMoist.H"
#include "ERF_EOS.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Compute saturation ratio */
static void saturation_ratio (MultiFab& a_sr, /*!< Saturation ratio */
                              const MultiFab& a_qv, /*!< Vapour mixing ratio */
                              const MultiFab& a_temperature, /*!< Temperature */
                              const MultiFab& a_pressure, /*!< Pressure */
                              const MaterialProperties& a_species_mat /*!< Species material properties */ )
{
    a_species_mat.computeSaturationVapFrac(a_sr, a_temperature, a_pressure);
    for (MFIter mfi(a_sr, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        bx.grow( a_sr.nGrowVect() );
        const Array4<Real>& sr_arr = a_sr.array(mfi);
        const Array4<Real const>& qv_arr = a_qv.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });
    }
    a_sr.FillBoundary();
}

/*! Compute phase change for a timestep: shrink/grow the super-droplet particles
    for a timestep, depending on the ambient flow conditions (saturation ratio,
    saturation pressure, and temperature). Update the Eulerial vapour and condensate
    mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange ( const Real& a_dt, /*!< Timestep */
                                       const Vector<std::unique_ptr<MultiFab>>& a_z /*!< terrain */ )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange()");
    auto dt_s = a_dt / static_cast<Real>(m_num_substeps_phase_change);

    for (int substep = 0; substep < m_num_substeps_phase_change; substep++) {

        // freezing/melting (solid <--> liquid) -  water
        phaseChange_SL_w(dt_s, a_z);

        // evaporation and condensation (vapour <--> liquid) - water
        phaseChange_LV_w(dt_s, a_z);

        // evaporation and condensation (vapour <--> liquid) - other species
        for (int is = m_istart_sp; is < m_num_species; is++) {
            phaseChange_LV_s(is, dt_s, a_z);
        }

        // Update pressure and temperature
        const auto& gvec = (*m_mic_fab_vars[MicVar_SD::temperature]).nGrowVect();
        for (MFIter mfi(*m_mic_fab_vars[MicVar_SD::temperature], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            Box bx = mfi.tilebox();
            bx.grow(gvec);

            const Array4<Real>& t_arr = m_mic_fab_vars[MicVar_SD::temperature]->array(mfi);
            const Array4<Real>& p_arr = m_mic_fab_vars[MicVar_SD::pressure]->array(mfi);

            const Array4<Real const>& rho_arr = m_mic_fab_vars[MicVar_SD::rho]->const_array(mfi);
            const Array4<Real const>& theta_arr = m_mic_fab_vars[MicVar_SD::theta]->const_array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);

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
            if (is == m_idx_i) { continue; } // skip ice

            auto& species = m_species[is];
            auto& species_mat = m_super_droplets->getSpeciesMaterial(species);
            AMREX_ALWAYS_ASSERT(species!=Species::Name::ice);

            MultiFab *qv_ptr, *sr_ptr;
            if (is == m_idx_w) {
                qv_ptr = m_mic_fab_vars[MicVar_SD::q_v].get();
                sr_ptr = m_mic_fab_vars[MicVar_SD::rh].get();
            } else {
                qv_ptr = m_mic_fab_vars[s_qv_idx(is,m_istart_sp)].get();
                sr_ptr = m_mic_fab_vars[s_sr_idx(is,m_istart_sp)].get();
            }

            saturation_ratio( (*sr_ptr),
                              (*qv_ptr),
                              (*m_mic_fab_vars[MicVar_SD::temperature]),
                              (*m_mic_fab_vars[MicVar_SD::pressure]),
                              species_mat );

        }

    }
}

/*! Condensation/evaporation (liquid <--> vapour) for water: shrink/grow the
 * super-droplet particles for a timestep, depending on the ambient flow conditions
 * (saturation ratio,  saturation pressure, and temperature). Update the Eulerian vapour
 * and condensate mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange_LV_w ( const Real& a_dt, /*!< Timestep */
                                            const Vector<std::unique_ptr<MultiFab>>& a_z /*!< terrain */ )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange_LV_w()");

    auto& species = m_species[m_idx_w];
    auto& species_mat = m_super_droplets->getSpeciesMaterial(species);
    bool is_water = species_mat.m_is_water;

    AMREX_ALWAYS_ASSERT(species == Species::Name::H2O);
    AMREX_ALWAYS_ASSERT(is_water);

    auto* qv_ptr = m_mic_fab_vars[MicVar_SD::q_v].get();
    auto* qc_ptr = m_mic_fab_vars[MicVar_SD::q_c].get();
    auto* qr_ptr = m_mic_fab_vars[MicVar_SD::q_r].get();
    auto* qt_ptr = m_mic_fab_vars[MicVar_SD::q_t].get();
    auto* sr_ptr = m_mic_fab_vars[MicVar_SD::rh].get();

    // Compute saturation pressure
    MultiFab mf_sat_pressure(   m_mic_fab_vars[MicVar_SD::pressure]->boxArray(),
                                m_mic_fab_vars[MicVar_SD::pressure]->DistributionMap(),
                                1,
                                m_mic_fab_vars[MicVar_SD::pressure]->nGrowVect() );
    species_mat.computeSaturationPressure( mf_sat_pressure,
                                           (*m_mic_fab_vars[MicVar_SD::temperature]) );
    mf_sat_pressure.FillBoundary();

    // Compute saturation ratio
    saturation_ratio( (*sr_ptr),
                      (*qv_ptr),
                      (*m_mic_fab_vars[MicVar_SD::temperature]),
                      (*m_mic_fab_vars[MicVar_SD::pressure]),
                      species_mat );

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
    m_super_droplets->MassChange_LV ( 0,
                                      a_dt,
                                      m_species[m_idx_w],
                                      (*m_mic_fab_vars[MicVar_SD::temperature]),
                                      (*m_mic_fab_vars[MicVar_SD::pressure]),
                                      mf_sat_pressure,
                                      (*sr_ptr),
                                      a_z,
                                      is_water );

    // Compute new condensate mixing ratio
    computeQcQrWater();

    // Update vapour mixing ratio
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->const_array(mfi);
        auto qv_arr = qv_ptr->array(mfi);
        auto qc_arr = qc_ptr->array(mfi);
        auto qr_arr = qr_ptr->array(mfi);

        auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->array(mfi);
        auto T_arr = m_mic_fab_vars[MicVar_SD::temperature]->const_array(mfi);
        auto dqc_arr = m_mic_fab_vars[MicVar_SD::dqcdt]->array(mfi);

        auto fac_cond = species_mat.m_lat_vap / m_Cp;

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
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
void SuperDropletsMoist::phaseChange_LV_s ( const int a_idx, /*!< Species index */
                                            const Real& a_dt, /*!< Timestep */
                                            const Vector<std::unique_ptr<MultiFab>>& a_z /*!< terrain */ )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange_LV_s()");

    auto& species = m_species[a_idx];
    auto& species_mat = m_super_droplets->getSpeciesMaterial(species);
    bool is_water = species_mat.m_is_water; // may be true for dummy water species used for debugging

    AMREX_ALWAYS_ASSERT((species!=Species::Name::ice) || (species!=Species::Name::H2O));

    auto* qv_ptr = m_mic_fab_vars[s_qv_idx(a_idx,m_istart_sp)].get();
    auto* qc_ptr = m_mic_fab_vars[s_qc_idx(a_idx,m_istart_sp)].get();
    auto* qt_ptr = m_mic_fab_vars[s_qt_idx(a_idx,m_istart_sp)].get();
    auto* sr_ptr = m_mic_fab_vars[s_sr_idx(a_idx,m_istart_sp)].get();

    // Compute saturation pressure
    MultiFab mf_sat_pressure(   m_mic_fab_vars[MicVar_SD::pressure]->boxArray(),
                                m_mic_fab_vars[MicVar_SD::pressure]->DistributionMap(),
                                1,
                                m_mic_fab_vars[MicVar_SD::pressure]->nGrowVect() );
    species_mat.computeSaturationPressure( mf_sat_pressure,
                                           (*m_mic_fab_vars[MicVar_SD::temperature]) );
    mf_sat_pressure.FillBoundary();

    // Compute saturation ratio
    saturation_ratio( (*sr_ptr),
                      (*qv_ptr),
                      (*m_mic_fab_vars[MicVar_SD::temperature]),
                      (*m_mic_fab_vars[MicVar_SD::pressure]),
                      species_mat );

    // Compute total species content qt
    computeQtSpecies(a_idx);

    // Compute super-droplets mass change
    m_super_droplets->MassChange_LV ( 0,
                                      a_dt,
                                      m_species[a_idx],
                                      (*m_mic_fab_vars[MicVar_SD::temperature]),
                                      (*m_mic_fab_vars[MicVar_SD::pressure]),
                                      mf_sat_pressure,
                                      (*sr_ptr),
                                      a_z,
                                      is_water );

    // Compute new condensate mixing ratio
    computeQcSpecies(a_idx);

    // Update vapour mixing ratio
    for ( MFIter mfi(*qv_ptr); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(qv_ptr->nGrowVect());

        auto qt_arr = qt_ptr->const_array(mfi);
        auto qv_arr = qv_ptr->array(mfi);
        auto qc_arr = qc_ptr->array(mfi);

        auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->array(mfi);
        auto T_arr = m_mic_fab_vars[MicVar_SD::temperature]->const_array(mfi);
        auto L_by_Cp = species_mat.m_lat_vap / m_Cp;

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
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
void SuperDropletsMoist::phaseChange_SL_w ( const Real& a_dt, /*!< Timestep */
                                            const Vector<std::unique_ptr<MultiFab>>& a_z /*!< terrain */ )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange_SL_w()");

    auto& species = m_species[m_idx_w];
    auto& species_mat = m_super_droplets->getSpeciesMaterial(species);
    bool is_water = species_mat.m_is_water;

    AMREX_ALWAYS_ASSERT(species == Species::Name::H2O);
    AMREX_ALWAYS_ASSERT(is_water);

    auto* qv_ptr = m_mic_fab_vars[MicVar_SD::q_v].get();
    auto* qc_ptr = m_mic_fab_vars[MicVar_SD::q_c].get();
    auto* qr_ptr = m_mic_fab_vars[MicVar_SD::q_r].get();
    auto* qi_ptr = m_mic_fab_vars[MicVar_SD::q_i].get();
    auto* qs_ptr = m_mic_fab_vars[MicVar_SD::q_s].get();
    auto* qg_ptr = m_mic_fab_vars[MicVar_SD::q_g].get();
    auto* qt_ptr = m_mic_fab_vars[MicVar_SD::q_t].get();
    auto* sr_ptr = m_mic_fab_vars[MicVar_SD::rh].get();

    // Compute saturation ratio
    saturation_ratio( (*sr_ptr),
                      (*qv_ptr),
                      (*m_mic_fab_vars[MicVar_SD::temperature]),
                      (*m_mic_fab_vars[MicVar_SD::pressure]),
                      species_mat );

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

    // Compute super-droplets mass change
    m_super_droplets->MassChange_SL ( 0,
                                      a_dt,
                                      (*m_mic_fab_vars[MicVar_SD::temperature]),
                                      (*sr_ptr),
                                      a_z );

    // Compute new cloud/rain mixing ratios
    computeQcQrWater();

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

        auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->array(mfi);
        auto T_arr = m_mic_fab_vars[MicVar_SD::temperature]->const_array(mfi);

        auto fac_cond = species_mat.m_lat_fus / m_Cp;

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            auto q_liquid = qc_arr(i,j,k) + qr_arr(i,j,k);

            auto q_solid  = qi_arr(i,j,k) + qs_arr(i,j,k) + qg_arr(i,j,k);
            auto q_liquid_old = qt_arr(i,j,k) - q_solid;

            auto dq_liquid = q_liquid_old - q_liquid;

            auto theta_over_T = theta_arr(i,j,k)/T_arr(i,j,k);
            theta_arr(i,j,k) += theta_over_T * fac_cond * dq_liquid;
        });
    }

    // Now compute new ice/snow/graupel mixing ratio
    computeQiQgQsWater();
}

#endif
