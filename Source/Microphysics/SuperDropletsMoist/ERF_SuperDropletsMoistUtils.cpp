#include "ERF_SuperDropletsMoist.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! \brief Ensure m_mic_fab_vars[lev] matches the geometry of a_cons_vars
 *
 * \param[in] a_lev AMR level
 * \param[in] a_cons_vars MultiFab providing target BoxArray and DistributionMap
 */
void SuperDropletsMoist::EnsureMicFabVars (const int a_lev, const MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::EnsureMicFabVars()");

    if (a_lev >= static_cast<int>(m_mic_fab_vars.size())) {
        m_mic_fab_vars.resize(a_lev + 1);
    }

    const int num_vars = MicVar_SD::NumVars
                        + (m_num_species-1) * MicVar_SD_Species::NumVars
                        + m_num_aerosols * MicVar_SD_Aerosols::NumVars;

    bool need_recreate = false;
    if (m_mic_fab_vars[a_lev].empty()) {
        need_recreate = true;
    } else if (static_cast<int>(m_mic_fab_vars[a_lev].size()) <= MicVar_SD::rho ||
               !m_mic_fab_vars[a_lev][MicVar_SD::rho] ||
               m_mic_fab_vars[a_lev][MicVar_SD::rho]->boxArray() != a_cons_vars.boxArray() ||
               m_mic_fab_vars[a_lev][MicVar_SD::rho]->DistributionMap() != a_cons_vars.DistributionMap()) {
        need_recreate = true;
    }

    if (need_recreate) {
        m_mic_fab_vars[a_lev].resize(num_vars);
        for (int i = 0; i < num_vars; i++) {
            m_mic_fab_vars[a_lev][i] = std::make_shared<MultiFab>(
                a_cons_vars.boxArray(), a_cons_vars.DistributionMap(),
                1, a_cons_vars.nGrowVect());
            m_mic_fab_vars[a_lev][i]->setVal(0.0);
        }
    }
}

/*! \brief Initialize m_mic_fab_vars for a new AMR level
 *
 * \param[in] a_lev AMR level
 * \param[in] a_cons_vars Conserved variables MultiFab
 */
void SuperDropletsMoist::InitLevel (const int a_lev, const MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::InitLevel()");
    EnsureMicFabVars(a_lev, a_cons_vars);
}

/*! \brief Build a coarse-cell mask: 1 where the coarse cell is exposed (not
 *  covered by a finer level), 0 where covered. If no finer level exists or the
 *  fine BoxArray is empty, the returned mask is all 1s.
 *
 * \param[in] a_lev AMR level whose grid is used as the coarse layout
 */
iMultiFab SuperDropletsMoist::buildFineMask (const int a_lev) const
{
    BL_PROFILE("SuperDropletsMoist::buildFineMask()");

    AMREX_ALWAYS_ASSERT(a_lev >= 0 && a_lev < static_cast<int>(m_mic_fab_vars.size()));
    const auto& cba = m_mic_fab_vars[a_lev][MicVar_SD::rho]->boxArray();
    const auto& cdm = m_mic_fab_vars[a_lev][MicVar_SD::rho]->DistributionMap();
    // Match the mic fabs' ghost width so the mask is valid over the same grown
    // region the phaseChange / Copy_Micro_to_State kernels iterate.
    const IntVect ng = m_mic_fab_vars[a_lev][MicVar_SD::rho]->nGrowVect();

    const int next_lev = a_lev + 1;
    const bool has_fine = (m_super_droplets &&
                           next_lev < m_super_droplets->numLevels() &&
                           !m_super_droplets->ParticleBoxArray(next_lev).empty());
    if (!has_fine) {
        iMultiFab mask(cba, cdm, 1, ng);
        mask.setVal(1);
        return mask;
    }

    const auto& fba = m_super_droplets->ParticleBoxArray(next_lev);
    const IntVect rr = m_super_droplets->GetParGDB()->refRatio(a_lev);
    const Periodicity period = m_super_droplets->Geom(a_lev).periodicity();
    return makeFineMask(cba, cdm, ng, fba, rr, period, 1, 0);
}

/*! \brief Copy moisture model variables from conserved state to member MultiFabs
 *
 * \param[in] a_cons_vars MultiFab containing the conserved state variables
 */
void SuperDropletsMoist::Copy_State_to_Micro (const MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::Copy_State_to_Micro()");

    const int lev = m_current_lev;
    EnsureMicFabVars(lev, a_cons_vars);

    const auto& gvec = a_cons_vars.nGrowVect();

    // Read all moisture (qv, qc, qr) from state to keep the micro vars
    // consistent with dycore-advected values; avoids q_t inconsistency on
    // AMR coarse/fine boundaries since q_c/q_r in micro fabs would otherwise
    // lag the dycore by one step.
    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        bx.grow(gvec);
        auto states_arr = a_cons_vars.const_array(mfi);

        // state variables
        {
            auto rho_arr = m_mic_fab_vars[lev][MicVar_SD::rho]->array(mfi);
            auto theta_arr = m_mic_fab_vars[lev][MicVar_SD::theta]->array(mfi);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                rho_arr(i,j,k) = states_arr(i,j,k,Rho_comp);
                theta_arr(i,j,k) = states_arr(i,j,k,RhoTheta_comp)/states_arr(i,j,k,Rho_comp);
            });
        }

        // water
        {
            auto q_t_arr = m_mic_fab_vars[lev][MicVar_SD::q_t]->array(mfi);
            auto q_v_arr = m_mic_fab_vars[lev][MicVar_SD::q_v]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[lev][MicVar_SD::q_c]->array(mfi);
            auto q_r_arr = m_mic_fab_vars[lev][MicVar_SD::q_r]->array(mfi);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real inv_rho = Real(1.0) / states_arr(i,j,k,Rho_comp);
                q_v_arr(i,j,k) = states_arr(i,j,k,RhoQ1_comp) * inv_rho;
                q_c_arr(i,j,k) = states_arr(i,j,k,RhoQ2_comp) * inv_rho;
                q_r_arr(i,j,k) = states_arr(i,j,k,RhoQ3_comp) * inv_rho;
                q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k) + q_r_arr(i,j,k);
            });
        }

        // other species
        for (int is = 1; is < m_num_species; is++) {
            auto q_v_arr = m_mic_fab_vars[lev][s_qv_idx(is)]->array(mfi);
            auto q_t_arr = m_mic_fab_vars[lev][s_qt_idx(is)]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[lev][s_qc_idx(is)]->const_array(mfi);
            auto qv_comp = q_qv_idx(is);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                q_v_arr(i,j,k) = states_arr(i,j,k,qv_comp) / states_arr(i,j,k,Rho_comp);
                q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k);
            });
        }
    }

    // Compute pressure and temperature
    for (MFIter mfi(*m_mic_fab_vars[lev][MicVar_SD::temperature], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        const Array4<Real      >& t_arr  = m_mic_fab_vars[lev][MicVar_SD::temperature]->array(mfi);
        const Array4<Real      >& p_arr  = m_mic_fab_vars[lev][MicVar_SD::pressure]->array(mfi);
        const Array4<Real const>& S_arr  = a_cons_vars.const_array(mfi);
        const Array4<Real const>& qv_arr = m_mic_fab_vars[lev][MicVar_SD::q_v]->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            t_arr(i,j,k,0) = getTgivenRandRTh(S_arr(i,j,k,Rho_comp),S_arr(i,j,k,RhoTheta_comp),qv_arr(i,j,k,0));
            p_arr(i,j,k,0) = getPgivenRTh(S_arr(i,j,k,RhoTheta_comp),qv_arr(i,j,k,0));
        });
    }

    AMREX_ASSERT( !m_mic_fab_vars[lev][MicVar_SD::pressure]->contains_nan() );
    AMREX_ASSERT( !m_mic_fab_vars[lev][MicVar_SD::temperature]->contains_nan() );

    // water
    {
        // Get vapour material properties object for water
        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(Species::Name::H2O);

        // Compute saturation ratio
        vapour_mat.computeSaturationVapFrac( (*m_mic_fab_vars[lev][MicVar_SD::rh]),
                                             (*m_mic_fab_vars[lev][MicVar_SD::temperature]),
                                             (*m_mic_fab_vars[lev][MicVar_SD::pressure]) );

        for (   MFIter mfi((*m_mic_fab_vars[lev][MicVar_SD::rh]),
                TilingIfNotGPU()); mfi.isValid();
                ++mfi ) {

            Box bx = mfi.tilebox();
            bx.grow( m_mic_fab_vars[lev][MicVar_SD::rh]->nGrowVect() );

            const Array4<Real>& sr_arr = m_mic_fab_vars[lev][MicVar_SD::rh]->array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[lev][MicVar_SD::q_v]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });

        }
    }

    // other species
    for (int is = 1; is < m_num_species; is++) {
        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(m_species[is]);
        vapour_mat.computeSaturationVapFrac((*m_mic_fab_vars[lev][s_sr_idx(is)]),
                                            (*m_mic_fab_vars[lev][MicVar_SD::temperature]),
                                            (*m_mic_fab_vars[lev][MicVar_SD::pressure]) );
        for (   MFIter mfi((*m_mic_fab_vars[lev][s_sr_idx(is)]),
                TilingIfNotGPU()); mfi.isValid();
                ++mfi ) {

            Box bx = mfi.tilebox();
            bx.grow( m_mic_fab_vars[lev][s_sr_idx(is)]->nGrowVect() );

            const Array4<Real>& sr_arr = m_mic_fab_vars[lev][s_sr_idx(is)]->array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[lev][s_qv_idx(is)]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });
        }
    }

    const Geometry& geom = (m_super_droplets && lev < m_super_droplets->numLevels())
                         ? m_super_droplets->Geom(lev)
                         : m_geom;
    for (auto i(0); i < MicVar_SD::NumVars; i++) {
        m_mic_fab_vars[lev][i]->FillBoundary(geom.periodicity());
    }
}

/*! \brief Copy moisture model variables to the conserved state vector
 *
 * \param[out] a_cons_vars MultiFab containing conserved state variables to be updated
 */
void SuperDropletsMoist::Copy_Micro_to_State (MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::Copy_Micro_to_state()");

    const int lev = m_current_lev;
    const auto& gvec = a_cons_vars.nGrowVect();

    iMultiFab fine_mask = buildFineMask(lev);

    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto states_arr = a_cons_vars.array(mfi);
        auto mask_arr = fine_mask.const_array(mfi);

        // state variables
        {
            auto theta_arr = m_mic_fab_vars[lev][MicVar_SD::theta]->const_array(mfi);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (mask_arr(i,j,k) == 0) return;
                states_arr(i,j,k,RhoTheta_comp) = states_arr(i,j,k,Rho_comp)*theta_arr(i,j,k);
            });
        }

        // water
        {
            auto q_v_arr = m_mic_fab_vars[lev][MicVar_SD::q_v]->const_array(mfi);
            auto q_c_arr = m_mic_fab_vars[lev][MicVar_SD::q_c]->const_array(mfi);
            auto q_r_arr = m_mic_fab_vars[lev][MicVar_SD::q_r]->const_array(mfi);

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (mask_arr(i,j,k) == 0) return;
                states_arr(i,j,k,RhoQ1_comp) = states_arr(i,j,k,Rho_comp)*q_v_arr(i,j,k);
                states_arr(i,j,k,RhoQ2_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
                states_arr(i,j,k,RhoQ3_comp) = states_arr(i,j,k,Rho_comp)*q_r_arr(i,j,k);
            });
        }

        // other species
        for (int is = 1; is < m_num_species; is++) {
            auto q_v_arr = m_mic_fab_vars[lev][s_qv_idx(is)]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[lev][s_qc_idx(is)]->array(mfi);
            auto qv_comp = q_qv_idx(is);
            auto qc_comp = q_qc_idx(is);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (mask_arr(i,j,k) == 0) return;
                states_arr(i,j,k,qv_comp) = states_arr(i,j,k,Rho_comp)*q_v_arr(i,j,k);
                states_arr(i,j,k,qc_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
            });
        }
    }

    const Geometry& geom = (m_super_droplets && lev < m_super_droplets->numLevels())
                         ? m_super_droplets->Geom(lev)
                         : m_geom;
    a_cons_vars.FillBoundary(geom.periodicity());
}

/*! \brief Update microphysics variables from conserved state */
void SuperDropletsMoist::Update_Micro_Vars (MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::Update_Micro_Vars()");
    Copy_State_to_Micro(a_cons_vars);
}

/*! \brief Average down moisture multifabs from the finest level down to 0
 *
 * For each level pair (lev, lev-1) from finest_level down to 1, restrict
 * fine-level micro multifabs onto coarse cells covered by the fine level.
 * This keeps coarse-level micro variables consistent with the particle deposit
 * on finer levels for tagging, plotting, and pre-regrid synchronization.
 *
 * \param[in] finest_level Top AMR level to restrict from
 */
void SuperDropletsMoist::AverageDownMicroVars (const int finest_level)
{
    BL_PROFILE("SuperDropletsMoist::AverageDownMicroVars()");

    for (int lev = finest_level; lev >= 1; --lev) {
        if (lev >= static_cast<int>(m_mic_fab_vars.size())) continue;
        if (m_mic_fab_vars[lev].empty() || m_mic_fab_vars[lev-1].empty()) continue;

        const IntVect rr = (m_super_droplets && lev <= m_super_droplets->numLevels())
                         ? m_super_droplets->GetParGDB()->refRatio(lev-1)
                         : IntVect(2);

        const int n = std::min<int>(static_cast<int>(m_mic_fab_vars[lev].size()),
                                    static_cast<int>(m_mic_fab_vars[lev-1].size()));
        for (int v = 0; v < n; ++v) {
            if (!m_mic_fab_vars[lev][v] || !m_mic_fab_vars[lev-1][v]) continue;
            average_down( *m_mic_fab_vars[lev][v],
                          *m_mic_fab_vars[lev-1][v],
                          0,
                          m_mic_fab_vars[lev][v]->nComp(),
                          rr );
        }
    }
}

/*! \brief Compute derived quantities and update state variables */
void SuperDropletsMoist::Update_State_Vars (MultiFab& a_cons_vars, const MultiFab& a_z_phys_nd)
{
    BL_PROFILE("SuperDropletsMoist::Update_State_Vars()");
    computeQcQrWater(a_z_phys_nd);
    computeQtWater();
    rainAccumulation(a_z_phys_nd);

    computeQcSpecies(a_z_phys_nd);
    computeQtSpecies();
    speciesAccumulation(a_z_phys_nd);
    aerosolAccumulation(a_z_phys_nd);

    if (!m_kinematic_mode) { Copy_Micro_to_State(a_cons_vars); }
}

/*! \brief Convert a density field to a mixing ratio field */
void SuperDropletsMoist::densityToRatio (MultiFab& a_var,
                                         const int a_comp)
{
    BL_PROFILE("SuperDropletsMoist::densityToRatio()");

    const int lev = m_current_lev;
    const auto& gvec = a_var.nGrowVect();

    for ( MFIter mfi(a_var); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto rho_arr = m_mic_fab_vars[lev][MicVar_SD::rho]->const_array(mfi);
        auto fab_arr = a_var.array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { fab_arr(i,j,k,a_comp) /= rho_arr(i,j,k); });
    }

    const Geometry& geom = (m_super_droplets && lev < m_super_droplets->numLevels())
                         ? m_super_droplets->Geom(lev)
                         : m_geom;
    a_var.FillBoundary(geom.periodicity());
}

/*! \brief Convert a mixing ratio field to a density field */
void SuperDropletsMoist::ratioToDensity (MultiFab& a_var,
                                         const int a_comp)
{
    BL_PROFILE("SuperDropletsMoist::ratioToDensity()");

    const int lev = m_current_lev;
    const auto& gvec = a_var.nGrowVect();

    for ( MFIter mfi(a_var); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto rho_arr = m_mic_fab_vars[lev][MicVar_SD::rho]->const_array(mfi);
        auto fab_arr = a_var.array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { fab_arr(i,j,k,a_comp) *= rho_arr(i,j,k); });
    }

    const Geometry& geom = (m_super_droplets && lev < m_super_droplets->numLevels())
                         ? m_super_droplets->Geom(lev)
                         : m_geom;
    a_var.FillBoundary(geom.periodicity());
}

/*! \brief Compute cloud and rain water mixing ratios from superdroplets */
void SuperDropletsMoist::computeQcQrWater (const MultiFab& a_z_phys_nd)
{
    BL_PROFILE("SuperDropletsMoist::computeQcQrWater()");

    const int lev = m_current_lev;

    m_super_droplets->cloudRainDensity( *(m_mic_fab_vars[lev][MicVar_SD::q_c]),
                                        a_z_phys_nd,
                                        lev,
                                        0,
                                        m_r_rain );
    m_super_droplets->cloudRainDensity( *(m_mic_fab_vars[lev][MicVar_SD::q_r]),
                                        a_z_phys_nd,
                                        lev,
                                        m_r_rain,
                                        one );

    if (m_dimensionality == SDMSimulationDim::one_d_z) {
        for ( MFIter mfi(*m_mic_fab_vars[lev][MicVar_SD::q_c]); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            int imin = bx.smallEnd(0);
            int jmin = bx.smallEnd(1);
            auto q_c_arr = m_mic_fab_vars[lev][MicVar_SD::q_c]->array(mfi);
            auto q_r_arr = m_mic_fab_vars[lev][MicVar_SD::q_r]->array(mfi);

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                q_c_arr(i,j,k) = q_c_arr(imin,jmin,k);
                q_r_arr(i,j,k) = q_r_arr(imin,jmin,k);
            });
        }
    }

    densityToRatio(*(m_mic_fab_vars[lev][MicVar_SD::q_c]));
    densityToRatio(*(m_mic_fab_vars[lev][MicVar_SD::q_r]));
}

/*! compute qt (total) for water */
void SuperDropletsMoist::computeQtWater ()
{
    BL_PROFILE("SuperDropletsMoist::computeQtWater()");

    const int lev = m_current_lev;

    for ( MFIter mfi(*m_mic_fab_vars[lev][MicVar_SD::q_t]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(m_mic_fab_vars[lev][MicVar_SD::q_t]->nGrowVect());

        auto q_c_arr = m_mic_fab_vars[lev][MicVar_SD::q_c]->const_array(mfi);
        auto q_r_arr = m_mic_fab_vars[lev][MicVar_SD::q_r]->const_array(mfi);
        auto q_v_arr = m_mic_fab_vars[lev][MicVar_SD::q_v]->const_array(mfi);
        auto q_t_arr = m_mic_fab_vars[lev][MicVar_SD::q_t]->array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k) + q_r_arr(i,j,k); });
    }
}

/*! Compute rain accumulation */
void SuperDropletsMoist::rainAccumulation (const MultiFab& a_z_phys_nd)
{
    BL_PROFILE("SuperDropletsMoist::rainAccumulation()");

    const int lev = m_current_lev;
    const Geometry& geom = (m_super_droplets && lev < m_super_droplets->numLevels())
                         ? m_super_droplets->Geom(lev)
                         : m_geom;
    auto domain = geom.Domain();
    int k_lo = domain.smallEnd(2);
    auto dt = m_dt;

    auto& vapour_mat = m_super_droplets->getSpeciesMaterial(Species::Name::H2O);
    auto mat_density = vapour_mat.m_density;

    MultiFab mf_zflux( m_mic_fab_vars[lev][MicVar_SD::rain_accum]->boxArray(),
                       m_mic_fab_vars[lev][MicVar_SD::rain_accum]->DistributionMap(),
                       1,
                       m_mic_fab_vars[lev][MicVar_SD::rain_accum]->nGrowVect() );
    m_super_droplets->speciesMassFlux(mf_zflux, a_z_phys_nd, lev, m_idx_w, 2);

    for ( MFIter mfi((*m_mic_fab_vars[lev][MicVar_SD::rain_accum]),TilingIfNotGPU());
          mfi.isValid(); ++mfi ) {
        Box bx = mfi.tilebox();
        const Array4<Real const>& zflux_arr = mf_zflux.const_array(mfi);
        const Array4<Real>& rain_accum_arr = m_mic_fab_vars[lev][MicVar_SD::rain_accum]->array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (k == k_lo) {
                auto rain_accum = std::max(Real(0), -zflux_arr(i,j,k)*dt/mat_density);
                rain_accum_arr(i,j,k) += (rain_accum * Real(1000.0) /* [m] -> [mm] */);
            }
        });
    }

}

/*! compute condensate mixing ratio */
void SuperDropletsMoist::computeQcSpecies (const int a_i, const MultiFab& a_z_phys_nd)
{
    BL_PROFILE("SuperDropletsMoist::computeQcSpecies()");

    const int lev = m_current_lev;

    m_super_droplets->speciesMassDensity( *(m_mic_fab_vars[lev][s_qc_idx(a_i)]),
                                          a_z_phys_nd,
                                          lev,
                                          a_i );
    if (m_dimensionality == SDMSimulationDim::one_d_z) {
        for ( MFIter mfi(*m_mic_fab_vars[lev][s_qc_idx(a_i)]); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            int imin = bx.smallEnd(0);
            int jmin = bx.smallEnd(1);
            auto q_c_arr = m_mic_fab_vars[lev][s_qc_idx(a_i)]->array(mfi);

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            { q_c_arr(i,j,k) = q_c_arr(imin,jmin,k); });
        }
    }
    densityToRatio(*(m_mic_fab_vars[lev][s_qc_idx(a_i)]));
}

/*! compute qt (total) */
void SuperDropletsMoist::computeQtSpecies (const int a_i)
{
    BL_PROFILE("SuperDropletsMoist::computeQtSpecies()");

    const int lev = m_current_lev;

    for ( MFIter mfi(*m_mic_fab_vars[lev][s_qt_idx(a_i)]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(m_mic_fab_vars[lev][MicVar_SD::q_t]->nGrowVect());

        auto q_c_arr = m_mic_fab_vars[lev][s_qc_idx(a_i)]->const_array(mfi);
        auto q_v_arr = m_mic_fab_vars[lev][s_qv_idx(a_i)]->const_array(mfi);
        auto q_t_arr = m_mic_fab_vars[lev][s_qt_idx(a_i)]->array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k); });
    }
}

/*! Compute ground accumulation for non-water species */
void SuperDropletsMoist::speciesAccumulation (const MultiFab& a_z_phys_nd)
{
    BL_PROFILE("SuperDropletsMoist::speciesAccumulation()");

    const int lev = m_current_lev;
    const Geometry& geom = (m_super_droplets && lev < m_super_droplets->numLevels())
                         ? m_super_droplets->Geom(lev)
                         : m_geom;
    auto domain = geom.Domain();
    const auto dx = geom.CellSizeArray();
    int k_lo = domain.smallEnd(2);
    auto dt = m_dt;

    for (int is = 1; is < m_num_species; is++) {
        MultiFab mf_zflux( m_mic_fab_vars[lev][s_sr_idx(is)]->boxArray(),
                           m_mic_fab_vars[lev][s_sr_idx(is)]->DistributionMap(),
                           1,
                           m_mic_fab_vars[lev][s_sr_idx(is)]->nGrowVect() );
        m_super_droplets->speciesMassFlux(mf_zflux, a_z_phys_nd, lev, is, 2);

        for ( MFIter mfi((*m_mic_fab_vars[lev][MicVar_SD::rain_accum]),TilingIfNotGPU());
              mfi.isValid(); ++mfi ) {
            Box bx = mfi.tilebox();
            const Array4<Real const>& zflux_arr = mf_zflux.const_array(mfi);
            const Array4<Real>& accum_arr = m_mic_fab_vars[lev][s_accum_idx(is)]->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (k == k_lo) {
                    auto accum = std::max(Real(0), -zflux_arr(i,j,k)*dt*dx[0]*dx[1]);
                    accum_arr(i,j,k) += accum;
                }
            });
        }
    }
}

/*! Compute ground accumulation for aerosol species */
void SuperDropletsMoist::aerosolAccumulation (const MultiFab& a_z_phys_nd)
{
    BL_PROFILE("SuperDropletsMoist::aerosolAccumulation()");

    const int lev = m_current_lev;
    const Geometry& geom = (m_super_droplets && lev < m_super_droplets->numLevels())
                         ? m_super_droplets->Geom(lev)
                         : m_geom;
    auto domain = geom.Domain();
    const auto dx = geom.CellSizeArray();
    int k_lo = domain.smallEnd(2);
    auto dt = m_dt;

    for (int ia = 0; ia < m_num_aerosols; ia++) {
        MultiFab mf_zflux( m_mic_fab_vars[lev][MicVar_SD::rain_accum]->boxArray(),
                           m_mic_fab_vars[lev][MicVar_SD::rain_accum]->DistributionMap(),
                           1,
                           m_mic_fab_vars[lev][MicVar_SD::rain_accum]->nGrowVect() );
        m_super_droplets->aerosolMassFlux(mf_zflux, a_z_phys_nd, lev, ia, 2);

        for ( MFIter mfi((*m_mic_fab_vars[lev][MicVar_SD::rain_accum]),TilingIfNotGPU());
              mfi.isValid(); ++mfi ) {
            Box bx = mfi.tilebox();
            const Array4<Real const>& zflux_arr = mf_zflux.const_array(mfi);
            const Array4<Real>& accum_arr = m_mic_fab_vars[lev][a_accum_idx(m_num_species,ia)]->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (k == k_lo) {
                    auto accum = std::max(Real(0), -zflux_arr(i,j,k)*dt*dx[0]*dx[1]);
                    accum_arr(i,j,k) += accum;
                }
            });
        }
    }
}

#endif
