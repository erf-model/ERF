#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_TableData.H>
#include <AMReX_GpuContainers.H>

#include <ERF_NumericalDiffusion.H>
#include <ERF_SrcHeaders.H>
#include <ERF_TI_slow_headers.H>
#include <ERF_MOSTStress.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in]  level level of resolution
 * @param[in]  nrk   which RK stage
 * @param[in]  dt    slow time step
 * @param[in]  S_data current solution
 * @param[in]  S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[in] source source terms for conserved variables
 * @param[in]  geom   Container for geometric information
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in] mapfac map factors
 * @param[in] dptr_rhotheta_src  custom temperature source term
 * @param[in] dptr_rhoqt_src  custom moisture source term
 * @param[in] dptr_wbar_sub  subsidence source term
 * @param[in] d_rayleigh_ptrs_at_lev  Vector of {strength of Rayleigh damping, reference value of theta} used to define Rayleigh damping
 * @param[in] d_sinesq_at_lev  sin( (pi/2) (z-z_t)/(damping depth)) at cell centers
 */

void make_sources (int level,
                   int /*nrk*/,
                   Real dt,
                   Real time,
                   const Vector<MultiFab>& S_data,
                   const  MultiFab & S_prim,
                          MultiFab & source,
                   const  MultiFab & base_state,
                   const  MultiFab*  z_phys_cc,
                   const  MultiFab & xvel,
                   const  MultiFab & yvel,
                   const MultiFab* qheating_rates,
                          MultiFab* terrain_blank,
                   const Geometry geom,
                   const SolverChoice& solverChoice,
                   Vector<std::unique_ptr<MultiFab>>& mapfac,
                   const Real* dptr_rhotheta_src,
                   const Real* dptr_rhoqt_src,
                   const Real* dptr_wbar_sub,
                   const Vector<Real*> d_rayleigh_ptrs_at_lev,
                   const Real* d_sinesq_at_lev,
                   InputSoundingData& input_sounding_data,
                   TurbulentPerturbation& turbPert,
                   bool is_slow_step)
{
    BL_PROFILE_REGION("erf_make_sources()");

    // *****************************************************************************
    // Initialize source to zero since we re-compute it every RK stage
    // *****************************************************************************
    if (is_slow_step) {
        source.setVal(0.);
    } else {
        source.setVal(0.0,Rho_comp,2);
    }

    const bool l_use_ndiff      = solverChoice.use_num_diff;

    TurbChoice tc = solverChoice.turbChoice[level];
    const bool l_use_KE  = tc.use_tke;
    const bool l_diff_KE = tc.diffuse_tke_3D;

    const Box& domain = geom.Domain();

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> dx    = geom.CellSizeArray();

    MultiFab r_hse (base_state, make_alias, BaseState::r0_comp , 1);

    Real* thetabar = d_rayleigh_ptrs_at_lev[Rayleigh::thetabar];

    // flags to apply certain source terms in substep call only
    bool use_Rayleigh_fast = ( (solverChoice.dampingChoice.rayleigh_damping_type == RayleighDampingType::FastExplicit) ||
                               (solverChoice.dampingChoice.rayleigh_damping_type == RayleighDampingType::FastImplicit) );
    bool use_ImmersedForcing_fast = solverChoice.immersed_forcing_substep;

    // flag for a moisture model
    bool has_moisture = (solverChoice.moisture_type != MoistureType::None);

    // *****************************************************************************
    // Planar averages for subsidence terms
    // *****************************************************************************
    Table1D<Real>      dptr_r_plane, dptr_t_plane, dptr_qv_plane, dptr_qc_plane;
    TableData<Real, 1>  r_plane_tab,  t_plane_tab,  qv_plane_tab,  qc_plane_tab;
    bool compute_averages = false;
    compute_averages = compute_averages ||
        ( is_slow_step && (dptr_wbar_sub || solverChoice.nudging_from_input_sounding) );

    if (compute_averages)
    {
        //
        // The call to "compute_averages" currently does all the components in one call
        // We can then extract each component separately with the "line_average" call
        //
        // We need just one ghost cell in the vertical
        //
        IntVect ng_c(S_data[IntVars::cons].nGrowVect()); ng_c[2] = 1;
        //
        // With no moisture we only (rho) and (rho theta); with moisture we also do qv and qc
        // We use the alias here to control ncomp inside the PlaneAverage
        //
        int ncomp = (!has_moisture) ? 2 : RhoQ2_comp+1;
        MultiFab cons(S_data[IntVars::cons], make_alias, 0, ncomp);

        PlaneAverage cons_ave(&cons, geom, solverChoice.ave_plane, ng_c);
        cons_ave.compute_averages(ZDir(), cons_ave.field());

        int ncell = cons_ave.ncell_line();

        Gpu::HostVector<    Real> r_plane_h(ncell);
        Gpu::DeviceVector<  Real> r_plane_d(ncell);

        Gpu::HostVector<    Real> t_plane_h(ncell);
        Gpu::DeviceVector<  Real> t_plane_d(ncell);

        cons_ave.line_average(Rho_comp     , r_plane_h);
        cons_ave.line_average(RhoTheta_comp, t_plane_h);

        Gpu::copyAsync(Gpu::hostToDevice, r_plane_h.begin(), r_plane_h.end(), r_plane_d.begin());
        Gpu::copyAsync(Gpu::hostToDevice, t_plane_h.begin(), t_plane_h.end(), t_plane_d.begin());

        Real* dptr_r = r_plane_d.data();
        Real* dptr_t = t_plane_d.data();

        Box tdomain  = domain; tdomain.grow(2,ng_c[2]);
        r_plane_tab.resize({tdomain.smallEnd(2)}, {tdomain.bigEnd(2)});
        t_plane_tab.resize({tdomain.smallEnd(2)}, {tdomain.bigEnd(2)});

        int offset = ng_c[2];

        dptr_r_plane = r_plane_tab.table();
        dptr_t_plane = t_plane_tab.table();
        ParallelFor(ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
        {
            dptr_r_plane(k-offset) = dptr_r[k];
            dptr_t_plane(k-offset) = dptr_t[k];
        });

        if (has_moisture)
        {
            Gpu::HostVector<  Real> qv_plane_h(ncell), qc_plane_h(ncell);
            Gpu::DeviceVector<Real> qv_plane_d(ncell), qc_plane_d(ncell);

            // Water vapor
            cons_ave.line_average(RhoQ1_comp, qv_plane_h);
            Gpu::copyAsync(Gpu::hostToDevice, qv_plane_h.begin(), qv_plane_h.end(), qv_plane_d.begin());

            // Cloud water
            cons_ave.line_average(RhoQ2_comp, qc_plane_h);
            Gpu::copyAsync(Gpu::hostToDevice, qc_plane_h.begin(), qc_plane_h.end(), qc_plane_d.begin());

            Real* dptr_qv = qv_plane_d.data();
            Real* dptr_qc = qc_plane_d.data();

            qv_plane_tab.resize({tdomain.smallEnd(2)}, {tdomain.bigEnd(2)});
            qc_plane_tab.resize({tdomain.smallEnd(2)}, {tdomain.bigEnd(2)});

            dptr_qv_plane = qv_plane_tab.table();
            dptr_qc_plane = qc_plane_tab.table();
            ParallelFor(ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
            {
                dptr_qv_plane(k-offset) = dptr_qv[k];
                dptr_qc_plane(k-offset) = dptr_qc[k];
            });
        }
    }

    // *****************************************************************************
    // Radiation flux vector for four stream approximation
    // *****************************************************************************
    // NOTE: The fluxes live on w-faces
    int klo = domain.smallEnd(0);
    int khi = domain.bigEnd(2);
    int nk  = khi - klo + 2;
    Gpu::DeviceVector<Real> radiation_flux(nk,0.0);
    Gpu::DeviceVector<Real> q_integral(nk,0.0);
    Real* rad_flux = radiation_flux.data();
    Real* q_int    = q_integral.data();

    // *****************************************************************************
    // Define source term for cell-centered conserved variables, from
    //    1. user-defined source terms for (rho theta) and (rho q_t)
    //    2. radiation           for (rho theta)
    //    3. Rayleigh damping    for (rho theta)
    //    4. custom forcing      for (rho theta) and (rho Q1)
    //    5. custom subsidence   for (rho theta) and (rho Q1)
    //    6. numerical diffusion for (rho theta)
    //    7. sponging
    //    8. turbulent perturbation
    //    9. nudging towards input sounding values (only for theta)
    //   10a. Immersed forcing for terrain
    //   10b. Immersed forcing for buildings
    //   11. Four stream radiation source for (rho theta)
    // *****************************************************************************

    // ***********************************************************************************************
    // Add remaining source terms
    // ***********************************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
    for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box bx  = mfi.tilebox();

        const Array4<const Real>& cell_data  = S_data[IntVars::cons].array(mfi);
        const Array4<const Real>& cell_prim  = S_prim.array(mfi);
        const Array4<Real>      & cell_src   = source.array(mfi);

        const Array4<const Real>& r0 = r_hse.const_array(mfi);

        const Array4<const Real>& z_cc_arr = z_phys_cc->const_array(mfi);

        const Array4<const Real>& t_blank_arr = (terrain_blank) ? terrain_blank->const_array(mfi) :
                                                               Array4<const Real>{};


        // *************************************************************************************
        // 2. Add radiation source terms to (rho theta)
        // *************************************************************************************
        if (solverChoice.rad_type != RadiationType::None && is_slow_step) {
            auto const& qheating_arr = qheating_rates->const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Short-wavelength and long-wavelength radiation source terms
                cell_src(i,j,k,RhoTheta_comp) += cell_data(i,j,k,Rho_comp) * ( qheating_arr(i,j,k,0) + qheating_arr(i,j,k,1) );
            });
        }


        // *************************************************************************************
        // 3. Add Rayleigh damping for (rho theta)
        // *************************************************************************************
        Real dampcoef = solverChoice.dampingChoice.rayleigh_dampcoef;

        if (solverChoice.dampingChoice.rayleigh_damp_T) {
            if ((is_slow_step && !use_Rayleigh_fast) || (!is_slow_step && use_Rayleigh_fast)) {
                int n  = RhoTheta_comp;
                int nr = Rho_comp;
                int np = PrimTheta_comp;

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real theta = cell_prim(i,j,k,np);
                    Real sinesq = d_sinesq_at_lev[k];
                    cell_src(i, j, k, n) -= dampcoef*sinesq * (theta - thetabar[k]) * cell_data(i,j,k,nr);
                });
            }
        }

        // *************************************************************************************
        // 4. Add custom forcing for (rho theta)
        // *************************************************************************************
        if (solverChoice.custom_rhotheta_forcing && is_slow_step) {
            const int n = RhoTheta_comp;
            if (solverChoice.custom_forcing_prim_vars) {
                const int nr = Rho_comp;
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_src(i, j, k, n) += cell_data(i,j,k,nr) * dptr_rhotheta_src[k];
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_src(i, j, k, n) += dptr_rhotheta_src[k];
                });
            }
        }

        // *************************************************************************************
        // 4. Add custom forcing for RhoQ1
        // *************************************************************************************
        if (solverChoice.custom_moisture_forcing && is_slow_step) {
            const int n = RhoQ1_comp;
            if (solverChoice.custom_forcing_prim_vars) {
                const int nr = Rho_comp;
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_src(i, j, k, n) += cell_data(i,j,k,nr) * dptr_rhoqt_src[k];
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_src(i, j, k, n) += dptr_rhoqt_src[k];
                });
            }
        }

        // *************************************************************************************
        // 5. Add custom subsidence for (rho theta)
        // *************************************************************************************
        if (solverChoice.custom_w_subsidence && is_slow_step) {
            const int n = RhoTheta_comp;
            if (solverChoice.custom_forcing_prim_vars) {
                const int nr = Rho_comp;
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real dzInv = (z_cc_arr) ? 1.0/ (z_cc_arr(i,j,k+1) - z_cc_arr(i,j,k-1)) : 0.5*dxInv[2];
                    Real T_hi = dptr_t_plane(k+1) / dptr_r_plane(k+1);
                    Real T_lo = dptr_t_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_cc = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    cell_src(i, j, k, n) -= cell_data(i,j,k,nr) * wbar_cc * (T_hi - T_lo) * dzInv;
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real dzInv = (z_cc_arr) ? 1.0/ (z_cc_arr(i,j,k+1) - z_cc_arr(i,j,k-1)) : 0.5*dxInv[2];
                    Real T_hi = dptr_t_plane(k+1) / dptr_r_plane(k+1);
                    Real T_lo = dptr_t_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_cc = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    cell_src(i, j, k, n) -= wbar_cc * (T_hi - T_lo) * dzInv;
                });
            }
        }

        // *************************************************************************************
        // 5. Add custom subsidence for RhoQ1 and RhoQ2
        // *************************************************************************************
        if (solverChoice.custom_w_subsidence && (solverChoice.moisture_type != MoistureType::None) && is_slow_step) {
            const int nv = RhoQ1_comp;
            if (solverChoice.custom_forcing_prim_vars) {
                const int nr = Rho_comp;
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real dzInv = (z_cc_arr) ? 1.0/ (z_cc_arr(i,j,k+1) - z_cc_arr(i,j,k-1)) : 0.5*dxInv[2];
                    Real Qv_hi = dptr_qv_plane(k+1) / dptr_r_plane(k+1);
                    Real Qv_lo = dptr_qv_plane(k-1) / dptr_r_plane(k-1);
                    Real Qc_hi = dptr_qc_plane(k+1) / dptr_r_plane(k+1);
                    Real Qc_lo = dptr_qc_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_cc = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    cell_src(i, j, k, nv  ) -= cell_data(i,j,k,nr) * wbar_cc * (Qv_hi - Qv_lo) * dzInv;
                    cell_src(i, j, k, nv+1) -= cell_data(i,j,k,nr) * wbar_cc * (Qc_hi - Qc_lo) * dzInv;
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real dzInv = (z_cc_arr) ? 1.0/ (z_cc_arr(i,j,k+1) - z_cc_arr(i,j,k-1)) : 0.5*dxInv[2];
                    Real Qv_hi = dptr_qv_plane(k+1) / dptr_r_plane(k+1);
                    Real Qv_lo = dptr_qv_plane(k-1) / dptr_r_plane(k-1);
                    Real Qc_hi = dptr_qc_plane(k+1) / dptr_r_plane(k+1);
                    Real Qc_lo = dptr_qc_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_cc = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    cell_src(i, j, k, nv  ) -= wbar_cc * (Qv_hi - Qv_lo) * dzInv;
                    cell_src(i, j, k, nv+1) -= wbar_cc * (Qc_hi - Qc_lo) * dzInv;
                });
            }
        }

        // *************************************************************************************
        // 6. Add numerical diffuion for rho and (rho theta)
        // *************************************************************************************
        if (l_use_ndiff && is_slow_step) {
            int sc;
            int nc;

            const Array4<const Real>& mf_mx   = mapfac[MapFacType::m_x]->const_array(mfi);
            const Array4<const Real>& mf_my   = mapfac[MapFacType::m_y]->const_array(mfi);

            // Rho is a special case
            NumericalDiffusion_Scal(bx, sc=0, nc=1, dt, solverChoice.num_diff_coeff,
                                    cell_data, cell_data, cell_src, mf_mx, mf_my);

            // Other scalars proceed as normal
            NumericalDiffusion_Scal(bx, sc=1, nc=1, dt, solverChoice.num_diff_coeff,
                                    cell_prim, cell_data, cell_src, mf_mx, mf_my);


            if (l_use_KE && l_diff_KE) {
                NumericalDiffusion_Scal(bx, sc=RhoKE_comp, nc=1, dt, solverChoice.num_diff_coeff,
                                        cell_prim, cell_data, cell_src, mf_mx, mf_my);
            }

            NumericalDiffusion_Scal(bx, sc=RhoScalar_comp, nc=NSCALARS, dt, solverChoice.num_diff_coeff,
                                    cell_prim, cell_data, cell_src, mf_mx, mf_my);
        }

        // *************************************************************************************
        // 7. Add sponging
        // *************************************************************************************
        if(!(solverChoice.spongeChoice.sponge_type == "input_sponge") && is_slow_step){
            ApplySpongeZoneBCsForCC(solverChoice.spongeChoice, geom, bx, cell_src, cell_data, r0, z_cc_arr);
        }

        // *************************************************************************************
        // 8. Add perturbation
        // *************************************************************************************
        if (solverChoice.pert_type == PerturbationType::Source && is_slow_step) {
            auto m_ixtype = S_data[IntVars::cons].boxArray().ixType(); // Conserved term
            const amrex::Array4<const amrex::Real>& pert_cell = turbPert.pb_cell[level].const_array(mfi);
            turbPert.apply_tpi(level, bx, RhoTheta_comp, m_ixtype, cell_src, pert_cell); // Applied as source term
        }

        // *************************************************************************************
        // 9. Add nudging towards value specified in input sounding
        // *************************************************************************************
        if (solverChoice.nudging_from_input_sounding && is_slow_step)
        {
            int itime_n    = 0;
            int itime_np1  = 0;
            Real coeff_n   = Real(1.0);
            Real coeff_np1 = Real(0.0);

            Real tau_inv = Real(1.0) / input_sounding_data.tau_nudging;

            int n_sounding_times = input_sounding_data.input_sounding_time.size();

            for (int nt = 1; nt < n_sounding_times; nt++) {
                if (time > input_sounding_data.input_sounding_time[nt]) itime_n = nt;
            }
            if (itime_n == n_sounding_times-1) {
                itime_np1 = itime_n;
            } else {
                itime_np1 = itime_n+1;
                coeff_np1 = (time                                               - input_sounding_data.input_sounding_time[itime_n]) /
                            (input_sounding_data.input_sounding_time[itime_np1] - input_sounding_data.input_sounding_time[itime_n]);
                coeff_n   = Real(1.0) - coeff_np1;
            }

            const Real* theta_inp_sound_n   = input_sounding_data.theta_inp_sound_d[itime_n].dataPtr();
            const Real* theta_inp_sound_np1 = input_sounding_data.theta_inp_sound_d[itime_np1].dataPtr();

            const int n  = RhoTheta_comp;
            const int nr = Rho_comp;

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real nudge = (coeff_n*theta_inp_sound_n[k] + coeff_np1*theta_inp_sound_np1[k]) - (dptr_t_plane(k)/dptr_r_plane(k));
                nudge *= tau_inv;
                cell_src(i, j, k, n) += cell_data(i, j, k, nr) * nudge;
            });
        }

        // *************************************************************************************
        // 10a. Add immersed source terms for terrain
        // *************************************************************************************
        if (solverChoice.terrain_type == TerrainType::ImmersedForcing &&
           ((is_slow_step && !use_ImmersedForcing_fast) || (!is_slow_step && use_ImmersedForcing_fast)))
        {
            const Array4<const Real>& u = xvel.array(mfi);
            const Array4<const Real>& v = yvel.array(mfi);

            // geometric properties
            const Real* dx_arr = geom.CellSize();
            const Real dx_x = dx_arr[0];
            const Real dx_y = dx_arr[1];
            const Real dx_z = dx_arr[2];

            const Real alpha_h          = solverChoice.if_Cd_scalar;
            const Real drag_coefficient = alpha_h / std::pow(dx_x*dx_y*dx_z, 1./3.);
            const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
            const Real U_s              = 1.0; // unit velocity scale

            // MOST parameters
            similarity_funs sfuns;
            const Real ggg                = CONST_GRAV;
            const Real kappa              = KAPPA;
            const Real z0                 = solverChoice.if_z0;
            const Real tflux              = solverChoice.if_surf_temp_flux;
            const Real init_surf_temp     = solverChoice.if_init_surf_temp;
            const Real surf_heating_rate  = solverChoice.if_surf_heating_rate;
            const Real Olen_in            = solverChoice.if_Olen_in;

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real t_blank       = t_blank_arr(i, j, k);
                const Real t_blank_above = t_blank_arr(i, j, k+1);
                const Real ux_cc_2r = 0.5 * (u(i  ,j  ,k+1) + u(i+1,j  ,k+1));
                const Real uy_cc_2r = 0.5 * (v(i  ,j  ,k+1) + v(i  ,j+1,k+1));
                const Real h_windspeed2r  = std::sqrt(ux_cc_2r * ux_cc_2r + uy_cc_2r * uy_cc_2r);

                const Real theta          = cell_data(i,j,k  ,RhoTheta_comp) / cell_data(i,j,k  ,Rho_comp);
                const Real theta_neighbor = cell_data(i,j,k+1,RhoTheta_comp) / cell_data(i,j,k+1,Rho_comp);

                // SURFACE TEMP AND HEATING/COOLING RATE
                if (init_surf_temp > 0.0) {
                    if (t_blank > 0 && (t_blank_above == 0.0)) { // force to MOST value
                        const Real surf_temp    = init_surf_temp + surf_heating_rate*time/3600;
                        const Real bc_forcing_rt_srf = -(cell_data(i,j,k-1,Rho_comp) * surf_temp - cell_data(i,j,k-1,RhoTheta_comp));
                        cell_src(i, j, k-1, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf; // k-1
                    }
                }

                // SURFACE HEAT FLUX
                if (tflux != 1e-8){
                    if (t_blank > 0 && (t_blank_above == 0.0)) { // force to MOST value
                        Real psi_m           = 0.0;
                        Real psi_h           = 0.0;
                        Real psi_h_neighbor  = 0.0;
                        Real ustar = h_windspeed2r * kappa / (std::log((1.5) * dx_z / z0) - psi_m);
                        const Real Olen  = -ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny);
                        const Real zeta          = (0.5) * dx_z / Olen;
                        const Real zeta_neighbor = (1.5) * dx_z / Olen;

                        // similarity functions
                        psi_m          = sfuns.calc_psi_m(zeta);
                        psi_h          = sfuns.calc_psi_h(zeta);
                        psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                        ustar = h_windspeed2r * kappa / (std::log((1.5) * dx_z / z0) - psi_m);

                        // prevent some unphysical math
                        if (!(ustar > 0.0 && !std::isnan(ustar))) { ustar = 0.0; }
                        if (!(ustar < 2.0 && !std::isnan(ustar))) { ustar = 2.0; }
                        if (psi_h_neighbor > std::log(1.5 * dx_z / z0)) { psi_h_neighbor = std::log(1.5 * dx_z / z0); }
                        if (psi_h > std::log(0.5 * dx_z / z0)) { psi_h = std::log(0.5 * dx_z / z0); }

                        // We do not know the actual temperature so use cell above
                        const Real thetastar    = theta * ustar * ustar / (kappa * ggg * Olen);
                        const Real surf_temp    = theta_neighbor - thetastar / kappa * (std::log((1.5) * dx_z / z0) - psi_h_neighbor);
                        const Real tTarget      = surf_temp + thetastar / kappa * (std::log((0.5) * dx_z / z0) - psi_h);

                        const Real bc_forcing_rt = -(cell_data(i,j,k,Rho_comp) * tTarget - cell_data(i,j,k,RhoTheta_comp));
                        cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
                    }
                }

                // OBUKHOV LENGTH
                if (Olen_in != 1e-8){
                    if (t_blank > 0 && (t_blank_above == 0.0)) { // force to MOST value
                        const Real Olen  = Olen_in;
                        const Real zeta          = (0.5) * dx_z / Olen;
                        const Real zeta_neighbor = (1.5) * dx_z / Olen;

                        // similarity functions
                        const Real psi_m          = sfuns.calc_psi_m(zeta);
                        const Real psi_h          = sfuns.calc_psi_h(zeta);
                        const Real psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                        const Real ustar = h_windspeed2r * kappa / (std::log((1.5) * dx_z / z0) - psi_m);

                        // We do not know the actual temperature so use cell above
                        const Real thetastar    = theta * ustar * ustar / (kappa * ggg * Olen);
                        const Real surf_temp    = theta_neighbor - thetastar / kappa * (std::log((1.5) * dx_z / z0) - psi_h_neighbor);
                        const Real tTarget      = surf_temp + thetastar / kappa * (std::log((0.5) * dx_z / z0) - psi_h);

                        const Real bc_forcing_rt = -(cell_data(i,j,k,Rho_comp) * tTarget - cell_data(i,j,k,RhoTheta_comp));
                        cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
                    }
                }

            });
        }

        // *************************************************************************************
        // 10b. Add immersed source terms for buildings
        // *************************************************************************************
        if ((solverChoice.buildings_type == BuildingsType::ImmersedForcing ) &&
           ((is_slow_step && !use_ImmersedForcing_fast) || (!is_slow_step && use_ImmersedForcing_fast)))
        {
            // geometric properties
            const Real* dx_arr = geom.CellSize();
            const Real dx_x = dx_arr[0];
            const Real dx_y = dx_arr[1];

            const Real alpha_h          = solverChoice.if_Cd_scalar;
            const Real U_s              = 1.0; // unit velocity scale
            const Real min_t_blank      = 0.005;

            const Real init_surf_temp     = solverChoice.if_init_surf_temp;
            const Real surf_heating_rate  = solverChoice.if_surf_heating_rate;

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real t_blank       = t_blank_arr(i, j, k);
                Real t_blank_below = t_blank_arr(i, j, k-1);
                Real t_blank_above = t_blank_arr(i, j, k+1);
                Real t_blank_north  = t_blank_arr(i  , j+1, k);
                Real t_blank_south  = t_blank_arr(i  , j-1, k);
                Real t_blank_east   = t_blank_arr(i+1, j  , k);
                Real t_blank_west   = t_blank_arr(i-1, j  , k);
                if (t_blank < min_t_blank) { t_blank = 0.0; } // deal with situations where very small volfrac exist
                if (t_blank_below < min_t_blank) { t_blank_below = 0.0; }
                if (t_blank_north < min_t_blank) { t_blank_north = 0.0; }
                if (t_blank_south < min_t_blank) { t_blank_south = 0.0; }
                if (t_blank_east < min_t_blank) { t_blank_east = 0.0; }
                if (t_blank_west < min_t_blank) { t_blank_west = 0.0; }

                Real dx_z    = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx[2];
                Real drag_coefficient = alpha_h / std::pow(dx_x*dx_y*dx_z, 1./3.);

                // SURFACE TEMP AND HEATING/COOLING RATE
                if (init_surf_temp > 0.0) {
                    const Real surf_temp    = init_surf_temp + surf_heating_rate*time/3600;
                    if (t_blank > 0 && (t_blank_above == 0.0) && (t_blank_below == 1.0)) { // building roof
                        const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                        cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;

                    } else if (((t_blank > 0 && t_blank < t_blank_west && t_blank_east == 0.0) ||
                                (t_blank > 0 && t_blank < t_blank_east && t_blank_west == 0.0) ||
                                (t_blank > 0 && t_blank < t_blank_north && t_blank_south == 0.0) ||
                                (t_blank > 0 && t_blank < t_blank_south && t_blank_north == 0.0))) {
                        // this should enter for just building walls
                        // walls are currently separated to allow for flexible in the future to heat walls differently

                        // south face
                        if ((t_blank < t_blank_north) && (t_blank_north == 1.0)) {
                            const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                            cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                        }

                        // north face
                        if ((t_blank < t_blank_south) && (t_blank_south == 1.0)) {
                            const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                            cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                        }

                        // west face
                        if ((t_blank < t_blank_east) && (t_blank_east == 1.0)) {
                            const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                            cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                        }

                        // east face
                        if ((t_blank < t_blank_west) && (t_blank_west == 1.0)) {
                            const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                            cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                        }

                    }
                }
            });
        }

        // *************************************************************************************
        // 11. Add 4 stream radiation src to RhoTheta
        // *************************************************************************************
        if (solverChoice.four_stream_radiation && has_moisture && is_slow_step)
        {
            AMREX_ALWAYS_ASSERT((bx.smallEnd(2) == klo) && (bx.bigEnd(2) == khi));
            Real D    = 3.75e-6; // [s^-1]
            Real F0   = 70;      // [W/m^2]
            Real F1   = 22;      // [W/m^2]
            Real krad = 85;      // [m^2 kg^-1]
            Real qt_i = 0.008;

            Box xybx = makeSlab(bx,2,klo);
            ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
            {
                // Inclusive scan at w-faces for the Q integral (also find "i" values)
                q_int[0] = 0.0;
                Real zi   = 0.5 * (z_cc_arr(i,j,khi) + z_cc_arr(i,j,khi-1));
                Real rhoi = 0.5 * (cell_data(i,j,khi,Rho_comp) + cell_data(i,j,khi-1,Rho_comp));
                for (int k(klo+1); k<=khi+1; ++k) {
                    int lk    = k - klo;
                    // Average to w-faces when looping w-faces
                    Real dz    = (z_cc_arr) ? 0.5 * (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-2)) : dx[2];
                    q_int[lk]  = q_int[lk-1] + krad * cell_data(i,j,k-1,Rho_comp) * cell_data(i,j,k-1,RhoQ2_comp) * dz;
                    Real qt_hi = cell_data(i,j,k  ,RhoQ1_comp) + cell_data(i,j,k  ,RhoQ2_comp);
                    Real qt_lo = cell_data(i,j,k-1,RhoQ1_comp) + cell_data(i,j,k-1,RhoQ2_comp);
                    if ( (qt_lo > qt_i) && (qt_hi < qt_i) ) {
                        zi   = 0.5 * (z_cc_arr(i,j,k) + z_cc_arr(i,j,k-1));
                        rhoi = 0.5 * (cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp));
                    }
                }

                // Decompose the integral to get the fluxes at w-faces
                Real q_int_inf = q_int[khi+1];
                for (int k(klo); k<=khi+1; ++k) {
                    int lk       = k - klo;
                    Real z       = 0.5 * (z_cc_arr(i,j,k) + z_cc_arr(i,j,k-1));
                    rad_flux[lk] = F1*std::exp(-q_int[lk]) + F0*std::exp(-(q_int_inf - q_int[lk]));
                    if (z > zi) {
                      rad_flux[lk] += rhoi * Cp_d * D * ( std::pow(z-zi,4./3.)/4. + zi*std::pow(z-zi,1./3.) ) ;
                    }
                }

                // Compute the radiative heating source
                for (int k(klo); k<=khi; ++k) {
                    int lk       = k - klo;
                    // Average to w-faces when looping CC
                    Real dzInv   = (z_cc_arr) ? 1.0/ (0.5 * (z_cc_arr(i,j,k+1) - z_cc_arr(i,j,k-1))) : dxInv[2];
                    // NOTE: Fnet  = Up - Dn (all fluxes are up here)
                    //       dT/dt = dF/dz * (1/(-rho*Cp))
                    Real dTdt    = (rad_flux[lk+1] - rad_flux[lk]) * dzInv / (-cell_data(i,j,k,Rho_comp)*Cp_d);
                    Real qv      = cell_data(i,j,k,RhoQ1_comp)/cell_data(i,j,k,Rho_comp);
                    Real iexner  = 1./getExnergivenRTh(cell_data(i,j,k,RhoTheta_comp), R_d/Cp_d, qv);
                    // Convert dT/dt to dTheta/dt and multiply rho
                    cell_src(i,j,k,RhoTheta_comp) += cell_data(i,j,k,Rho_comp) * dTdt * iexner;
                }
            });
        }
    } // mfi
    } // OMP
}
