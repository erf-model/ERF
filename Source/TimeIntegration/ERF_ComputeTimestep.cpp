#include <AMReX_Reduce.H>
#include "ERF_Constants.H"

#include <ERF_EOS.H>
#include <ERF_Utils.H>
#include <ERF_TimestepUtils.H>
#include <ERF_TerrainMetrics.H>
#include "Diffusion/ERF_TerrainDiffusionLimits.H"
#include <ERF.H>
#include "Diffusion/ERF_CloudChamberWallFlux.H"
#include "TimeIntegration/ERF_CloudChamberWallDtGuard.H"
#include "AuxiliaryState/ERF_AuxiliaryMappedTransport.H"
#include "Microphysics/SBM/ERF_SBMTransport.H"

#include <cmath>
#include <limits>
#include <sstream>

using namespace amrex;

/**
 * Function that calls estTimeStep for each level
 *
 */
void
ERF::ComputeDt (int step, double cur_time_d)
{
    Vector<double> dt_tmp(finest_level+1);

    // Explicit eddy-diffusion rates per level (1/s), from the diffusivities of the last step
    const bool do_diffusive_check = (diffusive_dt_check || diffusive_dt_limit);
    Vector<Real> rate_mom (finest_level+1, zero);
    Vector<Real> rate_scal(finest_level+1, zero);

    // Track exactly the current levels, so a level removed and later re-created is reported again
    slope_report_grids.resize(finest_level+1);
    slope_report_alpha.resize(finest_level+1, zero);

    const bool fitted = (solverChoice.terrain_type == TerrainType::StaticFittedMesh ||
                         solverChoice.terrain_type == TerrainType::MovingFittedMesh);
    const bool moving = (solverChoice.terrain_type == TerrainType::MovingFittedMesh);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // The terrain slope factor is reported only when one of the options it informs is on:
        // the diffusive check or limit, or a Smagorinsky2D limit at this level.  It is printed
        // whenever a level has a new set of grids (at initialisation, on the first step after a
        // restart, and on the coarse step after a regrid), and on a moving terrain also whenever
        // its largest value changes by more than 1 %.
        const TurbChoice& tc = solverChoice.turbChoice[lev];
        const bool slope_report = fitted && tc.use_kturb && tc.kh_kv_can_differ() &&
            (do_diffusive_check || tc.smag2d_slope_limiter || tc.smag2d_kh_cap > zero);
        if (slope_report) {
            const bool new_grids = (slope_report_grids[lev] != grids[lev]);
            if (SlopeReportEvaluate(new_grids, moving)) {
                ReportTerrainSlopeFactor(lev, new_grids);
                slope_report_grids[lev] = grids[lev];
            }
        }

        dt_tmp[lev] = estTimeStep(lev, dt_mri_ratio[lev]);

        if (do_diffusive_check) {
            ComputeDiffusiveRates(lev, rate_mom[lev], rate_scal[lev]);
            const Real rate = amrex::max(rate_mom[lev], rate_scal[lev]);
            if (diffusive_dt_limit && fixed_dt[lev] <= zero) {
                if (rate > zero) {
                    const double dt_diff = static_cast<double>(diffusive_cfl) / static_cast<double>(rate);
                    if (verbose) {
                        Print() << "Diffusive dt at level " << lev << ":  " << dt_diff
                                << " (erf.diffusive_cfl = " << diffusive_cfl << ")" << std::endl;
                    }
                    dt_tmp[lev] = std::min(dt_tmp[lev], dt_diff);
                } else if (diffusive_first_call && !restart_chkfile.empty() && istep[lev] > 0) {
                    // The eddy diffusivities are not checkpointed, so on the first step after a
                    // restart they are not known yet: do not let dt grow past the checkpointed
                    // step (which the limit bounded if the run that wrote it used the limit).
                    dt_tmp[lev] = std::min(dt_tmp[lev], dt[lev]);
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    double dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = amrex::min(dt_tmp[lev], static_cast<double>(change_max*dt[lev]));
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, static_cast<double>(n_factor*dt_tmp[lev]));

    }
    // Limit level 0 time step if requested
    if (step == 0) {
        dt_0 *= init_shrink;
        if (verbose && init_shrink != one) {
            Print() << "Timestep 0: shrink level 0 initial dt by " << init_shrink << std::endl;
        }
    }
    //
    // Limit dt by the value of stop_time.
    // Recall that stop_time is total time, but t_new is elapsed time,
    //     so we must add start_time to t_new
    //
    const double eps = 1.e-3*dt_0;
    if (cur_time_d + dt_0 > (stop_time - start_time) - eps) {
        dt_0 = (stop_time - start_time) - cur_time_d;
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }

    diffusive_first_call = false;

    // Warn when the step just chosen exceeds the explicit diffusive limit.  The warning repeats
    // on a level only when the Fourier number has grown by half since it was last reported.
    if (diffusive_dt_check) {
        if (static_cast<int>(diffusive_fourier_warned.size()) < finest_level+1) {
            diffusive_fourier_warned.resize(finest_level+1, zero);
        }
        for (int lev = 0; lev <= finest_level; ++lev) {
            const Real F_mom  = static_cast<Real>(dt[lev]) * rate_mom[lev];
            const Real F_scal = static_cast<Real>(dt[lev]) * rate_scal[lev];
            const Real F      = amrex::max(F_mom, F_scal);
            if (verbose > 1) {
                Print() << "Diffusive Fourier number at level " << lev << ": " << F
                        << " (momentum " << F_mom << ", scalars " << F_scal << ")" << std::endl;
            }
            if (F > diffusive_cfl && F > Real(1.5) * diffusive_fourier_warned[lev]) {
                Print() << "WARNING: explicit eddy diffusion at level " << lev
                        << " has Fourier number dt * rate = " << F
                        << " (momentum " << F_mom << ", scalars " << F_scal
                        << ") > erf.diffusive_cfl = " << diffusive_cfl
                        << " with dt = " << dt[lev] << ".\n"
                        << "         The rate is a conservative estimate of the horizontal, terrain-metric and"
                        << " explicit vertical eddy diffusion; the run may go unstable.  Consider"
                        << " erf.diffusive_dt_limit = true (adaptive dt), a smaller erf.fixed_dt, or for"
                        << " Smagorinsky2D erf.smag2d_slope_limiter / erf.smag2d_kh_cap."
                        << "  Repeated only if it grows by half." << std::endl;
                diffusive_fourier_warned[lev] = F;
            }
        }
    }
}

/**
 * Largest explicit eddy-diffusion rates on a level, for the diffusive time-step check: the
 * maximum over valid cells of CellDiffusiveRates (ERF_TerrainDiffusionLimits.H), which takes
 * the largest diffusivities and the smallest density over each cell's 3x3x3 neighbourhood,
 * built from the eddy diffusivities of the last step (zero
 * before the first step and on the first step after a restart, since they are not
 * checkpointed; a regrid happens before the advance that recomputes them).  Molecular
 * diffusion is not included.
 *
 * @param[in]  lev       level
 * @param[out] rate_mom  momentum rate (1/s), Mom_h and Mom_v
 * @param[out] rate_scal scalar rate (1/s): the largest over theta, and moisture, turbulent
 *                       kinetic energy and the advected scalar where they are carried
 */
void
ERF::ComputeDiffusiveRates (int lev, Real& rate_mom, Real& rate_scal) const
{
    rate_mom  = zero;
    rate_scal = zero;
    if (!eddyDiffs_lev[lev]) { return; }

    const auto dxinv       = geom[lev].InvCellSizeArray();
    const bool variable_dz = (SolverChoice::mesh_type != MeshType::ConstantDz);

    // Explicit fraction of the vertical diffusion (ExplicitVerticalFraction): 0 where the
    // partly implicit stages keep it bounded, 1 otherwise.  Compressible levels use the
    // three-stage scheme; anelastic MidPoint is two stages; anelastic RK2 turns the implicit
    // solve off (its factors are zero).
    const auto& fac = solverChoice.vert_implicit_fac[lev];
    const int nstages = (solverChoice.anelastic[lev] &&
                         solverChoice.anelastic_type[lev] == AnelasticType::MidPoint) ? 2 : 3;
    const Real e_uv = ExplicitVerticalFraction(nstages, fac, solverChoice.implicit_momentum_diffusion);
#ifdef ERF_IMPLICIT_W
    const Real e_w  = e_uv;
#else
    const Real e_w  = one;  // w's vertical diffusion is always explicit in this build
#endif
    const Real e_th = ExplicitVerticalFraction(nstages, fac, solverChoice.implicit_thermal_diffusion);
    const Real e_ke = ExplicitVerticalFraction(nstages, fac, solverChoice.implicit_ke_diffusion);
    // The implicit solve treats only the first moisture variable (qv); any other moist species
    // shares the Q diffusivities and diffuses explicitly, and so do the advected scalars.
    // "Only qv" is read from the condensate indices: every scheme that carries number or bin
    // variables also carries qc.
    const auto& mi = solverChoice.moisture_indices;
    const bool only_qv = (mi.qc < 0 && mi.qi < 0 && mi.qr < 0 && mi.qs < 0 && mi.qg < 0);
    const Real e_q  = only_qv ? ExplicitVerticalFraction(nstages, fac, solverChoice.implicit_moisture_diffusion)
                              : one;
    DiffusiveRateSettings set;
    set.variable_dz = variable_dz;
    set.has_q       = (solverChoice.moisture_type != MoistureType::None);
    set.has_ke      = solverChoice.turbChoice[lev].use_tke;
    set.has_scalar  = solverChoice.transport_scalar;
    set.e_uv = e_uv; set.e_w = e_w; set.e_th = e_th; set.e_q = e_q; set.e_ke = e_ke;

    const MultiFab& K = *eddyDiffs_lev[lev];
    const MultiFab& S = vars_new[lev][Vars::cons];
    // CellDiffusiveRates reads the 3x3x3 neighbourhood of each valid cell
    AMREX_ALWAYS_ASSERT(K.nGrowVect().min() >= 1 && S.nGrowVect().min() >= 1);

    ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
    ReduceData<Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<const Real> s_arr = S.const_array(mfi);
        const Array4<const Real> k_arr = K.const_array(mfi);
        const Array4<const Real> z_nd  = z_phys_nd[lev]->const_array(mfi);
        const Array4<const Real> mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const Array4<const Real> mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);

        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real rm = zero, rs = zero;
            CellDiffusiveRates(i, j, k, s_arr, k_arr, z_nd, mf_mx(i,j,0), mf_my(i,j,0), dxinv, set, rm, rs);
            return {rm, rs};
        });
    }

    ReduceTuple hv = reduce_data.value(reduce_op);
    rate_mom  = amrex::get<0>(hv);
    rate_scal = amrex::get<1>(hv);
    ParallelDescriptor::ReduceRealMax(rate_mom);
    ParallelDescriptor::ReduceRealMax(rate_scal);
}

/**
 * Print the largest terrain slope factor alpha = h dx/dz (TerrainSlopeFactor) on a level, and
 * how many cells have alpha > 1.  ComputeDt calls it on terrain-fitted meshes where the closure
 * lets K_h and K_v differ (the explicit terrain-metric diffusion K_h h^2 d2/dz2 then limits the
 * time step), only when the diffusive check or limit or a Smagorinsky2D limit is on.
 *
 * @param[in] lev   level
 * @param[in] force print even if the largest alpha is within 1 % of the last report (new grids)
 */
void
ERF::ReportTerrainSlopeFactor (int lev, bool force)
{
    AMREX_ALWAYS_ASSERT(z_phys_nd[lev]);
    const TurbChoice& tc = solverChoice.turbChoice[lev];

    const MultiFab& S = vars_new[lev][Vars::cons];

    ReduceOps<ReduceOpMax, ReduceOpSum> reduce_op;
    ReduceData<Real, Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<const Real> z_nd = z_phys_nd[lev]->const_array(mfi);
        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            const Real alpha = TerrainSlopeFactor(ComputeTerrainCellDrops(i,j,k,z_nd));
            return {alpha, (alpha > one) ? Long(1) : Long(0)};
        });
    }

    ReduceTuple hv = reduce_data.value(reduce_op);
    Real alpha_max = amrex::get<0>(hv);
    Long n_steep   = amrex::get<1>(hv);
    ParallelDescriptor::ReduceRealMax(alpha_max);

    if (!SlopeReportPrint(force, alpha_max, slope_report_alpha[lev])) { return; }
    ParallelDescriptor::ReduceLongSum(n_steep);   // needed only for the printed line
    slope_report_alpha[lev] = alpha_max;

    Print() << "Terrain slope factor alpha = h dx/dz at level " << lev << ": max " << alpha_max
            << ", alpha > 1 in " << n_steep << " of " << grids[lev].numPts() << " cells" << std::endl;
    if (alpha_max > one) {
        Print() << "    the explicit terrain-metric diffusion K h^2 d2/dz2 limits dt to about"
                << " dz^2 / (4 K_m h^2) for momentum and dz^2 / (2 K_s h^2) for scalars";
        if (tc.les_type == LESType::Smagorinsky && tc.smag2d) {
            Print() << ",\n    i.e. 1 / (4 Cs^2 |S| alpha^2) and Pr_t / (2 Cs^2 |S| alpha^2) for Smagorinsky2D"
                    << " without the WRF slope limiter (";
            if (tc.smag2d_slope_limiter) {
                Print() << "on here: K_h is divided by alpha or alpha^2)";
            } else {
                Print() << "off here)";
            }
        }
        Print() << std::endl;
    }
}

/**
 * Function that calls estTimeStep for each level
 *
 * @param[in] level level of refinement (coarsest level i 0)
 * @param[out] dt_fast_ratio ratio of slow to fast time step
 */
double
ERF::estTimeStep (int level, long& dt_fast_ratio) const
{
    BL_PROFILE("ERF::estTimeStep()");

    // Terrain aware (T) and terrain unaware (N) time step estimates.
    double estdt_comp_T = bogus_large_value;
    double estdt_comp_N = bogus_large_value;
    double estdt_lowM_T = bogus_large_value;
    double estdt_lowM_N = bogus_large_value;

    // We intentionally use the level 0 domain to compute whether to use this direction in the dt calculation
    const int nxc = geom[0].Domain().length(0);
    const int nyc = geom[0].Domain().length(1);

    auto const dxinv = geom[level].InvCellSizeArray();
    auto dxinv_EB = dxinv; dxinv_EB[2] = one / dz_min[level];

    MultiFab const& S_new = vars_new[level][Vars::cons];

    // Keep the thermodynamic samples alongside the cell-centered velocity so
    // the wall-rate reduction can call the same pointwise MOST evaluator as
    // production wall transfer without allocating a second global temporary.
    MultiFab ccvel_N(grids[level],dmap[level],7,0);
    MultiFab ccvel_T(grids[level],dmap[level],3,0);

    int klo = geom[level].Domain().smallEnd(2);
    int khi = geom[level].Domain().bigEnd(2);
    MultiFab omega(convert(grids[level],IntVect(0,0,1)),dmap[level],1,0);
    for (MFIter mfi(vars_new[level][Vars::zvel]); mfi.isValid(); ++mfi)
    {
        Box vbx = mfi.validbox();

        const Array4<      Real>& omega_arr = omega.array(mfi);
        const Array4<const Real>& u_arr     = vars_new[level][IntVars::xmom].const_array(mfi);
        const Array4<const Real>& v_arr     = vars_new[level][IntVars::ymom].const_array(mfi);
        const Array4<const Real>& w_arr     = vars_new[level][IntVars::zmom].const_array(mfi);

        const Array4<const Real>& z_nd_arr  = z_phys_nd[level]->const_array(mfi);

        const Array4<const Real>& mf_ux     = mapfac[level][MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy     = mapfac[level][MapFacType::v_y]->const_array(mfi);

        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k==klo || k==(khi+1)) {
                omega_arr(i,j,k) = zero;
            } else {
                omega_arr(i,j,k) = OmegaFromW(i,j,k,w_arr(i,j,k),
                                              u_arr,v_arr,mf_ux,mf_vy,
                                              z_nd_arr,dxinv);
            }
        });
    }

    average_face_to_cellcenter(ccvel_T,0,
                               Array<const MultiFab*,3>{&vars_new[level][Vars::xvel],
                                                        &vars_new[level][Vars::yvel],
                                                        &omega});
    average_face_to_cellcenter(ccvel_N,0,
                               Array<const MultiFab*,3>{&vars_new[level][Vars::xvel],
                                                        &vars_new[level][Vars::yvel],
                                                        &vars_new[level][Vars::zvel],});

    const bool chamber_cloudy = cloud_chamber_config.active &&
        cloud_chamber_config.cloudy;
    const Real rdOcp = solverChoice.rdOcp;
    const MultiFab& chamber_base_state = base_state[level];
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const Array4<const Real> state = S_new.const_array(mfi);
        const Array4<const Real> base = chamber_base_state.const_array(mfi);
        const Array4<Real> velocity = ccvel_N.array(mfi);
        const Box bx = mfi.validbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real rho = state(i,j,k,Rho_comp);
            velocity(i,j,k,3) = state(i,j,k,RhoTheta_comp) / rho;
            velocity(i,j,k,4) = chamber_cloudy ?
                state(i,j,k,RhoQ1_comp) / rho : Real(0.0);
            velocity(i,j,k,5) = base(i,j,k,BaseState::p0_comp);
            velocity(i,j,k,6) = rdOcp;
        });
    }

    bool l_substepping = (solverChoice.substepping_type[level] == SubsteppingType::Implicit);
    int  l_anelastic   = solverChoice.anelastic[level];

    bool l_comp_substepping_diag = (verbose && l_substepping && !l_anelastic && solverChoice.substepping_diag);

    Real estdt_comp_inv_N, estdt_comp_inv_T;
    Real estdt_lowM_inv_N, estdt_lowM_inv_T;
    Real estdt_vert_comp_inv, estdt_vert_lowM_inv;

    const MultiFab& z_nd_mf = *z_phys_nd[level];

    if (l_substepping && (nxc==1) && (nyc==1)) {
        // SCM -- should not depend on dx or dy; force minimum number of substeps
        estdt_comp_inv_T = std::numeric_limits<Real>::min();
        estdt_comp_inv_N = estdt_comp_inv_T;
    }
    else if (solverChoice.terrain_type == TerrainType::EB)
    {
        const eb_& eb_lev = get_eb(level);
        const MultiFab& detJ = (eb_lev.get_const_factory())->getVolFrac();

        estdt_comp_inv_N = ReduceMax(S_new, ccvel_N, detJ, 0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                   Array4<Real const> const& s,
                                   Array4<Real const> const& u,
                                   Array4<Real const> const& vf) -> Real
        {
           Real new_comp_dt = -bogus_large_value;
           amrex::Loop(b, [=,&new_comp_dt] (int i, int j, int k) noexcept
           {
               if (vf(i,j,k) > zero)
               {
                   const Real rho      = s(i, j, k, Rho_comp);
                   const Real rhotheta = s(i, j, k, RhoTheta_comp);

                   // NOTE: even when moisture is present,
                   //       we only use the partial pressure of the dry air
                   //       to compute the soundspeed
                   Real pressure = getPgivenRTh(rhotheta);
                   Real c = std::sqrt(Gamma * pressure / rho);

                   // If we are doing implicit acoustic substepping, then the z-direction does not contribute
                   //    to the computation of the time step
                   if (l_substepping) {
                       if ((nxc > 1) && (nyc==1)) {
                           // 2-D in x-z
                           new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,0))+c)*dxinv_EB[0]), new_comp_dt);
                       } else if ((nyc > 1) && (nxc==1)) {
                           // 2-D in y-z
                           new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,1))+c)*dxinv_EB[1]), new_comp_dt);
                       } else {
                           // 3-D
                           new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,0))+c)*dxinv_EB[0]),
                                                    ((amrex::Math::abs(u(i,j,k,1))+c)*dxinv_EB[1]), new_comp_dt);
                       }

                   // If we are not doing implicit acoustic substepping, then the z-direction contributes
                   //    to the computation of the time step
                   } else {
                       if (nxc > 1 && nyc > 1) {
                           new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,0))+c)*dxinv_EB[0]),
                                                    ((amrex::Math::abs(u(i,j,k,1))+c)*dxinv_EB[1]),
                                                    ((amrex::Math::abs(u(i,j,k,2))+c)*dxinv_EB[2]), new_comp_dt);
                       } else if (nxc > 1) {
                           new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,0))+c)*dxinv_EB[0]),
                                                    ((amrex::Math::abs(u(i,j,k,2))+c)*dxinv_EB[2]), new_comp_dt);
                       } else if (nyc > 1) {
                           new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,1))+c)*dxinv_EB[1]),
                                                    ((amrex::Math::abs(u(i,j,k,2))+c)*dxinv_EB[2]), new_comp_dt);
                       } else {
                           new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,2))+c)*dxinv_EB[2]), new_comp_dt);
                       }

                   }
               }
           });
           return new_comp_dt;
       });

        // The metric terms do not exist for EB, so the terrain aware
        // estimate is identical to the terrain unaware estimate
        estdt_comp_inv_T = estdt_comp_inv_N;

    } else {
        // One pass over the data returning both the terrain aware (T) estimate,
        // and the terrain unaware (N) estimate
        ReduceOps<ReduceOpMax,ReduceOpMax> reduce_op;
        ReduceData<Real,Real>              reduce_data(reduce_op);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            const Array4<const Real>& s    = S_new.const_array(mfi);
            const Array4<const Real>& u_T  = ccvel_T.const_array(mfi);
            const Array4<const Real>& u_N  = ccvel_N.const_array(mfi);
            const Array4<const Real>& z_nd = z_nd_mf.const_array(mfi);

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real,Real>
            {
                const Real rho      = s(i,j,k,Rho_comp);
                const Real rhotheta = s(i,j,k,RhoTheta_comp);

                // NOTE: even when moisture is present,
                //       we only use the partial pressure of the dry air
                //       to compute the soundspeed
                Real pressure = getPgivenRTh(rhotheta);
                Real c = std::sqrt(Gamma * pressure / rho);

                Real inv_dt_T = Compute_InvDt_Compressible(i,j,k, c,
                                    u_T(i,j,k,0), u_T(i,j,k,1), u_T(i,j,k,2),
                                    z_nd, dxinv, true, l_substepping, nxc, nyc);

                Real inv_dt_N = Compute_InvDt_Compressible(i,j,k, c,
                                    u_N(i,j,k,0), u_N(i,j,k,1), u_N(i,j,k,2),
                                    z_nd, dxinv, false, l_substepping, nxc, nyc);

                return {inv_dt_T, inv_dt_N};
            });
        }

        GpuTuple<Real,Real> hv = reduce_data.value(reduce_op);
        estdt_comp_inv_T       = amrex::get<0>(hv);
        estdt_comp_inv_N       = amrex::get<1>(hv);
    } // not EB

    {
        Real comp_inv[2] = {estdt_comp_inv_T, estdt_comp_inv_N};
        ParallelDescriptor::ReduceRealMax(comp_inv,2);
        estdt_comp_inv_T = comp_inv[0];
        estdt_comp_inv_N = comp_inv[1];
    }

    // Globally empty level -> ReduceMax = lowest(); treat level as non-constraining.
    estdt_comp_T = (estdt_comp_inv_T > zero) ? (cfl / estdt_comp_inv_T) : bogus_large_value;
    estdt_comp_N = (estdt_comp_inv_N > zero) ? (cfl / estdt_comp_inv_N) : bogus_large_value;

    //
    // Anelastic (low Mach) estimate -- purely advective, again terrain aware
    // and terrain unaware
    //
    if (solverChoice.terrain_type == TerrainType::EB)
    {
        estdt_lowM_inv_N = ReduceMax(ccvel_N, z_nd_mf, 0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                   Array4<Real const> const& u,
                                   Array4<Real const> const& z_nd) -> Real
        {
            Real new_lm_dt = -bogus_large_value;
            Loop(b, [=,&new_lm_dt] (int i, int j, int k) noexcept
            {
                Real inv_dt_lowM_N = Compute_InvDt_Anelastic(i,j,k,
                                             u(i,j,k,0), u(i,j,k,1), u(i,j,k,2),
                                             z_nd, dxinv_EB, false);
                new_lm_dt = amrex::max(inv_dt_lowM_N, new_lm_dt);
            });
            return new_lm_dt;
        });

        // The metric terms do not exist for EB
        estdt_lowM_inv_T = estdt_lowM_inv_N;

    } else {

        ReduceOps<ReduceOpMax,ReduceOpMax> reduce_op;
        ReduceData<Real,Real>              reduce_data(reduce_op);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ccvel_T,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            const Array4<const Real>& u_T  = ccvel_T.const_array(mfi);
            const Array4<const Real>& u_N  = ccvel_N.const_array(mfi);
            const Array4<const Real>& z_nd = z_nd_mf.const_array(mfi);

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real,Real>
            {
                Real inv_dt_T = Compute_InvDt_Anelastic(i,j,k,
                                    u_T(i,j,k,0), u_T(i,j,k,1), u_T(i,j,k,2),
                                    z_nd, dxinv, true);

                Real inv_dt_N = Compute_InvDt_Anelastic(i,j,k,
                                    u_N(i,j,k,0), u_N(i,j,k,1), u_N(i,j,k,2),
                                    z_nd, dxinv, false);

                return {inv_dt_T, inv_dt_N};
            });
        }

        GpuTuple<Real,Real> hv = reduce_data.value(reduce_op);
        estdt_lowM_inv_T       = amrex::get<0>(hv);
        estdt_lowM_inv_N       = amrex::get<1>(hv);
    }

    {
        Real lowM_inv[2] = {estdt_lowM_inv_T, estdt_lowM_inv_N};
        ParallelDescriptor::ReduceRealMax(lowM_inv,2);
        estdt_lowM_inv_T = lowM_inv[0];
        estdt_lowM_inv_N = lowM_inv[1];
    }

    if (estdt_lowM_inv_T > zero) { estdt_lowM_T = cfl / estdt_lowM_inv_T; }
    if (estdt_lowM_inv_N > zero) { estdt_lowM_N = cfl / estdt_lowM_inv_N; }

    // The host estimator uses a directional maximum, whereas M3's donor
    // positivity condition is based on the sum of all mapped outgoing faces.
    // Rebuild the current-time dry-air carrier with the same momentum and
    // static-terrain conventions as AdvectionSrcForRho.
    if (solverChoice.moisture_type == MoistureType::SBM) {
        const BoxArray& cell_ba = S_new.boxArray();
        const DistributionMapping& cell_dm = S_new.DistributionMap();
        const BoxArray xface_ba = convert(
            cell_ba, IntVect::TheDimensionVector(0));
        const BoxArray yface_ba = convert(
            cell_ba, IntVect::TheDimensionVector(1));
        const BoxArray zface_ba = convert(
            cell_ba, IntVect::TheDimensionVector(2));

        MultiFab rho_u(xface_ba, cell_dm, 1, 2);
        MultiFab rho_v(yface_ba, cell_dm, 1, 2);
        MultiFab rho_w(zface_ba, cell_dm, 1, 2);
        rho_u.setVal(Real(0.0));
        rho_v.setVal(Real(0.0));
        rho_w.setVal(Real(0.0));

        const IntVect valid_faces(0);
        const auto& current_u = vars_new[level][Vars::xvel];
        const auto& current_v = vars_new[level][Vars::yvel];
        const auto& current_w = vars_new[level][Vars::zvel];
        VelocityToMomentum(current_u, valid_faces, current_v, valid_faces,
                           current_w, valid_faces, S_new, rho_u, rho_v, rho_w,
                           geom[level].Domain(), domain_bcs_type, nullptr);
        rho_u.FillBoundary(geom[level].periodicity());
        rho_v.FillBoundary(geom[level].periodicity());
        rho_w.FillBoundary(geom[level].periodicity());

        erf_auxiliary::MappedFaceFluxRate mapped_carrier;
        mapped_carrier.define(cell_ba, cell_dm, 1, 0);
        std::string sbm_diagnostic;
        if (l_anelastic) {
            // Projection restores rho0*w before the anelastic stage seam, so
            // the donor estimate must use the same three momentum carriers.
            // Terrain Omega is only an internal projection representation.
            if (!erf_auxiliary::CopyNativeMappedDryAirCarrierFluxRate(
                    mapped_carrier, rho_u, rho_v, rho_w, sbm_diagnostic)) {
                amrex::Abort("SBM anelastic native carrier: " + sbm_diagnostic);
            }
        } else {
            MultiFab sbm_vertical_carrier(zface_ba, cell_dm, 1, 0);
            const bool terrain_fitted =
                solverChoice.mesh_type == MeshType::VariableDz;
            const int zlo = geom[level].Domain().smallEnd(2);
            const int zhi = geom[level].Domain().bigEnd(2);
            const auto inv_dx = geom[level].InvCellSizeArray();
            const MultiFab& z_nd = *z_phys_nd[level];
            const MultiFab& mf_ux = *mapfac[level][MapFacType::u_x];
            const MultiFab& mf_vy = *mapfac[level][MapFacType::v_y];
            for (MFIter mfi(sbm_vertical_carrier, TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const Box bx = mfi.tilebox();
                const auto ru = rho_u.const_array(mfi);
                const auto rv = rho_v.const_array(mfi);
                const auto rw = rho_w.const_array(mfi);
                const auto ux = mf_ux.const_array(mfi);
                const auto vy = mf_vy.const_array(mfi);
                const auto z = z_nd.const_array(mfi);
                const auto out = sbm_vertical_carrier.array(mfi);
                const bool fitted = terrain_fitted;
                const int bottom = zlo;
                const int top = zhi + 1;
                const auto dx = inv_dx;
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    if (!fitted) {
                        out(i, j, k, 0) = rw(i, j, k, 0);
                    } else if (k == bottom) {
                        out(i, j, k, 0) = Real(0.0);
                    } else if (k == top) {
                        out(i, j, k, 0) = rw(i, j, k, 0);
                    } else {
                        out(i, j, k, 0) = OmegaFromW(i, j, k, rw(i, j, k, 0),
                                                      ru, rv, ux, vy, z, dx);
                    }
                });
            }
            const bool carrier_ok =
                erf_auxiliary::BuildMappedDryAirCarrierFluxRate(
                    mapped_carrier, rho_u, rho_v, sbm_vertical_carrier,
                    *ax[level], *ay[level], *az[level],
                    *mapfac[level][MapFacType::u_y],
                    *mapfac[level][MapFacType::v_x],
                    *mapfac[level][MapFacType::m_x],
                    *mapfac[level][MapFacType::m_y], sbm_diagnostic);
            if (!carrier_ok) {
                Abort("SBM M3 current mapped carrier: " + sbm_diagnostic);
            }
        }

        if (sbm_transport == nullptr || !sbm_transport->is_defined(level) ||
            !sbm_transport->measure_is_ready(level)) {
            Abort("SBM M3 donor timestep estimate requires a ready static measure");
        }
        const MultiFab& measure = sbm_transport->static_measure(level);

        Real max_sbm_outgoing_rate = Real(0.0);
        const bool rate_ok = erf_auxiliary::ComputeMaxMappedOutgoingRate(
            mapped_carrier, measure, S_new, Rho_comp, dxinv,
            max_sbm_outgoing_rate, sbm_diagnostic);
        if (!rate_ok) {
            Abort("SBM M3 current donor-rate estimate: " + sbm_diagnostic);
        }
        if (max_sbm_outgoing_rate > Real(0.0)) {
            const double adaptive_sbm_dt =
                static_cast<double>(amrex::min(cfl, Real(1.0)) /
                                    max_sbm_outgoing_rate);
            estdt_comp_T = std::min(estdt_comp_T, adaptive_sbm_dt);
            estdt_comp_N = std::min(estdt_comp_N, adaptive_sbm_dt);
            estdt_lowM_T = std::min(estdt_lowM_T, adaptive_sbm_dt);
            estdt_lowM_N = std::min(estdt_lowM_N, adaptive_sbm_dt);

            double hard_sbm_dt = std::numeric_limits<double>::infinity();
            if (erf_auxiliary::FixedDtExceedsMappedDonorLimit(
                    static_cast<double>(fixed_dt[level]),
                    max_sbm_outgoing_rate, hard_sbm_dt)) {
                std::ostringstream message;
                message.precision(17);
                message << "SBM M3 fixed timestep exceeds donor positivity "
                           "limit: fixed_dt=" << fixed_dt[level]
                        << " hard_limit=" << hard_sbm_dt
                        << " max_outgoing_rate=" << max_sbm_outgoing_rate;
                Abort(message.str());
            }
            if (verbose) {
                Print() << "SBM mapped donor dt at level " << level << ": "
                        << adaptive_sbm_dt << " (hard fixed-dt limit "
                        << hard_sbm_dt << ")" << std::endl;
            }
        }
    }

     Real max_wall_rate = Real(0.0);
     Real estdt_wall = bogus_large_value;
     if (cloud_chamber_config.active &&
         cloud_chamber_config.physical_initialization &&
         cloud_chamber_config.has_wall_rate_channel()) {
         const auto walls = cloud_chamber_config.wall_boundary();
         const Box domain = geom[level].Domain();
         max_wall_rate = ReduceMax(ccvel_N, 0,
         [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                    Array4<Real const> const& velocity) -> Real
         {
             Real rate = Real(0.0);
             amrex::Loop(b, [=,&rate] (int i, int j, int k) noexcept
             {
                 // A free staggered component can receive tangential traction
                 // from every active perpendicular wall at an edge/corner.
                 // Fixed and neutral coefficients use the exact local
                 // velocity-parallel Jacobian row-sum factor.  MOST evaluates
                 // its current coefficients once per wall/cell and uses that
                 // state as a frozen-coefficient local rate estimate; it is
                 // not a nonlinear MOST Jacobian bound.
                 Real momentum_rate = Real(0.0);
                 amrex::GpuArray<Real, AMREX_SPACEDIM> low_momentum_rates{};
                 amrex::GpuArray<Real, AMREX_SPACEDIM> high_momentum_rates{};
                 // Evaluate each encountered wall/cell state once.  The
                 // component row sums below only compose these retained
                 // per-face rates; they must not repeat the MOST solve.
                 for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                     const bool low = (dir == 0 ? i == domain.smallEnd(0) :
                                       (dir == 1 ? j == domain.smallEnd(1) :
                                                    k == domain.smallEnd(2)));
                     const bool high = (dir == 0 ? i == domain.bigEnd(0) :
                                        (dir == 1 ? j == domain.bigEnd(1) :
                                                     k == domain.bigEnd(2)));
                     if (low) {
                         const auto& wall = walls[2*dir];
                         if (erf_cloud_chamber_wall_flux::
                             wall_rate_requires_tangential_speed(wall)) {
                             const Real U_t =
                                 erf_cloud_chamber_wall_flux::
                                 tangential_speed_cell_centered(
                                     dir, i, j, k, velocity, wall);
                             const auto runtime =
                                 erf_cloud_chamber_wall_flux::most_wall_coefficients(
                                     wall, velocity(i,j,k,3), velocity(i,j,k,4),
                                     velocity(i,j,k,5), velocity(i,j,k,6), U_t,
                                     Real(0.5) / dxinv[dir],
                                     dir == 2 ? 1 : 0);
                             rate = amrex::max(rate,
                                 erf_cloud_chamber_wall_flux::wall_rate_for_face(
                                     wall, U_t, dxinv[dir], runtime));
                             low_momentum_rates[dir] =
                                 erf_cloud_chamber_wall_flux::momentum_rate_for_face(
                                     wall, U_t, dxinv[dir], runtime);
                         }
                     }
                     if (high) {
                         const auto& wall = walls[2*dir+1];
                         if (erf_cloud_chamber_wall_flux::
                             wall_rate_requires_tangential_speed(wall)) {
                             const Real U_t =
                                 erf_cloud_chamber_wall_flux::
                                 tangential_speed_cell_centered(
                                     dir, i, j, k, velocity, wall);
                             const auto runtime =
                                 erf_cloud_chamber_wall_flux::most_wall_coefficients(
                                     wall, velocity(i,j,k,3), velocity(i,j,k,4),
                                     velocity(i,j,k,5), velocity(i,j,k,6), U_t,
                                     Real(0.5) / dxinv[dir],
                                     dir == 2 ? -1 : 0);
                             rate = amrex::max(rate,
                                 erf_cloud_chamber_wall_flux::wall_rate_for_face(
                                     wall, U_t, dxinv[dir], runtime));
                             high_momentum_rates[dir] =
                                 erf_cloud_chamber_wall_flux::momentum_rate_for_face(
                                     wall, U_t, dxinv[dir], runtime);
                         }
                     }
                 }
                 for (int component = 0; component < AMREX_SPACEDIM; ++component) {
                     const Real component_rate =
                         erf_cloud_chamber_wall_flux::momentum_row_sum_rate(
                             component, low_momentum_rates, high_momentum_rates);
                     momentum_rate = amrex::max(momentum_rate, component_rate);
                 }
                 rate = amrex::max(rate, momentum_rate);
             });
             return rate;
         });
         ParallelDescriptor::ReduceRealMax(max_wall_rate);
         if (max_wall_rate > Real(0.0)) {
             estdt_wall = erf_cloud_chamber_wall_flux::wall_dt_from_max_rate(max_wall_rate);
         }
     }
     // Wall-rate kernels and reductions use amrex::Real, while ERF's host
     // timestep estimates are double even in ERF_PRECISION=SINGLE builds.
     const double estdt_wall_host = static_cast<double>(estdt_wall);
     estdt_comp_T = std::min(estdt_comp_T, estdt_wall_host);
     estdt_comp_N = std::min(estdt_comp_N, estdt_wall_host);
     estdt_lowM_T = std::min(estdt_lowM_T, estdt_wall_host);
     estdt_lowM_N = std::min(estdt_lowM_N, estdt_wall_host);

     const double fixed_dt_level = static_cast<double>(fixed_dt[level]);
     erf_cloud_chamber_wall_dt_guard::enforce_fixed_dt_limit(
         level, fixed_dt_level, estdt_wall_host, max_wall_rate);

     // Additional vertical diagnostics
     if (l_comp_substepping_diag) {
         estdt_vert_comp_inv = ReduceMax(S_new, ccvel_T, z_nd_mf, 0,
         [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                    Array4<Real const> const& s,
                                    Array4<Real const> const& u,
                                    Array4<Real const> const& z_nd) -> Real
         {
             Real new_comp_dt = -bogus_large_value;
             amrex::Loop(b, [=,&new_comp_dt] (int i, int j, int k) noexcept
             {
                 {
                     const Real rho      = s(i, j, k, Rho_comp);
                     const Real rhotheta = s(i, j, k, RhoTheta_comp);

                     // NOTE: even when moisture is present,
                     //       we only use the partial pressure of the dry air
                     //       to compute the soundspeed
                     Real pressure = getPgivenRTh(rhotheta);
                     Real c = std::sqrt(Gamma * pressure / rho);

                     Real h_zeta  = Compute_h_zeta_AtCellCenter(i,j,k,dxinv,z_nd);
                     Real idz_loc = dxinv[2] / h_zeta;

                     // Look at z-direction only
                     new_comp_dt = amrex::max((amrex::Math::abs(u(i,j,k,2)) + c) * idz_loc, new_comp_dt);
                 }
             });
             return new_comp_dt;
         });

         estdt_vert_lowM_inv = ReduceMax(ccvel_T, z_nd_mf, 0,
         [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                    Array4<Real const> const& u,
                                    Array4<Real const> const& z_nd) -> Real
         {
             Real new_lowM_dt = -bogus_large_value;
             amrex::Loop(b, [=,&new_lowM_dt] (int i, int j, int k) noexcept
             {
                 Real h_zeta  = Compute_h_zeta_AtCellCenter(i,j,k,dxinv,z_nd);
                 Real idz_loc = dxinv[2] / h_zeta;
                 new_lowM_dt = amrex::max((amrex::Math::abs(u(i,j,k,2))) * idz_loc, new_lowM_dt);
             });
             return new_lowM_dt;
         });

         ParallelDescriptor::ReduceRealMax(estdt_vert_comp_inv);
         ParallelDescriptor::ReduceRealMax(estdt_vert_lowM_inv);
     }

     if (verbose) {
         // Terrain aware   (T): includes h_xi, h_eta and h_zeta
         // Terrain unaware (N): no metric terms, dxinv[2] used for the vertical spacing
         if (fixed_dt[level] <= zero) {
             Print() << "Using cfl = " << cfl << " and dx/dy/dz_min = " <<
               one/dxinv[0] << " " << one/dxinv[1] << " " << dz_min[level] << std::endl;
             Print() << "Compressible dt at level " << level << ":  "
                     << estdt_comp_T << " (terrain aware)  "
                     << estdt_comp_N << " (terrain unaware)" << std::endl;
             if (estdt_lowM_inv_T > 0.0_rt) {
                 Print() << "Anelastic   dt at level " << level << ":  "
                         << estdt_lowM_T << " (terrain aware)  "
                         << estdt_lowM_N << " (terrain unaware)" << std::endl;
             } else {
                 Print() << "Anelastic dt at level " << level << ": undefined " << std::endl;
             }
         }

         if (fixed_dt[level] > zero) {
             Print() << "Based on cfl of one " << std::endl;
             Print() << "Compressible dt at level " << level << " would be:  "
                     << estdt_comp_T/cfl << " (terrain aware)  "
                     << estdt_comp_N/cfl << " (terrain unaware)" << std::endl;
             if (estdt_lowM_inv_T > zero) {
                 Print() << "Anelastic    dt at level " << level << " would be:  "
                         << estdt_lowM_T/cfl << " (terrain aware)  "
                         << estdt_lowM_N/cfl << " (terrain unaware)" << std::endl;
             } else {
                 Print() << "Anelastic    dt at level " << level << " would be undefined " << std::endl;
             }
             Print() << "Fixed dt at level " << level << "       is:  " << fixed_dt[level] << std::endl;
             if (fixed_fast_dt[level] > zero) {
                 Print() << "Fixed fast dt at level " << level << "       is:  " << fixed_fast_dt[level] << std::endl;
             }
         }
     }

     if (solverChoice.substepping_type[level] != SubsteppingType::None) {
         if (fixed_dt[level] > zero && fixed_fast_dt[level] > zero) {
             dt_fast_ratio = static_cast<long>( fixed_dt[level] / fixed_fast_dt[level] );
             if (dt_fast_ratio < 1) {
                 Abort("Invalid fixed_fast_dt: must be <= fixed_dt so mri_dt_ratio >= 1");
             }
         } else if (fixed_dt[level] > zero) {
             // Max CFL_c = one for substeps by default, but we enforce a min of 4 substeps
             const double dt_sub_max = static_cast<double>(estdt_comp_T * (third/cfl) * sub_cfl);
             //
             // Check dt_sub_max BEFORE dividing by it, and the quotient before casting it.
             // Converting a non-finite or wildly out-of-range double to long is undefined
             // behaviour and in practice traps, and the division itself raises on a zero
             // denominator under erf's FPE trapping.  A level handed corrupted state
             // arrives here with an enormous sound speed, hence a vanishing estdt_comp_T
             // and a vanishing dt_sub_max; without this the first symptom is a bare signal
             // inside the time-step estimate rather than a message naming the cause.
             //
             auto bad_substep_ratio = [&] (const char* what, double value)
             {
                 std::ostringstream message;
                 message.precision(17);
                 message << "estTimeStep: cannot form the acoustic substep ratio at level "
                         << level << ": " << what << " = " << value
                         << " (fixed_dt=" << fixed_dt[level]
                         << ", estdt_comp_T=" << estdt_comp_T
                         << ", cfl=" << cfl << ", substepping_cfl=" << sub_cfl << ")."
                         << " The compressible time-step estimate came back zero or"
                            " non-finite, so the state this level was handed is not usable.";
                 Abort(message.str());
             };

             //
             // +inf is deliberately allowed through: an SCM level and a globally empty
             // level both report an unconstrained acoustic dt on purpose, the quotient
             // below is then zero, and the ratio correctly floors at the minimum of 4.
             // Only a zero, negative or NaN denominator is a real failure.
             //
             if (std::isnan(dt_sub_max) || (dt_sub_max <= 0.0)) {
                 bad_substep_ratio("dt_sub_max", dt_sub_max);
             }
             const double sub_ratio = std::max(fixed_dt[level]/dt_sub_max, 4.0);
             if (!std::isfinite(sub_ratio) || (sub_ratio > 1.0e9)) {
                 bad_substep_ratio("mri_dt_ratio", sub_ratio);
             }
             dt_fast_ratio = static_cast<long>( sub_ratio );
         } else {
             // auto dt_sub_max = (estdt_comp_T/cfl * sub_cfl);
             // dt_fast_ratio = static_cast<long>( std::max(estdt_comp_T/dt_sub_max,Real(4.)) );
         dt_fast_ratio = static_cast<long>( std::max((cfl/third) / sub_cfl, Real(4.)) );
         }

         // Force time step ratio to be an even value
         if (solverChoice.force_stage1_single_substep) {
             if ( dt_fast_ratio%2 != 0) dt_fast_ratio += 1;
         } else {
             if ( dt_fast_ratio%6 != 0) {
                 Print() << "mri_dt_ratio = " << dt_fast_ratio
                         << " not divisible by 6 for N/3 substeps in stage 1" << std::endl;
                 dt_fast_ratio = static_cast<int>(std::ceil(dt_fast_ratio/Real(6.0)) * 6);
             }
         }

         if (verbose) {
             Print() << "smallest even ratio is: " << dt_fast_ratio << std::endl;
         }
     } // if substepping

     // Print out some extra diagnostics -- dt calcs are repeated so as to not
     // disrupt the overall code flow...
     if (l_comp_substepping_diag) {
         double dt_diag = (fixed_dt[level] > zero) ? fixed_dt[level] : static_cast<double>(estdt_comp_T);
         int  ns      = (fixed_mri_dt_ratio > zero) ? fixed_mri_dt_ratio : dt_fast_ratio;

         // horizontal acoustic CFL must be < 1 (fully explicit)
         // vertical   acoustic CFL may  be > 1
         Print() << "effective horiz,vert acoustic CFL with " << ns << " substeps : "
            << (dt_diag / ns) * estdt_comp_inv_T << " "
            << (dt_diag / ns) * estdt_vert_comp_inv << std::endl;

         // vertical advective CFL should be < 1, otherwise w-damping may be needed
         Print() << "effective vert advective CFL : "
            << dt_diag * estdt_vert_lowM_inv << std::endl;
     }

     if (fixed_dt[level] > zero) {
         return fixed_dt[level];
     } else {
         // Anelastic (substepping is not allowed)
         if (l_anelastic) {

            // Make sure that timestep is less than the dt_max
            estdt_lowM_T = std::min(estdt_lowM_T, dt_max);

            // On the first timestep enforce dt_max_initial
            if (istep[level] == 0) {
                return std::min(dt_max_initial, estdt_lowM_T);
            } else {
                return estdt_lowM_T;
            }


         // Compressible with or without substepping
         } else {
             return estdt_comp_T;
         }
     }
}
