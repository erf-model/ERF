/**
 * \file ERF_WriteScalarProfiles.cpp
 */
#include <iomanip>
#include <limits>
#include "ERF_Constants.H"

#include "ERF.H"
#include "ERF_Derive.H"
#include "ERF_MoistUtils.H"
#include "ERF_MicrophysicsConstants.H"

using namespace amrex;

/**
 * Computes the integrated quantities on the grid such as the
 * total scalar and total mass quantities. Prints and writes to output file.
 *
 * @param time Current time
 */
void
ERF::sum_integrated_quantities (double time)
{
    BL_PROFILE("ERF::sum_integrated_quantities()");

    if (verbose <= 0)
      return;

    // Single level sum
    Real mass_sl;

    // Multilevel sums
    Real mass_ml = zero;
    Real rhth_ml = zero;
    Real scal_ml = zero;
    Real mois_ml = zero;
#ifdef ERF_USE_FIRE
    Real smoke_sl = zero;
    Real smoke_ml = zero;
#endif

    bool local = true;

    auto& mfx0 = *mapfac[0][MapFacType::m_x];
    auto& mfy0 = *mapfac[0][MapFacType::m_y];
    auto&  dJ0 = *detJ_cc[0];

    mass_sl = volWgtSumMF(0,vars_new[0][Vars::cons],Rho_comp,dJ0,mfx0,mfy0,false,local);

    for (int lev = 0; lev <= finest_level; lev++) {
        auto& mfx = *mapfac[lev][MapFacType::m_x];
        auto& mfy = *mapfac[lev][MapFacType::m_y];
        auto&  dJ = *detJ_cc[lev];
        mass_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons],Rho_comp,dJ,mfx,mfy,true);
    }

    Real rhth_sl = volWgtSumMF(0,vars_new[0][Vars::cons], RhoTheta_comp,dJ0,mfx0,mfy0,false);
    Real scal_sl = volWgtSumMF(0,vars_new[0][Vars::cons],RhoScalar_comp,dJ0,mfx0,mfy0,false);
    Real mois_sl = zero;
    if (solverChoice.moisture_type != MoistureType::None) {
        int n_qstate_into_total = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
        for (int qoff(0); qoff<n_qstate_into_total; ++qoff) {
            mois_sl += volWgtSumMF(0,vars_new[0][Vars::cons],RhoQ1_comp+qoff,dJ0,mfx0,mfy0,false);
        }
    }
#ifdef ERF_USE_FIRE
    if (m_fire_coupler) {
        smoke_sl = volWgtSumMF(0, vars_new[0][Vars::cons], RhoSmoke_comp, dJ0, mfx0, mfy0, false);
    }
#endif

    for (int lev = 0; lev <= finest_level; lev++) {
        auto& mfx = *mapfac[lev][MapFacType::m_x];
        auto& mfy = *mapfac[lev][MapFacType::m_y];
        auto&  dJ = *detJ_cc[lev];
        rhth_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons], RhoTheta_comp,dJ,mfx,mfy,true);
        scal_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons],RhoScalar_comp,dJ,mfx,mfy,true);
        if (solverChoice.moisture_type != MoistureType::None) {
            int n_qstate_into_total = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
            for (int qoff(0); qoff<n_qstate_into_total; ++qoff) {
                mois_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons],RhoQ1_comp+qoff,dJ,mfx,mfy,true);
            }
        }
#ifdef ERF_USE_FIRE
        if (m_fire_coupler) {
            smoke_ml += volWgtSumMF(lev, vars_new[lev][Vars::cons], RhoSmoke_comp, dJ, mfx, mfy, true);
        }
#endif
    }

    Real h_avg_ustar = zero;
    Real h_avg_tstar = zero;
    Real h_avg_olen  = zero;
    Orientation zlo(Direction::z, Orientation::low);
    if ((m_SurfaceLayer[zlo] != nullptr) && (NumDataLogs() > 0)) {
        // Sum each surface cell once: on grids split in z the planar arrays hold one
        // duplicate box per stacked 3D box
        h_avg_ustar = m_SurfaceLayer[zlo]->surface_sum(0, *m_SurfaceLayer[zlo]->get_u_star(0));
        h_avg_tstar = m_SurfaceLayer[zlo]->surface_sum(0, *m_SurfaceLayer[zlo]->get_t_star(0));
        h_avg_olen  = m_SurfaceLayer[zlo]->surface_sum(0, *m_SurfaceLayer[zlo]->get_olen(0));

        // Divide by the total number of cells we are averaging over
        const Box& domain = geom[0].Domain();
        const Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
        h_avg_ustar /= area_z;
        h_avg_tstar /= area_z;
        h_avg_olen  /= area_z;
    }

    // On pseudo-incompressible levels, the departure of (rho theta) from (rho0 theta0), which the
    // constraint div(rho0 theta0 u) = 0 holds at zero without heating or diffusion: max over the
    // valid cells (with fluid) of all such levels of |rho theta - rho0 theta0| / (rho0 theta0)
    Real pi_drift = -one;
    for (int lev = 0; lev <= finest_level; lev++) {
        if (solverChoice.compressibility_type[lev] != CompressibilityType::PseudoIncompressible) { continue; }
        auto const& s_arr  = vars_new[lev][Vars::cons].const_arrays();
        auto const& b_arr  = base_state[lev].const_arrays();
        auto const& dJ_arr = detJ_cc[lev]->const_arrays();
        Real lev_drift = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                                   vars_new[lev][Vars::cons], IntVect(0),
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            if (dJ_arr[b](i,j,k) <= zero) { return { zero }; }
            const Real rt0 = b_arr[b](i,j,k,BaseState::r0_comp) * b_arr[b](i,j,k,BaseState::th0_comp);
            return { std::abs(s_arr[b](i,j,k,RhoTheta_comp) - rt0) / rt0 };
        });
        pi_drift = std::max(pi_drift, lev_drift);
    }
    ParallelDescriptor::ReduceRealMax(pi_drift);

#ifdef ERF_USE_FIRE
    const int nfoo = 10;
    Real foo[nfoo] = {mass_sl,rhth_sl,scal_sl,mois_sl,smoke_sl,mass_ml,rhth_ml,scal_ml,mois_ml,smoke_ml};
#else
    const int nfoo = 8;
    Real foo[nfoo] = {mass_sl,rhth_sl,scal_sl,mois_sl,mass_ml,rhth_ml,scal_ml,mois_ml};
#endif
    Real zero_d = zero;
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
    ParallelDescriptor::ReduceRealSum(
        foo, nfoo, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
        int i = 0;
        mass_sl = foo[i++];
        rhth_sl = foo[i++];
        scal_sl = foo[i++];
        mois_sl = foo[i++];
#ifdef ERF_USE_FIRE
        smoke_sl = foo[i++];
#endif
        mass_ml = foo[i++];
        rhth_ml = foo[i++];
        scal_ml = foo[i++];
        mois_ml = foo[i++];
#ifdef ERF_USE_FIRE
        smoke_ml = foo[i++];
#endif

        // Print the integrals with every digit, so that a run can be checked for conservation
        // from its output (Tests/RunConservation.cmake)
        const int sumprecision = std::numeric_limits<Real>::max_digits10;

        Print() << '\n';
        Print() << "TIME= " << std::setw(datwidth) << std::setprecision(timeprecision) << std::left << time << '\n';
        if (finest_level ==  0) {
#if 1
           Print() << std::setprecision(sumprecision) << " MASS       = " << mass_sl << '\n';
#else
           Print() << std::setprecision(sumprecision) << " PERT MASS  = " << mass_sl << '\n';
#endif
           Print() << std::setprecision(sumprecision) << " RHO THETA  = " << rhth_sl << '\n';
           if (solverChoice.transport_scalar) { Print() << std::setprecision(sumprecision) << " RHO SCALAR = " << scal_sl << '\n'; }
           if (solverChoice.moisture_type != MoistureType::None) { Print() << std::setprecision(sumprecision) << " RHO QTOTAL = " << mois_sl << '\n'; }
#ifdef ERF_USE_FIRE
           if (m_fire_coupler) { Print() << std::setprecision(sumprecision) << " RHO SMOKE  = " << smoke_sl << '\n'; }
#endif
        } else {
#if 1
           Print() << std::setprecision(sumprecision) << " MASS       SL/ML = " << mass_sl << " " << mass_ml << '\n';
#else
           Print() << std::setprecision(sumprecision) << " PERT MASS  SL/ML = " << mass_sl << " " << mass_ml << '\n';
#endif
           Print() << std::setprecision(sumprecision) << " RHO THETA  SL/ML = " << rhth_sl << " " << rhth_ml << '\n';
           if (solverChoice.transport_scalar) { Print() << std::setprecision(sumprecision) << " RHO SCALAR SL/ML = " << scal_sl << " " << scal_ml << '\n'; }
           if (solverChoice.moisture_type != MoistureType::None) { Print() << std::setprecision(sumprecision) << " RHO QTOTAL SL/ML = " << mois_sl << " " << mois_ml << '\n'; }
#ifdef ERF_USE_FIRE
           if (m_fire_coupler) { Print() << std::setprecision(sumprecision) << " RHO SMOKE  SL/ML = " << smoke_sl << " " << smoke_ml << '\n'; }
#endif
        }
        if (pi_drift >= zero) {
           Print() << std::setprecision(sumprecision) << " PI DRIFT   = " << pi_drift
                   << "  (max |rho theta - rho0 theta0| / (rho0 theta0))" << '\n';
        }

        // The first data log only holds scalars
        if (NumDataLogs() > 0)
        {
            int n_d = 0;
            std::ostream& data_log1 = DataLog(n_d);
            if (data_log1.good()) {
                if (time == zero_d) {
                    data_log1 << std::setw(datwidth) << "          time";
                    data_log1 << std::setw(datwidth) << "          u_star";
                    data_log1 << std::setw(datwidth) << "          t_star";
                    data_log1 << std::setw(datwidth) << "          olen";
                    data_log1 << std::endl;
                } // time = 0

              // Write the quantities at this time
              data_log1 << std::setw(datwidth) << std::setprecision(timeprecision) << time;
              data_log1 << std::setw(datwidth) << std::setprecision(datprecision)  << h_avg_ustar;
              data_log1 << std::setw(datwidth) << std::setprecision(datprecision)  << h_avg_tstar;
              data_log1 << std::setw(datwidth) << std::setprecision(datprecision)  << h_avg_olen;
              data_log1 << std::endl;
            } // if good
        } // loop over i
      } // if IOProcessor
#ifdef AMREX_LAZY
    });
#endif

    if (sample_point_sampler) {
        run_line_sampler(*sample_point_sampler);
    }
    if (sample_line_sampler) {
        run_line_sampler(*sample_line_sampler);
    }
}

void
ERF::sum_derived_quantities (double time)
{
    if (verbose <= 0 || NumDerDataLogs() <= 0) return;

    int lev = 0;

    AMREX_ALWAYS_ASSERT(lev == 0);

    auto& mfx0 = *mapfac[0][MapFacType::m_x];
    auto& mfy0 = *mapfac[0][MapFacType::m_y];
    auto&  dJ0 = *detJ_cc[0];

    // ************************************************************************
    // WARNING: we are not filling ghost cells other than periodic outside the domain
    // ************************************************************************

    MultiFab mf_cc_vel(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(1,1,1));
    mf_cc_vel.setVal(0.); // We just do this to avoid uninitialized values

    // Average all three components of velocity (on faces) to the cell center
    average_face_to_cellcenter(mf_cc_vel,0,
                               Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],
                                                        &vars_new[lev][Vars::yvel],
                                                        &vars_new[lev][Vars::zvel]});
    mf_cc_vel.FillBoundary(geom[lev].periodicity());

    if (!geom[lev].isPeriodic(0) || !geom[lev].isPeriodic(1) || !geom[lev].isPeriodic(2)) {
        amrex::Warning("Ghost cells outside non-periodic physical boundaries are not filled -- vel set to 0 there");
    }

    MultiFab r_wted_magvelsq(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));
    MultiFab unwted_magvelsq(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));
    MultiFab     enstrophysq(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(1,1,1));
    MultiFab        theta_mf(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(unwted_magvelsq, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto& src_fab = mf_cc_vel[mfi];

        auto& dest1_fab = unwted_magvelsq[mfi];
        // NOTE: we send in src_fab where we should
        derived::erf_dermagvelsq(bx, dest1_fab, 0, 1, src_fab, (*z_phys_cc[lev])[mfi], Geom(lev),
                                 static_cast<Real>(t_new[0]), nullptr, lev);

        auto& dest2_fab = enstrophysq[mfi];
        derived::erf_derenstrophysq(bx, dest2_fab, 0, 1, src_fab, (*z_phys_cc[lev])[mfi], Geom(lev),
                                    static_cast<Real>(t_new[0]), nullptr, lev);
    }

    // Copy the MF holding 1/2(u^2 + v^2 + w^2) into the MF that will hold 1/2 rho (u^2 + v^2 + w^2)d
    MultiFab::Copy(r_wted_magvelsq, unwted_magvelsq, 0, 0, 1, 0);

    // Multiply the MF holding 1/2(u^2 + v^2 + w^2) by rho to get  1/2 rho (u^2 + v^2 + w^2)
    MultiFab::Multiply(r_wted_magvelsq, vars_new[lev][Vars::cons], 0, 0, 1, 0);

    // Copy the MF holding (rho theta) into "theta_mf"
    MultiFab::Copy(theta_mf, vars_new[lev][Vars::cons], RhoTheta_comp, 0, 1, 0);

    // Divide (rho theta) by rho to get theta in the MF "theta_mf"
    MultiFab::Divide(theta_mf, vars_new[lev][Vars::cons], Rho_comp, 0, 1, 0);

    Real  unwted_avg = volWgtSumMF(lev, unwted_magvelsq, 0, dJ0, mfx0, mfy0, false);
    Real  r_wted_avg = volWgtSumMF(lev, r_wted_magvelsq, 0, dJ0, mfx0, mfy0, false);
    Real enstrsq_avg = volWgtSumMF(lev, enstrophysq,     0, dJ0, mfx0, mfy0, false);
    Real   theta_avg = volWgtSumMF(lev, theta_mf,        0, dJ0, mfx0, mfy0, false);

    // Get volume including terrain (consistent with volWgtSumMF routine)
    MultiFab volume(grids[lev], dmap[lev], 1, 0);
    auto const& dx = geom[lev].CellSizeArray();
    Real cell_vol  = dx[0]*dx[1]*dx[2];
    volume.setVal(cell_vol);
    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
        MultiFab::Multiply(volume, *detJ_cc[lev], 0, 0, 1, 0);
    }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(volume, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx  = mfi.tilebox();
        auto dst        = volume.array(mfi);
        const auto& mfx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const auto& mfy = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst(i,j,k) /= (mfx(i,j,0)*mfy(i,j,0));
        });
    }
    Real vol = volume.sum();

     unwted_avg /= vol;
     r_wted_avg /= vol;
    enstrsq_avg /= vol;
      theta_avg /= vol;

    const int nfoo = 4;
    Real foo[nfoo] = {unwted_avg,r_wted_avg,enstrsq_avg,theta_avg};
    Real zero_d = zero;
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
    ParallelDescriptor::ReduceRealSum(
        foo, nfoo, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
        int i = 0;
        unwted_avg  = foo[i++];
        r_wted_avg  = foo[i++];
        enstrsq_avg = foo[i++];
          theta_avg = foo[i++];

        std::ostream& data_log_der = DerDataLog(0);

        if (time == zero_d) {
            data_log_der << std::setw(datwidth) << "          time";
            data_log_der << std::setw(datwidth) << "        ke_den";
            data_log_der << std::setw(datwidth) << "         velsq";
            data_log_der << std::setw(datwidth) << "     enstrophy";
            data_log_der << std::setw(datwidth) << "    int_energy";
            data_log_der << std::endl;
        }
        data_log_der << std::setw(datwidth) << std::setprecision(timeprecision) << time;
        data_log_der << std::setw(datwidth) << std::setprecision(datprecision)  <<  unwted_avg;
        data_log_der << std::setw(datwidth) << std::setprecision(datprecision)  <<  r_wted_avg;
        data_log_der << std::setw(datwidth) << std::setprecision(datprecision)  << enstrsq_avg;
        data_log_der << std::setw(datwidth) << std::setprecision(datprecision)  <<   theta_avg;
        data_log_der << std::endl;

      } // if IOProcessor
#ifdef AMREX_LAZY
    });
#endif
}

void
ERF::sum_energy_quantities (double time)
{
    if ( (verbose <= 0) || (tot_e_datalog.size() < 1) ) { return; }

    int lev = 0;

    auto& mfx0 = *mapfac[0][MapFacType::m_x];
    auto& mfy0 = *mapfac[0][MapFacType::m_y];
    auto&  dJ0 = *detJ_cc[0];

    AMREX_ALWAYS_ASSERT(lev == 0);

    bool local = true;

    // ************************************************************************
    // WARNING: we are not filling ghost cells other than periodic outside the domain
    // ************************************************************************

    MultiFab mf_cc_vel(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(1,1,1));
    mf_cc_vel.setVal(0.); // We just do this to avoid uninitialized values

    // Average all three components of velocity (on faces) to the cell center
    average_face_to_cellcenter(mf_cc_vel,0,
                               Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],
                                                        &vars_new[lev][Vars::yvel],
                                                        &vars_new[lev][Vars::zvel]});
    mf_cc_vel.FillBoundary(geom[lev].periodicity());

    if (!geom[lev].isPeriodic(0) || !geom[lev].isPeriodic(1) || !geom[lev].isPeriodic(2)) {
        amrex::Warning("Ghost cells outside non-periodic physical boundaries are not filled -- vel set to 0 there");
    }

    MultiFab tot_mass  (grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));
    MultiFab tot_energy(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));

    auto const& dx = geom[lev].CellSizeArray();
    auto const  prob_lo = geom[lev].ProbLoArray();
    bool is_moist = (solverChoice.moisture_type != MoistureType::None);
    const MoistureComponentIndices moist_idx = solverChoice.moisture_indices;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tot_mass, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const Array4<Real>& cc_vel_arr     = mf_cc_vel.array(mfi);
        const Array4<Real>& tot_mass_arr   = tot_mass.array(mfi);
        const Array4<Real>& tot_energy_arr = tot_energy.array(mfi);
        const Array4<const Real>& cons_arr = vars_new[lev][Vars::cons].const_array(mfi);
        const Array4<const Real>& z_arr    = (z_phys_cc[lev]) ? z_phys_cc[lev]->const_array(mfi) :
                                                                Array4<const Real>{};
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Mixing ratios of vapour, of all liquid (cloud and rain) and of all ice (cloud
            // ice, snow, graupel); the state holds rho_d times each species
            Real Qv = zero, Ql = zero, Qi = zero;
            if (is_moist) { GetMoistureVars(i, j, k, cons_arr, Qv, Ql, Qi, moist_idx); }
            Real Qt   = Qv + Ql + Qi;
            Real Rhod = cons_arr(i,j,k,Rho_comp);
            Real Rhot = Rhod * (one + Qt);
            Real Temp = getTgivenRandRTh(Rhod, cons_arr(i,j,k,RhoTheta_comp), Qv);
            Real TKE  = myhalf * ( cc_vel_arr(i,j,k,0)*cc_vel_arr(i,j,k,0)
                              + cc_vel_arr(i,j,k,1)*cc_vel_arr(i,j,k,1)
                              + cc_vel_arr(i,j,k,2)*cc_vel_arr(i,j,k,2) );
            // Height of the cell centre
            Real zval = (z_arr) ? z_arr(i,j,k) : prob_lo[2] + (Real(k) + myhalf)*dx[2];

            // Internal energy with constant heat capacities (Romps 2008, J. Atmos. Sci. 65,
            // 3779), measured from liquid water at the triple point: vapour carries
            // E0v = L_v - R_v T_trip, ice -E0s with E0s = L_f. The latent heats this
            // implies follow Kirchhoff's law, e.g. L_v(T) = L_v + (c_pv - c_l)(T - T_trip).
            Real Cv   = Cp_d - R_d;
            Real Cvv  = Cp_v - R_v;
            Real E0v  = L_v - R_v*h2otrip;
            Real E0s  = lfus;

            tot_mass_arr(i,j,k)   = Rhot;
            tot_energy_arr(i,j,k) = Rhod * ( (Cv + Cvv*Qv + Cp_l*Ql + Cp_i*Qi)*(Temp - h2otrip)
                                           + Qv*E0v - Qi*E0s
                                           + (one + Qt)*TKE + (one + Qt)*CONST_GRAV*zval );

        });

    }

    Real  tot_mass_avg   = volWgtSumMF(lev, tot_mass  , 0, dJ0, mfx0, mfy0, false, local);
    Real  tot_energy_avg = volWgtSumMF(lev, tot_energy, 0, dJ0, mfx0, mfy0, false, local);

    // Get volume including terrain (consistent with volWgtSumMF routine)
    MultiFab volume(grids[lev], dmap[lev], 1, 0);
    Real cell_vol  = dx[0]*dx[1]*dx[2];
    volume.setVal(cell_vol);
    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
        MultiFab::Multiply(volume, *detJ_cc[lev], 0, 0, 1, 0);
    }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(volume, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx  = mfi.tilebox();
        auto dst        = volume.array(mfi);
        const auto& mfx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const auto& mfy = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst(i,j,k) /= (mfx(i,j,0)*mfy(i,j,0));
        });
    }
    Real vol = volume.sum();

    // Divide by the volume
    tot_mass_avg   /= vol;
    tot_energy_avg /= vol;

    const int nfoo = 2;
    Real foo[nfoo] = {tot_mass_avg,tot_energy_avg};
    Real zero_d = zero;
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
    ParallelDescriptor::ReduceRealSum(
        foo, nfoo, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
        int i = 0;
        tot_mass_avg   = foo[i++];
        tot_energy_avg = foo[i++];

        std::ostream& data_log_energy = *tot_e_datalog[0];

        if (time == zero_d) {
            data_log_energy << std::setw(datwidth) << "          time";
            data_log_energy << std::setw(datwidth) << "      tot_mass";
            data_log_energy << std::setw(datwidth) << "    tot_energy";
            data_log_energy << std::endl;
        }
        data_log_energy << std::setw(datwidth) << std::setprecision(timeprecision) << time;
        data_log_energy << std::setw(datwidth) << std::setprecision(datprecision)  << tot_mass_avg;
        data_log_energy << std::setw(datwidth) << std::setprecision(datprecision)  << tot_energy_avg;
        data_log_energy << std::endl;

      } // if IOProcessor
#ifdef AMREX_LAZY
    });
#endif
}

Real
ERF::cloud_fraction (double /*time*/)
{
    BL_PROFILE("ERF::cloud_fraction()");

    int lev = 0;
    // This holds all of qc
    MultiFab qc(vars_new[lev][Vars::cons],make_alias,RhoQ2_comp,1);

    int direction = 2; // z-direction
    Box const& domain = geom[lev].Domain();

    auto const& qc_arr = qc.const_arrays();

    // qc_2d is an BaseFab<int> holding the max value over the column
    auto qc_2d = ReduceToPlane<ReduceOpMax,int>(direction, domain, qc,
         [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> int
         {
             if (qc_arr[box_no](i,j,k) > 0) {
                 return 1;
             } else {
                 return 0;
             }
         });

    auto* p = qc_2d.dataPtr();

    Long numpts = qc_2d.numPts();

    AMREX_ASSERT(numpts < Long(std::numeric_limits<int>::max));

#if 1
    if (ParallelDescriptor::UseGpuAwareMpi()) {
        ParallelDescriptor::ReduceIntMax(p,static_cast<int>(numpts));
    } else {
        Gpu::PinnedVector<int> hv(numpts);
        Gpu::copyAsync(Gpu::deviceToHost, p, p+numpts, hv.data());
        Gpu::streamSynchronize();
        ParallelDescriptor::ReduceIntMax(hv.data(),static_cast<int>(numpts));
        Gpu::copyAsync(Gpu::hostToDevice, hv.data(), hv.data()+numpts, p);
        // We must synchronize before hv goes out of scope and its pinned memory
        //    is returned to the arena, since the copy above may still be reading it
        Gpu::streamSynchronize();
    }

    // Sum over component 0
    Long num_cloudy = qc_2d.template sum<RunOn::Device>(0);

#else
    //
    // We need this if we allow domain decomposition in the vertical
    //    but for now we leave it commented out
    //
    Long num_cloudy = Reduce::Sum<Long>(numpts,
         [=] AMREX_GPU_DEVICE (Long i) -> Long {
             if (p[i] == 1) {
                 return 1;
             } else {
                 return 0;
             }
         });
    ParallelDescriptor::ReduceLongSum(num_cloudy);
#endif

    Real num_total = Real(qc_2d.box().d_numPts());

    Real cloud_frac = num_cloudy / num_total;

    return cloud_frac;
}

/**
 * The volume fraction of each level, which line and plane sampling write as volfrac and use
 * to blank the cells covered by immersed terrain (erf.terrain_type = EB or ImmersedForcing);
 * nullptr otherwise, and on the levels up to max_level that do not exist, whose EB factory
 * is not built
 */
Vector<const MultiFab*>
ERF::sampling_volfrac () const
{
    Vector<const MultiFab*> volfrac_lev(vars_new.size(), nullptr);
    if (solverChoice.terrain_type == TerrainType::EB ||
        solverChoice.terrain_type == TerrainType::ImmersedForcing) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            volfrac_lev[lev] = &EBFactory(lev).getVolFrac();
        }
    }
    return volfrac_lev;
}

/**
 * The ground of each level that exists, the bottom of its mesh or, with immersed terrain, the
 * terrain surface (station_ground), for a line sampler that needs the ground; nullptr otherwise
 */
Vector<const MultiFab*>
ERF::sampling_ground (const LineSampler& sampler)
{
    Vector<const MultiFab*> ground_lev(vars_new.size(), nullptr);
    if (sampler.needs_ground()) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            ground_lev[lev] = &station_ground(lev);
        }
    }
    return ground_lev;
}

/**
 * Sample the lines of a LineSampler and write them, passing it the stresses, the SGS heat
 * flux, whether each level carries a prognostic TKE, the heights of each level and the ground
 * on level 0.  Used for erf.do_line_sampling, for erf.sample_line, which is a LineSampler
 * of full-depth columns on level 0, and for erf.sample_point, a LineSampler of single cells.
 *
 * @param sampler the line sampler
 */
void
ERF::run_line_sampler (LineSampler& sampler)
{
    const int nlev = static_cast<int>(vars_new.size());
    Vector<Vector<MultiFab*>> tau_lev(nlev);
    Vector<MultiFab*> hfx3_lev(nlev, nullptr);
    Vector<int> prognostic_tke_available(nlev, 0);
    Vector<MultiFab*> z_cc_lev(nlev, nullptr);
    Vector<MultiFab*> z_nd_lev(nlev, nullptr);
    for (int lev = 0; lev < nlev; ++lev) {
        prognostic_tke_available[lev] = solverChoice.turbChoice[lev].use_tke ? 1 : 0;
        if (lev < static_cast<int>(z_phys_cc.size())) { z_cc_lev[lev] = z_phys_cc[lev].get(); }
        if (lev < static_cast<int>(z_phys_nd.size())) { z_nd_lev[lev] = z_phys_nd[lev].get(); }
        if (lev < static_cast<int>(Tau.size())) {
            for (const auto& tau : Tau[lev]) {
                tau_lev[lev].push_back(tau.get());
            }
        }
        if (lev < static_cast<int>(SFS_hfx3_lev.size())) {
            hfx3_lev[lev] = SFS_hfx3_lev[lev].get();
        }
    }
    // With immersed terrain the ground is the terrain surface, not the bottom of the mesh
    const MultiFab* ground_nd = sampler.needs_ground() ? &station_ground(0) : nullptr;
    const Vector<const MultiFab*> ground_lev = sampling_ground(sampler);
    const Vector<const MultiFab*> volfrac_lev = sampling_volfrac();
    // The line coordinates of a level made or regridded since they were last written
    sampler.write_coords(z_phys_cc, z_phys_nd, geom, ground_nd, ground_lev, volfrac_lev);
    sampler.get_sample_data(geom, vars_new, tau_lev, hfx3_lev, prognostic_tke_available,
                            z_cc_lev, z_nd_lev, ground_nd, volfrac_lev, ground_lev);
    sampler.write_sample_data(t_new, istep, ref_ratio, geom);
}

/**
 * Helper function which uses the current step number, time, and timestep to
 * determine whether it is time to take an action specified at every interval
 * of timesteps.
 *
 * @param nstep Timestep number
 * @param time Current time
 * @param dtlev Timestep for the current level
 * @param action_interval Interval in number of timesteps for taking action
 * @param action_per Interval in simulation time for taking action
 */
bool
ERF::is_it_time_for_action (int nstep, double time, double dtlev, int action_interval, double action_per)
{
  bool int_test = (action_interval > 0 && nstep % action_interval == 0);

  bool per_test = false;
  if (action_per > zero) {
    const int num_per_old = static_cast<int>(amrex::Math::floor((time - dtlev) / action_per));
    const int num_per_new = static_cast<int>(amrex::Math::floor((time) / action_per));

    if (num_per_old != num_per_new) {
      per_test = true;
    }
  }

  return int_test || per_test;
}
